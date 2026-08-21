from __future__ import annotations

import math
import random
import statistics

import torch
from torch.utils.data import DataLoader

from .config import ModelConfig
from .training import move_batch


def _pearson(predictions: torch.Tensor, labels: torch.Tensor) -> float | None:
    if predictions.numel() < 2:
        return None
    x = predictions - predictions.mean()
    y = labels - labels.mean()
    denominator = torch.sqrt(x.square().sum() * y.square().sum())
    if denominator <= 0:
        return None
    return float((x * y).sum() / denominator)


def regression_metrics(predictions: torch.Tensor, labels: torch.Tensor,
                       standard_deviation: torch.Tensor | None = None) -> dict:
    if predictions.numel() == 0:
        return {"count": 0}
    error = predictions - labels
    metrics = {
        "count": int(predictions.numel()),
        "mae": float(error.abs().mean()),
        "rmse": float(torch.sqrt(error.square().mean())),
        "pearson": _pearson(predictions, labels),
    }
    if standard_deviation is not None:
        absolute = error.abs()
        metrics["coverage_1sigma"] = float((absolute <= standard_deviation).float().mean())
        metrics["coverage_2sigma"] = float((absolute <= 2 * standard_deviation).float().mean())
    return metrics


def scalar_uncertainty_scale(
    predictions: torch.Tensor, standard_deviation: torch.Tensor, labels: torch.Tensor,
) -> float:
    if predictions.shape != labels.shape or standard_deviation.shape != labels.shape:
        raise ValueError("predictions, deviations, and labels must have identical shapes")
    if predictions.numel() == 0:
        raise ValueError("uncertainty calibration requires observations")
    if (not torch.isfinite(predictions).all() or not torch.isfinite(labels).all()
            or not torch.isfinite(standard_deviation).all()
            or (standard_deviation <= 0).any()):
        raise ValueError("uncertainty calibration inputs must be finite with positive deviations")
    standardized_square = ((labels - predictions) / standard_deviation).square()
    return float(torch.sqrt(standardized_square.mean()))


def uncertainty_metrics(
    predictions: torch.Tensor, standard_deviation: torch.Tensor, labels: torch.Tensor,
) -> dict:
    if predictions.shape != labels.shape or standard_deviation.shape != labels.shape:
        raise ValueError("predictions, deviations, and labels must have identical shapes")
    error = predictions - labels
    nll = 0.5 * (
        (error / standard_deviation).square()
        + 2.0 * torch.log(standard_deviation)
        + math.log(2.0 * math.pi)
    )
    return {
        "count": int(predictions.numel()),
        "coverage_1sigma": float((error.abs() <= standard_deviation).float().mean()),
        "coverage_2sigma": float((error.abs() <= 2.0 * standard_deviation).float().mean()),
        "gaussian_nll": float(nll.mean()),
    }


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end - 1) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def _fractional_top_membership(values: list[float], count: int) -> list[float]:
    membership = [0.0] * len(values)
    order = sorted(range(len(values)), key=values.__getitem__, reverse=True)
    start = 0
    remaining = count
    while start < len(order) and remaining > 0:
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        group = order[start:end]
        fraction = min(1.0, remaining / len(group))
        for index in group:
            membership[index] = fraction
        remaining -= min(remaining, len(group))
        start = end
    return membership


def campaign_ranking_metrics(predictions: list[float], observations: list[float], *,
                             seed: int = 42, max_pairs: int = 100_000) -> dict:
    if len(predictions) != len(observations):
        raise ValueError("predictions and observations must have equal length")
    if not predictions:
        return {"count": 0}
    predicted_ranks = torch.tensor(_average_ranks(predictions), dtype=torch.float64)
    observed_ranks = torch.tensor(_average_ranks(observations), dtype=torch.float64)
    spearman = _pearson(predicted_ranks, observed_ranks)
    count = len(predictions)
    top_count = max(1, math.ceil(count * 0.1))
    predicted_top = _fractional_top_membership(predictions, top_count)
    observed_top = _fractional_top_membership(observations, top_count)
    overlap = sum(predicted * observed for predicted, observed in zip(predicted_top, observed_top))

    all_pair_count = count * (count - 1) // 2
    rng = random.Random(seed)
    if all_pair_count <= max_pairs:
        pairs = ((left, right) for left in range(count) for right in range(left + 1, count))
    else:
        sampled = set()
        while len(sampled) < max_pairs:
            left, right = rng.sample(range(count), 2)
            sampled.add((min(left, right), max(left, right)))
        pairs = iter(sampled)
    concordant = 0.0
    comparable = 0
    for left, right in pairs:
        observed_delta = observations[left] - observations[right]
        predicted_delta = predictions[left] - predictions[right]
        if observed_delta == 0:
            continue
        comparable += 1
        if predicted_delta == 0:
            concordant += 0.5
        else:
            concordant += (observed_delta > 0) == (predicted_delta > 0)
    return {
        "count": count,
        "spearman": spearman,
        "pairwise_accuracy": concordant / comparable if comparable else None,
        "pairwise_comparisons": comparable,
        "top_10pct_count": top_count,
        "top_10pct_overlap": overlap,
        "top_10pct_recall": overlap / top_count,
        "top_10pct_enrichment": overlap * count / (top_count * top_count),
    }


def summarize_campaign_runs(runs: list[dict], *, bootstrap_samples: int = 10_000,
                            seed: int = 42) -> dict:
    if not runs:
        raise ValueError("at least one campaign evaluation run is required")
    if bootstrap_samples < 1:
        raise ValueError("bootstrap_samples must be positive")
    selection = runs[0].get("selection_identity")
    campaigns = tuple(sorted(runs[0].get("campaign_metrics", {})))
    if not campaigns:
        raise ValueError("campaign evaluation runs must contain campaign_metrics")
    for run in runs:
        if run.get("selection_identity") != selection:
            raise ValueError("campaign evaluation runs use different frozen selections")
        if tuple(sorted(run.get("campaign_metrics", {}))) != campaigns:
            raise ValueError("campaign evaluation runs contain different campaigns")
    metric_names = (
        "spearman", "pairwise_accuracy", "top_10pct_recall", "top_10pct_enrichment")
    per_seed = {name: [run["macro_metrics"][name] for run in runs] for name in metric_names}
    rng = random.Random(seed)
    summary = {}
    for name in metric_names:
        defined_seed = [value for value in per_seed[name] if value is not None]
        campaign_values = {
            campaign: [run["campaign_metrics"][campaign][name] for run in runs]
            for campaign in campaigns
        }
        defined_campaigns = [
            campaign for campaign, values in campaign_values.items()
            if any(value is not None for value in values)
        ]
        campaign_means = {
            campaign: statistics.mean(
                value for value in values if value is not None)
            for campaign, values in campaign_values.items()
            if any(value is not None for value in values)
        }
        if defined_campaigns:
            bootstrap = []
            for _ in range(bootstrap_samples):
                sample = [
                    defined_campaigns[rng.randrange(len(defined_campaigns))]
                    for _ in defined_campaigns
                ]
                bootstrap.append(statistics.mean(campaign_means[c] for c in sample))
            bootstrap.sort()
            ci = [bootstrap[int(0.025 * (bootstrap_samples - 1))],
                  bootstrap[int(0.975 * (bootstrap_samples - 1))]]
        else:
            ci = None
        summary[name] = {
            "seed_values": per_seed[name],
            "seed_mean": statistics.mean(defined_seed) if defined_seed else None,
            "seed_sample_stddev": (
                statistics.stdev(defined_seed) if len(defined_seed) > 1 else None),
            "campaign_bootstrap_mean": (
                statistics.mean(campaign_means.values()) if campaign_means else None),
            "campaign_bootstrap_95pct_ci": ci,
        }
    return {
        "run_count": len(runs),
        "campaign_count": len(campaigns),
        "selection_identity": selection,
        "bootstrap_samples": bootstrap_samples,
        "bootstrap_seed": seed,
        "metrics": summary,
    }


def summarize_paired_regression_runs(
    task_query: dict[int, dict], late_concat: dict[int, dict], *,
    metric_names: tuple[str, ...] = (
        "rmse", "mae", "pearson", "coverage_1sigma", "coverage_2sigma"),
    left_name: str = "task_query",
    right_name: str = "late_concat",
    difference_key: str = "late_minus_task_by_seed",
) -> dict:
    seeds = sorted(task_query)
    if not seeds or seeds != sorted(late_concat):
        raise ValueError("fusion runs must contain the same non-empty seed set")
    if not metric_names:
        raise ValueError("at least one paired metric is required")
    for seed in seeds:
        if any(name not in task_query[seed] or name not in late_concat[seed]
               for name in metric_names):
            raise ValueError(f"fusion run seed {seed} is missing required metrics")

    modes = {left_name: task_query, right_name: late_concat}
    summaries = {}
    for mode, runs in modes.items():
        summaries[mode] = {}
        for name in metric_names:
            values = [runs[seed][name] for seed in seeds]
            summaries[mode][name] = {
                "seed_values": values,
                "mean": statistics.mean(values),
                "sample_stddev": statistics.stdev(values) if len(values) > 1 else None,
            }
    paired = {}
    for name in metric_names:
        values = [late_concat[seed][name] - task_query[seed][name] for seed in seeds]
        paired[name] = {
            difference_key: values,
            "mean": statistics.mean(values),
            "sample_stddev": statistics.stdev(values) if len(values) > 1 else None,
        }
    return {"seeds": seeds, "modes": summaries, "paired_differences": paired}


def summarize_regression_seed_runs(runs: dict[str, dict[int, dict]]) -> dict:
    if not runs:
        raise ValueError("at least one regression baseline is required")
    metric_names = ("rmse", "mae", "pearson")
    expected_seeds = None
    output = {}
    for name, seed_runs in sorted(runs.items()):
        seeds = sorted(seed_runs)
        if not seeds:
            raise ValueError(f"regression baseline has no seeds: {name}")
        if expected_seeds is None:
            expected_seeds = seeds
        elif seeds != expected_seeds:
            raise ValueError("regression baselines must contain the same seed set")
        output[name] = {}
        for metric_name in metric_names:
            if any(metric_name not in seed_runs[seed] for seed in seeds):
                raise ValueError(f"regression baseline {name} is missing {metric_name}")
            values = [seed_runs[seed][metric_name] for seed in seeds]
            output[name][metric_name] = {
                "seed_values": values,
                "mean": statistics.mean(values),
                "sample_stddev": statistics.stdev(values) if len(values) > 1 else None,
            }
    return {"seeds": expected_seeds, "baselines": output}


@torch.no_grad()
def evaluate_model(model, loader: DataLoader, device: torch.device,
                   config: ModelConfig) -> dict:
    model.eval()
    predictions = []
    deviations = []
    labels = []
    masks = []
    for batch in loader:
        batch = move_batch(batch, device)
        output = model(batch)
        predictions.append(output["mean"].cpu())
        deviations.append(output["standard_deviation"].cpu())
        labels.append(batch["labels"].cpu())
        masks.append(batch["label_mask"].cpu())
    prediction_tensor = torch.cat(predictions)
    deviation_tensor = torch.cat(deviations)
    label_tensor = torch.cat(labels)
    mask_tensor = torch.cat(masks)
    return {
        task: regression_metrics(
            prediction_tensor[:, index][mask_tensor[:, index]],
            label_tensor[:, index][mask_tensor[:, index]],
            deviation_tensor[:, index][mask_tensor[:, index]],
        )
        for index, task in enumerate(config.task_names)
    }


def mean_baseline_metrics(train_records, evaluation_records,
                          task_names: tuple[str, ...]) -> dict:
    output = {}
    for task in task_names:
        train_values = [getattr(record, task) for record in train_records
                        if getattr(record, task) is not None]
        eval_values = [getattr(record, task) for record in evaluation_records
                       if getattr(record, task) is not None]
        if not train_values or not eval_values:
            output[task] = {"count": 0}
            continue
        mean = sum(train_values) / len(train_values)
        predictions = torch.full((len(eval_values),), mean)
        output[task] = regression_metrics(predictions, torch.tensor(eval_values))
    return output
