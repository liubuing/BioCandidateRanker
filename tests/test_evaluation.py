import unittest

import torch

from biocandidate.evaluation import (
    scalar_uncertainty_scale,
    summarize_campaign_runs,
    summarize_paired_regression_runs,
    summarize_regression_seed_runs,
    uncertainty_metrics,
)


def run(value, selection="same"):
    return {
        "selection_identity": selection,
        "macro_metrics": {
            "spearman": value,
            "pairwise_accuracy": value,
            "top_10pct_recall": value,
            "top_10pct_enrichment": value,
        },
        "campaign_metrics": {
            campaign: {
                "spearman": value + offset,
                "pairwise_accuracy": value + offset,
                "top_10pct_recall": value + offset,
                "top_10pct_enrichment": value + offset,
            }
            for campaign, offset in (("a", 0.0), ("b", 0.2))
        },
    }


class CampaignRunSummaryTest(unittest.TestCase):
    def test_summarizes_seed_and_campaign_uncertainty_deterministically(self):
        summary = summarize_campaign_runs(
            [run(0.1), run(0.3)], bootstrap_samples=100, seed=7)
        metric = summary["metrics"]["spearman"]
        self.assertAlmostEqual(metric["seed_mean"], 0.2)
        self.assertAlmostEqual(metric["campaign_bootstrap_mean"], 0.3)
        self.assertEqual(summary["campaign_count"], 2)
        self.assertEqual(
            summary,
            summarize_campaign_runs([run(0.1), run(0.3)], bootstrap_samples=100, seed=7),
        )

    def test_rejects_mixed_frozen_selections(self):
        with self.assertRaisesRegex(ValueError, "different frozen selections"):
            summarize_campaign_runs([run(0.1), run(0.2, selection="other")])

    def test_summarizes_paired_fusion_differences(self):
        names = ("rmse", "mae", "pearson", "coverage_1sigma", "coverage_2sigma")
        task = {
            seed: {name: float(seed + index) for index, name in enumerate(names)}
            for seed in (1, 2, 3)
        }
        late = {
            seed: {name: value + 0.25 for name, value in metrics.items()}
            for seed, metrics in task.items()
        }
        summary = summarize_paired_regression_runs(task, late)
        self.assertEqual(summary["seeds"], [1, 2, 3])
        self.assertAlmostEqual(summary["paired_differences"]["rmse"]["mean"], 0.25)

        point_summary = summarize_paired_regression_runs(
            task, late, metric_names=("rmse", "mae", "pearson"))
        self.assertNotIn("coverage_1sigma", point_summary["modes"]["task_query"])

    def test_summarizes_regression_baseline_seeds(self):
        runs = {
            "a": {
                seed: {"rmse": float(seed), "mae": 1.0, "pearson": 0.2}
                for seed in (1, 2, 3)
            },
            "b": {
                seed: {"rmse": float(seed + 1), "mae": 2.0, "pearson": 0.3}
                for seed in (1, 2, 3)
            },
        }
        summary = summarize_regression_seed_runs(runs)
        self.assertEqual(summary["seeds"], [1, 2, 3])
        self.assertAlmostEqual(summary["baselines"]["a"]["rmse"]["mean"], 2.0)

    def test_scalar_uncertainty_calibration_matches_standardized_rmse(self):
        predictions = torch.tensor([0.0, 0.0])
        deviations = torch.tensor([1.0, 2.0])
        labels = torch.tensor([2.0, 4.0])
        scale = scalar_uncertainty_scale(predictions, deviations, labels)
        self.assertAlmostEqual(scale, 2.0)
        metrics = uncertainty_metrics(predictions, deviations * scale, labels)
        self.assertAlmostEqual(metrics["coverage_1sigma"], 1.0)


if __name__ == "__main__":
    unittest.main()
