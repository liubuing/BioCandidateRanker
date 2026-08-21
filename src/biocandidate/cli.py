from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

import torch
from torch.utils.data import DataLoader

from .calibration import (
    AffineVarianceGaussianCalibrator,
    IdentityGaussianCalibrator,
    ScalarGaussianCalibrator,
    SplitConformalCalibrator,
    conformal_calibrator_from_dict,
    load_calibration_artifact,
    save_calibration_artifact,
    select_grouped_gaussian_calibrator,
)
from .config import ModelConfig
from .baselines import amino_acid_composition, morgan_fingerprints
from .feature_mlp import save_feature_checkpoint, train_feature_mlp
from .data import (
    EnzymeSubstrateCollator,
    ENZENGDB_ARCHIVE_SHA256,
    ENZENGDB_ARCHIVE_SIZE_BYTES,
    ENZENGDB_DOI,
    ENZENGDB_LICENSE,
    LUNZER_ARTICLE_DOI,
    LUNZER_DATASET_DOI,
    LUNZER_DATASET_LICENSE,
    LUNZER_SOURCE_SHA256,
    LUNZER_SOURCE_SIZE_BYTES,
    CANONICAL_KINETIC_UNITS,
    campaign_rank_records,
    campaign_representative_records,
    FileManifest,
    build_manifest,
    apply_split_manifest,
    assign_homology_splits,
    read_unikp_json,
    read_normalized_kinetics_jsonl,
    read_lunzer_tsv,
    read_enzengdb_directory,
    read_mmseqs_clusters,
    run_mmseqs_easy_cluster,
    run_mmseqs_easy_search,
    sequence_id,
    split_records,
    verify_manifest,
    write_split_manifest,
    write_unique_fasta,
    write_cold_split_manifest,
    aggregate_pair_measurements,
)
from .model.ranker import BioCandidateRanker
from .evaluation import (
    campaign_ranking_metrics,
    evaluate_model,
    mean_baseline_metrics,
    regression_metrics,
    scalar_uncertainty_scale,
    summarize_campaign_runs,
    summarize_paired_regression_runs,
    summarize_regression_seed_runs,
    uncertainty_metrics,
)
from .training import (
    CampaignBatchSampler,
    CampaignGroupedBatchSampler,
    RecordDataset,
    create_lr_scheduler,
    run_epoch,
    run_epoch_accum,
    run_ranking_epoch,
    save_checkpoint,
)
from .training import load_checkpoint


DATA_FORMATS = ("unikp", "normalized-kinetics")
KINETIC_TASKS = tuple(CANONICAL_KINETIC_UNITS)


def _manifest_dict(manifest: FileManifest) -> dict:
    return asdict(manifest)


def _binary_identity(path: str | Path) -> dict:
    payload = Path(path).read_bytes()
    return {"sha256": hashlib.sha256(payload).hexdigest(), "size_bytes": len(payload)}


def _calibration_identities(checkpoint: str | Path, task: str) -> dict:
    return {"checkpoint": _binary_identity(checkpoint), "task": task}


def audit_command(args: argparse.Namespace) -> None:
    manifest = build_manifest(args.data)
    result = read_unikp_json(args.data)
    output = {
        "source": str(Path(args.data).resolve()),
        "identity": _manifest_dict(manifest),
        "adapter": {
            "total_rows": result.audit.total_rows,
            "accepted_rows": result.audit.accepted_rows,
            "rejected_rows": len(result.audit.rejected),
            "duplicate_rows": len(result.audit.duplicates),
            "conflicting_duplicate_rows": result.audit.conflict_count,
        },
        "limitations": [
            "Source rows lack record-level accession provenance.",
            "No frozen independent benchmark is bundled.",
            "Exact protein-cold split is not a homology-cold split.",
        ],
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def audit_lunzer_command(args: argparse.Namespace) -> None:
    source_identity = build_manifest(args.data)
    if (source_identity.sha256 != LUNZER_SOURCE_SHA256
            or source_identity.size_bytes != LUNZER_SOURCE_SIZE_BYTES):
        raise ValueError("Lunzer source does not match the frozen Dryad file identity")
    observations = read_lunzer_tsv(args.data)
    records = [item.record for item in observations]
    overlap = None
    if args.unikp_data or args.unikp_manifest:
        if not args.unikp_data or not args.unikp_manifest:
            raise ValueError("--unikp-data and --unikp-manifest must be supplied together")
        manifest_payload = json.loads(Path(args.unikp_manifest).read_text(encoding="utf-8"))
        verify_manifest(args.unikp_data, FileManifest(**manifest_payload["identity"]))
        reference = read_unikp_json(args.unikp_data).records
        sequences = {record.sequence for record in reference}
        substrates = {record.substrate_smiles for record in reference}
        pairs = {(record.sequence, record.substrate_smiles) for record in reference}
        overlap = {
            "exact_sequence_records": sum(record.sequence in sequences for record in records),
            "exact_substrate_records": sum(record.substrate_smiles in substrates for record in records),
            "exact_sequence_substrate_pairs": sum(
                (record.sequence, record.substrate_smiles) in pairs for record in records),
            "publication_overlap": "unavailable_no_record_level_unikp_citations",
        }
    output = {
        "dataset": "Lunzer ancient adaptive landscape",
        "dataset_doi": LUNZER_DATASET_DOI,
        "article_doi": LUNZER_ARTICLE_DOI,
        "license": LUNZER_DATASET_LICENSE,
        "source_identity": _manifest_dict(source_identity),
        "records": len(observations),
        "genotypes": len({item.genotype for item in observations}),
        "cofactor_counts": dict(Counter(item.cofactor for item in observations)),
        "candidate_ids": [item.record.candidate_id for item in observations],
        "unikp_exact_overlap": overlap,
        "frozen_protocol": {
            "status": "confirmatory_no_model_selection",
            "primary_metrics": [
                "within_cofactor_spearman",
                "within_cofactor_pairwise_accuracy",
                "within_cofactor_top_10pct_enrichment",
            ],
            "endpoint": "derived_ln_kcat_equals_ln_km_plus_ln_kcat_over_km",
            "random_expectations": {
                "spearman": 0.0,
                "pairwise_accuracy": 0.5,
                "top_10pct_enrichment": 1.0,
            },
        },
        "limitations": [
            "All variants share one E. coli IMDH scaffold; this is not homology-cold.",
            "The deposited values are model-fit kinetics rather than independent raw replicates.",
            "An unresolved source-unit constant prevents absolute RMSE claims; ranking is invariant.",
            "The UniKP corpus lacks record-level citations, so publication overlap is unavailable.",
        ],
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def _select_enzengdb_observations(observations, *, min_campaign_size: int,
                                  max_records_per_campaign: int) -> tuple[list, dict]:
    campaigns = defaultdict(list)
    for observation in observations:
        campaigns[observation.campaign_id].append(observation)
    selected = []
    counts = {}
    for campaign_id, values in sorted(campaigns.items()):
        if len(values) < min_campaign_size:
            continue
        values = sorted(
            values,
            key=lambda item: hashlib.sha256(
                item.record.candidate_id.encode("utf-8")).hexdigest(),
        )
        if max_records_per_campaign > 0:
            values = values[:max_records_per_campaign]
        selected.extend(values)
        counts[campaign_id] = len(values)
    return selected, counts


def audit_enzengdb_command(args: argparse.Namespace) -> None:
    if args.exclude_homology_hits and not args.mmseqs:
        raise ValueError("--exclude-homology-hits requires --mmseqs")
    archive_identity = build_manifest(args.archive)
    if (archive_identity.sha256 != ENZENGDB_ARCHIVE_SHA256
            or archive_identity.size_bytes != ENZENGDB_ARCHIVE_SIZE_BYTES):
        raise ValueError("EnzEngDB archive does not match the frozen Zenodo v1 identity")
    result = read_enzengdb_directory(args.experiments)
    selected, campaign_counts = _select_enzengdb_observations(
        result.observations,
        min_campaign_size=args.min_campaign_size,
        max_records_per_campaign=args.max_records_per_campaign,
    )
    rejection_reasons = Counter(item.reason for item in result.audit.rejected)
    exact_overlap = None
    homology_overlap = None
    if args.unikp_data or args.unikp_manifest:
        if not args.unikp_data or not args.unikp_manifest:
            raise ValueError("--unikp-data and --unikp-manifest must be supplied together")
        manifest_payload = json.loads(Path(args.unikp_manifest).read_text(encoding="utf-8"))
        verify_manifest(args.unikp_data, FileManifest(**manifest_payload["identity"]))
        unikp_records = read_unikp_json(args.unikp_data).records
        sequences = {record.sequence for record in unikp_records}
        substrates = {record.substrate_smiles for record in unikp_records}
        pairs = {(record.sequence, record.substrate_smiles) for record in unikp_records}
        exact_overlap = {
            "sequence": sum(item.record.sequence in sequences for item in selected),
            "substrate": sum(item.record.substrate_smiles in substrates for item in selected),
            "sequence_substrate_pair": sum(
                (item.record.sequence, item.record.substrate_smiles) in pairs for item in selected),
            "unikp_publication_overlap": "unavailable_no_record_level_citations",
        }
        if args.mmseqs:
            homology_dir = Path(args.homology_output_dir)
            query_fasta = homology_dir / "enzengdb_selected.fasta"
            target_fasta = homology_dir / "unikp_reference.fasta"
            write_unique_fasta([item.record for item in selected], query_fasta)
            write_unique_fasta(unikp_records, target_fasta)
            hit_path = run_mmseqs_easy_search(
                args.mmseqs, query_fasta, target_fasta, homology_dir,
                min_identity=args.min_identity, coverage=args.coverage, threads=args.threads)
            hit_queries = set()
            with hit_path.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if line.strip():
                        hit_queries.add(line.split("\t", 1)[0])
            homology_overlap = {
                "method": "MMseqs2 easy-search",
                "min_identity": args.min_identity,
                "coverage": args.coverage,
                "coverage_mode": 0,
                "unique_external_sequences": len({
                    item.record.sequence for item in selected}),
                "unique_sequences_with_hit": len(hit_queries),
                "selected_records_with_hit": sum(
                    sequence_id(item.record.sequence) in hit_queries for item in selected),
                "hits_path": str(hit_path.resolve()),
                "hits_identity": _manifest_dict(build_manifest(hit_path)),
            }
            if args.exclude_homology_hits:
                selected = [
                    item for item in selected
                    if sequence_id(item.record.sequence) not in hit_queries
                ]
                selected, campaign_counts = _select_enzengdb_observations(
                    selected,
                    min_campaign_size=args.min_campaign_size,
                    max_records_per_campaign=0,
                )
                homology_overlap["retained_homology_cold_records"] = len(selected)
                homology_overlap["retained_homology_cold_campaigns"] = len(campaign_counts)
                exact_overlap = {
                    "sequence": sum(item.record.sequence in sequences for item in selected),
                    "substrate": sum(item.record.substrate_smiles in substrates for item in selected),
                    "sequence_substrate_pair": sum(
                        (item.record.sequence, item.record.substrate_smiles) in pairs
                        for item in selected),
                    "unikp_publication_overlap": "unavailable_no_record_level_citations",
                }
    selection = {
        "format_version": 1,
        "dataset": "EnzEngDB v1",
        "dataset_doi": ENZENGDB_DOI,
        "license": ENZENGDB_LICENSE,
        "source_archive_identity": _manifest_dict(archive_identity),
        "selection_parameters": {
            "min_campaign_size": args.min_campaign_size,
            "max_records_per_campaign": args.max_records_per_campaign,
            "exclude_homology_hits": args.exclude_homology_hits,
            "homology_threshold": {
                "min_identity": args.min_identity,
                "coverage": args.coverage,
                "coverage_mode": 0,
            } if args.exclude_homology_hits else None,
            "sampling": "lowest SHA256(candidate_id), independent of fitness_value",
        },
        "campaign_counts": campaign_counts,
        "candidate_ids": [item.record.candidate_id for item in selected],
    }
    if args.selection_output:
        Path(args.selection_output).write_text(
            json.dumps(selection, indent=2, sort_keys=True), encoding="utf-8")
    output = {
        "dataset": selection["dataset"],
        "dataset_doi": ENZENGDB_DOI,
        "license": ENZENGDB_LICENSE,
        "source_archive": str(Path(args.archive).resolve()),
        "source_archive_identity": selection["source_archive_identity"],
        "adapter": {
            "source_files": result.audit.source_files,
            "total_rows": result.audit.total_rows,
            "accepted_rows": result.audit.accepted_rows,
            "rejected_rows": len(result.audit.rejected),
            "accepted_campaigns": result.audit.campaigns,
            "rejection_reasons": dict(sorted(rejection_reasons.items())),
        },
        "frozen_selection": {
            "records": len(selected),
            "campaigns": len(campaign_counts),
            "campaign_counts": campaign_counts,
            "manifest": str(Path(args.selection_output).resolve()) if args.selection_output else None,
        },
        "exact_unikp_overlap": exact_overlap,
        "homology_unikp_overlap": homology_overlap,
        "limitations": [
            "fitness_value has campaign-specific units and is valid only for within-campaign ranking.",
            "Largest concrete carbon-containing reactant is selected deterministically.",
            "Publication overlap cannot be audited because UniKP lacks row-level citations.",
            ("Homology overlap was not requested." if homology_overlap is None else
             "MMseqs2 hits were removed from the frozen selection."
             if args.exclude_homology_hits else
             "Homology hits are reported, not silently removed from this transfer benchmark."),
        ],
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def evaluate_enzengdb_command(args: argparse.Namespace) -> None:
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    if args.task not in model.config.task_names:
        raise ValueError(f"checkpoint does not define task: {args.task}")
    validation_tasks = checkpoint.get("metrics", {}).get("validation_tasks", {})
    if int(validation_tasks.get(args.task, {}).get("count", 0)) == 0:
        raise ValueError(f"checkpoint task has no observed validation labels: {args.task}")
    selection = json.loads(Path(args.selection_manifest).read_text(encoding="utf-8"))
    archive_identity = build_manifest(args.archive)
    if selection.get("source_archive_identity") != _manifest_dict(archive_identity):
        raise ValueError("selection manifest does not match the supplied EnzEngDB archive")
    result = read_enzengdb_directory(args.experiments)
    by_id = {item.record.candidate_id: item for item in result.observations}
    requested_ids = selection.get("candidate_ids")
    if not isinstance(requested_ids, list) or not requested_ids:
        raise ValueError("selection manifest has no candidate_ids")
    missing = [candidate_id for candidate_id in requested_ids if candidate_id not in by_id]
    if missing:
        raise ValueError(f"selection contains {len(missing)} candidates absent from source data")
    selected = [by_id[candidate_id] for candidate_id in requested_ids]
    records = [item.record for item in selected]
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
    )
    loader = DataLoader(
        RecordDataset(records), batch_size=args.batch_size, shuffle=False,
        collate_fn=collator, num_workers=0,
    )
    task_index = model.config.task_names.index(args.task)
    predictions = []
    model.eval()
    with torch.no_grad():
        for batch in loader:
            batch = {name: value.to(device) for name, value in batch.items()}
            predictions.extend(model(batch)["mean"][:, task_index].cpu().tolist())
    campaigns = defaultdict(lambda: {"predictions": [], "observations": []})
    for observation, prediction in zip(selected, predictions):
        campaigns[observation.campaign_id]["predictions"].append(prediction)
        campaigns[observation.campaign_id]["observations"].append(observation.fitness_value)
    campaign_metrics = {}
    for campaign_id, values in sorted(campaigns.items()):
        campaign_seed = int(hashlib.sha256(campaign_id.encode("utf-8")).hexdigest()[:8], 16)
        campaign_metrics[campaign_id] = campaign_ranking_metrics(
            values["predictions"], values["observations"], seed=campaign_seed)
    metric_names = ("spearman", "pairwise_accuracy", "top_10pct_recall", "top_10pct_enrichment")
    macro = {}
    for name in metric_names:
        values = [metrics[name] for metrics in campaign_metrics.values()
                  if metrics.get(name) is not None]
        macro[name] = sum(values) / len(values) if values else None
    if args.task == "activity_rank":
        claim_boundary = (
            "Campaign-rank model evaluated on the frozen EnzEngDB selection. Metrics are "
            "within-campaign ranking measures, not absolute activity or kcat validation."
        )
    else:
        claim_boundary = (
            "Zero-shot cross-endpoint campaign ranking. The checkpoint was trained on kcat, not "
            "EnzEngDB fitness; this is not external kcat validation or activity-head validation."
        )
    output = {
        "dataset": "EnzEngDB v1",
        "dataset_doi": ENZENGDB_DOI,
        "license": ENZENGDB_LICENSE,
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "selection_identity": _manifest_dict(build_manifest(args.selection_manifest)),
        "task_used_as_proxy": args.task,
        "records": len(selected),
        "campaigns": len(campaign_metrics),
        "macro_metrics": macro,
        "campaign_metrics": campaign_metrics,
        "random_expectations": {
            "spearman": 0.0,
            "pairwise_accuracy": 0.5,
            "top_10pct_enrichment": 1.0,
        },
        "claim_boundary": claim_boundary,
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    print(text)


def evaluate_lunzer_command(args: argparse.Namespace) -> None:
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    if args.task not in model.config.task_names:
        raise ValueError(f"checkpoint does not define task: {args.task}")
    validation_tasks = checkpoint.get("metrics", {}).get("validation_tasks", {})
    if int(validation_tasks.get(args.task, {}).get("count", 0)) == 0:
        raise ValueError(f"checkpoint task has no observed validation labels: {args.task}")
    selection = json.loads(Path(args.selection_manifest).read_text(encoding="utf-8"))
    source_identity = build_manifest(args.data)
    if selection.get("source_identity") != _manifest_dict(source_identity):
        raise ValueError("selection manifest does not match the supplied Lunzer source")
    observations = read_lunzer_tsv(args.data)
    by_id = {item.record.candidate_id: item for item in observations}
    requested_ids = selection.get("candidate_ids")
    if (not isinstance(requested_ids, list) or len(requested_ids) != len(by_id)
            or set(requested_ids) != set(by_id)):
        raise ValueError("selection manifest does not contain the complete frozen landscape")
    selected = [by_id[candidate_id] for candidate_id in requested_ids]
    records = [item.record for item in selected]
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
    )
    loader = DataLoader(
        RecordDataset(records), batch_size=args.batch_size, shuffle=False,
        collate_fn=collator, num_workers=0,
    )
    task_index = model.config.task_names.index(args.task)
    predictions = []
    model.eval()
    with torch.no_grad():
        for batch in loader:
            batch = {name: value.to(device) for name, value in batch.items()}
            predictions.extend(model(batch)["mean"][:, task_index].cpu().tolist())
    cofactors = defaultdict(lambda: {"predictions": [], "observations": []})
    prediction_rows = []
    for observation, prediction in zip(selected, predictions):
        cofactors[observation.cofactor]["predictions"].append(prediction)
        cofactors[observation.cofactor]["observations"].append(observation.derived_ln_kcat)
        prediction_rows.append({
            "candidate_id": observation.record.candidate_id,
            "genotype": observation.genotype,
            "cofactor": observation.cofactor,
            "predicted_log10_kcat": prediction,
            "observed_derived_ln_kcat": observation.derived_ln_kcat,
        })
    cofactor_metrics = {
        cofactor: campaign_ranking_metrics(
            values["predictions"], values["observations"],
            seed=int(hashlib.sha256(cofactor.encode("ascii")).hexdigest()[:8], 16),
        )
        for cofactor, values in sorted(cofactors.items())
    }
    metric_names = (
        "spearman", "pairwise_accuracy", "top_10pct_recall", "top_10pct_enrichment")
    macro = {}
    for name in metric_names:
        values = [metrics[name] for metrics in cofactor_metrics.values()
                  if metrics.get(name) is not None]
        macro[name] = sum(values) / len(values) if values else None
    output = {
        "dataset": "Lunzer ancient adaptive landscape",
        "dataset_doi": LUNZER_DATASET_DOI,
        "article_doi": LUNZER_ARTICLE_DOI,
        "license": LUNZER_DATASET_LICENSE,
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "selection_identity": _manifest_dict(build_manifest(args.selection_manifest)),
        "task_used_as_proxy": args.task,
        "records": len(selected),
        "genotypes": len({item.genotype for item in selected}),
        "macro_metrics": macro,
        "cofactor_metrics": cofactor_metrics,
        "predictions": prediction_rows,
        "random_expectations": {
            "spearman": 0.0,
            "pairwise_accuracy": 0.5,
            "top_10pct_enrichment": 1.0,
        },
        "claim_boundary": (
            "Frozen within-cofactor mutant ranking on one E. coli IMDH scaffold. This is direct "
            "kinetic evidence but is neither homology-cold nor an absolute-RMSE validation."
        ),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def evaluate_lunzer_dlkcat_output_command(args: argparse.Namespace) -> None:
    selection = json.loads(Path(args.selection_manifest).read_text(encoding="utf-8"))
    source_identity = build_manifest(args.data)
    if selection.get("source_identity") != _manifest_dict(source_identity):
        raise ValueError("selection manifest does not match the supplied Lunzer source")
    observations = read_lunzer_tsv(args.data)
    by_id = {item.record.candidate_id: item for item in observations}
    requested_ids = selection.get("candidate_ids")
    if (not isinstance(requested_ids, list) or len(requested_ids) != len(by_id)
            or set(requested_ids) != set(by_id)):
        raise ValueError("selection manifest does not contain the complete frozen landscape")
    prediction_by_id = {}
    with Path(args.predictions).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Substrate Name", "Substrate SMILES", "Protein Sequence", "Kcat value (1/s)"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("DLKcat output is missing required columns")
        for row_number, row in enumerate(reader, start=2):
            candidate_id = row["Substrate Name"].strip()
            if candidate_id in prediction_by_id:
                raise ValueError(f"duplicate DLKcat candidate at row {row_number}: {candidate_id}")
            if candidate_id not in by_id:
                raise ValueError(f"unknown DLKcat candidate at row {row_number}: {candidate_id}")
            if row["Protein Sequence"].strip() != by_id[candidate_id].record.sequence:
                raise ValueError(f"DLKcat sequence mismatch at row {row_number}: {candidate_id}")
            try:
                value = float(row["Kcat value (1/s)"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid DLKcat prediction at row {row_number}") from exc
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"non-positive DLKcat prediction at row {row_number}")
            prediction_by_id[candidate_id] = math.log10(value)
    if set(prediction_by_id) != set(requested_ids):
        missing = set(requested_ids) - set(prediction_by_id)
        raise ValueError(f"DLKcat output is missing {len(missing)} frozen candidates")

    cofactors = defaultdict(lambda: {"predictions": [], "observations": []})
    prediction_rows = []
    for candidate_id in requested_ids:
        observation = by_id[candidate_id]
        prediction = prediction_by_id[candidate_id]
        cofactors[observation.cofactor]["predictions"].append(prediction)
        cofactors[observation.cofactor]["observations"].append(observation.derived_ln_kcat)
        prediction_rows.append({
            "candidate_id": candidate_id,
            "genotype": observation.genotype,
            "cofactor": observation.cofactor,
            "predicted_log10_kcat": prediction,
            "observed_derived_ln_kcat": observation.derived_ln_kcat,
        })
    cofactor_metrics = {
        cofactor: campaign_ranking_metrics(
            values["predictions"], values["observations"],
            seed=int(hashlib.sha256(cofactor.encode("ascii")).hexdigest()[:8], 16),
        )
        for cofactor, values in sorted(cofactors.items())
    }
    metric_names = (
        "spearman", "pairwise_accuracy", "top_10pct_recall", "top_10pct_enrichment")
    macro = {}
    for name in metric_names:
        values = [metrics[name] for metrics in cofactor_metrics.values()
                  if metrics.get(name) is not None]
        macro[name] = sum(values) / len(values) if values else None
    checkpoint_path = Path(args.dlkcat_checkpoint)
    checkpoint_bytes = checkpoint_path.read_bytes()
    output = {
        "dataset": "Lunzer ancient adaptive landscape",
        "dataset_doi": LUNZER_DATASET_DOI,
        "predictor": "DLKcat",
        "predictor_version": "upstream-7c15d0d4a7ac",
        "predictor_license": "GPL-3.0-only",
        "checkpoint_identity": {
            "sha256": hashlib.sha256(checkpoint_bytes).hexdigest(),
            "size_bytes": len(checkpoint_bytes),
        },
        "selection_identity": _manifest_dict(build_manifest(args.selection_manifest)),
        "prediction_output_identity": _manifest_dict(build_manifest(args.predictions)),
        "records": len(requested_ids),
        "genotypes": len({item.genotype for item in observations}),
        "macro_metrics": macro,
        "cofactor_metrics": cofactor_metrics,
        "predictions": prediction_rows,
        "random_expectations": {
            "spearman": 0.0,
            "pairwise_accuracy": 0.5,
            "top_10pct_enrichment": 1.0,
        },
        "claim_boundary": (
            "Established-predictor comparison on the frozen within-cofactor mutant landscape. "
            "DLKcat training-source overlap is not auditable at record level; this is neither "
            "homology-cold nor an absolute-RMSE validation."
        ),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def build_enzengdb_rank_split_command(args: argparse.Namespace) -> None:
    archive_identity = build_manifest(args.archive)
    if (archive_identity.sha256 != ENZENGDB_ARCHIVE_SHA256
            or archive_identity.size_bytes != ENZENGDB_ARCHIVE_SIZE_BYTES):
        raise ValueError("EnzEngDB archive does not match the frozen Zenodo v1 identity")
    frozen_test = json.loads(Path(args.test_selection).read_text(encoding="utf-8"))
    if frozen_test.get("source_archive_identity") != _manifest_dict(archive_identity):
        raise ValueError("frozen test selection does not match the EnzEngDB archive")
    test_campaigns = set(frozen_test.get("campaign_counts", {}))
    if not test_campaigns:
        raise ValueError("frozen test selection has no campaign_counts")

    result = read_enzengdb_directory(args.experiments)
    representatives = campaign_representative_records(result.observations)
    output_dir = Path(args.output_dir)
    representative_fasta = output_dir / "campaign_representatives.fasta"
    write_unique_fasta(representatives, representative_fasta)
    default_cluster_path = output_dir / "mmseqs" / "proteins_cluster.tsv"
    if args.cluster_tsv:
        cluster_path = Path(args.cluster_tsv)
        if not cluster_path.is_file():
            raise FileNotFoundError(f"campaign cluster TSV not found: {cluster_path}")
    elif default_cluster_path.is_file():
        cluster_path = default_cluster_path
    else:
        cluster_path = run_mmseqs_easy_cluster(
            args.mmseqs, representative_fasta, output_dir / "mmseqs",
            min_identity=args.min_identity, coverage=args.coverage, threads=args.threads)
    clusters = read_mmseqs_clusters(cluster_path)
    campaign_families = {
        record.campaign_group: clusters[sequence_id(record.sequence)]
        for record in representatives
    }
    test_families = {campaign_families[campaign] for campaign in test_campaigns}
    eligible_families = set(campaign_families.values()) - test_families
    if len(eligible_families) < 2:
        raise ValueError("fewer than two campaign families remain after frozen-test exclusion")
    eligible = [
        item for item in result.observations
        if campaign_families[item.campaign_id] in eligible_families
    ]
    selected, _ = _select_enzengdb_observations(
        eligible,
        min_campaign_size=args.min_campaign_size,
        max_records_per_campaign=args.max_records_per_campaign,
    )
    family_weights = Counter(campaign_families[item.campaign_id] for item in selected)
    if len(family_weights) < 2:
        raise ValueError("fewer than two eligible campaign families pass the row-count gate")
    target_validation_rows = round(len(selected) * args.validation_fraction)
    family_order = sorted(
        family_weights,
        key=lambda family: (
            -family_weights[family],
            hashlib.sha256(f"{args.seed}:{family}".encode("utf-8")).hexdigest(),
        ),
    )
    validation_families = set()
    validation_rows = 0
    for family in family_order:
        proposed = validation_rows + family_weights[family]
        if abs(proposed - target_validation_rows) < abs(validation_rows - target_validation_rows):
            validation_families.add(family)
            validation_rows = proposed
    if not validation_families:
        family = min(family_weights, key=lambda item: abs(family_weights[item] - target_validation_rows))
        validation_families.add(family)
    if validation_families == set(family_weights):
        validation_families.remove(min(validation_families, key=family_weights.__getitem__))

    selected_campaigns = {item.campaign_id for item in selected}
    campaign_partitions = {}
    for campaign, family in campaign_families.items():
        if family in test_families:
            campaign_partitions[campaign] = "frozen_test_family_excluded"
        elif campaign not in selected_campaigns:
            campaign_partitions[campaign] = "row_count_gate_excluded"
        elif family in validation_families:
            campaign_partitions[campaign] = "validation"
        else:
            campaign_partitions[campaign] = "train"
    partition_ids = {
        partition: [
            item.record.candidate_id for item in selected
            if campaign_partitions[item.campaign_id] == partition
        ]
        for partition in ("train", "validation")
    }
    if not partition_ids["train"] or not partition_ids["validation"]:
        raise ValueError("campaign-family split produced an empty train or validation partition")
    payload = {
        "format_version": 1,
        "dataset": "EnzEngDB v1",
        "dataset_doi": ENZENGDB_DOI,
        "source_archive_identity": _manifest_dict(archive_identity),
        "frozen_test_selection_identity": _manifest_dict(build_manifest(args.test_selection)),
        "parameters": {
            "method": "MMseqs2 campaign-representative family split",
            "cluster_identity": _manifest_dict(build_manifest(cluster_path)),
            "min_identity": args.min_identity,
            "coverage": args.coverage,
            "coverage_mode": 0,
            "seed": args.seed,
            "validation_fraction": args.validation_fraction,
            "validation_assignment": "family-weighted greedy by capped record count",
            "min_campaign_size": args.min_campaign_size,
            "max_records_per_campaign": args.max_records_per_campaign,
            "selection": "lowest SHA256(candidate_id), independent of fitness_value",
        },
        "campaign_families": campaign_families,
        "campaign_partitions": campaign_partitions,
        "family_counts": {
            "total": len(set(campaign_families.values())),
            "frozen_test_excluded": len(test_families),
            "row_count_eligible": len(family_weights),
            "train": len(set(family_weights) - validation_families),
            "validation": len(validation_families),
        },
        "partition_candidate_ids": partition_ids,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / "rank_split.json"
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "path": str(target.resolve()),
        "records": {name: len(rows) for name, rows in partition_ids.items()},
        "campaigns": dict(Counter(campaign_partitions.values())),
        "families": payload["family_counts"],
    }, indent=2, sort_keys=True))


def _initialize_activity_ranker(base_checkpoint: str, device: torch.device,
                                *, use_protein: bool, use_molecule: bool,
                                use_context: bool) -> tuple[BioCandidateRanker, dict]:
    base_model, base_payload = load_checkpoint(base_checkpoint, device)
    values = base_model.config.to_dict()
    values.update({
        "task_names": ("activity_rank",),
        "use_protein": use_protein,
        "use_molecule": use_molecule,
        "use_context": use_context,
        "shared_task_query": False,
    })
    model = BioCandidateRanker(ModelConfig.from_dict(values)).to(device)
    source_state = base_model.state_dict()
    target_state = model.state_dict()
    transferred = []
    for name, value in source_state.items():
        if name.startswith("heads.") or name == "task_queries":
            continue
        if name in target_state and target_state[name].shape == value.shape:
            target_state[name] = value.detach().clone()
            transferred.append(name)
    base_task_index = base_model.config.task_names.index("log10_kcat")
    target_state["task_queries"] = source_state["task_queries"][
        base_task_index:base_task_index + 1].detach().clone()
    model.load_state_dict(target_state)
    return model, {
        "checkpoint": str(Path(base_checkpoint).resolve()),
        "checkpoint_sha256": hashlib.sha256(Path(base_checkpoint).read_bytes()).hexdigest(),
        "source_task_query": "log10_kcat",
        "transferred_tensor_count": len(transferred) + 1,
        "ranking_head": "randomly_initialized",
        "optimizer_state": "not_transferred",
    }


@torch.no_grad()
def _evaluate_campaign_ranker(model, records, device, batch_size: int) -> dict:
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
        task_names=model.config.task_names,
    )
    loader = DataLoader(
        RecordDataset(records), batch_size=batch_size, shuffle=False,
        collate_fn=collator, num_workers=0)
    model.eval()
    predictions = []
    for batch in loader:
        batch = {name: value.to(device) for name, value in batch.items()}
        predictions.extend(model(batch)["mean"][:, 0].cpu().tolist())
    campaigns = defaultdict(lambda: {"predictions": [], "labels": []})
    for record, prediction in zip(records, predictions):
        campaigns[record.campaign_group]["predictions"].append(prediction)
        campaigns[record.campaign_group]["labels"].append(record.activity_rank)
    metrics = {}
    for campaign, values in sorted(campaigns.items()):
        seed = int(hashlib.sha256(campaign.encode("utf-8")).hexdigest()[:8], 16)
        metrics[campaign] = campaign_ranking_metrics(
            values["predictions"], values["labels"], seed=seed)
    names = ("spearman", "pairwise_accuracy", "top_10pct_recall", "top_10pct_enrichment")
    macro = {
        name: sum(item[name] for item in metrics.values() if item[name] is not None)
        / sum(item[name] is not None for item in metrics.values())
        for name in names
    }
    return {"count": len(records), "campaign_count": len(metrics), "macro": macro}


def train_enzengdb_ranker_command(args: argparse.Namespace) -> None:
    torch.manual_seed(args.seed)
    split = json.loads(Path(args.split_manifest).read_text(encoding="utf-8"))
    archive_identity = build_manifest(args.archive)
    if split.get("source_archive_identity") != _manifest_dict(archive_identity):
        raise ValueError("rank split does not match the supplied EnzEngDB archive")
    result = read_enzengdb_directory(args.experiments)
    by_id = {item.record.candidate_id: item for item in result.observations}
    partitions = split.get("partition_candidate_ids", {})
    observations = {}
    for partition in ("train", "validation"):
        identifiers = partitions.get(partition, [])
        missing = [identifier for identifier in identifiers if identifier not in by_id]
        if missing:
            raise ValueError(f"rank split {partition} has {len(missing)} missing candidates")
        observations[partition] = [by_id[identifier] for identifier in identifiers]
    train_records = campaign_rank_records(observations["train"])
    validation_records = campaign_rank_records(observations["validation"])
    device = torch.device(args.device)
    model, initialization = _initialize_activity_ranker(
        args.initialize_from, device,
        use_protein=not args.disable_protein,
        use_molecule=not args.disable_molecule,
        use_context=not args.disable_context,
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
        task_names=model.config.task_names,
    )
    train_sampler = CampaignBatchSampler(
        train_records, args.batch_size, shuffle=True, seed=args.seed)
    validation_sampler = CampaignGroupedBatchSampler(validation_records)
    train_loader = DataLoader(
        RecordDataset(train_records), batch_sampler=train_sampler,
        collate_fn=collator, num_workers=0)
    validation_loader = DataLoader(
        RecordDataset(validation_records), batch_sampler=validation_sampler,
        collate_fn=collator, num_workers=0)
    output_dir = Path(args.output_dir)
    data_manifest = {
        "source_archive": _manifest_dict(archive_identity),
        "split": {
            "path": str(Path(args.split_manifest).resolve()),
            "identity": _manifest_dict(build_manifest(args.split_manifest)),
        },
        "initialization": initialization,
    }
    best_loss = float("inf")
    for epoch in range(1, args.epochs + 1):
        train_sampler.set_epoch(epoch)
        train_metrics = run_ranking_epoch(model, train_loader, device, optimizer)
        validation_metrics = run_ranking_epoch(model, validation_loader, device)
        ranking = _evaluate_campaign_ranker(
            model, validation_records, device, args.batch_size)
        metrics = {
            "epoch": epoch,
            "train_loss": train_metrics.loss,
            "train_pairs": train_metrics.pair_count,
            "validation_loss": validation_metrics.loss,
            "validation_pairs": validation_metrics.pair_count,
            "validation_tasks": {"activity_rank": ranking},
        }
        print(json.dumps(metrics, sort_keys=True))
        save_checkpoint(
            output_dir / "latest.pt", model, optimizer, epoch=epoch,
            data_manifest=data_manifest, metrics=metrics)
        if validation_metrics.loss < best_loss:
            best_loss = validation_metrics.loss
            save_checkpoint(
                output_dir / "best.pt", model, optimizer, epoch=epoch,
                data_manifest=data_manifest, metrics=metrics)


def summarize_enzengdb_runs_command(args: argparse.Namespace) -> None:
    runs = [json.loads(Path(path).read_text(encoding="utf-8")) for path in args.inputs]
    summary = summarize_campaign_runs(
        runs, bootstrap_samples=args.bootstrap_samples, seed=args.seed)
    summary.update({
        "inputs": [str(Path(path).resolve()) for path in args.inputs],
        "claim_boundary": (
            "Seed variation and campaign-bootstrap uncertainty for within-campaign ranking. "
            "Intervals do not convert heterogeneous fitness endpoints into absolute activity."
        ),
    })
    text = json.dumps(summary, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def train_command(args: argparse.Namespace) -> None:
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    manifest_payload = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    identity = manifest_payload.get("identity", manifest_payload)
    manifest = FileManifest(**identity)
    verify_manifest(args.data, manifest)
    if args.data_format == "unikp":
        read_result = read_unikp_json(args.data)
        if read_result.audit.conflict_count and args.conflict_policy == "error":
            raise ValueError(
                f"found {read_result.audit.conflict_count} conflicting enzyme-substrate labels; "
                "curate them or pass --conflict-policy keep for development-only runs")
        source_records = read_result.records
        if args.conflict_policy == "median":
            if args.split_manifest:
                source_records = apply_split_manifest(
                    source_records, args.split_manifest, source_identity=identity)
            source_records = aggregate_pair_measurements(source_records)
    else:
        source_records = read_normalized_kinetics_jsonl(args.data).records
    task_names = tuple(args.tasks or (
        ("log10_kcat",) if args.data_format == "unikp" else KINETIC_TASKS
    ))
    if args.data_format == "unikp" and task_names != ("log10_kcat",):
        raise ValueError("UniKP supports only the log10_kcat task")
    if args.split_manifest:
        records = (
            source_records
            if args.data_format == "unikp" and args.conflict_policy == "median"
            else apply_split_manifest(source_records, args.split_manifest, source_identity=identity)
        )
    else:
        records = split_records(source_records, args.split_strategy, seed=args.seed)
    if args.max_records:
        records = records[:args.max_records]
    train_records = [record for record in records if record.split == "train"]
    validation_records = [record for record in records if record.split == "validation"]
    if not train_records or not validation_records:
        raise ValueError("split produced an empty train or validation partition")
    for task in task_names:
        for partition, partition_records in (
            ("train", train_records), ("validation", validation_records),
        ):
            if not any(getattr(record, task) is not None for record in partition_records):
                raise ValueError(f"task {task} has zero {partition} observations")

    config = ModelConfig(
        d_model=args.d_model,
        num_heads=args.num_heads,
        protein_layers=args.protein_layers,
        molecule_layers=args.molecule_layers,
        fusion_layers=args.fusion_layers,
        protein_chunk_size=args.chunk_size,
        context_buckets=args.context_buckets,
        use_protein=not args.disable_protein,
        use_molecule=not args.disable_molecule,
        use_context=not args.disable_context,
        shared_task_query=args.shared_task_query,
        fusion_mode=args.fusion_mode,
        protein_encoder=args.protein_encoder,
        uncertainty_mode=args.uncertainty_mode,
        esm2_model_name=args.esm2_model,
        esm2_frozen=not args.esm2_unfreeze,
        task_names=task_names,
    )
    device = torch.device(args.device)
    start_epoch = 1
    if args.resume:
        model, resume_payload = load_checkpoint(args.resume, device)
        if model.config.to_dict() != config.to_dict():
            raise ValueError("resume checkpoint model config does not match CLI config")
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
        optimizer.load_state_dict(resume_payload["optimizer_state_dict"])
        start_epoch = int(resume_payload["epoch"]) + 1
        torch.set_rng_state(resume_payload["torch_rng_state"].cpu())
        if torch.cuda.is_available() and resume_payload.get("cuda_rng_state_all"):
            torch.cuda.set_rng_state_all(
                [state.cpu() for state in resume_payload["cuda_rng_state_all"]])
    else:
        model = BioCandidateRanker(config).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    collator = EnzymeSubstrateCollator(
        context_buckets=config.context_buckets,
        task_names=config.task_names,
    )
    generator = torch.Generator().manual_seed(args.seed)
    if args.resume and resume_payload.get("data_loader_generator_state") is not None:
        generator.set_state(resume_payload["data_loader_generator_state"].cpu())
    train_loader = DataLoader(
        RecordDataset(train_records), batch_size=args.batch_size, shuffle=True,
        collate_fn=collator, num_workers=0, generator=generator)
    validation_loader = DataLoader(
        RecordDataset(validation_records), batch_size=args.batch_size, shuffle=False,
        collate_fn=collator, num_workers=0)

    best_validation = float("inf")
    output_dir = Path(args.output_dir)
    scheduler = None
    if args.cosine_schedule:
        scheduler = create_lr_scheduler(
            optimizer, warmup_epochs=args.warmup_epochs, total_epochs=args.epochs)
    run_manifest = {"source": manifest_payload}
    if args.split_manifest:
        run_manifest["split"] = {
            "path": str(Path(args.split_manifest).resolve()),
            "identity": _manifest_dict(build_manifest(args.split_manifest)),
        }
    else:
        run_manifest["split"] = {
            "strategy": args.split_strategy,
            "seed": args.seed,
            "warning": "development split; not homology-cold",
        }
    if args.resume and resume_payload.get("data_manifest") != run_manifest:
        raise ValueError("resume checkpoint data identity or split does not match this run")
    baseline = mean_baseline_metrics(train_records, validation_records, config.task_names)
    print(json.dumps({"mean_baseline": baseline}, sort_keys=True))
    epochs_without_improvement = 0
    for epoch in range(start_epoch, args.epochs + 1):
        if args.grad_accum_steps > 1:
            train_metrics = run_epoch_accum(
                model, train_loader, device, optimizer,
                accum_steps=args.grad_accum_steps)
        else:
            train_metrics = run_epoch(model, train_loader, device, optimizer)
        validation_metrics = run_epoch(model, validation_loader, device)
        metrics = {
            "epoch": epoch,
            "train_loss": train_metrics.loss,
            "validation_loss": validation_metrics.loss,
            "train_observed_labels": train_metrics.observed_labels,
            "validation_observed_labels": validation_metrics.observed_labels,
            "validation_tasks": evaluate_model(model, validation_loader, device, config),
            "mean_baseline": baseline,
        }
        print(json.dumps(metrics, sort_keys=True))
        save_checkpoint(
            output_dir / "latest.pt", model, optimizer, epoch=epoch,
            data_manifest=run_manifest, metrics=metrics,
            data_loader_generator_state=generator.get_state())
        if validation_metrics.loss < best_validation:
            best_validation = validation_metrics.loss
            epochs_without_improvement = 0
            save_checkpoint(
                output_dir / "best.pt", model, optimizer, epoch=epoch,
                data_manifest=run_manifest, metrics=metrics,
                data_loader_generator_state=generator.get_state())
        else:
            epochs_without_improvement += 1
        if scheduler is not None:
            scheduler.step()
        if args.patience > 0 and epochs_without_improvement >= args.patience:
            print(json.dumps({"early_stop": epoch, "patience": args.patience}))
            break


def homology_command(args: argparse.Namespace) -> None:
    manifest_payload = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    verify_manifest(args.data, FileManifest(**manifest_payload["identity"]))
    records = read_unikp_json(args.data).records
    output_dir = Path(args.output_dir)
    fasta_path = output_dir / "unique_proteins.fasta"
    write_unique_fasta(records, fasta_path)
    if args.cluster_tsv:
        cluster_path = Path(args.cluster_tsv)
        if not cluster_path.is_file():
            raise FileNotFoundError(f"cluster TSV not found: {cluster_path}")
    else:
        if not args.mmseqs:
            raise ValueError("--mmseqs is required unless --cluster-tsv is provided")
        cluster_path = run_mmseqs_easy_cluster(
            args.mmseqs, fasta_path, output_dir / "mmseqs",
            min_identity=args.min_identity, coverage=args.coverage, threads=args.threads)
    clusters = read_mmseqs_clusters(cluster_path)
    split_records_result = assign_homology_splits(records, clusters, seed=args.seed)
    split_path = output_dir / "homology_split.json"
    write_split_manifest(
        split_records_result, clusters, split_path,
        parameters={
            "method": "MMseqs2 linclust",
            "mmseqs_backend": args.mmseqs,
            "cluster_tsv": str(cluster_path.resolve()),
            "min_identity": args.min_identity,
            "coverage": args.coverage,
            "coverage_mode": 0,
            "seed": args.seed,
            "source_identity": manifest_payload["identity"],
        },
    )
    print(str(split_path.resolve()))


def cold_split_command(args: argparse.Namespace) -> None:
    manifest_payload = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    verify_manifest(args.data, FileManifest(**manifest_payload["identity"]))
    records = read_unikp_json(args.data).records
    assigned = split_records(records, args.strategy, seed=args.seed)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = write_cold_split_manifest(
        assigned, args.strategy, str(target), seed=args.seed,
        source_identity=manifest_payload["identity"])
    print(json.dumps({
        "path": str(target.resolve()),
        "split_counts": payload["split_counts"],
        "audit": payload["audit"],
    }, sort_keys=True))


def predict_command(args: argparse.Namespace) -> None:
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    calibrator = None
    calibration_task = None
    calibration_payload = None
    conformal_intervals = {}
    if args.calibration_artifact:
        calibrator, calibration_payload = load_calibration_artifact(args.calibration_artifact)
        identities = calibration_payload["identities"]
        calibration_task = identities.get("task")
        if calibration_task not in model.config.task_names:
            raise ValueError("calibration artifact task is not defined by the checkpoint")
        expected = _calibration_identities(args.checkpoint, calibration_task)
        if identities != expected:
            raise ValueError("calibration artifact identity mismatch")
        conformal_intervals = {
            name: conformal_calibrator_from_dict(payload)
            for name, payload in calibration_payload.get("conformal_intervals", {}).items()
        }
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    rows = payload.get("candidates", payload) if isinstance(payload, dict) else payload
    if not isinstance(rows, list):
        raise ValueError("candidate input must be a JSON list or a candidates list")
    from .data import EnzymeSubstrateRecord, FBAFeatureMetadata
    fba_metadata = None
    fba_features_path = getattr(args, "fba_features", None)
    if fba_features_path:
        fba_payload = json.loads(Path(fba_features_path).read_text(encoding="utf-8"))
        try:
            fba_metadata = FBAFeatureMetadata(
                schema_version=fba_payload["schema_version"],
                feature_ids=tuple(fba_payload["feature_ids"]),
                model_id=fba_payload["model_id"],
                solver_id=fba_payload["solver_id"],
                objective_id=fba_payload["objective_id"],
                condition_id=fba_payload["condition_id"],
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid FBA feature metadata in {fba_features_path}: {exc}") from exc
        if len(fba_metadata.feature_ids) != model.config.fba_context_dim:
            raise ValueError(
                f"FBA feature metadata width {len(fba_metadata.feature_ids)} does not match "
                f"checkpoint configured width {model.config.fba_context_dim}")
    records = []
    for index, row in enumerate(rows):
        fba_context = tuple(row.get("fba_context", ()))
        if fba_context and fba_metadata is None:
            raise ValueError(
                "candidate supplies fba_context but --fba-features metadata was not provided")
        records.append(EnzymeSubstrateRecord(
            candidate_id=str(row.get("candidate_id", index)),
            sequence=row["sequence"],
            substrate_smiles=row["substrate_smiles"],
            substrate_name=str(row.get("substrate_name", "")),
            organism=str(row.get("organism", "")),
            ec=str(row.get("ec", "")),
            enzyme_type=str(row.get("enzyme_type", "unknown")),
            reaction=str(row.get("reaction", "")),
            fba_context=fba_context,
            fba_feature_metadata=fba_metadata,
            source_dataset=str(Path(args.input).resolve()),
            source_row=index,
        ))
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim)
    batch = {name: value.to(device) for name, value in collator(records).items()}
    model.eval()
    with torch.no_grad():
        output = model(batch)
    predictions = []
    validation_tasks = checkpoint.get("metrics", {}).get("validation_tasks", {})
    for row, record in enumerate(records):
        tasks = {}
        for index, task in enumerate(model.config.task_names):
            observed = int(validation_tasks.get(task, {}).get("count", 0))
            if observed == 0:
                tasks[task] = {
                    "status": "untrained_no_labels",
                    "mean": None,
                    "standard_deviation": None,
                }
            else:
                raw_deviation = output["standard_deviation"][row, index]
                calibrated_deviation = None
                if task == calibration_task:
                    calibrated_mean, calibrated = calibrator.calibrate(
                        output["mean"][row, index].reshape(1), raw_deviation.reshape(1))
                    calibrated_deviation = float(calibrated[0])
                    intervals = {}
                    for name, conformal in conformal_intervals.items():
                        bounds = conformal.interval(
                            calibrated_mean,
                            calibrated if conformal.normalized else None,
                        )
                        intervals[name] = {
                            "lower": float(bounds[0][0]),
                            "upper": float(bounds[1][0]),
                        }
                else:
                    intervals = None
                tasks[task] = {
                    "status": "trained",
                    "mean": float(output["mean"][row, index]),
                    "standard_deviation": float(raw_deviation),
                    "raw_standard_deviation": float(raw_deviation),
                    "calibrated_standard_deviation": calibrated_deviation,
                    "conformal_intervals": intervals,
                }
        predictions.append({
            "candidate_id": record.candidate_id,
            "tasks": tasks,
        })
    result = {
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "training_data_manifest": checkpoint["data_manifest"],
        "input_identity": _manifest_dict(build_manifest(args.input)),
        "calibration": None if calibration_payload is None else {
            "artifact": str(Path(args.calibration_artifact).resolve()),
            "artifact_sha256": calibration_payload["artifact_sha256"],
            "method": calibration_payload["calibrator"]["kind"],
            "task": calibration_task,
            "conformal_intervals": sorted(conformal_intervals),
        },
        "predictions": predictions,
        "warning": "Computational ranking only; not experimental validation.",
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def evaluate_command(args: argparse.Namespace) -> None:
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    source_manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    verify_manifest(args.data, FileManifest(**source_manifest["identity"]))
    if args.data_format == "unikp":
        records = read_unikp_json(args.data).records
        if args.tasks and tuple(args.tasks) != ("log10_kcat",):
            raise ValueError("UniKP supports only the log10_kcat task")
        adapter_tasks = ("log10_kcat",)
    else:
        records = read_normalized_kinetics_jsonl(args.data).records
        adapter_tasks = KINETIC_TASKS
    task_names = tuple(args.tasks or (
        task for task in model.config.task_names if task in adapter_tasks
    ))
    if not task_names:
        raise ValueError("data adapter and checkpoint have no tasks in common")
    unsupported = [task for task in task_names if task not in model.config.task_names]
    if unsupported:
        raise ValueError(
            f"requested tasks are not supported by checkpoint: {', '.join(unsupported)}")
    if args.data_format == "normalized-kinetics":
        unsupported = [task for task in task_names if task not in KINETIC_TASKS]
        if unsupported:
            raise ValueError(
                f"requested tasks are not supported by normalized kinetics: "
                f"{', '.join(unsupported)}")
    records = apply_split_manifest(
        records, args.split_manifest, source_identity=source_manifest["identity"])
    if args.data_format == "unikp" and args.conflict_policy == "median":
        records = aggregate_pair_measurements(records)
    train_records = [record for record in records if record.split == "train"]
    evaluation_records = [record for record in records if record.split == args.partition]
    if not evaluation_records:
        raise ValueError(f"partition is empty: {args.partition}")
    run_manifest = {
        "source": source_manifest,
        "split": {
            "path": str(Path(args.split_manifest).resolve()),
            "identity": _manifest_dict(build_manifest(args.split_manifest)),
        },
    }
    if checkpoint.get("data_manifest") != run_manifest:
        raise ValueError("checkpoint was not trained with the supplied source and split identities")
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
        task_names=model.config.task_names)
    loader = DataLoader(
        RecordDataset(evaluation_records), batch_size=args.batch_size,
        shuffle=False, collate_fn=collator, num_workers=0)
    result = {
        "partition": args.partition,
        "record_count": len(evaluation_records),
        "checkpoint": str(Path(args.checkpoint).resolve()),
        "source_identity": source_manifest["identity"],
        "split_identity": run_manifest["split"]["identity"],
        "model_metrics": {
            task: metrics
            for task, metrics in evaluate_model(model, loader, device, model.config).items()
            if task in task_names
        },
        "mean_baseline": mean_baseline_metrics(
            train_records, evaluation_records, task_names),
        "model_parameters": sum(parameter.numel() for parameter in model.parameters()),
        "claim_boundary": "Internal homology-cold development evaluation; not external validation.",
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def evaluate_classical_baselines_command(args: argparse.Namespace) -> None:
    try:
        import numpy as np
        from sklearn.linear_model import Ridge
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise RuntimeError("classical baselines require the baseline optional dependency") from exc

    source_manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    verify_manifest(args.data, FileManifest(**source_manifest["identity"]))
    records = apply_split_manifest(
        read_unikp_json(args.data).records, args.split_manifest,
        source_identity=source_manifest["identity"])
    records = aggregate_pair_measurements(records)
    partitions = {
        name: [record for record in records if record.split == name]
        for name in ("train", "validation", "test")
    }
    if any(not values for values in partitions.values()):
        raise ValueError("classical baseline requires non-empty train, validation, and test")
    labels = {
        name: np.asarray([record.log10_kcat for record in values], dtype=np.float64)
        for name, values in partitions.items()
    }
    feature_builders = {
        "amino_acid_composition_ridge": amino_acid_composition,
        "morgan_r2_2048_ridge": morgan_fingerprints,
    }
    alpha_grid = (0.1, 1.0, 10.0, 100.0)
    results = {}
    for model_name, build_features in feature_builders.items():
        features = {name: build_features(values) for name, values in partitions.items()}
        validation_runs = []
        for alpha in alpha_grid:
            model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
            model.fit(features["train"], labels["train"])
            predictions = model.predict(features["validation"])
            rmse = float(np.sqrt(np.mean((predictions - labels["validation"]) ** 2)))
            validation_runs.append({"alpha": alpha, "rmse": rmse})
        selected = min(validation_runs, key=lambda item: (item["rmse"], item["alpha"]))
        model = make_pipeline(StandardScaler(), Ridge(alpha=selected["alpha"]))
        model.fit(features["train"], labels["train"])
        test_predictions = torch.from_numpy(model.predict(features["test"]))
        test_labels = torch.from_numpy(labels["test"])
        results[model_name] = {
            "feature_dimension": int(features["train"].shape[1]),
            "alpha_grid": list(alpha_grid),
            "validation_runs": validation_runs,
            "selected_alpha": selected["alpha"],
            "test_metrics": regression_metrics(test_predictions, test_labels),
        }
    output = {
        "source_identity": source_manifest["identity"],
        "split_identity": _manifest_dict(build_manifest(args.split_manifest)),
        "conflict_policy": "median_log10_by_sequence_smiles",
        "partition_counts": {name: len(values) for name, values in partitions.items()},
        "selection_policy": "fixed alpha grid selected by validation RMSE; train-only fit",
        "baselines": results,
        "claim_boundary": "Classical internal homology-cold baselines; not external validation.",
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def summarize_fusion_runs_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    seeds = protocol.get("seeds")
    if not isinstance(seeds, list) or not seeds or not all(isinstance(seed, int) for seed in seeds):
        raise ValueError("fusion protocol must define integer seeds")
    expected_source = protocol.get("source_identity")
    expected_split = protocol.get("split_identity")
    root = Path(args.runs_root)
    mode_directories = {"task_query": "task-query", "late_concat": "late-concat"}
    run_metrics = {mode: {} for mode in mode_directories}
    run_identities = {mode: {} for mode in mode_directories}
    for mode, directory_name in mode_directories.items():
        for seed in seeds:
            run_directory = root / f"{directory_name}-seed{seed}"
            metrics_path = run_directory / "test_metrics.json"
            checkpoint_path = run_directory / "best.pt"
            payload = json.loads(metrics_path.read_text(encoding="utf-8"))
            if payload.get("source_identity") != expected_source:
                raise ValueError(f"{mode} seed {seed} source identity mismatch")
            if payload.get("split_identity") != expected_split:
                raise ValueError(f"{mode} seed {seed} split identity mismatch")
            model, checkpoint = load_checkpoint(checkpoint_path, torch.device("cpu"))
            if model.config.fusion_mode != mode:
                raise ValueError(f"{mode} seed {seed} checkpoint fusion mode mismatch")
            checkpoint_source = checkpoint.get("data_manifest", {}).get("source", {}).get("identity")
            if checkpoint_source != expected_source:
                raise ValueError(f"{mode} seed {seed} checkpoint source mismatch")
            checkpoint_split = checkpoint.get("data_manifest", {}).get("split", {}).get("identity")
            if checkpoint_split != expected_split:
                raise ValueError(f"{mode} seed {seed} checkpoint split mismatch")
            run_metrics[mode][seed] = payload.get("model_metrics", {}).get("log10_kcat", {})
            checkpoint_bytes = checkpoint_path.read_bytes()
            run_identities[mode][str(seed)] = {
                "selected_epoch": checkpoint.get("epoch"),
                "model_parameters": payload.get("model_parameters"),
                "test_metrics_identity": _binary_identity(metrics_path),
                "checkpoint_sha256": hashlib.sha256(checkpoint_bytes).hexdigest(),
                "checkpoint_size_bytes": len(checkpoint_bytes),
            }
    summary = summarize_paired_regression_runs(
        run_metrics["task_query"], run_metrics["late_concat"])
    output = {
        "protocol_identity": _binary_identity(args.protocol),
        "protocol_claim_boundary": protocol.get("claim_boundary"),
        "run_identities": run_identities,
        "summary": summary,
        "interpretation": (
            "Post-hoc seed-stability analysis on an already observed internal test. "
            "Signed paired differences are late_concat minus task_query."
        ),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def train_feature_mlp_command(args: argparse.Namespace) -> None:
    import numpy as np

    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    if args.seed not in protocol.get("seeds", []):
        raise ValueError(f"seed is not registered in the feature MLP protocol: {args.seed}")
    feature_config = protocol.get("features", {}).get(args.feature)
    if not isinstance(feature_config, dict):
        raise ValueError(f"feature is not registered in the feature MLP protocol: {args.feature}")
    source_manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    if source_manifest.get("identity") != protocol.get("source_identity"):
        raise ValueError("feature MLP protocol source identity mismatch")
    verify_manifest(args.data, FileManifest(**source_manifest["identity"]))
    split_identity = _manifest_dict(build_manifest(args.split_manifest))
    if split_identity != protocol.get("split_identity"):
        raise ValueError("feature MLP protocol split identity mismatch")
    records = apply_split_manifest(
        read_unikp_json(args.data).records, args.split_manifest,
        source_identity=source_manifest["identity"])
    records = aggregate_pair_measurements(records)
    partitions = {
        name: [record for record in records if record.split == name]
        for name in ("train", "validation", "test")
    }
    if any(not values for values in partitions.values()):
        raise ValueError("feature MLP requires non-empty train, validation, and test")
    if args.feature == "amino_acid_composition":
        feature_builder = amino_acid_composition
    elif args.feature == "morgan_r2_2048":
        feature_builder = morgan_fingerprints
    else:
        raise ValueError(f"unsupported feature: {args.feature}")
    features = {name: feature_builder(values) for name, values in partitions.items()}
    if features["train"].shape[1] != feature_config.get("dimension"):
        raise ValueError("feature dimension does not match the frozen protocol")
    labels = {
        name: np.asarray([record.log10_kcat for record in values], dtype=np.float32)
        for name, values in partitions.items()
    }
    training = protocol["training"]
    model, metrics, feature_mean, feature_stddev = train_feature_mlp(
        features,
        labels,
        seed=args.seed,
        epochs=training["epochs"],
        batch_size=training["batch_size"],
        learning_rate=training["learning_rate"],
        weight_decay=training["weight_decay"],
        hidden_dimensions=tuple(training["hidden_dimensions"]),
        dropout=training["dropout"],
        device=torch.device(args.device),
    )
    output_directory = Path(args.output_dir)
    checkpoint_path = output_directory / "best.pt"
    checkpoint = {
        "format_version": 1,
        "feature": args.feature,
        "seed": args.seed,
        "model_state_dict": model.state_dict(),
        "input_dimension": features["train"].shape[1],
        "hidden_dimensions": tuple(training["hidden_dimensions"]),
        "dropout": training["dropout"],
        "feature_mean": feature_mean,
        "feature_stddev": feature_stddev,
        "protocol_identity": _binary_identity(args.protocol),
        "source_identity": source_manifest["identity"],
        "split_identity": split_identity,
        "selected_epoch": metrics["selected_epoch"],
    }
    save_feature_checkpoint(checkpoint, checkpoint_path)
    checkpoint_identity = _binary_identity(checkpoint_path)
    output = {
        "feature": args.feature,
        "seed": args.seed,
        "protocol_identity": _binary_identity(args.protocol),
        "source_identity": source_manifest["identity"],
        "split_identity": split_identity,
        "partition_counts": {name: len(values) for name, values in partitions.items()},
        "model_parameters": sum(parameter.numel() for parameter in model.parameters()),
        "checkpoint_identity": checkpoint_identity,
        **metrics,
        "claim_boundary": protocol.get("claim_boundary"),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "metrics.json").write_text(text, encoding="utf-8")
    print(text)


def summarize_feature_mlp_runs_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    seeds = protocol.get("seeds")
    if not isinstance(seeds, list) or not seeds:
        raise ValueError("feature MLP protocol must define seeds")
    features = {
        "amino_acid_composition": "amino-acid-composition",
        "morgan_r2_2048": "morgan-r2-2048",
    }
    root = Path(args.runs_root)
    run_metrics = {feature: {} for feature in features}
    run_identities = {feature: {} for feature in features}
    protocol_identity = _binary_identity(args.protocol)
    for feature, directory_name in features.items():
        for seed in seeds:
            run_directory = root / f"{directory_name}-seed{seed}"
            metrics_path = run_directory / "metrics.json"
            checkpoint_path = run_directory / "best.pt"
            payload = json.loads(metrics_path.read_text(encoding="utf-8"))
            if payload.get("feature") != feature or payload.get("seed") != seed:
                raise ValueError(f"feature MLP run identity mismatch: {feature} seed {seed}")
            if payload.get("protocol_identity") != protocol_identity:
                raise ValueError(f"feature MLP protocol mismatch: {feature} seed {seed}")
            if payload.get("source_identity") != protocol.get("source_identity"):
                raise ValueError(f"feature MLP source mismatch: {feature} seed {seed}")
            if payload.get("split_identity") != protocol.get("split_identity"):
                raise ValueError(f"feature MLP split mismatch: {feature} seed {seed}")
            checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
            if (checkpoint.get("feature") != feature or checkpoint.get("seed") != seed
                    or checkpoint.get("protocol_identity") != protocol_identity):
                raise ValueError(f"feature MLP checkpoint mismatch: {feature} seed {seed}")
            run_metrics[feature][seed] = payload["test_metrics"]
            run_identities[feature][str(seed)] = {
                "selected_epoch": payload["selected_epoch"],
                "selected_validation_rmse": payload["selected_validation_rmse"],
                "model_parameters": payload["model_parameters"],
                "metrics_identity": _binary_identity(metrics_path),
                "checkpoint_identity": _binary_identity(checkpoint_path),
            }
    output = {
        "protocol_identity": protocol_identity,
        "protocol_claim_boundary": protocol.get("claim_boundary"),
        "run_identities": run_identities,
        "summary": summarize_regression_seed_runs(run_metrics),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def summarize_internal_ablations_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    expected_source = protocol.get("source_identity")
    expected_split = protocol.get("split_identity")
    baseline_path = Path(args.baseline_metrics)
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if (baseline.get("source_identity") != expected_source
            or baseline.get("split_identity") != expected_split):
        raise ValueError("ablation baseline identity mismatch")
    baseline_metrics = baseline["model_metrics"]["log10_kcat"]
    variants = {
        "no_context": ("no-context", lambda config: not config.use_context),
        "shared_task_query": (
            "shared-task-query", lambda config: config.shared_task_query),
        "global_mean_protein": (
            "global-mean-protein", lambda config: config.protein_encoder == "global_mean"),
    }
    results = {}
    for name, (directory_name, config_check) in variants.items():
        run_directory = Path(args.runs_root) / directory_name
        metrics_path = run_directory / "test_metrics.json"
        checkpoint_path = run_directory / "best.pt"
        payload = json.loads(metrics_path.read_text(encoding="utf-8"))
        if (payload.get("source_identity") != expected_source
                or payload.get("split_identity") != expected_split):
            raise ValueError(f"ablation identity mismatch: {name}")
        model, checkpoint = load_checkpoint(checkpoint_path, torch.device("cpu"))
        if not config_check(model.config):
            raise ValueError(f"ablation checkpoint config mismatch: {name}")
        metrics = payload["model_metrics"]["log10_kcat"]
        results[name] = {
            "metrics": metrics,
            "delta_from_baseline": {
                metric: metrics[metric] - baseline_metrics[metric]
                for metric in ("rmse", "mae", "pearson")
            },
            "selected_epoch": checkpoint["epoch"],
            "model_parameters": payload["model_parameters"],
            "metrics_identity": _binary_identity(metrics_path),
            "checkpoint_identity": _binary_identity(checkpoint_path),
        }
    output = {
        "protocol_identity": _binary_identity(args.protocol),
        "protocol_claim_boundary": protocol.get("claim_boundary"),
        "baseline": {
            "metrics": baseline_metrics,
            "model_parameters": baseline.get("model_parameters"),
            "metrics_identity": _binary_identity(baseline_path),
        },
        "ablations": results,
        "non_identifiable_ablations": protocol.get("non_identifiable_ablations"),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


def summarize_uncertainty_ablation_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    seeds = protocol.get("seeds")
    if not isinstance(seeds, list) or not seeds:
        raise ValueError("uncertainty protocol must define seeds")
    expected_source = protocol.get("source_identity")
    expected_split = protocol.get("split_identity")
    runs = {"heteroscedastic": {}, "fixed_variance_mse": {}}
    identities = {mode: {} for mode in runs}
    roots = {
        "heteroscedastic": Path(args.baseline_runs_root),
        "fixed_variance_mse": Path(args.ablation_runs_root),
    }
    directory_templates = {
        "heteroscedastic": "task-query-seed{seed}",
        "fixed_variance_mse": "fixed-variance-seed{seed}",
    }
    for mode in runs:
        for seed in seeds:
            run_directory = roots[mode] / directory_templates[mode].format(seed=seed)
            metrics_path = run_directory / "test_metrics.json"
            checkpoint_path = run_directory / "best.pt"
            payload = json.loads(metrics_path.read_text(encoding="utf-8"))
            if (payload.get("source_identity") != expected_source
                    or payload.get("split_identity") != expected_split):
                raise ValueError(f"uncertainty run identity mismatch: {mode} seed {seed}")
            model, checkpoint = load_checkpoint(checkpoint_path, torch.device("cpu"))
            if model.config.uncertainty_mode != mode:
                raise ValueError(f"uncertainty checkpoint mode mismatch: {mode} seed {seed}")
            runs[mode][seed] = payload["model_metrics"]["log10_kcat"]
            identities[mode][str(seed)] = {
                "selected_epoch": checkpoint["epoch"],
                "model_parameters": payload["model_parameters"],
                "metrics_identity": _binary_identity(metrics_path),
                "checkpoint_identity": _binary_identity(checkpoint_path),
            }
    summary = summarize_paired_regression_runs(
        runs["heteroscedastic"], runs["fixed_variance_mse"],
        metric_names=("rmse", "mae", "pearson"),
        left_name="heteroscedastic",
        right_name="fixed_variance_mse",
        difference_key="fixed_minus_heteroscedastic_by_seed",
    )
    output = {
        "protocol_identity": _binary_identity(args.protocol),
        "protocol_claim_boundary": protocol.get("claim_boundary"),
        "run_identities": identities,
        "summary": summary,
        "paired_difference_definition": "fixed_variance_mse minus heteroscedastic",
        "excluded_metrics": protocol.get("excluded_metrics"),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


@torch.no_grad()
def _collect_task_outputs(model, records, task: str, batch_size: int, device: torch.device):
    if task not in model.config.task_names:
        raise ValueError(f"checkpoint does not define task: {task}")
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
    )
    loader = DataLoader(
        RecordDataset(records), batch_size=batch_size, shuffle=False,
        collate_fn=collator, num_workers=0)
    task_index = model.config.task_names.index(task)
    predictions = []
    deviations = []
    labels = []
    model.eval()
    for batch in loader:
        batch = {name: value.to(device) for name, value in batch.items()}
        output = model(batch)
        observed = batch["label_mask"][:, task_index]
        predictions.append(output["mean"][:, task_index][observed].cpu())
        deviations.append(output["standard_deviation"][:, task_index][observed].cpu())
        labels.append(batch["labels"][:, task_index][observed].cpu())
    return torch.cat(predictions), torch.cat(deviations), torch.cat(labels)


def fit_calibration_command(args: argparse.Namespace) -> None:
    source_manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    verify_manifest(args.data, FileManifest(**source_manifest["identity"]))
    split_identity = _manifest_dict(build_manifest(args.split_manifest))
    run_manifest = {
        "source": source_manifest,
        "split": {
            "path": str(Path(args.split_manifest).resolve()),
            "identity": split_identity,
        },
    }
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    if checkpoint.get("data_manifest") != run_manifest:
        raise ValueError("checkpoint was not trained with the supplied source and split identities")
    if args.task not in model.config.task_names:
        raise ValueError(f"checkpoint does not define task: {args.task}")

    split_payload = json.loads(Path(args.split_manifest).read_text(encoding="utf-8"))
    if split_payload.get("format_version") != 1 or not isinstance(split_payload.get("rows"), list):
        raise ValueError("unsupported split manifest format")
    validation_rows = {
        int(row["source_row"]): row["cluster"]
        for row in split_payload["rows"] if row.get("split") == "validation"
    }
    if not validation_rows:
        raise ValueError("validation partition is empty")
    records = read_unikp_json(args.data, source_rows=set(validation_rows)).records
    records = apply_split_manifest(
        records, args.split_manifest, source_identity=source_manifest["identity"])
    if args.conflict_policy == "median":
        records = aggregate_pair_measurements(records)
    validation_records = [record for record in records if record.split == "validation"]
    if not validation_records:
        raise ValueError("validation partition is empty")
    mean, deviation, labels = _collect_task_outputs(
        model, validation_records, args.task, args.batch_size, device)
    if labels.numel() == 0:
        raise ValueError(f"validation partition has no observed labels for task: {args.task}")

    selection = None
    if args.method == "grouped-cv":
        groups = []
        for record in validation_records:
            source_rows = record.source_rows or (record.source_row,)
            record_groups = {validation_rows[source_row] for source_row in source_rows}
            if len(record_groups) != 1:
                raise ValueError("aggregated validation record crosses frozen MMseqs groups")
            groups.append(record_groups.pop())
        selected = select_grouped_gaussian_calibrator(
            mean, deviation, labels, groups,
            n_splits=args.cv_folds, seed=args.cv_seed,
            nll_tolerance=args.nll_tolerance, crps_tolerance=args.crps_tolerance,
        )
        calibrator = selected.calibrator
        selection = selected.to_dict()
        selection["group_manifest_identity"] = split_identity
    else:
        calibrator_types = {
            "identity": IdentityGaussianCalibrator,
            "scalar": ScalarGaussianCalibrator,
            "affine-variance": AffineVarianceGaussianCalibrator,
        }
        calibrator = calibrator_types[args.method].fit(mean, deviation, labels)
    conformal = ()
    if args.method == "grouped-cv":
        calibrated_mean, calibrated_deviation = calibrator.calibrate(mean, deviation)
        conformal = tuple(
            SplitConformalCalibrator.fit(
                calibrated_mean,
                labels,
                alpha=alpha,
                standard_deviation=calibrated_deviation if normalized else None,
                normalized=normalized,
            )
            for alpha in (0.1, 0.05)
            for normalized in (False, True)
        )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    artifact = save_calibration_artifact(
        output_path,
        calibrator,
        identities=_calibration_identities(args.checkpoint, args.task),
        validation_count=int(labels.numel()),
        conformal_calibrators=conformal,
        selection=selection,
    )
    print(json.dumps(artifact, indent=2, sort_keys=True))


def calibrate_uncertainty_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    source_manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    if source_manifest.get("identity") != protocol.get("source_identity"):
        raise ValueError("uncertainty calibration source identity mismatch")
    verify_manifest(args.data, FileManifest(**source_manifest["identity"]))
    split_identity = _manifest_dict(build_manifest(args.split_manifest))
    if split_identity != protocol.get("split_identity"):
        raise ValueError("uncertainty calibration split identity mismatch")
    device = torch.device(args.device)
    model, checkpoint = load_checkpoint(args.checkpoint, device)
    if model.config.uncertainty_mode != "heteroscedastic":
        raise ValueError("uncertainty calibration requires a heteroscedastic checkpoint")
    checkpoint_manifest = checkpoint.get("data_manifest", {})
    if (checkpoint_manifest.get("source", {}).get("identity") != protocol["source_identity"]
            or checkpoint_manifest.get("split", {}).get("identity") != split_identity):
        raise ValueError("uncertainty calibration checkpoint identity mismatch")
    records = apply_split_manifest(
        read_unikp_json(args.data).records, args.split_manifest,
        source_identity=source_manifest["identity"])
    records = aggregate_pair_measurements(records)
    validation_records = [record for record in records if record.split == "validation"]
    test_records = [record for record in records if record.split == "test"]
    validation = _collect_task_outputs(
        model, validation_records, protocol["task"], args.batch_size, device)
    test = _collect_task_outputs(
        model, test_records, protocol["task"], args.batch_size, device)
    scale = scalar_uncertainty_scale(*validation)
    before = uncertainty_metrics(*test)
    after = uncertainty_metrics(test[0], test[1] * scale, test[2])
    output = {
        "protocol_identity": _binary_identity(args.protocol),
        "checkpoint_identity": _binary_identity(args.checkpoint),
        "selected_epoch": checkpoint["epoch"],
        "task": protocol["task"],
        "validation_records": len(validation_records),
        "test_records": len(test_records),
        "standard_deviation_scale": scale,
        "test_before": before,
        "test_after": after,
        "point_metrics": regression_metrics(test[0], test[2]),
        "nominal_gaussian_coverage": protocol["nominal_gaussian_coverage"],
        "claim_boundary": protocol.get("claim_boundary"),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    print(text)


def summarize_uncertainty_calibration_command(args: argparse.Namespace) -> None:
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    seeds = protocol.get("seeds")
    if not isinstance(seeds, list) or not seeds:
        raise ValueError("uncertainty calibration protocol must define seeds")
    protocol_identity = _binary_identity(args.protocol)
    runs = {}
    for seed in seeds:
        path = Path(args.runs_root) / f"seed{seed}.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("protocol_identity") != protocol_identity:
            raise ValueError(f"uncertainty calibration protocol mismatch: seed {seed}")
        runs[seed] = payload
    metrics = ("coverage_1sigma", "coverage_2sigma", "gaussian_nll")
    summary = {"scale": {}, "before": {}, "after": {}, "after_minus_before": {}}
    scales = [runs[seed]["standard_deviation_scale"] for seed in seeds]
    summary["scale"] = {
        "seed_values": scales,
        "mean": statistics.mean(scales),
        "sample_stddev": statistics.stdev(scales) if len(scales) > 1 else None,
    }
    for metric in metrics:
        before = [runs[seed]["test_before"][metric] for seed in seeds]
        after = [runs[seed]["test_after"][metric] for seed in seeds]
        difference = [right - left for left, right in zip(before, after)]
        for name, values in (
            ("before", before), ("after", after), ("after_minus_before", difference),
        ):
            summary[name][metric] = {
                "seed_values": values,
                "mean": statistics.mean(values),
                "sample_stddev": statistics.stdev(values) if len(values) > 1 else None,
            }
    output = {
        "protocol_identity": protocol_identity,
        "protocol_claim_boundary": protocol.get("claim_boundary"),
        "nominal_gaussian_coverage": protocol.get("nominal_gaussian_coverage"),
        "run_identities": {
            str(seed): _binary_identity(Path(args.runs_root) / f"seed{seed}.json")
            for seed in seeds
        },
        "summary": summary,
        "interpretation": (
            "Validation-only scalar calibration is successful only if held-out test NLL and "
            "coverage improve; no test labels were used to fit scales."
        ),
    }
    text = json.dumps(output, indent=2, sort_keys=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    print(text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="BioCandidateRanker research CLI")
    commands = parser.add_subparsers(dest="command", required=True)
    audit = commands.add_parser("audit-data")
    audit.add_argument("--data", required=True)
    audit.add_argument("--output")
    audit.set_defaults(func=audit_command)

    audit_lunzer = commands.add_parser("audit-lunzer")
    audit_lunzer.add_argument("--data", required=True)
    audit_lunzer.add_argument("--unikp-data")
    audit_lunzer.add_argument("--unikp-manifest")
    audit_lunzer.add_argument("--output")
    audit_lunzer.set_defaults(func=audit_lunzer_command)

    audit_enzengdb = commands.add_parser("audit-enzengdb")
    audit_enzengdb.add_argument("--archive", required=True)
    audit_enzengdb.add_argument("--experiments", required=True)
    audit_enzengdb.add_argument("--unikp-data")
    audit_enzengdb.add_argument("--unikp-manifest")
    audit_enzengdb.add_argument("--min-campaign-size", type=int, default=20)
    audit_enzengdb.add_argument("--max-records-per-campaign", type=int, default=2000)
    audit_enzengdb.add_argument("--selection-output")
    audit_enzengdb.add_argument("--output")
    audit_enzengdb.add_argument("--mmseqs")
    audit_enzengdb.add_argument(
        "--homology-output-dir", default="artifacts/external/enzengdb-v1/homology")
    audit_enzengdb.add_argument("--min-identity", type=float, default=0.3)
    audit_enzengdb.add_argument("--coverage", type=float, default=0.8)
    audit_enzengdb.add_argument("--threads", type=int, default=4)
    audit_enzengdb.add_argument("--exclude-homology-hits", action="store_true")
    audit_enzengdb.set_defaults(func=audit_enzengdb_command)

    train = commands.add_parser("train")
    train.add_argument("--data", required=True)
    train.add_argument("--data-format", choices=DATA_FORMATS, default="unikp")
    train.add_argument("--tasks", nargs="+", choices=KINETIC_TASKS)
    train.add_argument("--manifest", required=True)
    train.add_argument("--output-dir", default="artifacts/training")
    train.add_argument("--split-strategy", choices=(
        "protein_cold", "scaffold_cold", "double_cold"), default="protein_cold")
    train.add_argument("--split-manifest")
    train.add_argument("--epochs", type=int, default=10)
    train.add_argument("--batch-size", type=int, default=32)
    train.add_argument("--lr", type=float, default=2e-4)
    train.add_argument("--seed", type=int, default=42)
    train.add_argument("--max-records", type=int)
    train.add_argument("--conflict-policy", choices=("error", "keep", "median"), default="error")
    train.add_argument("--resume")
    train.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    train.add_argument("--d-model", type=int, default=128)
    train.add_argument("--num-heads", type=int, default=8)
    train.add_argument("--protein-layers", type=int, default=3)
    train.add_argument("--molecule-layers", type=int, default=4)
    train.add_argument("--fusion-layers", type=int, default=2)
    train.add_argument("--chunk-size", type=int, default=128)
    train.add_argument("--context-buckets", type=int, default=4096)
    train.add_argument("--disable-protein", action="store_true")
    train.add_argument("--disable-molecule", action="store_true")
    train.add_argument("--disable-context", action="store_true")
    train.add_argument("--shared-task-query", action="store_true")
    train.add_argument(
        "--fusion-mode", choices=("task_query", "late_concat"), default="task_query")
    train.add_argument(
        "--protein-encoder", choices=("chunk_transformer", "global_mean", "esm2"),
        default="chunk_transformer")
    train.add_argument(
        "--uncertainty-mode", choices=("heteroscedastic", "fixed_variance_mse"),
        default="heteroscedastic")
    train.add_argument("--esm2-model", default="esm2_t6_8M_UR50D",
                       help="ESM-2 pretrained model name (used with --protein-encoder esm2)")
    train.add_argument("--esm2-unfreeze", action="store_true",
                       help="Fine-tune ESM-2 backbone (default: frozen)")
    train.add_argument("--warmup-epochs", type=int, default=0,
                       help="Linear LR warmup epochs before cosine decay")
    train.add_argument("--cosine-schedule", action="store_true",
                       help="Enable cosine annealing LR schedule")
    train.add_argument("--grad-accum-steps", type=int, default=1,
                       help="Gradient accumulation steps (effective batch = batch_size * steps)")
    train.add_argument("--patience", type=int, default=0,
                       help="Early stopping patience (0 = disabled)")
    train.set_defaults(func=train_command)

    homology = commands.add_parser("build-homology-split")
    homology.add_argument("--data", required=True)
    homology.add_argument("--manifest", required=True)
    homology.add_argument("--mmseqs")
    homology.add_argument("--cluster-tsv")
    homology.add_argument("--output-dir", default="artifacts/homology")
    homology.add_argument("--min-identity", type=float, default=0.3)
    homology.add_argument("--coverage", type=float, default=0.8)
    homology.add_argument("--threads", type=int, default=4)
    homology.add_argument("--seed", type=int, default=42)
    homology.set_defaults(func=homology_command)

    rank_split = commands.add_parser("build-enzengdb-rank-split")
    rank_split.add_argument("--archive", required=True)
    rank_split.add_argument("--experiments", required=True)
    rank_split.add_argument("--test-selection", required=True)
    rank_split.add_argument("--mmseqs", required=True)
    rank_split.add_argument("--cluster-tsv")
    rank_split.add_argument("--output-dir", default="artifacts/external/enzengdb-v1/rank-split")
    rank_split.add_argument("--min-identity", type=float, default=0.3)
    rank_split.add_argument("--coverage", type=float, default=0.8)
    rank_split.add_argument("--threads", type=int, default=4)
    rank_split.add_argument("--validation-fraction", type=float, default=0.15)
    rank_split.add_argument("--min-campaign-size", type=int, default=20)
    rank_split.add_argument("--max-records-per-campaign", type=int, default=500)
    rank_split.add_argument("--seed", type=int, default=42)
    rank_split.set_defaults(func=build_enzengdb_rank_split_command)

    train_ranker = commands.add_parser("train-enzengdb-ranker")
    train_ranker.add_argument("--archive", required=True)
    train_ranker.add_argument("--experiments", required=True)
    train_ranker.add_argument("--split-manifest", required=True)
    train_ranker.add_argument("--initialize-from", required=True)
    train_ranker.add_argument("--output-dir", default="artifacts/enzengdb-ranker")
    train_ranker.add_argument("--epochs", type=int, default=3)
    train_ranker.add_argument("--batch-size", type=int, default=64)
    train_ranker.add_argument("--lr", type=float, default=1e-4)
    train_ranker.add_argument("--seed", type=int, default=42)
    train_ranker.add_argument("--disable-protein", action="store_true")
    train_ranker.add_argument("--disable-molecule", action="store_true")
    train_ranker.add_argument("--disable-context", action="store_true")
    train_ranker.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    train_ranker.set_defaults(func=train_enzengdb_ranker_command)

    cold_split = commands.add_parser("build-cold-split")
    cold_split.add_argument("--data", required=True)
    cold_split.add_argument("--manifest", required=True)
    cold_split.add_argument("--strategy", choices=(
        "protein_cold", "scaffold_cold", "double_cold"), required=True)
    cold_split.add_argument("--output", required=True)
    cold_split.add_argument("--seed", type=int, default=42)
    cold_split.set_defaults(func=cold_split_command)

    predict = commands.add_parser("predict")
    predict.add_argument("--checkpoint", required=True)
    predict.add_argument("--input", required=True)
    predict.add_argument("--output", required=True)
    predict.add_argument("--calibration-artifact")
    predict.add_argument(
        "--fba-features",
        help="JSON with FBAFeatureMetadata (schema_version, feature_ids, model_id, solver_id, "
             "objective_id, condition_id); required when candidates supply fba_context")
    predict.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    predict.set_defaults(func=predict_command)

    fit_calibration = commands.add_parser("fit-calibration")
    fit_calibration.add_argument("--checkpoint", required=True)
    fit_calibration.add_argument("--data", required=True)
    fit_calibration.add_argument("--manifest", required=True)
    fit_calibration.add_argument("--split-manifest", required=True)
    fit_calibration.add_argument("--task", required=True)
    fit_calibration.add_argument(
        "--method", choices=("identity", "scalar", "affine-variance", "grouped-cv"),
        required=True)
    fit_calibration.add_argument("--cv-folds", type=int, default=5)
    fit_calibration.add_argument("--cv-seed", type=int, default=0)
    fit_calibration.add_argument("--nll-tolerance", type=float, default=1e-12)
    fit_calibration.add_argument("--crps-tolerance", type=float, default=1e-12)
    fit_calibration.add_argument("--conflict-policy", choices=("keep", "median"), default="median")
    fit_calibration.add_argument("--batch-size", type=int, default=64)
    fit_calibration.add_argument("--output", required=True)
    fit_calibration.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    fit_calibration.set_defaults(func=fit_calibration_command)

    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("--checkpoint", required=True)
    evaluate.add_argument("--data", required=True)
    evaluate.add_argument("--data-format", choices=DATA_FORMATS, default="unikp")
    evaluate.add_argument("--tasks", nargs="+", choices=KINETIC_TASKS)
    evaluate.add_argument("--manifest", required=True)
    evaluate.add_argument("--split-manifest", required=True)
    evaluate.add_argument("--partition", choices=("validation", "test"), default="test")
    evaluate.add_argument("--conflict-policy", choices=("keep", "median"), default="median")
    evaluate.add_argument("--batch-size", type=int, default=64)
    evaluate.add_argument("--output", required=True)
    evaluate.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    evaluate.set_defaults(func=evaluate_command)

    classical = commands.add_parser("evaluate-classical-baselines")
    classical.add_argument("--data", required=True)
    classical.add_argument("--manifest", required=True)
    classical.add_argument("--split-manifest", required=True)
    classical.add_argument("--output", required=True)
    classical.set_defaults(func=evaluate_classical_baselines_command)

    summarize_fusion = commands.add_parser("summarize-fusion-runs")
    summarize_fusion.add_argument("--protocol", required=True)
    summarize_fusion.add_argument("--runs-root", required=True)
    summarize_fusion.add_argument("--output", required=True)
    summarize_fusion.set_defaults(func=summarize_fusion_runs_command)

    feature_mlp = commands.add_parser("train-feature-mlp")
    feature_mlp.add_argument("--protocol", required=True)
    feature_mlp.add_argument("--feature", required=True)
    feature_mlp.add_argument("--seed", type=int, required=True)
    feature_mlp.add_argument("--data", required=True)
    feature_mlp.add_argument("--manifest", required=True)
    feature_mlp.add_argument("--split-manifest", required=True)
    feature_mlp.add_argument("--output-dir", required=True)
    feature_mlp.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    feature_mlp.set_defaults(func=train_feature_mlp_command)

    summarize_feature_mlp = commands.add_parser("summarize-feature-mlp-runs")
    summarize_feature_mlp.add_argument("--protocol", required=True)
    summarize_feature_mlp.add_argument("--runs-root", required=True)
    summarize_feature_mlp.add_argument("--output", required=True)
    summarize_feature_mlp.set_defaults(func=summarize_feature_mlp_runs_command)

    summarize_ablations = commands.add_parser("summarize-internal-ablations")
    summarize_ablations.add_argument("--protocol", required=True)
    summarize_ablations.add_argument("--baseline-metrics", required=True)
    summarize_ablations.add_argument("--runs-root", required=True)
    summarize_ablations.add_argument("--output", required=True)
    summarize_ablations.set_defaults(func=summarize_internal_ablations_command)

    summarize_uncertainty = commands.add_parser("summarize-uncertainty-ablation")
    summarize_uncertainty.add_argument("--protocol", required=True)
    summarize_uncertainty.add_argument("--baseline-runs-root", required=True)
    summarize_uncertainty.add_argument("--ablation-runs-root", required=True)
    summarize_uncertainty.add_argument("--output", required=True)
    summarize_uncertainty.set_defaults(func=summarize_uncertainty_ablation_command)

    calibrate_uncertainty = commands.add_parser("calibrate-uncertainty")
    calibrate_uncertainty.add_argument("--protocol", required=True)
    calibrate_uncertainty.add_argument("--checkpoint", required=True)
    calibrate_uncertainty.add_argument("--data", required=True)
    calibrate_uncertainty.add_argument("--manifest", required=True)
    calibrate_uncertainty.add_argument("--split-manifest", required=True)
    calibrate_uncertainty.add_argument("--batch-size", type=int, default=64)
    calibrate_uncertainty.add_argument("--output", required=True)
    calibrate_uncertainty.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    calibrate_uncertainty.set_defaults(func=calibrate_uncertainty_command)

    summarize_calibration = commands.add_parser("summarize-uncertainty-calibration")
    summarize_calibration.add_argument("--protocol", required=True)
    summarize_calibration.add_argument("--runs-root", required=True)
    summarize_calibration.add_argument("--output", required=True)
    summarize_calibration.set_defaults(func=summarize_uncertainty_calibration_command)

    evaluate_enzengdb = commands.add_parser("evaluate-enzengdb")
    evaluate_enzengdb.add_argument("--checkpoint", required=True)
    evaluate_enzengdb.add_argument("--archive", required=True)
    evaluate_enzengdb.add_argument("--experiments", required=True)
    evaluate_enzengdb.add_argument("--selection-manifest", required=True)
    evaluate_enzengdb.add_argument("--task", default="log10_kcat")
    evaluate_enzengdb.add_argument("--batch-size", type=int, default=64)
    evaluate_enzengdb.add_argument("--output", required=True)
    evaluate_enzengdb.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    evaluate_enzengdb.set_defaults(func=evaluate_enzengdb_command)

    evaluate_lunzer = commands.add_parser("evaluate-lunzer")
    evaluate_lunzer.add_argument("--checkpoint", required=True)
    evaluate_lunzer.add_argument("--data", required=True)
    evaluate_lunzer.add_argument("--selection-manifest", required=True)
    evaluate_lunzer.add_argument("--task", default="log10_kcat")
    evaluate_lunzer.add_argument("--batch-size", type=int, default=64)
    evaluate_lunzer.add_argument("--output", required=True)
    evaluate_lunzer.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu")
    evaluate_lunzer.set_defaults(func=evaluate_lunzer_command)

    evaluate_lunzer_dlkcat = commands.add_parser("evaluate-lunzer-dlkcat-output")
    evaluate_lunzer_dlkcat.add_argument("--data", required=True)
    evaluate_lunzer_dlkcat.add_argument("--selection-manifest", required=True)
    evaluate_lunzer_dlkcat.add_argument("--predictions", required=True)
    evaluate_lunzer_dlkcat.add_argument("--dlkcat-checkpoint", required=True)
    evaluate_lunzer_dlkcat.add_argument("--output", required=True)
    evaluate_lunzer_dlkcat.set_defaults(func=evaluate_lunzer_dlkcat_output_command)

    summarize_runs = commands.add_parser("summarize-enzengdb-runs")
    summarize_runs.add_argument("--inputs", nargs="+", required=True)
    summarize_runs.add_argument("--bootstrap-samples", type=int, default=10_000)
    summarize_runs.add_argument("--seed", type=int, default=42)
    summarize_runs.add_argument("--output", required=True)
    summarize_runs.set_defaults(func=summarize_enzengdb_runs_command)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
