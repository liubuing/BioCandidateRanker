from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import platform
from dataclasses import asdict
from pathlib import Path

import torch
from torch.utils.data import DataLoader

from .config import ModelConfig
from .data import (
    CANONICAL_LOG10_FLUX_UNIT,
    EnzymeSubstrateCollator,
    EnzymeSubstrateRecord,
    FBAFeatureMetadata,
    FluxLabelMetadata,
    FluxLabelType,
    build_manifest,
)
from .model.ranker import BioCandidateRanker
from .training import RecordDataset, run_epoch, save_checkpoint


FORMAT_VERSION = 1
TARGET_REACTION_ID = "r_2111"
FEATURE_IDS = ("medium_lower_bound:r_1714", "medium_lower_bound:r_1992")
GLUCOSE_LEVELS = tuple(round(0.15 + 0.075 * index, 3) for index in range(12))
SPLITS = ("train",) * 8 + ("validation",) * 2 + ("test",) * 2
CLAIM_BOUNDARY = (
    "Targets are simulated FBA biomass fluxes under declared model assumptions. "
    "They are not experimental flux measurements or publication validation."
)
OUTPUT_NAMES = (
    "blocker.json", "manifest.json", "simulated_flux.jsonl", "smoke.pt",
    "smoke_metrics.json", "split.json",
)


def _canonical_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _binary_identity(path: Path) -> dict:
    return {"sha256": _file_hash(path), "size_bytes": path.stat().st_size}


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _clear_pipeline_outputs(output: Path) -> None:
    for name in OUTPUT_NAMES:
        path = output / name
        if path.is_file():
            path.unlink()


def _identities(model_path: Path, model, cobra_version: str, glpk_version: str) -> dict:
    objective = {
        "direction": "max",
        "reaction_coefficients": {TARGET_REACTION_ID: 1.0},
    }
    medium = {
        "schema": "simulated-flux-medium-v1",
        "glucose_exchange": "r_1714",
        "glucose_lower_bounds": [-value for value in GLUCOSE_LEVELS],
        "oxygen_exchange": "r_1992",
        "oxygen_lower_bound": -2.0,
        "all_other_bounds": "unchanged_from_model_asset",
    }
    solver = {
        "cobra_version": cobra_version,
        "glpk_version": glpk_version,
        "interface": model.solver.interface.__name__,
        "python": platform.python_version(),
    }
    return {
        "model": {
            "sha256": _file_hash(model_path),
            "size_bytes": model_path.stat().st_size,
        },
        "solver": {"sha256": _canonical_hash(solver), "specification": solver},
        "objective": {"sha256": _canonical_hash(objective), "specification": objective},
        "medium": {"sha256": _canonical_hash(medium), "specification": medium},
    }


def _record(row: dict) -> EnzymeSubstrateRecord:
    fba = row["fba"]
    label = row["label"]
    return EnzymeSubstrateRecord(
        # Neutral carriers satisfy the multimodal interface; both modalities are disabled.
        sequence="X",
        substrate_smiles="O",
        organism="Saccharomyces cerevisiae (model context only)",
        ec="not_applicable",
        enzyme_type="simulated_fba_condition",
        candidate_id=row["condition_id"],
        reaction=label["target_reaction_id"],
        fba_context=tuple(row["features"]),
        fba_feature_metadata=FBAFeatureMetadata(
            schema_version=fba["schema_version"],
            feature_ids=tuple(fba["feature_ids"]),
            model_id=fba["model_id"],
            solver_id=fba["solver_id"],
            objective_id=fba["objective_id"],
            condition_id=row["condition_id"],
        ),
        log10_flux=row["log10_flux"],
        flux_label_metadata=FluxLabelMetadata(
            label_type=label["type"],
            target_reaction_id=label["target_reaction_id"],
            canonical_unit=label["canonical_unit"],
            provenance_json=json.dumps(label["provenance"]),
        ),
        source_dataset="Yeast-MetaTwin simulated FBA smoke v1",
        source_row=row["source_row"],
        split=row["split"],
    )


def read_simulated_flux_jsonl(path: str | Path) -> tuple[EnzymeSubstrateRecord, ...]:
    records = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("format_version") != FORMAT_VERSION:
                raise ValueError(f"line {line_number}: unsupported simulated-flux format")
            records.append(_record(row))
    if not records:
        raise ValueError("simulated-flux dataset is empty")
    return tuple(records)


def generate_dataset(model_path: str | Path, output_dir: str | Path) -> dict:
    source = Path(model_path).resolve()
    if not source.is_file():
        raise RuntimeError(f"Yeast-MetaTwin model asset is not a file: {source}")
    try:
        import cobra
        import swiglpk
    except ImportError as exc:
        raise RuntimeError("missing required simulation dependency: cobra with GLPK") from exc

    initial_hash = _file_hash(source)
    model = cobra.io.load_yaml_model(source)
    model.solver = "glpk"
    for reaction_id in ("r_1714", "r_1992", TARGET_REACTION_ID):
        if reaction_id not in model.reactions:
            raise RuntimeError(f"required reaction is absent from model: {reaction_id}")
    model.objective = TARGET_REACTION_ID
    identities = _identities(
        source, model, importlib.metadata.version("cobra"), swiglpk.glp_version())
    model_id = "sha256:" + identities["model"]["sha256"]
    solver_id = "sha256:" + identities["solver"]["sha256"]
    objective_id = "sha256:" + identities["objective"]["sha256"]
    medium_hash = identities["medium"]["sha256"]

    rows = []
    for index, (glucose, split) in enumerate(zip(GLUCOSE_LEVELS, SPLITS), 1):
        with model:
            model.reactions.get_by_id("r_1714").lower_bound = -glucose
            model.reactions.get_by_id("r_1992").lower_bound = -2.0
            solution = model.optimize()
        if solution.status != "optimal":
            raise RuntimeError(f"condition {index} did not solve optimally: {solution.status}")
        flux = float(solution.fluxes[TARGET_REACTION_ID])
        if not math.isfinite(flux) or flux <= 0:
            raise RuntimeError(f"condition {index} produced nonpositive target flux: {flux}")
        condition = {"glucose_lower_bound": -glucose, "oxygen_lower_bound": -2.0}
        condition_id = "sha256:" + _canonical_hash(condition)
        provenance = {
            "claim_boundary": CLAIM_BOUNDARY,
            "condition": condition,
            "medium_sha256": medium_hash,
            "model_sha256": identities["model"]["sha256"],
            "objective_sha256": identities["objective"]["sha256"],
            "solver_sha256": identities["solver"]["sha256"],
        }
        row = {
            "format_version": FORMAT_VERSION,
            "source_row": index,
            "condition_id": condition_id,
            "split": split,
            "features": [-glucose, -2.0],
            "fba": {
                "schema_version": 1,
                "feature_ids": list(FEATURE_IDS),
                "model_id": model_id,
                "solver_id": solver_id,
                "objective_id": objective_id,
            },
            "log10_flux": math.log10(flux),
            "label": {
                "type": FluxLabelType.SIMULATED.value,
                "target_reaction_id": TARGET_REACTION_ID,
                "canonical_unit": CANONICAL_LOG10_FLUX_UNIT,
                "provenance": provenance,
            },
        }
        _record(row)
        rows.append(row)

    if _file_hash(source) != initial_hash:
        raise RuntimeError("read-only source model changed during simulation")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    dataset_path = output / "simulated_flux.jsonl"
    dataset_path.write_text(
        "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in rows),
        encoding="utf-8",
    )
    dataset_identity = asdict(build_manifest(dataset_path))
    split = {
        "format_version": 1,
        "dataset_identity": dataset_identity,
        "method": "fixed_label_independent_condition_order_v1",
        "counts": {name: SPLITS.count(name) for name in sorted(set(SPLITS))},
        "rows": [{"condition_id": row["condition_id"], "split": row["split"]} for row in rows],
    }
    _write_json(output / "split.json", split)
    manifest = {
        "format_version": 1,
        "task": "simulated_flux",
        "target": "log10_flux",
        "target_reaction_id": TARGET_REACTION_ID,
        "canonical_unit": CANONICAL_LOG10_FLUX_UNIT,
        "claim_boundary": CLAIM_BOUNDARY,
        "source_model_path": str(source),
        "source_access": "read_only",
        "identities": identities,
        "dataset_identity": dataset_identity,
        "split_identity": asdict(build_manifest(output / "split.json")),
        "leakage_control": {
            "target_reaction_excluded_from_feature_ids": TARGET_REACTION_ID not in FEATURE_IDS,
            "features_are_declared_medium_bounds_only": True,
            "solved_fluxes_in_features": False,
        },
    }
    _write_json(output / "manifest.json", manifest)
    return manifest


def train_cpu_smoke(output_dir: str | Path, *, seed: int = 7) -> dict:
    output = Path(output_dir)
    records = read_simulated_flux_jsonl(output / "simulated_flux.jsonl")
    partitions = {name: [record for record in records if record.split == name]
                  for name in ("train", "validation", "test")}
    if any(not values for values in partitions.values()):
        raise RuntimeError("simulated-flux split contains an empty partition")
    torch.manual_seed(seed)
    config = ModelConfig(
        d_model=16,
        num_heads=2,
        protein_layers=1,
        molecule_layers=1,
        fusion_layers=1,
        protein_chunk_size=8,
        context_buckets=32,
        fba_context_dim=len(FEATURE_IDS),
        dropout=0.0,
        use_protein=False,
        use_molecule=False,
        use_context=True,
        uncertainty_mode="fixed_variance_mse",
        task_names=("log10_flux",),
    )
    device = torch.device("cpu")
    model = BioCandidateRanker(config).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    collator = EnzymeSubstrateCollator(
        context_buckets=config.context_buckets,
        fba_context_dim=config.fba_context_dim,
        task_names=config.task_names,
    )
    loaders = {
        name: DataLoader(RecordDataset(values), batch_size=4, shuffle=False, collate_fn=collator)
        for name, values in partitions.items()
    }
    train_metrics = run_epoch(model, loaders["train"], device, optimizer)
    validation_metrics = run_epoch(model, loaders["validation"], device)
    test_metrics = run_epoch(model, loaders["test"], device)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    metrics = {
        "status": "smoke_passed",
        "device": "cpu",
        "epochs": 1,
        "seed": seed,
        "task": "simulated_flux",
        "claim_boundary": CLAIM_BOUNDARY,
        "train_loss": train_metrics.loss,
        "validation_loss": validation_metrics.loss,
        "test_loss": test_metrics.loss,
        "partition_observations": {name: len(values) for name, values in partitions.items()},
    }
    save_checkpoint(
        output / "smoke.pt",
        model,
        optimizer,
        epoch=1,
        data_manifest={
            "manifest_identity": _binary_identity(output / "manifest.json"),
            "dataset_identity": manifest["dataset_identity"],
            "split_identity": manifest["split_identity"],
            "simulation_identities": manifest["identities"],
            "claim_boundary": CLAIM_BOUNDARY,
        },
        metrics=metrics,
    )
    _write_json(output / "smoke_metrics.json", metrics)
    return metrics


def run_pipeline(model_path: str | Path, output_dir: str | Path, *, seed: int = 7) -> dict:
    output = Path(output_dir)
    _clear_pipeline_outputs(output)
    try:
        manifest = generate_dataset(model_path, output)
        metrics = train_cpu_smoke(output, seed=seed)
    except Exception as exc:
        _clear_pipeline_outputs(output)
        output.mkdir(parents=True, exist_ok=True)
        blocker = {
            "status": "blocked",
            "task": "simulated_flux",
            "claim_boundary": CLAIM_BOUNDARY,
            "model_path": str(Path(model_path).resolve()),
            "reason": f"{type(exc).__name__}: {exc}",
        }
        _write_json(output / "blocker.json", blocker)
        raise RuntimeError(blocker["reason"]) from exc
    return {"manifest": manifest, "metrics": metrics}
