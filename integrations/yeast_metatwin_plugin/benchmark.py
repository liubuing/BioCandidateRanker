"""Audit BioCandidateRanker benchmark independence and inference readiness.

Mirrors the structure of 04_prediction_plugins/UniKP/benchmark.py so the
platform's registry and gates can consume it unchanged.

Honest status, and why
    The benchmark gate is BLOCKED. BioCandidateRanker's internal homology-cold
    test is not an independent benchmark: its labels come from a development
    corpus that the project cannot redistribute, so no third party can rebuild
    the evaluation from public artifacts. The project's own reassessment reaches
    the same conclusion, and its paired family intervals against CatPred,
    CataPro and UniKP all cross zero, so no advantage is claimed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
PLUGIN = "BioCandidateRanker"
PLUGIN_VERSION = "software-implementation-v6+esm2-t6-seed7"
CAPABILITY = "kcat_prediction"
SOURCE_REPOSITORY = "https://github.com/liubuing/BioCandidateRanker"
RELEASE_ID = "software-implementation-v6"
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

# The adapter lives in plugin_runtime/; the package directory is its parent's
# sibling, and the checkpoint sits under this plugin's assets/.
CHECKPOINT_RELATIVE = "assets/esm2-t6-opt-seed7.best.pt"
ADAPTER_RELATIVE = "../plugin_runtime/biocandidate.py"

ASSETS = (CHECKPOINT_RELATIVE, "benchmark.py", ADAPTER_RELATIVE)

# The checkpoint only has log10_kcat trained; the adapter must not emit the
# other declared tasks, and the manifest states that so nobody assumes four.
TRAINED_TASKS = ("log10_kcat",)
UNTRAINED_TASKS = ("log10_km", "log10_activity", "log10_flux")

DEVELOPMENT_CORPUS = {
    "capability": CAPABILITY,
    "path": "(external) UniKP/DLKcat kcat development corpus",
    "source_role": "training_corpus",
    "overlap_evidence": (
        "all shippable BioCandidateRanker checkpoints were trained on this corpus; "
        "the frozen homology-cold split (configs/software_implementation_release_v6.json) "
        "partitions it into 13,157 / 1,640 / 1,646 rows"
    ),
    "exists": True,
    "size_bytes": 12132719,
    "sha256": "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719",
    "source_repository": SOURCE_REPOSITORY,
    "license": "unknown; redistribution not verified (configs/unikp_source_manifest.json)",
    "license_verified": False,
    "training_overlap": True,
    "independent_benchmark_eligible": False,
}

BLOCKERS = (
    "the internal homology-cold test is an already-observed split of the training corpus, "
    "not an independent benchmark",
    "the development corpus is not redistributable, so no third party can rebuild the "
    "evaluation from public artifacts",
    "paired family-bootstrap intervals against CatPred, CataPro and UniKP all cross zero, "
    "so no performance advantage is established",
    "no independently sourced, versioned holdout dataset is present",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def asset_record(relative_path: str) -> dict[str, Any]:
    path = (ROOT / relative_path).resolve()
    if not path.is_file():
        return {"path": relative_path, "exists": False, "size_bytes": None, "sha256": None}
    return {
        "path": relative_path,
        "exists": True,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def package_versions() -> tuple[dict[str, str | None], list[str]]:
    """Report the versions the adapter needs, plus anything missing."""
    versions: dict[str, str | None] = {"python": platform.python_version()}
    missing: list[str] = []
    for name in ("torch", "numpy", "rdkit", "esm"):
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", "unknown")
        except ImportError:
            versions[name] = None
            missing.append(name)
    return versions, missing


def build_manifests(
    smoke_status: str = "not_run",
    smoke_detail: str = "No successful end-to-end fixed-input inference result was supplied",
) -> tuple[dict[str, Any], dict[str, Any]]:
    generated = datetime.now(timezone.utc).isoformat()
    assets = tuple(asset_record(path) for path in ASSETS)
    versions, missing_packages = package_versions()

    blockers = list(BLOCKERS)
    if missing_packages:
        blockers.append("missing runtime packages: " + ", ".join(missing_packages))
    if smoke_status != "passed":
        blockers.append(f"inference smoke is {smoke_status}: {smoke_detail}")

    checkpoint = next((row for row in assets if row["path"] == CHECKPOINT_RELATIVE), None)
    if checkpoint and not checkpoint["exists"]:
        blockers.append(f"checkpoint asset missing: {CHECKPOINT_RELATIVE}")

    benchmark = {
        "schema_version": "1.0.0",
        "generated_at": generated,
        "plugin": PLUGIN,
        "status": "blocked",
        "independent": False,
        "metrics_published": False,
        "metrics": {},
        "blockers": list(BLOCKERS),
        "source": {
            "repository": SOURCE_REPOSITORY,
            "release_id": RELEASE_ID,
            "license": "MIT",
            "paper_target": "PLOS Computational Biology (not submitted)",
        },
        "datasets": [dict(DEVELOPMENT_CORPUS)],
    }

    readiness = {
        "schema_version": "1.0.0",
        "generated_at": generated,
        "plugin": PLUGIN,
        "status": "blocked",
        "inference_gate": "blocked",
        "inference_smoke": {"status": smoke_status, "detail": smoke_detail},
        "benchmark_gate": "blocked",
        "blockers": blockers,
        "runtime": {
            "python": versions.pop("python"),
            "executable": sys.executable,
            "packages": versions,
        },
        "assets": list(assets),
        "capabilities": (CAPABILITY,),
        "trained_tasks": list(TRAINED_TASKS),
        "untrained_declared_tasks": list(UNTRAINED_TASKS),
        "output_contract": {
            # BioCandidateRanker ships no applicability model, so this field is
            # always reported as not_available rather than asserted.
            "ood_field": "applicability.ood",
            "ood_status": "not_implemented",
            "truncation_fields": ["transforms[].truncated", "transforms[].strategy"],
            "uncertainty_field": "uncertainty.value",
            "uncertainty_method": "heteroscedastic_gaussian_log10_1sigma",
            "uncertainty_calibrated": False,
            # No tree/ensemble members: both tree fields stay null by contract.
            "tree_member_count_field": None,
            "tree_member_interval_field": None,
            "tree_interval_calibrated": False,
        },
        "benchmark_manifest": "benchmark_manifest.json",
        "notes": (
            "Loaded through a torch.serialization.safe_globals allowlist because the "
            "platform pins torch >= 2.6 and the checkpoints predate it. Requires "
            "BIOCANDIDATE_RANKER_ROOT (or ranker_root.txt) to locate the package."
        ),
    }

    benchmark_path = ROOT / "benchmark_manifest.json"
    if benchmark_path.is_file():
        readiness["benchmark_manifest_sha256"] = sha256_file(benchmark_path)
    return benchmark, readiness


def write_manifests(
    benchmark_path: Path,
    readiness_path: Path,
    smoke_status: str,
    smoke_detail: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    benchmark, readiness = build_manifests(smoke_status, smoke_detail)
    # newline="\n" keeps the manifests LF on Windows; Path.write_text would
    # otherwise translate to os.linesep and diverge from a fresh git checkout.
    benchmark_path.write_text(
        json.dumps(benchmark, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    readiness["benchmark_manifest_sha256"] = sha256_file(benchmark_path)
    readiness_path.write_text(
        json.dumps(readiness, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return benchmark, readiness


def validate_manifest(manifest: dict[str, Any]) -> None:
    if manifest["status"] == "blocked" and manifest.get("metrics"):
        raise ValueError("blocked manifests cannot publish metrics")
    for record in manifest.get("datasets", []) + manifest.get("assets", []):
        digest = record.get("sha256")
        if record.get("exists") and not (isinstance(digest, str) and SHA256_PATTERN.fullmatch(digest)):
            raise ValueError(f"missing or invalid sha256 for {record['path']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit BioCandidateRanker benchmark independence and inference readiness."
    )
    parser.add_argument("--benchmark-manifest", type=Path, default=ROOT / "benchmark_manifest.json")
    parser.add_argument("--readiness-manifest", type=Path, default=ROOT / "readiness_manifest.json")
    parser.add_argument("--inference-smoke-status", choices=("passed", "blocked", "not_run"), default="not_run")
    parser.add_argument("--inference-smoke-detail",
                        default="No successful end-to-end fixed-input inference result was supplied")
    args = parser.parse_args(argv)

    benchmark, readiness = write_manifests(
        args.benchmark_manifest, args.readiness_manifest,
        args.inference_smoke_status, args.inference_smoke_detail,
    )
    validate_manifest(benchmark)
    validate_manifest(readiness)
    print(json.dumps({
        "benchmark": benchmark["status"],
        "readiness": readiness["status"],
        "blockers": readiness["blockers"],
    }, indent=2))
    return 0 if benchmark["status"] == readiness["status"] == "ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
