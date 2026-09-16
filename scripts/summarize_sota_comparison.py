"""Aggregate Mode B SOTA comparison results into one summary receipt.

Scans per-seed test_metrics.json files under each predictor's run directory (disk is the
source of truth; driver summary.json files may cover only a subset of seeds after
incremental runs), recomputes the BioCandidateRanker ESM-2 reference numbers from its
checkpoint metrics, and writes a manuscript-ready comparison table with the frozen split
identity. Governed by configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOTA_ROOT = ROOT / "artifacts" / "external" / "sota-homology-cold"
PROTOCOL_SEEDS = (7, 42, 123)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_sd(values: list[float]) -> dict:
    return {
        "mean": statistics.mean(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def collect_seed_metrics(run_dir: Path) -> dict[int, dict]:
    metrics: dict[int, dict] = {}
    for path in sorted(run_dir.glob("seed*/test_metrics.json")):
        seed = int(path.parent.name.removeprefix("seed"))
        metrics[seed] = json.loads(path.read_text(encoding="utf-8"))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=SOTA_ROOT / "summary.json")
    args = parser.parse_args()

    rows: list[dict] = []

    esm2 = {}
    for path in sorted(ROOT.glob("artifacts/esm2-t6-opt-seed*/test_metrics.json")):
        seed = int(path.parent.name.removeprefix("esm2-t6-opt-seed"))
        esm2[seed] = json.loads(path.read_text(encoding="utf-8"))["model_metrics"]["log10_kcat"]
    rows.append({
        "model": "BioCandidateRanker ESM-2 (ours, reference)",
        "kind": "internal",
        "seeds": sorted(esm2),
        "rmse": mean_sd([esm2[seed]["rmse"] for seed in sorted(esm2)]),
        "mae": mean_sd([esm2[seed]["mae"] for seed in sorted(esm2)]),
        "pearson": mean_sd([esm2[seed]["pearson"] for seed in sorted(esm2)]),
        "source": "artifacts/esm2-t6-opt-seed{7,42,123}/test_metrics.json",
    })

    for predictor, run_dir in (
        ("UniKP Mode B (retrained)", SOTA_ROOT / "unikp-mode-b"),
        ("DLKcat Mode B (retrained)", SOTA_ROOT / "dlkcat-mode-b"),
    ):
        seeds_metrics = collect_seed_metrics(run_dir)
        if not seeds_metrics:
            print(f"no seed results, skipped: {run_dir}")
            continue
        seeds = sorted(seeds_metrics)
        rows.append({
            "model": predictor,
            "kind": "sota-mode-b",
            "seeds": seeds,
            "rmse": mean_sd([seeds_metrics[seed]["rmse"] for seed in seeds]),
            "mae": mean_sd([seeds_metrics[seed]["mae"] for seed in seeds]),
            "pearson": mean_sd([seeds_metrics[seed]["pearson"] for seed in seeds]),
            "source": f"{run_dir.relative_to(ROOT)}/seed*/test_metrics.json",
            "protocol_complete": seeds == list(PROTOCOL_SEEDS),
        })

    data_manifest = SOTA_ROOT / "data-manifest.json"
    receipt = {
        "schema_version": 1,
        "protocol_id": "sota-homology-cold-comparison-v1",
        "generated_on": "2026-09-16",
        "claim_boundary": "Internal homology-cold comparison; Mode B retraining uses published "
                          "hyperparameters and the frozen train partition only. No Mode A numbers "
                          "are included (see overlap-audit.json: 100% sequence overlap with the "
                          "published training corpus makes Mode A in-corpus recall).",
        "split_identity": json.loads(data_manifest.read_text(encoding="utf-8"))["corpus_identity"],
        "data_manifest_sha256": sha256(data_manifest),
        "overlap_audit_sha256": sha256(SOTA_ROOT / "overlap-audit.json"),
        "rows": rows,
    }
    args.out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        print(f"{row['model']}: RMSE {row['rmse']['mean']:.4f}±{row['rmse']['sd']:.4f} "
              f"MAE {row['mae']['mean']:.4f}±{row['mae']['sd']:.4f} "
              f"Pearson {row['pearson']['mean']:.4f}±{row['pearson']['sd']:.4f} "
              f"seeds {row['seeds']}")


if __name__ == "__main__":
    main()
