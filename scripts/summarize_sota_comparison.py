"""Aggregate Mode B SOTA comparison results into one summary receipt.

Reads per-predictor summary.json files produced by the Mode B drivers, recomputes the
BioCandidateRanker ESM-2 reference numbers from its checkpoint metrics, and writes a
manuscript-ready comparison table with the frozen split identity. Governed by
configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOTA_ROOT = ROOT / "artifacts" / "external" / "sota-homology-cold"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_sd(values: list[float]) -> dict:
    return {
        "mean": statistics.mean(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=SOTA_ROOT / "summary.json")
    args = parser.parse_args()

    rows: list[dict] = []

    esm2_metrics = []
    for seed in (7, 42, 123):
        payload = json.loads(
            (ROOT / f"artifacts/esm2-t6-opt-seed{seed}/test_metrics.json").read_text(encoding="utf-8")
        )
        esm2_metrics.append(payload["model_metrics"]["log10_kcat"])
    rows.append({
        "model": "BioCandidateRanker ESM-2 (ours, reference)",
        "kind": "internal",
        "seeds": [7, 42, 123],
        "rmse": mean_sd([m["rmse"] for m in esm2_metrics]),
        "mae": mean_sd([m["mae"] for m in esm2_metrics]),
        "pearson": mean_sd([m["pearson"] for m in esm2_metrics]),
        "source": "artifacts/esm2-t6-opt-seed{7,42,123}/test_metrics.json",
    })

    for predictor, path in (
        ("UniKP Mode B (retrained)", SOTA_ROOT / "unikp-mode-b" / "summary.json"),
        ("DLKcat Mode B (retrained)", SOTA_ROOT / "dlkcat-mode-b" / "summary.json"),
    ):
        if not path.is_file():
            print(f"missing, skipped: {path}")
            continue
        summary = json.loads(path.read_text(encoding="utf-8"))
        seeds = sorted(int(seed) for seed in summary["seeds"])
        rows.append({
            "model": predictor,
            "kind": "sota-mode-b",
            "seeds": seeds,
            "rmse": mean_sd([summary["seeds"][str(seed)]["rmse"] for seed in seeds]),
            "mae": mean_sd([summary["seeds"][str(seed)]["mae"] for seed in seeds]),
            "pearson": mean_sd([summary["seeds"][str(seed)]["pearson"] for seed in seeds]),
            "source": str(path.relative_to(ROOT)),
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
              f"Pearson {row['pearson']['mean']:.4f}±{row['pearson']['sd']:.4f}")


if __name__ == "__main__":
    main()
