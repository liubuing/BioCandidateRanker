"""Verify and summarize the frozen ESM-2 modality diagnostic results."""

from __future__ import annotations

import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = json.loads(
    (ROOT / "configs" / "esm2_modality_diagnostic_protocol.json").read_text(encoding="utf-8")
)
OUTPUT = ROOT / "artifacts" / "esm2-modality-diagnostic"
METRICS = ("rmse", "mae", "pearson")
VARIANTS = ("protein-only", "protein-molecule", "full")


def metric_path(variant: str, seed: int) -> Path:
    if variant == "full":
        return ROOT / "artifacts" / f"esm2-t6-opt-seed{seed}" / "test_metrics.json"
    return OUTPUT / f"{variant}-seed{seed}" / "test_metrics.json"


def main() -> None:
    values: dict[str, dict[int, dict[str, float]]] = {}
    for variant in VARIANTS:
        values[variant] = {}
        for seed in PROTOCOL["seeds"]:
            path = metric_path(variant, seed)
            payload = json.loads(path.read_text(encoding="utf-8"))
            for part in ("source", "split"):
                observed = payload[f"{part}_identity"]
                expected = PROTOCOL[f"{part}_identity"]
                for key in ("sha256", "size_bytes", "row_count"):
                    if observed[key] != expected[key]:
                        raise ValueError(f"{path}: {part} {key} mismatch")
            if payload["partition"] != "test" or payload["record_count"] != 1646:
                raise ValueError(f"{path}: unexpected partition or record count")
            values[variant][seed] = {
                name: payload["model_metrics"]["log10_kcat"][name] for name in METRICS
            }

    result = {
        "protocol": "configs/esm2_modality_diagnostic_protocol.json",
        "claim_boundary": PROTOCOL["claim_boundary"],
        "seeds": PROTOCOL["seeds"],
        "variants": {},
        "paired_differences": {},
    }
    for variant in VARIANTS:
        rows = values[variant]
        result["variants"][variant] = {
            "per_seed": rows,
            "summary": {
                name: {
                    "mean": statistics.mean(rows[seed][name] for seed in PROTOCOL["seeds"]),
                    "sd": statistics.stdev(rows[seed][name] for seed in PROTOCOL["seeds"]),
                }
                for name in METRICS
            },
        }
    for left, right in (("protein-only", "protein-molecule"),
                        ("protein-molecule", "full"), ("protein-only", "full")):
        result["paired_differences"][f"{right}_minus_{left}"] = {}
        for name in METRICS:
            differences = {
                seed: values[right][seed][name] - values[left][seed][name]
                for seed in PROTOCOL["seeds"]
            }
            result["paired_differences"][f"{right}_minus_{left}"][name] = {
                "per_seed": differences,
                "mean": statistics.mean(differences.values()),
                "sd": statistics.stdev(differences.values()),
            }
    target = OUTPUT / "summary.json"
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(target)


if __name__ == "__main__":
    main()
