"""Verify and summarize the frozen ESM-2 modality diagnostic results.

Reads each variant/seed ``test_metrics.json`` (protein-only and
protein-molecule variants under ``artifacts/esm2-modality-diagnostic``; full
models under their per-seed artifact directories), verifies source/split
identities against the diagnostic protocol plus the 1,646-row test partition,
and writes ``summary.json`` with per-variant means, SDs, and paired
variant-to-variant differences. Governance role: fail-closed verification and
aggregation only -- any identity or partition mismatch aborts before the
summary is written, and no model is re-run.

中文：校验并汇总冻结的 ESM-2 模态诊断结果：逐变体/种子读 test_metrics.json，
核对身份后聚合指标并写模态汇总。
"""

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
    """Locate the frozen test-metrics file for one variant/seed pair.
中文：定位一个变体/种子对的冻结测试指标文件。
"""
    if variant == "full":
        return ROOT / "artifacts" / f"esm2-t6-opt-seed{seed}" / "test_metrics.json"
    return OUTPUT / f"{variant}-seed{seed}" / "test_metrics.json"


def main() -> None:
    """Verify identities, aggregate metrics, and write the modality summary.
中文：核对身份、聚合指标并写模态汇总。
"""
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
