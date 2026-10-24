"""Freeze the F1a prospective prediction before the assay is run.

The royal jelly acid programme's design matrix defines F1a as a direct
enzymology comparison of HFD1, ALD2, ALD3 and ALD6 on 10-oxodecanoate. That
makes it the only planned experiment that satisfies all three requirements for
a prospective test of this project's central question: one substrate, one assay
condition, several enzyme candidates, and predictions that already exist with
recorded input hashes.

This script does not predict anything new. It reads the stored task-2
predictions, records their identity, and writes a manifest that fixes the
predicted ordering now, so a later comparison cannot be adjusted after the fact.

Output: configs/f1a_prospective_freeze.json

Reads the stored task-2 prediction JSON, its input spec, and candidate sequence
files under ``artifacts/10h2da-yeast-prediction``; writes the freeze manifest
with SHA-256 identity for every input, the pre-registered decision rule, and a
discriminating-power check. Governance role: pre-registration only -- nothing
is predicted here and the manifest exists so the comparison cannot be adjusted
after measurements are seen.

中文：在测定运行前冻结 F1a 前瞻预测（HFD1/ALD2/ALD3/ALD6 对 10-氧癸酸的直接
酶学比较），使预测顺序事后不可调整。
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

PREDICTIONS = ROOT / "artifacts/10h2da-yeast-prediction/task2_10hda_terminal_oxidation-predictions.json"
INPUT_SPEC = ROOT / "artifacts/10h2da-yeast-prediction/task2_10hda_terminal_oxidation.json"
SEQUENCE_DIR = ROOT / "artifacts/10h2da-yeast-prediction/sequences"
OUTPUT = ROOT / "configs/f1a_prospective_freeze.json"

# The four enzymes F1a names, by UniProt accession and the candidate id used in
# the stored prediction, so the manifest resolves to exact prediction entries.
F1A_ENZYMES = (
    {"uniprot": "Q04458", "gene": "HFD1", "role": "long-chain fatty aldehyde dehydrogenase (category evidence)"},
    {"uniprot": "P47771", "gene": "ALD2", "role": "broad aldehyde dehydrogenase"},
    {"uniprot": "P54114", "gene": "ALD3", "role": "broad aldehyde dehydrogenase"},
    {"uniprot": "P54115", "gene": "ALD6", "role": "broad aldehyde dehydrogenase"},
)

SUBSTRATE = {
    "name": "10-oxodecanoate",
    "note": "the F1a substrate, i.e. the second oxidation step of the 10-HDA route",
}

ASSAY = {
    "design": "purified enzyme or competent lysate, kcat and Km on 10-oxodecanoate",
    "detection": "NADH coupled at 340 nm; the coupling enzyme must be shown not to be rate limiting",
    "condition_matching": (
        "pH, temperature, buffer and enzyme concentration recorded per enzyme so the "
        "four are compared under matched conditions"
    ),
}

CLAIM_BOUNDARY = (
    "Prospective prediction freeze only. The predictions were produced before this "
    "manifest and are uncalibrated; the assay has not been performed. Agreement or "
    "disagreement will be reported as recorded, without selecting the comparison "
    "rule after seeing the measurements."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_freeze() -> dict[str, Any]:
    """Assemble the freeze manifest with predictions, ranks, and decision rule.
中文：组装带预测、排名与决策规则的冻结清单。
"""
    if not PREDICTIONS.is_file():
        raise SystemExit(f"missing stored predictions: {PREDICTIONS}")

    payload = json.loads(PREDICTIONS.read_text(encoding="utf-8"))
    by_candidate = {
        entry["candidate_id"]: entry for entry in payload.get("predictions", [])
    }

    frozen: list[dict[str, Any]] = []
    for enzyme in F1A_ENZYMES:
        matches = [
            candidate_id
            for candidate_id in by_candidate
            if candidate_id.startswith(f"{enzyme['gene']}_")
        ]
        if not matches:
            raise SystemExit(
                f"no stored prediction for {enzyme['gene']} ({enzyme['uniprot']}); "
                f"available: {sorted(by_candidate)[:8]}"
            )
        # one candidate per gene in this task's candidate set
        candidate_id = sorted(matches)[0]
        task = by_candidate[candidate_id]["tasks"]["log10_kcat"]
        if task.get("status") != "trained":
            raise SystemExit(f"{candidate_id} carries no trained kcat task")
        frozen.append(
            {
                "uniprot": enzyme["uniprot"],
                "gene": enzyme["gene"],
                "category_role": enzyme["role"],
                "candidate_id": candidate_id,
                "predicted_log10_kcat": task["mean"],
                "predicted_kcat_s-1": 10 ** task["mean"],
                "predicted_standard_deviation_log10": task["standard_deviation"],
            }
        )

    # The predicted ordering is what F1a tests. Sorting descending by prediction
    # and recording the rank here is the whole point of freezing early.
    for rank, entry in enumerate(
        sorted(frozen, key=lambda item: -item["predicted_log10_kcat"]), start=1
    ):
        entry["predicted_rank"] = rank

    estimates = [entry["predicted_log10_kcat"] for entry in frozen]
    spread = max(estimates) - min(estimates)
    deviations = [entry["predicted_standard_deviation_log10"] for entry in frozen]

    return {
        "schema_version": 1,
        "purpose": "F1a prospective prediction freeze for the royal jelly acid programme",
        "frozen_on": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "frozen_before_assay",
        "substrate": SUBSTRATE,
        "assay_plan": ASSAY,
        "predictions": sorted(frozen, key=lambda item: item["predicted_rank"]),
        "prediction_identity": {
            "predictions_path": PREDICTIONS.relative_to(ROOT).as_posix(),
            "predictions_sha256": sha256(PREDICTIONS),
            "input_spec_path": INPUT_SPEC.relative_to(ROOT).as_posix()
            if INPUT_SPEC.is_file()
            else None,
            "input_spec_sha256": sha256(INPUT_SPEC) if INPUT_SPEC.is_file() else None,
            "sequence_files": sorted(
                {
                    path.name: sha256(path)
                    for path in SEQUENCE_DIR.glob("*")
                    if path.is_file()
                }.items()
            )
            if SEQUENCE_DIR.is_dir()
            else [],
        },
        "decision_rule": {
            "primary": "Spearman rank correlation between predicted and measured kcat across the four enzymes",
            "secondary": "whether the measured best enzyme matches the predicted rank-1 enzyme",
            "pre_registered": True,
            "note": (
                "Fixed now, before any measurement. If the primary statistic is "
                "computed differently later, that must be reported as a deviation."
            ),
        },
        "discriminating_power_check": {
            "predicted_spread_log10": spread,
            "max_standard_deviation_log10": max(deviations),
            "spread_exceeds_uncertainty": spread > max(deviations),
            "note": (
                "Recorded before the assay because it bounds what the test can show. "
                "If the predicted spread is smaller than the model's own standard "
                "deviation, the ordering is not statistically distinguishable and a "
                "measured disagreement would not refute the model."
            ),
        },
        "claim_boundary": CLAIM_BOUNDARY,
    }


def main() -> int:
    """Write the freeze manifest and print the recorded ordering summary.
中文：写冻结清单并打印记录的排序摘要。
"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    freeze = build_freeze()
    target = args.output if args.output.is_absolute() else ROOT / args.output
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")

    check = freeze["discriminating_power_check"]
    print(json.dumps({
        "written": target.relative_to(ROOT).as_posix(),
        "frozen_enzymes": [entry["gene"] for entry in freeze["predictions"]],
        "predicted_order": [
            f"{entry['gene']} ({entry['predicted_kcat_s-1']:.2f})"
            for entry in freeze["predictions"]
        ],
        "spread_log10": round(check["predicted_spread_log10"], 3),
        "max_sd_log10": round(check["max_standard_deviation_log10"], 3),
        "spread_exceeds_uncertainty": check["spread_exceeds_uncertainty"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
