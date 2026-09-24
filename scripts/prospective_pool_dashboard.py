"""Summarize the frozen prospective-pool readiness gate and the options around it.

Reads the live readiness audit and the frozen release thresholds, then reports
per-gate shortfalls, per-family saturation, and what each resolution path would
require. This is a decision aid: it never scores records or generates predictions.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

DEFAULT_AUDIT = Path("artifacts/external/temporal-global-family-audit/readiness-audit.json")
DEFAULT_RELEASE_MANIFEST = Path("configs/software_implementation_release_v6.json")

GATES = (
    (
        "records_after_global_family_cap",
        "records_after_family_cap",
        "required_records",
        "capped records",
    ),
    (
        "global_families",
        "global_mmseqs_families",
        "required_global_mmseqs_families",
        "global MMseqs families",
    ),
    (
        "substrates_after_global_family_cap",
        "unique_substrates",
        "required_unique_substrates",
        "unique substrates",
    ),
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def family_table(audit: dict[str, Any], cap: int) -> list[dict[str, Any]]:
    counts = audit.get("family_record_counts_before_cap", {})
    rows = []
    for family, records in sorted(counts.items()):
        rows.append(
            {
                "family": family,
                "records": records,
                "at_cap": records >= cap,
                "headroom": max(cap - records, 0),
            }
        )
    return rows


def evaluate_gates(audit: dict[str, Any], thresholds: dict[str, Any]) -> list[dict[str, Any]]:
    """Compare live audit counts against frozen thresholds.

    The audit is the measurement of record; the release manifest is only a
    fallback for counts the audit does not carry.
    """
    counts = audit.get("counts", {})
    gates = []
    for audit_key, manifest_key, required_key, label in GATES:
        observed = counts.get(audit_key)
        source = "audit"
        if observed is None:
            observed = thresholds.get(manifest_key)
            source = "release_manifest"
        required = thresholds.get(required_key)
        if observed is None or required is None:
            continue
        gates.append(
            {
                "gate": label,
                "observed": observed,
                "required": required,
                "shortfall": max(required - observed, 0),
                "satisfied": observed >= required,
                "observed_source": source,
            }
        )
    return gates


def source_projection(audit: dict[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Estimate how many additional accepted sources the unmet gates imply."""
    counts = audit.get("counts", {})
    sources = counts.get("accepted_sources") or 0
    records = counts.get("records_after_global_family_cap") or 0
    per_source = records / sources if sources else None

    record_gate = next((gate for gate in gates if gate["gate"] == "capped records"), None)
    family_gate = next((gate for gate in gates if gate["gate"] == "global MMseqs families"), None)
    records_needed = record_gate["shortfall"] if record_gate else 0
    families_needed = family_gate["shortfall"] if family_gate else 0

    return {
        "accepted_sources": sources,
        "records_per_accepted_source": round(per_source, 2) if per_source else None,
        "records_needed": records_needed,
        "families_needed": families_needed,
        "sources_if_records_drive_acquisition": (
            round(records_needed / per_source, 1) if per_source else None
        ),
        "note": (
            "New families can only come from new sources, so the binding requirement is "
            "whichever of the two is larger."
        ),
    }


def build_options(
    gates: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    projection: dict[str, Any],
) -> list[dict[str, Any]]:
    unmet = [gate for gate in gates if not gate["satisfied"]]
    total_headroom = sum(row["headroom"] for row in rows)
    at_cap = [row["family"] for row in rows if row["at_cap"]]

    return [
        {
            "id": "acquire",
            "label": "Continue acquisition until the frozen gate passes",
            "unmet_gates": [gate["gate"] for gate in unmet],
            "requires": (
                f"{projection['records_needed']} additional accepted records and "
                f"{projection['families_needed']} additional global families; "
                f"roughly {projection['sources_if_records_drive_acquisition']} more accepted "
                "sources at the observed yield"
            ),
            "cost": "data-acquisition effort with an open-ended calendar; no compute needed",
            "keeps_frozen_evidence_intact": True,
        },
        {
            "id": "relax",
            "label": "Amend the gate under a new dated protocol and disclose it",
            "unmet_gates": [gate["gate"] for gate in unmet],
            "requires": (
                "a superseding protocol that states the relaxed thresholds and the reason, "
                "without editing the frozen temporal_absolute_kinetics_protocol.json"
            ),
            "cost": "editorial; reviewers may treat the pool as underpowered",
            "keeps_frozen_evidence_intact": True,
        },
        {
            "id": "reframe",
            "label": "Reframe the contribution onto endpoint-transfer analysis",
            "unmet_gates": [gate["gate"] for gate in unmet],
            "requires": (
                "make H1/H2/H3 the primary claims and report the prospective pool as an "
                "incomplete, unscored resource"
            ),
            "cost": "larger manuscript rewrite; independent validation stays absent",
            "keeps_frozen_evidence_intact": True,
        },
    ], {
        "total_family_headroom_before_cap": total_headroom,
        "families_at_cap": at_cap,
        "headroom_note": (
            "Headroom is theoretical: every additional record must still pass every "
            "inclusion gate, and records are capped per family."
        ),
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = ["# Prospective pool readiness dashboard", ""]
    lines.append(f"Source audit: `{report['audit_path']}`")
    lines.append(f"Threshold source: `{report['release_manifest_path']}`")
    lines.append("")
    lines.append("## Gate status")
    lines.append("")
    lines.append("| Gate | Observed | Required | Shortfall | Status |")
    lines.append("|---|---:|---:|---:|---|")
    for gate in report["gates"]:
        status = "pass" if gate["satisfied"] else "FAIL"
        lines.append(
            f"| {gate['gate']} | {gate['observed']} | {gate['required']} | "
            f"{gate['shortfall']} | {status} |"
        )
    lines.append("")
    lines.append(f"Frozen gate passes: **{report['readiness_gate_passes']}**")
    lines.append(f"Predictions permitted: **{report['predictions_permitted']}**")
    lines.append("")

    lines.append("## Family saturation")
    lines.append("")
    lines.append(
        f"Cap per family: {report['family_cap']}. Total theoretical headroom: "
        f"{report['saturation']['total_family_headroom_before_cap']}."
    )
    lines.append("")
    lines.append("| Family | Records | At cap | Headroom |")
    lines.append("|---|---:|---|---:|")
    for row in report["families"]:
        lines.append(
            f"| {row['family']} | {row['records']} | {'yes' if row['at_cap'] else 'no'} | "
            f"{row['headroom']} |"
        )
    lines.append("")
    lines.append(report["saturation"]["headroom_note"])
    lines.append("")

    lines.append("## Acquisition projection")
    lines.append("")
    for key, value in report["projection"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")

    lines.append("## Options")
    lines.append("")
    for option in report["options"]:
        lines.append(f"### {option['label']}")
        lines.append("")
        lines.append(f"- id: `{option['id']}`")
        lines.append(f"- unresolved gates: {', '.join(option['unmet_gates']) or 'none'}")
        lines.append(f"- requires: {option['requires']}")
        lines.append(f"- cost: {option['cost']}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--release-manifest", type=Path, default=DEFAULT_RELEASE_MANIFEST)
    parser.add_argument("--markdown", type=Path, default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    audit_path = args.audit if args.audit.is_absolute() else root / args.audit
    manifest_path = args.release_manifest if args.release_manifest.is_absolute() else root / args.release_manifest

    if not audit_path.is_file():
        print(json.dumps({"valid": False, "issues": [f"missing audit: {audit_path}"]}, indent=2))
        return 1
    if not manifest_path.is_file():
        print(json.dumps({"valid": False, "issues": [f"missing manifest: {manifest_path}"]}, indent=2))
        return 1

    audit = load_json(audit_path)
    release = load_json(manifest_path)
    thresholds = release.get("temporal_pool", {})
    cap = audit.get("parameters", {}).get("maximum_records_per_family", 20)

    gates = evaluate_gates(audit, thresholds)
    families = family_table(audit, cap)
    projection = source_projection(audit, gates)
    options, saturation = build_options(gates, families, projection)

    report = {
        "schema_version": 1,
        "audit_path": args.audit.as_posix(),
        "release_manifest_path": args.release_manifest.as_posix(),
        "generated_from_audit_date": audit.get("generated_on"),
        "readiness_gate_passes": audit.get("readiness_gate_passes"),
        "predictions_permitted": thresholds.get("predictions_permitted"),
        "family_cap": cap,
        "gates": gates,
        "families": families,
        "saturation": saturation,
        "projection": projection,
        "options": options,
        "claim_boundary": (
            "Decision aid only. No records were scored, selected by prediction, or used to "
            "train anything."
        ),
    }

    if args.markdown:
        markdown_path = args.markdown if args.markdown.is_absolute() else root / args.markdown
        markdown_path.parent.mkdir(parents=True, exist_ok=True)
        markdown_path.write_text(render_markdown(report), encoding="utf-8")
        report["markdown"] = args.markdown.as_posix()

    if args.json_out:
        json_path = args.json_out if args.json_out.is_absolute() else root / args.json_out
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
