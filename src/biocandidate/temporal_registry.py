"""Fail-closed registry for the prospective temporal curation pool.

Run with ``python -m biocandidate.temporal_registry``. Family counts are omitted by
default because source-local MMseqs clusters are not necessarily globally distinct.
The opt-in ``source-local-additive`` policy sums source-local clusters, but refuses to
report a count when a multicomponent source does not declare catalytic (rather than
component) families.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any


ACCEPTED_STATUS = "accepted_homology_cold_pool"
PENDING_PREFIX = "pending"
FAMILY_POLICIES = ("none", "source-local-additive")
TERMINAL_ZERO_DISPOSITION_WORDS = ("block", "complete", "exclude", "not_run")


def finalized_accepted_sources(root: Path, protocol_path: Path) -> list[dict[str, Any]]:
    """Return accepted rows only from sources passing the registry's fail-closed checks."""
    report = build_registry(root, protocol_path)
    finalized = {
        source["source_id"]
        for source in report["sources"]
        if source["accepted_evidence_final"] and source["csv_accepted_records"]
    }
    sources = []
    for source_id in sorted(finalized):
        source_dir = root / source_id
        csv_path = source_dir / "candidate_records.csv"
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            rows = [
                row
                for row in csv.DictReader(handle)
                if (row.get("status_at_normalization") or "").strip() == ACCEPTED_STATUS
            ]
        sources.append(
            {
                "source_id": source_id,
                "source_dir": source_dir,
                "candidate_csv": csv_path,
                "audit_path": source_dir / "homology-audit.json",
                "rows": rows,
            }
        )
    return sources


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _candidate_count(value: Any) -> int | None:
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, dict) and isinstance(value.get("count"), int):
        return value["count"]
    return None


def _audit_candidate_count(audit: dict[str, Any]) -> int | None:
    count = _candidate_count(audit.get("candidate_records"))
    if count is not None:
        return count
    for key in ("saturation_qualified_candidate_records", "candidate_record_count"):
        value = audit.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            return value
    return None


def _accepted_count(audit: dict[str, Any]) -> int | None:
    for key in (
        "accepted_records",
        "accepted_records_after_family_cap",
        "accepted_records_after_homology",
    ):
        value = audit.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            return value
    return None


def _is_finalized_zero_candidate_source(audit: dict[str, Any] | None) -> bool:
    if audit is None:
        return False
    if _accepted_count(audit) == 0:
        return True
    status = audit.get("status")
    return (
        _audit_candidate_count(audit) == 0
        and isinstance(status, str)
        and any(word in status.lower() for word in TERMINAL_ZERO_DISPOSITION_WORDS)
    ) or (
        audit.get("homology_cold_claimed") is False
        and isinstance(status, str)
        and any(word in status.lower() for word in TERMINAL_ZERO_DISPOSITION_WORDS)
    )


def _target_sha256(audit: dict[str, Any]) -> str | None:
    value = audit.get("development_target_sha256")
    if isinstance(value, str):
        return value
    target = audit.get("development_target")
    if isinstance(target, dict) and isinstance(target.get("sha256"), str):
        return target["sha256"]
    return None


def _has_hash(audit: dict[str, Any], *names: str) -> bool:
    return any(isinstance(audit.get(name), str) and len(audit[name]) == 64 for name in names)


def _has_query_hash(audit: dict[str, Any]) -> bool:
    if _has_hash(audit, "construct_sequences_sha256", "homology_queries_sha256"):
        return True
    variants = audit.get("variant_sequences")
    return isinstance(variants, dict) and _has_hash(variants, "sha256")


def _substrate_id(row: dict[str, str]) -> str | None:
    # The curated schemas use one of these columns for the varied, measured substrate.
    for name in (
        "substrate_pubchem_cid",
        "variable_substrate_pubchem_cid",
        "substrate_1_pubchem_cid",
    ):
        value = (row.get(name) or "").strip()
        if value:
            return f"pubchem:{value}"
    return None


def _source_family_count(audit: dict[str, Any]) -> tuple[int | None, str | None]:
    if not _has_hash(audit, "family_cluster_sha256"):
        return None, "audit lacks source-local family-cluster SHA256"
    if "complex_homology_rule" in audit or "candidate_component_families" in audit:
        value = audit.get("accepted_catalytic_families")
        if not isinstance(value, int) or isinstance(value, bool):
            return None, "multicomponent audit lacks accepted_catalytic_families"
        return value, None
    for key in ("accepted_catalytic_families", "candidate_mmseqs_families"):
        value = audit.get(key)
        if isinstance(value, int) and not isinstance(value, bool):
            return value, None
    return None, "audit lacks a source-local catalytic family count"


def build_registry(
    root: Path,
    protocol_path: Path,
    *,
    family_policy: str = "none",
) -> dict[str, Any]:
    """Reconcile source CSVs and homology audits into a readiness report."""
    if family_policy not in FAMILY_POLICIES:
        raise ValueError(f"family_policy must be one of {FAMILY_POLICIES}")

    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    leakage = protocol["leakage_gates"]
    gates = protocol["diversity_and_size_gate"]
    registry_policy = protocol.get("registry_readiness", {})
    require_zero_candidate_evidence = registry_policy.get(
        "finalized_zero_candidate_sources_must_pass_evidence_checks", False
    )
    source_dirs = sorted(
        path for path in root.iterdir() if path.is_dir() and (
            (path / "candidate_records.csv").exists() or (path / "homology-audit.json").exists()
        )
    )

    sources: list[dict[str, Any]] = []
    accepted_rows: list[dict[str, str]] = []
    blockers: list[str] = []
    family_total = 0
    family_unavailable: list[str] = []
    source_summary = {
        "accepted_pool_sources": 0,
        "pending_sources": 0,
        "finalized_zero_candidate_sources": 0,
        "unresolved_zero_candidate_sources": 0,
        "nonblocking_source_issues": 0,
    }

    for source_dir in source_dirs:
        source_id = source_dir.name
        csv_path = source_dir / "candidate_records.csv"
        audit_path = source_dir / "homology-audit.json"
        issues: list[str] = []
        rows: list[dict[str, str]] = []
        audit: dict[str, Any] | None = None

        if csv_path.exists():
            with csv_path.open(newline="", encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            if rows and "status_at_normalization" not in rows[0]:
                issues.append("candidate CSV lacks status_at_normalization")
        else:
            issues.append("missing candidate_records.csv")

        if audit_path.exists():
            try:
                loaded = json.loads(audit_path.read_text(encoding="utf-8"))
                if not isinstance(loaded, dict):
                    raise ValueError("top-level value is not an object")
                audit = loaded
            except (json.JSONDecodeError, ValueError) as error:
                issues.append(f"invalid homology-audit.json: {error}")
        else:
            issues.append("missing homology-audit.json")

        statuses: dict[str, int] = {}
        for row in rows:
            status = (row.get("status_at_normalization") or "").strip()
            statuses[status] = statuses.get(status, 0) + 1
        csv_accepted = statuses.get(ACCEPTED_STATUS, 0)
        pending = sum(count for status, count in statuses.items() if status.startswith(PENDING_PREFIX))

        if audit is not None:
            if audit.get("source_id") != source_id:
                issues.append("audit source_id does not match directory")
            if _audit_candidate_count(audit) != len(rows):
                issues.append("audit candidate count does not match CSV")
            audit_accepted = _accepted_count(audit)
            if audit_accepted is None:
                issues.append("audit lacks an accepted-record count")
            elif audit_accepted != csv_accepted:
                issues.append(
                    f"audit accepted count ({audit_accepted}) does not match CSV ({csv_accepted})"
                )
            if pending and audit_accepted:
                issues.append(f"audit claims accepted records while {pending} CSV rows are pending")
            if audit.get("min_identity") != leakage["mmseqs_min_identity"]:
                issues.append("audit MMseqs identity differs from frozen protocol")
            if audit.get("coverage") != leakage["mmseqs_coverage"]:
                issues.append("audit MMseqs coverage differs from frozen protocol")
            if audit.get("coverage_mode") != leakage["mmseqs_coverage_mode"]:
                issues.append("audit MMseqs coverage mode differs from frozen protocol")
            if not audit.get("method") or not audit.get("mmseqs_version"):
                issues.append("audit lacks MMseqs method/version")
            if not _target_sha256(audit):
                issues.append("audit lacks development-target SHA256")
            if not _has_query_hash(audit):
                issues.append("audit lacks candidate/query sequence SHA256")
            if not _has_hash(audit, "homology_hits_sha256"):
                hits = audit.get("homology_hits")
                if not (isinstance(hits, dict) and _has_hash(hits, "sha256")):
                    issues.append("audit lacks homology-hit evidence SHA256")
            candidate_meta = audit.get("candidate_records")
            if (
                csv_path.exists()
                and isinstance(candidate_meta, dict)
                and candidate_meta.get("sha256") != _sha256(csv_path)
            ):
                issues.append("candidate CSV SHA256 differs from audit")

        finalized_zero_candidate = (
            not csv_accepted
            and not pending
            and _is_finalized_zero_candidate_source(audit)
        )
        blocks_accepted_pool_readiness = bool(
            csv_accepted
            or pending
            or not finalized_zero_candidate
            or require_zero_candidate_evidence
        )
        accepted_evidence_final = not issues and not pending
        if csv_accepted:
            source_summary["accepted_pool_sources"] += 1
        if pending:
            source_summary["pending_sources"] += 1
        if finalized_zero_candidate:
            source_summary["finalized_zero_candidate_sources"] += 1
        elif not csv_accepted and not pending:
            source_summary["unresolved_zero_candidate_sources"] += 1

        if csv_accepted and not accepted_evidence_final:
            blockers.append(f"{source_id}: accepted records lack finalized homology evidence")
        if pending:
            blockers.append(f"{source_id}: {pending} candidate records remain pending")
        if issues and blocks_accepted_pool_readiness:
            blockers.extend(f"{source_id}: {issue}" for issue in issues)
        elif issues:
            source_summary["nonblocking_source_issues"] += len(issues)
        if accepted_evidence_final:
            accepted_rows.extend(
                row
                for row in rows
                if (row.get("status_at_normalization") or "").strip() == ACCEPTED_STATUS
            )

        if family_policy == "source-local-additive" and csv_accepted and accepted_evidence_final:
            assert audit is not None
            count, reason = _source_family_count(audit)
            if reason:
                family_unavailable.append(f"{source_id}: {reason}")
            else:
                family_total += count or 0

        sources.append(
            {
                "source_id": source_id,
                "candidate_records": len(rows),
                "csv_accepted_records": csv_accepted,
                "pending_records": pending,
                "accepted_evidence_final": accepted_evidence_final,
                "disposition": (
                    "accepted_pool_source"
                    if csv_accepted
                    else "pending_source"
                    if pending
                    else "finalized_zero_candidate_source"
                    if finalized_zero_candidate
                    else "unresolved_zero_candidate_source"
                ),
                "blocks_accepted_pool_readiness": blocks_accepted_pool_readiness,
                "issues": issues,
            }
        )

    substrate_ids = {_substrate_id(row) for row in accepted_rows}
    if None in substrate_ids:
        missing = sum(_substrate_id(row) is None for row in accepted_rows)
        blockers.append(f"{missing} accepted records lack a supported substrate structure identifier")
        substrate_ids.discard(None)

    family_count: int | None = None
    family_note = "not reported; select an explicit family policy"
    if family_policy == "source-local-additive":
        if family_unavailable:
            family_note = "; ".join(family_unavailable)
            blockers.append(f"family count unavailable: {family_note}")
        else:
            family_count = family_total
            family_note = (
                "sum of source-local MMseqs catalytic-family counts; assumes families from "
                "different sources are distinct"
            )

    counts = {
        "accepted_records": len(accepted_rows),
        "globally_unique_substrates": len(substrate_ids),
        "families": family_count,
    }
    if counts["accepted_records"] < gates["minimum_records"]:
        blockers.append(
            f"accepted records {counts['accepted_records']} < required {gates['minimum_records']}"
        )
    if counts["globally_unique_substrates"] < gates["minimum_unique_substrates"]:
        blockers.append(
            "globally unique substrates "
            f"{counts['globally_unique_substrates']} < required {gates['minimum_unique_substrates']}"
        )
    if family_count is None:
        blockers.append("global family readiness gate cannot be evaluated")
    elif family_count < gates["minimum_mmseqs_families"]:
        blockers.append(f"families {family_count} < required {gates['minimum_mmseqs_families']}")

    return {
        "schema_version": 2,
        "generated_on": date.today().isoformat(),
        "root": root.as_posix(),
        "protocol": protocol_path.as_posix(),
        "protocol_sha256": _sha256(protocol_path),
        "accepted_status": ACCEPTED_STATUS,
        "substrate_identity_policy": "PubChem CID of the varied, measured substrate",
        "family_policy": family_policy,
        "family_policy_note": family_note,
        "registry_readiness_policy": {
            "finalized_zero_candidate_sources_must_pass_evidence_checks": (
                require_zero_candidate_evidence
            )
        },
        "source_summary": source_summary,
        "counts": counts,
        "readiness_gate_passes": not blockers,
        "blockers": blockers,
        "sources": sources,
        "claim_boundary": "Curation pool only; no model predictions were generated.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("artifacts/external/temporal-absolute-kinetics"),
    )
    parser.add_argument(
        "--protocol", type=Path, default=Path("configs/temporal_absolute_kinetics_protocol.json")
    )
    parser.add_argument("--family-policy", choices=FAMILY_POLICIES, default="none")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = build_registry(args.root, args.protocol, family_policy=args.family_policy)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="ascii")
    else:
        print(rendered, end="")
    return 0 if report["readiness_gate_passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
