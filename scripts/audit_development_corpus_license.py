"""Audit the declared license and redistribution state of every project corpus.

Reads configs/data_license_inventory.json, checks that every declared evidence
path exists, cross-checks the inventory against the manifests that independently
record license state, and reports which corpora block a public release.

The script asserts nothing about license terms on its own; it verifies that the
declarations are complete, internally consistent, and still supported.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

DEFAULT_INVENTORY = Path("configs/data_license_inventory.json")
DEFAULT_RELEASE_MANIFEST = Path("configs/software_implementation_release_v6.json")

# Declarations this project records elsewhere; the inventory must agree with them.
CROSS_CHECKS = (
    {
        "corpus_id": "unikp-kcat-development-corpus",
        "path": Path("configs/unikp_source_manifest.json"),
        "expected_key": "license_status",
        "expected_value": "redistribution_not_verified",
        "maps_to_field": "redistribution",
        "explanation": "unikp_source_manifest.json records redistribution_not_verified",
    },
)

REQUIRED_FIELDS = ("id", "role", "license_name", "license_basis", "redistribution", "evidence")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(
    root: Path,
    inventory: dict[str, Any],
    release_manifest: dict[str, Any] | None,
    cross_checks: tuple[dict[str, Any], ...] | None = None,
) -> dict[str, Any]:
    checks = CROSS_CHECKS if cross_checks is None else cross_checks
    issues: list[str] = []
    allowed = set(inventory.get("policy", {}).get("allowed_redistribution_states", []))

    corpora: list[dict[str, Any]] = []
    seen_ids: set[str] = set()

    for entry in inventory.get("corpora", []):
        corpus_id = entry.get("id", "<unnamed>")
        entry_issues: list[str] = []

        for field in REQUIRED_FIELDS:
            if not entry.get(field):
                entry_issues.append(f"missing required field: {field}")

        if corpus_id in seen_ids:
            entry_issues.append(f"duplicate corpus id: {corpus_id}")
        seen_ids.add(corpus_id)

        redistribution = entry.get("redistribution")
        if allowed and redistribution not in allowed:
            entry_issues.append(
                f"redistribution {redistribution!r} is not in the allowed set {sorted(allowed)}"
            )

        if not entry.get("evidence"):
            entry_issues.append("no evidence path recorded, so the claim cannot be checked")

        missing_evidence = [
            path for path in entry.get("evidence", []) if not (root / path).exists()
        ]
        for path in missing_evidence:
            entry_issues.append(f"missing evidence path: {path}")

        shipping = entry.get("redistribute_in_release", False)
        blocked = redistribution in {"not_verified", "blocked_pending_permission"}
        if shipping and blocked:
            entry_issues.append(
                "declared shippable in a release while redistribution is unresolved"
            )
        if blocked and not entry.get("blocking_reason"):
            entry_issues.append("blocked corpus has no blocking_reason")

        corpora.append(
            {
                "id": corpus_id,
                "role": entry.get("role"),
                "license_name": entry.get("license_name"),
                "redistribution": redistribution,
                "redistribute_in_release": shipping,
                "blocking_reason": entry.get("blocking_reason"),
                "evidence_count": len(entry.get("evidence", [])),
                "missing_evidence": missing_evidence,
                "issues": entry_issues,
            }
        )
        issues.extend(f"{corpus_id}: {issue}" for issue in entry_issues)

    for check in checks:
        path = root / check["path"]
        if not path.is_file():
            issues.append(f"cross-check source missing: {check['path'].as_posix()}")
            continue
        recorded = load_json(path).get(check["expected_key"])
        if recorded != check["expected_value"]:
            issues.append(
                f"{check['path'].as_posix()} records {check['expected_key']}={recorded!r}, "
                f"but the inventory assumes {check['expected_value']!r}"
            )
        else:
            match = next(
                (item for item in corpora if item["id"] == check["corpus_id"]), None
            )
            if match and match["redistribution"] != "not_verified":
                issues.append(
                    f"inventory marks {check['corpus_id']} as {match['redistribution']!r} "
                    f"while {check['path'].as_posix()} says {check['expected_value']!r}"
                )

    shippable = [item["id"] for item in corpora if item["redistribute_in_release"]]
    blocking = [item["id"] for item in corpora if item["redistribution"] in
                {"not_verified", "blocked_pending_permission"}]

    release_blockers: list[str] = []
    if "unikp-kcat-development-corpus" in blocking:
        release_blockers.append(
            "the training corpus cannot ship, so a public archive must be code-and-protocols only"
        )
    if release_manifest is not None:
        release_blockers.append(
            "release manifest already excludes the training corpus; the archive ships "
            f"{len(release_manifest.get('frozen_files', []))} frozen evidence files only"
        )

    return {
        "schema_version": 1,
        "inventory_path": DEFAULT_INVENTORY.as_posix(),
        "corpus_count": len(corpora),
        "valid": not issues,
        "issues": issues,
        "corpora": corpora,
        "shippable_corpora": shippable,
        "blocking_corpora": blocking,
        "release_blockers": release_blockers,
        "claim_boundary": (
            "Declaration audit only. It confirms the project's recorded license basis is "
            "complete and consistent; it does not grant or verify any third-party right."
        ),
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = ["# Development corpus license audit", ""]
    lines.append(f"Inventory: `{report['inventory_path']}`")
    lines.append("")
    lines.append("| Corpus | Role | License | Redistribution | In release |")
    lines.append("|---|---|---|---|---|")
    for item in report["corpora"]:
        lines.append(
            f"| {item['id']} | {item['role']} | {item['license_name']} | "
            f"{item['redistribution']} | {'yes' if item['redistribute_in_release'] else 'no'} |"
        )
    lines.append("")
    lines.append(f"Consistent: **{report['valid']}**")
    lines.append("")
    if report["blocking_corpora"]:
        lines.append("## Corpora blocking a public release")
        lines.append("")
        for item in report["corpora"]:
            if item["id"] in report["blocking_corpora"]:
                lines.append(f"- `{item['id']}`: {item['blocking_reason']}")
        lines.append("")
    if report["issues"]:
        lines.append("## Issues")
        lines.append("")
        for issue in report["issues"]:
            lines.append(f"- {issue}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--release-manifest", type=Path, default=DEFAULT_RELEASE_MANIFEST)
    parser.add_argument("--markdown", type=Path, default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit non-zero when any declaration is incomplete or inconsistent",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    inventory_path = args.inventory if args.inventory.is_absolute() else root / args.inventory
    manifest_path = args.release_manifest if args.release_manifest.is_absolute() else root / args.release_manifest

    if not inventory_path.is_file():
        print(json.dumps({"valid": False, "issues": [f"missing inventory: {args.inventory}"]}, indent=2))
        return 1

    release_manifest = load_json(manifest_path) if manifest_path.is_file() else None
    report = audit(root, load_json(inventory_path), release_manifest)

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
    if args.strict and not report["valid"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
