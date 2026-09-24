"""Draft the manuscript data availability statement from the license audit.

Consumes the JSON produced by scripts/audit_development_corpus_license.py and
emits a statement that separates what ships in the release, what a reader must
obtain from the original deposit, and what remains unresolved.

The generated text is a draft for author review. It never asserts a right the
audit did not record.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SHIPPABLE_STATES = {
    "permitted_with_attribution": "redistributed in the release with attribution",
    "permitted_unconditionally": "redistributed in the release without further condition",
}
OBTAIN_STATES = {
    "not_verified": "not redistributed; obtainable only from the original source",
    "blocked_pending_permission": "not redistributed; terms unresolved",
    "not_applicable_reference_only": "not redistributed; referenced as an external tool",
}


def group_corpora(report: dict[str, Any]) -> tuple[list[dict], list[dict], list[dict]]:
    ships, obtain, unresolved = [], [], []
    for corpus in report.get("corpora", []):
        state = corpus.get("redistribution")
        if state in SHIPPABLE_STATES and corpus.get("redistribute_in_release"):
            ships.append(corpus)
        elif state == "permitted_with_attribution":
            # Licensed, but the project deliberately does not ship it (for example an
            # unscored pool or a release-scope decision).
            obtain.append(corpus)
        elif state in OBTAIN_STATES:
            obtain.append(corpus)
        else:
            unresolved.append(corpus)
    unresolved.extend(
        corpus for corpus in report.get("corpora", [])
        if corpus.get("redistribution") not in SHIPPABLE_STATES
        and corpus.get("redistribution") not in OBTAIN_STATES
    )
    return ships, obtain, unresolved


def build_statement(report: dict[str, Any], release_id: str, archive_note: str | None) -> dict[str, Any]:
    ships, obtain, unresolved = group_corpora(report)
    valid = report.get("valid", False)

    paragraphs: list[str] = []
    paragraphs.append(
        f"Code, frozen evaluation protocols, and all generated result artifacts for "
        f"{release_id} are deposited in a public archive. The archive contains no "
        f"third-party measurement data."
    )

    if ships:
        names = ", ".join(item["id"] for item in ships)
        paragraphs.append(
            f"Evaluated benchmark data that may be redistributed ({names}) is included "
            f"under its original license with attribution."
        )

    if obtain:
        lines = []
        for item in obtain:
            lines.append(f"{item['id']} ({item['license_name']}): {OBTAIN_STATES.get(item['redistribution'], 'see inventory')}")
        paragraphs.append(
            "The following data are not redistributed and must be obtained from the "
            "original deposit: " + "; ".join(lines) + "."
        )

    training_corpus = next(
        (item for item in obtain if item["id"] == "unikp-kcat-development-corpus"), None
    )
    if training_corpus:
        paragraphs.append(
            "The kcat development corpus used for training is not redistributed. "
            f"{training_corpus['blocking_reason']} Scripts that reconstruct the frozen "
            "split from a locally obtained copy are provided so the evaluation can be "
            "reproduced without redistributing the corpus."
        )

    if not valid:
        paragraphs.append(
            "WARNING: the license audit reported unresolved issues "
            f"({len(report.get('issues', []))} recorded). Resolve them before submission."
        )

    if archive_note:
        paragraphs.append(archive_note)

    return {
        "schema_version": 1,
        "release_id": release_id,
        "audit_valid": valid,
        "shipped_corpora": [item["id"] for item in ships],
        "obtain_elsewhere_corpora": [item["id"] for item in obtain],
        "unresolved_corpora": [item["id"] for item in unresolved],
        "draft_statement": "\n\n".join(paragraphs),
        "review_required": [
            "Confirm the archive DOI and repository URL before submission.",
            "Confirm redistribution terms with the data provider where the audit says not_verified.",
            "Confirm the Author Summary and data statement length limits for the target journal.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--audit",
        type=Path,
        default=Path("artifacts/license-audit/development-corpus-license-audit.json"),
    )
    parser.add_argument("--release-id", default="software-implementation-v6")
    parser.add_argument("--archive-note", default=None)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    audit_path = args.audit if args.audit.is_absolute() else root / args.audit
    if not audit_path.is_file():
        print(
            json.dumps(
                {
                    "valid": False,
                    "issues": [
                        f"missing audit: {args.audit.as_posix()}",
                        "run: python scripts/audit_development_corpus_license.py "
                        "--json-out artifacts/license-audit/development-corpus-license-audit.json",
                    ],
                },
                indent=2,
            )
        )
        return 1

    report = json.loads(audit_path.read_text(encoding="utf-8"))
    statement = build_statement(report, args.release_id, args.archive_note)

    if args.out:
        out_path = args.out if args.out.is_absolute() else root / args.out
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            "# Data availability statement (draft)\n\n"
            + statement["draft_statement"]
            + "\n\n## Review checklist\n\n"
            + "\n".join(f"- {item}" for item in statement["review_required"])
            + "\n",
            encoding="utf-8",
        )
        statement["written_to"] = args.out.as_posix()

    print(json.dumps(statement, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
