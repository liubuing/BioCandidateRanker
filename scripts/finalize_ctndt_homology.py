"""Synchronize CtNDT candidate statuses with its already-frozen homology audit."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "biostudies-S-EPMC12645435"
)
CSV_PATH = SOURCE / "candidate_records.csv"
AUDIT_PATH = SOURCE / "homology-audit.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    candidate = audit.get("candidate_records", {})
    expected_hash = candidate.get("sha256")
    if expected_hash != sha256(CSV_PATH):
        raise ValueError("CtNDT candidate CSV differs from the frozen pre-finalization audit")
    if audit.get("homology_hits", {}).get("count") != 0:
        raise ValueError("CtNDT homology audit is not a zero-hit audit")
    if audit.get("accepted_records") != candidate.get("count"):
        raise ValueError("CtNDT audit does not accept every candidate row")

    with CSV_PATH.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames
    if not fieldnames or len(rows) != candidate["count"]:
        raise ValueError("CtNDT candidate CSV shape differs from the frozen audit")
    if any(row["status_at_normalization"] != "pending_homology" for row in rows):
        raise ValueError("CtNDT rows are not in the expected pending state")
    for row in rows:
        row["status_at_normalization"] = "accepted_homology_cold_pool"
    with CSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    candidate["size_bytes"] = CSV_PATH.stat().st_size
    candidate["sha256"] = sha256(CSV_PATH)
    AUDIT_PATH.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
