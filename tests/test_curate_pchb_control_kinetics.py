import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "europepmc-PMC12329711"
)


def test_unresolved_control_constructs_are_terminally_excluded() -> None:
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    unresolved = [row for row in rows if row["sequence_id"] in {
        "control-mtcm", "control-mtcm-v73e", "control-mtcm-r49l", "control-ecm"
    }]
    assert len(unresolved) == 4
    assert {row["status_at_normalization"] for row in unresolved} == {
        "excluded_construct_unresolved"
    }
    assert not any(row["status_at_normalization"].startswith("pending") for row in rows)
    assert all("not the exact assayed" in row["construct_sequence_scope"] for row in rows)

    evidence = json.loads((SOURCE / "construct-resolution.json").read_text(encoding="ascii"))
    assert evidence["exact_assayed_chains_reconstructed"] == 0
    assert evidence["terminal"] is True
    assert len(evidence["rows"]) == 4
    assert all("CATALYTIC_ACCESSION_NOT_SUBSTITUTED" in row["evidence_codes"] for row in evidence["rows"])

    audit = json.loads((SOURCE / "homology-audit.json").read_text(encoding="ascii"))
    assert audit["readiness_gate_passes"] is True
    assert audit["blockers"] == []
    assert audit["accepted_records"] == 0
