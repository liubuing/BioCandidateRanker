import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "europepmc-PMC13210114"
)


def test_rmla_package_fails_closed_without_recalculated_labels() -> None:
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle)) == []

    with (SOURCE / "exclusions.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["variable_substrate"], row["reported_kcat"]) for row in rows] == [
        ("Glc-1-P", "33.685"),
        ("dTTP", "171.4"),
    ]
    assert all(row["candidate_label_created"] == "False" for row in rows)
    assert all(
        "REPORTED_KCAT_VMAX_ENZYME_CONCENTRATION_INCONSISTENT"
        in row["exclusion_codes"]
        for row in rows
    )

    audit = json.loads((SOURCE / "homology-audit.json").read_text(encoding="ascii"))
    assert audit["mmseqs_version"] == "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
    assert audit["homology_hit_alignments"] == 3
    assert audit["accepted_records"] == 0

    construct = json.loads(
        (SOURCE / "construct-resolution.json").read_text(encoding="ascii")
    )
    assert construct["status"] == "unresolved_fail_closed"
    assert construct["native_sequence"]["accession"] == "AYM85572.1"
