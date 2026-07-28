import csv
import json
from pathlib import Path


SOURCE = Path("artifacts/external/temporal-absolute-kinetics/europepmc-PMC10620931")


def test_all_reported_rows_are_retained_but_none_are_admitted():
    with (SOURCE / "excluded_records.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 11
    assert {row["source_table"] for row in rows} == {"Table 1", "Table 2", "Table 3"}
    assert all(row["status_at_normalization"] == "excluded_fail_closed" for row in rows)
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle)) == []


def test_family_specific_construct_and_homology_decisions():
    construct = json.loads((SOURCE / "construct-resolution.json").read_text(encoding="ascii"))
    audit = json.loads((SOURCE / "homology-audit.json").read_text(encoding="ascii"))
    assert construct["blc23o"]["status"] == "resolved_exact"
    assert construct["ro12ctd"]["status"] == "unresolved_fail_closed"
    assert construct["trypsin"]["status"] == "unresolved_fail_closed"
    assert audit["families"]["blc23o"]["hits"] == 0
    assert audit["families"]["ro12ctd"]["hits"] == 1
    assert audit["families"]["bovine_trypsin"]["hits"] == 7
    assert audit["accepted_records"] == 0


def test_pubchem_and_saturation_are_row_level():
    provenance = json.loads((SOURCE / "provenance.json").read_text(encoding="ascii"))
    saturation = json.loads((SOURCE / "saturation-audit.json").read_text(encoding="ascii"))
    assert provenance["pubchem"] == {
        "3-methylcatechol": {"cid": 340},
        "L-BApNA hydrochloride": {"cid": 16219022},
    }
    assert len(saturation["rows"]) == 11
    assert sum(row["varied_substrate_decision"] == "fail" for row in saturation["rows"]) == 4
    assert sum(row["oxygen_decision"] == "unresolved" for row in saturation["rows"]) == 8
