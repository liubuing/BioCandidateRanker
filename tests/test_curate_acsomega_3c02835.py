"""Pins the ACS Omega 3c02835 curation verdict: all eleven reported rows are
retained as evidence yet fail closed (zero admitted records), with per-family
construct and homology decisions recorded row by row.

中文：钉住 ACS Omega 3c02835 策展裁定：11 行全部保留为证据但零入选，逐行记录构建体与同源决策。
"""

import csv
import json
from pathlib import Path

import pytest


SOURCE = Path("artifacts/external/temporal-absolute-kinetics/europepmc-PMC10620931")

pytestmark = pytest.mark.skipif(
    not (SOURCE / "excluded_records.csv").is_file(),
    reason="local curation evidence under gitignored artifacts/ is required",
)


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
