import json
import math
import tempfile
import unittest
from pathlib import Path

from biocandidate.data import (
    CANONICAL_KINETIC_UNITS,
    EvidenceTier,
    read_normalized_kinetics_jsonl,
)


def row(labels):
    return {
        "schema_version": 1,
        "sequence": "ACDE",
        "substrate_smiles": "CCO",
        "organism": "E. coli",
        "ec": "1.1.1.1",
        "enzyme_type": "wildtype",
        "source_dataset": "direct-study",
        "labels": labels,
    }


def label(value, unit, **extra):
    return {
        "value": value,
        "unit": unit,
        "evidence": "direct",
        "provenance": {"doi": "10.1000/example", "table": "Table 1"},
        **extra,
    }


class NormalizedKineticsAdapterTest(unittest.TestCase):
    def read_rows(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "kinetics.jsonl"
            path.write_text(
                "".join(json.dumps(item) + "\n" for item in rows), encoding="utf-8")
            return read_normalized_kinetics_jsonl(path)

    def test_converts_km_to_log10_molar_and_retains_metadata(self):
        result = self.read_rows([row({"log10_km": label(250, "uM")})])
        record = result.records[0]
        self.assertAlmostEqual(record.log10_km, math.log10(250e-6))
        self.assertIsNone(record.log10_kcat)
        metadata = record.kinetic_label_metadata[0]
        self.assertEqual(metadata.canonical_unit, CANONICAL_KINETIC_UNITS["log10_km"])
        self.assertEqual(metadata.evidence_tier, EvidenceTier.DIRECT)
        self.assertEqual(metadata.provenance["table"], "Table 1")

    def test_converts_efficiency_units(self):
        rows = [row({"log10_kcat_per_km": label(2.5, "mM^-1 s^-1")})]
        record = self.read_rows(rows).records[0]
        self.assertAlmostEqual(record.log10_kcat_per_km, math.log10(2500))

    def test_accepts_canonical_log_values_and_mixed_sparse_labels(self):
        rows = [
            row({"log10_kcat": label(1.25, "log10(s^-1)")}),
            row({
                "log10_km": label(-4.0, "log10(M)"),
                "log10_kcat_per_km": label(5.0, "log10(M^-1 s^-1)"),
            }),
        ]
        records = self.read_rows(rows).records
        self.assertEqual(records[0].log10_kcat, 1.25)
        self.assertIsNone(records[0].log10_km)
        self.assertEqual(records[1].log10_km, -4.0)
        self.assertEqual(records[1].log10_kcat_per_km, 5.0)

    def test_rejects_unknown_or_dimensionally_wrong_units(self):
        for task, unit in (
            ("log10_km", "mg/mL"),
            ("log10_kcat_per_km", "s^-1"),
            ("log10_kcat", "rpm"),
        ):
            with self.subTest(task=task, unit=unit):
                with self.assertRaisesRegex(ValueError, "unsupported unit"):
                    self.read_rows([row({task: label(1.0, unit)})])

    def test_rejects_nonpositive_linear_values(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            self.read_rows([row({"log10_km": label(0, "M")})])
        with self.assertRaisesRegex(ValueError, "must be numeric"):
            self.read_rows([row({"log10_kcat": label(True, "s^-1")})])

    def test_rejects_unresolved_activity_or_flux_labels(self):
        for task in ("log10_activity", "log10_flux", "activity"):
            with self.subTest(task=task):
                with self.assertRaisesRegex(ValueError, "unsupported kinetic task"):
                    self.read_rows([row({task: label(1.0, "log10(a.u.)")})])

    def test_requires_direct_evidence_and_provenance(self):
        indirect = label(1.0, "s^-1", evidence="inferred")
        with self.assertRaisesRegex(ValueError, "evidence must be 'direct'"):
            self.read_rows([row({"log10_kcat": indirect})])
        missing_provenance = label(1.0, "s^-1")
        missing_provenance["provenance"] = {}
        with self.assertRaisesRegex(ValueError, "provenance"):
            self.read_rows([row({"log10_kcat": missing_provenance})])

    def test_rejects_unsupported_schema_version(self):
        payload = row({"log10_kcat": label(1.0, "s^-1")})
        payload["schema_version"] = 2
        with self.assertRaisesRegex(ValueError, "schema_version"):
            self.read_rows([payload])
        payload["schema_version"] = True
        with self.assertRaisesRegex(ValueError, "schema_version"):
            self.read_rows([payload])


if __name__ == "__main__":
    unittest.main()
