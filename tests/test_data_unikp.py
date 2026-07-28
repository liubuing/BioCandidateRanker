import json
import math
import tempfile
import unittest
from pathlib import Path

from biocandidate.data import aggregate_pair_measurements, read_unikp_json


class UniKPAdapterTest(unittest.TestCase):
    def test_filters_tracks_rows_logs_values_and_audits_duplicates(self):
        rows = [
            {
                "Sequence": "AAAA",
                "Substrate": "ethanol",
                "Smiles": "CCO",
                "Organism": "E. coli",
                "ECNumber": "1.1.1.1",
                "Type": "wildtype",
                "Value": 100,
            },
            {"Sequence": "BBBB", "Substrate": "mixture", "Smiles": "CC.O", "Value": 2},
            {"Sequence": "CCCC", "Substrate": "amine", "Smiles": "CCN", "Value": 0},
            {
                "Sequence": "AAAA",
                "Substrate": "ethanol",
                "Smiles": "CCO",
                "Organism": "E. coli",
                "ECNumber": "1.1.1.1",
                "Type": "wildtype",
                "Value": 10,
            },
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unikp.json"
            path.write_text(json.dumps(rows), encoding="utf-8")
            before = path.read_bytes()
            result = read_unikp_json(path)
            self.assertEqual(path.read_bytes(), before)

        self.assertEqual([record.source_row for record in result.records], [0, 3])
        self.assertEqual(result.records[0].log10_kcat, 2.0)
        self.assertEqual(result.records[0].substrate_name, "ethanol")
        self.assertTrue(math.isclose(result.records[1].log10_kcat, 1.0))
        self.assertEqual(result.audit.total_rows, 4)
        self.assertEqual(result.audit.accepted_rows, 2)
        self.assertEqual(len(result.audit.rejected), 2)
        self.assertEqual(result.audit.duplicates[0].first_source_row, 0)
        self.assertTrue(result.audit.duplicates[0].conflicting_value)
        self.assertEqual(result.audit.conflict_count, 1)
        aggregated = aggregate_pair_measurements(result.records)
        self.assertEqual(len(aggregated), 1)
        self.assertEqual(aggregated[0].replicate_count, 2)
        self.assertEqual(aggregated[0].source_rows, (0, 3))
        self.assertTrue(math.isclose(aggregated[0].log10_kcat, 1.5))

    def test_rejects_non_record_json(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text('{"not_data": 1}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "array of record objects"):
                read_unikp_json(path)


if __name__ == "__main__":
    unittest.main()
