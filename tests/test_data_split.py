import unittest

from biocandidate.data import (
    EnzymeSubstrateRecord,
    RDKitUnavailableError,
    split_records,
    write_cold_split_manifest,
)


def record(sequence, smiles, row):
    return EnzymeSubstrateRecord(
        sequence=sequence,
        substrate_smiles=smiles,
        organism="",
        ec="",
        enzyme_type="unknown",
        source_dataset="tiny",
        source_row=row,
    )


class SplitTest(unittest.TestCase):
    def test_protein_cold_is_deterministic_and_has_no_protein_leakage(self):
        records = [
            record("AAAA", "CCO", 0),
            record("AAAA", "CCN", 1),
            record("CCCC", "CCC", 2),
            record("DDDD", "CCCl", 3),
            record("EEEE", "CCBr", 4),
        ]
        first = split_records(records, "protein_cold", ratios=(0.6, 0.2, 0.2), seed=9)
        second = split_records(records, "protein_cold", ratios=(0.6, 0.2, 0.2), seed=9)
        self.assertEqual(first, second)
        splits_by_protein = {}
        for item in first:
            splits_by_protein.setdefault(item.sequence, set()).add(item.split)
        self.assertTrue(all(len(splits) == 1 for splits in splits_by_protein.values()))
        self.assertIsNone(records[0].split)

    def test_double_cold_keeps_connected_bipartite_component_together(self):
        records = [
            record("AAAA", "c1ccccc1O", 0),
            record("AAAA", "CCO", 1),
            record("BBBB", "c1ccccc1N", 2),
            record("CCCC", "C1CCCCC1", 3),
            record("DDDD", "CCN", 4),
        ]
        try:
            assigned = split_records(records, "double_cold", ratios=(0.6, 0.2, 0.2), seed=2)
        except RDKitUnavailableError as exc:
            self.skipTest(str(exc))
        self.assertEqual(assigned[0].split, assigned[1].split)
        self.assertEqual(assigned[0].split, assigned[2].split)

    def test_rejects_unknown_strategy(self):
        with self.assertRaisesRegex(ValueError, "strategy"):
            split_records([record("AAAA", "CCO", 0)], "random")

    def test_manifest_audits_double_cold_overlap(self):
        records = [
            record("AAAA", "c1ccccc1O", 0),
            record("AAAA", "CCO", 1),
            record("BBBB", "c1ccccc1N", 2),
            record("CCCC", "C1CCCCC1", 3),
        ]
        try:
            assigned = split_records(records, "double_cold", ratios=(0.5, 0.25, 0.25))
        except RDKitUnavailableError as exc:
            self.skipTest(str(exc))
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            payload = write_cold_split_manifest(
                assigned, "double_cold", str(Path(directory) / "split.json"),
                seed=0, source_identity={"sha256": "a" * 64})
        self.assertEqual(payload["audit"]["protein_cross_split_count"], 0)
        self.assertEqual(payload["audit"]["scaffold_cross_split_count"], 0)


if __name__ == "__main__":
    unittest.main()
