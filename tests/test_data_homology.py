import tempfile
import unittest
from pathlib import Path

from biocandidate.data import (
    EnzymeSubstrateRecord,
    apply_split_manifest,
    assign_homology_splits,
    read_mmseqs_clusters,
    sequence_id,
    write_split_manifest,
    write_unique_fasta,
)


def record(sequence, row):
    return EnzymeSubstrateRecord(
        sequence=sequence,
        substrate_smiles="CCO",
        organism="host",
        ec="1.1.1.1",
        enzyme_type="wildtype",
        source_dataset="test",
        source_row=row,
    )


class HomologySplitTest(unittest.TestCase):
    def test_cluster_never_crosses_split_and_manifest_round_trips(self):
        records = [record("AAAA", 0), record("AAAT", 1), record("CCCC", 2), record("DDDD", 3)]
        first = sequence_id("AAAA")
        clusters = {
            first: first,
            sequence_id("AAAT"): first,
            sequence_id("CCCC"): sequence_id("CCCC"),
            sequence_id("DDDD"): sequence_id("DDDD"),
        }
        assigned = assign_homology_splits(records, clusters, ratios=(0.5, 0.25, 0.25))
        self.assertEqual(assigned[0].split, assigned[1].split)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "split.json"
            write_split_manifest(assigned, clusters, path, parameters={"test": True})
            restored = apply_split_manifest(records, path)
        self.assertEqual([item.split for item in restored], [item.split for item in assigned])

    def test_fasta_and_cluster_parser(self):
        records = [record("AAAA", 0), record("AAAA", 1), record("CCCC", 2)]
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "sequences.fasta"
            sequences = write_unique_fasta(records, fasta)
            self.assertEqual(len(sequences), 2)
            cluster_file = Path(directory) / "clusters.tsv"
            identifiers = sorted(sequences)
            cluster_file.write_text(
                f"{identifiers[0]}\t{identifiers[0]}\n{identifiers[0]}\t{identifiers[1]}\n",
                encoding="utf-8",
            )
            parsed = read_mmseqs_clusters(cluster_file)
        self.assertEqual(len(parsed), 2)


if __name__ == "__main__":
    unittest.main()
