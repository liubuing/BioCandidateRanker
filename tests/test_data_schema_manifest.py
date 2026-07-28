import hashlib
import json
import tempfile
import unittest
import zipfile
from dataclasses import FrozenInstanceError
from pathlib import Path

from biocandidate.data import (
    EnzymeSubstrateRecord,
    EvidenceTier,
    FileManifest,
    ManifestVerificationError,
    build_manifest,
    verify_manifest,
)


class SchemaTest(unittest.TestCase):
    def test_record_is_canonical_and_immutable(self):
        record = EnzymeSubstrateRecord(
            sequence=" ac d\n",
            substrate_smiles=" CCO ",
            organism=" E. coli ",
            ec="1.1.1.1",
            enzyme_type="wildtype",
            evidence_tier="direct",
            source_dataset=" tiny ",
            source_row=2,
        )
        self.assertEqual(record.sequence, "ACD")
        self.assertEqual(record.substrate_smiles, "CCO")
        self.assertEqual(record.evidence_tier, EvidenceTier.DIRECT)
        with self.assertRaises(FrozenInstanceError):
            record.split = "train"

    def test_activity_rank_is_optional_and_bounded(self):
        record = EnzymeSubstrateRecord(
            sequence="ACD",
            substrate_smiles="CCO",
            organism="",
            ec="",
            enzyme_type="variant",
            campaign_group=" campaign-a ",
            activity_rank=0.75,
        )
        self.assertEqual(record.campaign_group, "campaign-a")
        self.assertEqual(record.activity_rank, 0.75)
        with self.assertRaisesRegex(ValueError, "between 0 and 1"):
            EnzymeSubstrateRecord(
                sequence="ACD",
                substrate_smiles="CCO",
                organism="",
                ec="",
                enzyme_type="variant",
                activity_rank=1.1,
            )


class ManifestTest(unittest.TestCase):
    def test_verifies_all_fields_without_changing_source(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.json"
            path.write_text(json.dumps([{"x": 1}, {"x": 2}]), encoding="utf-8")
            before = path.read_bytes()
            manifest = FileManifest(
                sha256=hashlib.sha256(before).hexdigest(),
                size_bytes=len(before),
                row_count=2,
            )
            self.assertIsNone(verify_manifest(path, manifest))
            self.assertEqual(path.read_bytes(), before)

    def test_fails_closed_on_any_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.json"
            path.write_text("[]", encoding="utf-8")
            manifest = FileManifest("0" * 64, 2, 0)
            with self.assertRaises(ManifestVerificationError):
                verify_manifest(path, manifest)

    def test_counts_candidate_container(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidates.json"
            path.write_text(
                json.dumps({"candidates": [{"id": 1}, {"id": 2}]}), encoding="utf-8")
            manifest = build_manifest(path)
        self.assertEqual(manifest.row_count, 2)

    def test_counts_csv_data_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.csv"
            path.write_text("id,value\n1,a\n2,b\n", encoding="utf-8")
            manifest = build_manifest(path)
        self.assertEqual(manifest.row_count, 2)

    def test_counts_experiment_rows_in_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "data.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("Data/experiments/a.csv", "id,value\n1,a\n2,b\n")
                archive.writestr("Data/not_an_experiment.csv", "id\nignored\n")
            manifest = build_manifest(path)
        self.assertEqual(manifest.row_count, 2)


if __name__ == "__main__":
    unittest.main()
