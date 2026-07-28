import csv
import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biocandidate.global_temporal_audit import build_global_audit


ACCEPTED = "accepted_homology_cold_pool"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_protocol(path: Path, *, cap=20, records=1, families=1, substrates=1):
    path.write_text(
        json.dumps(
            {
                "leakage_gates": {
                    "mmseqs_min_identity": 0.3,
                    "mmseqs_coverage": 0.8,
                    "mmseqs_coverage_mode": 0,
                },
                "diversity_and_size_gate": {
                    "minimum_records": records,
                    "minimum_mmseqs_families": families,
                    "minimum_unique_substrates": substrates,
                    "maximum_records_per_family": cap,
                },
            }
        ),
        encoding="utf-8",
    )


def write_source(root: Path, name: str, sequence_ids: list[str], *, duplicate_fasta=False):
    source = root / name
    source.mkdir()
    candidate = source / "candidate_records.csv"
    with candidate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "candidate_id",
                "sequence_id",
                "substrate_pubchem_cid",
                "status_at_normalization",
            ],
        )
        writer.writeheader()
        for index, sequence_id in enumerate(sequence_ids):
            writer.writerow(
                {
                    "candidate_id": f"candidate-{index}",
                    "sequence_id": sequence_id,
                    "substrate_pubchem_cid": str(100 + index),
                    "status_at_normalization": ACCEPTED,
                }
            )
    fasta = source / "construct_sequences.fasta"
    fasta_text = ">seq-a\nMAAA\n>seq-b\nMCCC\n"
    if duplicate_fasta:
        fasta_text += ">seq-a\nMDDD\n"
    fasta.write_text(fasta_text, encoding="ascii")
    audit = {
        "source_id": name,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": "test-version",
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": "a" * 64,
        "candidate_records": len(sequence_ids),
        "construct_sequences_sha256": sha256(fasta),
        "homology_hits_sha256": "b" * 64,
        "accepted_records": len(sequence_ids),
    }
    (source / "homology-audit.json").write_text(json.dumps(audit), encoding="utf-8")
    return source


def write_clusters(path: Path, representatives: list[str]):
    path.write_text(
        "".join(
            f"{representative}\tgseq-{index:06d}\n"
            for index, representative in enumerate(representatives, 1)
        ),
        encoding="ascii",
    )


class GlobalTemporalAuditTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = self.base / "pool"
        self.root.mkdir()
        self.protocol = self.base / "protocol.json"
        write_protocol(self.protocol)

    def tearDown(self):
        self.temporary.cleanup()

    def test_exact_mapping_writes_fasta_manifest_and_counts(self):
        write_source(self.root, "europepmc-PMC12284513", ["seq-a", "seq-b"])
        clusters = self.base / "clusters.tsv"
        write_clusters(clusters, ["gseq-000001", "gseq-000002"])
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None,
            cluster_tsv=clusters,
        )
        self.assertEqual(
            report["counts"],
            {
                "accepted_sources": 1,
                "excluded_or_unresolved_sources": 0,
                "mapped_records_before_cap": 2,
                "mapped_component_sequences": 2,
                "global_families": 2,
                "records_after_global_family_cap": 2,
                "substrates_after_global_family_cap": 2,
            },
        )
        self.assertTrue(report["readiness_gate_passes"])
        manifest = list(csv.DictReader((self.base / "output/global-accepted-mapping.csv").open()))
        self.assertEqual([row["source_sequence_id"] for row in manifest], ["seq-a", "seq-b"])
        self.assertEqual([row["global_family_id"] for row in manifest], ["family-0001", "family-0002"])

    def test_global_cap_is_stable_and_label_independent(self):
        write_protocol(self.protocol, cap=2, records=2)
        write_source(self.root, "europepmc-PMC12284513", ["seq-a", "seq-a", "seq-a"])
        clusters = self.base / "clusters.tsv"
        write_clusters(clusters, ["gseq-000001"] * 3)
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None,
            cluster_tsv=clusters,
        )
        self.assertEqual(report["counts"]["global_families"], 1)
        self.assertEqual(report["counts"]["records_after_global_family_cap"], 2)
        manifest = list(csv.DictReader((self.base / "output/global-accepted-mapping.csv").open()))
        self.assertEqual(
            [row["retained_after_family_cap"] for row in manifest], ["True", "True", "False"]
        )
        self.assertFalse(report["parameters"]["cap_uses_labels"])

    def test_ambiguous_duplicate_fasta_id_excludes_entire_source(self):
        write_source(
            self.root, "europepmc-PMC12284513", ["seq-a"], duplicate_fasta=True
        )
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None
        )
        self.assertEqual(report["counts"]["mapped_records_before_cap"], 0)
        self.assertEqual(report["counts"]["excluded_or_unresolved_sources"], 1)
        self.assertTrue(any("duplicate" in blocker for blocker in report["blockers"]))

    def test_undeclared_source_schema_is_not_guessed(self):
        write_source(self.root, "unknown-source", ["seq-a"])
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None
        )
        self.assertEqual(report["counts"]["accepted_sources"], 0)
        self.assertIn("no explicit source mapping schema", report["sources"][0]["blockers"])

    def test_construct_fasta_requires_exact_audit_hash(self):
        source = write_source(self.root, "europepmc-PMC12284513", ["seq-a"])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit.pop("construct_sequences_sha256")
        audit["homology_queries_sha256"] = "c" * 64
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None
        )
        self.assertEqual(report["counts"]["mapped_records_before_cap"], 0)
        self.assertTrue(any("construct_sequences_sha256" in item for item in report["blockers"]))

    def test_cluster_tsv_must_cover_every_component_exactly(self):
        write_source(self.root, "europepmc-PMC12284513", ["seq-a", "seq-b"])
        clusters = self.base / "clusters.tsv"
        write_clusters(clusters, ["gseq-000001"])
        report = build_global_audit(
            self.root, self.protocol, self.base / "output", mmseqs=None,
            cluster_tsv=clusters,
        )
        self.assertIsNone(report["counts"]["global_families"])
        self.assertTrue(any("membership mismatch" in blocker for blocker in report["blockers"]))

    @patch("biocandidate.global_temporal_audit.run_mmseqs_easy_cluster")
    @patch("biocandidate.global_temporal_audit.subprocess.run")
    def test_mmseqs_rerun_replaces_managed_workspace(self, run, cluster):
        write_source(self.root, "europepmc-PMC12284513", ["seq-a"])
        output = self.base / "output"
        stale = output / "mmseqs" / "stale-db"
        stale.parent.mkdir(parents=True)
        stale.write_text("stale", encoding="ascii")
        cluster_tsv = self.base / "clusters.tsv"
        write_clusters(cluster_tsv, ["gseq-000001"])
        run.return_value = subprocess.CompletedProcess([], 0, "5d152c612b6ad2a56f657b7a02c127eceaea2a75\n", "")
        cluster.return_value = cluster_tsv
        report = build_global_audit(self.root, self.protocol, output, mmseqs="mmseqs")
        self.assertFalse(stale.exists())
        self.assertEqual(report["counts"]["global_families"], 1)


if __name__ == "__main__":
    unittest.main()
