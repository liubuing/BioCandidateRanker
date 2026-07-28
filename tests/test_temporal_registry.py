import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from biocandidate.temporal_registry import build_registry, main


def write_protocol(path: Path, *, minimum_records=1, minimum_substrates=1, minimum_families=1):
    path.write_text(
        json.dumps(
            {
                "leakage_gates": {
                    "mmseqs_min_identity": 0.3,
                    "mmseqs_coverage": 0.8,
                    "mmseqs_coverage_mode": 0,
                },
                "diversity_and_size_gate": {
                    "minimum_records": minimum_records,
                    "minimum_unique_substrates": minimum_substrates,
                    "minimum_mmseqs_families": minimum_families,
                },
            }
        ),
        encoding="utf-8",
    )


def write_source(root: Path, name: str, statuses: list[str], *, complex_source=False):
    source = root / name
    source.mkdir()
    candidate = source / "candidate_records.csv"
    with candidate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["candidate_id", "substrate_pubchem_cid", "status_at_normalization"],
        )
        writer.writeheader()
        for index, status in enumerate(statuses):
            writer.writerow(
                {
                    "candidate_id": f"row-{index}",
                    "substrate_pubchem_cid": str(100 + index),
                    "status_at_normalization": status,
                }
            )
    empty_sha = hashlib.sha256(b"").hexdigest()
    audit = {
        "source_id": name,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": "test-version",
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": "a" * 64,
        "candidate_records": len(statuses),
        "construct_sequences_sha256": "b" * 64,
        "homology_hits_sha256": empty_sha,
        "accepted_records": statuses.count("accepted_homology_cold_pool"),
        "candidate_mmseqs_families": 1,
        "family_cluster_sha256": "c" * 64,
    }
    if complex_source:
        audit["complex_homology_rule"] = "exclude if either component hits"
        audit["candidate_component_families"] = 2
    (source / "homology-audit.json").write_text(json.dumps(audit), encoding="utf-8")
    return source


class TemporalRegistryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = self.base / "pool"
        self.root.mkdir()
        self.protocol = self.base / "protocol.json"
        write_protocol(self.protocol)

    def tearDown(self):
        self.temporary.cleanup()

    def test_final_evidence_counts_records_and_global_substrates(self):
        write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        write_source(self.root, "source-b", ["accepted_homology_cold_pool"])
        report = build_registry(
            self.root, self.protocol, family_policy="source-local-additive"
        )
        self.assertEqual(report["counts"], {
            "accepted_records": 2,
            "globally_unique_substrates": 1,
            "families": 2,
        })
        self.assertTrue(report["readiness_gate_passes"])

    def test_accepted_row_without_audit_fails_closed(self):
        source = write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        (source / "homology-audit.json").unlink()
        report = build_registry(self.root, self.protocol)
        self.assertEqual(report["counts"]["accepted_records"], 0)
        self.assertTrue(any("lack finalized homology evidence" in item for item in report["blockers"]))

    def test_pending_csv_disagrees_with_accepted_audit(self):
        source = write_source(self.root, "source-a", ["pending_homology"])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["accepted_records"] = 1
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        report = build_registry(self.root, self.protocol)
        issues = report["sources"][0]["issues"]
        self.assertTrue(any("pending" in issue for issue in issues))
        self.assertFalse(report["sources"][0]["accepted_evidence_final"])
        self.assertEqual(report["counts"]["accepted_records"], 0)

    def test_pending_rows_fail_closed_even_when_audit_reports_zero_accepted(self):
        write_source(self.root, "source-a", ["pending_homology"])
        report = build_registry(self.root, self.protocol)
        source = report["sources"][0]
        self.assertEqual(source["disposition"], "pending_source")
        self.assertTrue(source["blocks_accepted_pool_readiness"])
        self.assertFalse(source["accepted_evidence_final"])
        self.assertTrue(any("remain pending" in item for item in report["blockers"]))

    def test_finalized_zero_candidate_issues_are_reported_but_do_not_block(self):
        source = write_source(self.root, "excluded-source", [])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["candidate_records"] = 4
        audit.pop("construct_sequences_sha256")
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        write_source(self.root, "accepted-source", ["accepted_homology_cold_pool"])

        report = build_registry(
            self.root, self.protocol, family_policy="source-local-additive"
        )
        excluded = next(
            item for item in report["sources"] if item["source_id"] == "excluded-source"
        )
        self.assertEqual(excluded["disposition"], "finalized_zero_candidate_source")
        self.assertFalse(excluded["blocks_accepted_pool_readiness"])
        self.assertIn("audit candidate count does not match CSV", excluded["issues"])
        self.assertIn("audit lacks candidate/query sequence SHA256", excluded["issues"])
        self.assertFalse(any(item.startswith("excluded-source:") for item in report["blockers"]))
        self.assertEqual(report["source_summary"]["nonblocking_source_issues"], 2)
        self.assertTrue(report["readiness_gate_passes"])

    def test_finalized_zero_candidate_source_does_not_require_candidate_csv(self):
        source = write_source(self.root, "excluded-source", [])
        (source / "candidate_records.csv").unlink()

        report = build_registry(self.root, self.protocol)
        excluded = report["sources"][0]
        self.assertEqual(excluded["disposition"], "finalized_zero_candidate_source")
        self.assertFalse(excluded["blocks_accepted_pool_readiness"])
        self.assertIn("missing candidate_records.csv", excluded["issues"])
        self.assertFalse(any(item.startswith("excluded-source:") for item in report["blockers"]))

    def test_protocol_can_require_finalized_zero_candidate_evidence(self):
        source = write_source(self.root, "excluded-source", [])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["candidate_records"] = 1
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        protocol = json.loads(self.protocol.read_text(encoding="utf-8"))
        protocol["registry_readiness"] = {
            "finalized_zero_candidate_sources_must_pass_evidence_checks": True
        }
        self.protocol.write_text(json.dumps(protocol), encoding="utf-8")

        report = build_registry(self.root, self.protocol)
        self.assertTrue(report["sources"][0]["blocks_accepted_pool_readiness"])
        self.assertTrue(any(item.startswith("excluded-source:") for item in report["blockers"]))

    def test_malformed_accepted_source_still_fails_closed(self):
        source = write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit.pop("construct_sequences_sha256")
        audit_path.write_text(json.dumps(audit), encoding="utf-8")

        report = build_registry(self.root, self.protocol)
        self.assertTrue(report["sources"][0]["blocks_accepted_pool_readiness"])
        self.assertEqual(report["counts"]["accepted_records"], 0)
        self.assertTrue(any("lack finalized homology evidence" in item for item in report["blockers"]))

    def test_multicomponent_family_ambiguity_is_not_guessed(self):
        write_source(
            self.root, "source-a", ["accepted_homology_cold_pool"], complex_source=True
        )
        report = build_registry(
            self.root, self.protocol, family_policy="source-local-additive"
        )
        self.assertIsNone(report["counts"]["families"])
        self.assertIn("multicomponent", report["family_policy_note"])
        self.assertFalse(report["readiness_gate_passes"])

    def test_stale_candidate_hash_is_rejected(self):
        source = write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["candidate_records"] = {"count": 1, "sha256": "d" * 64}
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        report = build_registry(self.root, self.protocol)
        self.assertEqual(report["counts"]["accepted_records"], 0)
        self.assertIn("candidate CSV SHA256 differs from audit", report["sources"][0]["issues"])

    def test_source_specific_candidate_count_is_supported(self):
        source = write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        audit_path = source / "homology-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit.pop("candidate_records")
        audit["saturation_qualified_candidate_records"] = 1
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
        report = build_registry(self.root, self.protocol)
        self.assertTrue(report["sources"][0]["accepted_evidence_final"])

    def test_cli_writes_report_and_returns_nonzero_when_blocked(self):
        write_source(self.root, "source-a", ["accepted_homology_cold_pool"])
        output = self.base / "report.json"
        result = main([
            "--root", str(self.root), "--protocol", str(self.protocol), "--output", str(output)
        ])
        self.assertEqual(result, 1)
        self.assertTrue(output.exists())


if __name__ == "__main__":
    unittest.main()
