import unittest
import json
from unittest.mock import patch

import torch

from biocandidate.data import (
    AminoAcidTokenizer,
    EnzymeSubstrateCollator,
    EnzymeSubstrateRecord,
    FBAFeatureMetadata,
    FluxLabelMetadata,
    FluxLabelType,
    MolecularGraph,
    RDKitUnavailableError,
    smiles_to_sparse_graph,
)


def record(sequence, smiles, row, **labels):
    return EnzymeSubstrateRecord(
        sequence=sequence,
        substrate_smiles=smiles,
        organism="E. coli",
        ec="1.1.1.1",
        enzyme_type="wildtype",
        source_dataset="tiny",
        source_row=row,
        **labels,
    )


class TokenizerTest(unittest.TestCase):
    def test_padding_and_unknown_residue(self):
        tokens, mask = AminoAcidTokenizer().pad(["AC", "AXZ"])
        self.assertEqual(tokens.shape, (2, 3))
        self.assertEqual(tokens[0, 2].item(), 0)
        self.assertEqual(tokens[1, 2].item(), 1)
        self.assertEqual(mask.tolist(), [[True, True, False], [True, True, True]])


class GraphCollatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            smiles_to_sparse_graph("CCO")
        except RDKitUnavailableError as exc:
            raise unittest.SkipTest(str(exc))

    def test_sparse_graph_has_directed_edges(self):
        graph = smiles_to_sparse_graph("CCO")
        self.assertEqual(graph.atom_features.shape, (3, 6))
        self.assertEqual(graph.edge_index.shape, (2, 4))
        self.assertEqual(graph.edge_features.shape, (4, 3))

    def test_collates_graph_context_and_masked_multitask_labels(self):
        records = [
            record("ACD", "CCO", 0, log10_kcat=1.5),
            record("AC", "O", 1, log10_km=-3.0, log10_kcat_per_km=4.0),
            record("ACE", "CC", 2, log10_activity=3.0),
        ]
        collator = EnzymeSubstrateCollator(context_buckets=64)
        batch = collator(records)
        again = collator(records)
        self.assertEqual(batch["sequence_tokens"].shape, (3, 3))
        self.assertEqual(batch["atom_features"].shape[0], 6)
        self.assertEqual(batch["graph_batch"].tolist(), [0, 0, 0, 1, 2, 2])
        self.assertEqual(batch["labels"].shape, (3, 5))
        self.assertEqual(batch["edge_features"].shape, (6, 3))
        self.assertEqual(
            batch["label_mask"].tolist(),
            [
                [True, False, False, False, False],
                [False, True, True, False, False],
                [False, False, False, True, False],
            ],
        )
        self.assertTrue(torch.equal(batch["context_ids"], again["context_ids"]))

    def test_empty_and_exact_width_fba_context(self):
        metadata = FBAFeatureMetadata(
            schema_version=1,
            feature_ids=("EX_glc", "BIOMASS"),
            model_id="model-sha256:abc",
            solver_id="highs-1.7",
            objective_id="BIOMASS",
            condition_id="minimal-aerobic-v1",
        )
        records = [
            record("ACD", "CCO", 0),
            record("ACE", "CC", 1, fba_context=(1.5, -0.25), fba_feature_metadata=metadata),
        ]
        batch = EnzymeSubstrateCollator(context_buckets=64, fba_context_dim=2)(records)
        self.assertEqual(batch["fba_context"].tolist(), [[0.0, 0.0], [1.5, -0.25]])
        self.assertEqual(batch["fba_context_mask"].tolist(), [False, True])

    def test_rejects_fba_context_width_mismatch_instead_of_padding_or_truncating(self):
        metadata = FBAFeatureMetadata(
            schema_version=1,
            feature_ids=("R1", "R2"),
            model_id="model-v1",
            solver_id="solver-v1",
            objective_id="objective-v1",
            condition_id="condition-v1",
        )
        item = record(
            "ACD", "CCO", 0, fba_context=(1.0, 2.0), fba_feature_metadata=metadata)
        with self.assertRaisesRegex(ValueError, "width 2.*configured width 3"):
            EnzymeSubstrateCollator(context_buckets=64, fba_context_dim=3)([item])


class FBASchemaTest(unittest.TestCase):
    def test_nonempty_context_requires_versioned_identity_metadata(self):
        with self.assertRaisesRegex(ValueError, "requires FBAFeatureMetadata"):
            record("ACD", "CCO", 0, fba_context=(1.0,))
        with self.assertRaisesRegex(ValueError, "schema_version"):
            FBAFeatureMetadata(
                schema_version=2,
                feature_ids=("R1",),
                model_id="model",
                solver_id="solver",
                objective_id="objective",
                condition_id="condition",
            )

    def test_record_rejects_declared_feature_width_mismatch(self):
        metadata = FBAFeatureMetadata(
            schema_version=1,
            feature_ids=("R1", "R2"),
            model_id="model",
            solver_id="solver",
            objective_id="objective",
            condition_id="condition",
        )
        with self.assertRaisesRegex(ValueError, "feature_ids exactly"):
            record("ACD", "CCO", 0, fba_context=(1.0,), fba_feature_metadata=metadata)

    def test_flux_requires_explicit_semantics_and_rejects_target_leakage(self):
        fba_metadata = FBAFeatureMetadata(
            schema_version=1,
            feature_ids=("EX_glc", "TARGET_RXN"),
            model_id="model",
            solver_id="solver",
            objective_id="biomass",
            condition_id="condition",
        )
        with self.assertRaisesRegex(ValueError, "requires FluxLabelMetadata"):
            record("ACD", "CCO", 0, log10_flux=1.0)
        flux_metadata = FluxLabelMetadata(
            label_type=FluxLabelType.SIMULATED,
            target_reaction_id="target_rxn",
            canonical_unit="log10(mmol gDW^-1 h^-1)",
            provenance_json=json.dumps({"simulation_run_id": "run-123"}),
        )
        with self.assertRaisesRegex(ValueError, "target reaction"):
            record(
                "ACD",
                "CCO",
                0,
                fba_context=(1.0, 2.0),
                fba_feature_metadata=fba_metadata,
                log10_flux=1.0,
                flux_label_metadata=flux_metadata,
            )

    def test_flux_metadata_distinguishes_experimental_measurements(self):
        metadata = FluxLabelMetadata(
            label_type="experimental",
            target_reaction_id="EX_product",
            canonical_unit="log10(mmol gDW^-1 h^-1)",
            provenance_json=json.dumps({"assay_id": "assay-7"}),
        )
        item = record(
            "ACD", "CCO", 0, log10_flux=-0.2, flux_label_metadata=metadata)
        self.assertEqual(item.flux_label_metadata.label_type, FluxLabelType.EXPERIMENTAL)
        self.assertEqual(item.flux_label_metadata.provenance, {"assay_id": "assay-7"})

    def test_flux_metadata_rejects_noncanonical_units(self):
        with self.assertRaisesRegex(ValueError, "canonical_unit"):
            FluxLabelMetadata(
                label_type="experimental",
                target_reaction_id="EX_product",
                canonical_unit="mmol gDW^-1 h^-1",
                provenance_json=json.dumps({"assay_id": "assay-7"}),
            )


class CampaignCollatorTest(unittest.TestCase):
    @patch("biocandidate.data.collate.smiles_to_sparse_graph")
    def test_custom_activity_rank_task_and_stable_campaign_ids(self, graph_mock):
        graph_mock.return_value = MolecularGraph(
            atom_features=torch.zeros((1, 6), dtype=torch.long),
            edge_index=torch.zeros((2, 0), dtype=torch.long),
            edge_features=torch.zeros((0, 3), dtype=torch.long),
        )
        records = [
            record("ACD", "CCO", 0, activity_rank=0.25, campaign_group="alpha"),
            record("AC", "O", 1, activity_rank=0.75, campaign_group="alpha"),
            record("ACE", "CC", 2, activity_rank=0.5, campaign_group="beta"),
        ]
        collator = EnzymeSubstrateCollator(context_buckets=64, task_names=("activity_rank",))
        batch = collator(records)
        again = collator(records)
        self.assertEqual(batch["labels"].shape, (3, 1))
        self.assertEqual(batch["labels"].flatten().tolist(), [0.25, 0.75, 0.5])
        self.assertEqual(batch["campaign_ids"].dtype, torch.int64)
        self.assertEqual(batch["campaign_ids"][0], batch["campaign_ids"][1])
        self.assertNotEqual(batch["campaign_ids"][0], batch["campaign_ids"][2])
        self.assertTrue(torch.equal(batch["campaign_ids"], again["campaign_ids"]))


if __name__ == "__main__":
    unittest.main()
