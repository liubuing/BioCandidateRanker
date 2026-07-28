import unittest

import torch

from biocandidate import BioCandidateRanker, ModelConfig
from biocandidate.data import (
    EnzymeSubstrateCollator,
    EnzymeSubstrateRecord,
    FBAFeatureMetadata,
    RDKitUnavailableError,
    smiles_to_sparse_graph,
)
from biocandidate.model import (
    masked_multitask_gaussian_loss,
    masked_multitask_mse_loss,
    pairwise_logistic_ranking_loss,
)
from biocandidate.baselines import amino_acid_composition, morgan_fingerprints


class ModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            smiles_to_sparse_graph("CCO")
        except RDKitUnavailableError as exc:
            raise unittest.SkipTest(str(exc))

    def test_forward_backward_and_uncertainty(self):
        records = [
            EnzymeSubstrateRecord(
                sequence="ACDEFGHIK" * 20,
                substrate_smiles="CCO",
                organism="Escherichia coli",
                ec="1.1.1.1",
                enzyme_type="wildtype",
                log10_kcat=1.2,
                reaction="ethanol oxidation",
                fba_context=(0.1,) * 8,
                fba_feature_metadata=FBAFeatureMetadata(
                    schema_version=1,
                    feature_ids=tuple(f"R{i}" for i in range(8)),
                    model_id="ecoli-core-v1",
                    solver_id="highs-1.7",
                    objective_id="biomass",
                    condition_id="aerobic-glucose",
                ),
                source_dataset="test",
                source_row=0,
            ),
            EnzymeSubstrateRecord(
                sequence="MNPQRSTVWY" * 33,
                substrate_smiles="CC(=O)O",
                organism="Saccharomyces cerevisiae",
                ec="2.7.1.1",
                enzyme_type="mutant",
                log10_kcat=0.4,
                source_dataset="test",
                source_row=1,
            ),
        ]
        batch = EnzymeSubstrateCollator(context_buckets=64)(records)
        config = ModelConfig(
            d_model=32,
            num_heads=4,
            protein_layers=1,
            protein_chunk_size=32,
            molecule_layers=2,
            fusion_layers=1,
            context_buckets=64,
            dropout=0.0,
        )
        model = BioCandidateRanker(config)
        output = model(batch)
        self.assertEqual(output["mean"].shape, (2, 5))
        self.assertTrue(torch.isfinite(output["standard_deviation"]).all())
        loss, per_task = masked_multitask_gaussian_loss(
            output["mean"], output["log_variance"], batch["labels"],
            batch["label_mask"], batch["evidence_weight"])
        self.assertEqual(per_task.shape, (5,))
        loss.backward()
        self.assertTrue(torch.isfinite(loss))
        self.assertTrue(any(parameter.grad is not None for parameter in model.parameters()))

    def test_long_sequence_is_chunked(self):
        config = ModelConfig(
            d_model=16, num_heads=4, protein_layers=1,
            protein_chunk_size=64, molecule_layers=1, fusion_layers=1,
            context_buckets=32, dropout=0.0)
        model = BioCandidateRanker(config)
        tokens = torch.ones(2, 1025, dtype=torch.long)
        mask = torch.ones_like(tokens, dtype=torch.bool)
        chunks, chunk_mask = model.protein(tokens, mask)
        self.assertEqual(chunks.shape, (2, 17, 16))
        self.assertTrue(chunk_mask.all())

    def test_rejects_model_without_modalities(self):
        with self.assertRaisesRegex(ValueError, "modality"):
            ModelConfig(use_protein=False, use_molecule=False, use_context=False)

    def test_rejects_nonpositive_fba_context_dimension(self):
        with self.assertRaisesRegex(ValueError, "fba_context_dim"):
            ModelConfig(fba_context_dim=0)

    def test_late_fusion_forward_does_not_instantiate_query_fusion(self):
        records = [EnzymeSubstrateRecord(
            sequence="ACDEFGHIK", substrate_smiles="CCO", organism="", ec="",
            enzyme_type="wildtype", log10_kcat=1.0,
        )]
        batch = EnzymeSubstrateCollator(context_buckets=32)(records)
        config = ModelConfig(
            d_model=16, num_heads=4, protein_layers=1, protein_chunk_size=16,
            molecule_layers=1, fusion_layers=1, context_buckets=32, dropout=0.0,
            fusion_mode="late_concat",
        )
        model = BioCandidateRanker(config)
        output = model(batch)
        self.assertEqual(output["mean"].shape, (1, 5))
        self.assertFalse(hasattr(model, "cross_attention"))
        self.assertTrue(torch.isfinite(output["standard_deviation"]).all())

    def test_global_mean_protein_encoder_has_one_token(self):
        config = ModelConfig(
            d_model=16, num_heads=4, protein_layers=1, protein_chunk_size=4,
            molecule_layers=1, fusion_layers=1, context_buckets=32, dropout=0.0,
            protein_encoder="global_mean",
        )
        model = BioCandidateRanker(config)
        tokens = torch.tensor([[1, 2, 3, 0], [4, 5, 0, 0]])
        mask = tokens != 0
        encoded, encoded_mask = model.protein(tokens, mask)
        self.assertEqual(encoded.shape, (2, 1, 16))
        self.assertTrue(encoded_mask.all())

    def test_fixed_variance_model_has_unit_standard_deviation(self):
        records = [EnzymeSubstrateRecord(
            sequence="ACDE", substrate_smiles="CCO", organism="", ec="",
            enzyme_type="wildtype", log10_kcat=1.0,
        )]
        batch = EnzymeSubstrateCollator(context_buckets=32)(records)
        config = ModelConfig(
            d_model=16, num_heads=4, protein_layers=1, protein_chunk_size=8,
            molecule_layers=1, fusion_layers=1, context_buckets=32, dropout=0.0,
            uncertainty_mode="fixed_variance_mse",
        )
        output = BioCandidateRanker(config)(batch)
        self.assertTrue(torch.equal(output["log_variance"], torch.zeros_like(output["mean"])))
        self.assertTrue(torch.equal(
            output["standard_deviation"], torch.ones_like(output["mean"])))

    def test_masked_mse_ignores_missing_labels(self):
        mean = torch.tensor([[1.0, 10.0]], requires_grad=True)
        labels = torch.tensor([[3.0, 0.0]])
        mask = torch.tensor([[True, False]])
        loss, per_task, numerator, denominator = masked_multitask_mse_loss(
            mean, labels, mask, return_components=True)
        self.assertAlmostEqual(float(loss.detach()), 4.0)
        self.assertAlmostEqual(float(per_task[0].detach()), 4.0)
        self.assertAlmostEqual(float(numerator.detach()), 4.0)
        self.assertAlmostEqual(float(denominator.detach()), 1.0)
        loss.backward()
        self.assertEqual(mean.grad.tolist(), [[-4.0, 0.0]])


class PairwiseRankingLossTest(unittest.TestCase):
    def test_correct_direction_has_lower_loss_and_backpropagates(self):
        labels = torch.tensor([0.0, 1.0, 0.5, 0.5])
        campaign_ids = torch.tensor([1, 1, 2, 2], dtype=torch.int64)
        correct = torch.tensor([-2.0, 2.0, 1.0, -1.0], requires_grad=True)
        reversed_scores = -correct.detach()
        loss, pair_count = pairwise_logistic_ranking_loss(correct, labels, campaign_ids)
        reverse_loss, reverse_count = pairwise_logistic_ranking_loss(
            reversed_scores, labels, campaign_ids)
        self.assertEqual(pair_count, 1)
        self.assertEqual(reverse_count, 1)
        self.assertLess(loss.item(), reverse_loss.item())
        loss.backward()
        self.assertIsNotNone(correct.grad)
        self.assertTrue(torch.isfinite(correct.grad).all())

    def test_no_pairs_returns_differentiable_zero(self):
        scores = torch.tensor([1.0, 2.0], requires_grad=True)
        loss, pair_count = pairwise_logistic_ranking_loss(
            scores, torch.tensor([0.5, 0.5]), torch.tensor([1, 1], dtype=torch.int64))
        self.assertEqual(pair_count, 0)
        loss.backward()
        self.assertEqual(scores.grad.tolist(), [0.0, 0.0])


class ClassicalFeatureTest(unittest.TestCase):
    def test_composition_and_morgan_features_are_stable(self):
        records = [EnzymeSubstrateRecord(
            sequence="ACAA", substrate_smiles="CCO", organism="", ec="",
            enzyme_type="wildtype",
        )]
        composition = amino_acid_composition(records)
        fingerprints = morgan_fingerprints(records, bits=64)
        self.assertEqual(composition.shape, (1, 20))
        self.assertAlmostEqual(float(composition.sum()), 1.0)
        self.assertEqual(fingerprints.shape, (1, 64))
        self.assertGreater(int(fingerprints.sum()), 0)


if __name__ == "__main__":
    unittest.main()
