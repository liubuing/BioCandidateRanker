import tempfile
import unittest
from pathlib import Path

import torch
from torch import nn

from biocandidate import BioCandidateRanker, ModelConfig
from biocandidate.data import EnzymeSubstrateRecord
from biocandidate.training import (
    CampaignBatchSampler,
    load_checkpoint,
    run_epoch,
    run_ranking_epoch,
    save_checkpoint,
)


class CheckpointTest(unittest.TestCase):
    def test_checkpoint_preserves_config_and_data_identity(self):
        config = ModelConfig(
            d_model=16, num_heads=4, protein_layers=1,
            molecule_layers=1, fusion_layers=1)
        model = BioCandidateRanker(config)
        optimizer = torch.optim.AdamW(model.parameters())
        manifest = {"identity": {"sha256": "a" * 64, "size_bytes": 10, "row_count": 2}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.pt"
            save_checkpoint(
                path, model, optimizer, epoch=2,
                data_manifest=manifest, metrics={"loss": 1.0})
            restored, payload = load_checkpoint(path, torch.device("cpu"))
        self.assertEqual(restored.config.d_model, 16)
        self.assertEqual(payload["data_manifest"], manifest)
        self.assertEqual(payload["epoch"], 2)


class CampaignBatchSamplerTest(unittest.TestCase):
    def test_batches_never_mix_campaigns_and_shuffle_is_reproducible(self):
        records = [
            EnzymeSubstrateRecord(
                sequence="ACD",
                substrate_smiles="CCO",
                organism="",
                ec="",
                enzyme_type="variant",
                candidate_id=str(index),
                campaign_group=campaign,
            )
            for index, campaign in enumerate(("a", "b", "a", "b", "a", "c"))
        ]
        sampler = CampaignBatchSampler(records, 2, shuffle=True, seed=17)
        batches = list(sampler)
        for batch in batches:
            self.assertEqual(len({records[index].campaign_group for index in batch}), 1)
        self.assertEqual(batches, list(CampaignBatchSampler(records, 2, shuffle=True, seed=17)))
        sampler.set_epoch(1)
        self.assertNotEqual(batches, list(sampler))


class TinyRankingModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.weight = nn.Parameter(torch.tensor(0.0))

    def forward(self, batch):
        return {"mean": batch["features"] * self.weight}


class TinyRegressionModel(nn.Module):
    def __init__(self, uncertainty_mode):
        super().__init__()
        self.config = ModelConfig(uncertainty_mode=uncertainty_mode)

    def forward(self, batch):
        mean = batch["features"]
        return {"mean": mean, "log_variance": torch.zeros_like(mean)}


class RegressionEpochTest(unittest.TestCase):
    def test_evaluation_loss_is_invariant_to_batch_partitioning(self):
        full_batch = {
            "features": torch.zeros(3, 2),
            "labels": torch.tensor([[1.0, 100.0], [3.0, 4.0], [5.0, 6.0]]),
            "label_mask": torch.tensor([
                [True, False], [True, True], [False, True],
            ]),
            "evidence_weight": torch.tensor([1.0, 2.0, 0.5]),
        }
        partitioned = [
            {name: value[:1] for name, value in full_batch.items()},
            {name: value[1:] for name, value in full_batch.items()},
        ]

        for uncertainty_mode in ("fixed_variance_mse", "heteroscedastic"):
            with self.subTest(uncertainty_mode=uncertainty_mode):
                model = TinyRegressionModel(uncertainty_mode)
                combined_metrics = run_epoch(model, [full_batch], torch.device("cpu"))
                partitioned_metrics = run_epoch(model, partitioned, torch.device("cpu"))
                self.assertAlmostEqual(combined_metrics.loss, partitioned_metrics.loss)
                self.assertEqual(combined_metrics.observed_labels, 4)
                self.assertEqual(partitioned_metrics.observed_labels, 4)

    def test_fails_when_whole_epoch_has_no_observed_labels(self):
        batch = {
            "features": torch.zeros(2, 2),
            "labels": torch.zeros(2, 2),
            "label_mask": torch.zeros(2, 2, dtype=torch.bool),
            "evidence_weight": torch.ones(2),
        }
        model = TinyRegressionModel("fixed_variance_mse")
        with self.assertRaisesRegex(ValueError, "no observed labels"):
            run_epoch(model, [batch], torch.device("cpu"))


class RankingEpochTest(unittest.TestCase):
    def test_skips_empty_pair_batch_and_trains_on_comparable_pairs(self):
        loader = [
            {
                "features": torch.tensor([[1.0], [2.0]]),
                "labels": torch.tensor([[0.5], [0.5]]),
                "label_mask": torch.tensor([[True], [True]]),
                "campaign_ids": torch.tensor([1, 1], dtype=torch.int64),
            },
            {
                "features": torch.tensor([[1.0], [2.0], [3.0]]),
                "labels": torch.tensor([[0.0], [1.0], [0.0]]),
                "label_mask": torch.tensor([[True], [True], [True]]),
                "campaign_ids": torch.tensor([2, 2, 3], dtype=torch.int64),
            },
        ]
        model = TinyRankingModel()
        optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
        metrics = run_ranking_epoch(model, loader, torch.device("cpu"), optimizer)
        self.assertEqual(metrics.pair_count, 1)
        self.assertAlmostEqual(metrics.loss, torch.log(torch.tensor(2.0)).item())
        self.assertGreater(model.weight.item(), 0.0)

    def test_fails_when_loader_has_no_comparable_pairs(self):
        loader = [{
            "features": torch.tensor([[1.0], [2.0]]),
            "labels": torch.tensor([[0.5], [0.5]]),
            "label_mask": torch.tensor([[True], [True]]),
            "campaign_ids": torch.tensor([1, 1], dtype=torch.int64),
        }]
        with self.assertRaisesRegex(ValueError, "no comparable"):
            run_ranking_epoch(TinyRankingModel(), loader, torch.device("cpu"))


if __name__ == "__main__":
    unittest.main()
