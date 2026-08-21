"""Tests for ESM-2 protein encoder and LR scheduler."""
from __future__ import annotations

import math
import unittest

import torch

from biocandidate.config import ModelConfig
from biocandidate.training import create_lr_scheduler


def _esm2_available() -> bool:
    try:
        import esm  # noqa: F401
        return True
    except ImportError:
        return False


class TestModelConfigESM2(unittest.TestCase):
    def test_esm2_encoder_accepted(self):
        config = ModelConfig(protein_encoder="esm2", d_model=64, num_heads=4)
        self.assertEqual(config.protein_encoder, "esm2")
        self.assertEqual(config.esm2_model_name, "esm2_t6_8M_UR50D")
        self.assertTrue(config.esm2_frozen)

    def test_esm2_config_roundtrip(self):
        config = ModelConfig(
            protein_encoder="esm2",
            esm2_model_name="esm2_t12_35M_UR50D",
            esm2_frozen=False,
            d_model=64,
            num_heads=4,
        )
        restored = ModelConfig.from_dict(config.to_dict())
        self.assertEqual(restored.protein_encoder, "esm2")
        self.assertEqual(restored.esm2_model_name, "esm2_t12_35M_UR50D")
        self.assertFalse(restored.esm2_frozen)

    def test_invalid_encoder_rejected(self):
        with self.assertRaises(ValueError):
            ModelConfig(protein_encoder="bert")

    def test_default_encoder_unchanged(self):
        config = ModelConfig()
        self.assertEqual(config.protein_encoder, "chunk_transformer")


class TestLRScheduler(unittest.TestCase):
    def test_no_warmup_cosine(self):
        model = torch.nn.Linear(4, 2)
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
        scheduler = create_lr_scheduler(optimizer, warmup_epochs=0, total_epochs=10)
        lrs = []
        for _ in range(10):
            lrs.append(optimizer.param_groups[0]["lr"])
            scheduler.step()
        # First LR should be base LR (no warmup)
        self.assertAlmostEqual(lrs[0], 1e-3, places=6)
        # Last LR should be near min_lr_ratio * base
        self.assertLess(lrs[-1], lrs[0] * 0.1)

    def test_warmup_ramps_up(self):
        model = torch.nn.Linear(4, 2)
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
        scheduler = create_lr_scheduler(optimizer, warmup_epochs=5, total_epochs=20)
        lrs = []
        for _ in range(20):
            lrs.append(optimizer.param_groups[0]["lr"])
            scheduler.step()
        # LR should increase during warmup
        self.assertLess(lrs[0], lrs[4])
        # Peak should be at end of warmup
        self.assertAlmostEqual(lrs[4], 1e-3, places=6)
        # Should decay after warmup
        self.assertLess(lrs[-1], lrs[5])

    def test_min_lr_floor(self):
        model = torch.nn.Linear(4, 2)
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
        scheduler = create_lr_scheduler(
            optimizer, warmup_epochs=0, total_epochs=100, min_lr_ratio=0.05)
        for _ in range(200):
            scheduler.step()
        final_lr = optimizer.param_groups[0]["lr"]
        self.assertGreaterEqual(final_lr, 1e-3 * 0.05 - 1e-8)


@unittest.skipUnless(
    _esm2_available(), "fair-esm not installed")
class TestESM2Encoder(unittest.TestCase):
    def test_output_shape(self):
        from biocandidate.model.encoders import ESM2ProteinEncoder
        encoder = ESM2ProteinEncoder(d_model=64, chunk_size=32, dropout=0.0)
        tokens = torch.randint(2, 22, (2, 50))
        mask = torch.ones(2, 50, dtype=torch.bool)
        mask[1, 40:] = False
        output, output_mask = encoder(tokens, mask)
        expected_chunks = math.ceil(50 / 32)
        self.assertEqual(output.shape, (2, expected_chunks, 64))
        self.assertEqual(output_mask.shape, (2, expected_chunks))

    def test_padding_zeros_and_stability(self):
        from biocandidate.model.encoders import ESM2ProteinEncoder
        encoder = ESM2ProteinEncoder(d_model=64, chunk_size=64, dropout=0.0)
        encoder.eval()
        tokens = torch.randint(2, 22, (1, 30))
        mask = torch.ones(1, 30, dtype=torch.bool)
        out1, mask1 = encoder(tokens, mask)
        # Pad to 64 with zeros — same chunk count (1)
        padded_tokens = torch.cat((tokens, torch.zeros(1, 34, dtype=torch.long)), dim=1)
        padded_mask = torch.cat((mask, torch.zeros(1, 34, dtype=torch.bool)), dim=1)
        out2, mask2 = encoder(padded_tokens, padded_mask)
        # Both produce 1 chunk
        self.assertEqual(out1.shape, out2.shape)
        # EOS now sits immediately after the valid residues, so representations (and
        # therefore chunk-pooled outputs) must be identical regardless of trailing padding.
        self.assertTrue(torch.allclose(out1, out2, atol=1e-4),
                        f"padding-dependent output: {out1.flatten()[:4]} vs {out2.flatten()[:4]}")
        # Padded positions still masked out of chunk pooling.

    def test_output_is_invariant_to_batch_composition(self):
        from biocandidate.model.encoders import ESM2ProteinEncoder
        encoder = ESM2ProteinEncoder(d_model=64, chunk_size=4, dropout=0.0)
        encoder.eval()
        torch.manual_seed(0)
        width = 20  # collator pads all rows of a batch to the same width
        token_ids = torch.randint(2, 22, (1, width))
        short_len = 6
        short_mask = torch.ones(1, width, dtype=torch.bool)
        short_mask[0, short_len:] = False
        solo, _ = encoder(token_ids, short_mask)

        # Same 6-residue sequence batched beside a genuinely longer 15-residue sequence.
        batch = torch.cat((token_ids, torch.randint(2, 22, (1, width))), dim=0)
        long_mask = torch.ones(2, width, dtype=torch.bool)
        long_mask[0, short_len:] = False
        long_mask[1, 15:] = False
        joined, _ = encoder(batch, long_mask)

        # EOS sits immediately after each row's own residues, so the short sequence's
        # framing and therefore its chunk-pooled embedding are unchanged by batch neighbors.
        self.assertEqual(solo.shape[1], joined.shape[1])
        self.assertTrue(
            torch.allclose(solo.flatten(), joined[0].flatten(), atol=1e-4),
            "short-sequence embedding changed when batched beside a longer sequence")

    def test_padded_chunks_are_zero(self):
        from biocandidate.model.encoders import ESM2ProteinEncoder
        encoder = ESM2ProteinEncoder(d_model=64, chunk_size=32, dropout=0.0)
        encoder.eval()
        # 40 tokens -> 2 chunks; second chunk has only 8 valid positions
        tokens = torch.randint(2, 22, (1, 40))
        mask = torch.ones(1, 40, dtype=torch.bool)
        mask[0, 35:] = False  # only 3 valid in second chunk
        out, out_mask = encoder(tokens, mask)
        self.assertEqual(out.shape[1], 2)
        self.assertTrue(out_mask[0, 0].item())
        self.assertTrue(out_mask[0, 1].item())
        # A fully-padded chunk should produce zeros
        tokens2 = torch.randint(2, 22, (1, 64))
        mask2 = torch.zeros(1, 64, dtype=torch.bool)
        mask2[0, :10] = True  # only first 10 valid, second chunk fully padded
        out2, out_mask2 = encoder(tokens2, mask2)
        self.assertTrue(out_mask2[0, 0].item())
        self.assertFalse(out_mask2[0, 1].item())
        torch.testing.assert_close(
            out2[0, 1], torch.zeros(64), atol=1e-6, rtol=0.0)

    def test_frozen_backbone_no_grad(self):
        from biocandidate.model.encoders import ESM2ProteinEncoder
        encoder = ESM2ProteinEncoder(d_model=64, chunk_size=32, dropout=0.0,
                                     esm2_frozen=True)
        for param in encoder.backbone.parameters():
            self.assertFalse(param.requires_grad)
        # Projection should still be trainable
        trainable = [p for p in encoder.projection.parameters() if p.requires_grad]
        self.assertGreater(len(trainable), 0)

    def test_full_model_with_esm2(self):
        from biocandidate.model.ranker import BioCandidateRanker
        config = ModelConfig(
            d_model=64, num_heads=4, protein_encoder="esm2",
            protein_chunk_size=32, task_names=("log10_kcat",))
        model = BioCandidateRanker(config)
        batch = {
            "sequence_tokens": torch.randint(2, 22, (2, 40)),
            "sequence_mask": torch.ones(2, 40, dtype=torch.bool),
            "atom_features": torch.randint(0, 10, (6, 6)),
            "edge_index": torch.tensor([[0, 1, 2, 3, 4, 5], [1, 0, 3, 2, 5, 4]]),
            "edge_features": torch.randint(0, 4, (6, 3)),
            "graph_batch": torch.tensor([0, 0, 0, 1, 1, 1]),
            "context_ids": torch.randint(1, 100, (2, 4)),
            "fba_context": torch.zeros(2, 8),
            "fba_context_mask": torch.zeros(2, dtype=torch.bool),
        }
        output = model(batch)
        self.assertEqual(output["mean"].shape, (2, 1))
        self.assertTrue(torch.isfinite(output["mean"]).all())


if __name__ == "__main__":
    unittest.main()
