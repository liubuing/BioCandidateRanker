# Supplementary Material

**Pretrained Protein Representations Unlock Multimodal Fusion Gains for Enzyme Kinetic Prediction**

[Authors anonymized for review]

---

## S1: Frozen Split Identities

- Source SHA256: `13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719`
- Split SHA256: `ae1793afea078b4a2d8b834e9b9cd7c836d43508e26185ed87d26b8b82f96300`
- MMseqs2 version: `5d152c612b6ad2a56f657b7a02c127eceaea2a75`
- Clustering: 30% identity, 80% coverage, coverage mode 0
- Clusters: 2,204 (13,157 train / 1,640 validation / 1,646 test after median aggregation)

---

## S2: ESM-2 Model Configuration

- Model: esm2_t6_8M_UR50D (6 layers, 320-dim, 8M parameters, frozen)
- Projection: Linear(320, 256) + SiLU + Dropout(0.1) + LayerNorm(256)
- Chunk pooling: 128-residue chunks, masked mean
- Sequence truncation: 512 residues maximum
- Training: AdamW, lr=3e-4, cosine schedule, 2-epoch warmup, patience 4
- Effective batch size: 32 (micro-batch 4, accumulation 8)
- Model parameters (trainable): 12,946,716

---

## S3: Architecture Ablation Results (Seed 42, From-Scratch Model)

**Supplementary Table S1.** Single-seed ablation results.

| Ablation | Test RMSE | Delta RMSE | Test Pearson |
|---|---|---|---|
| Task-query baseline | 1.4638 | — | 0.2731 |
| Remove reaction context | 1.5046 | +0.0409 | 0.1908 |
| Shared task query | 1.4911 | +0.0273 | 0.2520 |
| Global-mean protein encoder | 1.4910 | +0.0272 | 0.2278 |

---

## S4: Per-Seed Test Results (ESM-2 Multimodal)

**Supplementary Table S2.** Per-seed test metrics.

| Seed | Test RMSE | Test MAE | Test Pearson | Best Epoch |
|---|---|---|---|---|
| 7 | 1.3880 | 1.0683 | 0.3906 | 2 |
| 42 | 1.3754 | 1.0606 | 0.3968 | 1 |
| 123 | 1.4114 | 1.0835 | 0.3556 | 3 |
| **Mean** | **1.3916** | **1.0708** | **0.3810** | — |
| **SD** | **0.0183** | **0.0117** | **0.0218** | — |

---

## S5: EnzEngDB Selection Protocol

- Archive: 30,050,527 bytes, SHA256 `8013ad81586db2187162aada0709c1cabc7e7d69f03dd5c199776aaf000dd6ea`
- Accepted rows: 245,945 across 160 campaigns
- Label-independent cap: 2,000 rows per campaign by SHA256 of candidate ID
- Homology filtering: MMseqs2 30%/80% against all UniKP sequences
- Final selection: 6,423 rows, 51 campaigns, zero exact sequence overlap
- Selection SHA256: `8e76932fc76b884ea8cf127934683b43b9e896b057b817488a9dc5533cf1fc4b`

---

## S6: Claim Boundary

This work establishes an internal homology-cold predictive gain. It does not establish: calibrated uncertainty, external kcat validation, activity or flux prediction, candidate ranking utility, or publication-grade efficacy for enzyme engineering decisions. The prospective independent benchmark remains incomplete (192/300 records, 25/30 families). No model predictions have been generated for the prospective pool.
