# Supplementary Material

**Pretrained Protein Representations Improve Enzyme Kinetic Prediction under Homology-Cold Evaluation**

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
- Model parameters (total, including frozen backbone): 12,946,716

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
| **SD** | **0.0183** | **0.0117** | **0.0222** | — |

---

## S4a: Post-Hoc ESM-2 Modality Diagnostics

The internal test was observed before this comparison. The protocol was frozen before the six new runs. All variants use the same ESM-2 backbone, model width, split, target, training schedule, and seeds. Disabled modalities are masked in the task-query fusion model. Results are diagnostic, not independent confirmation. Source and split identities were checked against the full-model reference before summary.

**Supplementary Table S3.** Per-seed test results and three-seed summaries (1,646 records).

| Variant | Seed | RMSE | MAE | Pearson |
|---|---:|---:|---:|---:|
| Protein only | 7 | 1.5013 | 1.1597 | 0.2171 |
| Protein only | 42 | 1.5094 | 1.1682 | 0.2586 |
| Protein only | 123 | 1.4969 | 1.1690 | 0.1766 |
| Protein only | Mean ± SD | 1.5025 ± 0.0064 | 1.1656 ± 0.0051 | 0.2174 ± 0.0410 |
| Protein + substrate | 7 | 1.4571 | 1.1287 | 0.2734 |
| Protein + substrate | 42 | 1.4278 | 1.1028 | 0.3109 |
| Protein + substrate | 123 | 1.4513 | 1.1286 | 0.2529 |
| Protein + substrate | Mean ± SD | 1.4454 ± 0.0156 | 1.1200 ± 0.0149 | 0.2791 ± 0.0294 |
| Full model | 7 | 1.3880 | 1.0683 | 0.3906 |
| Full model | 42 | 1.3754 | 1.0606 | 0.3968 |
| Full model | 123 | 1.4114 | 1.0835 | 0.3556 |
| Full model | Mean ± SD | 1.3916 ± 0.0183 | 1.0708 ± 0.0117 | 0.3810 ± 0.0222 |

Paired mean RMSE changes are −0.0571 ± 0.0213 for protein + substrate versus protein only, and −0.0538 ± 0.0147 for full versus protein + substrate (sample SD across three seeds). Both directions are consistent across seeds. Machine-readable values and input identities are in `artifacts/esm2-modality-diagnostic/summary.json`.

---

## S5: EnzEngDB Selection Protocol

- Archive: 30,050,527 bytes, SHA256 `8013ad81586db2187162aada0709c1cabc7e7d69f03dd5c199776aaf000dd6ea`
- Accepted rows: 245,945 across 160 campaigns
- Label-independent cap: 2,000 rows per campaign by SHA256 of candidate ID
- Homology filtering: MMseqs2 30%/80% against all UniKP sequences
- Final selection: 6,423 rows, 51 campaigns, zero exact sequence overlap
- Selection SHA256: `8e76932fc76b884ea8cf127934683b43b9e896b057b817488a9dc5533cf1fc4b`

---

## S5a: Comparison With Retrained Modern Baselines

All predictors were retrained on the frozen training partition and evaluated once on the 1,646-record test partition. UniKP and DLKcat follow their frozen protocols; CataPro uses the official head and default optimizer over cached ProtT5/MolT5 features with fixed-validation early stopping; CatPred is an adapted pipeline (non-standard residues mapped to a zero embedding, memory-mapped feature cache, per-member seed revision) around frozen 650M-parameter ESM-2 features with official ensemble inference (three group seeds x 10 ensemble members), not a bit-level reproduction of the upstream run.

**Supplementary Table S4.** Three-seed summaries (mean +/- SD over seeds).

| Model | RMSE | MAE | Pearson |
|---|---:|---:|---:|
| BioCandidateRanker ESM-2 (historical implementation) | 1.3916 +/- 0.0183 | 1.0708 +/- 0.0117 | 0.3810 +/- 0.0222 |
| BioCandidateRanker (current-code replay of historical weights) | 1.3956 +/- 0.0174 | 1.0744 +/- 0.0113 | 0.3749 +/- 0.0223 |
| CataPro Mode B | 1.3876 +/- 0.0195 | 1.0650 +/- 0.0105 | 0.4320 +/- 0.0004 |
| CatPred Mode B (adapted) | 1.3893 +/- 0.0175 | 1.0744 +/- 0.0144 | 0.3897 +/- 0.0279 |
| UniKP Mode B | 1.4093 +/- 0.0094 | 1.0785 +/- 0.0041 | 0.3777 +/- 0.0222 |
| DLKcat Mode B | 1.7138 +/- 0.1415 | 1.3351 +/- 0.1251 | 0.1463 +/- 0.0764 |

**Supplementary Table S5.** Exploratory family-level paired bootstrap, RMSE difference (ours minus baseline), 95% interval over 5,000 family resamples (219 families, seed 20260918; no multiple-comparison correction).

| Baseline | Reference implementation | RMSE difference | 95% interval |
|---|---|---:|---|
| CatPred | historical | +0.0023 | [-0.0555, +0.0634] |
| CataPro | current-code replay | +0.0080 | [-0.0731, +0.0967] |
| UniKP | current-code replay | -0.0137 | [-0.0717, +0.0511] |

Every interval crosses zero: neither superiority nor equivalence is established. Machine-readable receipts: `artifacts/catapro-mode-b/internal-test/`, `artifacts/catpred-mode-b/internal-test/`, `artifacts/paired-family-analysis/current-code/paired-bootstrap.json`; frozen protocols `configs/catapro_internal_test_freeze.json`, `configs/catpred_internal_test_protocol.json`.

---

## S6: Claim Boundary

This work establishes an internal homology-cold predictive gain and post-hoc modality diagnostic. It does not establish: calibrated uncertainty, external kcat validation, activity or flux prediction, candidate ranking utility, or publication-grade efficacy for enzyme engineering decisions. The prospective independent benchmark remains incomplete (209/300 records, 28/30 families). No model predictions have been generated for the prospective pool.
