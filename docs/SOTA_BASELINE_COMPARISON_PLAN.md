# SOTA Baseline Comparison Plan (Pre-Submission Experiment)

Date: 2026-09-16. Status: proposed, not yet frozen. This plan closes the most likely
reviewer objection against `manuscript/Bioinformatics_Main_Final_EN.md`: Table 1 contains
only internal baselines, with no published predictor (UniKP, DLKcat, TurNuP, CatPred)
evaluated on the same frozen homology-cold split.

## Goal

Add one row group to Table 1: published kcat predictors evaluated on the frozen
homology-cold test partition (1,646 records), under two clearly separated modes:

- **Mode A — as-published checkpoints (zero-shot).** Answers "how do released models do
  on our cold split?" Expected confound: their training corpora are BRENDA-derived and
  may overlap our test rows; this is a *leakage-favorable upper reference*, not a fair
  comparison. Must be accompanied by an exact-overlap audit (see Audit below).
- **Mode B — retrained on our frozen train partition.** Answers "how do the published
  architectures fare when trained under our split discipline?" This is the fair
  comparison and the primary result. Mode A is secondary context.

## Protocol (freeze before any run)

1. **Inputs.** Protein sequence + substrate SMILES for every row of the frozen split
   (`artifacts/homology-final/homology_split.json`, SHA256-bound; test rows never used
   for any hyperparameter decision).
2. **Labels.** Median-aggregated log10(kcat) on the same conflict policy as all
   reported runs. No re-aggregation.
3. **Metrics.** RMSE / MAE / Pearson on the test partition, identical evaluator code
   path as existing runs (`src/biocandidate/evaluation.py`).
4. **No peeking.** Validation partition may be used for predictor-internal calibration
   (e.g., isotonic scaling of outputs) at most once, before any test evaluation; protocol
   receipt frozen first, following `configs/temporal_absolute_kinetics_protocol.json`
   conventions.
5. **Identity binding.** External weights/archives recorded with SHA256 + version pin;
   any download failure writes a `blocker.json` (fail-closed), no substitution.

## Model worklist

| Predictor | License / source | Mode A | Mode B | Notes |
|---|---|---|---|---|
| DLKcat (Li et al. 2022) | GPL-3.0-only, runnable locally | yes | yes | Precedent exists: already ran for the IMDH landscape; ingestion/validation pattern is `evaluate-lunzer` in `src/biocandidate/cli.py` (line ~529) — add a sibling `evaluate-homology-dlkcat` command. |
| UniKP (Yu et al. 2023) | code+weights on Zenodo (10.5281/zenodo.10115498); verify license before use | yes | yes | Light architecture (pretrained embeddings + tree model); Mode B retraining on 13K rows is cheap. |
| TurNuP (Kroll et al. 2023) | check availability/license | optional | optional | Only if trivially runnable; skip if it drags the schedule. |
| CatPred (Boorla et al. 2025) | open code (maranasgroup/CatPred) | optional | optional | Newest competitor (kcat/Km/Ki); at minimum discuss in Related Work, run if time allows. |

## Overlap audit (mandatory for Mode A)

For each predictor, obtain its published training corpus (DLKcat/UniKP corpora are
public), then compute against our frozen test partition:

- exact protein-sequence overlap counts;
- exact (sequence, SMILES) pair overlap counts;
- MMseqs2 30%/80% homology-hit counts (same pins as existing audits).

Report as a supplementary table next to the Mode A numbers so the leakage asymmetry is
explicit. This audit reuses the exact-overlap machinery already built for the temporal
curation pool.

## Compute budget

Local RTX 5060 Laptop (8 GB) suffices: DLKcat training on ~13K rows is small; UniKP Mode
B is minutes-scale. Estimate 3–5 working days for both modes of DLKcat + UniKP including
the overlap audit, plus 1 day to regenerate the manuscript table and figures.

## Deliverables

1. Extended Table 1 (Mode B primary, Mode A footnoted with overlap-audit reference).
2. New supplementary section: per-predictor overlap audit + ingestion receipts.
3. Updated Discussion 4.1/4.2 paragraph interpreting the fair comparison.
4. Machine-readable run manifests under `artifacts/external/sota-homology-cold/` with
   SHA256 identities of weights, inputs, and outputs.

## Explicit non-goals

- No re-architecture, no tuning of BioCandidateRanker against the test partition.
- No external benchmark claims beyond the internal homology-cold split (the prospective
  192/300 pool remains untouched and unscored).
- No use of the test partition for any selection in either mode.
