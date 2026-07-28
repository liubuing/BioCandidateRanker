# Acceptance Criteria

The target is measured by gates, not a subjective percentage.

## G1: Code Correctness

- All unit and integration tests pass on CPU and CUDA.
- Padding invariance, missing-label masking, sparse graph batching, long-sequence
  chunking, checkpoint restoration, and manifest failure behavior are tested.
- No NaN or infinite prediction, loss, gradient, or uncertainty in stress tests.
- A checkpoint reproduces predictions from its saved configuration and data identity.

## G2: Stable Training

- A real-data smoke run completes forward, backward, validation, and atomic checkpointing.
- Full training is resumable and deterministic within documented tolerance.
- Validation reports a train-only mean baseline and task-level RMSE, MAE, Pearson,
  and uncertainty coverage.
- Conflicting source labels require an explicit policy; they cannot pass silently.

## G3: Scalable Architecture

- Protein attention complexity is `O(L * chunk_size)` rather than `O(L^2)`.
- Molecular message passing uses sparse edges and never materializes `[B, L, L, D]`.
- Benchmark peak memory and throughput at protein lengths 256, 512, 1024, and 2048.
- The fusion sequence contains compressed modality tokens and task queries only.

## G4: Biological Data Readiness

- Every source has URL/version/license/retrieval date/SHA256/row count.
- Record-level provenance is retained where upstream provides it.
- Experimental, curated, inferred, predicted, and simulated evidence remain distinct.
- Protein homology-cold and substrate scaffold-cold test sets are frozen.
- Duplicate and conflicting enzyme-substrate measurements are audited and resolved by
  a documented policy.
- Absolute-kinetics labels use canonical units and direct record-level provenance;
  activity and flux cannot enter through the absolute-kinetics adapter.
- A nonempty FBA vector has frozen ordered feature, model, solver, objective, and condition
  identities; a flux target has explicit simulated/experimental semantics and cannot leak
  through its FBA inputs.

## G5: Publication Readiness

- Compare against mean, sequence-only, molecule-only, late-fusion, and an established
  kinetic predictor baseline.
- Run modality, task-query, FBA-context, and uncertainty ablations.
- Report at least protein-homology-cold and scaffold-cold results.
- Use an external dataset that does not overlap model training data.
- Calibrate uncertainty and report confidence intervals across at least three seeds.
- Do not claim host-level benefit without enzyme-constrained FBA and independent
  experimental evidence.
- Do not describe Km, activity, or flux as trained unless governed labels for that exact
  task and the corresponding checkpoint/split identities are present.
- A new prospective ranking claim requires the separate pre-label readiness and checkpoint
  freeze, complete attempted-candidate rosters, development/closed-test family nonoverlap,
  blind prediction deposit, and one final campaign-level evaluation receipt.

## Current Gate State

| Gate | State |
|---|---|
| G1 code correctness | Passed for current scope: CPU tests, CUDA train/infer, resume, failure gates |
| G2 stable training | Passed for internal baseline: full train, resume, holdout metrics, atomic checkpoints |
| G3 scalability | Passed for current scope: sparse graph and 256-2048 residue GPU benchmark |
| G4 biological readiness | Partial: the frozen global audit has 192 accepted records, 25 global MMseqs families, and 54 primary/varied-substrate PubChem CIDs after the label-independent family cap; the substrate gate passes, but the 300-record and 30-family gates fail; recent NQO, Thermus CMP-kinase/pyrophosphatase, and phosphoglycerate-kinase sources add compliant exact-construct rows, while other discovery candidates remain blocked by exact-construct, substrate, saturation/endpoint, license, or homology gates; row-level UniKP citations remain unavailable |
| G5 publication readiness | Not ready: record and family readiness remain blocked; EnzEngDB and IMDH rankings fail and remain closed; the prospective ranking gate is implemented but has no frozen campaign evidence, blind prediction deposit, or one-time evaluation receipt; DLKcat overlap is not auditable; Morgan MLP nearly matches multimodal models; fusion/ablation evidence is post-hoc; validation-only grouped-CV calibration selects identity and provides conformal intervals but has no new independent external calibration result; governed absolute Km and activity training evidence is absent, while the flux smoke is explicitly simulated rather than experimental validation |

The software implementation itself is closed as `software-implementation-v1`. G4/G5
scientific gaps are explicitly `blocked_pending_external_data`, not incomplete software.
Their recovery conditions are frozen in `configs/external_data_blockers.json`. Broad
literature discovery is no longer the primary completion strategy.

No component is described as 90% complete until every mandatory item in its gate passes.
Implementation of interfaces does not complete the associated scientific validation gap;
the project does not claim that six scientific gaps are complete.
