# Current Results

Date: 2026-07-27

## Software Release State

`software-implementation-v1` is closed at the file identities in
`configs/software_implementation_release_v1.json`. The normalized kinetics path, frozen
homology/global audit, grouped validation-only calibration, conformal intervals,
prospective-ranking governance, and simulated-flux smoke are implemented and tested.

The following scientific workstreams are `blocked_pending_external_data`:

- independent absolute-kinetics benchmark expansion;
- governed Km training;
- governed absolute-activity training;
- experimental-flux training;
- prospective candidate-ranking evaluation.

Broad literature discovery is no longer the primary workstream. Targeted requests,
experimental collaboration, and the post-arrival workflow are documented in
`docs/TARGETED_DATA_REQUESTS.md`, `docs/EXPERIMENTAL_COLLABORATION_SPEC.md`, and
`docs/DATA_ARRIVAL_RUNBOOK.md`.

The targeted evidence phase has now been executed and is recorded in
`artifacts/external/evidence-phase-receipt-2026-07-27.json`. It found no qualifying public
108-row/five-family absolute-kinetics batch, no independently custodied unobserved campaign,
and no currently governable Km or experimental-flux corpus. SABIO-RK exposes substantial
row-level Km metadata, but exact assay chains and permission for the intended training and
artifact use remain unresolved. The targeted yeast flux table lacks units, condition and
replicate metadata, uncertainty, methods, citations, and complete model mappings. No model
was trained from either source, and no closed test was reopened.

## Data Audit

- Raw rows: 17,010
- Accepted after upstream-compatible filtering: 16,838
- Rejected: 172
- Duplicate sequence-SMILES pairs: 395
- Conflicting duplicate measurements: 215
- Conflict policy for reported runs: median on the log10 scale
- Source SHA256: `13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719`
- Source size: 12,132,719 bytes

The corpus has corpus-level provenance only, redistribution permission is not verified,
and overlap with existing kinetic predictors prevents independent benchmark claims.

## Prospective Temporal Curation

The prospective protocol was frozen before curation at
`configs/temporal_absolute_kinetics_protocol.json`. It requires a primary article dated
2023 or later, CC0 or CC BY 4.0 measurement data, record-level provenance, zero exact
sequence-substrate overlap, and no UniKP MMseqs2 hit at 30% identity and 80% coverage.
No model predictions have been generated for this pool.

The strict registry audit at
`artifacts/external/temporal-absolute-kinetics/registry-readiness-audit.json` supersedes
earlier pool summaries. It reports 192 finalized records with status
`accepted_homology_cold_pool` and 54 globally unique substrates, where substrate identity
is the PubChem CID of the primary or varied, measured substrate. The frozen global audit maps all 192
records from 18 sources to 203 exact row/component sequences and groups them into 25 global
families. The label-independent 20-record family cap retains all 192 records and 54 registry
substrate identities. These counts
must not be reported as a final benchmark size.

CM/CDT homology is now complete for source `europepmc-PMC10520331`. The source contributes
17 saturation-qualified accepted records, comprising 10 isolated CM and 7 isolated CDT
measurements across 12 accepted exact sequences. MMseqs2 found zero development-corpus
hits. Frozen `linclust` grouped all 13 exact assayed constructs into seven source-level
clusters using MMseqs2 version `5d152c612b6ad2a56f657b7a02c127eceaea2a75`, 30% minimum
identity, 80% coverage, and coverage mode 0. The family-cluster SHA256 is
`19a029ee1e80a645c9253e9aa24a552ee8c45330249860b8e55321bf20363dc5`.

The frozen global family count is 25. The record gate is blocked at 192 < 300 and the
family gate is blocked at 25 < 30; the 50-substrate gate passes with 54 identities. The fourth- and
fifth-round machine-readable discovery registries screen 13 and at least 15 additional
absent-by-source families, respectively. No source in either round clears every exact
construct, substrate, saturation/endpoint, license, and homology gate. Fifth-round source
`biostudies-S-EPMC11742315` was fully acquired and curated: eight unique finite ShPepV
rows have exact `QKQ29470.1`/`A0A657M1L9` plus `ENLYFQG-His10` constructs and defined
PubChem substrates. They remain exclusion evidence because the ShPatB-DTNB coupling step
was not shown non-rate-limiting across the fitted series. The pool therefore remains
blocked, is not an external benchmark, and has not been scored by a model.

PMC12676663 (`10.1021/acs.biochem.5c00561`) contributes no accepted records. Two direct
reciprocal-saturation measurements were curated for phosphite varied at 1.0 mM NAD+ and
NAD+ varied at 20 mM phosphite. The exact 339-aa assay construct is the retained `GSH`
thrombin scar plus O69054 carrying all 17 reported substitutions; its calculated mass of
36,764.37 Da agrees with the reported 36,764.4 Da and observed 36,764.1 Da. Pinned WSL
MMseqs2 found 22 development-corpus hits, including a 98.8%-identity hit over 98.8% query
coverage, so both otherwise eligible `kcat` rows are preserved as homology exclusions.
The fragment-activation rate constants are not direct accepted endpoints. No model
predictions were generated.

BioStudies S-EPMC12781114 (`10.1021/acs.biochem.5c00559`) contributes two compliant direct
reaction-level `kcat` rows: 11 +/- 1 s-1 for NQO-WT and 5.4 +/- 0.2 s-1 for NQO-P78G. Both
are global ping-pong bi-bi fits whose `kcat` is defined at saturated NADH and coenzyme Q0;
the values are not duplicated into substrate-specific labels. The exact pET20b(+) products
are Q9I4V0 residues 1-328 followed directly by His6 (334 aa), with P78G at native residue
78 and no pelB signal or linker. PubChem mappings are CID 439153 (NADH) and CID 69068
(coenzyme Q0). Exact checks against all 17,010 development rows and pinned MMseqs2 searches
both found zero overlap/hits. The explicitly unreliable WT NADH Km/specificity values and
pre-steady-state half-reaction constants remain machine-readable exclusions. No model
predictions were generated.

PMC12894115 (`10.1007/s12010-025-05449-0`) contributes no accepted records. The
supplement discloses an unaccessioned 516-aa SnFDHal ORF, but not the complete His6-tagged
assay construct; the reported `NdeI/XhoI` cloning also conflicts with the SI `EcoRI`
reverse primer. All three Table 1 `kcat` values come from 30-minute coupled endpoint
assays, are mass-balance-incompatible with the stated 5.5 uM enzyme concentration, and
lack saturation evidence for Fre, FAD/reduced-flavin supply, NADH, chloride, and oxygen.
Pinned WSL MMseqs2 found five development-corpus hits for the disclosed ORF, with a best
hit of 39.0% identity over 99.8% query coverage. The rows are retained only as audited
exclusions; no model predictions were generated.

## Governed Task Interfaces

The strict normalized-kinetics v1 JSONL adapter accepts direct experimental `kcat`, `Km`,
and `kcat/Km` labels with per-label provenance. Supported source units are normalized to
`log10(s^-1)`, `log10(M)`, and `log10(M^-1 s^-1)`. The `train` and `evaluate` CLI commands
accept `--data-format normalized-kinetics` and explicit subsets of `log10_kcat`,
`log10_km`, and `log10_kcat_per_km`; training fails when a requested task has no train or
validation observations. The adapter explicitly rejects activity and flux labels so that
unresolved endpoints cannot enter the absolute-kinetics path.

The version-2 calibration workflow performs five-fold model selection across identity,
scalar, and affine-variance Gaussian calibrators using frozen validation MMseqs groups only.
Gaussian NLL is primary, CRPS is secondary, and deterministic ties prefer the simpler model.
For the current `log10_kcat` checkpoint, identity won with validation NLL `1.67734` and CRPS
`0.72184`; scalar and affine variance were worse, so no rescaling is forced. The checkpoint-
bound artifact at `artifacts/uncertainty-calibration-v2/grouped-cv.json` also freezes finite-
sample normalized and unnormalized 90% and 95% conformal intervals. Tests prove that test
rows are not requested and MMseqs groups do not cross folds. This validates the calibration
procedure and its current identity selection; it does not create an independent external
calibration claim.

Nonempty FBA context follows schema v1: an ordered, unique feature-ID list must exactly
match vector width and include model, solver, objective, and condition identities. Flux
labels require an explicit `simulated` or `experimental` type, target reaction, canonical
`log10(mmol gDW^-1 h^-1)` unit, and nonempty provenance. The target reaction is forbidden
from the FBA feature IDs to prevent direct target leakage. Current UniKP FBA masks are all
false. No governed absolute Km, flux, or absolute-activity training result is reported;
the EnzEngDB campaign-percentile ranker is only a relative activity-ranking task.

## Frozen Splits

### Protein Homology Cold

- MMseqs2 version: `5d152c612b6ad2a56f657b7a02c127eceaea2a75`
- Identity threshold: 30%
- Coverage threshold: 80%
- Clusters: 2,204
- Raw accepted rows: 13,470 train / 1,684 validation / 1,684 test
- After median aggregation in the full run: 13,157 train / 1,640 validation / 1,646 test

### Scaffold Cold

- Rows: 13,470 train / 1,684 validation / 1,684 test
- Cross-split scaffolds: 0
- Cross-split proteins: 1,580

### Double Cold Stress Split

- Rows: 16,268 train / 285 validation / 285 test
- Cross-split proteins: 0
- Cross-split scaffolds: 0
- Limitation: a large bipartite component makes the split severely imbalanced.

## Full Model

Architecture: 64 hidden dimensions, two local protein layers, two sparse molecular
message-passing layers, one task-query fusion layer, 536,520 parameters.

Best checkpoint was selected on validation loss after epoch 1 of a three-epoch internal
baseline run.

| Model | Test RMSE | Test MAE | Test Pearson |
|---|---:|---:|---:|
| Train mean | 1.5101 | 1.1795 | 0.0000 |
| Amino-acid composition ridge | 1.5262 | 1.1701 | 0.0891 |
| Morgan radius-2 2,048-bit ridge | 1.6883 | 1.2813 | 0.1717 |
| Protein only | 1.5417 | 1.1821 | 0.1222 |
| Molecule only | 1.5370 | 1.1936 | 0.1085 |
| Task-query multimodal | 1.4638 | **1.1222** | 0.2731 |
| Late-concatenation multimodal | **1.4556** | 1.1287 | **0.3061** |

The seed-42 late-concatenation model used the same split, conflict policy, three-epoch
budget, encoder dimensions, and optimizer settings as the task-query run. It has 485,960
parameters versus 536,520. Its small single-seed RMSE and Pearson advantage did not
establish that late fusion was superior, so a fixed three-seed robustness comparison was
run for the two already evaluated architectures.

### Three-Seed Fusion Stability

The post-hoc robustness protocol was frozen before the six runs at seeds 7, 42, and 123.
This internal test had already been observed, so the analysis measures seed stability and
is not a new confirmatory architecture-selection experiment.

| Fusion mode | RMSE mean +/- SD | MAE mean +/- SD | Pearson mean +/- SD |
|---|---:|---:|---:|
| Task query | **1.4699 +/- 0.0064** | **1.1268 +/- 0.0063** | 0.2471 +/- 0.0282 |
| Late concatenation | 1.4801 +/- 0.0230 | 1.1396 +/- 0.0114 | **0.2821 +/- 0.0278** |

Paired `late_concat - task_query` differences were `+0.0101 +/- 0.0168` RMSE,
`+0.0128 +/- 0.0056` MAE, and `+0.0350 +/- 0.0357` Pearson. Late concatenation had lower
RMSE only for seed 42, while task query had lower RMSE for seeds 7 and 123. Late
concatenation had higher Pearson for all three seeds, although the seed-123 difference
was only 0.0003. Neither fusion mode dominates across endpoints, and the current evidence
does not support a robust task-query advantage.

The two ridge baselines used train-only standardized features and selected
`alpha` from `{0.1, 1, 10, 100}` by validation RMSE before test evaluation. Both selected
the grid boundary at 100. The grid must not be expanded based on the observed test result;
these are fixed-budget classical references rather than fully optimized MLP baselines.

### Feature MLP Baselines

The planned amino-acid-composition and Morgan MLP baselines were run at seeds 7, 42,
and 123 under a protocol frozen before training. Both use train-only standardization, a
`128 -> 64` MLP, dropout 0.1, AdamW, 40 epochs, and lowest validation RMSE selection.
Because the internal test had already been observed, these runs complete the baseline
matrix but are not new confirmatory evidence.

| Feature MLP | RMSE mean +/- SD | MAE mean +/- SD | Pearson mean +/- SD |
|---|---:|---:|---:|
| Amino-acid composition | 1.5777 +/- 0.0121 | 1.2185 +/- 0.0104 | 0.0428 +/- 0.0064 |
| Morgan radius-2 2,048-bit | **1.4756 +/- 0.0024** | **1.1469 +/- 0.0022** | **0.2566 +/- 0.0188** |

The Morgan MLP is a strong single-modality reference. Its mean RMSE is only 0.0057 worse
than the three-seed task-query mean and 0.0045 better than the late-concatenation mean.
The current data therefore do not demonstrate a practically large multimodal advantage.

### Internal Architecture Ablations

Three informative seed-42 ablations used the same task-query training budget. Positive
RMSE deltas indicate degradation relative to task-query seed 42 (`1.4638` RMSE):

| Ablation | Test RMSE | Delta RMSE | Test Pearson |
|---|---:|---:|---:|
| Remove context | 1.5046 | +0.0409 | 0.1908 |
| Shared task query | 1.4911 | +0.0273 | 0.2520 |
| Global-mean protein encoder | 1.4910 | +0.0272 | 0.2278 |

These single-seed, post-hoc diagnostics are consistent with contributions from context,
task-specific queries, and the chunk Transformer, but do not establish causal or
confirmatory superiority. FBA removal is not identifiable because all UniKP FBA masks are
false. Evidence-weight removal is also not identifiable because every record has the same
`unknown` tier and constant weight.

### Uncertainty Objective And Calibration

The heteroscedastic Gaussian NLL objective was compared with a fixed-variance mean-only
head trained by masked MSE at the same three seeds. Fixed-variance coverage is not
interpreted because its unit standard deviation is only an interface placeholder.

| Objective | RMSE mean +/- SD | MAE mean +/- SD | Pearson mean +/- SD |
|---|---:|---:|---:|
| Heteroscedastic Gaussian NLL | **1.4699 +/- 0.0064** | **1.1268 +/- 0.0063** | 0.2471 +/- 0.0282 |
| Fixed-variance MSE | 1.5023 +/- 0.0042 | 1.1524 +/- 0.0161 | **0.2586 +/- 0.0038** |

Paired `fixed - heteroscedastic` RMSE was positive for every seed, with mean
`+0.0324 +/- 0.0053`; MAE increased by `+0.0256 +/- 0.0195`. The heteroscedastic objective
therefore provides a stable point-error benefit under this fixed budget, although it does
not improve mean Pearson.

Each heteroscedastic checkpoint then received one validation-only standard-deviation
scale fitted by closed-form Gaussian NLL minimization. Mean scale was
`1.0017 +/- 0.0514`. Held-out test results were:

| Metric | Before calibration | After calibration | Nominal |
|---|---:|---:|---:|
| 1-sigma coverage | 0.6614 | 0.6594 | 0.6827 |
| 2-sigma coverage | 0.9176 | 0.9188 | 0.9545 |
| Gaussian NLL | 1.8247 | 1.8236 | not applicable |

NLL change was inconsistent by seed (`+0.0130`, `+0.0018`, `-0.0180`) and mean coverage
did not move materially toward nominal values. Validation-only scalar calibration is thus
a negative result and does not establish calibrated uncertainty. No alternative
calibration method is selected using this observed test.

Full-model uncertainty coverage was 0.6659 within one predicted standard deviation and
0.9228 within two. These values are not yet calibrated to nominal 68%/95% targets.

## Scaling

RTX 5060 Laptop GPU, batch size 4, hidden size 64, median of 10 warmed iterations:

| Protein length | Peak allocated memory |
|---:|---:|
| 256 | 20.7 MB |
| 512 | 31.1 MB |
| 1024 | 49.3 MB |
| 2048 | 87.4 MB |

All outputs were finite. Memory growth is consistent with chunked local attention and
sparse molecular message passing rather than full-sequence quadratic attention.

## External EnzEngDB Transfer

EnzEngDB v1 was retrieved from Zenodo DOI `10.5281/zenodo.17310823` under CC BY 4.0.
The official `Data.zip` identity is:

- Size: 30,050,527 bytes
- SHA256: `8013ad81586db2187162aada0709c1cabc7e7d69f03dd5c199776aaf000dd6ea`
- Experiment CSV rows: 462,092
- Strictly accepted rows: 245,945 across 160 campaigns

Rows with non-canonical protein sequences, missing/non-finite fitness, non-concrete
wildcard reactions, or unusable reaction SMILES were rejected. Selection retained at
most 2,000 rows per campaign by lowest SHA256 of candidate ID, independent of fitness.

MMseqs2 search against all UniKP development sequences used 30% identity, 80% coverage,
and coverage mode 0. Homology hits were removed, then the minimum 20 rows per campaign
gate was reapplied. The frozen result contains 6,423 rows across 51 campaigns, with zero
exact sequence overlap and zero exact sequence-substrate overlap. Its selection SHA256
is `8e76932fc76b884ea8cf127934683b43b9e896b057b817488a9dc5533cf1fc4b`.

The internal kcat checkpoint was evaluated as a zero-shot proxy for within-campaign
fitness ranking:

| Metric, macro over campaigns | Result | Random expectation |
|---|---:|---:|
| Spearman | 0.0710 | 0.0000 |
| Pairwise accuracy | 0.5237 | 0.5000 |
| Top-10% enrichment | 0.5943 | 1.0000 |
| Top-10% recall | 0.0648 | approximately 0.1000 |

This is a negative/weak transfer result. The current checkpoint is not suitable for
ranking EnzEngDB engineering variants. EnzEngDB fitness is not kcat, and the activity
head had no training labels, so this is neither external kcat validation nor activity
validation.

### Campaign-Aware Activity Ranking

A dedicated one-task ranker was initialized from the kcat encoder while its ranking head
was initialized from scratch. Fitness values were converted to average-tie percentile
ranks within each campaign. Training batches contained one campaign only and optimized
same-campaign pairwise logistic loss.

Campaign representatives were clustered with MMseqs2 at 30% identity and 80% coverage.
All 13 families containing any frozen-test campaign were excluded. After the 20-row
campaign gate and 500-row cap, the family-cold development split contained:

- Train: 7,345 records, 24 campaigns, 1 large family
- Validation: 441 records, 12 campaigns, 8 families
- Frozen test: 6,423 records, 51 campaigns

The single-family training concentration is an important limitation of the available
data. Three full-model seeds were selected independently by lowest family-cold validation
loss and evaluated on the same frozen test selection:

| Metric | Three-seed mean +/- SD | Campaign bootstrap 95% CI | Random |
|---|---:|---:|---:|
| Spearman | -0.0256 +/- 0.0795 | [-0.0549, 0.0033] | 0.0000 |
| Pairwise accuracy | 0.4910 +/- 0.0275 | [0.4808, 0.5014] | 0.5000 |
| Top-10% enrichment | 0.5973 +/- 0.1363 | [0.4443, 0.7585] | 1.0000 |
| Top-10% recall | 0.0665 +/- 0.0157 | [0.0491, 0.0845] | about 0.1000 |

Seed 7 diagnostic comparisons on the frozen test were:

| Model | Spearman | Pairwise accuracy | Top-10% enrichment |
|---|---:|---:|---:|
| Unmodified kcat proxy | 0.0710 | 0.5247 | 0.5859 |
| Full activity ranker | 0.0647 | 0.5223 | 0.4465 |
| Protein-only ranker | 0.0645 | 0.5202 | 0.6731 |
| Molecule-only ranker | -0.0564 | 0.4816 | 0.8054 |

No trained ranker exceeds random top-decile enrichment, and the three-seed confidence
intervals do not support generalization. Activity ranking therefore remains a negative
result and must not be used for candidate prioritization.

## Frozen IMDH Mutation Landscape

The CC0 Dryad dataset `10.5061/dryad.7nd70`, associated with article
`10.1126/science.1115649`, contains a complete 512-genotype landscape at E. coli
3-isopropylmalate dehydrogenase positions 236, 289, 290, 296, 337, and 341. Each genotype
has fitted `ln(Km)` and `ln(kcat/Km)` values for NAD and NADP, yielding 1,024 derived
`ln(kcat)` observations. Because the deposited file does not resolve the source-unit
constant, only within-cofactor ranking metrics were frozen; absolute RMSE is not reported.

- Source size: 38,761 bytes
- Source SHA256: `7bfd71d78235ade06a9e117b8229ad99b1fe7adceb0a5ed72d708c3e925c2a69`
- Exact UniKP sequence, substrate, and sequence-substrate overlap: 0 / 0 / 0 records
- MMseqs2 hits at 30% identity and 80% coverage: 512 / 512 unique sequences
- Best-hit identity range: 67.7% to 69.4%

The selection and metrics were frozen before either predictor was evaluated. Results are:

| Predictor | Spearman | Pairwise accuracy | Top-10% enrichment |
|---|---:|---:|---:|
| BioCandidateRanker full kcat checkpoint | -0.1659 | 0.4452 | 0.4734 |
| DLKcat `upstream-7c15d0d4a7ac` | -0.2586 | 0.4168 | 0.2840 |
| Random expectation | 0.0000 | 0.5000 | 1.0000 |

BioCandidateRanker scored NAD at Spearman -0.0538 and NADP at -0.2779. DLKcat scored
NAD at -0.2159 and NADP at -0.3012. Neither model recovers the experimental mutation
landscape. DLKcat is a legally runnable GPL-3.0-only comparison, but its record-level
training overlap cannot be audited. This is direct kinetic evidence from a familiar
homologous family, not a homology-cold external validation. The landscape is now closed
to model selection and hyperparameter tuning.

## Interpretation

The implementation is operational and the multimodal signal survives the internal
homology-cold test, but the gain is modest. Three-seed results show lower average RMSE and
MAE for task-query fusion but higher Pearson for late concatenation, with inconsistent
paired RMSE direction. A Morgan MLP nearly matches both multimodal models. The existing
kcat model fails both the
homology-cold EnzEngDB cross-endpoint ranking test and the homologous direct-kinetics IMDH
mutation landscape. DLKcat also fails the latter, so the result is not evidence that the
new model is publication-ready. The next scientific requirement remains a genuinely
independent absolute-kinetics benchmark with auditable training-source overlap, not
further tuning against either frozen external result.
