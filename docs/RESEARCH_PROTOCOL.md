# Research Protocol

## Primary Question

Does task-query fusion of protein sequence, substrate graph, reaction metadata, and host
context improve enzyme kinetic prediction under cold-start evaluation compared with
single-modality and late-fusion baselines?

## Primary Endpoint

`log10(kcat)` RMSE on a frozen protein-homology-cold test set.

## Secondary Endpoints

- MAE and Pearson correlation;
- substrate scaffold-cold RMSE;
- predictive interval coverage;
- inference memory and throughput by protein length;
- future Km, absolute activity, and host flux tasks when governed labels and frozen splits
  exist;
- zero-shot within-campaign EnzEngDB fitness ranking, reported as Spearman, pairwise
  accuracy, and top-decile enrichment rather than cross-campaign RMSE.

## Dataset Policy

The local 17,010-row UniKP/DLKcat-derived corpus is a development corpus. It has
training overlap with existing predictors, lacks record-level accessions, and is not an
independent benchmark. Six non-positive values and multi-component SMILES are excluded
using the upstream preprocessing policy. Duplicate and conflicting values are reported.

Direct absolute-kinetics datasets may use normalized-kinetics schema v1. Each JSONL label
must be one of `log10_kcat`, `log10_km`, or `log10_kcat_per_km`, carry direct evidence and
nonempty record-level provenance, and resolve to `log10(s^-1)`, `log10(M)`, or
`log10(M^-1 s^-1)`. Supported positive linear source units are converted by the adapter.
Activity and flux labels are prohibited in this adapter. CLI training and evaluation must
declare `--data-format normalized-kinetics` and the requested tasks; each requested task
must have observations in train and validation before training begins.

The prospective temporal registry currently has 192 finalized accepted records and 54
globally unique primary/varied-substrate PubChem CIDs. The frozen global audit maps all 192 records
from 18 sources to 203 exact row/component sequences and groups them into 25 families under
MMseqs2 version `5d152c612b6ad2a56f657b7a02c127eceaea2a75`, 30% identity, 80% coverage,
and coverage mode 0. The label-independent 20-record family cap retains all 192 records.
The substrate gate passes; the 300-record and 30-family gates fail. The fourth
discovery round records 13 additional absent-by-annotation families, but none clears every
exact-construct, substrate, saturation, license, and homology requirement. Readiness remains
blocked pending eligible new data, and model predictions remain prohibited.

The first MMseqs2 development split uses 30% identity and 80% coverage. It is stored
with the source hash and clustering parameters. EnzEngDB v1 is a CC BY 4.0 external
sequence-function source. Its frozen transfer set removes all sequences with MMseqs2
hits to UniKP at the same 30% identity and 80% coverage thresholds. EnzEngDB fitness is
campaign-specific and cannot be interpreted as an absolute kinetic endpoint.

The EnzEngDB selection is capped by a label-independent hash before homology filtering.
Metrics are computed within each campaign and macro-averaged so large deep-mutational
scans cannot dominate. Future activity training must split whole campaign families and
must not evaluate on campaigns used for optimization or model selection.

The Lunzer IMDH landscape is a CC0 direct-kinetics mutation benchmark. Its 512 complete
six-site genotypes are reconstructed on UniProt P30125 and paired with NAD and NADP.
Observed ranking uses `ln(Km) + ln(kcat/Km)` within each cofactor. This removes dependence
on an unresolved constant source-unit offset. The frozen primary metrics are Spearman,
pairwise accuracy, and top-decile enrichment, macro-averaged over the two cofactors.
All 512 variants have MMseqs2 hits to UniKP at 30% identity and 80% coverage, so this is a
local mutation-sensitivity test, not a homology-cold benchmark.

## Planned Baselines

1. Train-set mean.
2. Protein amino-acid composition MLP.
3. Morgan fingerprint MLP.
4. Independent protein/molecule encoders with late concatenation.
5. Full task-query model.
6. An established kinetic predictor evaluated under a legally and scientifically valid
   non-overlapping protocol.

Implemented fixed-budget references additionally include amino-acid-composition ridge,
Morgan radius-2 2,048-bit ridge, and a late-concatenation model with the same encoders and
training budget as task-query fusion. A post-hoc paired comparison at seeds 7, 42, and
123 finds lower mean RMSE/MAE for task-query fusion and higher mean Pearson for late
concatenation. Paired RMSE direction is inconsistent. Because the internal test was
already observed before this protocol, these results are seed-stability evidence rather
than a new confirmatory test and do not establish either architecture as superior.

The composition and Morgan MLP baselines use train-only standardization, fixed hidden
dimensions 128 and 64, and seeds 7, 42, and 123. The Morgan MLP nearly matches the
multimodal RMSE, limiting claims about multimodal effect size. Seed-42 no-context,
shared-query, and global-mean-protein ablations are diagnostic only because the internal
test was already observed. FBA and evidence-weight ablations are not identifiable on the
current corpus; all FBA masks are false and all evidence weights are constant.

The uncertainty-objective ablation compares the existing heteroscedastic Gaussian NLL
with a mean-only fixed-variance MSE model at seeds 7, 42, and 123. Heteroscedastic training
has lower RMSE and MAE at every paired seed. A single standard-deviation scale is then fit
on validation only using the Gaussian NLL closed form and applied unchanged to test.
Calibration does not consistently improve test NLL or coverage, so uncertainty remains
uncalibrated. No further calibration method may be selected using this internal test.

Reusable calibration artifacts follow a validation-only workflow. `fit-calibration`
accepts an explicit identity, scalar, or affine-variance method, fits only labeled
validation rows, and writes the fit partition, validation count, exact checkpoint/task
identities, and an integrity SHA256. `predict --calibration-artifact` must reject an
identity mismatch and report raw and calibrated standard deviations separately. Test
labels are evaluation-only and must never fit the artifact. Saving or applying an artifact
does not establish calibration, and calibration does not turn an untrained task into a
trained one.

## FBA And Flux Contract

FBA context is optional. Every nonempty schema-v1 vector must carry ordered unique feature
IDs whose count exactly matches tensor width, plus nonempty model, solver, objective, and
condition identities. Empty context carries no FBA metadata. FBA values are simulated
context under the named assumptions, not experimental truth.

A supplied flux target must use the canonical `log10(mmol gDW^-1 h^-1)` unit and metadata
that identifies it as simulated or experimental, names the target reaction, and retains
nonempty provenance. The target reaction must not appear among the FBA input feature IDs.
Simulated and experimental flux labels remain distinct. No flux head may be described as
trained without governed labels satisfying this contract and a checkpoint that records
observed training and validation labels.

## Ablations

- Remove organism and EC context.
- Remove molecular graph.
- Replace protein chunks with global mean embedding.
- Replace task queries with one shared pooled vector.
- Remove FBA context.
- Remove evidence weighting.
- Remove heteroscedastic uncertainty.

## Claim Boundary

The model ranks computational candidates. It does not establish catalytic activity,
metabolic benefit, safety, or experimental validity. FBA demonstrates feasibility under
model assumptions only. Predictions from UniKP, DLKcat, DisorderFlow, or other models
are auxiliary evidence and never experimental ground truth.

Implemented adapters, heads, masks, calibration artifacts, and FBA/flux schemas are
engineering capabilities, not completed scientific validations. In particular, Km,
absolute activity, and flux cannot be claimed trained without governed labels for the
exact task. The current evidence does not support a claim that six scientific gaps are
complete.

The current EnzEngDB run uses a kcat-trained head as a zero-shot proxy. Its weak ranking
result cannot support candidate selection and is not described as external kcat or
activity validation.

The follow-up activity ranker uses campaign-local percentile labels and pairwise loss.
All families containing frozen-test campaigns are excluded before training. The resulting
development data are concentrated in one large training family, and three-seed frozen
test intervals overlap or underperform random expectations. No further architecture or
hyperparameter choices may be made using this frozen test; a new governed benchmark is
required for subsequent confirmatory claims.

The Lunzer selection was frozen before evaluating BioCandidateRanker or DLKcat. Both
predictors underperform random ranking expectations. No architecture, feature, training,
or hyperparameter decision may now use this landscape. DLKcat provides an established
predictor comparison under GPL-3.0-only, but its record-level training-source overlap is
not auditable and therefore does not satisfy the independent absolute-kinetics gate.

A future candidate-ranking claim must use the separate prospective gate in
`docs/PROSPECTIVE_RANKING_BENCHMARK.md`. EnzEngDB and IMDH cannot be reused for model
selection. Before any new campaign labels are available, the full attempted-candidate
roster, family nonoverlap against all development and closed tests, model checkpoints, and
campaign metrics must be frozen. Predictions are deposited blind and labels are joined
once for campaign-level reporting; inactive and censored candidates remain auditable.

### Acquisition Failures

An HTTP 404 is not retried during the same curation cycle. The endpoint is skipped and
recorded in `artifacts/external/temporal-absolute-kinetics/unfinished-acquisitions.json`
with its source impact and evidence path. Other transport failures are listed separately
and are not mislabeled as HTTP 404 responses. Scientific exclusions are not acquisition
failures and remain in their source-specific audit artifacts.

### Closed Software Baseline

The software baseline is frozen in `configs/software_implementation_release_v1.json` using
file-level SHA256 identities. Scientific work that requires unavailable external evidence is
marked `blocked_pending_external_data` in `configs/external_data_blockers.json`. Broad
literature discovery is not the primary workstream after this closure. Use targeted data
requests or experimental collaboration, then execute `docs/DATA_ARRIVAL_RUNBOOK.md`.

For an isolated acquisition, schema, or source issue, attempt one bounded correction. If it
remains unresolved, record its evidence and impact, skip the affected rows, and continue.
Fail the overall workflow only when accepted-data identity, frozen evidence identity, or
test correctness is affected.
