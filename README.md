# BioCandidateRanker

Host-aware multimodal enzyme candidate scoring research system.

The project is deliberately isolated from the source projects. External datasets and
tools are opened read-only through explicit paths; generated manifests, splits,
checkpoints, and reports stay under this repository.

## Software Implementation Release

The current software implementation is formally closed as `software-implementation-v3`.
Its file-level SHA256 freeze and claim boundary are in
`configs/software_implementation_release_v3.json`; v1 and v2 remain as historical records
under `configs/software_implementation_release_v*.json`. Verify the closed baseline with:

```powershell
python scripts\verify_software_release.py --manifest configs\software_implementation_release_v3.json
```

Absolute-kinetics benchmark expansion, governed Km training, governed absolute-activity
training, experimental-flux training, and prospective candidate-ranking evaluation are
`blocked_pending_external_data`. The machine-readable states and recovery conditions are in
`configs/external_data_blockers.json`. Broad literature discovery is no longer the primary
workstream. New work proceeds through `docs/TARGETED_DATA_REQUESTS.md` and
`docs/EXPERIMENTAL_COLLABORATION_SPEC.md`; once data arrive, follow
`docs/DATA_ARRIVAL_RUNBOOK.md`.

An isolated source or endpoint failure receives one bounded correction attempt. If it
remains unresolved, it is recorded and skipped rather than blocking unrelated work. Errors
that alter accepted-data identity, frozen evidence, or test correctness still fail closed.

## Current Scope

The first trained task predicts `log10(kcat)` from:

- protein sequence, encoded by local chunked attention;
- substrate molecular graph, encoded by sparse message passing;
- organism, EC, enzyme type, and reaction context;
- optional FBA context features;
- task-specific query tokens with predictive uncertainty.

Nonempty FBA context uses the strict v1 feature contract: the ordered feature IDs must
exactly match the configured tensor width and carry model, solver, objective, and condition
identities. A supplied flux target requires metadata declaring it as simulated or
experimental and cannot name a reaction included in its FBA inputs. Empty context remains
supported. These semantics reserve a governed task interface; they do not imply that an
experimental flux model has been trained.

### Simulated Flux Smoke

An isolated reproducibility pipeline can generate a small **simulated** biomass-flux
task from a read-only Yeast-MetaTwin YAML asset and run one CPU training epoch. Inputs are
only declared glucose and oxygen medium bounds; the solved target reaction is excluded
from FBA features. This is an executable engineering smoke test under model assumptions,
not experimental flux and not publication validation.

```powershell
$env:PYTHONPATH="src"
python scripts\run_simulated_flux_smoke.py `
  --model "C:\biological\Metabolic model prediction\Yeast-MetaTwin\Data\model\Yeast-MetaTwin.yml" `
  --output-dir artifacts\simulated-flux-smoke
```

The output contains normalized `log10(mmol gDW^-1 h^-1)` labels, a frozen split,
`manifest.json`, one CPU smoke checkpoint, and metrics. The manifest and checkpoint bind
SHA256 identities for the model bytes, solver specification, objective specification,
medium protocol, dataset, and split. Solved fluxes are never input features. Missing
COBRA/GLPK support, a missing or incompatible model, failed/nonpositive solutions, source
mutation, or contract failure removes partial outputs and writes only `blocker.json`.

The architecture also supports `log10(Km)`, `log10(kcat/Km)`, activity, and flux task
interfaces. Missing labels are masked and do not contribute to training. A head being
present is not evidence that it was trained: Km, absolute activity, and flux may be
described as trained only when the checkpoint records governed labels for that task. The
campaign-percentile EnzEngDB ranker is a separate relative-ranking task, not a trained
absolute-activity head.

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
```

## Data Audit

Audit the external UniKP development corpus without modifying it:

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.cli audit-data `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --output artifacts\unikp_audit.json
```

The generated manifest records SHA256, byte size, row count, accepted/rejected rows,
duplicates, conflicts, and scientific limitations. Training fails if the source no
longer matches this identity.

### Normalized Absolute Kinetics

The strict v1 JSONL adapter accepts sparse direct experimental labels for `log10_kcat`,
`log10_km`, and `log10_kcat_per_km`. It converts positive linear `kcat`, `Km`, and
`kcat/Km` values from supported units to `log10(s^-1)`, `log10(M)`, and
`log10(M^-1 s^-1)`, respectively. Every label requires `evidence: "direct"` and a nonempty
JSON provenance object; activity and flux labels are rejected by this adapter.

Train or evaluate selected normalized tasks through the same manifest- and split-bound
CLI:

```powershell
python -m biocandidate.cli train `
  --data artifacts\normalized-kinetics.jsonl `
  --data-format normalized-kinetics `
  --tasks log10_kcat log10_km log10_kcat_per_km `
  --manifest artifacts\normalized-kinetics-manifest.json `
  --split-manifest artifacts\normalized-kinetics-split.json `
  --output-dir artifacts\normalized-training
```

Training preflight rejects any requested task with zero observations in train or
validation. This interface does not establish that a Km checkpoint has been trained; that
claim requires a governed labeled dataset, frozen split, and corresponding checkpoint.

## Prospective Temporal Benchmark

The frozen curation policy is documented in `docs/TEMPORAL_BENCHMARK_PROTOCOL.md` and
`configs/temporal_absolute_kinetics_protocol.json`. Candidate metadata, acquisition
hashes, and exclusion evidence are under
`artifacts/external/temporal-absolute-kinetics`.

The strict registry currently sees 192 accepted records and 54 primary/varied-substrate PubChem
CIDs across finalized sources before exact global sequence mapping. The global audit maps
all 192 records from 18 sources to 203 exact row/component sequences; they form 25 global
families and contain 54 substrates. No family exceeds 20 records, so the label-independent
cap retains all 192 mapped records. The eight held PMC12329711 controls resolve to three
homology exclusions, one nonfinite exclusion, and four exact-tag-mapping blockers; they add
no accepted records. MMseqs2 uses frozen version
`5d152c612b6ad2a56f657b7a02c127eceaea2a75` at 30% identity, 80% coverage, and coverage
mode 0.

The fourth and fifth discovery rounds are frozen in machine-readable registries under
`artifacts/external/temporal-absolute-kinetics`. The fifth round adds 10 previously
unscreened BioStudies-backed sources spanning at least 15 families. Its NQO source contributes
two direct global-fit `kcat` rows for exact Q9I4V0-His6 WT/P78G products after zero exact
sequence-substrate overlaps and zero frozen-MMseqs development hits. Its fully acquired ShPepV source has exact tagged constructs, defined PubChem
substrates, and pinned-MMseqs evidence, but remains excluded because the coupled ShPatB-DTNB
assay was not demonstrated non-rate-limiting across the fitted substrate series.
The ENPP1 source resolves the glycosylated construct-770 mature fusion, both direct Table 2
`kcat` rows, two-step cGAMP endpoint semantics, PubChem structures, and saturation. It adds
no accepted rows because pinned MMseqs2 finds an 80.0%-identity development hit covering
100% of the exact ENPP1 catalytic component; zero exact sequence-substrate overlaps were found.

The separate PMC12676663 audit resolves the exact `GSH`-scarred 17X-PsPTDH construct and
two direct reciprocal-saturation `kcat` rows, but frozen MMseqs2 finds 22 UniKP development
hits. The rows and complete hit evidence are retained as homology exclusions under
`artifacts/external/temporal-absolute-kinetics/europepmc-PMC12676663`; fragment-activation
rate constants are excluded as non-`kcat` endpoints.

This remains a curation pool, not a benchmark. The exact mapped subset passes the
50-substrate gate with 54 identities but fails the 300-record and 30-family gates; unresolved rows fail closed.
No model predictions may be generated until every gate passes and the final record list and
checkpoint are frozen.

### Prospective Candidate Ranking

A separate fail-closed readiness gate for a future prospective campaign benchmark is
documented in `docs/PROSPECTIVE_RANKING_BENCHMARK.md`, with a data-free example protocol at
`configs/prospective_candidate_ranking_protocol.example.json`. It requires complete
candidate rosters including inactive and censored entries, family nonoverlap with all
development and closed tests, pre-label checkpoint identities, label-independent
selection, a blind prediction deposit, and one final campaign-level evaluation. It does
not reopen or extend the closed EnzEngDB or IMDH tests and creates no benchmark data or
scores unless a supplied evidence package passes readiness.

## Development Training

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.cli train `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --split-strategy protein_cold `
  --conflict-policy keep `
  --output-dir artifacts\training
```

`--conflict-policy keep` is explicitly development-only. The source has no frozen
independent benchmark and exact protein-cold splitting is not homology-cold splitting.
Neither setting is acceptable for a publication claim.

## Homology-Cold Split

MMseqs2 can run through WSL while all generated files remain in this project:

```powershell
python -m biocandidate.cli build-homology-split `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --mmseqs "wsl:Ubuntu-24.04:/home/<user>/.local/bin/mmseqs" `
  --output-dir artifacts\homology
```

## Global Temporal Family Audit

The temporal readiness audit consumes only sources finalized by `temporal_registry`, maps
accepted rows through explicit source-schema columns to exact construct FASTA IDs, and
fails closed for a missing, ambiguous, or unhashed mapping. It globally clusters every
mapped component with the frozen MMseqs2 30% identity, 80% coverage, coverage-mode 0
settings. Multicomponent records share a family when any component clusters connect them;
the 20-record family cap then follows source ID and CSV row order without using labels.

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.global_temporal_audit `
  --mmseqs "mmseqs" `
  --output-dir artifacts/external/temporal-global-family-audit
```

The output directory contains `global-accepted-components.fasta`, the exact row/component
mapping in `global-accepted-mapping.csv`, MMseqs cluster evidence, and
`readiness-audit.json`. Unresolved sources are reported but contribute no records,
families, or substrates. This command performs no model scoring.

Train against the frozen result with `--split-manifest
artifacts\homology\homology_split.json`. The recorded development split contains 2,204
MMseqs2 clusters and 13,470/1,684/1,684 train/validation/test rows.

## Scaling Benchmark

```powershell
$env:PYTHONPATH="src"
python scripts\benchmark_scaling.py --device cuda --output artifacts\scaling.json
```

On the current RTX 5060 Laptop GPU, batch size 4 and 64 hidden dimensions used about
20.7, 31.1, 49.3, and 87.4 MB peak allocated memory at lengths 256, 512, 1024, and
2048 respectively. This is an engineering measurement, not a cross-model benchmark.

## Candidate Scoring

```powershell
python -m biocandidate.cli predict `
  --checkpoint artifacts\training\best.pt `
  --input configs\candidate_example.json `
  --output artifacts\candidate_predictions.json
```

The output retains the checkpoint training-data manifest, input file identity,
task-specific means, uncertainty estimates, and a non-validation warning.

Fit a reusable uncertainty calibration artifact on validation labels only, then bind it to
prediction explicitly:

```powershell
python -m biocandidate.cli fit-calibration `
  --checkpoint artifacts\full-homology-baseline\best.pt `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --split-manifest artifacts\homology-final\homology_split.json `
  --task log10_kcat `
  --method scalar `
  --output artifacts\calibration\log10-kcat-scalar.json

python -m biocandidate.cli predict `
  --checkpoint artifacts\full-homology-baseline\best.pt `
  --calibration-artifact artifacts\calibration\log10-kcat-scalar.json `
  --input configs\candidate_example.json `
  --output artifacts\candidate_predictions.json
```

`fit-calibration` accepts explicit `identity`, `scalar`, or `affine-variance` methods and
uses only the frozen validation partition. The artifact records `fit_partition:
"validation"`, validation count, task and checkpoint identities, and an integrity SHA256.
Prediction rejects a checkpoint identity mismatch and reports raw and calibrated standard
deviations separately. A calibration artifact changes uncertainty reporting only; it does
not train a task or justify a calibration claim without held-out evidence. The current
validation-only scalar result remains negative and uncertainty remains uncalibrated.

Training supports controlled ablations with `--disable-protein`, `--disable-molecule`,
`--disable-context`, and `--shared-task-query`. A run with every modality disabled is
rejected.

Evaluate a frozen holdout without changing the checkpoint:

```powershell
python -m biocandidate.cli evaluate `
  --checkpoint artifacts\full-homology-baseline\best.pt `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --split-manifest artifacts\homology-final\homology_split.json `
  --partition test `
  --output artifacts\full-homology-baseline\test_metrics.json
```

Evaluation refuses a source or split identity that differs from the checkpoint.

Run the same-encoder late-concatenation baseline by adding
`--fusion-mode late_concat` to `train`, then evaluate its selected checkpoint with the
same command above. Fixed-budget classical references are produced with:

```powershell
python -m biocandidate.cli evaluate-classical-baselines `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --split-manifest artifacts\homology-final\homology_split.json `
  --output artifacts\classical-homology-baselines.json
```

This optional command requires `pip install -e ".[baseline]"`. Alpha is selected only on
validation from the fixed `{0.1, 1, 10, 100}` grid.

The paired three-seed fusion protocol and summary are reproduced with:

```powershell
python -m biocandidate.cli summarize-fusion-runs `
  --protocol configs\fusion_three_seed_protocol.json `
  --runs-root artifacts\fusion-three-seed `
  --output artifacts\fusion-three-seed\summary.json
```

Task-query has lower mean RMSE/MAE, while late concatenation has higher mean Pearson;
paired RMSE direction varies by seed. This post-hoc internal analysis does not establish
either fusion mechanism as superior.

Feature MLP runs are governed by `configs/feature_mlp_three_seed_protocol.json`:

```powershell
python -m biocandidate.cli train-feature-mlp `
  --protocol configs\feature_mlp_three_seed_protocol.json `
  --feature morgan_r2_2048 --seed 42 `
  --data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --manifest artifacts\unikp_audit.json `
  --split-manifest artifacts\homology-final\homology_split.json `
  --output-dir artifacts\feature-mlp-three-seed\morgan-r2-2048-seed42

python -m biocandidate.cli summarize-feature-mlp-runs `
  --protocol configs\feature_mlp_three_seed_protocol.json `
  --runs-root artifacts\feature-mlp-three-seed `
  --output artifacts\feature-mlp-three-seed\summary.json
```

The Morgan MLP reaches `1.4756 +/- 0.0024` test RMSE, close to both multimodal models.
Informative internal ablations are summarized in
`artifacts/internal-ablations/summary.json`; FBA and evidence-weight removals are
non-identifiable on UniKP because their inputs/weights are constant.

The three-seed uncertainty objective and validation-only calibration results are stored
in `artifacts/uncertainty-ablation/summary.json` and
`artifacts/uncertainty-calibration/summary.json`. Heteroscedastic Gaussian NLL improves
point RMSE over fixed-variance MSE, but one-parameter validation scaling does not
consistently improve held-out NLL or move coverage to nominal values. Uncertainty outputs
must therefore still be described as uncalibrated.

## External Campaign Benchmark

EnzEngDB v1 is retrieved from Zenodo DOI `10.5281/zenodo.17310823` under CC BY 4.0.
The official archive is verified as 30,050,527 bytes with SHA256
`8013ad81586db2187162aada0709c1cabc7e7d69f03dd5c199776aaf000dd6ea`.

Freeze a label-independent, UniKP-homology-cold campaign selection:

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.cli audit-enzengdb `
  --archive artifacts\external\enzengdb-v1\Data.zip `
  --experiments artifacts\external\enzengdb-v1\extracted\Data\experiments `
  --unikp-data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --unikp-manifest configs\unikp_source_manifest.json `
  --mmseqs "wsl:Ubuntu-24.04:/home/<user>/.local/bin/mmseqs" `
  --exclude-homology-hits `
  --selection-output artifacts\external\enzengdb-v1\selection-homology-cold.json `
  --output artifacts\external\enzengdb-v1\audit-homology-cold.json
```

Run the frozen zero-shot transfer test:

```powershell
python -m biocandidate.cli evaluate-enzengdb `
  --checkpoint artifacts\full-homology-baseline\best.pt `
  --archive artifacts\external\enzengdb-v1\Data.zip `
  --experiments artifacts\external\enzengdb-v1\extracted\Data\experiments `
  --selection-manifest artifacts\external\enzengdb-v1\selection-homology-cold.json `
  --task log10_kcat `
  --output artifacts\external\enzengdb-v1\zero-shot-kcat-proxy.json
```

EnzEngDB `fitness_value` has campaign-specific meaning. Metrics are calculated only
within campaigns and macro-averaged. This command tests cross-endpoint ranking transfer;
it is not external kcat validation and does not validate the untrained activity head.

Campaign-aware activity ranking is available as a separate protocol. The split excludes
every family containing a frozen-test campaign, then assigns complete MMseqs2 campaign
families to train or validation:

```powershell
python -m biocandidate.cli build-enzengdb-rank-split `
  --archive artifacts\external\enzengdb-v1\Data.zip `
  --experiments artifacts\external\enzengdb-v1\extracted\Data\experiments `
  --test-selection artifacts\external\enzengdb-v1\selection-homology-cold.json `
  --mmseqs "wsl:Ubuntu-24.04:/home/<user>/.local/bin/mmseqs" `
  --output-dir artifacts\external\enzengdb-v1\rank-split

python -m biocandidate.cli train-enzengdb-ranker `
  --archive artifacts\external\enzengdb-v1\Data.zip `
  --experiments artifacts\external\enzengdb-v1\extracted\Data\experiments `
  --split-manifest artifacts\external\enzengdb-v1\rank-split\rank_split.json `
  --initialize-from artifacts\full-homology-baseline\best.pt `
  --output-dir artifacts\external\enzengdb-v1\ranker-seed42
```

Labels are average-tie percentiles within each campaign. Batches never mix campaigns,
and optimization uses same-campaign pairwise logistic loss. The three-seed result is
summarized in `artifacts/external/enzengdb-v1/ranker-three-seed-summary.json`.

## Direct-Kinetics Mutation Benchmark

The CC0 Lunzer IMDH landscape (`10.5061/dryad.7nd70`) contains 512 complete six-site
genotypes measured against NAD and NADP. Freeze the full selection and exact-overlap audit
before evaluation:

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.cli audit-lunzer `
  --data "artifacts\external\absolute-kinetics-screen\dryad-4964723\Lunzer Curated Biochemistry.xls" `
  --unikp-data "<external-path>\Kcat_combination_0918_wildtype_mutant.json" `
  --unikp-manifest configs\unikp_source_manifest.json `
  --output artifacts\external\absolute-kinetics-screen\dryad-4964723\selection.json

python -m biocandidate.cli evaluate-lunzer `
  --checkpoint artifacts\full-homology-baseline\best.pt `
  --data "artifacts\external\absolute-kinetics-screen\dryad-4964723\Lunzer Curated Biochemistry.xls" `
  --selection-manifest artifacts\external\absolute-kinetics-screen\dryad-4964723\selection.json `
  --output artifacts\external\absolute-kinetics-screen\dryad-4964723\biocandidate-full-evaluation.json
```

The deposited endpoint is reconstructed as `ln(Km) + ln(kcat/Km)`. Only ranking within
each cofactor is interpreted because the source-unit constant is unresolved. Every
variant has a UniKP homolog at the 30% identity/80% coverage gate, so this is a mutation
sensitivity test rather than homology-cold validation. The frozen BioCandidateRanker and
DLKcat results both underperform random ranking expectations and cannot support candidate
selection. See `docs/CURRENT_RESULTS.md` for identities and metrics.

## Tests

```powershell
$env:PYTHONPATH="src"
python -m unittest discover -s tests -p "test_*.py" -v
```

## Boundaries

- DisorderFlow data is not used.
- FBA outputs are simulation context, not experimental truth.
- UniKP/DLKcat predictions are never promoted to labels.
- Existing biological project directories are not modified.
- Current outputs are research artifacts, not wet-lab or clinical recommendations.

See `docs/ACCEPTANCE_CRITERIA.md` and `docs/RESEARCH_PROTOCOL.md` before interpreting
metrics.

Frozen internal and external transfer results are reported in `docs/CURRENT_RESULTS.md`.
The internal RMSE gain is modest, and the external ranking result is not strong enough
for candidate selection or a publication-grade efficacy claim.
