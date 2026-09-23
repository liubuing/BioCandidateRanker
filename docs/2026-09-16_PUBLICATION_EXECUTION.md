# Publication execution ledger (2026-09-16)

## Paper scope and claim boundary

The next manuscript is a `log10(kcat)` prediction and evaluation paper. The current
homology-cold result is an internal comparison on a frozen UniKP-associated corpus.
Retrained UniKP is close to BioCandidateRanker (RMSE 1.4093 vs. 1.3916); this supports
parity or a small possible improvement, not established superiority. EnzEngDB fitness
and IMDH relative mutation effects are different endpoints and show poor transfer.
No reliable candidate-ranking, absolute Km/activity, or experimental-flux claim is made.

## Workstream A: manuscript repair

- [x] Incorporate three-seed UniKP and DLKcat Mode B rows from
  `artifacts/external/sota-homology-cold/summary.json` in both main manuscript sources.
- [x] Correct the known fatally misattributed citations in both main manuscripts.
- [x] State that the ESM-2 protein-only control is missing and that the existing
  fusion-specific claim is unresolved.
- [ ] Check every in-text citation and DOI mechanically; inspect both supplementary
  manuscripts and regenerate any DOCX/PDF exports from corrected sources.
- [ ] Add actual author/contact/affiliation information and a stable public code archive
  at submission. Establish a data-availability statement consistent with the source
  corpus's unverified redistribution permission.

## Workstream B: representation and modality controls

The existing homology test and SOTA comparison have been observed. New runs on that
test are **post-hoc diagnostic**; they cannot be presented as a fresh confirmatory
selection experiment. Freeze the following comparison before training:

Execution update (2026-09-16): the original checkpoints were inspected and the model
configuration, initial learning rate (0.0003), source and split identities, and three
seeds were recovered. Both source and split files are accessible and their SHA256
hashes match the checkpoints. The post-hoc protocol is frozen at
`configs/esm2_modality_diagnostic_protocol.json`. The first protein-only seed 7 run
has begun and produced its first checkpoint; it is not yet a completed/evaluated result.
The sequential runner for the remaining runs is
`scripts/run_esm2_modality_diagnostic.py`. A hidden continuation process is waiting
for the first run to exit; it will verify the first run's exit status and checkpoints,
then run/evaluate all six protocol cells sequentially. Its progress and errors are
recorded in `artifacts/esm2-modality-diagnostic/continuation.log` and
`continuation.err`; each subsequent training/evaluation gets a per-run log.

Progress update (2026-09-17): the first protein-only seed 7 run exited normally after
epoch 5 and was evaluated once on the frozen internal test. It scored RMSE 1.5013,
MAE 1.1597, Pearson 0.2171; the corresponding previously evaluated full-model seed 7
scored RMSE 1.3880, MAE 1.0683, Pearson 0.3906. This is one post-hoc seed only and
is not a final fusion estimate. The sequential runner has started protein-only seed 42;
its first epoch completed without error. All remaining cells are pending.

Completion update (2026-09-17): all six new runs finished normally and each checkpoint
was evaluated once on the frozen internal test. `scripts/summarize_esm2_modality_diagnostic.py`
verified 1,646 test records and matching source/split identities for every variant and
wrote `artifacts/esm2-modality-diagnostic/summary.json`.

| ESM-2 variant | RMSE mean ± SD | MAE mean ± SD | Pearson mean ± SD |
|---|---:|---:|---:|
| Protein only | 1.5025 ± 0.0064 | 1.1656 ± 0.0051 | 0.2174 ± 0.0410 |
| Protein + substrate | 1.4454 ± 0.0156 | 1.1200 ± 0.0149 | 0.2791 ± 0.0294 |
| Full reference | 1.3916 ± 0.0183 | 1.0708 ± 0.0117 | 0.3810 ± 0.0222 |

The paired mean RMSE changes are −0.0571 ± 0.0213 when adding substrate and
−0.0538 ± 0.0147 when adding context; all six paired seed differences have the
same favorable direction. These results support a post-hoc internal modality
diagnostic, not prospective or independent generalization. The English and Chinese
main and supplementary manuscripts now include the results and claim boundary.

1. Frozen ESM-2 protein-only (`use_protein=true`, molecule/context false).
2. Frozen ESM-2 protein + molecule (context false).
3. Frozen ESM-2 full model (existing reference; do not retune it).

For 1 and 2, use the same source identity, MMseqs2 split, target, duplicate policy,
ESM-2 t6 backbone, 256-dimensional head, optimization schedule, seeds 7/42/123,
checkpoint criterion, and metric implementation as the existing full model. Select
checkpoints by validation only. Record train/validation/test row identities, run flags,
source SHA256, per-seed metrics, and predictions. Analyze paired per-seed RMSE, MAE,
and Pearson; do not claim a fusion benefit from an unpaired comparison. Run the same
predeclared variants on scaffold-cold only after freezing its complete protocol.

The CLI already exposes `--disable-molecule`, `--disable-context`,
`--protein-encoder esm2`, `--split-strategy`, and `--split-manifest`. The recovered
checkpoints contain architecture and source identities. Training schedule details not
serialized in the checkpoints are taken from the manuscript and explicitly declared
in the new protocol before any diagnostic test evaluation.

## Workstream C: independent same-endpoint evaluation

The temporal pool has 192 admitted records, 25 global families, and 54 substrates.
The frozen gate requires at least 300 records, 30 families, and 50 substrates, so
108 records and five families remain. Work from `docs/TARGETED_DATA_REQUESTS.md`:
request the Purdue Cdc14 deposition, then missing constructs or saturation evidence
from already audited zero-hit families, and seek an independent five-family assay panel.
Each delivery must pass license, exact-construct, substrate, assay, provenance, temporal,
exact-overlap, and MMseqs checks. The independent custodian retains labels until the
final record list and checkpoints are frozen. Do not score the current 192-record pool.

## Decision gate

- **Independent pool passes and scores well:** submit a restrained methods/evaluation
  article with same-endpoint external performance and a transparent baseline comparison.
- **Pool remains blocked:** prepare a narrower benchmarking/negative-results article.
  Explicitly state the absence of an independent same-endpoint test; select a venue whose
  requirements permit that evidence level. Do not relax the frozen pool criteria.

No external requests have been sent. Sending them requires a named sender,
institutional affiliation, and contact address.

Execution update (2026-09-17): [the kcat acquisition plan](2026-09-17_KCAT_DATA_ACQUISITION.md)
now ranks four source-specific evidence requests and an independent multi-family assay
invitation. The author-facing drafts are in `docs/external-requests/`. Purdue remains
inaccessible from this environment; the other three source audits identify exact
construct, substrate-structure, or O2-saturation gaps. Known rows from the phage and
multicopper sources can contribute at most 24 before the remaining acceptance gates,
so an independent panel remains necessary. All five drafts are unsent; no new records
were admitted and no external pool was scored.

Planning update (2026-09-17): the Chinese [wet-lab validation plan](2026-09-17_WET_LAB_VALIDATION_PLAN_CN.md)
and its Word report specify an eight-family prospective kcat cohort, staged pilot,
staffing, costs, blinded comparison, and an optional engineering-ranking campaign.
These are proposed experiments; none has been run. Unpublished new measurements must
be evaluated under a separately preregistered prospective protocol and cannot be
automatically counted in the existing publication-gated temporal pool. If new families
alone supply the 108-record deficit, the 20-per-family cap requires at least six families.

Target update (2026-09-17): the user selected **PLOS Computational Biology** as the
primary publication target. The [new research plan](PLOS_COMPUTATIONAL_BIOLOGY_PLAN_CN.md)
supersedes broad venue targeting and prioritizes biological insight, modern baseline
comparisons, and independent verification of transfer to candidate ranking. This is
a planning update, not a new scientific result or a submission-readiness declaration.
Frozen evaluation gates remain unchanged. Next: CatPred/CataPro feasibility and
endpoint/provenance audits before defining additional experiments.

Execution update (2026-09-18): CataPro fixed-split retraining is complete for seeds
7/42/123. Frozen single-pass internal test: RMSE 1.3876 ± 0.0195, MAE 1.0650 ± 0.0105,
Pearson 0.4320 ± 0.0004. Checkpoint/prediction identities and saved-metric recomputation
passed; seven focused tests passed. These results do not establish BioCandidateRanker
superiority. See `docs/CATAPRO_INTERNAL_COMPARISON_CN.md` and
`docs/CATAPRO_MODE_B_RUNBOOK_CN.md`. No independent temporal pool was scored.
