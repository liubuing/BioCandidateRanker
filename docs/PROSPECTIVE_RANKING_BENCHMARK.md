# Prospective Candidate-Ranking Benchmark Gate

This is a governance interface for a future, prospectively run candidate-ranking
benchmark. It is separate from the already closed EnzEngDB and IMDH evaluations. The
implementation creates audit and freeze manifests only; it does not supply candidates,
experimental labels, model predictions, or a readiness claim.

## Required Evidence

The pre-label gate requires five inputs:

- a protocol with configurable minimum campaigns, distinct families, and candidates per
  campaign, plus frozen campaign-level ranking metrics;
- a candidate CSV containing exactly identified campaign and family membership and an
  `active`, `inactive`, or `censored` disposition for every attempted candidate;
- a campaign registry declaring exact roster counts, label-independent selection, roster
  completeness, and that labels were unavailable at freeze time;
- a family audit proving exact roster-family coverage and zero overlap with identified,
  SHA256-bound development and closed-test family sets;
- a model freeze listing unique model IDs and byte-verified checkpoint identities, with an
  explicit declaration that checkpoints were frozen before labels were available.

Label-like columns are prohibited in the candidate roster. Inactive and censored entries
must not be filtered out. The example at
`configs/prospective_candidate_ranking_protocol.example.json` defines protocol gates but
intentionally contains no data and cannot establish readiness.

## Lifecycle

First create a new pre-label readiness receipt. Output files are exclusive-create and are
never overwritten:

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.prospective_benchmark audit-freeze `
  --protocol protocol.json --candidates candidates.csv `
  --campaign-registry campaigns.json --family-audit family-audit.json `
  --model-freeze models.json --output readiness-freeze.json
```

A blocked audit still writes its issues and exits nonzero. Do not generate or deposit
predictions unless `ready` is `true`.

After a passing freeze, generate predictions without opening labels. Deposit a CSV with
exactly `model_id,candidate_id,prediction`; every frozen model-candidate pair is required:

```powershell
python -m biocandidate.prospective_benchmark deposit-freeze `
  --readiness readiness-freeze.json --protocol protocol.json `
  --candidates candidates.csv --campaign-registry campaigns.json `
  --family-audit family-audit.json --model-freeze models.json `
  --predictions blind-predictions.csv --output prediction-deposit.json
```

Only after the blind deposit is frozen may the independent label custodian expose a label
CSV with exactly `candidate_id,disposition,label`. It must contain every candidate;
censored rows may have an empty label, while non-censored rows require a finite numeric
label. Perform the one final join:

```powershell
python -m biocandidate.prospective_benchmark evaluate-once `
  --readiness readiness-freeze.json --deposit prediction-deposit.json `
  --candidates candidates.csv --predictions blind-predictions.csv `
  --labels custodian-labels.csv --output final-evaluation-receipt.json
```

The final receipt reports each metric for each campaign and the macro mean across
campaigns. Evaluation exclusively creates `<prediction-deposit>.consumed.json` before
metric computation, so the deposit cannot be evaluated again under another output path;
the marker remains if final receipt writing fails. Operational governance must make the
deposit directory append-only and deny the modeling team access to labels until deposit
acceptance. Final results are confirmatory and cannot be used to tune, select, or replace
a model.
