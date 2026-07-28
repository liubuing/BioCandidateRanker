# Project State Archive: 2026-07-27

## Status

The software implementation is closed as `software-implementation-v1`. Scientific work is
paused at `blocked_pending_external_data`.

The independent absolute-kinetics curation pool currently contains:

- 192 of 300 required records;
- 25 of 30 required global MMseqs families;
- 54 of 50 required substrates;
- 18 accepted sources;
- 203 mapped component sequences.

The pool is not a benchmark yet. It has not been scored by the model, and predictions remain
prohibited until the complete 300/30/50 gate passes and final records and checkpoints are
frozen.

## Completed

- Strict normalized `kcat`, `Km`, and `kcat/Km` ingestion.
- Frozen source, split, model, and artifact identities.
- Pinned MMseqs homology and global family audit.
- Validation-only grouped uncertainty calibration and conformal intervals.
- Explicit FBA/flux contracts and simulated-flux smoke training.
- Prospective candidate-ranking freeze, blind-deposit, and one-time evaluation software.
- Targeted public acquisition audits for absolute kinetics, Km, experimental flux, and a
  genuinely unobserved campaign.
- External request drafts, response tracker, acceptance rules, and campaign handoff package.

## Not Completed

- The additional 108 eligible absolute-kinetics records and five global families.
- A frozen and evaluated independent external benchmark.
- A governed non-`kcat` checkpoint trained from row-level experimental provenance.
- A real blind candidate campaign with an independent label custodian.

These are external evidence gaps, not missing software features.

## External Requests

No messages have been sent. The request queue is
`artifacts/external/external-request-queue.json` and records `messages_sent: 0`.

Verified recipients:

- Purdue Cdc14: Mark C. Hall, `mchall@purdue.edu`.
- SABIO-RK: `sabiork@h-its.org`.
- yeast-GEM: Eduard Kerkhoven, `eduardk@chalmers.se`, with the official GitHub issue form
  as fallback.

The supplied sender name is `chy`. Sender organization and sender email are still missing.
No prospective campaign collaborator or independent custodian has been identified.

## Key Evidence

- Software release: `configs/software_implementation_release_v1.json`
- External blockers: `configs/external_data_blockers.json`
- Combined evidence receipt:
  `artifacts/external/evidence-phase-receipt-2026-07-27.json`
- Global readiness:
  `artifacts/external/temporal-global-family-audit/readiness-audit.json`
- Km blocker: `artifacts/external/governed-km-sabio-rk/blocker.json`
- Experimental flux blocker: `artifacts/experimental-yeast-flux/blocker.json`
- Campaign dependency receipt:
  `artifacts/external/prospective-candidate-ranking/external-dependency-receipt.json`
- Campaign collaborator handoff:
  `artifacts/external/prospective-candidate-ranking/handoff`
- Data-arrival instructions: `docs/DATA_ARRIVAL_RUNBOOK.md`

## Last Verification

```text
139 passed, 8 subtests passed
Ruff: all checks passed
software release verification: valid, no identity issues
```

## Resume

Run these commands before changing anything:

```powershell
python scripts\verify_software_release.py
python -m pytest -q
python -m ruff check .
wsl -d Ubuntu-24.04 mmseqs version
```

Expected pinned MMseqs version:

```text
5d152c612b6ad2a56f657b7a02c127eceaea2a75
```

Re-check external licenses, contacts, APIs, and download endpoints because they may change
during a long pause. Do not overwrite the v1 release manifest. New governed data or code
should produce a new release manifest.

When new experimental data arrive, follow `docs/DATA_ARRIVAL_RUNBOOK.md`. Attempt one
bounded correction for an isolated source problem; if it remains unresolved, record and
skip it. Do not lower scientific gates to force completion.
