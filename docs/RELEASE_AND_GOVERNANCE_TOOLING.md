# Release and governance tooling

These scripts implement the submission-facing checks that the manuscripts assume.
Each one exits non-zero when a claim would not survive review.

## Why these exist

The manuscripts state that a public archived release is "required before submission"
and that "redistribution permission for the development corpus remains unresolved".
Both statements were aspirations without tooling. This set makes them executable:
the archive can now be built and verified repeatably, and the corpus terms are
recorded in one auditable place.

## Tools

### Build a deposit archive

```bash
python scripts/build_release_archive.py
```

Writes `artifacts/release-archive/<release-id>.zip` plus a sidecar manifest.

The archive holds every Git-tracked file and every frozen evidence artifact named by
`configs/software_implementation_release_v6.json`. Entries are sorted with fixed
timestamps, so rebuilding the same commit yields a byte-identical archive — the
SHA256 is stable enough to cite in a manuscript.

The builder refuses to produce an archive when:

- a frozen artifact belongs to a corpus that may not ship (`--skip-license-check` overrides
  this for local inspection only);
- `--fail-on-local-paths` is set and a member embeds an absolute local path.

### Verify an archive

```bash
python scripts/verify_release_archive.py artifacts/release-archive/software-implementation-v6.zip
```

Recomputes every member hash against the embedded `RELEASE_MANIFEST.json`, checks that
nothing is unlisted or absent, confirms frozen evidence still matches the release
manifest, and confirms the archive was built from the current manifest snapshot.

### Corpus license audit

```bash
python scripts/audit_development_corpus_license.py --strict \
  --json-out artifacts/license-audit/development-corpus-license-audit.json \
  --markdown artifacts/license-audit/development-corpus-license-audit.md
```

Reads `configs/data_license_inventory.json` and checks that every claim is complete,
evidenced by a path that exists, and consistent with the manifests that record license
state elsewhere. `--strict` fails on any inconsistency.

The inventory also maps frozen artifact paths to corpora, which is what lets the
archive builder decide whether a file may ship.

### Data availability statement

```bash
python scripts/generate_data_availability_statement.py \
  --out artifacts/license-audit/data-availability-statement.md
```

Drafts the manuscript statement from the audit. It separates what ships, what a reader
must obtain from the original deposit, and what stays unresolved. The draft still needs
author review; it never asserts a right the audit did not record.

### Prospective pool dashboard

```bash
python scripts/prospective_pool_dashboard.py --markdown artifacts/prospective-pool/dashboard.md
```

Compares live audit counts against the frozen thresholds and reports per-gate shortfalls,
per-family saturation, and what each resolution path would cost. This supports the
decision recorded in `docs/PLOS_COMPUTATIONAL_BIOLOGY_PLAN_CN.md`, and never scores
records or generates predictions.

### External request status

```bash
python scripts/external_request_status.py --markdown artifacts/external-requests/status.md
```

Reconciles the three records of outreach: `artifacts/external/external-request-queue.json`,
`artifacts/external/external-response-tracker.csv`, and the drafts under
`docs/external-requests/`. Reports drafts with no queue entry, queue entries whose draft
is missing, and genuine send-state contradictions. Equivalent labels across the two
vocabularies are not reported as conflicts.

### Figure 1

```bash
python scripts/figure1_data_task_map.py
```

Regenerates `manuscript/figures/figure1-data-task-map.svg` and its legend. One lane per
corpus, so labels cannot collide: position encodes sequence novelty, colour encodes the
endpoint, area encodes record count, and a dashed outline marks a corpus that has never
been scored.

### Interactive readiness dashboard

```bash
# Refresh the four JSON inputs, then build the page.
python scripts/audit_development_corpus_license.py \
  --json-out artifacts/license-audit/development-corpus-license-audit.json
python scripts/prospective_pool_dashboard.py \
  --json-out artifacts/prospective-pool/dashboard.json
python scripts/external_request_status.py \
  --json-out artifacts/external-requests/status.json
python scripts/build_release_archive.py
python scripts/build_readiness_dashboard.py
```

Writes `artifacts/dashboard/readiness-dashboard.html`: one self-contained page with five
tabs (Overview, Prospective pool, Corpus licenses, External requests, Archive). No network
access and no build step, so it opens from disk or from a local static server.

What is interactive:

- tab navigation, with the current tab reflected in the URL hash so back and forward work;
- the corpus table sorts by any column, filters by ship status, and searches across
  corpus, role, licence, and redistribution state;
- clicking a corpus row opens its blocking reason and evidence count;
- each overview blocker expands to list the affected files and the exact command that
  clears it.

The page embeds absolute local paths when they exist, so it is generated under the
gitignored `artifacts/` tree. Pass `--redact-paths` to build a version safe to share:

```bash
python scripts/build_readiness_dashboard.py --redact-paths \
  --output artifacts/dashboard/readiness-dashboard-shareable.html
```

The headline verdict is derived, not authored: the page reports "Not submission ready"
with the count of blocking items, and each blocker is computed from a specific record
rather than a hand-maintained list.

## Known open items surfaced by these tools

- **No external request has been sent.** Four drafts have no queue entry at all, and two
  still carry a placeholder recipient.
- **18 archived files embed absolute local paths**, including `/home/liubuing/...` in
  `artifacts/external/temporal-global-family-audit/readiness-audit.json`. Because the Git
  author and the upstream repository account are both `liubuing`, this is a direct
  de-anonymization vector for a double-blind submission. Redact before deposit.
- **The training corpus cannot be redistributed**, so the deposit is code, protocols, and
  derived artifacts only. `artifacts/full-homology-baseline/best.pt` is a checkpoint
  trained on that corpus; it ships as a declared derivative.
- **The prospective pool is 108 records and 5 families short** of its frozen gate, and
  two families are already at the 20-record cap.

## Order of operations before submission

1. `audit_development_corpus_license.py --strict` — resolve any inconsistency.
2. `prospective_pool_dashboard.py` — record the decision on the pool gate.
3. `external_request_status.py` — record what outreach has actually happened.
4. `figure1_data_task_map.py` — regenerate Figure 1.
5. `build_release_archive.py --fail-on-local-paths` — redact first, then build.
6. `verify_release_archive.py` — verify the exact file you will deposit.
7. `generate_data_availability_statement.py` — draft the statement from the final audit.
8. `build_readiness_dashboard.py` — render the state for a human reader.

Steps 1 to 3 write the JSON the dashboard reads, so run them first or the page will
report missing sources.

## Current dashboard verdict

At the last build the page reports **Not submission ready — 6 blocking items**:
unsatisfied local-path redaction, the unmet prospective-pool gate, zero external
requests sent, and three unwritten submission components (Author Summary, cover letter,
placeholder fields).
