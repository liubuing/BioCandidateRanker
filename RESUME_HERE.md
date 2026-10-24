# Resume Here

The project was archived on 2026-07-27.

Start with:

1. `docs/PROJECT_STATE_ARCHIVE_2026-07-27.md`
2. `artifacts/project-state-archive-2026-07-27.json`
3. `configs/software_implementation_release_v11.json`
4. `configs/external_data_blockers.json`
5. `docs/DATA_ARRIVAL_RUNBOOK.md`

Then run:

```powershell
python scripts\verify_software_release.py --manifest configs\software_implementation_release_v11.json
python -m pytest -q
python -m ruff check .
```

Current scientific state: `blocked_pending_external_data`.

Current external benchmark pool: 209 accepted records / 28 global families / 66
substrates (201 / 28 / 65 after the eight documented quality holds). Do not
generate predictions for this pool until all frozen readiness gates pass; the
raw pool is 91 records and two families short, and the carry-forward pool is 99
records and two families short.

External requests are drafts only; zero messages were sent. Sender name `chy` was supplied,
but organization and email were not. No real prospective campaign collaborator or
independent label custodian was supplied.
