# Resume Here

The project was archived on 2026-07-27.

Start with:

1. `docs/PROJECT_STATE_ARCHIVE_2026-07-27.md`
2. `artifacts/project-state-archive-2026-07-27.json`
3. `configs/software_implementation_release_v1.json`
4. `configs/external_data_blockers.json`
5. `docs/DATA_ARRIVAL_RUNBOOK.md`

Then run:

```powershell
python scripts\verify_software_release.py
python -m pytest -q
python -m ruff check .
```

Current scientific state: `blocked_pending_external_data`.

Current external benchmark pool: 192 records, 25 global families, 54 substrates. Do not
generate predictions for this pool until all frozen readiness gates pass.

External requests are drafts only; zero messages were sent. Sender name `chy` was supplied,
but organization and email were not. No real prospective campaign collaborator or
independent label custodian was supplied.
