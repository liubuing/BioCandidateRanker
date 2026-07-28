# Data Arrival Runbook

This runbook starts only after a targeted external delivery arrives. Work on one source at
a time; an isolated failure is corrected once, then recorded and skipped so other sources
continue.

## 1. Preserve The Delivery

1. Copy the unmodified delivery into a source-specific external evidence directory.
2. Record source URL/DOI, version, license, retrieval timestamp, byte size, and SHA256.
3. Keep originals read-only. Write normalized files separately.
4. Reject predictions, inferred values, screenshots without source rows, and unknown units.

## 2. Absolute Kinetics

1. Convert direct `kcat`, `Km`, and `kcat/Km` rows to the v1 normalized JSONL contract.
2. Audit exact constructs, PubChem structures, source cells, units, saturation, and formulae.
3. Search exact constructs against the frozen UniKP FASTA with pinned MMseqs2.
4. Regenerate the strict registry and global family audit.
5. Stop before predictions unless records are at least 300, families at least 30, and
   substrates at least 50.

```powershell
$env:PYTHONPATH="src"
python -m biocandidate.temporal_registry `
  --output artifacts\external\temporal-absolute-kinetics\registry-readiness-audit.json

python -m biocandidate.global_temporal_audit `
  --mmseqs "mmseqs" `
  --output-dir artifacts\external\temporal-global-family-audit
```

When every gate passes, freeze the final manifest and checkpoint identities before running
any model against the pool.

## 3. Km Training

1. Require original linear Km values and units plus row-level provenance.
2. Convert to `log10(M)` with `normalized-kinetics`; do not infer the unit from value range.
3. Audit duplicates and conflicts, then build an MMseqs family-aware split.
4. Train only `log10_km`; bind the source and split identities into the checkpoint.

```powershell
python -m biocandidate.cli train `
  --data artifacts\governed-km.jsonl `
  --data-format normalized-kinetics `
  --tasks log10_km `
  --manifest artifacts\governed-km-manifest.json `
  --split-manifest artifacts\governed-km-split.json `
  --output-dir artifacts\governed-km-training
```

Fit calibration from validation outputs only. Keep any independent test labels closed.

## 4. Absolute Activity Training

Do not reuse the absolute-kinetics adapter until a dedicated activity adapter has been
implemented for the delivered canonical unit. First freeze the endpoint definition,
denominator, assay compatibility rules, provenance schema, and family-aware split. Add
adapter tests before training. EnzEngDB remains campaign ranking only.

## 5. Experimental Flux Training

1. Validate `FluxLabelMetadata` as `experimental` and bind model, solver, objective,
   condition, method, and source identities.
2. Map every target to the frozen reaction namespace and record mapping confidence.
3. Exclude the target reaction and solved fluxes from FBA input features.
4. Freeze condition/strain splits, then train a separate experimental-flux checkpoint.
5. Never merge the simulated smoke labels into the experimental target.

## 6. Prospective Candidate Ranking

Follow `docs/PROSPECTIVE_RANKING_BENCHMARK.md` exactly:

1. `audit-freeze` before labels;
2. generate and `deposit-freeze` complete blind predictions;
3. receive labels from the independent custodian;
4. run `evaluate-once` and preserve the consumption marker and final receipt.

## 7. Calibration And Verification

```powershell
python -m biocandidate.cli fit-calibration `
  --checkpoint <best.pt> --data <data> --manifest <manifest.json> `
  --split-manifest <split.json> --task <task> --method grouped-cv `
  --cv-folds 5 --output <calibration.json> --device cpu

python -m pytest -q
python -m ruff check .
python scripts\verify_software_release.py
```

The final release verifier checks the closed software baseline. New governed datasets and
checkpoints require a new release manifest rather than silently changing v1 identities.
