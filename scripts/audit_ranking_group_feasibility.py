"""Metadata-only inventory; no labels, predictions, or ranking selection."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
AUDIT = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted'
OUT = ROOT / 'artifacts/research-scope-20261005-v2'


def missing_condition(value):
    return value.strip().casefold() in {'', 'not reported', 'not specified', 'unknown', 'na', 'n/a', 'none', 'null'}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    audit = json.loads((AUDIT / 'readiness-audit.json').read_text())
    mapping = AUDIT / 'global-accepted-mapping.csv'
    assert sha(mapping) == audit['artifacts']['mapping_manifest_sha256']
    records = defaultdict(list)
    with mapping.open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            if row['retained_after_family_cap'] == 'True':
                records[(row['source_id'], row['candidate_id'])].append(row)
    assert len(records) == audit['counts']['records_after_global_family_cap']
    metadata = {}
    hashes = {}
    for source in sorted({s for s, _ in records}):
        path = BASE / source / 'candidate_records.csv'
        hashes[source] = sha(path)
        with path.open(encoding='utf-8-sig', newline='') as f:
            for row in csv.DictReader(f):
                key = (source, row.get('candidate_id') or row.get('source_record_id'))
                if key in records:
                    # Deliberately whitelist only metadata: never use numeric response columns.
                    metadata[key] = {field: row.get(field, '').strip() for field in
                                     ['assay_pH', 'assay_temperature_C', 'condition_id', 'reaction_direction']}
    assert records.keys() <= metadata.keys()
    groups = defaultdict(list)
    for key, components in records.items():
        m = metadata[key]
        signature = tuple(x['sequence_sha256'] for x in sorted(components, key=lambda x: int(x['component_index'])))
        group = (key[0], components[0]['substrate_id'], m['assay_pH'], m['assay_temperature_C'])
        groups[group].append(dict(candidate_id=key[1], sequence_signature=signature,
                                  multi_component=len(components) > 1,
                                  reported_condition_id=m['condition_id'], reported_direction=m['reaction_direction']))
    inventory = []
    for key, rows in sorted(groups.items()):
        inventory.append(dict(source_id=key[0], substrate_id=key[1], reported_pH=key[2],
                              reported_temperature_C=key[3], records=len(rows),
                              distinct_sequence_signatures=len({r['sequence_signature'] for r in rows}),
                              missing_pH_or_temperature=missing_condition(key[2]) or missing_condition(key[3]),
                              any_multi_component=any(r['multi_component'] for r in rows),
                              candidates=rows, status='metadata_bucket_not_validated_campaign'))
    prospective = [g for g in inventory if g['distinct_sequence_signatures'] >= 3]
    result = dict(date='2026-10-05', input_records=len(records), metadata_buckets=len(inventory),
                  buckets_with_at_least_three_sequence_signatures=len(prospective),
                  sources_in_these_buckets=len({g['source_id'] for g in prospective}),
                  these_buckets_with_missing_pH_or_temperature=sum(g['missing_pH_or_temperature'] for g in prospective),
                  confirmed_independent_campaigns=None,
                  claim_boundary='Feasibility inventory only. Equal source, substrate, pH and temperature do not establish equal buffer, cofactors, direction, stoichiometry or experimental campaign. Distinct sequence signatures are not necessarily distinct enzyme candidates. Three is an inventory display threshold, not a frozen evaluation gate.',
                  labels_used=False, predictions_run=False, mapping_sha256=sha(mapping),
                  source_csv_sha256=hashes, groups=inventory)
    (OUT / 'ranking-metadata-feasibility.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['groups', 'source_csv_sha256', 'claim_boundary']}))


if __name__ == '__main__':
    main()
