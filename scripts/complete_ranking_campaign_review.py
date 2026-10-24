"""Classify all remaining metadata buckets without scoring or changing the pool.

Closes the 18 remaining campaign buckets with a fixed per-source disposition
table grounded in published experimental-design metadata only (never response
magnitudes), re-reading each source's candidate rows against recorded SHA256s
and asserting the exact disposition counts. Writes
campaign-audit-complete-20261005/closeout.json; the pool stays unchanged,
evaluation_ready stays false everywhere, and the frozen protocol and prior
pool hashes are re-verified.

中文：本脚本在不打分、不改动池的情况下为全部剩余元数据桶分类：用一张仅基于已发表实验
设计元数据（绝不使用响应幅值）的固定逐源处置表收口剩余 18 个候选桶，按记录的 SHA256
复读各来源候选行并断言精确的处置计数。写出
campaign-audit-complete-20261005/closeout.json；池保持不变、处处 evaluation_ready=false，
冻结协议与先前池哈希被复验。治理角色：候选组评审收口审计收据。
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
OUT = BASE / 'campaign-audit-complete-20261005'


def sha(p):
    """Return the SHA256 hex digest of a file.

    中文：返回文件的 SHA256 十六进制摘要。
    """
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    """Apply the disposition table to every queued bucket and write the closeout.

    中文：把处置表应用到每个排队桶并写出收尾收据。
    """
    initial = json.loads((BASE / 'campaign-audit-20261005/closeout.json').read_text())
    inventory = json.loads((ROOT / 'artifacts/research-scope-20261005-v2/ranking-metadata-feasibility.json').read_text())
    groups = {(g['source_id'], g['substrate_id']): g for g in inventory['groups'] if g['distinct_sequence_signatures'] >= 3}
    # Decisions concern metadata and published experimental design, never response magnitude.
    decisions = {
        'biostudies-S-EPMC12412097': ('hold_metadata', 'Temperature and reaction-buffer identity unresolved; physical-state differences require separate interpretation.', 'provenance.json'),
        'biostudies-S-EPMC12645435': ('common_conditions_candidate', 'Inosine/adenine curves share chemical conditions; enzyme doses vary. Check normalization/initial-rate comparability; other substrates are not pooled.', 'raw/cb5c00120_si_001.pdf'),
        'biostudies-S-EPMC12667214': ('common_conditions_candidate', 'Separate CMP-series, ATP-series and pyrophosphate-series. CMPK reciprocal series share constructs and are not independent campaigns. Coupled reporter enzymes are not extra candidates.', 'raw/PMC12667214-fullText.xml'),
        'europepmc-PMC10520331': ('common_conditions_candidate', 'CM and CDT are distinct reaction/readout tasks. Common buffer and temperature documented; enzyme dose and construct architecture vary; two biological replicates share one study.', 'raw/PMC10520331-fullText.xml'),
        'europepmc-PMC11656708': ('common_conditions_candidate', 'Four compounds remain separate. Common CHES buffer/temperature documented in supplement; dose, substrate ranges and Zn loading vary. Preserve source/family cap and fit limitations.', 'raw/bi4c00350_si_001.pdf'),
        'europepmc-PMC11659886': ('common_conditions_candidate', 'Separate sulfated and nonsulfated substrate tasks. Shared fluorescence conditions documented for requested enzymes; concentrations vary and three parent enzymes/variants are correlated. Bt4394 different-pH comparators are not included.', 'raw/PMC11659886-fullText.xml'),
        'europepmc-PMC12284513': ('split_conditions', 'Split ethanol 2.5% and 10%; repeated WT records under different solvent conditions are not new candidates. Full tag mapping remains a caveat.', 'raw/PMC12284513-fullText.xml'),
        'europepmc-PMC12444516': ('hold_metadata', 'Two substrate buckets share nine parent proteins; kinetic temperature and full tag/linker provenance unresolved. Do not infer temperature from another paper.', 'raw/ao5c05616_si_001.pdf'),
        'europepmc-PMC12994016': ('hold_metadata', 'Kinetic temperature not reported; chromatography temperature is not a substitute. Preserve deposited/reference discrepancy and fitted-range caveats.', 'raw/PMC12994016.xml'),
        'zenodo-14055918': ('split_conditions', 'Split ATP preincubation versus none; reaction pH/temperature remain unresolved for the HPLC endpoint. Charcoal assay buffer is not automatically transferable.', 'raw/PMC11658767-fullText.xml')
    }
    reviewed = []
    for queued in initial['remaining_queue']:
        source, substrate = queued['source_id'], queued['substrate_id']
        g = groups[(source, substrate)]
        status, reason, evidence = decisions[source]
        ids = {x['candidate_id'] for x in g['candidates']}
        with (BASE / source / 'candidate_records.csv').open(encoding='utf-8-sig', newline='') as f:
            rows = [r for r in csv.DictReader(f) if r['candidate_id'] in ids]
        assert len(rows) == g['records'] == queued['records']
        split = []
        field = {'europepmc-PMC12284513': 'ethanol_percent_vv', 'zenodo-14055918': 'atp_preincubation'}.get(source)
        if field:
            for val in sorted({r[field] for r in rows}):
                members = [r for r in rows if r[field] == val]
                split.append(dict(condition_field=field, value=val, candidate_ids=[r['candidate_id'] for r in members],
                                  distinct_sequences=len({r['sequence_id'] for r in members})))
        path = BASE / source / evidence
        assert path.exists()
        reviewed.append(dict(source_id=source, substrate_id=substrate, candidate_ids=sorted(ids),
                             records=len(rows), disposition=status, reason=reason,
                             missing_pH_or_temperature=g['missing_pH_or_temperature'],
                             condition_splits=split, evidence_file=str(path.relative_to(ROOT)), evidence_sha256=sha(path),
                             source_csv_sha256=sha(BASE / source / 'candidate_records.csv'),
                             study_dependency_block=source, evaluation_ready=False))
    assert len(reviewed) == 18 and len({(r['source_id'], r['substrate_id']) for r in reviewed}) == 18
    counts = dict(Counter(r['disposition'] for r in reviewed))
    assert counts == {'hold_metadata': 4, 'common_conditions_candidate': 12, 'split_conditions': 2}
    pool = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json'
    prior = json.loads((BASE / 'targeted-curation-20261004/closeout.json').read_text())
    assert sha(pool) == prior['global_audit_sha256']
    protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    assert sha(protocol) == prior['protocol_sha256']
    result = dict(date='2026-10-05', remaining_queue_disposition_complete=True, reviewed_this_round=18,
                  reviewed_total=21, unresolved_experimental_details_remain=True, counts=counts,
                  common_condition_candidate_sources=len({r['source_id'] for r in reviewed if r['disposition']=='common_conditions_candidate'}),
                  corrected_inventory_missing_conditions=inventory['these_buckets_with_missing_pH_or_temperature'],
                  new_accepted_records=0, pool_counts=json.loads(pool.read_text())['counts'], evaluation_ready_groups=0,
                  model_predictions_run=False, selection_uses_response_values=False,
                  prior_pool_sha256=sha(pool), protocol_sha256=sha(protocol), groups=reviewed,
                  prior_three_group_receipt_sha256=sha(BASE / 'campaign-audit-20261005/closeout.json'))
    (OUT / 'closeout.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['groups', 'pool_counts']}))


if __name__ == '__main__':
    main()
