"""Reconcile all 209 historical pool records; no training, scoring or label edits.

Record-level carry-forward is distinct from final benchmark readiness. Existing
source audits supply primary-evidence decisions; this run verifies their inputs
and applies the subsequently documented campaign quality holds.
"""
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
GLOBAL = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted'
OUT = ROOT / 'artifacts/external/pool-quality-register-20261005'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def coverage(rows):
    return dict(records=len(rows), sources=len({r['source_id'] for r in rows}),
                families=len({r['global_family_id'] for r in rows}),
                substrates=len({r['substrate_id'] for r in rows}),
                unique_sequence_signatures=len({tuple(r['sequence_signature']) for r in rows}))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protocol_path = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    protocol = read(protocol_path)
    audit = read(GLOBAL / 'readiness-audit.json')
    inv = read(ROOT / 'artifacts/research-scope-20261005-v2/ranking-metadata-feasibility.json')
    first = read(BASE / 'campaign-audit-20261005/closeout.json')
    rest = read(BASE / 'campaign-audit-complete-20261005/closeout.json')
    detail = read(BASE / 'campaign-detail-20261005/closeout.json')
    assert sha(protocol_path) == detail['protocol_sha256']
    assert sha(GLOBAL / 'readiness-audit.json') == detail['pool_sha256']
    assert sha(GLOBAL / 'global-accepted-mapping.csv') == inv['mapping_sha256'] == audit['artifacts']['mapping_manifest_sha256']
    assert sha(GLOBAL / 'global-accepted-components.fasta') == audit['artifacts']['global_fasta_sha256']
    assert sha(Path(audit['mmseqs']['cluster_tsv'])) == audit['mmseqs']['cluster_tsv_sha256']

    mapping = defaultdict(list)
    with (GLOBAL / 'global-accepted-mapping.csv').open(encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f):
            if r['retained_after_family_cap'] == 'True':
                mapping[(r['source_id'], r['candidate_id'])].append(r)
    # Recheck mapped component content against global FASTA, not only receipt text.
    fasta = {}
    for block in (GLOBAL / 'global-accepted-components.fasta').read_text().split('>')[1:]:
        lines = block.splitlines()
        fasta[lines[0].split()[0]] = ''.join(lines[1:])
    for components in mapping.values():
        for c in components:
            assert hashlib.sha256(fasta[c['global_sequence_id']].encode()).hexdigest() == c['sequence_sha256']

    reviews = {}
    for source, g in first['campaign_reviews'].items():
        for cid in g['ids']:
            reviews[(source, cid)] = dict(status=g['status'], notes=g['limitations'])
    for g in rest['groups']:
        for cid in g['candidate_ids']:
            reviews[(g['source_id'], cid)] = dict(status=g['disposition'], notes=[g['reason']])
    details = {(g['source_id'], m['candidate_id']): (g, m)
               for g in detail['groups'] for m in g['members']}
    buckets = {}
    for i, g in enumerate(inv['groups'], 1):
        for c in g['candidates']:
            buckets[(g['source_id'], c['candidate_id'])] = (f'bucket-{i:03d}', g)

    source_register = {}
    rows = []
    evidence = {}
    for source in sorted({k[0] for k in mapping}):
        folder = BASE / source
        p = folder / 'candidate_records.csv'
        assert sha(p) == inv['source_csv_sha256'][source]
        hp = folder / 'homology-audit.json'
        h = read(hp)
        assert h['mmseqs_version'] == audit['mmseqs']['required_version']
        assert (h['min_identity'], h['coverage'], h['coverage_mode']) == (.3, .8, 0)
        pp = folder / 'provenance.json'
        prov = read(pp) if pp.exists() else {}
        licenses = [v for obj in (prov, h) for k, v in obj.items()
                    if 'license' in k and isinstance(v, str) and v in protocol['license_gate']['accepted']]
        assert licenses, source
        caveats = {k: v for k, v in prov.items() if any(x in k for x in ('caveat', 'anomal', 'quality_', 'inference_boundary'))}
        source_register[source] = dict(article_doi=prov.get('article_doi', h.get('article_doi')),
            accepted_license_evidence=sorted(set(licenses)), source_caveats=caveats,
            primary_review_basis='existing source provenance and homology audit; not a new full-text re-curation',
            historical_homology_audit_sha256=sha(hp))
        for ep in [p, hp] + ([pp] if pp.exists() else []):
            evidence[str(ep.relative_to(ROOT))] = sha(ep)
        with p.open(encoding='utf-8-sig', newline='') as f:
            all_rows = list(csv.DictReader(f))
        selected = {r.get('candidate_id') or r.get('source_record_id'): r for r in all_rows}
        assert len(selected) == len(all_rows), source
        for key in sorted(k for k in mapping if k[0] == source):
            cid = key[1]
            r = selected[cid]
            value = float(r['kcat_s-1'])
            basic_valid = math.isfinite(value) and value > 0
            assert basic_valid, (source, cid)
            comps = sorted(mapping[key], key=lambda x: int(x['component_index']))
            families = {c['global_family_id'] for c in comps}
            assert len(families) == 1
            bid, bucket = buckets[key]
            signature = [c['sequence_sha256'] for c in comps]
            condition = {}
            if source == 'europepmc-PMC12284513':
                condition['ethanol_percent_vv'] = r['ethanol_percent_vv']
            elif source == 'zenodo-14055918':
                condition['atp_preincubation'] = r['atp_preincubation']
            elif source == 'biostudies-S-EPMC10795190':
                condition = {k: r[k] for k in ['metal', 'metal_concentration_mM']}
            elif source == 'biostudies-S-EPMC11955521':
                condition = {k: r.get(k, '') for k in ['fixed_reaction_substrate', 'fixed_reaction_substrate_mM']}
            elif source == 'biostudies-S-EPMC13033379':
                condition = {k: r[k] for k in ['constant_cosubstrate', 'constant_cosubstrate_concentration', 'potassium_phosphate_mM']}
            review = reviews.get(key, dict(status='not_yet_campaign_reviewed', notes=[]))
            held = key in details and details[key][0]['findings']['evidence_gate'] == 'hold_source_consistency_review'
            reasons = ['historical source admission retained; current input integrity and finite positive endpoint checked']
            if held:
                reasons = ['CtNDT inosine campaign: repeat-count conflict and D62N figure/table inconsistency; group review hold, not eight proven wrong values']
            flags = ['development_publication_overlap_unavailable', 'final_readiness_gate_closed']
            if bucket['missing_pH_or_temperature']:
                flags.append('missing_condition_allowed_if_unreported_by_original_absolute_protocol')
            if caveats:
                flags.append('source_caveats_apply_see_source_register')
            if source == 'biostudies-S-EPMC12645435':
                flags.append('SD_column_name_not_verified_uncertainty_type')
            if source == 'europepmc-PMC12362431':
                flags.append('NADPH_readout_not_product_specific_fucose_endpoint')
            row = dict(source_id=source, candidate_id=cid, article_doi=source_register[source]['article_doi'],
                substrate_id=comps[0]['substrate_id'], global_family_id=next(iter(families)),
                sequence_signature=signature, component_count=len(comps),
                parent_block=details[key][1]['parent'] if key in details else None,
                parent_block_unknown=key not in details,
                original_kcat_s_1=r['kcat_s-1'], original_value_modified=False,
                absolute_status='pending' if held else 'carry_forward_pass', absolute_reasons=reasons,
                absolute_pass_scope='record-level original protocol only; not final evaluation permission',
                ranking_status='pending', campaign_review_status=review['status'], campaign_notes=review['notes'],
                metadata_bucket_id=bid, condition_partition=condition,
                ranking_group_key=bid+'|'+json.dumps(condition, sort_keys=True),
                reported_pH=bucket['reported_pH'], reported_temperature_C=bucket['reported_temperature_C'],
                flags=flags, prediction_eligible=False,
                evidence=[str(p.relative_to(ROOT)), str(hp.relative_to(ROOT))]+([str(pp.relative_to(ROOT))] if pp.exists() else []))
            if key in details:
                row['detail_findings'] = details[key][0]['findings']
                row['uncertainty_interpretation'] = details[key][1]['uncertainty_interpretation']
            rows.append(row)

    pair_groups = defaultdict(list)
    for r in rows:
        pair_groups[(tuple(r['sequence_signature']), r['substrate_id'])].append(r)
    repeated_pairs = []
    for members in pair_groups.values():
        if len(members) > 1:
            pid = f'repeated-pair-{len(repeated_pairs)+1:03d}'
            for r in members:
                r['repeated_sequence_substrate_block'] = pid
                r['flags'].append('repeated_sequence_substrate_not_independent_enzyme_candidates')
            repeated_pairs.append(dict(block_id=pid,
                cross_source=len({r['source_id'] for r in members}) > 1,
                records=[dict(source_id=r['source_id'], candidate_id=r['candidate_id'],
                              pH=r['reported_pH'], temperature=r['reported_temperature_C'],
                              condition_partition=r['condition_partition']) for r in members],
                disposition='retain condition-specific measurements; no automatic averaging or deletion'))
    grouped = defaultdict(list)
    for r in rows:
        grouped[r['ranking_group_key']].append(r)
    group_summary = []
    for group_id, members in grouped.items():
        nseq = len({tuple(r['sequence_signature']) for r in members})
        for r in members:
            r['group_distinct_sequence_signatures'] = nseq
            if nseq < 2:
                r['ranking_status'] = 'not_applicable_single_sequence_in_current_pool'
            elif r['absolute_status'] == 'pending':
                r['ranking_status'] = 'pending_source_consistency'
            elif r['campaign_review_status'] == 'hold_metadata' or 'missing_condition_allowed_if_unreported_by_original_absolute_protocol' in r['flags']:
                r['ranking_status'] = 'pending_conditions'
            else:
                r['ranking_status'] = 'pending_campaign_or_analysis_protocol'
        group_summary.append(dict(group_key=group_id, records=len(members), distinct_sequence_signatures=nseq,
            source_id=members[0]['source_id'], candidate_ids=[r['candidate_id'] for r in members],
            evaluation_ready=False, independent_campaign_confirmed=False))
    assert len(rows) == len(mapping) == 209
    assert len({(r['source_id'], r['candidate_id']) for r in rows}) == 209
    assert coverage(rows)['families'] == 28 and coverage(rows)['substrates'] == 66
    passed = [r for r in rows if r['absolute_status'] == 'carry_forward_pass']
    held = [r for r in rows if r['absolute_status'] == 'pending']
    counts = dict(historical_pool=coverage(rows), absolute_carry_forward=coverage(passed),
        absolute_pending=coverage(held), absolute_excluded=0,
        ranking_status_records=dict(Counter(r['ranking_status'] for r in rows)),
        condition_partitioned_metadata_groups=len(grouped),
        multi_sequence_candidate_groups=sum(g['distinct_sequence_signatures'] >= 2 for g in group_summary),
        multi_sequence_candidate_records=sum(g['records'] for g in group_summary if g['distinct_sequence_signatures'] >= 2),
        confirmed_independent_ranking_campaigns=0, final_evaluation_eligible_records=0)
    counts['repeated_sequence_substrate_blocks'] = len(repeated_pairs)
    counts['repeated_sequence_substrate_records'] = sum(len(g['records']) for g in repeated_pairs)
    counts['cross_source_exact_sequence_substrate_repeat_blocks'] = sum(g['cross_source'] for g in repeated_pairs)
    summary = dict(date='2026-10-05', counts=counts,
        record_gap_if_holds_unresolved=300-len(passed), family_gap_if_holds_unresolved=max(0,30-coverage(passed)['families']),
        classification_basis='fresh integrity checks plus inherited primary-source decisions and documented new holds; not 209 newly repeated full-text audits',
        parent_coverage_note='Known parent labels recorded only for 29 detailed records; sequence signatures are not parent counts. No invented total independent-parent count.',
        model_predictions_generated=False, protocol_changed=False, labels_modified=False,
        limitation='Historical cap retained. No automatic replacement of held rows or global reclustering; subset family counts use fixed original global families.',
        evidence_sha256=evidence, protocol_sha256=sha(protocol_path), original_pool_sha256=sha(GLOBAL/'readiness-audit.json'))
    for ep in [ROOT / 'artifacts/research-scope-20261005-v2/ranking-metadata-feasibility.json',
               BASE / 'campaign-audit-20261005/closeout.json', BASE / 'campaign-audit-complete-20261005/closeout.json',
               BASE / 'campaign-detail-20261005/closeout.json', GLOBAL / 'global-accepted-mapping.csv']:
        evidence[str(ep.relative_to(ROOT))] = sha(ep)
    for name, data in [('records.json', rows), ('groups.json', group_summary), ('sources.json', source_register), ('repeated-pairs.json', repeated_pairs), ('summary.json', summary)]:
        (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n',encoding='utf-8')
    labels = {'carry_forward_pass':'原协议记录资格保留', 'pending':'待核'}
    ranking_labels = {'not_applicable_single_sequence_in_current_pool':'当前池仅一种序列，不构成排序',
        'pending_source_consistency':'原文一致性待核', 'pending_conditions':'实验条件待核',
        'pending_campaign_or_analysis_protocol':'实验组/统计方案待核'}
    lines = ['# 209 条逐记录质量处置清单', '', '绝对值“保留”指沿用已审来源资格并通过本轮完整性检查，不代表正式测试已就绪。原文证据、限制和逐条原因见同目录 records.json 与 sources.json。排序组最低规模尚未冻结；只有一种序列的组在当前池内无法比较排序。', '',
             '| 来源 | 记录 ID | 家族 | 底物 | 绝对值资格 | 排序资格 |', '|---|---|---|---|---|---|']
    for r in rows:
        lines.append('| '+' | '.join([r['source_id'],r['candidate_id'],r['global_family_id'],r['substrate_id'],labels[r['absolute_status']],ranking_labels[r['ranking_status']]])+' |')
    (OUT / 'RECORDS_CN.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k in ['counts','record_gap_if_holds_unresolved','family_gap_if_holds_unresolved']},indent=2))


if __name__ == '__main__':
    main()
