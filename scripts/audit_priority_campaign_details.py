"""Save source-grounded assay metadata overlays; never modify accepted labels."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
OUT = BASE / 'campaign-detail-20261005'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    prior = json.loads((BASE / 'campaign-audit-complete-20261005/closeout.json').read_text())
    groups = []
    for source, substrate in [('biostudies-S-EPMC12645435', 'pubchem:135398641'),
                              ('europepmc-PMC11659886', 'pubchem:18740504')]:
        original = next(g for g in prior['groups'] if (g['source_id'], g['substrate_id']) == (source, substrate))
        p = BASE / source / 'candidate_records.csv'
        with p.open(encoding='utf-8-sig', newline='') as f:
            rows = [r for r in csv.DictReader(f) if r['candidate_id'] in original['candidate_ids']]
        assert len(rows) == original['records']
        assert sha(p) == original['source_csv_sha256']
        ndt = source.startswith('biostudies')
        members = []
        for r in rows:
            variant = r['variant'] if ndt else r['enzyme_variant']
            parent = 'CtNDT' if ndt else r['sequence_id'].rsplit('-', 1)[0]
            members.append(dict(candidate_id=r['candidate_id'], variant=variant, parent=parent,
                enzyme_concentration_uM=(0.8 if variant == 'Y7A' else 1.0) if ndt else float(r['enzyme_concentration_nM']) / 1000,
                uncertainty_interpretation='fit-derived error; SD not established; main/supp replicate conflict' if ndt else 'author +/- statistic undefined; do not treat as SD or SEM',
                prediction_eligible=False))
        if ndt:
            findings = dict(conditions='inosine varied, adenine 10 mM; pH 7, 65 C; CHES/MES/HEPES',
                time_course_evidence='SI Figure S4 PDF pages 50,52-54 shows eight variant time-course panels; graphical support, not a raw-data re-fit',
                calibration='SI Figure S4 caption: Ado calibration converts HPLC peak area to concentration',
                fitting='Michaelis-Menten in Prism 9; per-variant concentration documented in SI page 7 and Table S3 page 24',
                repeat_status='conflict: main Figure 3 says triplicate, SI page 7 says duplicate',
                error_status='main Figure 3 describes error propagated from fit standard errors; SI Table S3 does not define its +/- as SD',
                source_table_scale_flag='D62N Table S3 kcat 0.002 s^-1 at 1 uM enzyme implies Vmax 0.12 uM/min; Figure S4 page 54 plots about 0.75 uM/min at 1000 uM substrate. Visual inconsistency, no replacement value inferred.',
                evidence_gate='hold_source_consistency_review')
        else:
            findings = dict(conditions='4MU-6S-GlcNAc, pH 6, 25 C; Bis-tris propane/citrate/NaCl',
                initial_rate_evidence='main Michaelis-Menten Kinetics explicitly says continuous initial release rates',
                calibration='4MU standard curve measured in same buffer',
                fitting='main formula explicitly includes Et*kcat*x/(Km+x), Prism 6.01',
                repeat_status='triplicate documented; biological versus technical replication not specified in reviewed kinetics section',
                error_status='no SD/SEM definition located for reported kinetic +/-; do not borrow definitions from other experiments',
                parent_blocks=dict(Counter(m['parent'] for m in members)),
                range_flags=['SI Fig S2 panel b endpoint around 3000 uM but caption says 300 uM',
                             'SI Fig S8 panel e extends to about 1500 uM but caption says 300 uM',
                             'SI Fig S7 panel b extends to 500 uM but caption says 300 uM'],
                evidence_gate='initial_rate_and_normalization_supported_uncertainty_and_range_qa_pending')
        groups.append(dict(source_id=source, substrate_id=substrate, records=len(rows), members=members,
                           source_csv_sha256=sha(p), findings=findings, evaluation_ready=False))
    assert [g['records'] for g in groups] == [8, 21]
    assert groups[1]['findings']['parent_blocks'] == {'wg': 7, 'yg': 7, 'f3-orf26': 7}
    pool = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json'
    protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    assert sha(pool) == prior['prior_pool_sha256']
    assert sha(protocol) == prior['protocol_sha256']
    evidence = [BASE / 'biostudies-S-EPMC12645435/raw/cb5c00120_si_001.pdf',
                BASE / 'europepmc-PMC11659886/raw/bg4c00088_si_001.pdf',
                BASE / 'europepmc-PMC11659886/raw/PMC11659886-fullText.xml',
                OUT / 'PMC12645435.xml', OUT / 'amuc-fulltext.xml', OUT / 'amuc-supp.docx']
    result = dict(date='2026-10-05', audited_records=29, audited_groups=2, groups=groups,
        accepted_additions=0, evaluation_ready_groups=0, model_predictions_generated=False,
        labels_modified=False, label_magnitude_used_for_selection=False,
        numerical_consistency_inspected_for_quality_control=True,
        warning='Historical accepted-pool count is unchanged; it is not an assertion that all ranking metadata and source consistency issues are closed.',
        existing_pool_counts=prior['pool_counts'], pool_sha256=sha(pool), protocol_sha256=sha(protocol),
        evidence_sha256={str(p.relative_to(ROOT)): sha(p) for p in evidence})
    OUT.mkdir(exist_ok=True)
    (OUT / 'closeout.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: result[k] for k in ['audited_records','audited_groups','accepted_additions','evaluation_ready_groups','model_predictions_generated']}))


if __name__ == '__main__':
    main()
