"""Evidence-backed campaign review; retains all original pool decisions."""
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
OUT = BASE / 'campaign-audit-20261005'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    inventory = ROOT / 'artifacts/research-scope-20261005/ranking-metadata-feasibility.json'
    groups = [g for g in json.loads(inventory.read_text())['groups'] if g['distinct_sequence_signatures'] >= 3]
    reviews = {
        'biostudies-S-EPMC12751057': dict(substrate_id='pubchem:447206',
            status='shared_chemical_conditions_assay_design_review_pending',
            ids=['human-dnph1-001', 'human-dnph1-003', 'human-dnph1-004'],
            common_conditions='5hmdUMP hydrolysis; 100 mM HEPES pH 7.0; 25 C; continuous absorbance assay',
            limitations=['Enzyme concentrations, substrate ranges and acquisition durations differ across constructs.',
                        'Concentration normalization is inherent to kcat; differing enzyme dose is not by itself an exclusion. Initial-rate, depletion and fit comparability need source-specific review.',
                        'pH 8.5 TAPS WT measurement must remain separate.'],
            evidence_terms=['Spectrophotometric assay', '100 mM HEPES', '600 s']),
        'europepmc-PMC12838360': dict(substrate_id='pubchem:5962',
            status='reported_common_assay_conditions_verified_subset_only',
            ids=['racx-001', 'racx-017', 'racx-019', 'racx-020'],
            common_conditions='L-lysine racemization; 50 mM Tris-HCl pH 7.5; 37 C; 2.5 ug enzyme in 50 uL; 30 min endpoint; triplicate',
            limitations=['Only the four L-lysine records already retained by the existing family cap are listed; this is not a full campaign selection freeze.',
                        'Endpoint assay linearity is reported up to 60 min; variant-wise raw linearity and uncertainty definitions remain to be checked.',
                        'Reverse D-lysine reactions are not merged into this group.'],
            evidence_terms=['Racemase activity assays', '30-min incubation', 'Kinetic analyses']),
        'europepmc-PMC12362431': dict(substrate_id='pubchem:135398621',
            status='readout_and_product_specificity_review_pending',
            ids=['gfus-001', 'gfus-002', 'gfus-003', 'gfus-004'],
            common_conditions='GDP-keto sugar conversion with NADPH; 37 C; Tris pH 8.0; NADPH-consumption readout',
            limitations=['C116S changes product composition; NADPH consumption is not a product-specific GDP-L-fucose production score.',
                        'Enzyme doses differ; fixed cofactor concentrations and fit comparability require supplemental verification.',
                        'Do not remove a variant based on its response value or silently redefine the target product.'],
            evidence_terms=['altered product specificity', 'GDP-d-altrose', '0.24 and 0.90'])
    }
    evidence = {}
    for source, review in reviews.items():
        pmc = 'PMC' + source.rsplit('PMC', 1)[1]
        xml = BASE / source / 'raw' / (pmc + '-fullText.xml')
        tree = ET.parse(xml)
        paragraphs = [' '.join(''.join(p.itertext()).split()) for p in tree.iter('p')]
        evidence[source] = dict(xml_sha256=sha(xml), paragraphs=[t for t in paragraphs if any(s in t for s in review['evidence_terms'])])
        assert evidence[source]['paragraphs']
        match = [g for g in groups if g['source_id'] == source and g['substrate_id'] == review['substrate_id']]
        assert len(match) == 1
        assert set(review['ids']) == {r['candidate_id'] for r in match[0]['candidates']}
        review['source_csv_sha256'] = sha(BASE / source / 'candidate_records.csv')
        review['evaluation_eligible'] = False
    queue = [dict(source_id=g['source_id'], substrate_id=g['substrate_id'],
                  records=g['records'], status='not_yet_reviewed') for g in groups if g['source_id'] not in reviews]
    assert len(queue) == 18
    (OUT / 'campaign-source-evidence.json').write_text(json.dumps(evidence, indent=2, ensure_ascii=False), encoding='utf-8')
    protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    assert sha(protocol) == 'a9ebdc022d77925ebdf32abd016e4ab09fae6a1ab4421dde7bb675107742e36f'
    pool = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json'
    prior = json.loads((BASE / 'targeted-curation-20261004/closeout.json').read_text())
    assert sha(pool) == prior['global_audit_sha256']
    result = dict(date='2026-10-05', reviewed_groups=3, remaining_groups=18,
                  common_assay_subset_verified=1, evaluation_ready_groups=0,
                  campaign_reviews=reviews, remaining_queue=queue, new_accepted_records=0,
                  source_review=dict(pmc='PMC11213884', status='not_eligible_as_absolute_kcat_from_reviewed_material',
                                    evidence='Supplementary methods use fixed 1 mM acceptor and 2 mM donor, reporting initial specific activities and time courses. No saturating kcat fit verified. Do not equate specific activity times molecular weight with kcat.',
                                    followup='May inform a separate specific-activity task, not the frozen absolute-kcat pool.'),
                  model_predictions_run=False, selections_use_response_magnitude=False,
                  prior_audit_sha256=sha(pool), protocol_sha256=sha(protocol), inventory_sha256=sha(inventory),
                  artifact_sha256={p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'closeout.json'})
    (OUT / 'closeout.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['reviewed_groups', 'remaining_groups', 'common_assay_subset_verified', 'evaluation_ready_groups', 'new_accepted_records']}))


if __name__ == '__main__':
    main()
