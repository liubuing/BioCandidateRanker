"""Close second source-screen batch; keep screening separate from accepted additions."""
import hashlib,json
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
OUT=BASE/'public-batch2-20261004'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    batch=json.loads((OUT/'batch-register.json').read_text(encoding='utf-8'))
    seqs={r.id:str(r.seq) for r in SeqIO.parse(OUT/'reference-queries.fasta','fasta')}
    hits={k:[] for k in seqs}
    for line in (OUT/'development-hits.tsv').read_text().splitlines():
        r=line.split('\t')
        if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8:hits[r[0]].append(r)
    reasons={
        'PMC13396203':('priority_construct_review','Three direct rows in Table 3. Reference R75estA no hit; synthesized signal-peptide truncation and tag still need exact mapping.'),
        'PMC12896513':('pending_figure_and_sequence_review','Methods describe kinetic fitting; no extracted kcat table. Need original figures and exact constructs.'),
        'PMC12734981':('held_internal_consistency_review','Two direct paragraph kcat values; pNPB reported kcat 9.7 s-1 and Km 2.83 mM imply about 3428 M-1 s-1, whereas paragraph prints 59400. No label changes; investigate primary figure.'),
        'PMC12648627':('priority_mutant_provenance_review','Measured esterase kcat is for published CaCAmut, in min-1. Only deposited WT reference screened; no-hit WT is not a completed mutant screen.'),
        'PMC11118615':('priority_construct_review','Three direct Table 2 rows and no reference hit. Verify retained expression tag, conditions and full assayed sequence.'),
        'PMC11975207':('pending_supplement_and_sequence_review','Supplement Table 2 contains donor kinetics; identify actual varied donor, primary versus cited data and full construct.'),
        'PMC12994016':('priority_ltn3d_construct_review','PlLtn3D deposited sequence has no hit. Multi-substrate Tables 1/3/4/5 include other enzymes, cofactor series, WT repeats and ND entries; must not sum all rows as unique pairs. Units min-1.'),
        'PMC12302338':('held_sequence_download','Six direct Table 2 values; named PDB 5JGY FASTA endpoint returned 404. No homology result for this enzyme.'),
        'PMC10339645':('pending_derived_endpoint_evidence','Extracted Table 4 is prior-study comparison, not new direct kcat. New Vmax measurements require separate primary endpoint and enzyme concentration audit.'),
        'PMC11199699':('wt_homology_excluded_other_entries_unresolved','HhKDFDH deposited WT sequence has qualifying development hits. Other published homologs and variants were not screened; do not exclude or accept the entire paper based on one sequence.'),
        'PMC13358373':('pending_derived_endpoint_evidence','Main table gives Vmax and catalytic efficiencies, including unsaturated cases. Do not multiply rounded values to invent independent absolute kcat labels.'),
        'PMC10572201':('pending_sequence_and_pair_dedup','Eight table values include reaction-direction/co-substrate changes for repeated sequence-substrate pairs. Requires direct sequence check and protocol-consistent deduplication.'),
        'PMC11372007':('priority_mature_construct_review','Two direct apparent kcat values in Table 2; deposited precursor no hit. Verify secreted mature construct and glycosylation/fit caveats.')}
    reviews=[]
    for e in batch['entries']:
        if e.get('license_screen')!='compatible_article_license':
            state,note='held_license','Article license restricted or not resolved; no independently compatible measurement dataset verified.'
        else:state,note=reasons[e['pmcid']]
        reviews.append(dict(pmcid=e['pmcid'],doi=e['doi'],status=state,note=note))
    current=ROOT/'artifacts/external/temporal-global-family-audit-20261004-batch2/readiness-audit.json'
    previous=ROOT/'artifacts/external/temporal-global-family-audit-20261004/readiness-audit.json'
    a=json.loads(current.read_text());b=json.loads(previous.read_text())
    sourceaudit=BASE/'europepmc-PMC12774947/homology-audit.json'
    assert json.loads(sourceaudit.read_text())['accepted_records']==3
    oldclose=json.loads((BASE/'public-batch-20261004/batch-closeout.json').read_text())
    assert sha(previous)==oldclose['global_audit_sha256']
    previous_batch=json.loads((BASE/'public-batch-20261004/batch-register.json').read_text(encoding='utf-8'))
    for rel,h in previous_batch['protected_hashes'].items():assert sha(ROOT/rel)==h
    audit=dict(date='2026-10-04',new_fulltexts=batch['fulltexts'],compatible_article_licenses=batch['compatible_article_licenses'],
               new_reference_queries=len(seqs),reference_queries_with_hit=sum(bool(v) for v in hits.values()),reference_queries_without_hit=sum(not v for v in hits.values()),
               results={k:dict(length=len(seqs[k]),qualifying_hits=len(v),best_hit=max(v,key=lambda r:float(r[7])) if v else None) for k,v in hits.items()},
               source_reviews=reviews,accepted_additions=3,addition_origin='Previously pending AaCSE1 source, not the 20 new fulltexts',new_batch_source_acceptances=0,
               remaining_old_candidate_pairs=3,current_counts=a['counts'],record_delta=a['counts']['records_after_global_family_cap']-b['counts']['records_after_global_family_cap'],
               family_delta=a['counts']['global_families']-b['counts']['global_families'],substrate_delta=a['counts']['substrates_after_global_family_cap']-b['counts']['substrates_after_global_family_cap'],
               remaining_records=101,remaining_families=3,model_predictions_generated=False,
               mmseqs_version=json.loads(sourceaudit.read_text())['mmseqs_version'],min_identity=.3,coverage=.8,coverage_mode=0,
               global_audit=str(current.relative_to(ROOT)),global_audit_sha256=sha(current),previous_global_audit_sha256=sha(previous),
               reference_sha256=sha(ROOT/'artifacts/homology-final/unique_proteins.fasta'),source_audit_sha256=sha(sourceaudit),
               evidence_sha256={name:sha(OUT/name) for name in ['batch-register.json','sequence-register.json','reference-queries.fasta','development-hits.tsv','mmseqs.log','Novagen-pET15b.pdf']})
    assert audit['record_delta']==3 and audit['family_delta']==1 and audit['substrate_delta']==2
    (OUT/'batch-closeout.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if isinstance(v,(int,bool))},indent=2))
if __name__=='__main__':main()
