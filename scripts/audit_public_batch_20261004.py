"""Record manually checked source cells and component-screen dispositions."""
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/external/temporal-absolute-kinetics/public-batch-20261004'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((OUT/name).read_text(encoding='utf-8'))
def main():
    batch=read('batch-register.json')
    for rel,h in batch['protected_hashes'].items():assert sha(ROOT/rel)==h
    ref=ROOT/'artifacts/homology-final/unique_proteins.fasta'
    assert sha(ref)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    queries={r.id:str(r.seq) for r in SeqIO.parse(OUT/'component-queries.fasta','fasta')}
    hits={k:[] for k in queries}
    for line in (OUT/'development-hits.tsv').read_text().splitlines():
        r=line.split('\t')
        if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8:hits[r[0]].append(r)
    result={k:dict(length=len(queries[k]),qualifying_hits=len(v),best_hit=max(v,key=lambda r:float(r[7])) if v else None) for k,v in hits.items()}
    # Exact source values; no estimates read off plotted curves and no computed labels.
    cells=[
        ('PMC12774947','AaCSE1','caffeoyl-5-O-shikimic acid',.0066,'Table 1 / AaCSE1 / nonlinear fit','primary',None,None),
        ('PMC12774947','AaCSE1','caffeoyl-5-O-shikimic acid',.0062,'Table 1 / AaCSE1 / ECB-method','alternative_fit_not_independent',None,None),
        ('PMC12774947','AaCSE1','caffeoyl-5-O-quinic acid',.0159,'Table 1 / AaCSE1 / chlorogenic acid','primary',None,None),
        ('PMC12774947','AaCSE1','4-nitrophenyl butyrate',.5233,'Table 1 / AaCSE1 / 4-NPB','primary',None,None),
        ('PMC12774947','MeMAGL3','4-nitrophenyl butyrate',11.5,'Table 1 / MeMAGL3 / 4-NPB','apparent_parameter_requires_review',None,None),
        ('PMC12405362','CiCQE1','5-O-caffeoylquinic acid',.256,'Section 3 / Figure 6A accompanying paragraph','primary',None,None),
        ('PMC12405362','CiCQE3','5-O-caffeoylquinic acid',.466,'Section 3 / Figure 6A accompanying paragraph','primary',None,None),
    ]
    for name,values in [('LbXyn43B',[349.6,247.2,150.0]),('LbXyn43B_T274A',[317.4,199.8,127.6])]:
        for substrate,value in zip(['beta-1,4-D-xylobiose','beta-1,4-D-xylotriose','beta-1,4-D-xylotetraose'],values):
            cells.append(('PMC12432433',name,substrate,value,'Table 4 / Present work / '+name+' / '+substrate,'primary',37,7.0))
    for substrate,value in [('p-nitrophenyl acetate',16.6),('p-nitrophenyl butyrate',13.4),('p-nitrophenyl hexanoate',11.1),('p-nitrophenyl octanoate',5.0)]:
        cells.append(('PMC13498779','Ta0887',substrate,value,'Table 2 / '+substrate,'primary',None,8.0))
    source={x['pmcid']:x for x in batch['entries'] if x.get('pmcid')}
    rows=[]
    for pmc,name,substrate,value,location,kind,temp,ph in cells:
        assert source[pmc]['license_screen']=='compatible_article_license'
        row=dict(candidate_id=f'public-batch-20261004-{len(rows)+1:03d}',pmcid=pmc,doi=source[pmc]['doi'],
                 enzyme=name,substrate=substrate,kcat_s_1=value,source_location=location,measurement_role=kind,
                 assay_temperature_C=temp,assay_pH=ph,source_sha256=sha(OUT/(pmc+'.xml')),
                 status='excluded_component_development_homology' if hits[name] else 'pending_construct_substrate_conditions_and_global_audit')
        if kind=='alternative_fit_not_independent':row['status']='secondary_fit_do_not_count_independently'
        if kind=='apparent_parameter_requires_review':row['status']='held_apparent_fit_and_construct_review'
        rows.append(row)
    with (OUT/'staged-measurements.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    reviews=[]
    for e in batch['entries']:
        if 'fulltext_sha256' not in e:continue
        pmc=e['pmcid'];article=ET.parse(OUT/(pmc+'.xml')).getroot()
        if article.get('article-type')=='review-article':state='not_primary_measurement_article'
        elif e['license_screen']!='compatible_article_license':state='held_article_license_no_separate_dataset_license_verified'
        elif pmc in {'PMC12405362','PMC12774947','PMC13498779'}:state='priority_no_component_hit_pending_remaining_gates'
        elif pmc in {'PMC12432433','PMC13249509'}:state='tested_kinetic_components_have_development_hits'
        elif pmc=='PMC13533970':state='not_selected_for_benign_industrial_plant_batch'
        elif pmc=='PMC13202482':state='table_is_prior_study_kcat_over_km_not_new_absolute_kcat'
        elif pmc=='PMC13304128':state='three_direct_rows_identified_but_reference_accession_not_yet_verified_as_assayed_sequence'
        else:state='pending_endpoint_or_supplement_manual_review'
        reviews.append(dict(pmcid=pmc,doi=e['doi'],title=e['title'],status=state))
    receipt=dict(date='2026-10-04',metadata_records=batch['unique_metadata_records'],fulltexts=batch['fulltexts_acquired'],
                 compatible_article_licenses=batch['compatible_article_licenses'],component_queries=len(queries),
                 component_hit_queries=sum(bool(x) for x in hits.values()),component_no_hit_queries=sum(not x for x in hits.values()),
                 source_cell_records=len(rows),nonindependent_alternative_fits=1,homology_excluded_cell_records=sum(r['status']=='excluded_component_development_homology' for r in rows),
                 pending_unique_pairs=sum(not hits[r['enzyme']] and r['measurement_role']!='alternative_fit_not_independent' for r in rows),
                 accepted_additions=0,accepted_global_family_additions=0,predictions_generated=False,
                 mmseqs_version='5d152c612b6ad2a56f657b7a02c127eceaea2a75',
                 mmseqs_arguments=['easy-search','component-queries.fasta',str(ref),'development-hits.tsv','mmseqs-tmp','--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4','--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits'],
                 reference_sha256=sha(ref),results=result,source_reviews=reviews,
                 limitations=['No-hit catalytic components are not fully audited assayed constructs or independent global families.',
                              'Blank staged assay conditions mean not yet extracted; they do not assert that the paper failed to report conditions.',
                              'No structural identifiers, complete source gate audit or global family update has been completed for the pending cells.',
                              'Broad full-text queries retrieve reference mentions and reviews. Metadata counts are not a measurement count.',
                              'AaCSE1 and MeMAGL3 paper lengths 324/344 differ from deposited translated CDS lengths 323/343; stop codons present; keep source discrepancy visible.',
                              'MeMAGL3 is described as an apparent fit with inhibition/solubility limitations; requires separate quality decision.'],
                 protected_hashes=batch['protected_hashes'],artifacts={n:sha(OUT/n) for n in ['batch-register.json','sequence-register.json','component-queries.fasta','development-hits.tsv','mmseqs.log','staged-measurements.csv']})
    (OUT/'curation-audit.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if isinstance(v,(int,bool))},indent=2))
if __name__=='__main__':main()
