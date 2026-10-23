"""Component-level development homology screen for the plant CCD source.

Reads the staged PMC12292525 article XML and OeCCD1/InCCD1 FASTAs under
``artifacts/external/temporal-absolute-kinetics/ccd-prescreen-20261003``,
screens the catalytic components of the reported GST fusions against the
SHA-256-pinned development reference with frozen MMseqs2, and parses Table 2
kinetics into ``source-records.csv`` plus a ``unit-consistency.json`` audit and
``homology-audit.json`` receipt.

Accept/exclude logic: either component with a qualifying development hit
(>=30% identity, >=80% query and target coverage) is
``excluded_component_development_homology`` (a GST fusion must not hide
development homology by diluting whole-chain coverage); hit-free rows are
further ``held_unit_definition_conflict`` because the published nmol/min unit
definition and Table 2 kcat cannot both follow the stated conversion -- values
are audited, never corrected. Governance role: fail-closed prescreen with
``accepted_additions`` 0 and no predictions.
"""
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
OUT=BASE/'ccd-prescreen-20261003'
ARTICLE=BASE/'evidence-sprint-20261003/PMC12292525.xml'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
VERSION='5d152c612b6ad2a56f657b7a02c127eceaea2a75'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def linux(p):return '/mnt/d/'+p.relative_to(Path('D:/')).as_posix()

def main():
    """Screen CCD components, audit unit consistency, and write the homology audit."""
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in protected}
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    text=''.join(ET.parse(ARTICLE).getroot().itertext())
    queries={};mapping={'OeCCD1':('XP_022895849.1',542),'InCCD1':('XP_019150854.1',567)}
    for name,(accession,length) in mapping.items():
        assert accession in text
        file=OUT/(accession+'.fasta')
        s=''.join(file.read_text().splitlines()[1:]);assert len(s)==length
        queries[name]=s
    query=OUT/'catalytic-components.fasta'
    query.write_text(''.join(f'>{name}\n{s}\n' for name,s in queries.items()))
    exe=['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version=subprocess.check_output(exe+['version']).decode().strip();assert version==VERSION
    hits=OUT/'development-hits.tsv'
    cmd=exe+['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'search-tmp'),
        '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
        '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'mmseqs.log').open('wb') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    result={k:[] for k in queries}
    for line in hits.read_text().splitlines():
        r=line.split('\t')
        if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8:result[r[0]].append(r)
    rows=[]
    table=next(n for n in ET.parse(ARTICLE).iter('table-wrap') if n.get('id')=='biology-14-00752-t002')
    for tr in table.iter('tr'):
        c=[''.join(n.itertext()).strip() for n in tr]
        if not c or c[0] not in mapping:continue
        rows.append(dict(enzyme=c[0],article_doi='10.3390/biology14070752',source_location='Table 2 / '+c[0],
            substrate='beta-apo-8-prime-carotenal',accession=mapping[c[0]][0],
            km_mM=c[1],vmax_U_mg=c[2],kcat_s_1=c[3],
            status='excluded_component_development_homology' if result[c[0]] else 'pending_full_construct_review'))
    assert len(rows)==2
    assert '1 nmol' in text and '87.43' in text and '89.55' in text
    unit_checks=[]
    masses={'OeCCD1':87430.0,'InCCD1':89550.0}
    for row in rows:
        nmol_rate=float(row['vmax_U_mg'])*masses[row['enzyme']]/60000000
        unit_checks.append(dict(enzyme=row['enzyme'],mass_Da=masses[row['enzyme']],
            reported_kcat_s_1=float(row['kcat_s_1']),kcat_from_reported_nmol_unit=nmol_rate,
            kcat_if_U_meant_umol_per_min=nmol_rate*1000,
            reported_to_nmol_ratio=float(row['kcat_s_1'])/nmol_rate))
        if not result[row['enzyme']]:row['status']='held_unit_definition_conflict'
    (OUT/'unit-consistency.json').write_text(json.dumps(dict(
        definition_source='Section 2.7: one unit produces 1 nmol beta-ionone per minute',
        mass_source='Section 3.3: GST fusions 87.43 and 89.55 kDa',
        kinetics_source='Table 2: Vmax U/mg and kcat s-1',
        formula='kcat_s-1 = Vmax_nmol_min-1_mg-1 * molecular_mass_Da / 60000000',
        rows=unit_checks,interpretation='Published unit definition and Table 2 kcat cannot both follow the stated conversion. Possible unit-definition typo; no correction inferred or applied.',
        original_values_changed=False),indent=2)+'\n')
    with (OUT/'source-records.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt=dict(date='2026-10-03',doi='10.3390/biology14070752',
        source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12292525/fullTextXML',
        article_sha256=sha(ARTICLE),mmseqs_version=version,command=cmd,reference_sha256=sha(REF),
        query_scope='Catalytic CCD components of reported GST fusions; no claim that these are the full assayed fusion sequences.',
        interpretation='A qualifying catalytic-component hit is evidence against family-cold independence. GST fusion must not hide development homology by diluting whole-chain coverage.',
        results={k:dict(length=len(queries[k]),qualifying_hits=len(v),best_hit=max(v,key=lambda x:float(x[7])) if v else None) for k,v in result.items()},
        accepted_additions=0,excluded_rows=sum(x['status'].startswith('excluded') for x in rows),
        held_rows=sum(x['status']=='held_unit_definition_conflict' for x in rows),
        predictions_generated=False,protected_hashes=before,
        artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='homology-audit.json'})
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    (OUT/'homology-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt['results'],indent=2))

if __name__=='__main__':main()
