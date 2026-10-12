"""Reconstruct the reported standard-vector interpretation for an evidence audit.

This does not claim independent sequencing of the experimental plasmid.
"""
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
PRE=BASE/'esterase-prescreen-20261003'
OUT=BASE/'estsit01-curation-20261003'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
VECTOR_URL='https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/product/documents/104/753/tb036vm.pdf'
# Manufacturer TB036, page 1, coding strand; transcribed from the published map.
VECTOR='ATGGCTAGCATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGACAAGCTTGCGGCCGCACTCGAGCACCACCACCACCACCACTGA'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fasta(p):return ''.join(p.read_text().splitlines()[1:])
def linux(p):return '/mnt/d/'+p.relative_to(Path('D:/')).as_posix()

def main():
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in protected}
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    cds=fasta(OUT/'EstSIT01-cds.fasta');native=fasta(PRE/'ANG84115.1.fasta')
    assert len(cds)==1113 and str(Seq(cds).translate())==native+'*'
    # Original supplement Table S1 primer annealing regions.
    assert cds.startswith('ATGACCCTGTTCGACGGCATCACGTCT')
    reverse='GTCGGCGGAGCGGATGATGATCGCCTC'
    assert cds[:-3].endswith(str(Seq(reverse).reverse_complement()))
    prefix=VECTOR[:VECTOR.index('GAATTC')+6]
    suffix=VECTOR[VECTOR.index('AAGCTT'):]
    full_dna=prefix+cds[:-3]+suffix
    full=str(Seq(full_dna).translate())
    assert full.endswith('*') and '*' not in full[:-1]
    ntag=str(Seq(prefix).translate());ctag=str(Seq(suffix).translate()).rstrip('*')
    assert full[:-1]==ntag+native+ctag
    sequences={}
    for name in ['WT','E301P','E301P/G215A','E301P/G215A/G293A']:
        core=list(native)
        for change in ([] if name=='WT' else name.split('/')):
            ix=int(change[1:-1])-1;assert core[ix]==change[0];core[ix]=change[-1]
        sequences[name.replace('/','_')]=ntag+''.join(core)+ctag
    query=OUT/'article-inferred-constructs.fasta'
    query.write_text(''.join(f'>{k}\n{v}\n' for k,v in sequences.items()))
    (OUT/'vector-map-transcription.json').write_text(json.dumps(dict(source=VECTOR_URL,page=1,coding_region=VECTOR,method='manual transcription from manufacturer map; not a deposited experimental plasmid'),indent=2)+'\n')
    exe=['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version=subprocess.check_output(exe+['version']).decode().strip()
    assert version=='5d152c612b6ad2a56f657b7a02c127eceaea2a75'
    hits=OUT/'construct-hits.tsv'
    cmd=exe+['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'search-tmp'),
             '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
             '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'mmseqs.log').open('wb') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    result={k:[] for k in sequences}
    for line in hits.read_text().splitlines():
        r=line.split('\t');result[r[0]].append(r)
    refseqs={''.join(part.splitlines()[1:]) for part in REF.read_text().split('>') if part}
    compound=json.loads((OUT/'pubchem-11234411.json').read_text())['PropertyTable']['Properties'][0]
    assert compound['CID']==11234411 and compound['MolecularFormula']=='C21H22N2O5'
    table=next(n for n in ET.parse(PRE/'PMC11676285.xml').iter('table-wrap') if 'Kinetic parameters and specific activities' in ''.join(n.itertext()))
    consistency=[]
    for tr in table.iter('tr'):
        c=[''.join(n.itertext()).strip() for n in tr]
        if not c or c[0] not in ['Wild-type','E301P','E301P/G215A','E301P/G215A/G293A']:continue
        name='WT' if c[0]=='Wild-type' else c[0].replace('/','_')
        km=float(c[1].split('±')[0]);kcat=float(c[2].split('±')[0]);reported=float(c[4].split('±')[0])
        mass=ProteinAnalysis(sequences[name]).molecular_weight()
        expected=kcat*60000/mass*10/(10+km)
        consistency.append(dict(variant=name,molecular_mass_Da=mass,reported_specific_activity_U_mg=reported,
                                expected_U_mg_at_10mM=expected,ratio_expected_to_reported=expected/reported))
    (OUT/'measurement-consistency.json').write_text(json.dumps(dict(
        source='PMC11676285 Table 3, both footnotes: pH 8.0, 30 C, 0.05 M PBS; specific activity at 10 mM substrate',
        formula='U/mg = kcat_s-1 * 60000 / molecular_mass_Da * S_mM/(Km_mM+S_mM)',
        assumptions='Same preparation and active fraction; reported Michaelis-Menten model; mass normalization to inferred monomer. Differences may reflect assay/preparation or normalization details; this calculation alone does not establish an error.',
        rows=consistency,labels_modified=False),indent=2)+'\n')
    with (PRE/'candidate-measurements.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    for row in rows:
        row['substrate_pubchem_cid']=11234411
        row['substrate_smiles']=compound['SMILES']
        row['duplicate_group']='EstSIT01-WT-CID11234411' if row['variant']=='WT' else ''
        row['duplicate_resolution']='retain_earliest_publication_as_primary' if row['variant']=='WT' and row['source']=='PMC11180644' else ('secondary_repeat_not_independent' if row['variant']=='WT' else 'unique_variant')
        row['status']='held_for_measurement_consistency_review'
        row['construct_interpretation']='standard pET21a EcoRI/HindIII reconstruction; not independently deposited'
    with (OUT/'reviewed-measurements.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    audit=dict(date='2026-10-03',status='held_for_measurement_consistency_review',
        accepted_additions=0,measurements=5,independent_variant_substrate_pairs=4,
        n_terminal_extension=ntag,c_terminal_extension=ctag,inferred_length=len(full)-1,
        construct_evidence='Original Table S1 primer ends match CP015810.1:135837..136949; published EcoRI/HindIII and standard manufacturer vector map give an in-frame fusion.',
        inference_limit='Actual experimental plasmid sequence not deposited; standard unmodified vector interpretation only.',
        duplicate_policy='Within this source group, earliest publication is the primary WT record regardless of kcat; later WT remains a repeat and is not independently counted.',
        compound=compound,mmseqs_version=version,command=cmd,homology_hits=result,
        exact_construct_sequence_overlap={k:v in refseqs for k,v in sequences.items()},
        remaining=['Reconcile Table 3 turnover numbers with reported mass-specific activities under the same assay conditions; do not silently rescale labels.',
                   'Confirm substrate structure against original scheme (CID/name match acquired, scheme graphic not yet inspected).',
                   'Global family mapping/cap only after source-level evidence review passes.'],
        predictions_generated=False,protected_hashes=before,
        artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='construct-review.json'})
    (OUT/'construct-review.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(inferred_construct_length=len(full)-1,hit_counts={k:len(v) for k,v in result.items()},accepted_additions=0)))

if __name__=='__main__':main()
