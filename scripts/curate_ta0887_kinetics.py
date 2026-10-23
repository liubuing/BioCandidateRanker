"""Curate four primary Ta0887 esterase measurements with deposited processed construct.

Reads staged raw evidence for Europe PMC PMC13498779 (article XML, PDB 10XN
CIF/FASTA, UniProt Q9HJS7 FASTA, PubChem 4-nitrophenyl ester records) from
``artifacts/external/temporal-absolute-kinetics/public-batch-20261004`` into
``europepmc-PMC13498779/raw``, screens the exact processed construct
(GPH + Q9HJS7) against the SHA-256-pinned development reference with frozen
MMseqs2, and writes ``construct_sequences.fasta``, ``candidate_records.csv``,
``provenance.json``, and ``homology-audit.json``.

Accept/exclude logic: the four Table 2 p-NP ester kcat rows (acetate, butyrate,
hexanoate, octanoate) are accepted as ``accepted_homology_cold_pool`` only when
the construct has zero development-corpus hits, the article carries CC-BY-4.0,
and a Vmax/native-mass unit cross-check agrees within 1.5%; any failed marker
assertion aborts the run before outputs are written. Governance role:
fail-closed cold-pool curation -- readiness gate stays False, no model
predictions are generated.
"""
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from biocandidate.data.homology import run_mmseqs_easy_search, run_mmseqs_easy_cluster

BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
BATCH=BASE/'public-batch-20261004'
SOURCE=BASE/'europepmc-PMC13498779'
RAW=SOURCE/'raw'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
DOI='10.1002/2211-5463.70327'
VERSION='5d152c612b6ad2a56f657b7a02c127eceaea2a75'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    """Validate evidence, run the pinned homology screen, and emit accepted records."""
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json', ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in protected}
    RAW.mkdir(parents=True,exist_ok=True)
    names=['PMC13498779.xml','10XN.cif','10XN.fasta','Q9HJS7.fasta']
    compounds=['acetate','butyrate','hexanoate','octanoate']
    names += ['pubchem-4-nitrophenyl-'+n+'.json' for n in compounds]
    for n in names:
        dest=RAW/n
        if dest.exists():assert sha(dest)==sha(BATCH/n)
        else:shutil.copyfile(BATCH/n,dest)
    tree=ET.parse(RAW/'PMC13498779.xml');text=' '.join(tree.getroot().itertext())
    for evidence in [DOI,'10XN','Q9HJS7','remove the 6XHis tag','40','pH 8.0']:
        assert evidence in text,evidence
    licenses=[' '.join(n.itertext()) for n in tree.iter('license')]
    assert any('creativecommons.org/licenses/by/4.0' in x for x in licenses)
    cif=MMCIF2Dict(str(RAW/'10XN.cif'))
    assert DOI in cif['_citation.pdbx_database_id_DOI']
    native=str(SeqIO.read(RAW/'Q9HJS7.fasta','fasta').seq)
    records=list(SeqIO.parse(RAW/'10XN.fasta','fasta'));assert len(records)==1
    construct=str(records[0].seq)
    assert len(native)==235 and construct=='GPH'+native and len(construct)==238
    assert cif['_entity_poly.pdbx_seq_one_letter_code_can'][0].replace('\n','')==construct
    reference={str(r.seq) for r in SeqIO.parse(REF,'fasta')}
    assert native not in reference and construct not in reference
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    fasta=SOURCE/'construct_sequences.fasta'
    fasta.write_text('>ta0887-10xn-processed\n'+construct+'\n')
    exe='wsl:Ubuntu-24.04-D:/home/liubuing/.local/bin/mmseqs'
    observed=subprocess.check_output(['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs','version']).decode().strip()
    assert observed==VERSION
    hit=run_mmseqs_easy_search(exe,fasta,REF,SOURCE/'homology',threads=4)
    cluster=run_mmseqs_easy_cluster(exe,fasta,SOURCE/'homology/family-cluster',threads=4)
    assert not hit.read_text().strip(), 'Construct has a development hit: do not admit'
    component=json.loads((BATCH/'curation-audit.json').read_text(encoding='utf-8'))
    assert component['results']['Ta0887']['qualifying_hits']==0
    table=next(n for n in tree.iter('table-wrap') if n.get('id')=='feb470327-tbl-0002')
    table_text=' '.join(table.itertext())
    rows=[];checks=[]
    for i,(suffix,kcat,sd,km,vmax,cid) in enumerate(zip(compounds,[16.6,13.4,11.1,5.0],[6.1,.5,.5,.4],[16.1,2.5,1.0,1.7],[.056,.045,.037,.017],[13243,75834,567162,97416]),1):
        assert str(kcat) in table_text and str(vmax) in table_text
        data=json.loads((RAW/('pubchem-4-nitrophenyl-'+suffix+'.json')).read_text())['PropertyTable']['Properties'][0]
        assert data['CID']==cid
        mol=Chem.MolFromSmiles(data['SMILES']);assert mol is not None
        assert sum(a.GetAtomicNum()==6 for a in mol.GetAtoms())==[8,10,12,14][i-1]
        # Independent units check using author-reported native mass; retain original labels.
        implied=vmax*26700/(60*1.5)
        assert abs(implied-kcat)/kcat<.015
        checks.append(dict(substrate=suffix,reported_kcat=kcat,from_Vmax_native_mass=implied))
        rows.append(dict(candidate_id=f'ta0887-{i:03d}',article_doi=DOI,stable_record_url='https://europepmc.org/articles/PMC13498779',
                         source_file='raw/PMC13498779.xml',source_row='Table 2 / p-NP '+suffix,
                         organism='Thermoplasma acidophilum',enzyme_identity='Ta0887 esterase',sequence_accession='PDB 10XN entity 1; UniProt Q9HJS7',
                         sequence_id='ta0887-10xn-processed',construct='Deposited processed protein, GPH + Q9HJS7; N-terminal His tag removed',
                         variable_substrate='4-nitrophenyl '+suffix,substrate_pubchem_cid=cid,substrate_isomeric_smiles=data['SMILES'],
                         endpoint='kcat_s-1',**{'kcat_s-1':kcat},kcat_sd_s_1=sd,km_mM=km,assay_pH=8.0,assay_temperature_C=40,
                         replicates=3,error_type='SD',fit_method='nonlinear Michaelis-Menten',
                         quality_caveat='Acetate Km exceeds maximum assayed substrate concentration; fit has larger uncertainty; original reported value retained' if i==1 else '',
                         status_at_normalization='accepted_homology_cold_pool'))
    path=SOURCE/'candidate_records.csv'
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    provenance=dict(source_id=SOURCE.name,article_doi=DOI,article_published='2026-08-21',license='CC-BY-4.0',
                    license_scope='Article Table 2 author-generated measurements',
                    sequence_evidence='Article names PDB 10XN and describes His-tag cleavage before subsequent experiments. Its deposited entity polymer has GPH plus the full Q9HJS7 sequence; the coordinate-observed residues alone are not used.',
                    sequence_inference_boundary='Assayed processed form mapped to deposited 10XN entity from the same paper and purification workflow; complete expression plasmid not independently obtained.',
                    assay_evidence='Biochemical characterization: standard 40 C; kinetics explicitly pH 8.0, overriding default pH 7.0. Table 2 identifies standard assay conditions.',
                    substrate_evidence='Author specifies para-nitrophenyl C2/C4/C6/C8 esters; PubChem names, structures and carbon counts verified.',
                    unit_check=checks,unit_check_note='Published 26.7 kDa mass and 1.5 ug in 1 mL agree with kcat within 1.5%; no label conversion or rescaling.',
                    quality_caveats=['Acetate saturation is incomplete (Km 16.1 mM versus max concentration 10 mM); SD 6.1 s-1 retained. Protocol has no additional saturation threshold.',
                                     'Other primary-substrate conditions contain additives; this is condition-specific kinetics.',
                                     'Publication overlap with development unavailable because development corpus lacks row-level citations.'],
                    exact_sequence_substrate_overlap=0,exact_overlap_basis='Neither native nor full processed sequence appears in frozen complete development FASTA; therefore no exact sequence-substrate pair can overlap.',
                    raw_urls={'PMC13498779.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13498779/fullTextXML','10XN.cif':'https://files.rcsb.org/download/10XN.cif','10XN.fasta':'https://www.rcsb.org/fasta/entry/10XN/display','Q9HJS7.fasta':'https://rest.uniprot.org/uniprotkb/Q9HJS7.fasta'},
                    raw_sha256={n:sha(RAW/n) for n in names},model_predictions_run=False)
    (SOURCE/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    audit=dict(audited_on='2026-10-04',source_id=SOURCE.name,article_doi=DOI,method='MMseqs2 easy-search and linclust; native component independently screened',
               mmseqs_version=observed,min_identity=.3,coverage=.8,coverage_mode=0,development_target_sha256=sha(REF),
               candidate_records=dict(count=4,sha256=sha(path)),construct_sequences_sha256=sha(fasta),homology_hits_sha256=sha(hit),
               family_cluster_sha256=sha(cluster),candidate_mmseqs_families=1,accepted_records=4,accepted_unique_sequences=1,accepted_unique_substrates=4,
               homology_hit_sequences=0,exact_sequence_substrate_overlap=0,status='accepted_homology_cold_pool',
               provenance_sha256=sha(SOURCE/'provenance.json'),component_prescreen_sha256=sha(BATCH/'curation-audit.json'),
               readiness_gate_passes=False,claim_boundary='Four accepted curation-pool records; pending global audit. Not benchmark scores.')
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    (SOURCE/'homology-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(source=SOURCE.name,accepted_records=4,construct_length=len(construct),exact_overlap=0,development_hits=0)))
if __name__=='__main__':main()
