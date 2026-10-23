"""Curate three directly measured AaCSE1 kcat rows from DOI 10.1007/s00425-025-04915-7
(PMC12774947).

Evidence written under artifacts/external/temporal-absolute-kinetics/europepmc-PMC12774947:
candidate_records.csv (3 accepted rows), provenance.json (caveats and withheld fits),
homology-audit.json, and the reconstructed construct_sequences.fasta.

Acceptance requires a documented vector reconstruction (deposited PV962689.1 CDS, exact
primer-end matches, Novagen pET-15b leader) and a homology-cold pool: MMseqs2 easy-search
against the frozen development reference must return zero hits. The alternative ECB fit of
the same shikimate measurements is withheld as non-independent evidence.
"""
import csv, hashlib, json, shutil, subprocess, sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from biocandidate.data.homology import run_mmseqs_easy_search,run_mmseqs_easy_cluster
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
OLD=BASE/'public-batch-20261004';NEW=BASE/'public-batch2-20261004'
SOURCE=BASE/'europepmc-PMC12774947';RAW=SOURCE/'raw'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
DOI='10.1007/s00425-025-04915-7'
VERSION='5d152c612b6ad2a56f657b7a02c127eceaea2a75'
def sha(p):
    """Return the SHA-256 hex digest of the file at ``p``."""
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    """Reconstruct the AaCSE1 construct, audit homology, and write the evidence bundle
    fail-closed."""
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit-20261004/readiness-audit.json']
    before={str(p):sha(p) for p in protected}
    RAW.mkdir(parents=True,exist_ok=True)
    for parent,names in [(OLD,['PMC12774947.xml','PV962689.gb','AaCSE-supplement.pdf','AaCSE-supplement.txt','pubchem-4-nitrophenyl-butyrate.json']),
                         (NEW,['Novagen-pET15b.pdf','pubchem-327-97-9.json','pubchem-73263-62-4.json'])]:
        for name in names:
            dest=RAW/name
            if dest.exists():assert sha(dest)==sha(parent/name)
            else:shutil.copyfile(parent/name,dest)
    article=(RAW/'PMC12774947.xml').read_text(encoding='utf-8')
    supplement=(RAW/'AaCSE-supplement.txt').read_text(encoding='utf-8')
    for token in [DOI,'pET-15b','NdeI','XhoI','PV962689','creativecommons.org/licenses/by/4.0']:assert token in article
    rec=SeqIO.read(RAW/'PV962689.gb','genbank');dna=str(rec.seq)
    feature=next(x for x in rec.features if x.type=='CDS');core=feature.qualifiers['translation'][0]
    assert len(core)==323 and len(dna)==972 and str(Seq(dna).translate())==core+'*'
    fwd='ATGGCCTCCGCAGCGCCGCC';rev='TCACGACTTGATCTGGCTGTCATGCTTCTCC'
    assert fwd in supplement and rev in supplement and dna.startswith(fwd) and dna.endswith(str(Seq(rev).reverse_complement()))
    # Manufacturer TB045 coding strand upstream of the insert ATG, visually checked.
    leader_dna='ATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCAT'
    leader=str(Seq(leader_dna).translate());assert leader=='MGSSHHHHHHSSGLVPRGSH'
    construct=leader+core;assert len(construct)==343
    refseq={str(r.seq) for r in SeqIO.parse(REF,'fasta')}
    assert core not in refseq and construct not in refseq
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    fasta=SOURCE/'construct_sequences.fasta';fasta.write_text('>aacse1-pet15b-reconstructed\n'+construct+'\n')
    # Frozen WSL mmseqs binary: verify the pinned version, then hand the homology helpers a
    # "wsl:<distro>:<path>" handle so they translate Windows artifact paths to /mnt/... form.
    version=subprocess.check_output(['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs','version']).decode().strip();assert version==VERSION
    exe='wsl:Ubuntu-24.04-D:/home/liubuing/.local/bin/mmseqs'
    hits=run_mmseqs_easy_search(exe,fasta,REF,SOURCE/'homology',threads=4)
    clusters=run_mmseqs_easy_cluster(exe,fasta,SOURCE/'homology/family-cluster',threads=4)
    assert not hits.read_text().strip()
    component=json.loads((OLD/'curation-audit.json').read_text(encoding='utf-8'));assert component['results']['AaCSE1']['qualifying_hits']==0
    rows=[]
    specs=[('5-O-caffeoylshikimic acid','73263-62-4',5281762,.0066,.0005,36.4,6.5,'6 in main Table 1; 9 in Supplement Table S6'),
           ('chlorogenic acid','327-97-9',1794427,.0159,.0009,425.1,6.5,'9'),
           ('4-nitrophenyl butyrate',None,75834,.5233,None,662.1,7.0,'3')]
    for i,(name,cas,cid,kcat,sd,km,ph,n) in enumerate(specs,1):
        f=RAW/('pubchem-'+cas+'.json' if cas else 'pubchem-4-nitrophenyl-butyrate.json')
        c=json.loads(f.read_text())['PropertyTable']['Properties'][0];assert c['CID']==cid and Chem.MolFromSmiles(c['SMILES']) is not None
        if cas:assert cas in supplement
        assert str(kcat) in article
        rows.append(dict(candidate_id=f'aacse1-{i:03d}',article_doi=DOI,stable_record_url='https://europepmc.org/articles/PMC12774947',
                         source_file='raw/PMC12774947.xml; raw/AaCSE-supplement.pdf',source_row='Table 1 / AaCSE1 / '+name,
                         organism='Anthoceros agrestis',enzyme_identity='AaCSE1 caffeoylshikimate esterase',
                         sequence_accession='PV962689.1; XXS03338.1',sequence_id='aacse1-pet15b-reconstructed',
                         construct='Reconstructed from deposited CDS, verified terminal primers and Novagen pET-15b TB045; not a deposited full plasmid',
                         variable_substrate=name,substrate_pubchem_cid=cid,substrate_isomeric_smiles=c['SMILES'],endpoint='kcat_s-1',
                         **{'kcat_s-1':kcat},kcat_sd_s_1=sd,km_uM=km,assay_pH=ph,assay_temperature_C=25,replicates=n,
                         error_type='SD where reported',fit_method='Michaelis-Menten; alternative ECB fit not counted independently',
                         status_at_normalization='accepted_homology_cold_pool'))
    csvfile=SOURCE/'candidate_records.csv'
    with csvfile.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    provenance=dict(source_id=SOURCE.name,article_doi=DOI,article_published='2026-01-06',license='CC-BY-4.0',license_scope='Article Table 1 primary measurements',
                    sequence_evidence='Author-deposited complete CDS; Supplement Table S4 primers exactly match both ends; article states unmodified pET-15b, NdeI/XhoI, retained N-terminal His tag. Manufacturer TB045 defines translated leader. Reverse primer retains stop, so no C-terminal vector extension.',
                    sequence_evidence_level='documented reconstruction from accession and explicit cloning description; not independent sequencing of the physical experimental plasmid',
                    prior_hold_resolution='Earlier leader was only a hypothesis. Original manufacturer map plus exact primer-end checks now supports the documented reconstruction. Frozen protocol allows full sequence or unambiguous accession and variant; it does not require depositing the entire plasmid.',
                    manufacturer_document_url='https://mcb.berkeley.edu/labs/krantz/pdf/vectors/pet15b_vector.pdf',
                    paper_length_discrepancy='972 bases comprise 323 amino acids and a stop; deposited translation used. Paper prints 324 amino acids.',
                    computed_mass_Da=ProteinAnalysis(construct).molecular_weight(),paper_mass_Da=38040,
                    quality_caveats=['Shikimate saturation incomplete; main n=6 and supplemental n=9 disagree; both retained.',
                                     'pH is reported reaction-buffer pH, not pH 7 of substrate stock. All assays 25 C.',
                                     'Development row-level citation overlap is unavailable.'],
                    withheld=[dict(enzyme='MeMAGL3',reason='Apparent fit with inhibition/solubility limitations still pending quality review; not included in this source CSV.'),
                              dict(enzyme='AaCSE1',value=.0062,reason='Alternative ECB fit of same shikimate measurements; not independent.')],
                    raw_sha256={p.name:sha(p) for p in RAW.iterdir() if p.is_file()},exact_sequence_substrate_overlap=0,model_predictions_run=False)
    (SOURCE/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    audit=dict(audited_on='2026-10-04',source_id=SOURCE.name,article_doi=DOI,method='MMseqs2 easy-search and linclust; independent native-component prescreen',
               mmseqs_version=version,min_identity=.3,coverage=.8,coverage_mode=0,development_target_sha256=sha(REF),
               candidate_records=dict(count=3,sha256=sha(csvfile)),construct_sequences_sha256=sha(fasta),homology_hits_sha256=sha(hits),
               family_cluster_sha256=sha(clusters),candidate_mmseqs_families=1,accepted_records=3,accepted_unique_sequences=1,accepted_unique_substrates=3,
               exact_sequence_substrate_overlap=0,status='accepted_homology_cold_pool',provenance_sha256=sha(SOURCE/'provenance.json'),
               readiness_gate_passes=False,claim_boundary='Accepted curation pool only; requires global cap audit; no model scores.')
    assert before=={str(p):sha(p) for p in protected}
    (SOURCE/'homology-audit.json').write_text(json.dumps(audit,indent=2)+'\n');print('Accepted source records: 3')
if __name__=='__main__':main()
