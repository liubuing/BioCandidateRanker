"""Audit published PlLtn3D measurements; no predictions or protocol changes."""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from biocandidate.data.homology import run_mmseqs_easy_search, run_mmseqs_easy_cluster

BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
STAGE = BASE / 'targeted-curation-20261004'
SOURCE = BASE / 'europepmc-PMC12994016'
RAW = SOURCE / 'raw'
REF = ROOT / 'artifacts/homology-final/unique_proteins.fasta'
DOI = '10.1016/j.jbc.2026.111280'
VERSION = '5d152c612b6ad2a56f657b7a02c127eceaea2a75'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    protected = [ROOT / 'configs/temporal_absolute_kinetics_protocol.json',
                 ROOT / 'artifacts/external/temporal-global-family-audit-20261004-batch2/readiness-audit.json']
    before = {str(p): sha(p) for p in protected}
    RAW.mkdir(parents=True, exist_ok=True)
    inputs = [BASE / 'public-batch2-20261004/PMC12994016.xml', STAGE / '9X6I.cif',
              STAGE / 'A0A844HLS7.fasta', *sorted(STAGE.glob('pubchem-*.json'))]
    for p in inputs:
        dest = RAW / p.name
        if dest.exists():
            assert sha(dest) == sha(p)
        else:
            shutil.copyfile(p, dest)
    article = (RAW / 'PMC12994016.xml').read_text(encoding='utf-8')
    tree = ET.fromstring(article)
    for token in [DOI, 'creativecommons.org/licenses/by/4.0', 'MRGSHHHHHHG']:
        assert token in article
    date = tree.find('.//pub-date')
    published = '-'.join(date.findtext(x).zfill(2) for x in ['year', 'month', 'day'])
    assert published == '2026-02-12'
    cif = MMCIF2Dict(str(RAW / '9X6I.cif'))
    assert DOI in cif['_citation.pdbx_database_id_DOI']
    assert 'A0A844HLS7' in cif['_struct_ref.pdbx_db_accession']
    deposited = ''.join(cif['_entity_poly.pdbx_seq_one_letter_code_can'][0].split())
    native = str(SeqIO.read(RAW / 'A0A844HLS7.fasta', 'fasta').seq)
    assert len(native) == 244 and len(deposited) == 255
    assert deposited == 'MRGSHHHHHHG' + 'S' + native[1:]
    constructs = {'plltn3d-WT': deposited}
    for mutation in ['S136A', 'S138A', 'M186A', 'M238A']:
        index = int(mutation[1:-1]) + 10  # 11-residue leader; native positions are 1-based.
        assert deposited[index] == mutation[0]
        constructs['plltn3d-' + mutation] = deposited[:index] + mutation[-1] + deposited[index + 1:]
    fasta = SOURCE / 'construct_sequences.fasta'
    fasta.write_text(''.join(f'>{key}\n{seq}\n' for key, seq in constructs.items()))
    cores = SOURCE / 'core_sequences.fasta'
    cores.write_text(''.join(f'>{key}-core\n{seq[11:]}\n' for key, seq in constructs.items()) + f'>plltn3d-reference-native\n{native}\n')
    reference = {str(r.seq) for r in SeqIO.parse(REF, 'fasta')}
    assert not reference.intersection([*constructs.values(), *[s[11:] for s in constructs.values()], native])
    assert sha(REF) == '69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    version = subprocess.check_output(['wsl', '-d', 'Ubuntu-24.04-D', '--', '/home/liubuing/.local/bin/mmseqs', 'version']).decode().strip()
    assert version == VERSION
    exe = 'wsl:Ubuntu-24.04-D:/home/liubuing/.local/bin/mmseqs'
    hits = run_mmseqs_easy_search(exe, fasta, REF, SOURCE / 'homology', threads=4)
    corehits = run_mmseqs_easy_search(exe, cores, REF, SOURCE / 'homology/core', threads=4)
    clusters = run_mmseqs_easy_cluster(exe, fasta, SOURCE / 'homology/family-cluster', threads=4)
    assert not hits.read_text().strip() and not corehits.read_text().strip()
    # Table order fixes selection without using response magnitudes. NAD+ duplicate excluded.
    specs = [('WT', 'L-threonate', 'L-threonic-acid', 5460407, '198', '4', 6.50, 'tbl1', 0),
             ('WT', 'D-fuconate', 'D-fuconic-acid', 5459874, '489', '29', 6.09, 'tbl1', 2),
             ('WT', 'L-arabinonate', 'L-arabinonic-acid', 6602426, '614', '34', 9.98, 'tbl1', 3),
             ('WT', 'D-xylonate', 'D-xylonic-acid', 6602431, '104', '1', 1.96, 'tbl1', 5),
             ('WT', 'D-galactonate', 'D-galactonic-acid', 128869, '883', '18', 50.8, 'tbl1', 6),
             ('WT', 'D-glycerate', 'D-glyceric-acid', 439194, '7.52', '0.26', 6.01, 'tbl1', 7),
             ('S136A', 'L-threonate', 'L-threonic-acid', 5460407, '0.150', '0.002', 6.16, 'tbl4', 1),
             ('S138A', 'L-threonate', 'L-threonic-acid', 5460407, '200', '5', 6.36, 'tbl4', 2),
             ('M186A', 'L-threonate', 'L-threonic-acid', 5460407, '269', '8', 12.3, 'tbl4', 6),
             ('M238A', 'L-threonate', 'L-threonic-acid', 5460407, '513', '3', 125, 'tbl4', 8)]
    rows = []
    for i, (variant, name, chemname, cid, value, sd, km, table, rowindex) in enumerate(specs, 1):
        tab = tree.find(f'.//table-wrap[@id="{table}"]')
        cells = [''.join(c.itertext()) for c in tab.findall('.//tbody/tr')[rowindex].findall('td')]
        assert cells[-2] == f'{value} ± {sd}'
        compound = json.loads((RAW / f'pubchem-{chemname}.json').read_text())['PropertyTable']['Properties'][0]
        assert compound['CID'] == cid and Chem.MolFromSmiles(compound['SMILES']) is not None
        rows.append(dict(candidate_id=f'plltn3d-{i:03d}', article_doi=DOI,
                         stable_record_url='https://europepmc.org/articles/PMC12994016',
                         source_file='raw/PMC12994016.xml', source_row=f'{table} / body row {rowindex + 1} / kcat',
                         organism='Paracoccus litorisediminis', enzyme_identity='PlLtn3D GL300_RS07945',
                         sequence_accession='PDB 9X6I; UniProt A0A844HLS7', sequence_id='plltn3d-' + variant,
                         construct='Deposited full entity; N-terminal leader and reference M1S difference preserved; published variant positions use native numbering',
                         variant=variant, variable_substrate=name, substrate_pubchem_cid=cid,
                         substrate_isomeric_smiles=compound['SMILES'], structure_normalization='Neutral parent acid of named sugar-acid anion',
                         endpoint='kcat_s-1', **{'kcat_s-1': float(value) / 60}, kcat_sd_s_1=float(sd) / 60,
                         reported_kcat=value, reported_sd=sd, reported_unit='min^-1', conversion='divide value and SD by 60',
                         km_mM=km, assay_pH=10, assay_temperature_C='', temperature_evidence='Not reported for kinetic assay; chromatography temperature not substituted',
                         cofactor='NADP+', replicates=3, error_type='SD', fit_method='Michaelis-Menten',
                         quality_caveat='Some fitted Km exceed measured concentration range; see table footnotes',
                         status_at_normalization='accepted_homology_cold_pool'))
    csvfile = SOURCE / 'candidate_records.csv'
    with csvfile.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    provenance = dict(source_id=SOURCE.name, article_doi=DOI, article_published=published, license='CC-BY-4.0',
                      license_scope='Article primary Tables 1 and 4',
                      sequence_evidence='Same-paper PDB 9X6I complete entity, article states retained 11-residue N-terminal leader. Mutants were made from same plasmid and prepared by same procedures. Four published substitutions checked against native numbering.',
                      sequence_caveat='Deposited core starts S whereas UniProt starts M; deposited full sequence used. This is a deposited construct/reference discrepancy, not an invented correction or independently sequenced assay plasmid.',
                      selection_rule='Primary carbon-substrate measurements in table order; first cofactor condition for repeated WT pair; only numeric kcat; no model or magnitude-based selection.',
                      withheld=[{'item': 'Table 1 g-Gluconate', 'reason': 'Nonstandard stereochemical prefix conflicts with narrative D-gluconate; hold identity mapping.'},
                                {'item': 'Table 1 D-glyceraldehyde', 'reason': 'Reverse-reaction condition mapping remains unresolved.'},
                                {'item': 'Table 1 WT L-threonate NAD+', 'reason': 'Duplicate enzyme-substrate pair; first table condition retained.'},
                                {'item': 'Table 4 WT', 'reason': 'Repeated Table 1 measurement.'},
                                {'item': 'Table 4 ND variants', 'reason': 'No absolute kcat; ND is not zero.'},
                                {'item': 'Table 3 and other enzymes', 'reason': 'Varied-cofactor mapping or full construct not yet audited.'}],
                      quality_caveats=['Kinetic temperature not reported; 35 C is HPLC temperature, not used.',
                                       'Some Km exceed substrate range. Published values retained with caveat; no post-hoc threshold added.',
                                       'Development row-level citation overlap is unavailable.'],
                      raw_sha256={p.name: sha(p) for p in RAW.iterdir() if p.is_file()},
                      exact_sequence_substrate_overlap=0, model_predictions_run=False)
    (SOURCE / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    families = {line.split('\t')[0] for line in clusters.read_text().splitlines() if line}
    audit = dict(audited_on='2026-10-04', source_id=SOURCE.name, article_doi=DOI,
                 method='MMseqs2 easy-search and linclust; full constructs, tag-free cores and native reference checked',
                 mmseqs_version=version, min_identity=.3, coverage=.8, coverage_mode=0,
                 development_target_sha256=sha(REF), candidate_records=dict(count=len(rows), sha256=sha(csvfile)),
                 construct_sequences_sha256=sha(fasta), core_sequences_sha256=sha(cores),
                 homology_hits_sha256=sha(hits), core_homology_hits_sha256=sha(corehits), family_cluster_sha256=sha(clusters),
                 candidate_mmseqs_families=len(families), accepted_records=len(rows), accepted_unique_sequences=5,
                 accepted_unique_substrates=6, exact_sequence_substrate_overlap=0, status='accepted_homology_cold_pool',
                 provenance_sha256=sha(SOURCE / 'provenance.json'), readiness_gate_passes=False,
                 claim_boundary='Accepted curation pool only; global cap audit required; no model scores.')
    assert before == {str(p): sha(p) for p in protected}
    (SOURCE / 'homology-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({'accepted_source_records': len(rows), 'local_families': len(families), 'protected_hashes': before}))


if __name__ == '__main__':
    main()
