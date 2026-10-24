"""Reference/core-only prescreen. These outputs must not enter accepted data.

Reads the staged EstSit01 evidence under
``artifacts/external/temporal-absolute-kinetics/esterase-prescreen-20261003``
(two article XMLs, ANG84115.1/WP_053548180.1 FASTAs), builds the four
reference-core variant sequences (WT, E301P, E301P/G215A, +G293A), screens
them against the SHA-256-pinned development reference with frozen MMseqs2, and
parses the primary kinetic table into ``candidate-measurements.csv`` plus
``prescreen-receipt.json``.

Governance role: fail-closed prescreen -- every measurement is staged with
status ``pending_construct_and_substrate_audit`` and ``accepted_additions``
stays 0 (reference cores are not the assayed fusion proteins; cross-paper WT
duplication is unresolved); protected config/audit hashes are re-verified at
the end, and no predictions are generated.

中文：仅参考核心的预筛（输出不得进入接受数据）：对暂存的 EstSit01 证据
构建参考核心、暂挂测量并写收据。
"""
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'artifacts/external/temporal-absolute-kinetics/esterase-prescreen-20261003'
REF = ROOT/'artifacts/homology-final/unique_proteins.fasta'
VERSION = '5d152c612b6ad2a56f657b7a02c127eceaea2a75'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def seq(p):
    return ''.join(p.read_text().splitlines()[1:])

def linux(p):
    """Map a Windows path to its /mnt/d WSL equivalent.
中文：把 Windows 路径映射为其 /mnt/d WSL 等价路径。
"""
    return '/mnt/d/'+p.relative_to(Path('D:/')).as_posix()

def main():
    """Screen reference cores, stage pending measurements, and write the receipt.
中文：筛查参考核心、暂挂测量并写收据。
"""
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in protected}
    native=seq(OUT/'ANG84115.1.fasta')
    assert len(native)==370 and native==seq(OUT/'WP_053548180.1.fasta')
    variants=['WT','E301P','E301P/G215A','E301P/G215A/G293A']
    sequences={}
    for name in variants:
        s=list(native)
        for change in ([] if name=='WT' else name.split('/')):
            i=int(change[1:-1])-1
            assert s[i]==change[0]
            s[i]=change[-1]
        sequences[name.replace('/','_')]=''.join(s)
    query=OUT/'reference-cores.fasta'
    query.write_text(''.join(f'>{k}\n{v}\n' for k,v in sequences.items()),encoding='ascii')
    exe=['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version=subprocess.check_output(exe+['version']).decode().strip()
    assert version==VERSION
    hits=OUT/'core-hits.tsv'
    cmd=exe+['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'core-search-tmp'),
             '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
             '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'core-mmseqs.log').open('wb') as f:
        subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    result={name:[] for name in sequences}
    for line in hits.read_text().splitlines():
        r=line.split('\t');result[r[0]].append(r)
    development={''.join(part.splitlines()[1:]) for part in REF.read_text().split('>') if part}
    # Parse the primary table rather than transcribing the variant numbers.
    tree=ET.parse(OUT/'PMC11676285.xml')
    tables=[n for n in tree.iter('table-wrap') if 'Kinetic parameters and specific activities' in ''.join(n.itertext())]
    assert len(tables)==1
    entries=[]
    for row in tables[0].iter('tr'):
        cells=[''.join(n.itertext()).strip() for n in row]
        if cells and cells[0] in ['Wild-type']+variants[1:]:
            entries.append(cells)
    assert len(entries)==4
    candidates=[dict(source='PMC11180644',doi='10.1186/s40643-024-00776-2',variant='WT',
                     kcat_s_1='5.808',kcat_sd_s_1='',km_mM='0.147',
                     source_location='Results / kinetic analysis; Methods / Determination of enzyme properties')]
    for c in entries:
        candidates.append(dict(source='PMC11676285',doi='10.3390/molecules29245839',
                               variant='WT' if c[0]=='Wild-type' else c[0],
                               kcat_s_1=c[2].split('±')[0].strip(), kcat_sd_s_1=c[2].split('±')[1].strip(),
                               km_mM=c[1].split('±')[0].strip(),source_location='Table 3 / '+c[0]))
    for row in candidates:
        row.update(reference_accession='ANG84115.1',substrate='meso-dimethyl ester (biotin synthesis intermediate)',
                   temperature_C=30,pH=8.0,status='pending_construct_and_substrate_audit',
                   caveat='Reference cores only; not complete assayed fusion proteins. Cross-paper WT duplication unresolved.')
    with (OUT/'candidate-measurements.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(candidates[0]));w.writeheader();w.writerows(candidates)
    assert hashes=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    receipt=dict(date='2026-10-03',status='pending_not_accepted',sources=2,candidate_measurements=5,
                 core_variant_substrate_pairs=4,independent_families_not_established=True,
                 accepted_additions=0,predictions_generated=False,
                 accession_reconciliation='ANG84115.1 from original supplement Table S2 equals WP_053548180.1 cited in later paper; 370 residues.',
                 mmseqs_version=version,command=cmd,homology_hits=result,
                 exact_core_sequence_overlap={k:v in development for k,v in sequences.items()},
                 pending=['complete_assayed_vector_fusion_sequence','substrate_structure_identity',
                          'cross_paper_WT_duplicate_adjudication','final_construct_homology_and_exact_pair_audit',
                          'global_family_mapping_and_cap'],
                 other_screened_source=dict(pmc='PMC10816834',doi='10.3390/ijms25021272',
                    status='pending_measurement_consistency',reason='Table 2 units and row alignment need verification against enzyme concentration and Methods; no rows admitted.'),
                 protected_hashes=hashes,
                 artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='prescreen-receipt.json'},
                 provenance_urls={
                     'PMC11180644.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11180644/fullTextXML',
                     'PMC11676285.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11676285/fullTextXML',
                     'PMC10816834.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10816834/fullTextXML',
                     '40643_2024_776_MOESM2_ESM.docx':'https://static-content.springer.com/esm/art%3A10.1186%2Fs40643-024-00776-2/MediaObjects/40643_2024_776_MOESM2_ESM.docx',
                     'ANG84115.1.fasta':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=ANG84115.1&rettype=fasta&retmode=text',
                     'WP_053548180.1.fasta':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=WP_053548180.1&rettype=fasta&retmode=text'})
    (OUT/'prescreen-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(candidate_measurements=5,core_variants=4,homology_hits={k:len(v) for k,v in result.items()},accepted_additions=0)))

if __name__=='__main__':
    main()
