"""Figure-transcribed ancestor verification and exclusion-only homology prescreen."""
import hashlib
import json
import re
import subprocess
from pathlib import Path
from Bio.Seq import Seq

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
SOURCE=BASE/'source-review-20261004'
OUT=BASE/'sabath-verification-20261004'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
# Supplementary Fig S2, upper sequence, six 60-column lines plus 11 residues.
LINES=[
 'MEVQQVLHMNGGDGETSYANNSSLQKKVISMAKPIIEEAITDLYCSNTFPESLTIADLGC',
 'SSGPNTLLVVSEIIDTIHKKCRQLGRPPPEFQVFLNDLPGNDFNTIFRSLPSFYEKLKKE',
 'KGNGFGPCFIAGVPGSFYGRLFPSKSLHFVHSSYSLHWLSQVPQGLENEAGGVLNKGNIY',
 'IAKTSPPSVLKAYLEQFQRDFSLFLKSRSEEMVPGGRMVLTLLGRRSADPSSKECCYLWE',
 'LLAQALNDMVSEGLIEEEKLDSFNLPYYAPSAEEVKAVIEKEGSFTIDRLETFEVNWDAN',
 'NNDDDGNKSDKLKRGKNVAKYIRAVAEPMLASHFGEAIMDELFRRYAKIVAEHMSKEKTK',
 'YINLVISLTKK']
SITES={54:'T',88:'P',123:'N',157:'H',181:'I',238:'L',267:'Y',272:'A',301:'N',322:'I',329:'M',348:'K',361:'Y'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def linux(p):return '/mnt/d/'+p.relative_to(Path('D:/')).as_posix()

def main():
    OUT.mkdir(exist_ok=True)
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in protected}
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    assert [len(x) for x in LINES]==[60]*6+[11]
    native=''.join(LINES);assert len(native)==371
    assert all(native[i-1]==aa for i,aa in SITES.items())
    seqs={'AncSBXMT_figure':native,'AncSBXMT_Y267Q_figure':native[:266]+'Q'+native[267:]}
    # Independent local evidence from published mutagenesis primers, not designed primers.
    primer_text=(SOURCE/'sabath-page-8.txt').read_text(encoding='utf-8')
    checks=[]
    for name,dna in re.findall(r'([A-Z]\d+[A-Z][^\s]*)\s+([acgt]{15,})',primer_text):
        mutation=re.match(r'([A-Z])(\d+)([A-Z])',name)
        old,position,new=mutation.groups();position=int(position)
        assert native[position-1]==old
        mutant=native[:position-1]+new+native[position:]
        matches=[]
        for strand in [dna,str(Seq(dna).reverse_complement())]:
            for frame in range(3):
                usable=strand[frame:];usable=usable[:len(usable)//3*3]
                peptide=str(Seq(usable).translate())
                start=mutant.find(peptide)
                if len(peptide)>=5 and start>=0 and start<=position-1<start+len(peptide):
                    matches.append(dict(peptide=peptide,start_1based=start+1))
        checks.append(dict(primer_name=name,mutation=f'{old}{position}{new}',matches=matches))
    query=OUT/'ancestor-queries.fasta'
    query.write_text(''.join(f'>{k}\n{s}\n' for k,s in seqs.items()))
    exe=['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version=subprocess.check_output(exe+['version']).decode().strip()
    assert version=='5d152c612b6ad2a56f657b7a02c127eceaea2a75'
    hits=OUT/'development-hits.tsv'
    cmd=exe+['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'search-tmp'),
        '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
        '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'mmseqs.log').open('wb') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    result={k:[] for k in seqs}
    for line in hits.read_text().splitlines():
        r=line.split('\t')
        if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8:result[r[0]].append(r)
    audit=dict(date='2026-10-04',doi='10.1093/gbe/evae016',
        source_pdf_sha256=sha(SOURCE/'sabath-supplement.pdf'),figure_page=2,figure='S2A upper row',
        figure_render_sha256=sha(SOURCE/'sabath-page2-high.png'),transcribed_lines=LINES,
        residue_count=len(native),landmark_checks=SITES,primer_checks=checks,
        primer_matches=sum(bool(c['matches']) for c in checks),primer_count=len(checks),
        name_note='Figure upper-row label AncBSXMT and body AncSBXMT refer to the same ancestor here; no substitution with descendant AncSAMT.',
        limitation='Manual figure transcription reviewed against image and local primer evidence, not independent full-length machine-readable deposition. Query is catalytic ancestor, not complete pET15b fusion. Primer matches validate local segments only.',
        mmseqs_version=version,command=cmd,development_sha256=sha(REF),
        results={k:dict(qualifying_hits=len(v),best_hit=max(v,key=lambda r:float(r[7])) if v else None) for k,v in result.items()},
        source_disposition='exclude_family_cold_source' if all(result.values()) else 'pending_exact_construct_and_remaining_gates',
        accepted_additions=0,predictions_generated=False,protected_hashes=before,
        artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='verification-audit.json'})
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    (OUT/'verification-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(residues=len(native),primer_matches=audit['primer_matches'],primer_count=len(checks),results=audit['results'],disposition=audit['source_disposition']),indent=2))

if __name__=='__main__':main()
