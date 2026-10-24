"""Summarize source-level eligibility and audit G6PD reference variants.

Classifies Fordyce portal studies by publication year and direct-kcat
endpoints (portal metadata years are screening clues, not verified dates),
builds six in-silico G6PD reference variants, screens them with a pinned
MMseqs2 build via WSL, and records the portal's missing license and the
Sabath ancestral-sequence state. Writes source-review-20261004/{
g6pd-reference-variants.fasta, g6pd-variant-hits.tsv, variant-mmseqs.log,
review-audit.json}; protected hashes are re-verified before and after and
accepted_additions stays 0.

中文：本脚本汇总来源层面的资格并审计 G6PD 参考变异：按发表年份与直接 kcat 端点对
Fordyce 门户研究分类（门户元数据年份只是筛查线索，不是已核实日期），构建六条 in-silico
G6PD 参考变异并经 WSL 用钉住版本的 MMseqs2 筛查，同时记录门户缺失许可与 Sabath 祖先
序列状态。写出 source-review-20261004/{g6pd-reference-variants.fasta,
g6pd-variant-hits.tsv, variant-mmseqs.log, review-audit.json}；保护哈希在前后各复验一次，
accepted_additions 保持 0。治理角色：fail-closed 来源评审审计收据。
"""
import collections
import hashlib
import json
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/external/temporal-absolute-kinetics/source-review-20261004'
REF=ROOT/'artifacts/homology-final/unique_proteins.fasta'
def sha(p):
    """Return the SHA256 hex digest of a file.

    中文：返回文件的 SHA256 十六进制摘要。
    """
    return hashlib.sha256(p.read_bytes()).hexdigest()
def linux(p):
    """Map a Windows path to its WSL /mnt/d/ equivalent for the MMseqs invocation.

    中文：把 Windows 路径映射为 WSL 的 /mnt/d/ 等价路径，供 MMseqs 调用使用。
    """
    return '/mnt/d/'+p.relative_to(Path('D:/')).as_posix()

def main():
    """Screen portal studies and G6PD variants, then write the review audit.

    中文：筛查门户研究与 G6PD 变异，随后写出评审审计。
    """
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in protected}
    assert sha(REF)=='69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    data=json.loads((OUT/'fordyce-data.json').read_text())
    studies=[]
    for pub in sorted({r['p'] for r in data}):
        rows=[r for r in data if r['p']==pub]
        years=sorted({r['y'] for r in rows})
        direct=sum(r['mt']=='kcat' for r in rows)
        after=all(int(y)>=2023 for y in years)
        studies.append(dict(study_key=pub,years=years,portal_rows=len(rows),
            endpoints=dict(collections.Counter(r['mt'] for r in rows)),
            protein_systems=sorted({r['pw'] for r in rows}),direct_kcat_rows=direct,
            disposition=('time_gate_fail_portal_year' if not after else ('endpoint_not_direct_kcat' if not direct else 'pending_primary_source_and_dataset_license'))))
    native=''.join((OUT/'g6pd-reference.fasta').read_text().splitlines()[1:])
    assert len(native)==515
    variants=['WT','P409R','P409Q','P409S','G410A','G410C','G410D']
    sequences={}
    for v in variants:
        s=list(native)
        if v!='WT':
            i=int(v[1:-1])-1;assert s[i]==v[0];s[i]=v[-1]
        sequences[v]=''.join(s)
    query=OUT/'g6pd-reference-variants.fasta'
    query.write_text(''.join(f'>{k}\n{s}\n' for k,s in sequences.items()))
    exe=['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version=subprocess.check_output(exe+['version']).decode().strip()
    assert version=='5d152c612b6ad2a56f657b7a02c127eceaea2a75'
    hits=OUT/'g6pd-variant-hits.tsv'
    cmd=exe+['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'variant-tmp'),
             '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
             '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'variant-mmseqs.log').open('wb') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    results={v:[] for v in variants}
    for line in hits.read_text().splitlines():
        r=line.split('\t')
        if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8:results[r[0]].append(r)
    assert all(results.values())
    audit=dict(date='2026-10-04',portal_url='https://www.fordycelab.com/data',
        data_url='https://fordycelab.github.io/data-portal/data.json',
        portal_rows=len(data),study_count=len(studies),studies=studies,
        post_2022_direct_kcat_rows=sum(s['direct_kcat_rows'] for s in studies if all(int(y)>=2023 for y in s['years'])),
        portal_license='No license identified in repository metadata, root listing or README; absence is not permission. Separate original-study data licenses remain unresolved.',
        repository_commit=json.loads((OUT/'portal-commits.json').read_text())[0]['sha'],
        note='Portal metadata years and endpoint labels are screening clues, not verified publication dates, independent sample counts or dataset license grants.',
        g6pd=dict(doi='10.3390/ijms26178464',reference='UniProt P11413',
            status='excluded_reference_family_homology',query_scope='Reference human G6PD plus article-reported six variants; not proof of exact assayed plasmid sequence.',
            mmseqs_version=version,command=cmd,reference_sha256=sha(REF),
            results={v:dict(qualifying_hits=len(h),best_hit=max(h,key=lambda r:float(r[7]))) for v,h in results.items()}),
        sabath=dict(doi='10.1093/gbe/evae016',status='pending_ancestral_sequence_verification',
            supplement_url='https://www.ebi.ac.uk/biostudies/files/S-EPMC10853004/evae016_supplementary_data.pdf',
            evidence='Eight-page supplement acquired. Page 2 visually contains ancestral sequence alignment; text extraction omits alignment letters. Authoritative sequence is not absent; machine-readable/double-checked transcription remains required.',
            figshare='https://figshare.com/articles/dataset/Figure_2_alignment_for_ASR_fasta/25050875',
            figshare_access='API returned HTTP 403; no access-control bypass attempted',
            homology_run=False),
        accepted_additions=0,predictions_generated=False,protected_hashes=before,
        artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='review-audit.json'})
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    (OUT/'review-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(portal_rows=len(data),studies=len(studies),recent_direct_kcat=audit['post_2022_direct_kcat_rows'],g6pd_hit_variants=len(results),accepted=0)))

if __name__=='__main__':main()
