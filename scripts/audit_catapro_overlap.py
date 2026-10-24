"""Download a pinned public dataset for input overlap audit only (no predictions).

Fetches the CataPro kcat-data_0.4simi-10fold.csv blob named by
artifacts/pcb-readiness/CataPro/tree.json, verifies its Git blob SHA1, and
compares exact sequence / RDKit canonical SMILES pairs against every local
partition recorded in artifacts/external/sota-homology-cold/data-manifest.json
(each re-hashed before use). Writes artifacts/pcb-readiness/catapro-data-audit.json;
the receipt records predictions_generated=false and its interpretation is
deliberately scoped to input overlap, with no homology or license claim.

中文：本脚本仅为输入重叠审计下载一个钉住的公开数据集（不做预测）：按
artifacts/pcb-readiness/CataPro/tree.json 指定的 blob 获取 CataPro 的
kcat-data_0.4simi-10fold.csv，校验其 Git blob SHA1，并将精确的 序列/RDKit 规范 SMILES
对与 artifacts/external/sota-homology-cold/data-manifest.json 记录的每个本地分区（使用前
逐一重新哈希）比对。写出 artifacts/pcb-readiness/catapro-data-audit.json；收据记录
predictions_generated=false，解释范围刻意限定为输入重叠，不做同源性或许可声明。
治理角色：只读审计收据（重叠证据，不产生任何预测或标签）。
"""
from pathlib import Path
from collections import Counter
import base64,csv,hashlib,io,json,urllib.request
from rdkit import Chem, RDLogger

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'artifacts/pcb-readiness'
receipt=json.loads((OUT/'source-receipt.json').read_text())
repo=next(r for r in receipt['repositories'] if r['repository']=='zchwang/CataPro')
tree=json.loads((OUT/'CataPro/tree.json').read_text())
entry=next(x for x in tree['tree'] if x['path']=='datasets/kcat-data_0.4simi-10fold.csv')
target=OUT/'CataPro/kcat-data_0.4simi-10fold.csv'
if not target.exists():
    req=urllib.request.Request(entry['url'],headers={'User-Agent':'BioCandidateRanker-overlap-audit'})
    with urllib.request.urlopen(req,timeout=30) as r:blob=json.load(r)
    data=base64.b64decode(blob['content'])
    assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==entry['sha']
    target.write_bytes(data)
data=target.read_bytes()
assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==entry['sha']
reader=csv.DictReader(io.StringIO(data.decode('utf-8-sig')));rows=list(reader)
print('columns:',reader.fieldnames,'rows:',len(rows),flush=True)
result={'upstream_commit':repo['commit'],'path':entry['path'],'bytes':len(data),
        'sha256':hashlib.sha256(data).hexdigest(),'rows':len(rows),'columns':reader.fieldnames,
        'column_nonempty_counts':{k:sum(r.get(k) not in ('',None) for r in rows) for k in reader.fieldnames},
        'predictions_generated':False}
RDLogger.DisableLog('rdApp.*')
cache={}
for value in {r['Smiles'] for r in rows}:
    mol=Chem.MolFromSmiles(value)
    cache[value]=Chem.MolToSmiles(mol) if mol is not None else None
seqs={r['Sequence'].strip() for r in rows}
pairs={(r['Sequence'].strip(),cache[r['Smiles']]) for r in rows if cache[r['Smiles']] is not None}
result['upstream_invalid_smiles_rows']=sum(cache[r['Smiles']] is None for r in rows)
result['overlap_definition']='Exact sequence and RDKit canonical SMILES; row counts in each local partition, across union of all upstream folds. No homology or publication audit.'
result['fold_counts']=dict(Counter(r['fold'] for r in rows))
result['local_partition_overlap']={}
manifest=json.loads((ROOT/'artifacts/external/sota-homology-cold/data-manifest.json').read_text())
for name,part in manifest['partitions'].items():
    path=ROOT/part['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==part['sha256']
    with path.open(encoding='utf-8',newline='') as f: local=list(csv.DictReader(f))
    result['local_partition_overlap'][name]={
      'rows':len(local),'sequence_match_rows':sum(r['sequence'] in seqs for r in local),
      'sequence_smiles_match_rows':sum((r['sequence'],Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles']))) in pairs for r in local)}
result['interpretation']='Published ten-fold ensemble is not an independent baseline on the local internal test when inputs overlap its training-fold union. Use same-split retraining; unit labels do not resolve row-level primary measurement provenance.'
(OUT/'catapro-data-audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
