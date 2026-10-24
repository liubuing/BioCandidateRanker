"""One diagnostic test pass per frozen retrained checkpoint.

Reads the evaluation freeze and retrain protocol manifests under
``artifacts/current-code-retrain`` plus ``configs/current_code_retrain_protocol.json``,
rebuilds the 1,646-row internal test partition from the pinned UniKP source and
homology split, and writes per-seed prediction CSVs with receipts and
``summary.json`` under ``internal-test/``. Governance role: fail-closed,
previously-observed-test diagnostic only -- every SHA-256 pin is re-verified,
existing receipts must match the freeze and prediction hashes or the run
aborts, and no checkpoint or selection is changed.

中文：对每个冻结重训检查点做一次诊断性测试评估：校验评估冻结与重训协议
清单后，重建 1,646 行内部测试分区，每种子恰评估一次并写收据。
"""
from pathlib import Path
import csv,json,hashlib,sys,os
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/current-code-retrain'
OUT=BASE/'internal-test'
sys.path.insert(0,str(BASE/'source'))
os.environ['OMP_NUM_THREADS']='4'
import torch
from biocandidate.training import load_checkpoint
from biocandidate.cli import _collect_task_outputs
from biocandidate.data.unikp import read_unikp_json,aggregate_pair_measurements
from biocandidate.data.homology import apply_split_manifest
from biocandidate.evaluation import regression_metrics

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    """Verify freezes, run each seed checkpoint once over the internal test, write receipts.
中文：校验冻结后，每个种子检查点在内部测试上跑一次并写收据。
"""
    freeze=read(BASE/'evaluation-freeze.json')
    for e in freeze['files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
    protocol=read(ROOT/'configs/current_code_retrain_protocol.json')
    for e in protocol['frozen_files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
    os.environ['TORCH_HOME']=protocol['torch_home']
    assert sha(protocol['source_path'])==protocol['source_sha256']
    manifest=read(ROOT/'configs/unikp_source_manifest.json')
    records=read_unikp_json(protocol['source_path']).records
    records=apply_split_manifest(records,ROOT/'artifacts/homology-final/homology_split.json',source_identity=manifest['identity'])
    records=[r for r in aggregate_pair_measurements(records) if r.split=='test']
    with (ROOT/'artifacts/external/sota-homology-cold/data/test.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    assert len(records)==len(rows)==1646
    for r,row in zip(records,rows):
        assert r.sequence==row['sequence'] and r.substrate_smiles==row['smiles'] and abs(r.log10_kcat-float(row['log10_kcat']))<1e-10
    OUT.mkdir(exist_ok=True);results=[]
    for seed in [7,42,123]:
        dest=OUT/f'seed{seed}-predictions.csv';receipt=OUT/f'seed{seed}-receipt.json'
        if receipt.exists():
            old=read(receipt);assert sha(dest)==old['predictions_sha256'] and old['freeze_sha256']==sha(BASE/'evaluation-freeze.json');results.append(old);continue
        if dest.exists():raise FileExistsError(dest)
        with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):model,payload=load_checkpoint(BASE/f'seed{seed}/best.pt',torch.device('cuda'))
        p,sd,y=_collect_task_outputs(model,records,'log10_kcat',4,torch.device('cuda'))
        assert len(p)==1646 and torch.isfinite(p).all() and torch.isfinite(sd).all()
        with dest.open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['pair_id','prediction','label']);w.writerows((row['pair_id'],float(a),float(b)) for row,a,b in zip(rows,p,y))
        r={'seed':seed,'selected_epoch':payload['epoch'],'metrics':regression_metrics(p,y),'predictions_sha256':sha(dest),'freeze_sha256':sha(BASE/'evaluation-freeze.json')}
        receipt.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');results.append(r);print(json.dumps(r),flush=True)
        del model;torch.cuda.empty_cache()
    (OUT/'summary.json').write_text(json.dumps({'scope':'Current-code retrained model; previously observed internal test, diagnostic only','runs':results},indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
