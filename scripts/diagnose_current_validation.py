"""Read-only model diagnostics on validation data; no checkpoint selection changes.

Reads the frozen retrain protocol, the completed-training manifest, the three
per-seed checkpoints under ``artifacts/current-code-retrain``, and the
1,640-row homology-cold validation CSV; verifies every SHA-256 pin, runs each
frozen checkpoint over the validation partition, and writes per-seed prediction
CSVs plus ``report.json`` under ``validation-diagnostic/`` with stratified
metrics (length and canonicality strata), SD-coverage, and seed-disagreement
statistics. Governance role: descriptive diagnosis only -- checkpoints,
selected epochs, and test predictions are never modified or re-selected, and a
pre-existing ``report.json`` aborts the run (fail-closed).

中文：对验证数据做只读模型诊断；不改变任何检查点选择。读取冻结重训协议、
完成训练清单与三个种子检查点，在 1,640 行同源冷验证集上分层诊断。
"""
from pathlib import Path
import csv,json,hashlib,os,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'artifacts/current-code-retrain'
OUT=BASE/'validation-diagnostic';sys.path.insert(0,str(BASE/'source'))
os.environ['OMP_NUM_THREADS']='4'
import numpy as np
import torch
from biocandidate.training import load_checkpoint
from biocandidate.cli import _collect_task_outputs
from biocandidate.data.unikp import read_unikp_json,aggregate_pair_measurements
from biocandidate.data.homology import apply_split_manifest

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def metrics(p,y,sd):
    """Compute RMSE/MAE/bias, SD coverage, and constant-free Gaussian NLL for one stratum.
中文：计算单层的 RMSE/MAE/偏差、SD 覆盖与常数无关高斯 NLL。
"""
    e=p-y
    return {'n':len(y),'rmse':float(np.sqrt(np.mean(e**2))),'mae':float(np.mean(abs(e))),
        'bias':float(e.mean()),'mean_sd':float(sd.mean()),'coverage_1sigma':float(np.mean(abs(e)<=sd)),
        'coverage_2sigma':float(np.mean(abs(e)<=2*sd)),
        'gaussian_nll_without_constant':float(np.mean(.5*(e/sd)**2+np.log(sd)))}

def main():
    """Verify all frozen inputs, then diagnose the three seed checkpoints on validation.
中文：校验全部冻结输入后，对三个种子检查点做验证集诊断。
"""
    if (OUT/'report.json').exists():raise FileExistsError(OUT/'report.json')
    OUT.mkdir(exist_ok=True)
    protocol=read(ROOT/'configs/current_code_retrain_protocol.json');os.environ['TORCH_HOME']=protocol['torch_home']
    for e in protocol['frozen_files']:assert sha(ROOT/e['path'])==e['sha256']
    freeze=read(BASE/'completed-training-manifest.json')
    for r in freeze['runs']:assert sha(BASE/f"seed{r['seed']}/best.pt")==r['checkpoint_sha256']
    assert sha(protocol['source_path'])==protocol['source_sha256']
    records=read_unikp_json(protocol['source_path']).records
    records=apply_split_manifest(records,ROOT/'artifacts/homology-final/homology_split.json',source_identity=read(ROOT/'configs/unikp_source_manifest.json')['identity'])
    records=[r for r in aggregate_pair_measurements(records) if r.split=='validation']
    with (ROOT/'artifacts/external/sota-homology-cold/data/validation.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    assert len(rows)==len(records)==1640
    for r,row in zip(records,rows):assert r.sequence==row['sequence'] and r.substrate_smiles==row['smiles'] and abs(r.log10_kcat-float(row['log10_kcat']))<1e-10
    lengths=np.array([len(r.sequence) for r in records]);noncanonical=np.array([bool(set(r.sequence)-set('ACDEFGHIKLMNPQRSTVWY')) for r in records])
    groups={'all':np.ones(len(rows),dtype=bool),'length_le_128':lengths<=128,'length_129_512':(lengths>128)&(lengths<=512),'length_gt_512':lengths>512,'noncanonical':noncanonical,'canonical':~noncanonical}
    result={'scope':'Validation-only descriptive diagnosis; no changes to model, selected checkpoint or test predictions','script_sha256':sha(__file__),'validation_sha256':sha(ROOT/'artifacts/external/sota-homology-cold/data/validation.csv'),'seeds':{}}
    all_preds=[]
    for seed in [7,42,123]:
        with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):model,payload=load_checkpoint(BASE/f'seed{seed}/best.pt',torch.device('cuda'))
        pred,sd,label=_collect_task_outputs(model,records,'log10_kcat',4,torch.device('cuda'));p=pred.numpy().astype(float);s=sd.numpy().astype(float);y=label.numpy().astype(float)
        assert np.isfinite(p).all() and np.isfinite(s).all() and (s>0).all()
        dest=OUT/f'seed{seed}.csv'
        with dest.open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['pair_id','prediction','label','sd','sequence_length']);w.writerows((r['pair_id'],a,b,c,int(length)) for r,a,b,c,length in zip(rows,p,y,s,lengths))
        logs=[json.loads(line) for line in (BASE/f'seed{seed}/train.log').read_text().splitlines() if line.startswith('{')];epochs=[r for r in logs if 'epoch' in r]
        best_rmse=min(epochs,key=lambda r:r['validation_tasks']['log10_kcat']['rmse'])
        selected=next(r for r in epochs if r['epoch']==payload['epoch'])
        result['seeds'][str(seed)]={'selected_epoch':payload['epoch'],'lowest_logged_validation_rmse_epoch':best_rmse['epoch'],'selected_logged_rmse':selected['validation_tasks']['log10_kcat']['rmse'],'lowest_logged_validation_rmse':best_rmse['validation_tasks']['log10_kcat']['rmse'],'epochs':epochs,'strata':{k:metrics(p[m],y[m],s[m]) for k,m in groups.items() if m.any()},'prediction_sha256':sha(dest)}
        all_preds.append(p);print(seed,json.dumps(result['seeds'][str(seed)]['strata']['all']),flush=True)
        del model;torch.cuda.empty_cache()
    a=np.array(all_preds);result['seed_disagreement']={'mean_prediction_sd':float(a.std(axis=0,ddof=1).mean()),'p95_prediction_sd':float(np.quantile(a.std(axis=0,ddof=1),.95)),'prediction_correlations':np.corrcoef(a).tolist()}
    (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
