"""Frozen outer-family evaluation of two validation checkpoint selectors."""
from pathlib import Path
import csv,json,hashlib,os,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'artifacts/nested-family-selection';OUT=BASE/'outer-evaluation'
sys.path.insert(0,str(BASE/'source'));os.environ['OMP_NUM_THREADS']='4'
import torch
import numpy as np
from biocandidate.training import load_checkpoint
from biocandidate.cli import _collect_task_outputs
from biocandidate.data.unikp import read_unikp_json,aggregate_pair_measurements
from biocandidate.data.homology import apply_split_manifest

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write(p,d):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');os.replace(tmp,p)

def draws(p,y,s,groups,w):
    e=p-y
    row=np.column_stack([np.ones(len(y)),e**2,abs(e),p,y,p*p,y*y,p*y,.5*(e/s)**2+np.log(s)+.5*np.log(2*np.pi),abs(e)<=s,abs(e)<=2*s])
    stats=np.zeros((w.shape[1],row.shape[1]));np.add.at(stats,groups,row)
    n,sq,ae,x,z,xx,zz,xz,nll,c1,c2=(w@stats).T
    return np.column_stack([np.sqrt(sq/n),ae/n,(xz-x*z/n)/np.sqrt(np.maximum(0,xx-x*x/n)*np.maximum(0,zz-z*z/n)),nll/n,c1/n,c2/n])

def main():
    freeze=BASE/'outer-evaluation-freeze.json';f=read(freeze)
    for e in f['files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
    protocol=read(ROOT/'configs/nested_family_selection_protocol.json')
    for e in protocol['frozen_files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
    assert sha(protocol['source_path'])==protocol['source_sha256']
    os.environ['TORCH_HOME']=protocol['torch_home']
    split=read(BASE/'split.json');outer_rows={r['source_row'] for r in split['rows'] if r['split']=='test'}
    records=read_unikp_json(protocol['source_path']).records
    records=aggregate_pair_measurements(apply_split_manifest(records,BASE/'split.json',source_identity=read(BASE/'source-manifest.json')['identity']))
    records=[r for r in records if r.split=='test'];assert len(records)==2127
    mapping={r['source_row']:r['cluster'] for r in split['rows']};families=[mapping[r.source_row] for r in records]
    unique,groups=np.unique(families,return_inverse=True);assert len(unique)==266
    for r in records:assert set(r.source_rows)<=outer_rows
    w=np.vstack([np.ones(266),np.random.default_rng(20260928).multinomial(266,np.full(266,1/266),size=5000)])
    OUT.mkdir(exist_ok=True);all_stats={'nll':[],'rmse':[]};receipts=[]
    for seed in [7,42,123]:
        for rule,file in [('nll','best.pt'),('rmse','best_rmse.pt')]:
            dest=OUT/f'seed{seed}-{rule}.csv';receipt=OUT/f'seed{seed}-{rule}-receipt.json'
            if receipt.exists():
                r=read(receipt);assert r['freeze_sha256']==sha(freeze) and sha(dest)==r['predictions_sha256']
                with dest.open(newline='') as stream:rows=list(csv.DictReader(stream))
                assert [int(x['source_row']) for x in rows]==[x.source_row for x in records]
                p,y,s=[np.array([float(x[k]) for x in rows]) for k in ['prediction','label','sd']]
            else:
                if dest.exists():raise FileExistsError(dest)
                with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):model,payload=load_checkpoint(BASE/f'seed{seed}'/file,torch.device('cuda'))
                pred,sd,label=_collect_task_outputs(model,records,'log10_kcat',4,torch.device('cuda'))
                p,y,s=[v.numpy().astype(float) for v in [pred,label,sd]]
                assert len(p)==2127 and np.isfinite(p).all() and np.isfinite(s).all() and (s>0).all()
                with dest.open('w',encoding='utf-8',newline='') as stream:
                    writer=csv.writer(stream);writer.writerow(['source_row','family','prediction','label','sd']);writer.writerows(zip([r.source_row for r in records],families,p,y,s))
                r={'seed':seed,'rule':rule,'epoch':payload['epoch'],'checkpoint_sha256':sha(BASE/f'seed{seed}'/file),'predictions_sha256':sha(dest),'freeze_sha256':sha(freeze)}
                write(receipt,r);del model;torch.cuda.empty_cache()
            assert np.allclose(y,[r.log10_kcat for r in records],rtol=0,atol=1e-6)
            stat=draws(p,y,s,groups,w);assert np.isfinite(stat).all();all_stats[rule].append(stat);receipts.append(r)
            print(seed,rule,stat[0].tolist(),flush=True)
    names=['rmse','mae','pearson','gaussian_nll','coverage_1sigma','coverage_2sigma']
    mean={k:np.mean(v,axis=0) for k,v in all_stats.items()};delta=mean['rmse']-mean['nll']
    result={'scope':'Single nested family holdout within historical training data; not pristine external confirmation','rows':2127,'families':266,'bootstrap_replicates':5000,'bootstrap_seed':20260928,'rules':{k:{name:{'mean':float(mean[k][0,i]),'seed_sd':float(np.std([x[0,i] for x in all_stats[k]],ddof=1))} for i,name in enumerate(names)} for k in all_stats},'rmse_rule_minus_nll_rule':{name:{'difference':float(delta[0,i]),'ci95':np.quantile(delta[1:,i],[.025,.975]).tolist()} for i,name in enumerate(names)},'receipts':receipts,'limitations':'Fixed seeds and single split; exploratory secondary metrics without multiplicity adjustment; zero-crossing intervals do not prove equivalence'}
    if (OUT/'summary.json').exists():raise FileExistsError(OUT/'summary.json')
    write(OUT/'summary.json',result)

if __name__=='__main__':main()
