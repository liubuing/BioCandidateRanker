"""Frozen internal evaluation of the completed adapted CatPred ensembles (WSL)."""
from pathlib import Path
import csv, gzip, hashlib, json, os, sys, time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/catpred-mode-b'
OUT = BASE / 'internal-test'
PROTOCOL = ROOT / 'configs/catpred_internal_test_protocol.json'
os.environ.setdefault('CATPRED_CACHE_PATH', str(BASE/'cache'))
os.environ.setdefault('CATPRED_ESM_BATCH_SIZE', '2')
os.environ.setdefault('TORCH_HOME', '/mnt/d/DisorderFlowRuntime/cache/torch')
os.environ.setdefault('OMP_NUM_THREADS', '4')
sys.path[:0] = [str(BASE/'upstream'), str(ROOT/'src')]

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def write(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    os.replace(tmp, path)

def main():
    import numpy as np
    import torch
    protocol = json.loads(PROTOCOL.read_text())
    for e in protocol['frozen_files']:
        assert sha(ROOT/e['path']) == e['sha256'], e['path']
    manifest = json.loads((BASE/'full-recoverable/completed-ensemble-manifest.json').read_text())
    members = manifest['members']
    assert len(members)==30 and {(m['group_seed'],m['member']) for m in members}=={(s,i) for s in [7,42,123] for i in range(10)}
    for m in members: assert sha(ROOT/m['checkpoint']) == m['checkpoint_sha256']
    OUT.mkdir(exist_ok=True)
    with (ROOT/protocol['test_source']).open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    assert len(rows)==1646 and len({r['pair_id'] for r in rows})==1646
    # Prediction input deliberately excludes test labels.
    records = {}
    with (OUT/'inputs.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f); w.writerow(['pair_id','reactant_smiles','pdbpath','sequence'])
        for r in rows:
            key='seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]+'.pdb'
            records[key]={'name':key,'seq':r['sequence']}
            w.writerow([r['pair_id'],r['smiles'],key,r['sequence']])
    with gzip.open(OUT/'protein_records.json.gz','wt',encoding='utf-8') as f: json.dump(records,f)
    from catpred.data.esm_utils import get_many_esm_reprs, ESM_MAX_LENGTH, ESM_EMBED_DIM
    from catpred.data.cache_utils import cache_entry_path
    sequences = list(dict.fromkeys(r['seq'] for r in records.values()))
    features=[]
    for offset in range(0,len(sequences),32):
        batch=sequences[offset:offset+32]
        values=get_many_esm_reprs(batch,device='cpu',batch_size=2)
        for seq in batch:
            value=values[seq]
            assert tuple(value.shape)==(min(len(seq),ESM_MAX_LENGTH-1),ESM_EMBED_DIM)
            assert torch.isfinite(value).all()
            path=cache_entry_path(path='esm/proteins',cache_key=seq)
            features.append({'sequence_sha256':hashlib.sha256(seq.encode('ascii')).hexdigest(),'cache_path':str(path.relative_to(ROOT)),'sha256':sha(path)})
        del values
        print('Test features',len(features),'/',len(sequences),flush=True)
    write(OUT/'feature-manifest.json',{'features':features,'scope':'Frozen pretrained features only; no fitting on test data'})
    from catpred_memory_entry import install_adapters
    install_adapters()
    from catpred.args import PredictArgs
    from catpred.train.make_predictions import make_predictions
    from biocandidate.evaluation import regression_metrics
    reports=[]
    for seed in [7,42,123]:
        target=OUT/f'seed{seed}-predictions.csv'; receipt=OUT/f'seed{seed}-receipt.json'
        if receipt.exists():
            old=json.loads(receipt.read_text());assert old['protocol_sha256']==sha(PROTOCOL) and sha(target)==old['predictions_sha256']
            reports.append(old);continue
        raw=OUT/f'seed{seed}-upstream.csv'
        if raw.exists() or target.exists(): raise FileExistsError('Partial prediction output requires audit before retry')
        paths=[str(ROOT/m['checkpoint']) for m in members if m['group_seed']==seed]
        args=PredictArgs().parse_args(['--test_path',str(OUT/'inputs.csv'),'--preds_path',str(raw),'--protein_records_path',str(OUT/'protein_records.json.gz'),'--smiles_columns','reactant_smiles','--batch_size','16','--num_workers','0','--checkpoint_paths',*paths])
        print('Predicting ensemble',seed,flush=True)
        preds=np.asarray(make_predictions(args),dtype=float)
        assert preds.shape==(1646,1) and np.isfinite(preds).all()
        with target.open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['pair_id','label_log10_kcat','prediction_log10_kcat'])
            for r,v in zip(rows,preds[:,0]):w.writerow([r['pair_id'],r['log10_kcat'],float(v)])
        labels=torch.tensor([float(r['log10_kcat']) for r in rows],dtype=torch.float64)
        result={'seed':seed,'protocol_sha256':sha(PROTOCOL),'predictions_sha256':sha(target),'metrics':regression_metrics(torch.tensor(preds[:,0],dtype=torch.float64),labels)}
        write(receipt,result);reports.append(result)
        print(json.dumps(result),flush=True)
    write(OUT/'summary.json',{'scope':'Adapted CatPred Mode B; internal homology-cold test only','seeds':reports,'mean_sd':{k:{'mean':float(np.mean([r['metrics'][k] for r in reports])),'sd':float(np.std([r['metrics'][k] for r in reports],ddof=1))} for k in ['rmse','mae','pearson']}})

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    import fcntl
    lock=(OUT/'evaluation.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    write(OUT/'execution.json',{'status':'running','started_unix':time.time()})
    try:
        main()
    except Exception as exc:
        write(OUT/'execution.json',{'status':'failed','error':repr(exc)});raise
    write(OUT/'execution.json',{'status':'completed','finished_unix':time.time()})
