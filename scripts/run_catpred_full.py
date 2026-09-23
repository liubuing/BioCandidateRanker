"""Run a frozen CatPred development ensemble; test data are never supplied."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/catpred-mode-b'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seed',type=int,required=True)
    args=parser.parse_args()
    protocol_path=ROOT/'configs/catpred_mode_b_protocol.json'
    protocol=json.loads(protocol_path.read_text())
    assert args.seed in protocol['seeds']
    assert sha(__file__)==protocol['runner_sha256']
    for entry in protocol['frozen_files']:
        assert sha(ROOT/entry['path'])==entry['sha256'],entry['path']
    manifest=json.loads((BASE/'development-feature-manifest.json').read_text())
    for entry in manifest['features']:
        assert sha(ROOT/entry['cache_path'])==entry['sha256'],entry['cache_path']
    target=BASE/'full'/f'seed{args.seed}'
    if target.exists():raise FileExistsError('Existing run preserved: '+str(target))
    target.mkdir(parents=True)
    env=os.environ.copy()
    env.update(CATPRED_CACHE_PATH=str(BASE/'cache'),TORCH_HOME='/mnt/d/DisorderFlowRuntime/cache/torch',
        PYTHONPYCACHEPREFIX=str(BASE/'pycache'),MPLCONFIGDIR=str(BASE/'matplotlib'),
        OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',CATPRED_ESM_BATCH_SIZE='2')
    data=BASE/'development-data'
    command=[sys.executable,'-u',str(BASE/'upstream/train.py'),
        '--protein_records_path',str(data/'protein_records.json.gz'),
        '--data_path',str(data/'train.csv'),'--separate_val_path',str(data/'validation.csv'),
        '--separate_test_path',str(data/'validation.csv'),
        '--dataset_type','regression','--smiles_columns','reactant_smiles','--target_columns','log10kcat_max',
        '--extra_metrics','mae','mse','r2','--ensemble_size',str(protocol['ensemble_size']),
        '--seq_embed_dim','36','--seq_self_attn_nheads','6','--loss_function','mve',
        '--batch_size',str(protocol['batch_size']),'--epochs',str(protocol['epochs']),
        '--add_esm_feats','--seed',str(args.seed),'--pytorch_seed',str(args.seed),
        '--save_dir',str(target/'run')]
    started=time.time()
    receipt={'status':'running','seed':args.seed,'protocol_sha256':sha(protocol_path),
        'command':command,'started_unix':started,'evaluation_scope':'upstream test alias is validation only; internal test unopened'}
    status=target/'execution.json';status.write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt),flush=True)
    with (target/'training.log').open('w') as log:
        process=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt.update(status='completed' if process.returncode==0 else 'failed',returncode=process.returncode,elapsed_seconds=time.time()-started)
    receipt['checkpoints']=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in sorted(target.rglob('model.pt'))]
    if process.returncode==0 and len(receipt['checkpoints'])!=protocol['ensemble_size']:
        receipt['status']='failed';receipt['error']='Unexpected checkpoint count'
    status.write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
    if receipt['status']!='completed':raise RuntimeError('Training did not complete; inspect receipt')

if __name__=='__main__':main()
