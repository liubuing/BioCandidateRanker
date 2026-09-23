"""Single-lock queue of independently seeded members with epoch recovery."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,time,fcntl
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'artifacts/catpred-mode-b'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2));os.replace(tmp,path)

def main():
    root=BASE/'full-recoverable';root.mkdir(exist_ok=True)
    lock=(root/'queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    protocol_path=ROOT/'configs/catpred_recoverable_protocol.json';protocol=json.loads(protocol_path.read_text());identity=sha(protocol_path)
    assert sha(__file__)==protocol['runner_sha256']
    for entry in protocol['frozen_files']:assert sha(ROOT/entry['path'])==entry['sha256'],entry['path']
    features=json.loads((BASE/'development-feature-manifest.json').read_text())
    for i,entry in enumerate(features['features'],1):
        assert sha(ROOT/entry['cache_path'])==entry['sha256']
        if i%1000==0:print('Verified features',i,flush=True)
    retained=json.loads((BASE/'retained-member0.json').read_text());assert sha(ROOT/retained['checkpoint'])==retained['checkpoint_sha256']
    completed=[retained]
    data=BASE/'development-data'
    for group in protocol['seeds']:
        for member in range(10):
            if group==7 and member==0:continue
            folder=root/f'seed{group}'/f'member{member}';folder.mkdir(parents=True,exist_ok=True)
            receipt_path=folder/'execution.json'
            if receipt_path.exists():
                old=json.loads(receipt_path.read_text());assert old['protocol_sha256']==identity
                if old['status']=='completed':
                    assert sha(ROOT/old['checkpoint'])==old['checkpoint_sha256'];completed.append(old);continue
            initial_seed=group if member==0 else group*1000+member
            env=os.environ.copy();env.pop('CATPRED_TEST_STOP_EPOCH',None)
            env.update(CATPRED_RECOVERY_ID=f'{identity}:{group}:{member}',CATPRED_CACHE_PATH=str(BASE/'cache'),
                PYTHONPYCACHEPREFIX=str(BASE/'pycache'),MPLCONFIGDIR=str(BASE/'matplotlib'),
                TORCH_HOME='/mnt/d/DisorderFlowRuntime/cache/torch',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
            command=[sys.executable,'-u',str(ROOT/'scripts/catpred_recovery_entry.py'),
                '--protein_records_path',str(data/'protein_records.json.gz'),'--data_path',str(data/'train.csv'),
                '--separate_val_path',str(data/'validation.csv'),'--separate_test_path',str(data/'validation.csv'),
                '--dataset_type','regression','--smiles_columns','reactant_smiles','--target_columns','log10kcat_max',
                '--extra_metrics','mae','mse','r2','--ensemble_size','1','--seq_embed_dim','36','--seq_self_attn_nheads','6',
                '--loss_function','mve','--batch_size','16','--epochs','30','--add_esm_feats',
                '--seed',str(group),'--pytorch_seed',str(initial_seed),'--save_dir',str(folder/'run')]
            receipt={'status':'running','protocol_sha256':identity,'group_seed':group,'member':member,'initial_seed':initial_seed,'started_unix':time.time(),'scope':'validation only; independent member seeding amendment'}
            write(receipt_path,receipt);print('Starting/resuming',group,member,flush=True)
            with (folder/'training.log').open('a') as log:r=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
            receipt.update(status='completed' if r.returncode==0 else 'failed',returncode=r.returncode,elapsed_this_attempt=time.time()-receipt['started_unix'])
            checkpoint=folder/'run/fold_0/model_0/model.pt'
            if r.returncode==0:
                receipt.update(checkpoint=str(checkpoint.relative_to(ROOT)),checkpoint_sha256=sha(checkpoint));completed.append(receipt)
            write(receipt_path,receipt)
            if r.returncode:raise RuntimeError(f'Member failed {group}/{member}; state preserved')
    write(root/'completed-ensemble-manifest.json',{'protocol_sha256':identity,'members':completed,'scope':'Checkpoints complete, no held-out test evaluation yet'})

if __name__=='__main__':main()
