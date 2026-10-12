"""Run three frozen current-code training jobs, without test evaluation (Windows)."""
from pathlib import Path
import hashlib, json, os, subprocess, sys, time, msvcrt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/selection-rule-comparison'
PROTOCOL=ROOT/'configs/selection_rule_comparison_protocol.json'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def write(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8');os.replace(tmp,path)

def main():
    BASE.mkdir(exist_ok=True)
    lock=(BASE/'queue.lock').open('a+b');lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0)
    msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    p=json.loads(PROTOCOL.read_text(encoding='utf-8'))
    assert sha(__file__)==p['runner_sha256']
    for e in p['frozen_files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
    assert sha(p['source_path'])==p['source_sha256']
    for e in p['weights']:assert sha(e['path'])==e['sha256']
    env=os.environ.copy();env['PYTHONPATH']=str(BASE/'source');env['TORCH_HOME']=p['torch_home'];env['PYTHONUNBUFFERED']='1'
    env['OMP_NUM_THREADS']='4';env['MKL_NUM_THREADS']='4'
    completed=[]
    for seed in p['seeds']:
        folder=BASE/f'seed{seed}';folder.mkdir(exist_ok=True);receipt=folder/'execution.json'
        if receipt.exists():
            old=json.loads(receipt.read_text())
            if old['status']=='completed':
                assert old['protocol_sha256']==sha(PROTOCOL) and sha(folder/'best.pt')==old['checkpoint_sha256']
                assert sha(folder/'best_rmse.pt')==old['rmse_checkpoint_sha256']
                completed.append(old);continue
            raise RuntimeError(f'Incomplete run needs recovery audit: {folder}')
        for e in p['frozen_files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
        cmd=[sys.executable,'-u','-m','biocandidate.cli','train','--data',p['source_path'],
             '--manifest',str(ROOT/'artifacts/unikp_audit.json'),'--split-manifest',str(ROOT/'artifacts/homology-final/homology_split.json'),
             '--conflict-policy','median','--protein-encoder','esm2','--esm2-model','esm2_t6_8M_UR50D',
             '--d-model','256','--num-heads','8','--protein-layers','3','--molecule-layers','4','--fusion-layers','2',
             '--chunk-size','128','--context-buckets','4096','--fusion-mode','task_query','--uncertainty-mode','heteroscedastic',
             '--lr','0.0003','--batch-size','4','--grad-accum-steps','8','--epochs','10','--warmup-epochs','2','--patience','0',
             '--seed',str(seed),'--cosine-schedule','--device','cuda','--output-dir',str(folder)]
        r={'status':'running','seed':seed,'started_unix':time.time(),'protocol_sha256':sha(PROTOCOL),'command':cmd}
        write(receipt,r);print('Starting training',seed,flush=True)
        with (folder/'train.log').open('w',encoding='utf-8') as log:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        r.update(status='completed' if result.returncode==0 else 'failed',returncode=result.returncode,finished_unix=time.time())
        if result.returncode==0:
            r['checkpoint_sha256']=sha(folder/'best.pt')
            r['rmse_checkpoint_sha256']=sha(folder/'best_rmse.pt')
        write(receipt,r)
        if result.returncode:raise RuntimeError(f'Training failed: seed{seed}')
        completed.append(r)
    write(BASE/'completed-training-manifest.json',{'protocol_sha256':sha(PROTOCOL),'runs':completed,'scope':'Selection rule comparison training complete; no held-out evaluation performed'})

if __name__=='__main__':main()
