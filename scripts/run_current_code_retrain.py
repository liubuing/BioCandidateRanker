"""Run the three frozen current-code training jobs; no test evaluation (Windows).

Reads configs/current_code_retrain_protocol.json, which pins this runner's own
SHA256, every frozen input file, the training source tree, the pretrained
weights and the seed set (7, 42, 123). For each seed it re-verifies all
hashes, launches ``biocandidate.cli train`` as a subprocess with the frozen
configuration, and writes an execution.json receipt plus train.log under
artifacts/current-code-retrain/seed<seed>/; completed-training-manifest.json
is written only after every seed finishes. Governance: an exclusive msvcrt
file lock serialises runs; a completed receipt short-circuits that seed while
an incomplete receipt fails closed and demands a manual recovery audit; the
held-out test partition is never evaluated.
"""
from pathlib import Path
import hashlib, json, os, subprocess, sys, time, msvcrt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/current-code-retrain'
PROTOCOL=ROOT/'configs/current_code_retrain_protocol.json'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def write(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8');os.replace(tmp,path)

def main():
    """Verify the freeze and run one training job per protocol seed, fail-closed throughout."""
    BASE.mkdir(exist_ok=True)
    # Single-byte non-blocking lock: a second concurrent retrain on this
    # machine would race the per-seed receipts, so refuse to start.
    lock=(BASE/'queue.lock').open('a+b');lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0)
    msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
    p=json.loads(PROTOCOL.read_text(encoding='utf-8'))
    # The protocol pins this runner's own hash, so any edit invalidates the freeze.
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
                completed.append(old);continue
            # An interrupted run must be audited by hand; never silently restart it.
            raise RuntimeError(f'Incomplete run needs recovery audit: {folder}')
        for e in p['frozen_files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
        cmd=[sys.executable,'-u','-m','biocandidate.cli','train','--data',p['source_path'],
             '--manifest',str(ROOT/'artifacts/unikp_audit.json'),'--split-manifest',str(ROOT/'artifacts/homology-final/homology_split.json'),
             '--conflict-policy','median','--protein-encoder','esm2','--esm2-model','esm2_t6_8M_UR50D',
             '--d-model','256','--num-heads','8','--protein-layers','3','--molecule-layers','4','--fusion-layers','2',
             '--chunk-size','128','--context-buckets','4096','--fusion-mode','task_query','--uncertainty-mode','heteroscedastic',
             '--lr','0.0003','--batch-size','4','--grad-accum-steps','8','--epochs','10','--warmup-epochs','2','--patience','4',
             '--seed',str(seed),'--cosine-schedule','--device','cuda','--output-dir',str(folder)]
        r={'status':'running','seed':seed,'started_unix':time.time(),'protocol_sha256':sha(PROTOCOL),'command':cmd}
        write(receipt,r);print('Starting training',seed,flush=True)
        with (folder/'train.log').open('w',encoding='utf-8') as log:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        r.update(status='completed' if result.returncode==0 else 'failed',returncode=result.returncode,finished_unix=time.time())
        if result.returncode==0:r['checkpoint_sha256']=sha(folder/'best.pt')
        write(receipt,r)
        if result.returncode:raise RuntimeError(f'Training failed: seed{seed}')
        completed.append(r)
    write(BASE/'completed-training-manifest.json',{'protocol_sha256':sha(PROTOCOL),'runs':completed,'scope':'Current-code retraining complete; no held-out evaluation performed'})

if __name__=='__main__':main()
