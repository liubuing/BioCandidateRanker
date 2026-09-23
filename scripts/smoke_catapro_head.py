"""Exercise reviewed, pinned CataPro head on synthetic tensors only."""
from pathlib import Path
import hashlib
import importlib.util
import json
import time
import torch

root=Path(__file__).resolve().parents[1]
out=root/'artifacts/pcb-readiness'
source=out/'CataPro/source/training/kcat/model.py'
receipt=json.loads((out/'source-receipt.json').read_text())
repo=next(r for r in receipt['repositories'] if r['repository']=='zchwang/CataPro')
expected=next(f['sha256'] for f in repo['files'] if f['path']=='training/kcat/model.py')
assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
spec=importlib.util.spec_from_file_location('reviewed_catapro_head',source)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
torch.manual_seed(7)
device='cuda' if torch.cuda.is_available() else 'cpu'
model=module.KcatModel(rate=0.,device=device)
optimizer=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=.01)
if device=='cuda':torch.cuda.reset_peak_memory_stats();torch.cuda.synchronize()
started=time.perf_counter()
x=torch.randn(8,1024,device=device);s=torch.randn(8,935,device=device);y=torch.randn(8,1,device=device)
optimizer.zero_grad();pred=model(x,s);loss=torch.sqrt(torch.mean((pred-y)**2)+1e-6)
assert pred.shape==(8,1) and torch.isfinite(loss)
loss.backward();optimizer.step();model.eval()
with torch.no_grad():test=model(x,s)
assert torch.isfinite(test).all()
if device=='cuda':torch.cuda.synchronize()
result={'status':'pass','synthetic_only':True,'real_dataset_predictions':False,
        'upstream_commit':repo['commit'],'source_sha256':expected,'device':device,
        'parameters':sum(p.numel() for p in model.parameters()),
        'elapsed_seconds':time.perf_counter()-started,
        'peak_allocated_gpu_bytes':torch.cuda.max_memory_allocated() if device=='cuda' else None,
        'limitation':'Head forward/backward only; excludes ProtT5/MolT5 feature extraction, checkpoints, and full-data training.'}
(out/'catapro-head-smoke.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
