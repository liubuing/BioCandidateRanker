"""Verify mmap values and simultaneously map the complete development cache."""
from pathlib import Path
import json,os,time,hashlib

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/catpred-mode-b'
os.environ['CATPRED_CACHE_PATH']=str(BASE/'cache')

def main():
    import torch,psutil
    from catpred_memory_entry import install_memory_adapter
    import catpred.data.cache_utils as cache
    manifest=json.loads((BASE/'development-feature-manifest.json').read_text())
    entries=sorted(manifest['features'],key=lambda r:r['shape'][0])
    install_memory_adapter()
    probes=[]
    for entry in [entries[0],entries[len(entries)//2],entries[-1]]:
        path=ROOT/entry['cache_path']
        ordinary=torch.load(path,map_location='cpu',weights_only=True)
        mapped=cache.load_torch_artifact(path,purpose='esm cache entry',map_location='cpu',roots=[BASE/'cache'])
        assert torch.equal(ordinary,mapped)
        probes.append({'path':entry['cache_path'],'shape':entry['shape'],'bitwise_equal':True})
        del ordinary,mapped
    tensors=[];started=time.monotonic();process=psutil.Process()
    for i,entry in enumerate(entries,1):
        tensor=cache.load_torch_artifact(ROOT/entry['cache_path'],purpose='esm cache entry',map_location='cpu',roots=[BASE/'cache'])
        assert list(tensor.shape)==entry['shape']
        tensors.append(tensor)
        if i%1000==0:print(i,process.memory_info().rss,flush=True)
    receipt={'status':'passed','entries':len(tensors),'logical_feature_bytes':sum(t.numel()*t.element_size() for t in tensors),'resident_bytes_after_mapping':process.memory_info().rss,'map_elapsed_seconds':time.monotonic()-started,'probes':probes,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Mapping feasibility and exact representative values; training memory still to be measured.'}
    (BASE/'mmap-verification.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
