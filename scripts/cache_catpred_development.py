"""Bounded-memory official CatPred ESM features, development data only (WSL)."""
from pathlib import Path
import hashlib,json,gzip,time,sys,os

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/catpred-mode-b'
os.environ.setdefault('CATPRED_CACHE_PATH',str(BASE/'cache'))
os.environ.setdefault('CATPRED_ESM_BATCH_SIZE','2')
sys.path.insert(0,str(BASE/'upstream'))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    import torch
    from catpred.data.esm_utils import get_many_esm_reprs,ESM_MAX_LENGTH,ESM_EMBED_DIM
    from catpred.data.cache_utils import cache_entry_path
    data=BASE/'development-data'
    receipt=json.loads((data/'receipt.json').read_text())
    source=data/'protein_records.json.gz'
    assert sha(source)==receipt['records_sha256']
    with gzip.open(source,'rt') as f:records=json.load(f)
    sequences=list(dict.fromkeys(r['seq'] for r in records.values()))
    result=[];started=time.monotonic()
    for offset in range(0,len(sequences),32):
        batch=sequences[offset:offset+32]
        features=get_many_esm_reprs(batch,device='cpu',batch_size=2)
        for seq in batch:
            value=features[seq]
            assert tuple(value.shape)==(min(len(seq),ESM_MAX_LENGTH-1),ESM_EMBED_DIM)
            assert torch.isfinite(value).all()
            path=cache_entry_path(path='esm/proteins',cache_key=seq)
            result.append({'sequence_sha256':hashlib.sha256(seq.encode('ascii')).hexdigest(),'cache_path':str(path.relative_to(ROOT)),'sha256':sha(path),'shape':list(value.shape)})
        del features
        elapsed=time.monotonic()-started
        print(json.dumps({'completed':len(result),'total':len(sequences),'seconds':round(elapsed,2),'peak_gpu_allocated_bytes':torch.cuda.max_memory_allocated()}),flush=True)
    (BASE/'development-feature-manifest.json').write_text(json.dumps({'scope':'training and validation only; official ESM2 cache','input_sha256':sha(source),'script_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-started,'features':result},indent=2))

if __name__=='__main__':main()
