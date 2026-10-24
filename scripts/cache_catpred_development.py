"""Bounded-memory official CatPred ESM features, development data only (WSL).

Reads development-data/protein_records.json.gz (verified against receipt.json),
computes official CatPred ESM2 features for every distinct sequence in small
CPU batches through the WSL upstream tree, and asserts shape and finiteness on
each cache entry. Writes
artifacts/catpred-mode-b/development-feature-manifest.json, which records this
script's own SHA256 and the training-and-validation-only scope.

中文：本脚本经 WSL 上游树以有界内存计算官方 CatPred ESM 特征（仅限开发数据）：读取
development-data/protein_records.json.gz（对照 receipt.json 校验），用小 CPU 批次为每条
去重序列计算官方 CatPred ESM2 特征，并对每个缓存条目断言形状与有限性。写出
artifacts/catpred-mode-b/development-feature-manifest.json，其中记录本脚本自身 SHA256 与
“仅训练与验证”的范围。治理角色：范围受限的特征缓存清单收据（预处理证据，不产生预测）。
"""
from pathlib import Path
import hashlib,json,gzip,time,sys,os

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/catpred-mode-b'
os.environ.setdefault('CATPRED_CACHE_PATH',str(BASE/'cache'))
os.environ.setdefault('CATPRED_ESM_BATCH_SIZE','2')
sys.path.insert(0,str(BASE/'upstream'))

def sha(p):
    """Return the SHA256 hex digest of a file.

    中文：返回文件的 SHA256 十六进制摘要。
    """
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    """Cache ESM features for all development sequences and write the manifest.

    中文：为全部开发序列缓存 ESM 特征并写出清单。
    """
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
