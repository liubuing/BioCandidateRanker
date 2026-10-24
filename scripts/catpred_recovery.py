"""Atomic epoch-boundary state for single-member CatPred runs.

save() persists complete epoch-boundary state (model, optimizer, scheduler,
sampler and Python/NumPy/Torch/CUDA RNG states, plus the latest best
checkpoint bytes) atomically via a tmp file and os.replace; restore() replays
it so a rerun is deterministic. The saved identity must match the
CATPRED_RECOVERY_ID environment variable, and CATPRED_TEST_STOP_EPOCH makes
save() exit 75 to rehearse interruption at a chosen epoch.

中文：本模块为单成员 CatPred 运行提供原子化的 epoch 边界状态：save() 通过 tmp 文件与
os.replace 原子持久化完整的 epoch 边界状态（模型、优化器、调度器、采样器与
Python/NumPy/Torch/CUDA 的 RNG 状态，以及最新 best checkpoint 字节）；restore() 完整
重放该状态使重跑确定。保存的身份必须与 CATPRED_RECOVERY_ID 环境变量一致；
CATPRED_TEST_STOP_EPOCH 使 save() 以退出码 75 退出，用于在选定 epoch 演练中断。
治理角色：确定性恢复协议（供重启后精确续训）。
"""
from pathlib import Path
import os,random
import numpy as np
import torch

def pack(value):
    """Recursively convert NumPy arrays/scalars into JSON-safe structures.

    中文：递归地把 NumPy 数组/标量转换为 JSON 安全的结构。
    """
    if isinstance(value,np.ndarray):return {'__array__':value.tolist(),'dtype':str(value.dtype)}
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:pack(v) for k,v in value.items()}
    if isinstance(value,tuple):return tuple(pack(v) for v in value)
    if isinstance(value,list):return [pack(v) for v in value]
    return value

def unpack(value):
    """Inverse of pack: rebuild NumPy arrays and containers from saved structures.

    中文：pack 的逆操作：从已保存的结构重建 NumPy 数组与容器。
    """
    if isinstance(value,dict) and set(value)=={'__array__','dtype'}:return np.array(value['__array__'],dtype=value['dtype'])
    if isinstance(value,dict):return {k:unpack(v) for k,v in value.items()}
    if isinstance(value,tuple):return tuple(unpack(v) for v in value)
    if isinstance(value,list):return [unpack(v) for v in value]
    return value

def restore(folder,model,optimizer,scheduler,loader):
    """Load recovery.pt if present; return (start_epoch, best_score, best_epoch, n_iter).

    中文：若存在 recovery.pt 则加载并重放全部状态；返回
    (start_epoch, best_score, best_epoch, n_iter)。
    """
    path=Path(folder)/'recovery.pt'
    if not path.exists():return 0,float('inf'),0,0
    state=unpack(torch.load(path,map_location='cpu',weights_only=True))
    assert state['identity']==os.environ['CATPRED_RECOVERY_ID']
    model.load_state_dict(state['model']);optimizer.load_state_dict(state['optimizer'])
    scheduler.load_state_dict(unpack(state['scheduler']))
    loader._sampler._random.setstate(state['sampler_rng'])
    random.setstate(state['python_rng']);np.random.set_state(unpack(state['numpy_rng']))
    torch.set_rng_state(state['torch_rng'])
    if state['cuda_rng'] is not None:torch.cuda.set_rng_state_all(state['cuda_rng'])
    best=Path(folder)/'model.pt';tmp=best.with_suffix('.restore.tmp');tmp.write_bytes(state['best_checkpoint'].numpy().tobytes());os.replace(tmp,best)
    print('RECOVERY restored through epoch',state['epoch'],flush=True)
    return state['epoch']+1,state['best_score'],state['best_epoch'],state['n_iter']

def save(folder,model,optimizer,scheduler,loader,epoch,best_score,best_epoch,n_iter):
    """Atomically persist all epoch-boundary state after each epoch.

    中文：每个 epoch 结束后原子持久化全部 epoch 边界状态。
    """
    state={'identity':os.environ['CATPRED_RECOVERY_ID'],'epoch':epoch,'model':model.state_dict(),
        'optimizer':optimizer.state_dict(),'scheduler':pack(scheduler.state_dict()),
        'sampler_rng':loader._sampler._random.getstate(),'python_rng':random.getstate(),
        'numpy_rng':pack(np.random.get_state()),'torch_rng':torch.get_rng_state(),
        'cuda_rng':torch.cuda.get_rng_state_all() if next(model.parameters()).is_cuda else None,
        'best_score':best_score,'best_epoch':best_epoch,'n_iter':n_iter,
        'best_checkpoint':torch.frombuffer(bytearray((Path(folder)/'model.pt').read_bytes()),dtype=torch.uint8)}
    path=Path(folder)/'recovery.pt';tmp=path.with_suffix('.tmp');torch.save(pack(state),tmp);os.replace(tmp,path)
    print('RECOVERY saved epoch',epoch,flush=True)
    if os.getenv('CATPRED_TEST_STOP_EPOCH')==str(epoch):raise SystemExit(75)
