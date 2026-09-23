"""Atomic epoch-boundary state for single-member CatPred runs."""
from pathlib import Path
import os,random
import numpy as np
import torch

def pack(value):
    if isinstance(value,np.ndarray):return {'__array__':value.tolist(),'dtype':str(value.dtype)}
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:pack(v) for k,v in value.items()}
    if isinstance(value,tuple):return tuple(pack(v) for v in value)
    if isinstance(value,list):return [pack(v) for v in value]
    return value

def unpack(value):
    if isinstance(value,dict) and set(value)=={'__array__','dtype'}:return np.array(value['__array__'],dtype=value['dtype'])
    if isinstance(value,dict):return {k:unpack(v) for k,v in value.items()}
    if isinstance(value,tuple):return tuple(unpack(v) for v in value)
    if isinstance(value,list):return [unpack(v) for v in value]
    return value

def restore(folder,model,optimizer,scheduler,loader):
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
