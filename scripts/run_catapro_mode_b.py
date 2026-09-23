"""CataPro kcat head on frozen train/validation data; no test evaluation command."""
from __future__ import annotations
import argparse
import csv
import gc
import hashlib
import importlib.util
import json
import random
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import MACCSkeys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/catapro-mode-b'
PROTOCOL=ROOT/'configs/catapro_mode_b_protocol.json'
PROTT5=Path(r'C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment\04_prediction_plugins\UniKP\models\prot_t5_xl_uniref50')
MOLT5=OUT/'models/molt5-base-smiles2caption'
LEGACY=ROOT/'artifacts/external/sota-homology-cold/unikp-mode-b/embedding-cache'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))

def save_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8');temp.replace(path)

def json_hash(value):return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()

def load_partitions(limit=0):
    protocol=read_json(PROTOCOL)
    manifest=read_json(ROOT/protocol['data_manifest'])
    if sha(ROOT/protocol['data_manifest'])!=protocol['data_manifest_sha256']:
        raise ValueError('Frozen data manifest changed')
    partitions={}
    for name in ['train','validation']:
        entry=manifest['partitions'][name];p=ROOT/entry['path']
        if sha(p)!=entry['sha256']:raise ValueError(f'{name} content changed')
        with p.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
        if len(rows)!=entry['rows']:raise ValueError('Row count changed')
        partitions[name]=rows[:limit] if limit else rows
    return partitions

def sequence_text(sequence):
    sequence=sequence if len(sequence)<=1000 else sequence[:500]+sequence[-500:]
    return re.sub('[UZOB]','X',' '.join(sequence))

def model_identity(path):
    names=['config.json','pytorch_model.bin','spiece.model','tokenizer_config.json','special_tokens_map.json','tokenizer.json']
    return {n:sha(path/n) for n in names if (path/n).is_file()}

def source_module(relative):
    receipt=read_json(ROOT/'artifacts/pcb-readiness/source-receipt.json')
    repo=next(r for r in receipt['repositories'] if r['repository']=='zchwang/CataPro')
    if repo['commit']!=read_json(PROTOCOL)['upstream_commit']:raise ValueError('Upstream commit changed')
    item=next(f for f in repo['files'] if f['path']==relative)
    path=ROOT/'artifacts/pcb-readiness/CataPro/source'/relative
    if sha(path)!=item['sha256']:raise ValueError('Upstream code changed')
    spec=importlib.util.spec_from_file_location('catapro_'+Path(relative).stem,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def cache_features(limit=0,device='cuda'):
    from transformers import T5EncoderModel,T5Tokenizer
    started=time.perf_counter();rows=load_partitions(limit)
    mode='smoke' if limit else 'full';dest=OUT/mode;dest.mkdir(parents=True,exist_ok=True)
    prot_identity=model_identity(PROTT5)
    if prot_identity.get('pytorch_model.bin')!=read_json(PROTOCOL)['prott5_weight_sha256']:
        raise ValueError('ProtT5 weights differ from recorded UniKP encoder')
    mol_identity=model_identity(MOLT5)
    download=read_json(MOLT5/'download-receipt.json')
    if download['revision']!=read_json(PROTOCOL)['molt5_revision']:raise ValueError('MolT5 revision changed')
    for name,value in mol_identity.items():
        if value!=download['files'][name]['sha256']:raise ValueError('MolT5 file changed')
    key_path=LEGACY/'prott5-index.json';vec_path=LEGACY/'prott5-embeddings.npy'
    keys=read_json(key_path);vectors=np.load(vec_path,allow_pickle=False)
    if vectors.shape!=(len(keys),1024) or len(set(keys))!=len(keys):raise ValueError('Invalid ProtT5 cache')
    lookup={k:i for i,k in enumerate(keys)}
    selected=sorted({r['sequence'] for part in rows.values() for r in part})
    if any(k not in lookup for k in selected):raise ValueError('Missing cached protein')
    train_keys=sorted({r['sequence'] for r in load_partitions(0)['train']},key=lambda x:(len(x),x))
    probes=list(dict.fromkeys([train_keys[0],train_keys[len(train_keys)//2],train_keys[-1]]))
    check_identity={'prot_files':prot_identity,'keys_sha256':sha(key_path),'vectors_sha256':sha(vec_path),
                    'probe_sequence_hashes':[hashlib.sha256(s.encode()).hexdigest() for s in probes],
                    'protocol_sha256':sha(PROTOCOL)}
    gate_path=OUT/'prott5-reuse-check.json'
    if gate_path.exists() and read_json(gate_path).get('identity')==check_identity:
        gate=read_json(gate_path)
        if not gate.get('passed'):raise ValueError('Prior ProtT5 compatibility failed')
    else:
        # Select probes from training inputs only, without using labels.
        tok=T5Tokenizer.from_pretrained(PROTT5,do_lower_case=False,local_files_only=True)
        encoder=T5EncoderModel.from_pretrained(PROTT5,local_files_only=True).to(device).eval()
        errors=[]
        with torch.inference_mode():
            for seq in probes:
                batch=tok([sequence_text(seq)],add_special_tokens=True,padding=True,return_tensors='pt').to(device)
                hidden=encoder(**batch).last_hidden_state[0].cpu().numpy()
                length=int(batch['attention_mask'][0].sum())-1
                fresh=hidden[:length].astype(np.float64).mean(axis=0)
                errors.append(float(np.max(np.abs(fresh-vectors[lookup[seq]]))))
        del encoder,tok;gc.collect()
        if torch.cuda.is_available():torch.cuda.empty_cache()
        gate={'identity':check_identity,'probe_lengths':[len(s) for s in probes],
              'max_absolute_errors':errors,'tolerance':1e-4,'passed':max(errors)<=1e-4,
              'scope':'training-only representative probes; same preprocessing and weight hash; not exhaustive equality'}
        save_json(gate_path,gate)
        if not gate['passed']:raise ValueError('ProtT5 cache differs; re-extraction required')
    identity={'protocol_sha256':sha(PROTOCOL),'prott5_reuse_check_sha256':sha(gate_path),
              'molt5_files':mol_identity,'representation':'ProtT5 mean1024 + MolT5 mean768 + MACCS167; float32'}
    namespace=json_hash(identity);cache=OUT/'molecule-cache'/namespace;cache.mkdir(parents=True,exist_ok=True)
    smiles=sorted({r['smiles'] for part in rows.values() for r in part})
    tokenizer=None;encoder=None;mol_vectors={}
    for i,smile in enumerate(smiles):
        key=hashlib.sha256(smile.encode()).hexdigest();path=cache/(key+'.npy');meta=cache/(key+'.json')
        if path.exists() and meta.exists():
            info=read_json(meta)
            if info.get('smiles')!=smile or sha(path)!=info['sha256']:raise ValueError('Molecule cache mismatch')
            value=np.load(path,allow_pickle=False)
        else:
            if encoder is None:
                tokenizer=T5Tokenizer.from_pretrained(MOLT5,local_files_only=True)
                encoder=T5EncoderModel.from_pretrained(MOLT5,local_files_only=True).to(device).eval()
            ids=tokenizer(smile,return_tensors='pt').input_ids.to(device)
            with torch.inference_mode():value=encoder(input_ids=ids).last_hidden_state[0,:-1].mean(dim=0).cpu().numpy()
            mol=Chem.MolFromSmiles(smile)
            if mol is None:raise ValueError('Invalid molecule')
            value=np.concatenate([value,np.asarray(list(MACCSkeys.GenMACCSKeys(mol).ToBitString()),dtype=np.float32)]).astype(np.float32)
            np.save(path,value);save_json(meta,{'smiles':smile,'sha256':sha(path),'identity':namespace})
        if value.shape!=(935,) or not np.isfinite(value).all():raise ValueError('Invalid molecule feature')
        mol_vectors[smile]=value
        if (i+1)%50==0:print(f'MolT5 {i+1}/{len(smiles)}',flush=True)
    if encoder is not None:del encoder,tokenizer
    gc.collect()
    if torch.cuda.is_available():torch.cuda.empty_cache()
    manifest={'mode':mode,'limit':limit,'identity':identity,'partitions':{},'test_accessed':False}
    for name,part in rows.items():
        features=np.stack([np.concatenate([vectors[lookup[r['sequence']]],mol_vectors[r['smiles']]]) for r in part]).astype(np.float32)
        path=dest/f'{name}-features.npy';np.save(path,features)
        manifest['partitions'][name]={'row_ids':[r['pair_id'] for r in part],'path':str(path.relative_to(ROOT)),
                                      'sha256':sha(path),'shape':list(features.shape)}
    manifest['seconds']=time.perf_counter()-started
    save_json(dest/'feature-manifest.json',manifest)
    print(json.dumps({'status':'cached','mode':mode,'sizes':{k:len(v) for k,v in rows.items()},'seconds':manifest['seconds']}),flush=True)

def load_feature_data(mode):
    manifest=read_json(OUT/mode/'feature-manifest.json')
    if manifest['identity']['protocol_sha256']!=sha(PROTOCOL):raise ValueError('Protocol changed')
    rows=load_partitions(manifest['limit']);result={}
    for name,part in rows.items():
        item=manifest['partitions'][name];path=ROOT/item['path']
        if [r['pair_id'] for r in part]!=item['row_ids']:raise ValueError('Feature row ordering changed')
        if sha(path)!=item['sha256']:raise ValueError('Features changed')
        features=np.load(path,allow_pickle=False)
        if features.shape!=(len(part),1959) or not np.isfinite(features).all():raise ValueError('Bad feature dimensions')
        labels=np.asarray([float(r['log10_kcat']) for r in part],dtype=np.float32)
        if not np.isfinite(labels).all():raise ValueError('Invalid target')
        result[name]=(features,labels)
    return result,manifest

def train(mode,device):
    protocol=read_json(PROTOCOL);cfg=protocol['training'];data,feature_manifest=load_feature_data(mode)
    if mode=='full':
        check=read_json(OUT/'molecule-equivalence-check.json')
        if not check.get('passed') or check.get('feature_identity')!=feature_manifest['identity']:
            raise ValueError('Full training requires matching upstream molecule equivalence check')
    module=source_module('training/kcat/model.py')
    max_epochs=2 if mode=='smoke' else cfg['max_epochs'];seeds=[7] if mode=='smoke' else cfg['seeds']
    for seed in seeds:
        dest=OUT/mode/f'seed{seed}';dest.mkdir(parents=True,exist_ok=True)
        if (dest/'receipt.json').exists():
            old=read_json(dest/'receipt.json')
            if old['protocol_sha256']!=sha(PROTOCOL) or old['feature_manifest_sha256']!=sha(OUT/mode/'feature-manifest.json'):
                raise ValueError('Completed run has incompatible identities')
            if sha(dest/'best.pt')!=old['checkpoint_sha256']:raise ValueError('Completed checkpoint changed')
            print(f'Seed {seed} already complete',flush=True);continue
        random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
        if torch.cuda.is_available():torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.benchmark=False
        model=module.KcatModel(rate=cfg['dropout'],device=device)
        optimizer=torch.optim.Adam(model.parameters(),lr=cfg['learning_rate'],betas=(.9,.999),weight_decay=cfg['weight_decay'])
        tensors={k:(torch.from_numpy(x).to(device),torch.from_numpy(y).to(device)) for k,(x,y) in data.items()}
        best=float('inf');stale=0;logs=[];start=time.perf_counter()
        if device.startswith('cuda'):torch.cuda.reset_peak_memory_stats()
        for epoch in range(max_epochs):
            model.train();x,y=tensors['train'];order=torch.randperm(len(y),device=device)
            batches=list(order.split(cfg['batch_size']))
            if len(batches[-1])==1:
                if len(batches)<2:raise ValueError('BatchNorm requires at least two training rows')
                batches[-2]=torch.cat([batches[-2],batches[-1]]);batches.pop()
            for indices in batches:
                optimizer.zero_grad();p=model(x[indices,:1024],x[indices,1024:]).ravel()
                loss=torch.sqrt(torch.mean((p-y[indices])**2)+1e-6)
                if not torch.isfinite(loss):raise ValueError('Nonfinite loss')
                loss.backward();optimizer.step()
            model.eval();vx,vy=tensors['validation'];preds=[];losses=[]
            with torch.inference_mode():
                for s in range(0,len(vy),cfg['batch_size']):
                    p=model(vx[s:s+cfg['batch_size'],:1024],vx[s:s+cfg['batch_size'],1024:]).ravel()
                    preds.append(p);losses.append(float(torch.sqrt(torch.mean((p-vy[s:s+cfg['batch_size']])**2)+1e-6)))
            validation_loss=float(np.mean(losses));p=torch.cat(preds)
            rmse=float(torch.sqrt(torch.mean((p-vy)**2)));mae=float(torch.mean(torch.abs(p-vy)))
            if not np.isfinite(validation_loss):raise ValueError('Nonfinite validation')
            improved=validation_loss<best-cfg['min_delta']
            if improved:
                best=validation_loss;stale=1 if epoch==0 else 0;best_epoch=epoch
                torch.save(model.state_dict(),dest/'best.pt')
            else:stale+=1
            logs.append({'epoch':epoch,'validation_batch_mean_rmse_loss':validation_loss,'validation_rmse':rmse,'validation_mae':mae})
            save_json(dest/'history.json',logs)
            print(f'{mode} seed={seed} epoch={epoch} validation_rmse={rmse:.4f}',flush=True)
            if stale>=cfg['patience']:break
        result={'mode':mode,'seed':seed,'upstream_commit':protocol['upstream_commit'],'protocol_sha256':sha(PROTOCOL),
                'feature_manifest_sha256':sha(OUT/mode/'feature-manifest.json'),'checkpoint_sha256':sha(dest/'best.pt'),
                'best_epoch':best_epoch,'epochs':len(logs),'best_validation_loss':best,'seconds':time.perf_counter()-start,
                'rows':{k:len(v[1]) for k,v in data.items()},'test_accessed':False,
                'runner_sha256':sha(Path(__file__)),
                'peak_allocated_gpu_bytes':torch.cuda.max_memory_allocated() if device.startswith('cuda') else None,
                'torch_version':torch.__version__,'numpy_version':np.__version__}
        save_json(dest/'receipt.json',result)
    print('Training complete; test not evaluated.',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['cache','train'])
    p.add_argument('--limit',type=int,default=0);p.add_argument('--mode',choices=['smoke','full'],default='smoke')
    p.add_argument('--device',default='cuda');args=p.parse_args()
    if args.limit<0:raise ValueError('limit cannot be negative')
    if args.action=='cache':cache_features(args.limit,args.device)
    else:train(args.mode,args.device)

if __name__=='__main__':main()
