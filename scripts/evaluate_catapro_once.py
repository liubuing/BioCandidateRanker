"""Evaluate three pre-frozen CataPro checkpoints on the existing internal test once."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,sys,statistics
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import MACCSkeys
import run_catapro_mode_b as run

ROOT=run.ROOT
FREEZE=ROOT/'configs/catapro_internal_test_freeze.json'
DEST=run.OUT/'internal-test'
sys.path.insert(0,str(ROOT/'src'))
from biocandidate.evaluation import regression_metrics

def verify_files(entries,root=ROOT):
    for item in entries:
        if run.sha(root/item['path'])!=item['sha256']:
            raise ValueError(f"Frozen file changed: {item['path']}")

def features(rows,identity):
    from transformers import T5EncoderModel,T5Tokenizer
    gate=run.read_json(run.OUT/'prott5-reuse-check.json')
    if not gate['passed'] or run.sha(run.OUT/'prott5-reuse-check.json')!=identity['prott5_reuse_check_sha256']:
        raise ValueError('ProtT5 gate mismatch')
    keys=run.read_json(run.LEGACY/'prott5-index.json')
    if run.sha(run.LEGACY/'prott5-index.json')!=gate['identity']['keys_sha256'] or run.sha(run.LEGACY/'prott5-embeddings.npy')!=gate['identity']['vectors_sha256']:
        raise ValueError('ProtT5 cache changed')
    vectors=np.load(run.LEGACY/'prott5-embeddings.npy',allow_pickle=False);lookup={s:i for i,s in enumerate(keys)}
    if run.model_identity(run.MOLT5)!=identity['molt5_files']:raise ValueError('MolT5 changed')
    namespace=run.json_hash(identity);cache=run.OUT/'molecule-cache'/namespace
    molecules={};encoder=None;tokenizer=None
    for i,smile in enumerate(sorted({r['smiles'] for r in rows})):
        key=hashlib.sha256(smile.encode()).hexdigest();path=cache/(key+'.npy');meta=cache/(key+'.json')
        if path.exists() and meta.exists():
            receipt=run.read_json(meta)
            if receipt['smiles']!=smile or receipt['sha256']!=run.sha(path):raise ValueError('Molecule cache changed')
            value=np.load(path,allow_pickle=False)
        else:
            if encoder is None:
                tokenizer=T5Tokenizer.from_pretrained(run.MOLT5,local_files_only=True)
                encoder=T5EncoderModel.from_pretrained(run.MOLT5,local_files_only=True).to('cuda').eval()
            ids=tokenizer(smile,return_tensors='pt').input_ids.to('cuda')
            with torch.inference_mode():value=encoder(input_ids=ids).last_hidden_state[0,:-1].mean(dim=0).cpu().numpy()
            mol=Chem.MolFromSmiles(smile)
            if mol is None:raise ValueError('Invalid SMILES')
            value=np.concatenate([value,np.asarray(list(MACCSkeys.GenMACCSKeys(mol).ToBitString()),dtype=np.float32)]).astype(np.float32)
            np.save(path,value);run.save_json(meta,{'smiles':smile,'sha256':run.sha(path),'identity':namespace})
        if value.shape!=(935,) or not np.isfinite(value).all():raise ValueError('Invalid molecule features')
        molecules[smile]=value
    if encoder is not None:del encoder,tokenizer
    x=np.stack([np.concatenate([vectors[lookup[r['sequence']]],molecules[r['smiles']]]) for r in rows]).astype(np.float32)
    if x.shape!=(1646,1959) or not np.isfinite(x).all():raise ValueError('Test feature shape or values invalid')
    return x

def main():
    freeze=run.read_json(FREEZE);verify_files(freeze['files'])
    DEST.mkdir(parents=True,exist_ok=True)
    summary_path=DEST/'summary.json'
    if summary_path.exists():
        result=run.read_json(summary_path)
        if result['freeze_sha256']!=run.sha(FREEZE):raise ValueError('Existing evaluation has different freeze')
        verify_files(result['prediction_files'])
        print('Already evaluated; verified saved receipt, no model inference.',flush=True);return
    if (DEST/'started.json').exists():raise ValueError('Evaluation already started; audit partial outputs before recovery')
    full=run.read_json(run.OUT/'full/feature-manifest.json')
    if full['identity']['protocol_sha256']!=run.sha(run.PROTOCOL):raise ValueError('Training protocol changed')
    check=run.read_json(run.OUT/'molecule-equivalence-check.json')
    if not check['passed'] or check['feature_identity']!=full['identity']:raise ValueError('Feature equivalence gate failed')
    for seed in freeze['seeds']:
        receipt=run.read_json(run.OUT/f'full/seed{seed}/receipt.json')
        if receipt['test_accessed'] or receipt['mode']!='full':raise ValueError('Invalid training provenance')
        if receipt['checkpoint_sha256']!=run.sha(run.OUT/f'full/seed{seed}/best.pt'):raise ValueError('Checkpoint mismatch')
        if receipt['protocol_sha256']!=run.sha(run.PROTOCOL) or receipt['feature_manifest_sha256']!=run.sha(run.OUT/'full/feature-manifest.json'):
            raise ValueError('Training provenance mismatch')
    manifest=run.read_json(ROOT/freeze['data_manifest'])
    with (ROOT/manifest['partitions']['test']['path']).open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)!=freeze['test_rows']:raise ValueError('Unexpected row count')
    x=features(rows,full['identity']);np.save(DEST/'test-features.npy',x)
    run.save_json(DEST/'feature-receipt.json',{'identity':full['identity'],'row_ids':[r['pair_id'] for r in rows],
        'sha256':run.sha(DEST/'test-features.npy'),'frozen_encoders_only':True})
    module=run.source_module('training/kcat/model.py')
    # Exclusive marker prevents accidental repeat inference after outcomes are seen.
    with (DEST/'started.json').open('x',encoding='utf-8') as f:
        json.dump({'freeze_sha256':run.sha(FREEZE),'started_at':datetime.now(timezone.utc).isoformat()},f)
    xt=torch.from_numpy(x).to('cuda');y=torch.tensor([float(r['log10_kcat']) for r in rows])
    result={'freeze_sha256':run.sha(FREEZE),'claim_boundary':freeze['claim_boundary'],'seeds':{},'prediction_files':[]}
    for seed in freeze['seeds']:
        model=module.KcatModel(rate=0.,device='cuda')
        model.load_state_dict(torch.load(run.OUT/f'full/seed{seed}/best.pt',map_location='cuda',weights_only=True));model.eval()
        with torch.inference_mode():
            pred=torch.cat([model(xt[s:s+64,:1024],xt[s:s+64,1024:]).ravel().cpu() for s in range(0,len(rows),64)])
        if not torch.isfinite(pred).all():raise ValueError('Nonfinite prediction')
        metrics=regression_metrics(pred,y);result['seeds'][str(seed)]=metrics
        path=DEST/f'seed{seed}-predictions.csv'
        with path.open('w',encoding='utf-8',newline='') as f:
            writer=csv.writer(f);writer.writerow(['pair_id','prediction_log10_kcat','label_log10_kcat'])
            writer.writerows((r['pair_id'],float(p),float(v)) for r,p,v in zip(rows,pred,y))
        result['prediction_files'].append({'path':str(path.relative_to(ROOT)),'sha256':run.sha(path)})
        run.save_json(DEST/f'seed{seed}-metrics.json',metrics)
        print(seed,metrics,flush=True)
    result['mean_sd']={key:{'mean':statistics.mean(r[key] for r in result['seeds'].values()),
        'sd':statistics.stdev(r[key] for r in result['seeds'].values())} for key in ['rmse','mae','pearson']}
    result['completed_at']=datetime.now(timezone.utc).isoformat()
    run.save_json(summary_path,result);print(json.dumps(result['mean_sd'],indent=2))

if __name__=='__main__':main()
