"""Replay frozen ESM2 checkpoints solely to recover aligned predictions for inference."""
from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
import torch

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from biocandidate.training import load_checkpoint
from biocandidate.cli import _collect_task_outputs
from biocandidate.data.unikp import read_unikp_json,aggregate_pair_measurements
from biocandidate.data.homology import apply_split_manifest
from biocandidate.evaluation import regression_metrics

OUT=ROOT/'artifacts/paired-family-analysis'
SOURCE=Path(r'C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    protocol=read(ROOT/'configs/paired_family_analysis_protocol.json')
    for entry in protocol['frozen_files']:
        if sha(ROOT/entry['path'])!=entry['sha256']:raise ValueError('Frozen file changed: '+entry['path'])
    source_manifest=read(ROOT/'configs/unikp_source_manifest.json')
    assert sha(SOURCE)==source_manifest['identity']['sha256']
    split=read(ROOT/'artifacts/homology-final/homology_split.json')
    indices={r['source_row'] for r in split['rows'] if r['split']=='test'}
    records=read_unikp_json(SOURCE).records
    records=apply_split_manifest(records,ROOT/'artifacts/homology-final/homology_split.json',source_identity=source_manifest['identity'])
    records=[r for r in aggregate_pair_measurements(records) if r.split=='test']
    with (ROOT/'artifacts/external/sota-homology-cold/data/test.csv').open(encoding='utf-8',newline='') as f:reference=list(csv.DictReader(f))
    assert len(records)==len(reference)==1646
    for record,row in zip(records,reference):
        assert record.sequence==row['sequence'] and record.substrate_smiles==row['smiles']
        assert abs(record.log10_kcat-float(row['log10_kcat']))<1e-10
    OUT.mkdir(parents=True,exist_ok=True)
    device=torch.device('cuda')
    for seed in protocol['seeds']:
        target=OUT/f'reference-seed{seed}.csv';receipt_path=OUT/f'reference-seed{seed}-receipt.json'
        ckpt=ROOT/f'artifacts/esm2-t6-opt-seed{seed}/best.pt'
        if receipt_path.exists():
            receipt=read(receipt_path)
            assert receipt['predictions_sha256']==sha(target) and receipt['checkpoint_sha256']==sha(ckpt)
            continue
        with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):
            model,payload=load_checkpoint(ckpt,device)
        # Relocated checkout paths differ; authenticate source/split content, not old absolute paths.
        assert payload['data_manifest']['source']['identity']==source_manifest['identity']
        assert payload['data_manifest']['split']['identity']['sha256']==sha(ROOT/'artifacts/homology-final/homology_split.json')
        pred,_,labels=_collect_task_outputs(model,records,'log10_kcat',32,device)
        metrics=regression_metrics(pred,labels)
        prior=read(ROOT/f'artifacts/esm2-t6-opt-seed{seed}/test_metrics.json')['model_metrics']['log10_kcat']
        for key in ['rmse','mae','pearson']:
            if abs(metrics[key]-prior[key])>protocol['replay_metric_tolerance']:
                raise ValueError(f'Replay differs from frozen metrics: {seed} {key}')
        with target.open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['pair_id','prediction','label'])
            w.writerows((r['pair_id'],float(p),float(y)) for r,p,y in zip(reference,pred,labels))
        receipt={'checkpoint_sha256':sha(ckpt),'predictions_sha256':sha(target),'metrics':metrics,
                 'purpose':'post-hoc paired analysis; fixed checkpoint replay, no fitting or selection',
                 'protocol_sha256':sha(ROOT/'configs/paired_family_analysis_protocol.json')}
        receipt_path.write_text(json.dumps(receipt,indent=2),encoding='utf-8')
        print(seed,metrics,flush=True);del model;torch.cuda.empty_cache()

if __name__=='__main__':main()
