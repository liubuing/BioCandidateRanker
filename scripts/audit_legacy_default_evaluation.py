"""Reconstruct historical evaluation using the CLI's default batch size 64."""
import importlib.util,json,csv
from pathlib import Path
import torch
import recover_reference_predictions as r

def main():
    protocol=r.read(r.ROOT/'configs/legacy_default_replay_protocol.json')
    assert r.sha(__file__)==protocol['script_sha256']
    legacy=r.ROOT/protocol['legacy_source']
    assert r.sha(legacy)==protocol['legacy_sha256']
    original=r.read(r.ROOT/'configs/paired_family_analysis_protocol.json')
    for entry in original['frozen_files']:assert r.sha(r.ROOT/entry['path'])==entry['sha256']
    spec=importlib.util.spec_from_file_location('historical_encoders',legacy)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    from biocandidate.model.encoders import ESM2ProteinEncoder
    ESM2ProteinEncoder.forward=module.ESM2ProteinEncoder.forward
    manifest=r.read(r.ROOT/'configs/unikp_source_manifest.json')
    assert r.sha(r.SOURCE)==manifest['identity']['sha256']
    records=r.read_unikp_json(r.SOURCE).records
    records=r.apply_split_manifest(records,r.ROOT/'artifacts/homology-final/homology_split.json',source_identity=manifest['identity'])
    records=[x for x in r.aggregate_pair_measurements(records) if x.split=='test']
    with (r.ROOT/'artifacts/external/sota-homology-cold/data/test.csv').open(encoding='utf-8',newline='') as f:reference=list(csv.DictReader(f))
    assert len(reference)==len(records)
    for x,row in zip(records,reference):assert x.sequence==row['sequence'] and x.substrate_smiles==row['smiles'] and abs(x.log10_kcat-float(row['log10_kcat']))<1e-10
    out=r.ROOT/'artifacts/paired-family-analysis/legacy-default'
    out.mkdir(exist_ok=True)
    device=torch.device('cuda')
    for seed in protocol['seeds']:
        target=out/f'reference-seed{seed}.csv'
        if target.exists():raise FileExistsError(target)
        checkpoint=r.ROOT/f'artifacts/esm2-t6-opt-seed{seed}/best.pt'
        with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):model,_=r.load_checkpoint(checkpoint,device)
        p,_,y=r._collect_task_outputs(model,records,'log10_kcat',64,device)
        metrics=r.regression_metrics(p,y)
        prior=r.read(r.ROOT/f'artifacts/esm2-t6-opt-seed{seed}/test_metrics.json')['model_metrics']['log10_kcat']
        delta={k:metrics[k]-prior[k] for k in ['rmse','mae','pearson']}
        passed=all(abs(v)<=1e-5 for v in delta.values())
        receipt={'metrics':metrics,'historical_metrics':prior,'delta':delta,'passed':passed,'torch':torch.__version__,'batch_size':64,'encoder_forward_git':'57f3f36','checkpoint_sha256':r.sha(checkpoint),'protocol_sha256':r.sha(r.ROOT/'configs/legacy_default_replay_protocol.json')}
        print(seed,json.dumps(receipt),flush=True)
        if passed:
            with target.open('w',encoding='utf-8',newline='') as f:
                w=csv.writer(f);w.writerow(['pair_id','prediction','label']);w.writerows((row['pair_id'],float(a),float(b)) for row,a,b in zip(reference,p,y))
            receipt['predictions_sha256']=r.sha(target)
        (out/f'reference-seed{seed}-receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
        if not passed:break
        del model;torch.cuda.empty_cache()

if __name__=='__main__':main()
