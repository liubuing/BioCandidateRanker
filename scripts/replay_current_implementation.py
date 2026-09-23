"""Separate current-code audit; never replace historical predictions or metrics."""
import json
import sys
from pathlib import Path
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import recover_reference_predictions as replay

def main():
    protocol_path = replay.ROOT/'configs/current_implementation_replay_protocol.json'
    protocol = replay.read(protocol_path)
    assert replay.sha(__file__) == protocol['script_sha256']
    original = replay.read(replay.ROOT/'configs/paired_family_analysis_protocol.json')
    for entry in original['frozen_files']:
        assert replay.sha(replay.ROOT/entry['path']) == entry['sha256']
    for entry in protocol['source_files']:
        assert replay.sha(replay.ROOT/entry['path']) == entry['sha256']
    manifest = replay.read(replay.ROOT/'configs/unikp_source_manifest.json')
    assert replay.sha(replay.SOURCE) == manifest['identity']['sha256']
    records = replay.read_unikp_json(replay.SOURCE).records
    records = replay.apply_split_manifest(records,replay.ROOT/'artifacts/homology-final/homology_split.json',source_identity=manifest['identity'])
    records = [r for r in replay.aggregate_pair_measurements(records) if r.split=='test']
    import csv
    with (replay.ROOT/'artifacts/external/sota-homology-cold/data/test.csv').open(encoding='utf-8',newline='') as f:
        reference = list(csv.DictReader(f))
    assert len(reference)==len(records)==1646
    for r, row in zip(records, reference):
        assert r.sequence==row['sequence'] and r.substrate_smiles==row['smiles']
        assert abs(r.log10_kcat-float(row['log10_kcat']))<1e-10
    out = replay.ROOT/'artifacts/paired-family-analysis/current-code'
    out.mkdir(exist_ok=True)
    device = torch.device('cuda')
    for seed in protocol['seeds']:
        target = out/f'reference-seed{seed}.csv'
        if target.exists():
            raise FileExistsError(target)
        ckpt = replay.ROOT/f'artifacts/esm2-t6-opt-seed{seed}/best.pt'
        with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):
            model, payload = replay.load_checkpoint(ckpt,device)
        p,_,y = replay._collect_task_outputs(model,records,'log10_kcat',protocol['batch_size'],device)
        metrics = replay.regression_metrics(p,y)
        prior = replay.read(replay.ROOT/f'artifacts/esm2-t6-opt-seed{seed}/test_metrics.json')['model_metrics']['log10_kcat']
        with target.open('w',encoding='utf-8',newline='') as f:
            writer=csv.writer(f);writer.writerow(['pair_id','prediction','label'])
            writer.writerows((r['pair_id'],float(a),float(b)) for r,a,b in zip(reference,p,y))
        receipt={'metrics':metrics,'historical_metrics':prior,
            'difference_from_historical':{k:metrics[k]-prior[k] for k in ['rmse','mae','pearson']},
            'predictions_sha256':replay.sha(target),'checkpoint_sha256':replay.sha(ckpt),
            'protocol_sha256':replay.sha(protocol_path),'scope':protocol['scope']}
        (out/f'reference-seed{seed}-receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
        print(seed,receipt,flush=True)
        del model
        torch.cuda.empty_cache()

if __name__=='__main__':main()
