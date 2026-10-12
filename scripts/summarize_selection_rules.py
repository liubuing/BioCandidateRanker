"""Audit both selectors against the complete common validation history."""
from pathlib import Path
import json,hashlib
import torch
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/selection-rule-comparison'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    out=BASE/'selection-summary.json'
    if out.exists():raise FileExistsError(out)
    manifest=json.loads((BASE/'completed-training-manifest.json').read_text())
    assert len(manifest['runs'])==3
    result={'scope':'Validation-selection diagnostic only; biased by selection on the same validation set; no test inference','runs':[]}
    for run in manifest['runs']:
        seed=run['seed'];folder=BASE/f'seed{seed}'
        rows=[json.loads(x) for x in (folder/'train.log').read_text().splitlines() if x.startswith('{')]
        rows=[x for x in rows if 'epoch' in x];assert [x['epoch'] for x in rows]==list(range(1,11))
        chosen={'nll':min(rows,key=lambda x:x['validation_loss']),'rmse':min(rows,key=lambda x:x['validation_tasks']['log10_kcat']['rmse'])}
        entry={'seed':seed,'rules':{}}
        for rule,filename,key in [('nll','best.pt','checkpoint_sha256'),('rmse','best_rmse.pt','rmse_checkpoint_sha256')]:
            path=folder/filename;assert sha(path)==run[key]
            with torch.serialization.safe_globals([torch.torch_version.TorchVersion]):payload=torch.load(path,map_location='cpu',weights_only=True)
            assert payload['epoch']==chosen[rule]['epoch']
            entry['rules'][rule]={'epoch':payload['epoch'],'validation_loss':chosen[rule]['validation_loss'],'metrics':chosen[rule]['validation_tasks']['log10_kcat'],'checkpoint_sha256':sha(path)}
        result['runs'].append(entry)
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
