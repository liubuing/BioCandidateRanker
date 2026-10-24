"""Audit both selectors against the complete common validation history.

Reads the completed-training manifest and each seed's full epoch log under
``artifacts/selection-rule-comparison``, re-derives the epoch each rule (NLL
and RMSE) would select, confirms the saved checkpoints match those epochs and
their recorded hashes, and writes ``selection-summary.json``. Governance role:
validation-only diagnostic -- explicitly biased because selection uses the same
validation history it is measured on, no test inference occurs, and a
pre-existing summary aborts the run.

中文：对照完整共同验证历史审计两种选择器：读完成训练清单与逐种子全量
epoch 日志，重推 NLL 与 RMSE 各会选的 epoch，确认后写选择汇总。
"""
from pathlib import Path
import json,hashlib
import torch
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/selection-rule-comparison'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    """Re-derive both selection rules per seed and write the selection summary.
中文：逐种子重推两种选择规则并写选择汇总。
"""
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
