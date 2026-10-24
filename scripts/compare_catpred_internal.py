"""Exploratory paired-family comparison after frozen CatPred inference.

Reads the frozen internal test split, the homology family mapping, and the
per-seed prediction CSVs for CatPred, CataPro, UniKP and both
BioCandidateRanker replays, then computes family-stratified paired bootstrap
contrasts (seed 20260918, 5000 replicates) against the adapted CatPred run.
Writes artifacts/catpred-mode-b/internal-test/paired-comparison.json,
refusing to overwrite an existing report. Exploratory only: no multiplicity
adjustment and a CI crossing zero never establishes equivalence.

中文：本脚本在冻结的 CatPred 推理之后做探索性的配对家族比较：读取冻结内部测试划分、
同源家族映射，以及 CatPred、CataPro、UniKP 与两个 BioCandidateRanker 回放的逐种子预测
CSV，再针对适配后的 CatPred 运行计算家族分层的配对 bootstrap 对比（种子 20260918、
5000 次重抽样）。写出 artifacts/catpred-mode-b/internal-test/paired-comparison.json，
拒绝覆盖既有报告。治理角色：探索性对比收据——无多重性校正，跨零置信区间绝不确立
等效性，不构成确认性结论。
"""
from pathlib import Path
import json, hashlib
import numpy as np
from paired_family_analysis import rows, aligned, metrics_by_draw, sha

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/catpred-mode-b/internal-test'

def main():
    """Compute and write the paired family-bootstrap comparison report.

    中文：计算并写出配对家族 bootstrap 对比报告。
    """
    target=OUT/'paired-comparison.json'
    if target.exists(): raise FileExistsError(target)
    data=rows(ROOT/'artifacts/external/sota-homology-cold/data/test.csv')
    ids=[r['pair_id'] for r in data];labels=np.array([float(r['log10_kcat']) for r in data])
    split=json.loads((ROOT/'artifacts/homology-final/homology_split.json').read_text())
    mapping={r['sequence_id']:r['cluster'] for r in split['rows'] if r['split']=='test'}
    families=[mapping['seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]] for r in data]
    unique,groups=np.unique(families,return_inverse=True);g=len(unique)
    weights=np.vstack([np.ones(g),np.random.default_rng(20260918).multinomial(g,np.full(g,1/g),size=5000)])
    specs={
      'CatPred adapted':('artifacts/catpred-mode-b/internal-test/seed{seed}-predictions.csv','prediction_log10_kcat','label_log10_kcat'),
      'CataPro':('artifacts/catapro-mode-b/internal-test/seed{seed}-predictions.csv','prediction_log10_kcat','label_log10_kcat'),
      'UniKP':('artifacts/external/sota-homology-cold/unikp-mode-b/seed{seed}/predictions.csv','log10_kcat_predicted','log10_kcat_observed'),
      'BioCandidateRanker historical':('artifacts/paired-family-analysis/legacy-default/reference-seed{seed}.csv','prediction','label'),
      'BioCandidateRanker current-code replay':('artifacts/paired-family-analysis/current-code/reference-seed{seed}.csv','prediction','label')}
    metrics={};summaries={};inputs=[];names=['rmse','mae','pearson']
    for name,(pattern,pk,lk) in specs.items():
        ps=[]
        for seed in [7,42,123]:
            path=ROOT/pattern.format(seed=seed);ps.append(aligned(path,ids,labels,pk,lk));inputs.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path)})
        metrics[name]=metrics_by_draw(ps,labels,groups,weights)
        assert np.isfinite(metrics[name]).all()
        per=np.array([metrics_by_draw([v],labels,groups,np.ones((1,g)))[0] for v in ps])
        summaries[name]={k:{'mean':float(per[:,i].mean()),'sd':float(per[:,i].std(ddof=1))} for i,k in enumerate(names)}
    contrasts={}
    for name in specs:
        if name=='CatPred adapted':continue
        delta=metrics[name]-metrics['CatPred adapted']
        contrasts[name+' minus CatPred']={k:{'difference':float(delta[0,i]),'ci95':np.quantile(delta[1:,i],[.025,.975]).tolist()} for i,k in enumerate(names)}
    report={'rows':len(data),'families':g,'bootstrap_replicates':5000,'bootstrap_seed':20260918,'mean_sd':summaries,'contrasts':contrasts,'inputs':inputs,'scope':'Exploratory internal comparison; fixed seeds, family resampling, no multiplicity adjustment. Current-code replay is not retraining; CI crossing zero does not establish equivalence.'}
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))

if __name__=='__main__':main()
