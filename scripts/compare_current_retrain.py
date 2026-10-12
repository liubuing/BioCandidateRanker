"""Exploratory family-bootstrap comparisons of the retrained current model."""
from pathlib import Path
import json,hashlib
import numpy as np
from paired_family_analysis import rows,aligned,metrics_by_draw,sha
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/current-code-retrain/internal-test'

def main():
    target=OUT/'paired-comparison.json'
    if target.exists():raise FileExistsError(target)
    data=rows(ROOT/'artifacts/external/sota-homology-cold/data/test.csv');ids=[r['pair_id'] for r in data];y=np.array([float(r['log10_kcat']) for r in data])
    split=json.loads((ROOT/'artifacts/homology-final/homology_split.json').read_text());mapping={r['sequence_id']:r['cluster'] for r in split['rows'] if r['split']=='test'}
    families=[mapping['seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]] for r in data];unique,groups=np.unique(families,return_inverse=True);g=len(unique)
    weights=np.vstack([np.ones(g),np.random.default_rng(20260918).multinomial(g,np.full(g,1/g),size=5000)])
    specs={
      'current_retrained':('artifacts/current-code-retrain/internal-test/seed{seed}-predictions.csv','prediction','label'),
      'CatPred':('artifacts/catpred-mode-b/internal-test/seed{seed}-predictions.csv','prediction_log10_kcat','label_log10_kcat'),
      'CataPro':('artifacts/catapro-mode-b/internal-test/seed{seed}-predictions.csv','prediction_log10_kcat','label_log10_kcat'),
      'UniKP':('artifacts/external/sota-homology-cold/unikp-mode-b/seed{seed}/predictions.csv','log10_kcat_predicted','log10_kcat_observed'),
      'historical':('artifacts/paired-family-analysis/legacy-default/reference-seed{seed}.csv','prediction','label')}
    stats={};summary={};inputs=[];names=['rmse','mae','pearson']
    for name,(pattern,pk,lk) in specs.items():
        preds=[]
        for seed in [7,42,123]:
            f=ROOT/pattern.format(seed=seed);preds.append(aligned(f,ids,y,pk,lk));inputs.append({'path':f.relative_to(ROOT).as_posix(),'sha256':sha(f)})
        stats[name]=metrics_by_draw(preds,y,groups,weights);assert np.isfinite(stats[name]).all()
        per=np.array([metrics_by_draw([p],y,groups,np.ones((1,g)))[0] for p in preds])
        summary[name]={k:{'mean':float(per[:,i].mean()),'sd':float(per[:,i].std(ddof=1))} for i,k in enumerate(names)}
    contrasts={}
    for name in specs:
        if name=='current_retrained':continue
        d=stats['current_retrained']-stats[name];contrasts[name]={k:{'difference':float(d[0,i]),'ci95':np.quantile(d[1:,i],[.025,.975]).tolist()} for i,k in enumerate(names)}
    result={'rows':len(y),'families':g,'replicates':5000,'bootstrap_seed':20260918,'mean_sd':summary,'current_minus_comparator':contrasts,'inputs':inputs,'scope':'Exploratory internal diagnostic, already observed split; no multiplicity correction; zero-crossing CI is not equivalence'}
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='inputs'},indent=2))

if __name__=='__main__':main()
