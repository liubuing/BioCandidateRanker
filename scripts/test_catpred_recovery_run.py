"""Actual upstream CPU training: uninterrupted vs epoch-0 interruption/restart.

Runs the pinned upstream CatPred entry point three times on the smoke data
(continuous, interrupted after epoch 0, and resumed) under
``artifacts/catpred-mode-b/recovery-test-v3``, then compares the recovery
states bitwise (model tensors, sampler RNG, torch RNG, iteration, best score).
Writes ``result.json`` on success. Governance role: integration test of
crash-resume determinism only -- CPU scope, with CUDA nondeterminism explicitly
excluded from the claim.

中文：真实上游 CPU 训练：连续运行 vs epoch-0 中断重启的对比，断言恢复态
逐位一致。
"""
from pathlib import Path
import os,subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'artifacts/catpred-mode-b'

def main():
    """Run the interrupt/resume scenarios and assert bitwise-equal recovery states.
中文：跑中断/恢复场景并断言逐位一致的恢复状态。
"""
    import torch
    out=BASE/'recovery-test-v3';out.mkdir(exist_ok=True)
    env=os.environ.copy();env.update(CATPRED_RECOVERY_ID='recovery-test-v1',CATPRED_CACHE_PATH=str(BASE/'cache'),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    common=[sys.executable,'-u',str(ROOT/'scripts/catpred_recovery_entry.py'),'--protein_records_path',str(BASE/'smoke-data/protein_records.json.gz'),'--data_path',str(BASE/'smoke-data/train.csv'),'--separate_val_path',str(BASE/'smoke-data/validation.csv'),'--separate_test_path',str(BASE/'smoke-data/validation.csv'),'--dataset_type','regression','--smiles_columns','reactant_smiles','--target_columns','log10kcat_max','--ensemble_size','1','--seq_embed_dim','36','--seq_self_attn_nheads','6','--loss_function','mve','--batch_size','16','--epochs','3','--add_esm_feats','--seed','7','--pytorch_seed','7','--no_cuda']
    for name,directory,stop,expected in [('continuous','continuous',False,0),('interrupted','resumed',True,75),('resumed','resumed',False,0)]:
        runenv=env.copy()
        if stop:runenv['CATPRED_TEST_STOP_EPOCH']='0'
        with (out/f'{name}.log').open('w') as log:r=subprocess.run(common+['--save_dir',str(out/directory)],cwd=ROOT,env=runenv,stdout=log,stderr=subprocess.STDOUT)
        print(name,r.returncode,flush=True)
        if r.returncode!=expected:raise RuntimeError(name)
    states=[torch.load(out/d/'fold_0/model_0/recovery.pt',map_location='cpu',weights_only=True) for d in ['continuous','resumed']]
    assert states[0]['epoch']==states[1]['epoch']==2
    assert states[0]['n_iter']==states[1]['n_iter']
    assert states[0]['sampler_rng']==states[1]['sampler_rng']
    assert all(torch.equal(states[0]['model'][k],states[1]['model'][k]) for k in states[0]['model'])
    assert states[0]['best_score']==states[1]['best_score']
    assert torch.equal(states[0]['torch_rng'],states[1]['torch_rng'])
    receipt={'status':'passed','epochs':3,'interruption_after_epoch':0,'final_model_bitwise_equal':True,'sampler_rng_equal':True,'torch_rng_equal':True,'best_score_equal':True,'scope':'CPU actual upstream integration; CUDA nondeterministic reductions may still differ slightly.'}
    (out/'result.json').write_text(json.dumps(receipt,indent=2));print(receipt,flush=True)

if __name__=='__main__':main()
