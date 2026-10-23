"""Compare cached GPU molecule features against pinned upstream CPU extraction.

Reads the first 32 training rows via ``run_catapro_mode_b`` helpers (train
partition only, never test), recomputes MolT5 embeddings and MACCS keys with
the pinned upstream CPU code for three fixed positions, and compares them to
the cached GPU-extracted features. Writes ``molecule-equivalence-check.json``
with the max absolute error and a pass/fail verdict at the recorded 1e-4
tolerance; a failure raises and leaves the recorded verdict in place.
Governance role: cache-equivalence gate for the CataPro Mode B pipeline.
"""
import json
import numpy as np
import run_catapro_mode_b as runner

upstream=runner.source_module('inference/utils.py')
rows=runner.load_partitions(32)['train']
indices=[0,len(rows)//2,len(rows)-1]
smiles=[rows[i]['smiles'] for i in indices]
reference=upstream.get_molT5_embed(smiles,str(runner.MOLT5))
maccs=upstream.GetMACCSKeys(smiles)
data,manifest=runner.load_feature_data('smoke')
actual=data['train'][0][indices,1024:]
error=float(np.max(np.abs(reference-actual[:,:768])))
result={'upstream_commit':runner.read_json(runner.PROTOCOL)['upstream_commit'],
        'feature_identity':manifest['identity'],
        'inputs':'3 fixed positions in first 32 training rows; no test inputs',
        'molt5_max_absolute_error':error,'maccs_exact':bool(np.array_equal(maccs,actual[:,768:])),
        'tolerance':1e-4,'passed':bool(error<=1e-4 and np.array_equal(maccs,actual[:,768:]))}
runner.save_json(runner.OUT/'molecule-equivalence-check.json',result)
print(json.dumps(result,indent=2))
if not result['passed']:raise ValueError('Upstream molecule feature comparison failed')
