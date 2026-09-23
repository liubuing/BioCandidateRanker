"""Compare cached GPU molecule features against pinned upstream CPU extraction."""
import json
from pathlib import Path
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
