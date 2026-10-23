"""Apply the unchanged bootstrap estimator to authenticated historical replay.

Reads configs/recovered_historical_analysis_protocol.json and re-verifies the
SHA256 of this runner, the paired_family_analysis estimator, and every frozen
input listed there before delegating to paired_family_analysis.main(), which
writes the artifacts/paired-family-analysis/legacy-default outputs. The three
per-seed reference receipts must already exist, pass, and match the legacy
replay protocol hash. Any failed assertion aborts before the estimator runs;
this runner is itself hash-bound by the protocol's runner_sha256.
"""
from pathlib import Path
import json
import paired_family_analysis as analysis

def main():
    """Re-verify every frozen identity, then rerun the estimator unchanged."""
    plan_path=analysis.ROOT/'configs/recovered_historical_analysis_protocol.json'
    plan=json.loads(plan_path.read_text())
    assert analysis.sha(__file__)==plan['runner_sha256']
    assert analysis.sha(Path(analysis.__file__))==plan['estimator_sha256']
    for entry in plan['frozen_files']:
        assert analysis.sha(analysis.ROOT/entry['path'])==entry['sha256']
    analysis.OUT=analysis.OUT/'legacy-default'
    for seed in [7,42,123]:
        receipt=json.loads((analysis.OUT/f'reference-seed{seed}-receipt.json').read_text())
        assert receipt['passed'] and all(abs(v)<=1e-5 for v in receipt['delta'].values())
        assert receipt['protocol_sha256']==analysis.sha(analysis.ROOT/'configs/legacy_default_replay_protocol.json')
    analysis.main()

if __name__=='__main__':main()
