#!/usr/bin/env bash
set -euo pipefail
cd /mnt/d/biological/BioCandidateRanker
export PYTHONPYCACHEPREFIX="$PWD/artifacts/catpred-mode-b/pycache"
for catpred_seed in 7 42 123; do
  artifacts/catpred-mode-b/venv/bin/python -u scripts/run_catpred_full_v3_attempt2.py --seed "$catpred_seed"
done
