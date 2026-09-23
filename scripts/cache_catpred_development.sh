#!/usr/bin/env bash
set -euo pipefail
cd /mnt/d/biological/BioCandidateRanker
export CATPRED_CACHE_PATH="$PWD/artifacts/catpred-mode-b/cache"
export TORCH_HOME=/mnt/d/DisorderFlowRuntime/cache/torch
export PYTHONPYCACHEPREFIX="$PWD/artifacts/catpred-mode-b/pycache"
export MPLCONFIGDIR="$PWD/artifacts/catpred-mode-b/matplotlib"
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
artifacts/catpred-mode-b/venv/bin/python -u scripts/cache_catpred_development.py
