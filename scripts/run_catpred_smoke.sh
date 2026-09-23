#!/usr/bin/env bash
set -euo pipefail
cd /mnt/d/biological/BioCandidateRanker
if [[ -e artifacts/catpred-mode-b/smoke-run/fold_0/model_0/model.pt ]]; then
 echo "Existing smoke checkpoint preserved; choose a new reviewed output directory." >&2
 exit 1
fi
export CATPRED_CACHE_PATH="$PWD/artifacts/catpred-mode-b/cache"
export TORCH_HOME=/mnt/d/DisorderFlowRuntime/cache/torch
export PYTHONPYCACHEPREFIX="$PWD/artifacts/catpred-mode-b/pycache"
export MPLCONFIGDIR="$PWD/artifacts/catpred-mode-b/matplotlib"
export CATPRED_ESM_BATCH_SIZE=2
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
python_bin="$PWD/artifacts/catpred-mode-b/venv/bin/python"
"$python_bin" -m pip freeze > artifacts/catpred-mode-b/environment-freeze.txt
"$python_bin" -u artifacts/catpred-mode-b/upstream/train.py \
 --protein_records_path artifacts/catpred-mode-b/smoke-data/protein_records.json.gz \
 --data_path artifacts/catpred-mode-b/smoke-data/train.csv \
 --separate_val_path artifacts/catpred-mode-b/smoke-data/validation.csv \
 --separate_test_path artifacts/catpred-mode-b/smoke-data/validation.csv \
 --dataset_type regression --smiles_columns reactant_smiles --target_columns log10kcat_max \
 --extra_metrics mae mse r2 --ensemble_size 1 --seq_embed_dim 36 --seq_self_attn_nheads 6 \
 --loss_function mve --batch_size 16 --epochs 3 --add_esm_feats --seed 7 --pytorch_seed 7 \
 --save_dir artifacts/catpred-mode-b/smoke-run
