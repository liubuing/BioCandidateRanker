#!/usr/bin/env bash
set -euo pipefail
cd /mnt/d/biological/BioCandidateRanker
if [[ -e artifacts/catpred-mode-b/length-pilot-run ]]; then exit 1; fi
export CATPRED_CACHE_PATH="$PWD/artifacts/catpred-mode-b/cache"
export TORCH_HOME=/mnt/d/DisorderFlowRuntime/cache/torch
export PYTHONPYCACHEPREFIX="$PWD/artifacts/catpred-mode-b/pycache"
export MPLCONFIGDIR="$PWD/artifacts/catpred-mode-b/matplotlib"
export CATPRED_ESM_BATCH_SIZE=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
artifacts/catpred-mode-b/venv/bin/python -u artifacts/catpred-mode-b/upstream/train.py \
 --protein_records_path artifacts/catpred-mode-b/length-pilot-data/protein_records.json.gz \
 --data_path artifacts/catpred-mode-b/length-pilot-data/train.csv \
 --separate_val_path artifacts/catpred-mode-b/length-pilot-data/validation.csv \
 --separate_test_path artifacts/catpred-mode-b/length-pilot-data/validation.csv \
 --dataset_type regression --smiles_columns reactant_smiles --target_columns log10kcat_max \
 --extra_metrics mae mse r2 --ensemble_size 1 --seq_embed_dim 36 --seq_self_attn_nheads 6 \
 --loss_function mve --batch_size 16 --epochs 3 --add_esm_feats --seed 7 --pytorch_seed 7 \
 --save_dir artifacts/catpred-mode-b/length-pilot-run
