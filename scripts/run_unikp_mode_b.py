"""UniKP Mode B: retrain UniKP's ExtraTrees regressor on the frozen homology-cold split.

Uses the local UniKP deployment's frozen pretrained feature extractors exactly as
published (SMILES transformer + ProtT5-XL-UniRef50 mean pooling), computes embeddings
once over unique rows, and fits one ExtraTreesRegressor per frozen seed on the train
partition only. Test predictions and metrics are written per seed under the protocol's
artifacts root. Governed by configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pickle
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biocandidate.evaluation import regression_metrics  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_partition(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unikp-root", type=Path, required=True,
                        help="local UniKP deployment directory (code/, models/)")
    parser.add_argument("--data-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/data")
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/unikp-mode-b")
    parser.add_argument("--seeds", type=int, nargs="+", default=[7, 42, 123])
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--prot-t5-batch", type=int, default=8)
    args = parser.parse_args()

    code_dir = args.unikp_root / "code"
    models_dir = args.unikp_root / "models"
    sys.path.insert(0, str(code_dir))
    from build_vocab import WordVocab as _WordVocab  # noqa: E402
    from pretrain_trfm import TrfmSeq2seq  # noqa: E402
    from utils import split as split_smiles  # noqa: E402

    # vocab.pkl was written by build_vocab.py executed as __main__, so its pickle
    # resolves WordVocab against the loading process's __main__ module.
    import __main__  # noqa: E402

    __main__.WordVocab = _WordVocab
    WordVocab = _WordVocab

    partitions = {
        name: load_partition(args.data_dir / f"{name}.csv")
        for name in ("train", "validation", "test")
    }
    for name, rows in partitions.items():
        for row in rows:
            if float(row["log10_kcat"]) != float(row["log10_kcat"]):
                raise ValueError(f"non-finite label in {name}")

    unique_smiles = sorted({row["smiles"] for rows in partitions.values() for row in rows})
    unique_sequences = sorted({
        row["sequence"] for rows in partitions.values() for row in rows
    })
    print(f"unique SMILES: {len(unique_smiles)}, unique sequences: {len(unique_sequences)}")

    cache_dir = args.out_dir / "embedding-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    smiles_cache = cache_dir / "smiles-trfm-embeddings.npy"
    smiles_index = cache_dir / "smiles-index.json"
    if smiles_cache.is_file() and smiles_index.is_file():
        smiles_vec = np.load(smiles_cache)
        smiles_keys = json.loads(smiles_index.read_text(encoding="utf-8"))
        if smiles_keys != unique_smiles:
            raise ValueError("cached SMILES index does not match the partition data")
        print("reusing cached SMILES embeddings")
    else:
        vocab = WordVocab.load_vocab(str(models_dir / "vocab.pkl"))
        trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
        trfm.load_state_dict(torch.load(models_dir / "trfm_12_23000.pkl",
                                        map_location="cpu", weights_only=True))
        trfm.eval()

        pad_index, unk_index, sos_index, eos_index = 0, 1, 2, 3

        def get_inputs(sm: str) -> tuple[list[int], list[int]]:
            seq_len = 220
            sm = sm.split()
            if len(sm) > 218:
                sm = sm[:109] + sm[-109:]
            ids = [vocab.stoi.get(token, unk_index) for token in sm]
            ids = [sos_index] + ids + [eos_index]
            seg = [1] * len(ids)
            padding = [pad_index] * (seq_len - len(ids))
            ids.extend(padding)
            seg.extend(padding)
            return ids, seg

        tokenized = [split_smiles(sm) for sm in unique_smiles]
        x_id, _ = [], []
        for sm in tokenized:
            ids, seg = get_inputs(sm)
            x_id.append(ids)
            _.append(seg)
        x_id = torch.tensor(x_id)
        vectors = trfm.encode(torch.t(x_id))
        if vectors.shape != (len(unique_smiles), 1024):
            raise ValueError(f"unexpected SMILES feature shape {vectors.shape}")
        smiles_vec = vectors
        np.save(smiles_cache, smiles_vec)
        smiles_index.write_text(json.dumps(unique_smiles), encoding="utf-8")
        print("SMILES embeddings written")

    smiles_lookup = {smiles: i for i, smiles in enumerate(unique_smiles)}

    prot_cache = cache_dir / "prott5-embeddings.npy"
    prot_index = cache_dir / "prott5-index.json"
    if prot_cache.is_file() and prot_index.is_file():
        prot_vec = np.load(prot_cache)
        prot_keys = json.loads(prot_index.read_text(encoding="utf-8"))
        if prot_keys != unique_sequences:
            raise ValueError("cached sequence index does not match the partition data")
        print("reusing cached ProtT5 embeddings")
    else:
        from transformers import AutoTokenizer, T5EncoderModel

        model_path = models_dir / "prot_t5_xl_uniref50"
        tokenizer = AutoTokenizer.from_pretrained(str(model_path), do_lower_case=False)
        encoder = T5EncoderModel.from_pretrained(str(model_path)).to(args.device)
        encoder.eval()
        features: list[np.ndarray] = []
        started = time.time()
        with torch.no_grad():
            for start in range(0, len(unique_sequences), args.prot_t5_batch):
                batch = unique_sequences[start:start + args.prot_t5_batch]
                prepared = []
                for sequence in batch:
                    if len(sequence) > 1000:
                        sequence = sequence[:500] + sequence[-500:]
                    prepared.append(re.sub(r"[UZOB]", "X", " ".join(sequence)))
                encoded = tokenizer(prepared, add_special_tokens=True, padding=True)
                input_ids = torch.tensor(encoded["input_ids"]).to(args.device)
                attention_mask = torch.tensor(encoded["attention_mask"]).to(args.device)
                embedding = encoder(input_ids=input_ids,
                                    attention_mask=attention_mask).last_hidden_state.cpu().numpy()
                for row_index in range(len(embedding)):
                    seq_len = int(attention_mask[row_index].sum())
                    seq_emb = embedding[row_index][:seq_len - 1]
                    features.append(seq_emb.mean(axis=0))
                done = start + len(batch)
                if done % 160 < args.prot_t5_batch:
                    rate = done / max(time.time() - started, 1e-9)
                    print(f"ProtT5 {done}/{len(unique_sequences)} ({rate:.1f}/s)", flush=True)
        prot_vec = np.asarray(features, dtype=float)
        if prot_vec.shape != (len(unique_sequences), 1024):
            raise ValueError(f"unexpected sequence feature shape {prot_vec.shape}")
        np.save(prot_cache, prot_vec)
        prot_index.write_text(json.dumps(unique_sequences), encoding="utf-8")
        print("ProtT5 embeddings written")

    prot_lookup = {sequence: i for i, sequence in enumerate(unique_sequences)}

    def design(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray]:
        indices_s = [smiles_lookup[row["smiles"]] for row in rows]
        indices_p = [prot_lookup[row["sequence"]] for row in rows]
        x = np.concatenate([smiles_vec[indices_s], prot_vec[indices_p]], axis=1)
        y = np.asarray([float(row["log10_kcat"]) for row in rows])
        return x, y

    x_train, y_train = design(partitions["train"])
    x_test, y_test = design(partitions["test"])
    print(f"train design {x_train.shape}, test design {x_test.shape}")

    from sklearn.ensemble import ExtraTreesRegressor

    summary = {
        "predictor": "UniKP Mode B (ExtraTrees retrained on frozen homology-cold train split)",
        "protocol_id": "sota-homology-cold-comparison-v1",
        "feature_dimensions": int(x_train.shape[1]),
        "train_rows": int(x_train.shape[0]),
        "test_rows": int(x_test.shape[0]),
        "identity": {
            "trfm_12_23000.pkl": sha256(models_dir / "trfm_12_23000.pkl"),
            "vocab.pkl": sha256(models_dir / "vocab.pkl"),
            "prot_t5_pytorch_model.bin": sha256(models_dir / "prot_t5_xl_uniref50" / "pytorch_model.bin"),
            "data-manifest.json": sha256(ROOT / "artifacts/external/sota-homology-cold/data-manifest.json"),
        },
        "seeds": {},
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for seed in args.seeds:
        run_dir = args.out_dir / f"seed{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)
        model = ExtraTreesRegressor(random_state=seed)
        model.fit(x_train, y_train)
        predictions = model.predict(x_test)
        metrics = regression_metrics(
            torch.tensor(predictions, dtype=torch.float64),
            torch.tensor(y_test, dtype=torch.float64),
        )
        with (run_dir / "predictions.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["pair_id", "log10_kcat_observed", "log10_kcat_predicted"])
            for row, predicted in zip(partitions["test"], predictions):
                writer.writerow([row["pair_id"], row["log10_kcat"], repr(float(predicted))])
        with (run_dir / "model.pkl").open("wb") as handle:
            pickle.dump(model, handle)
        (run_dir / "test_metrics.json").write_text(
            json.dumps(metrics, indent=2) + "\n", encoding="utf-8"
        )
        summary["seeds"][str(seed)] = metrics
        print(f"seed {seed}: {metrics}")

    mean = {
        key: float(np.mean([summary["seeds"][str(seed)][key] for seed in args.seeds]))
        for key in ("rmse", "mae", "pearson")
    }
    sd = {
        key: float(np.std([summary["seeds"][str(seed)][key] for seed in args.seeds], ddof=1))
        for key in ("rmse", "mae", "pearson")
    }
    summary["three_seed_mean"] = mean
    summary["three_seed_sd"] = sd
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"mean": mean, "sd": sd}, indent=2))


if __name__ == "__main__":
    main()
