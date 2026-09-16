"""DLKcat Mode B: retrain the published DLKcat architecture on the frozen split.

Imports the model classes verbatim from the pinned upstream run_model.py and injects
the published hyperparameters as module globals; the upstream 80/10/10 random split is
replaced by the frozen homology-cold partitions. Per seed: torch/np seeds set, 50
upstream iterations with lr decay every 10, and exactly one test evaluation at the end
through the shared evaluator. Governed by
configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import pickle
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "artifacts" / "external" / "sota-homology-cold" / "upstream-dlkcat"
LOG10_OF_2 = math.log10(2)

# Published setting string:
# all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3
# --lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50
PUBLISHED = {
    "dim": 20, "layer_gnn": 3, "window": 11, "layer_cnn": 3, "layer_output": 3,
    "lr": 1e-3, "lr_decay": 0.5, "decay_interval": 10, "weight_decay": 1e-6,
    "iteration": 50,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_upstream_run_model() -> object:
    spec = importlib.util.spec_from_file_location("dlkcat_upstream_run_model", UPSTREAM / "run_model.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_partition(input_dir: Path, split: str, device: torch.device) -> list[tuple]:
    def load_tensor(name: str, dtype) -> list:
        return [dtype(d).to(device)
                for d in np.load(input_dir / split / f"{name}.npy", allow_pickle=True)]

    compounds = load_tensor("compounds", torch.LongTensor)
    adjacencies = load_tensor("adjacencies", torch.FloatTensor)
    proteins = load_tensor("proteins", torch.LongTensor)
    interactions = load_tensor("regression", torch.FloatTensor)
    return list(zip(compounds, adjacencies, proteins, interactions))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/dlkcat-input")
    parser.add_argument("--data-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/data")
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/dlkcat-mode-b")
    parser.add_argument("--seeds", type=int, nargs="+", default=[7, 42, 123])
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    upstream = load_upstream_run_model()
    device = torch.device(args.device)

    with (args.input_dir / "fingerprint_dict.pickle").open("rb") as handle:
        fingerprint_dict = pickle.load(handle)
    with (args.input_dir / "sequence_dict.pickle").open("rb") as handle:
        word_dict = pickle.load(handle)

    # Upstream classes read these module globals at construction and forward time.
    upstream.n_fingerprint = len(fingerprint_dict)
    upstream.n_word = len(word_dict)
    for name, value in PUBLISHED.items():
        setattr(upstream, name, value)
    upstream.device = device

    dataset_train = load_partition(args.input_dir, "train", device)
    dataset_dev = load_partition(args.input_dir, "validation", device)
    dataset_test = load_partition(args.input_dir, "test", device)
    print(f"train {len(dataset_train)}, validation {len(dataset_dev)}, test {len(dataset_test)}")

    with (args.data_dir / "test.csv").open(newline="", encoding="utf-8") as handle:
        test_rows = list(csv.DictReader(handle))
    if len(test_rows) != len(dataset_test):
        raise ValueError("test CSV row count does not match encoded test partition")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "predictor": "DLKcat Mode B (published architecture retrained on frozen homology-cold train split)",
        "protocol_id": "sota-homology-cold-comparison-v1",
        "hyperparameters": PUBLISHED,
        "label_encoding": "log2(kcat) internally; converted to log10(kcat) for metrics",
        "identity": {
            "run_model.py": sha256(UPSTREAM / "run_model.py"),
            "input-manifest.json": sha256(args.input_dir / "input-manifest.json"),
            "data-manifest.json": sha256(ROOT / "artifacts/external/sota-homology-cold/data-manifest.json"),
        },
        "seeds": {},
    }

    for seed in args.seeds:
        run_dir = args.out_dir / f"seed{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)
        torch.manual_seed(seed)
        np.random.seed(seed)
        model = upstream.KcatPrediction().to(device)
        trainer = upstream.Trainer(model)
        tester = upstream.Tester(model)
        log_path = run_dir / "training-log.txt"
        started = time.time()
        with log_path.open("w", encoding="utf-8") as log:
            log.write("Epoch\tTime(sec)\tLoss_train\tMAE_dev\tRMSE_dev\tR2_dev\n")
            for epoch in range(1, PUBLISHED["iteration"] + 1):
                if epoch % PUBLISHED["decay_interval"] == 0:
                    trainer.optimizer.param_groups[0]["lr"] *= PUBLISHED["lr_decay"]
                loss_train, _, _ = trainer.train(dataset_train)
                mae_dev, rmse_dev, r2_dev = tester.test(dataset_dev)
                log.write(
                    f"{epoch}\t{time.time() - started:.0f}\t{loss_train:.4f}\t"
                    f"{mae_dev:.4f}\t{rmse_dev:.4f}\t{r2_dev:.4f}\n"
                )
                log.flush()
                if epoch % 10 == 0:
                    print(f"seed {seed} epoch {epoch}: loss {loss_train:.4f} "
                          f"dev RMSE {rmse_dev:.4f}", flush=True)
        # One test evaluation per seed, on the final model, on the log10 scale.
        predicted_log2 = []
        with torch.no_grad():
            for data in dataset_test:
                _, predicted = model(data, train=False)
                predicted_log2.append(float(predicted[0]))
        predicted_log10 = np.asarray(predicted_log2) * LOG10_OF_2
        observed_log10 = np.asarray([float(row["log10_kcat"]) for row in test_rows])
        from biocandidate.evaluation import regression_metrics
        metrics = regression_metrics(
            torch.tensor(predicted_log10, dtype=torch.float64),
            torch.tensor(observed_log10, dtype=torch.float64),
        )
        torch.save(model.state_dict(), run_dir / "model.pt")
        with (run_dir / "predictions.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["pair_id", "log10_kcat_observed", "log10_kcat_predicted"])
            for row, predicted in zip(test_rows, predicted_log10):
                writer.writerow([row["pair_id"], row["log10_kcat"], repr(float(predicted))])
        (run_dir / "test_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
        summary["seeds"][str(seed)] = metrics
        print(f"seed {seed}: {metrics}", flush=True)

    seeds = args.seeds
    summary["three_seed_mean"] = {
        key: float(np.mean([summary["seeds"][str(s)][key] for s in seeds]))
        for key in ("rmse", "mae", "pearson")
    }
    summary["three_seed_sd"] = {
        key: float(np.std([summary["seeds"][str(s)][key] for s in seeds], ddof=1))
        for key in ("rmse", "mae", "pearson")
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"mean": summary["three_seed_mean"], "sd": summary["three_seed_sd"]}, indent=2))


if __name__ == "__main__":
    main()
