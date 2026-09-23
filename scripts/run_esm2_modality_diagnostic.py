"""Run the frozen post-hoc ESM-2 modality diagnostic sequentially on one GPU."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "configs" / "esm2_modality_diagnostic_protocol.json"
SOURCE = Path(
    r"C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment"
    r"\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json"
)
SPLIT = ROOT / "artifacts" / "homology-final" / "homology_split.json"
MANIFEST = ROOT / "artifacts" / "unikp_audit.json"
OUTPUT = ROOT / "artifacts" / "esm2-modality-diagnostic"


def identity(path: Path) -> dict[str, int | str]:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"size_bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def checked_identity(path: Path, expected: dict) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = identity(path)
    for key in ("size_bytes", "sha256"):
        if actual[key] != expected[key]:
            raise ValueError(f"{path}: {key} mismatch: {actual[key]}")


def call_logged(arguments: list[str], logfile: Path) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    with logfile.open("w", encoding="utf-8") as output:
        subprocess.run(
            [sys.executable, "-u", "-m", "biocandidate.cli", *arguments],
            cwd=ROOT,
            env=environment,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    checked_identity(SOURCE, protocol["source_identity"])
    checked_identity(SPLIT, protocol["split_identity"])
    training = protocol["training"]
    for variant in ("esm2_protein_only", "esm2_protein_molecule"):
        variant_flags = protocol["variants"][variant]
        for seed in protocol["seeds"]:
            name = variant.removeprefix("esm2_").replace("_", "-")
            directory = OUTPUT / f"{name}-seed{seed}"
            directory.mkdir(parents=True, exist_ok=True)
            metrics_path = directory / "test_metrics.json"
            if args.skip_existing and metrics_path.is_file():
                print(f"skip {variant} seed {seed}: evaluated", flush=True)
                continue
            checkpoint = directory / "best.pt"
            completion_marker = directory / "training-complete.json"
            if not completion_marker.is_file():
                if checkpoint.is_file():
                    raise RuntimeError(
                        f"{directory}: checkpoint exists without training-complete.json; "
                        "inspect the interrupted run before resuming")
                train_arguments = [
                    "train", "--data", str(SOURCE), "--manifest", str(MANIFEST),
                    "--split-manifest", str(SPLIT), "--conflict-policy", "median",
                    "--protein-encoder", "esm2", "--d-model", str(training["d_model"]),
                    "--num-heads", str(training["num_heads"]),
                    "--protein-layers", str(training["protein_layers"]),
                    "--molecule-layers", str(training["molecule_layers"]),
                    "--fusion-layers", str(training["fusion_layers"]),
                    "--chunk-size", str(training["protein_chunk_size"]),
                    "--context-buckets", str(training["context_buckets"]),
                    "--lr", str(training["learning_rate"]),
                    "--batch-size", str(training["batch_size"]),
                    "--grad-accum-steps", str(training["gradient_accumulation_steps"]),
                    "--epochs", str(training["epochs"]),
                    "--warmup-epochs", str(training["warmup_epochs"]),
                    "--patience", str(training["early_stopping_patience"]),
                    "--seed", str(seed), "--cosine-schedule", "--output-dir", str(directory),
                    *variant_flags,
                ]
                print(f"train {variant} seed {seed}", flush=True)
                call_logged(train_arguments, directory / "train.log")
                completion_marker.write_text(
                    json.dumps({"checkpoint": identity(checkpoint)}, indent=2) + "\n",
                    encoding="utf-8",
                )
            if not metrics_path.is_file():
                print(f"evaluate {variant} seed {seed}", flush=True)
                call_logged(
                    [
                        "evaluate", "--checkpoint", str(checkpoint), "--data", str(SOURCE),
                        "--manifest", str(MANIFEST), "--split-manifest", str(SPLIT),
                        "--partition", "test", "--conflict-policy", "median",
                        "--batch-size", str(training["batch_size"]),
                        "--output", str(metrics_path),
                    ],
                    directory / "evaluate.log",
                )


if __name__ == "__main__":
    main()
