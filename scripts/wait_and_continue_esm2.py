"""Continue the frozen modality diagnostic after the first interactive run exits.

Waits on a training PID passed via --pid, fails closed if it exited nonzero,
either checkpoint is missing, or the run stopped before the earliest plausible
patience stop (epoch >= 5); otherwise records ``training-complete.json`` with
checkpoint identity and re-invokes ``run_esm2_modality_diagnostic.py`` to run
the remaining variants. Governance role: supervised continuation helper only --
no training parameters are read or changed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import psutil
import torch


ROOT = Path(__file__).resolve().parents[1]
FIRST = ROOT / "artifacts" / "esm2-modality-diagnostic" / "protein-only-seed7"


def identity(path: Path) -> dict[str, int | str]:
    """Return size and SHA-256 identity for a produced checkpoint file."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"size_bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def main() -> None:
    """Wait for the first run, validate its outcome, then continue the diagnostic."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--pid", type=int, required=True)
    arguments = parser.parse_args()
    print(f"Waiting for training PID {arguments.pid}", flush=True)
    exit_code = psutil.Process(arguments.pid).wait(timeout=None)
    print(f"First run exited with code {exit_code}", flush=True)
    if exit_code != 0:
        raise RuntimeError("First training run failed; inspect its interactive output")
    checkpoint = FIRST / "best.pt"
    latest = FIRST / "latest.pt"
    if not checkpoint.is_file() or not latest.is_file():
        raise RuntimeError("First run exited without both checkpoints")
    payload = torch.load(latest, map_location="cpu", weights_only=False)
    if payload["epoch"] < 5:
        raise RuntimeError("First run ended before the earliest plausible patience stop")
    (FIRST / "training-complete.json").write_text(
        json.dumps({"checkpoint": identity(checkpoint), "last_epoch": payload["epoch"]}, indent=2)
        + "\n",
        encoding="utf-8",
    )
    subprocess.run(
        [sys.executable, "-u", str(ROOT / "scripts" / "run_esm2_modality_diagnostic.py")],
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
