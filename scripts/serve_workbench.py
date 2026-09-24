"""Serve an interactive workbench for the project's own capabilities.

The page can run real kcat predictions through the frozen checkpoints, browse
the checkpoints with their recorded validation metrics, and execute the
read-only governance scripts. It exists because the capabilities are otherwise
locked behind CLI invocations with JSON files.

Endpoints (all JSON unless noted):
  GET  /                the workbench page
  GET  /dashboard       the readiness dashboard, if it has been generated
  GET  /api/health      server and device status
  GET  /api/models      checkpoints with config and validation metrics
  POST /api/predict     {model, candidates[], calibration_artifact?}
  POST /api/governance  {action} from a fixed allowlist

Local-only by design: the server binds to 127.0.0.1 and the governance actions
are an allowlist, not a shell.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

INDEX_PATH = Path(__file__).resolve().parent / "workbench_static" / "index.html"
DEFAULT_CHECKPOINT_GLOBS = ("artifacts/*/best.pt", "artifacts/*/*/best.pt")

GOVERNANCE_ACTIONS: dict[str, list[str]] = {
    "license_audit": [sys.executable, "scripts/audit_development_corpus_license.py"],
    "pool_dashboard": [sys.executable, "scripts/prospective_pool_dashboard.py"],
    "request_status": [sys.executable, "scripts/external_request_status.py"],
    "figure1_regenerate": [sys.executable, "scripts/figure1_data_task_map.py"],
    "archive_verify": [
        sys.executable,
        "scripts/verify_release_archive.py",
        "artifacts/release-archive/software-implementation-v6.zip",
    ],
}
GOVERNANCE_TIMEOUT_SECONDS = 120

WARNING_TEXT = "Computational ranking only; not experimental validation."
CHECKPOINT_SUFFIX = ".pt"


def _safe_load_globals():
    """Allowlist the globals torch 2.6+ refuses by default.

    Checkpoints saved before torch 2.6 embed torch.torch_version.TorchVersion,
    and the feature-MLP checkpoints embed numpy object reconstruction, which
    weights_only loading rejects. training.load_checkpoint is part of the
    frozen release manifest, so the shim lives here instead. Only stdlib-free
    trusted types are allowed: no arbitrary class from the checkpoint resolves.
    """
    import numpy as np
    import torch
    from torch.torch_version import TorchVersion

    allowed = [TorchVersion, np.ndarray, np.dtype]
    for name in ("_reconstruct",):
        allowed.append(getattr(np._core.multiarray, name, None))
    for dtype_name in ("float16", "float32", "float64", "int8", "int16",
                       "int32", "int64", "uint8", "bool_"):
        allowed.append(getattr(np, f"{dtype_name}", None))
        allowed.append(getattr(np.dtypes, f"{dtype_name[0].upper()}{dtype_name[1:]}DType", None))
    allowed = [item for item in allowed if item is not None]
    return torch.serialization.safe_globals(allowed)


class ModelCache:
    """Load checkpoints on demand and keep a small pool of live models."""

    def __init__(self, root: Path, device: str, capacity: int = 2):
        self.root = root
        self.device = device
        self.capacity = capacity
        self._entries: dict[str, tuple[Any, dict, float]] = {}
        self._order: list[str] = []
        self.lock = threading.Lock()

    def get(self, relative: str):
        import torch
        from biocandidate.training import load_checkpoint

        path = self.resolve(relative)
        stamp = path.stat().st_mtime
        with self.lock:
            entry = self._entries.get(relative)
            if entry and entry[2] == stamp:
                self._order.remove(relative)
                self._order.append(relative)
                return entry[0], entry[1]
        with _safe_load_globals():
            model, checkpoint = load_checkpoint(path, torch.device(self.device))
        with self.lock:
            self._entries[relative] = (model, checkpoint, stamp)
            self._order.append(relative)
            while len(self._order) > self.capacity:
                self._entries.pop(self._order.pop(0), None)
        return model, checkpoint

    def resolve(self, relative: str) -> Path:
        """Refuse anything outside the repository or outside the checkpoint set."""
        candidate = self.resolve_data(relative)
        if candidate.suffix != CHECKPOINT_SUFFIX:
            raise ValueError("model path must point at a checkpoint file")
        return candidate

    def resolve_data(self, relative: str) -> Path:
        """Any repository file: checkpoints, calibration artifacts, and so on."""
        candidate = (self.root / relative).resolve()
        if not candidate.is_relative_to(self.root):
            raise ValueError("path must stay inside the repository")
        if not candidate.is_file():
            raise ValueError(f"file not found: {relative}")
        return candidate


def list_models(root: Path) -> list[dict[str, Any]]:
    """Summarize every checkpoint without loading it onto the accelerator."""
    import torch

    models: list[dict[str, Any]] = []
    for pattern in DEFAULT_CHECKPOINT_GLOBS:
        for path in sorted(root.glob(pattern)):
            try:
                with _safe_load_globals():
                    payload = torch.load(path, map_location="cpu", weights_only=True)
                if "model_config" not in payload:
                    # Feature-MLP baselines store a different layout; they are
                    # listed for completeness but are not prediction models.
                    models.append(
                        {
                            "path": path.relative_to(root).as_posix(),
                            "tasks": [],
                            "feature_baseline": payload.get("feature"),
                            "epoch": payload.get("selected_epoch"),
                            "size_bytes": path.stat().st_size,
                            "validation": {},
                        }
                    )
                    continue
                config = payload["model_config"]
                validation = payload.get("metrics", {}).get("validation_tasks", {})
                models.append(
                    {
                        "path": path.relative_to(root).as_posix(),
                        "tasks": list(config.get("task_names", [])),
                        "protein_encoder": config.get("protein_encoder"),
                        "esm2_model_name": config.get("esm2_model_name"),
                        "d_model": config.get("d_model"),
                        "fusion_mode": config.get("fusion_mode"),
                        "uncertainty_mode": config.get("uncertainty_mode"),
                        "epoch": payload.get("epoch"),
                        "size_bytes": path.stat().st_size,
                        "validation": {
                            task: {
                                key: validation[task][key]
                                for key in ("count", "mae", "rmse", "pearson")
                                if key in validation[task]
                            }
                            for task in config.get("task_names", [])
                            if task in validation
                        },
                    }
                )
            except Exception as error:  # noqa: BLE001 - listing must survive one bad file
                models.append(
                    {
                        "path": path.relative_to(root).as_posix(),
                        "error": f"{type(error).__name__}: {error}",
                    }
                )
    return models


def parse_candidates(payload: Any) -> list[dict[str, Any]]:
    """Validate and normalize the candidate list from a request body."""
    if isinstance(payload, dict):
        rows = payload.get("candidates")
    else:
        rows = payload
    if not isinstance(rows, list) or not rows:
        raise ValueError("request must carry a non-empty candidates list")
    normalized: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"candidate {index} is not an object")
        sequence = str(row.get("sequence", "")).strip()
        smiles = str(row.get("substrate_smiles", "")).strip()
        if not sequence:
            raise ValueError(f"candidate {index} is missing sequence")
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYXBZJUOacdefghiklmnpqrstvwyxbzjuo\s]+", sequence):
            raise ValueError(
                f"candidate {index} sequence contains non-residue characters"
            )
        if not smiles:
            raise ValueError(f"candidate {index} is missing substrate_smiles")
        normalized.append(
            {
                "candidate_id": str(row.get("candidate_id", index)),
                "sequence": sequence.upper().replace(" ", ""),
                "substrate_smiles": smiles,
                "substrate_name": str(row.get("substrate_name", "")),
                "organism": str(row.get("organism", "")),
                "ec": str(row.get("ec", "")),
                "enzyme_type": str(row.get("enzyme_type", "unknown")),
                "reaction": str(row.get("reaction", "")),
            }
        )
    return normalized


def run_prediction(
    cache: ModelCache,
    model_path: str,
    candidates: list[dict[str, Any]],
    calibration_artifact: str | None,
    fba_payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Mirror cli.predict_command's in-memory flow and result shape.

    cli.py is part of the frozen release manifest, so the shared logic is
    re-expressed here rather than refactored; a test pins the two together.
    An optional fba_payload attaches governed FBA context features produced by
    biocandidate.fba_context.
    """
    import torch
    from biocandidate.data import FBAFeatureMetadata, EnzymeSubstrateRecord

    model, checkpoint = cache.get(model_path)

    fba_metadata = None
    fba_context: tuple[float, ...] = ()
    if fba_payload:
        body = fba_payload.get("metadata", {})
        features = fba_payload.get("features")
        if not isinstance(features, list) or not features:
            raise ValueError("fba payload must carry a non-empty features list")
        fba_metadata = FBAFeatureMetadata(
            schema_version=int(body.get("schema_version", 1)),
            feature_ids=tuple(body.get("feature_ids", ())),
            model_id=str(body.get("model_id", "")),
            solver_id=str(body.get("solver_id", "")),
            objective_id=str(body.get("objective_id", "")),
            condition_id=str(body.get("condition_id", "")),
        )
        if len(fba_metadata.feature_ids) != model.config.fba_context_dim:
            raise ValueError(
                f"FBA feature width {len(fba_metadata.feature_ids)} does not match "
                f"checkpoint configured width {model.config.fba_context_dim}"
            )
        if len(features) != len(fba_metadata.feature_ids):
            raise ValueError(
                f"fba features length {len(features)} != feature_ids length "
                f"{len(fba_metadata.feature_ids)}"
            )
        fba_context = tuple(float(value) for value in features)

    records = [
        EnzymeSubstrateRecord(
            candidate_id=row["candidate_id"],
            sequence=row["sequence"],
            substrate_smiles=row["substrate_smiles"],
            substrate_name=row["substrate_name"],
            organism=row["organism"],
            ec=row["ec"],
            enzyme_type=row["enzyme_type"],
            reaction=row["reaction"],
            fba_context=fba_context,
            fba_feature_metadata=fba_metadata,
            source_dataset="workbench",
            source_row=index,
        )
        for index, row in enumerate(candidates)
    ]

    from biocandidate.data.collate import EnzymeSubstrateCollator

    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
    )
    batch = {
        name: value.to(cache.device) for name, value in collator(records).items()
    }
    model.eval()
    with torch.no_grad():
        output = model(batch)

    calibration_payload = None
    calibration_task = None
    calibrator = None
    conformal_intervals: dict[str, Any] = {}
    if calibration_artifact:
        from biocandidate.calibration import (
            conformal_calibrator_from_dict,
            load_calibration_artifact,
        )

        calibrator, calibration_payload = load_calibration_artifact(
            cache.resolve_data(calibration_artifact)
        )
        identities = calibration_payload["identities"]
        calibration_task = identities.get("task")
        if calibration_task not in model.config.task_names:
            raise ValueError("calibration artifact task is not defined by the checkpoint")
        checkpoint_path = cache.resolve(model_path)
        expected = {
            "checkpoint": {
                "sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
                "size_bytes": checkpoint_path.stat().st_size,
            },
            "task": calibration_task,
        }
        if identities != expected:
            raise ValueError("calibration artifact identity mismatch")
        conformal_intervals = {
            name: conformal_calibrator_from_dict(body)
            for name, body in calibration_payload.get("conformal_intervals", {}).items()
        }

    validation_tasks = checkpoint.get("metrics", {}).get("validation_tasks", {})
    predictions = []
    for row, record in enumerate(records):
        tasks: dict[str, Any] = {}
        for index, task in enumerate(model.config.task_names):
            observed = int(validation_tasks.get(task, {}).get("count", 0))
            if observed == 0:
                tasks[task] = {
                    "status": "untrained_no_labels",
                    "mean": None,
                    "standard_deviation": None,
                }
                continue
            mean_value = output["mean"][row, index]
            deviation_value = output["standard_deviation"][row, index]
            mean = float(mean_value)
            deviation = float(deviation_value)
            intervals: dict[str, Any] | None = None
            calibrated_deviation: float | None = None
            if calibrator is not None and task == calibration_task:
                calibrated_mean, calibrated = calibrator.calibrate(
                    mean_value.reshape(1), deviation_value.reshape(1)
                )
                calibrated_deviation = float(calibrated[0])
                intervals = {}
                for name, conformal in conformal_intervals.items():
                    bounds = conformal.interval(
                        calibrated_mean,
                        calibrated if conformal.normalized else None,
                    )
                    intervals[name] = {
                        "lower": None if not math.isfinite(float(bounds[0][0])) else float(bounds[0][0]),
                        "upper": None if not math.isfinite(float(bounds[1][0])) else float(bounds[1][0]),
                    }
            tasks[task] = {
                "status": "trained",
                "mean": mean,
                "standard_deviation": deviation,
                "raw_standard_deviation": deviation,
                "calibrated_standard_deviation": calibrated_deviation,
                "conformal_intervals": intervals,
            }
        predictions.append({"candidate_id": record.candidate_id, "tasks": tasks})

    result = {
        "checkpoint": model_path,
        "input_identity": {
            "row_count": len(candidates),
            "sha256": hashlib.sha256(
                json.dumps(candidates, sort_keys=True).encode("utf-8")
            ).hexdigest(),
        },
        "predictions": predictions,
        "warning": WARNING_TEXT,
    }
    if fba_metadata is not None:
        result["fba"] = {
            "attached": True,
            "metadata": {
                "model_id": fba_metadata.model_id,
                "objective_id": fba_metadata.objective_id,
                "condition_id": fba_metadata.condition_id,
                "feature_ids": list(fba_metadata.feature_ids),
            },
            "note": (
                "FBA context was attached, but no shipped checkpoint was trained "
                "with flux labels: the flux pathway is untrained, so this "
                "demonstrates the plumbing rather than a validated input."
            ),
        }
    return result


def run_governance(root: Path, action: str) -> dict[str, Any]:
    """Execute one allowlisted governance script and return its JSON output."""
    if action not in GOVERNANCE_ACTIONS:
        known = ", ".join(sorted(GOVERNANCE_ACTIONS))
        raise ValueError(f"unknown governance action {action!r}; known: {known}")
    command = GOVERNANCE_ACTIONS[action]
    completed = subprocess.run(
        command,
        cwd=root,
        capture_output=True,
        text=True,
        timeout=GOVERNANCE_TIMEOUT_SECONDS,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "action": action,
            "ok": False,
            "returncode": completed.returncode,
            "stderr": completed.stderr[-4000:],
        }
    try:
        parsed = json.loads(completed.stdout)
    except json.JSONDecodeError:
        parsed = {"raw": completed.stdout[-4000:]}
    return {"action": action, "ok": True, "result": parsed}


_FBA_MODEL_CACHE: dict[str, Any] = {}
_FBA_LOCK = threading.Lock()


def _fba_model(model_id: str, root: Path):
    """Load each registered GEM once; cobra models are reused via context managers."""
    from biocandidate import fba_context

    with _FBA_LOCK:
        if model_id not in _FBA_MODEL_CACHE:
            config = fba_context.load_model_config(model_id, root=root)
            resolved = root / config["path"]
            _FBA_MODEL_CACHE[model_id] = {
                "model": fba_context.load_model(model_id, resolved),
                "sha256": config["sha256"],
            }
        return _FBA_MODEL_CACHE[model_id]


def fba_status(root: Path) -> dict[str, Any]:
    from biocandidate import fba_context

    models = {}
    for model_id in fba_context.model_ids():
        entry = {}
        try:
            config = fba_context.load_model_config(model_id, root=root)
            resolved = root / config["path"]
            entry.update(
                {
                    "available": resolved.is_file(),
                    "sha256": config["sha256"],
                    "objective_id": fba_context.MODEL_REGISTRY[model_id]["objective_id"],
                    "citation": fba_context.MODEL_REGISTRY[model_id]["citation"],
                    "presets": fba_context.presets_for(model_id),
                }
            )
        except Exception as error:  # noqa: BLE001 - status must survive a bad config
            entry.update({"available": False, "error": str(error)})
        models[model_id] = entry
    return {
        "models": models,
        "claim_boundary": fba_context.CLAIM_BOUNDARY,
        "feature_width_note": (
            "Both models emit 8-wide vectors matching ModelConfig.fba_context_dim."
        ),
    }


def fba_run(root: Path, body: dict[str, Any]) -> dict[str, Any]:
    from biocandidate import fba_context

    model_id = str(body.get("model_id", "iML1515"))
    entry = _fba_model(model_id, root)
    return fba_context.run_fba_context(
        model_id=model_id,
        preset=str(body.get("preset", "baseline")),
        glucose=float(body.get("glucose", fba_context.DEFAULT_GLUCOSE)),
        oxygen=float(body.get("oxygen", fba_context.DEFAULT_OXYGEN)),
        growth_min=float(body.get("growth_min", 0.0)),
        with_fva=bool(body.get("with_fva", False)),
        model=entry["model"],
        model_sha256=entry["sha256"],
    )


class WorkbenchHandler(BaseHTTPRequestHandler):
    server_version = "BioCandidateRankerWorkbench/1.0"
    cache: ModelCache
    root: Path

    def log_message(self, fmt: str, *args: Any) -> None:  # keep the console quiet
        pass

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, status: int, payload: Any) -> None:
        body = json.dumps(payload, indent=2).encode("utf-8")
        self._send(status, body, "application/json; charset=utf-8")

    def do_GET(self) -> None:  # noqa: N802 - stdlib API
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            if INDEX_PATH.is_file():
                self._send(200, INDEX_PATH.read_bytes(), "text/html; charset=utf-8")
            else:
                self._send_json(404, {"error": f"missing {INDEX_PATH}"})
        elif path == "/api/health":
            self._send_json(
                200,
                {
                    "status": "ok",
                    "device": self.cache.device,
                    "loaded_models": len(self.cache._entries),
                },
            )
        elif path == "/api/models":
            try:
                self._send_json(200, {"models": list_models(self.root)})
            except Exception as error:  # noqa: BLE001
                self._send_json(500, {"error": str(error)})
        elif path == "/api/fba/status":
            try:
                self._send_json(200, fba_status(self.root))
            except Exception as error:  # noqa: BLE001
                self._send_json(500, {"error": str(error)})
        elif path == "/dashboard":
            target = self.root / "artifacts/dashboard/readiness-dashboard.html"
            if target.is_file():
                self._send(200, target.read_bytes(), "text/html; charset=utf-8")
            else:
                self._send_json(
                    404,
                    {
                        "error": "dashboard not generated yet",
                        "remedy": "python scripts/build_readiness_dashboard.py",
                    },
                )
        else:
            self._send_json(404, {"error": f"no route for {path}"})

    def do_POST(self) -> None:  # noqa: N802 - stdlib API
        path = self.path.split("?", 1)[0]
        length = int(self.headers.get("Content-Length") or 0)
        if length > 8_000_000:
            self._send_json(413, {"error": "request body too large"})
            return
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
        except json.JSONDecodeError as error:
            self._send_json(400, {"error": f"invalid JSON: {error}"})
            return

        if path == "/api/predict":
            try:
                candidates = parse_candidates(body)
                model_path = body.get("model")
                if not model_path:
                    raise ValueError("request must name a model checkpoint")
                result = run_prediction(
                    self.cache,
                    str(model_path),
                    candidates,
                    body.get("calibration_artifact"),
                    fba_payload=body.get("fba"),
                )
                self._send_json(200, result)
            except Exception as error:  # noqa: BLE001 - surfaced to the page
                self._send_json(400, {"error": f"{type(error).__name__}: {error}"})
        elif path == "/api/fba/run":
            try:
                self._send_json(200, fba_run(self.root, body))
            except Exception as error:  # noqa: BLE001
                self._send_json(400, {"error": f"{type(error).__name__}: {error}"})
        elif path == "/api/governance":
            try:
                self._send_json(200, run_governance(self.root, str(body.get("action"))))
            except Exception as error:  # noqa: BLE001
                self._send_json(400, {"error": str(error)})
        else:
            self._send_json(404, {"error": f"no route for {path}"})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8808)
    parser.add_argument(
        "--device",
        default=None,
        help="torch device (default: cuda when available, else cpu)",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    if args.device:
        device = args.device
    else:
        import torch

        device = "cuda" if torch.cuda.is_available() else "cpu"

    WorkbenchHandler.root = root
    WorkbenchHandler.cache = ModelCache(root, device)

    server = ThreadingHTTPServer((args.host, args.port), WorkbenchHandler)
    print(
        json.dumps(
            {
                "serving": f"http://{args.host}:{args.port}/",
                "root": str(root),
                "device": device,
                "dashboard": "http://"
                + f"{args.host}:{args.port}"
                + "/dashboard",
            },
            indent=2,
        )
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
