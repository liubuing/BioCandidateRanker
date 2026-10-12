"""BioCandidateRanker kcat prediction plugin for the Yeast-MetaTwin Platform.

Drop-in for 04_prediction_plugins/plugin_runtime/. Follows the contract in
.schema and the conventions of the existing UniKP adapter.

What the plugin predicts
    log10(kcat) for one enzyme-substrate pair, from the frozen ESM-2 multimodal
    checkpoint of the BioCandidateRanker release. Only the log10_kcat head is
    trained in that checkpoint; the other three declared tasks carry no labels
    and are never emitted.

Environment note (important)
    The platform's main environment pins torch 2.13, where
    biocandidate.training.load_checkpoint rejects pre-2.6 checkpoints under
    weights_only=True. training.py is part of that project's frozen
    software-implementation-v6 manifest and cannot be edited, so this adapter
    wraps the call in a torch.serialization.safe_globals allowlist - the same
    shim scripts/serve_workbench.py uses. With torch < 2.6 the shim is a no-op.

Context inputs
    The checkpoint has use_context=True with four context fields (organism, EC,
    enzyme type, reaction). Those fields are hashed into the forward pass, so a
    caller that supplies none of them queries a context the model never saw
    during training. The adapter reads them from PluginInput.context when
    present, falls back to the same defaults the project's own CLI uses, and
    reports which fields were supplied so the platform can judge how far a
    request sits from the training distribution.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

from .schema import (
    ApplicabilityReport,
    BenchmarkReport,
    InputTransform,
    PluginInput,
    PluginResult,
    Prediction,
    RunStatus,
    UncertaintyReport,
)


PLUGIN_VERSION = "software-implementation-v6+esm2-t6-seed7"
CAPABILITY = "kcat_prediction"
SEQUENCE_LIMIT = 512
TRUNCATION_STRATEGY = "head_first_512"

# biocandidate maps these to trained token ids; anything else becomes its
# single unknown token, which is a real information loss worth reporting.
CANONICAL_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")
TOLERATED_RESIDUES = frozenset("XBZUOJ")  # accepted, mapped to unknown
CONTEXT_FIELDS = ("organism", "ec", "enzyme_type", "reaction")
CONTEXT_DEFAULTS = {
    "organism": "",
    "ec": "",
    "enzyme_type": "unknown",
    "reaction": "",
}

# Assets are resolved relative to plugin_root, which the platform passes as the
# 04_prediction_plugins directory - the same convention the UniKP adapter uses
# (plugin_root / "UniKP" / "models"). Asset paths here therefore match the
# PluginSpec entries in catalog.py exactly.
DEFAULT_CHECKPOINT = "BioCandidateRanker/assets/esm2-t6-opt-seed7.best.pt"
CHECKPOINT_SHA256 = "5a85d9613cceac0d906dcf17a4440006ce38e899f03c7472db467223104ea0d7"
RANKER_ROOT_ENV = "BIOCANDIDATE_RANKER_ROOT"
RANKER_ROOT_FILE = "ranker_root.txt"
PLUGIN_DIR_NAME = "BioCandidateRanker"

UNTRAINED_TASKS = ("log10_km", "log10_activity", "log10_flux")

_CACHED_MODEL: dict[str, Any] = {}


def resolve_ranker_root(plugin_root: Path) -> Path:
    """Locate the BioCandidateRanker checkout that provides the package.

    Order: the BIOCANDIDATE_RANKER_ROOT environment variable, then the
    ranker_root.txt the installer writes into this plugin's directory, then an
    already-importable biocandidate package. plugin_root is the platform's
    04_prediction_plugins directory.
    """
    configured = os.environ.get(RANKER_ROOT_ENV)
    if configured:
        candidate = Path(configured).expanduser().resolve()
        if (candidate / "src/biocandidate").is_dir():
            return candidate
        raise ValueError(f"{RANKER_ROOT_ENV} does not contain src/biocandidate: {candidate}")

    for pointer in (plugin_root / PLUGIN_DIR_NAME / RANKER_ROOT_FILE,
                    plugin_root / RANKER_ROOT_FILE):
        if pointer.is_file():
            candidate = Path(pointer.read_text(encoding="utf-8").strip()).expanduser()
            if (candidate / "src/biocandidate").is_dir():
                return candidate
            raise ValueError(f"{pointer} does not contain src/biocandidate: {candidate}")

    try:
        import biocandidate  # noqa: F401

        return Path(biocandidate.__file__).resolve().parents[2]
    except ImportError as exc:  # pragma: no cover - configuration error
        raise ValueError(
            f"BioCandidateRanker is not importable and no pointer was found; set "
            f"{RANKER_ROOT_ENV} or write {PLUGIN_DIR_NAME}/{RANKER_ROOT_FILE} "
            f"under {plugin_root}"
        ) from exc


def _import_path(root: Path) -> None:
    source = str(root / "src")
    if source not in sys.path:
        sys.path.insert(0, source)


def _safe_globals():
    """Allowlist the globals torch >= 2.6 refuses when loading BCR checkpoints.

    Only trusted, stdlib-free types are allowed; no arbitrary class from the
    checkpoint resolves. Mirrors scripts/serve_workbench.py.
    """
    import numpy as np
    import torch
    from torch.torch_version import TorchVersion

    allowed: list[Any] = [TorchVersion, np.ndarray, np.dtype]
    allowed.append(getattr(np._core.multiarray, "_reconstruct", None))
    for name in ("float16", "float32", "float64", "int8", "int16", "int32",
                 "int64", "uint8", "bool_"):
        allowed.append(getattr(np, name, None))
        allowed.append(getattr(np.dtypes, f"{name[0].upper()}{name[1:]}DType", None))
    return torch.serialization.safe_globals([item for item in allowed if item is not None])


def prepare_input(
    request: PluginInput,
) -> tuple[str, dict[str, str], list[str], tuple[InputTransform, ...]]:
    """Validate and normalise the request; report every transform applied."""
    if request.capability != CAPABILITY:
        raise ValueError(f"unsupported BioCandidateRanker capability: {request.capability}")
    if request.sequence is None or request.substrate_smiles is None:
        raise ValueError("BioCandidateRanker requires sequence and substrate_smiles")

    sequence = re.sub(r"\s+", "", request.sequence.upper())
    if not sequence:
        raise ValueError("sequence is empty after whitespace removal")

    invalid = sorted(set(sequence) - CANONICAL_RESIDUES - TOLERATED_RESIDUES)
    if invalid:
        raise ValueError(f"sequence contains invalid residues: {''.join(invalid)}")

    transforms: list[InputTransform] = []
    effective_sequence = sequence
    if len(sequence) > SEQUENCE_LIMIT:
        effective_sequence = sequence[:SEQUENCE_LIMIT]
        transforms.append(
            InputTransform("sequence", len(sequence), len(effective_sequence), True, TRUNCATION_STRATEGY)
        )
    else:
        transforms.append(
            InputTransform("sequence", len(sequence), len(sequence), False, "none")
        )

    # Count positions, not distinct residue types: a 36-residue canonical
    # sequence contains about 20 distinct residues but zero substitutions.
    non_canonical = sum(1 for residue in sequence if residue not in CANONICAL_RESIDUES)
    if non_canonical:
        transforms.append(
            InputTransform(
                "sequence:non_canonical_residues",
                len(sequence),
                len(sequence),
                False,
                f"mapped_to_unknown_token({non_canonical})",
            )
        )

    context: dict[str, str] = {}
    supplied: list[str] = []
    raw_context = request.context or {}
    for field in CONTEXT_FIELDS:
        value = raw_context.get(field)
        if value is None or str(value).strip() == "":
            context[field] = CONTEXT_DEFAULTS[field]
        else:
            context[field] = str(value).strip()
            supplied.append(field)

    messages = []
    absent = [field for field in CONTEXT_FIELDS if field not in supplied]
    if absent:
        messages.append(
            "context fields absent and left at their project defaults: "
            + ", ".join(absent)
            + "; these fields are hashed into the forward pass, so an empty context "
              "is a distribution shift from the training corpus"
        )
    if non_canonical:
        messages.append(
            f"{non_canonical} non-canonical residue(s) mapped to the unknown token"
        )
    return effective_sequence, context, messages, tuple(transforms)


def _load_model(checkpoint_path: Path, ranker_root: Path):
    """Load and cache the BCR checkpoint for this process."""
    key = str(checkpoint_path)
    cached = _CACHED_MODEL.get(key)
    if cached is not None:
        return cached

    _import_path(ranker_root)
    import torch

    from biocandidate.training import load_checkpoint

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    with _safe_globals():
        model, checkpoint = load_checkpoint(checkpoint_path, device)
    _CACHED_MODEL[key] = (model, checkpoint, device)
    return _CACHED_MODEL[key]


def _predict_log10_kcat(
    sequence: str,
    smiles: str,
    context: dict[str, str],
    checkpoint_path: Path,
    ranker_root: Path,
) -> tuple[float, float, str]:
    """Return (mean, standard_deviation, unit) in log10 space."""
    model, checkpoint, device = _load_model(checkpoint_path, ranker_root)

    import torch

    from biocandidate.data import EnzymeSubstrateRecord
    from biocandidate.data.collate import EnzymeSubstrateCollator

    if "log10_kcat" not in model.config.task_names:
        raise ValueError("checkpoint does not declare a log10_kcat task")
    observed = (
        checkpoint.get("metrics", {})
        .get("validation_tasks", {})
        .get("log10_kcat", {})
        .get("count", 0)
    )
    if not observed:
        raise ValueError("log10_kcat head carries no training labels in this checkpoint")

    record = EnzymeSubstrateRecord(
        candidate_id="plugin-request",
        sequence=sequence,
        substrate_smiles=smiles,
        substrate_name="",
        organism=context["organism"],
        ec=context["ec"],
        enzyme_type=context["enzyme_type"],
        reaction=context["reaction"],
        fba_context=(),
        fba_feature_metadata=None,
        source_dataset="yeast-metatwin-platform-plugin",
        source_row=0,
    )
    collator = EnzymeSubstrateCollator(
        context_buckets=model.config.context_buckets,
        fba_context_dim=model.config.fba_context_dim,
    )
    batch = {name: value.to(device) for name, value in collator([record]).items()}
    model.eval()
    with torch.no_grad():
        output = model(batch)

    index = list(model.config.task_names).index("log10_kcat")
    return (
        float(output["mean"][0, index]),
        float(output["standard_deviation"][0, index]),
        "s^-1",
    )


def assess_applicability(
    transforms: tuple[InputTransform, ...], context_supplied: list[str]
) -> ApplicabilityReport:
    """BCR ships no out-of-distribution model, so domain status is not asserted.

    The reasons list what a caller can still act on: whether the sequence was
    truncated and how much of the context was supplied.
    """
    reasons = [
        "BioCandidateRanker has no applicability or OOD model; domain status is not asserted",
    ]
    truncated = [item for item in transforms if item.field == "sequence" and item.truncated]
    if truncated:
        reasons.append(
            f"sequence truncated {truncated[0].original_length} -> "
            f"{truncated[0].effective_length} residues ({truncated[0].strategy})"
        )
    missing = [field for field in CONTEXT_FIELDS if field not in context_supplied]
    if missing:
        reasons.append("context supplied for: " + (", ".join(context_supplied) or "none")
                       + "; defaulted: " + ", ".join(missing))
    else:
        reasons.append("all four context fields supplied")
    return ApplicabilityReport("not_available", None, None, "not_implemented", tuple(reasons))


def predict(request: PluginInput, plugin_root: Path) -> PluginResult:
    """Entry point matching the UniKP adapter signature."""
    try:
        sequence, context, messages, transforms = prepare_input(request)
    except Exception as exc:
        # Input problems are the caller's to fix, not a runtime failure.
        return PluginResult(
            request_id=request.request_id,
            plugin="BioCandidateRanker",
            plugin_version=PLUGIN_VERSION,
            capability=request.capability,
            status=RunStatus.UNSUPPORTED,
            messages=(f"{type(exc).__name__}: {exc}",),
        )

    try:
        ranker_root = resolve_ranker_root(plugin_root)
        checkpoint_path = (plugin_root / DEFAULT_CHECKPOINT).resolve()
        if not checkpoint_path.is_file():
            raise FileNotFoundError(
                f"checkpoint asset missing: {checkpoint_path}; run the plugin installer "
                "or copy esm2-t6-opt-seed7.best.pt into the plugin assets directory"
            )
        mean, deviation, unit = _predict_log10_kcat(
            sequence, request.substrate_smiles, context, checkpoint_path, ranker_root
        )
        context_supplied = [field for field in CONTEXT_FIELDS if context[field] != CONTEXT_DEFAULTS[field]]
        uncertainty = UncertaintyReport(
            "available",
            "heteroscedastic Gaussian head; interval is the uncalibrated 1-sigma band in "
            "log10 space (no matching calibration artifact ships with this release)",
            deviation,
            interval=(mean - deviation, mean + deviation),
            calibrated=False,
        )
        return PluginResult(
            request_id=request.request_id,
            plugin="BioCandidateRanker",
            plugin_version=PLUGIN_VERSION,
            capability=request.capability,
            status=RunStatus.READY,
            predictions=(
                Prediction("log10_kcat", mean, f"log10({unit})", "log10"),
                Prediction("kcat", math.pow(10.0, mean), unit),
            ),
            transforms=transforms,
            applicability=assess_applicability(transforms, context_supplied),
            uncertainty=uncertainty,
            benchmark=BenchmarkReport(),
            messages=tuple(messages) + (
                "Prediction is for prioritization and is not curated kinetic evidence.",
                "Only log10_kcat is trained in this checkpoint; "
                + ", ".join(UNTRAINED_TASKS)
                + " carry no labels and are not emitted.",
            ),
        )
    except Exception as exc:
        return PluginResult(
            request_id=request.request_id,
            plugin="BioCandidateRanker",
            plugin_version=PLUGIN_VERSION,
            capability=request.capability,
            status=RunStatus.ERROR,
            transforms=transforms,
            applicability=ApplicabilityReport(
                "not_available", None, None, "runtime failed before domain assessment"
            ),
            uncertainty=UncertaintyReport("not_available"),
            messages=(f"{type(exc).__name__}: {exc}",),
        )


def checkpoint_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
