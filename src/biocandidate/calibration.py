"""Validation-only uncertainty calibration and interval evaluation utilities.

Fitting functions in this module accept one labeled calibration set.  Applying a
fitted calibrator accepts predictions only, which keeps test labels out of the
calibration interface.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import torch


Tensor = torch.Tensor
_FORMAT_VERSION = 2
_MODEL_ORDER = ("identity", "scalar", "affine_variance")


def _vector(value: Tensor, name: str, *, positive: bool = False) -> Tensor:
    if not isinstance(value, torch.Tensor):
        raise TypeError(f"{name} must be a torch.Tensor")
    if value.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if value.numel() == 0:
        raise ValueError(f"{name} must not be empty")
    if value.dtype not in (torch.float32, torch.float64):
        raise TypeError(f"{name} must have float32 or float64 dtype")
    if not torch.isfinite(value).all():
        raise ValueError(f"{name} must contain only finite values")
    if positive and (value <= 0).any():
        raise ValueError(f"{name} must contain only positive values")
    return value


def _same_vectors(**values: Tensor) -> None:
    shapes = {tuple(value.shape) for value in values.values()}
    devices = {value.device for value in values.values()}
    if len(shapes) != 1:
        raise ValueError(f"{', '.join(values)} must have identical shapes")
    if len(devices) != 1:
        raise ValueError(f"{', '.join(values)} must be on the same device")


def _probability(value: float, name: str = "alpha") -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result) or not 0.0 < result < 1.0:
        raise ValueError(f"{name} must be strictly between zero and one")
    return result


@dataclass(frozen=True)
class IdentityGaussianCalibrator:
    """Leave Gaussian predictive means and deviations unchanged."""

    kind = "identity"

    @classmethod
    def fit(
        cls, mean: Tensor, standard_deviation: Tensor, labels: Tensor
    ) -> IdentityGaussianCalibrator:
        _validated_gaussian_sample(mean, standard_deviation, labels)
        return cls()

    def calibrate(self, mean: Tensor, standard_deviation: Tensor) -> tuple[Tensor, Tensor]:
        _validated_predictions(mean, standard_deviation)
        return mean, standard_deviation

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind}


@dataclass(frozen=True)
class ScalarGaussianCalibrator:
    """Multiply predictive standard deviations by one validation-fitted scale."""

    scale: float
    kind = "scalar"

    def __post_init__(self) -> None:
        if isinstance(self.scale, bool) or not isinstance(self.scale, (int, float)):
            raise TypeError("scale must be a real number")
        if not math.isfinite(self.scale) or self.scale <= 0.0:
            raise ValueError("scale must be finite and positive")

    @classmethod
    def fit(
        cls, mean: Tensor, standard_deviation: Tensor, labels: Tensor
    ) -> ScalarGaussianCalibrator:
        _validated_gaussian_sample(mean, standard_deviation, labels)
        scale = torch.sqrt(torch.mean(((labels - mean) / standard_deviation).square()))
        if not torch.isfinite(scale) or scale <= 0:
            raise ValueError("scalar calibration requires at least one nonzero residual")
        return cls(float(scale))

    def calibrate(self, mean: Tensor, standard_deviation: Tensor) -> tuple[Tensor, Tensor]:
        _validated_predictions(mean, standard_deviation)
        return mean, standard_deviation * self.scale

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "scale": self.scale}


@dataclass(frozen=True)
class AffineVarianceGaussianCalibrator:
    """Transform predictive variance as ``scale * variance + offset``."""

    scale: float
    offset: float
    kind = "affine_variance"

    def __post_init__(self) -> None:
        for value, name in ((self.scale, "scale"), (self.offset, "offset")):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{name} must be a real number")
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.scale == 0.0 and self.offset == 0.0:
            raise ValueError("scale and offset cannot both be zero")

    @classmethod
    def fit(
        cls,
        mean: Tensor,
        standard_deviation: Tensor,
        labels: Tensor,
        *,
        max_steps: int = 250,
    ) -> AffineVarianceGaussianCalibrator:
        _validated_gaussian_sample(mean, standard_deviation, labels)
        if isinstance(max_steps, bool) or not isinstance(max_steps, int):
            raise TypeError("max_steps must be an integer")
        if max_steps < 1:
            raise ValueError("max_steps must be positive")
        residual_square = (labels.detach().double() - mean.detach().double()).square()
        variance = standard_deviation.detach().double().square()
        if residual_square.max() == 0:
            raise ValueError("affine calibration requires at least one nonzero residual")

        # Softplus gives a smooth nonnegative parameterization for deterministic LBFGS.
        raw = torch.tensor([0.541324854612918, -4.0], dtype=torch.float64, requires_grad=True)
        optimizer = torch.optim.LBFGS(
            [raw],
            max_iter=max_steps,
            tolerance_grad=1e-12,
            tolerance_change=1e-14,
            line_search_fn="strong_wolfe",
        )

        def closure() -> Tensor:
            optimizer.zero_grad()
            parameters = torch.nn.functional.softplus(raw)
            calibrated_variance = parameters[0] * variance + parameters[1]
            loss = 0.5 * torch.mean(
                torch.log(calibrated_variance) + residual_square / calibrated_variance
            )
            loss.backward()
            return loss

        optimizer.step(closure)
        parameters = torch.nn.functional.softplus(raw.detach())
        return cls(float(parameters[0]), float(parameters[1]))

    def calibrate(self, mean: Tensor, standard_deviation: Tensor) -> tuple[Tensor, Tensor]:
        _validated_predictions(mean, standard_deviation)
        variance = self.scale * standard_deviation.square() + self.offset
        return mean, torch.sqrt(variance)

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "offset": self.offset, "scale": self.scale}


GaussianCalibrator = (
    IdentityGaussianCalibrator | ScalarGaussianCalibrator | AffineVarianceGaussianCalibrator
)


def _validated_predictions(mean: Tensor, standard_deviation: Tensor) -> None:
    _vector(mean, "mean")
    _vector(standard_deviation, "standard_deviation", positive=True)
    _same_vectors(mean=mean, standard_deviation=standard_deviation)


def _validated_gaussian_sample(mean: Tensor, standard_deviation: Tensor, labels: Tensor) -> None:
    _validated_predictions(mean, standard_deviation)
    _vector(labels, "labels")
    _same_vectors(mean=mean, standard_deviation=standard_deviation, labels=labels)


def gaussian_nll(mean: Tensor, standard_deviation: Tensor, labels: Tensor) -> Tensor:
    """Return mean Gaussian negative log likelihood, including its constant term."""
    _validated_gaussian_sample(mean, standard_deviation, labels)
    z = (labels - mean) / standard_deviation
    return torch.mean(
        0.5 * z.square() + torch.log(standard_deviation) + 0.5 * math.log(2 * math.pi)
    )


def gaussian_crps(mean: Tensor, standard_deviation: Tensor, labels: Tensor) -> Tensor:
    """Return the mean proper continuous ranked probability score for Gaussians."""
    _validated_gaussian_sample(mean, standard_deviation, labels)
    z = (labels - mean) / standard_deviation
    phi = torch.exp(-0.5 * z.square()) / math.sqrt(2.0 * math.pi)
    cdf = 0.5 * (1.0 + torch.erf(z / math.sqrt(2.0)))
    return torch.mean(
        standard_deviation * (z * (2.0 * cdf - 1.0) + 2.0 * phi - 1.0 / math.sqrt(math.pi))
    )


def interval_score(lower: Tensor, upper: Tensor, labels: Tensor, *, alpha: float) -> Tensor:
    """Return the mean proper central interval score at miscoverage ``alpha``."""
    alpha = _probability(alpha)
    _vector(lower, "lower")
    _vector(upper, "upper")
    _vector(labels, "labels")
    _same_vectors(lower=lower, upper=upper, labels=labels)
    if (lower > upper).any():
        raise ValueError("lower must not exceed upper")
    score = upper - lower
    score = score + (2.0 / alpha) * (lower - labels) * (labels < lower)
    score = score + (2.0 / alpha) * (labels - upper) * (labels > upper)
    return torch.mean(score)


@dataclass(frozen=True)
class SplitConformalCalibrator:
    """A finite-sample split-conformal symmetric interval calibrator."""

    radius: float
    alpha: float
    normalized: bool = False
    tie_policy: str = "higher"

    def __post_init__(self) -> None:
        _probability(self.alpha)
        if not isinstance(self.normalized, bool):
            raise TypeError("normalized must be a boolean")
        if self.tie_policy != "higher":
            raise ValueError("tie_policy must be 'higher'")
        if isinstance(self.radius, bool) or not isinstance(self.radius, (int, float)):
            raise TypeError("radius must be a real number")
        if math.isnan(self.radius) or self.radius < 0.0:
            raise ValueError("radius must be nonnegative")

    @classmethod
    def fit(
        cls,
        mean: Tensor,
        labels: Tensor,
        *,
        alpha: float,
        standard_deviation: Tensor | None = None,
        normalized: bool = False,
    ) -> SplitConformalCalibrator:
        alpha = _probability(alpha)
        _vector(mean, "mean")
        _vector(labels, "labels")
        _same_vectors(mean=mean, labels=labels)
        if not isinstance(normalized, bool):
            raise TypeError("normalized must be a boolean")
        scores = (labels - mean).abs()
        if normalized:
            if standard_deviation is None:
                raise ValueError("normalized conformal calibration requires standard_deviation")
            _vector(standard_deviation, "standard_deviation", positive=True)
            _same_vectors(mean=mean, labels=labels, standard_deviation=standard_deviation)
            scores = scores / standard_deviation
        elif standard_deviation is not None:
            raise ValueError("standard_deviation is only accepted when normalized=True")

        rank = math.ceil((scores.numel() + 1) * (1.0 - alpha))
        radius = math.inf if rank > scores.numel() else float(torch.sort(scores).values[rank - 1])
        return cls(radius=radius, alpha=alpha, normalized=normalized)

    def interval(
        self,
        mean: Tensor,
        standard_deviation: Tensor | None = None,
    ) -> tuple[Tensor, Tensor]:
        _vector(mean, "mean")
        if self.normalized:
            if standard_deviation is None:
                raise ValueError("normalized conformal intervals require standard_deviation")
            _vector(standard_deviation, "standard_deviation", positive=True)
            _same_vectors(mean=mean, standard_deviation=standard_deviation)
            width = standard_deviation * self.radius
        else:
            if standard_deviation is not None:
                raise ValueError("standard_deviation is only accepted for normalized intervals")
            width = torch.full_like(mean, self.radius)
        return mean - width, mean + width

    def to_dict(self) -> dict[str, Any]:
        return {
            "alpha": self.alpha,
            "normalized": self.normalized,
            "radius": None if math.isinf(self.radius) else self.radius,
            "tie_policy": self.tie_policy,
        }


def conformal_calibrator_from_dict(payload: Mapping[str, Any]) -> SplitConformalCalibrator:
    if not isinstance(payload, Mapping) or set(payload) != {
        "alpha", "normalized", "radius", "tie_policy",
    }:
        raise ValueError("invalid conformal calibrator payload")
    radius = payload["radius"]
    if radius is None:
        radius = math.inf
    return SplitConformalCalibrator(
        alpha=payload["alpha"],
        normalized=payload["normalized"],
        radius=radius,
        tie_policy=payload["tie_policy"],
    )


def deterministic_group_folds(
    groups: Sequence[str | int],
    *,
    n_splits: int,
    seed: int = 0,
) -> tuple[int, ...]:
    """Assign rows to balanced folds without splitting a group."""
    if isinstance(n_splits, bool) or not isinstance(n_splits, int):
        raise TypeError("n_splits must be an integer")
    if n_splits < 2:
        raise ValueError("n_splits must be at least two")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer")
    if isinstance(groups, (str, bytes)) or not isinstance(groups, Sequence):
        raise TypeError("groups must be a sequence")
    if not groups:
        raise ValueError("groups must not be empty")
    rows: dict[tuple[str, str], list[int]] = {}
    for index, group in enumerate(groups):
        if isinstance(group, bool) or not isinstance(group, (str, int)):
            raise TypeError("group identifiers must be strings or integers")
        key = (type(group).__name__, str(group))
        rows.setdefault(key, []).append(index)
    if n_splits > len(rows):
        raise ValueError("n_splits cannot exceed the number of distinct groups")

    ordered = sorted(
        rows.items(),
        key=lambda item: (
            -len(item[1]),
            hashlib.sha256(f"{seed}:{item[0][0]}:{item[0][1]}".encode("utf-8")).hexdigest(),
        ),
    )
    counts = [0] * n_splits
    assignment = [-1] * len(groups)
    for _, indices in ordered:
        fold = min(range(n_splits), key=lambda candidate: (counts[candidate], candidate))
        for index in indices:
            assignment[index] = fold
        counts[fold] += len(indices)
    return tuple(assignment)


def grouped_fold_indices(
    groups: Sequence[str | int],
    *,
    n_splits: int,
    seed: int = 0,
) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """Return deterministic ``(train, validation)`` indices for grouped folds."""
    assignment = deterministic_group_folds(groups, n_splits=n_splits, seed=seed)
    return tuple(
        (
            tuple(index for index, value in enumerate(assignment) if value != fold),
            tuple(index for index, value in enumerate(assignment) if value == fold),
        )
        for fold in range(n_splits)
    )


def _fit_gaussian_calibrator(
    kind: str, mean: Tensor, standard_deviation: Tensor, labels: Tensor
) -> GaussianCalibrator:
    calibrator_types = {
        "identity": IdentityGaussianCalibrator,
        "scalar": ScalarGaussianCalibrator,
        "affine_variance": AffineVarianceGaussianCalibrator,
    }
    if kind not in calibrator_types:
        raise ValueError(f"unsupported calibration model: {kind}")
    return calibrator_types[kind].fit(mean, standard_deviation, labels)


@dataclass(frozen=True)
class GroupedCalibrationSelection:
    """Validation-only grouped cross-validation selection result."""

    calibrator: GaussianCalibrator
    selected_kind: str
    scores: tuple[dict[str, Any], ...]
    fold_assignment: tuple[int, ...]
    n_splits: int
    seed: int
    nll_tolerance: float
    crps_tolerance: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_order_simplest_first": list(_MODEL_ORDER),
            "crps_tolerance": self.crps_tolerance,
            "fold_assignment": list(self.fold_assignment),
            "grouping": "frozen_mmseqs_group_id",
            "n_splits": self.n_splits,
            "nll_tolerance": self.nll_tolerance,
            "primary_metric": "gaussian_nll",
            "scores": [dict(score) for score in self.scores],
            "secondary_metric": "gaussian_crps",
            "seed": self.seed,
            "selected_kind": self.selected_kind,
            "tie_rule": "nll_then_crps_then_simpler_model_within_tolerances",
        }


def select_grouped_gaussian_calibrator(
    mean: Tensor,
    standard_deviation: Tensor,
    labels: Tensor,
    groups: Sequence[str | int],
    *,
    n_splits: int = 5,
    seed: int = 0,
    nll_tolerance: float = 1e-12,
    crps_tolerance: float = 1e-12,
) -> GroupedCalibrationSelection:
    """Select variance calibration by grouped CV on one validation partition."""
    _validated_gaussian_sample(mean, standard_deviation, labels)
    if len(groups) != labels.numel():
        raise ValueError("groups must have one identifier per validation observation")
    for value, name in ((nll_tolerance, "nll_tolerance"), (crps_tolerance, "crps_tolerance")):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number")
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")

    assignment = deterministic_group_folds(groups, n_splits=n_splits, seed=seed)
    folds = grouped_fold_indices(groups, n_splits=n_splits, seed=seed)
    scores = []
    for kind in _MODEL_ORDER:
        nll_sum = 0.0
        crps_sum = 0.0
        count = 0
        for train_indices, validation_indices in folds:
            train = torch.tensor(train_indices, dtype=torch.long, device=mean.device)
            validation = torch.tensor(validation_indices, dtype=torch.long, device=mean.device)
            calibrator = _fit_gaussian_calibrator(
                kind, mean[train], standard_deviation[train], labels[train]
            )
            fold_mean, fold_deviation = calibrator.calibrate(
                mean[validation], standard_deviation[validation]
            )
            fold_count = len(validation_indices)
            nll_sum += float(gaussian_nll(fold_mean, fold_deviation, labels[validation])) * fold_count
            crps_sum += float(gaussian_crps(fold_mean, fold_deviation, labels[validation])) * fold_count
            count += fold_count
        scores.append({"kind": kind, "gaussian_nll": nll_sum / count,
                       "gaussian_crps": crps_sum / count})

    selected = scores[0]
    for candidate in scores[1:]:
        nll_delta = candidate["gaussian_nll"] - selected["gaussian_nll"]
        if nll_delta < -nll_tolerance:
            selected = candidate
            continue
        if abs(nll_delta) <= nll_tolerance:
            crps_delta = candidate["gaussian_crps"] - selected["gaussian_crps"]
            if crps_delta < -crps_tolerance:
                selected = candidate
    calibrator = _fit_gaussian_calibrator(
        selected["kind"], mean, standard_deviation, labels
    )
    return GroupedCalibrationSelection(
        calibrator=calibrator,
        selected_kind=selected["kind"],
        scores=tuple(scores),
        fold_assignment=assignment,
        n_splits=n_splits,
        seed=seed,
        nll_tolerance=float(nll_tolerance),
        crps_tolerance=float(crps_tolerance),
    )


def calibrator_from_dict(payload: Mapping[str, Any]) -> GaussianCalibrator:
    if not isinstance(payload, Mapping):
        raise TypeError("calibrator payload must be a mapping")
    kind = payload.get("kind")
    expected_keys = {
        "identity": {"kind"},
        "scalar": {"kind", "scale"},
        "affine_variance": {"kind", "scale", "offset"},
    }
    if kind not in expected_keys or set(payload) != expected_keys[kind]:
        raise ValueError("invalid calibrator payload")
    if kind == "identity":
        return IdentityGaussianCalibrator()
    if kind == "scalar":
        return ScalarGaussianCalibrator(payload["scale"])
    return AffineVarianceGaussianCalibrator(payload["scale"], payload["offset"])


def _canonical_json(payload: Any) -> bytes:
    try:
        encoded = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError("artifact values must be finite JSON values") from exc
    return encoded


def save_calibration_artifact(
    path: str | os.PathLike[str],
    calibrator: GaussianCalibrator,
    *,
    identities: Mapping[str, Any],
    validation_count: int,
    conformal_calibrators: Sequence[SplitConformalCalibrator] = (),
    selection: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Serialize a calibrator bound to exact input/checkpoint identities."""
    if not isinstance(
        calibrator,
        (
            IdentityGaussianCalibrator,
            ScalarGaussianCalibrator,
            AffineVarianceGaussianCalibrator,
        ),
    ):
        raise TypeError("calibrator has an unsupported type")
    if not isinstance(identities, Mapping) or not identities:
        raise ValueError("identities must be a non-empty mapping")
    if isinstance(validation_count, bool) or not isinstance(validation_count, int):
        raise TypeError("validation_count must be an integer")
    if validation_count < 1:
        raise ValueError("validation_count must be positive")
    conformal = {}
    for interval in conformal_calibrators:
        if not isinstance(interval, SplitConformalCalibrator):
            raise TypeError("conformal_calibrators contains an unsupported type")
        key = f"{int(round((1.0 - interval.alpha) * 100))}_{'normalized' if interval.normalized else 'unnormalized'}"
        if key in conformal:
            raise ValueError(f"duplicate conformal interval: {key}")
        conformal[key] = interval.to_dict()
    extended = bool(conformal) or selection is not None
    body = {
        "format_version": _FORMAT_VERSION if extended else 1,
        "calibrator": calibrator.to_dict(),
        "identities": dict(identities),
        "fit_partition": "validation",
        "validation_count": validation_count,
    }
    if extended:
        if not conformal:
            raise ValueError("version 2 artifacts require conformal calibrators")
        body["conformal_intervals"] = conformal
        body["selection"] = None if selection is None else dict(selection)
    digest = hashlib.sha256(_canonical_json(body)).hexdigest()
    payload = {**body, "artifact_sha256": digest}
    destination = Path(path)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return payload


def load_calibration_artifact(
    path: str | os.PathLike[str],
    *,
    expected_identities: Mapping[str, Any] | None = None,
) -> tuple[GaussianCalibrator, dict[str, Any]]:
    """Load an artifact after integrity and optional exact identity checks."""
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    required = {
        "format_version",
        "calibrator",
        "identities",
        "fit_partition",
        "validation_count",
        "artifact_sha256",
    }
    if isinstance(payload, dict) and payload.get("format_version") == 2:
        required |= {"conformal_intervals", "selection"}
    if not isinstance(payload, dict) or set(payload) != required:
        raise ValueError("invalid calibration artifact schema")
    digest = payload.pop("artifact_sha256")
    actual = hashlib.sha256(_canonical_json(payload)).hexdigest()
    payload["artifact_sha256"] = digest
    if not isinstance(digest, str) or not hmac.compare_digest(digest, actual):
        raise ValueError("calibration artifact integrity check failed")
    if payload["format_version"] not in (1, _FORMAT_VERSION) or payload["fit_partition"] != "validation":
        raise ValueError("unsupported calibration artifact")
    if not isinstance(payload["identities"], dict) or not payload["identities"]:
        raise ValueError("calibration artifact identities must be a non-empty mapping")
    validation_count = payload["validation_count"]
    if isinstance(validation_count, bool) or not isinstance(validation_count, int):
        raise ValueError("calibration artifact validation_count must be an integer")
    if validation_count < 1:
        raise ValueError("calibration artifact validation_count must be positive")
    if expected_identities is not None:
        if not isinstance(expected_identities, Mapping):
            raise TypeError("expected_identities must be a mapping")
        if payload["identities"] != dict(expected_identities):
            raise ValueError("calibration artifact identity mismatch")
    calibrator = calibrator_from_dict(payload["calibrator"])
    if payload["format_version"] == 2:
        intervals = payload["conformal_intervals"]
        if not isinstance(intervals, dict) or not intervals:
            raise ValueError("calibration artifact conformal_intervals must be non-empty")
        for interval in intervals.values():
            conformal_calibrator_from_dict(interval)
        if payload["selection"] is not None and not isinstance(payload["selection"], dict):
            raise ValueError("calibration artifact selection must be an object or null")
    return calibrator, payload
