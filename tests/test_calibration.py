import json
import math

import pytest
import torch

from biocandidate.calibration import (
    AffineVarianceGaussianCalibrator,
    IdentityGaussianCalibrator,
    ScalarGaussianCalibrator,
    SplitConformalCalibrator,
    conformal_calibrator_from_dict,
    deterministic_group_folds,
    gaussian_crps,
    gaussian_nll,
    grouped_fold_indices,
    interval_score,
    load_calibration_artifact,
    save_calibration_artifact,
    select_grouped_gaussian_calibrator,
)


def test_gaussian_scores_and_interval_score_formulas():
    mean = torch.tensor([0.0], dtype=torch.float64)
    deviation = torch.tensor([2.0], dtype=torch.float64)
    label = torch.tensor([0.0], dtype=torch.float64)
    assert gaussian_nll(mean, deviation, label).item() == pytest.approx(
        math.log(2.0) + 0.5 * math.log(2.0 * math.pi)
    )
    assert gaussian_crps(mean, deviation, label).item() == pytest.approx(
        2.0 * (math.sqrt(2.0) - 1.0) / math.sqrt(math.pi)
    )
    assert interval_score(
        torch.tensor([-1.0]), torch.tensor([1.0]), torch.tensor([2.0]), alpha=0.1
    ).item() == pytest.approx(22.0)


def test_identity_scalar_and_affine_variance_formulas():
    mean = torch.tensor([1.0, 2.0])
    deviation = torch.tensor([2.0, 3.0])
    identity_mean, identity_sd = IdentityGaussianCalibrator().calibrate(mean, deviation)
    assert identity_mean is mean
    assert identity_sd is deviation

    scalar = ScalarGaussianCalibrator.fit(mean, deviation, torch.tensor([3.0, 8.0]))
    assert scalar.scale == pytest.approx(math.sqrt(2.5))
    assert torch.equal(scalar.calibrate(mean, deviation)[1], deviation * scalar.scale)

    affine = AffineVarianceGaussianCalibrator(scale=2.0, offset=3.0)
    assert torch.allclose(
        affine.calibrate(mean, deviation)[1], torch.sqrt(2.0 * deviation.square() + 3.0)
    )


def test_affine_variance_fit_improves_validation_nll():
    mean = torch.zeros(6, dtype=torch.float64)
    deviation = torch.tensor([0.5, 0.5, 1.0, 1.0, 2.0, 2.0], dtype=torch.float64)
    labels = torch.tensor([-1.0, 1.0, -1.5, 1.5, -2.5, 2.5], dtype=torch.float64)
    fitted = AffineVarianceGaussianCalibrator.fit(mean, deviation, labels)
    _, calibrated_deviation = fitted.calibrate(mean, deviation)
    assert gaussian_nll(mean, calibrated_deviation, labels) < gaussian_nll(mean, deviation, labels)


@pytest.mark.parametrize(
    "operation",
    [
        lambda: gaussian_nll(torch.tensor([]), torch.tensor([]), torch.tensor([])),
        lambda: gaussian_crps(torch.tensor([0]), torch.tensor([1]), torch.tensor([0])),
        lambda: ScalarGaussianCalibrator.fit(
            torch.tensor([0.0]), torch.tensor([0.0]), torch.tensor([1.0])
        ),
        lambda: interval_score(
            torch.tensor([1.0]), torch.tensor([0.0]), torch.tensor([0.5]), alpha=0.1
        ),
        lambda: SplitConformalCalibrator.fit(torch.tensor([0.0]), torch.tensor([1.0]), alpha=1.0),
    ],
)
def test_invalid_inputs_are_rejected(operation):
    with pytest.raises((TypeError, ValueError)):
        operation()


def test_conformal_finite_sample_rank_normalization_and_tie_policy():
    mean = torch.zeros(4)
    labels = torch.tensor([1.0, 1.0, 2.0, 9.0])
    conformal = SplitConformalCalibrator.fit(mean, labels, alpha=0.4)
    # ceil((4 + 1) * .6) = 3, with the conservative higher order statistic.
    assert conformal.radius == 2.0
    lower, upper = conformal.interval(torch.tensor([2.0]))
    assert lower.item() == 0.0
    assert upper.item() == 4.0
    assert bool((labels.abs() <= conformal.radius)[2])  # Boundary ties are covered.

    normalized = SplitConformalCalibrator.fit(
        mean,
        torch.tensor([1.0, 2.0, 6.0, 8.0]),
        alpha=0.4,
        standard_deviation=torch.tensor([1.0, 2.0, 2.0, 2.0]),
        normalized=True,
    )
    assert normalized.radius == 3.0
    assert normalized.interval(torch.tensor([10.0]), torch.tensor([2.0]))[0].item() == 4.0

    finite_sample = SplitConformalCalibrator.fit(mean, labels, alpha=0.1)
    assert math.isinf(finite_sample.radius)


def test_infinite_radius_round_trips_through_json_artifact_without_crashing():
    mean = torch.tensor([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], dtype=torch.float64)
    labels = torch.tensor([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], dtype=torch.float64)
    infinite_radius = SplitConformalCalibrator.fit(mean, labels, alpha=0.05)
    assert math.isinf(infinite_radius.radius)
    assert infinite_radius.to_dict()["radius"] is None

    serialized = json.dumps(infinite_radius.to_dict(), allow_nan=False)
    restored = conformal_calibrator_from_dict(json.loads(serialized))
    assert math.isinf(restored.radius)
    lower, upper = restored.interval(torch.tensor([4.0]))
    assert not math.isfinite(float(lower[0]))
    assert not math.isfinite(float(upper[0]))


def test_group_folds_are_deterministic_and_isolate_groups():
    groups = ["a", "a", "b", "c", "c", "d"]
    first = deterministic_group_folds(groups, n_splits=3, seed=17)
    assert first == deterministic_group_folds(groups, n_splits=3, seed=17)
    for group in set(groups):
        assert len({first[index] for index, value in enumerate(groups) if value == group}) == 1
    for train, validation in grouped_fold_indices(groups, n_splits=3, seed=17):
        assert set(train).isdisjoint(validation)
        assert {groups[index] for index in train}.isdisjoint(groups[index] for index in validation)


def test_grouped_selection_scores_only_held_out_groups_and_is_deterministic():
    mean = torch.zeros(8, dtype=torch.float64)
    deviation = torch.ones(8, dtype=torch.float64)
    labels = torch.tensor([-2.0, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0],
                          dtype=torch.float64)
    groups = ["a", "a", "b", "b", "c", "c", "d", "d"]
    result = select_grouped_gaussian_calibrator(
        mean, deviation, labels, groups, n_splits=2, seed=9,
    )
    repeat = select_grouped_gaussian_calibrator(
        mean, deviation, labels, groups, n_splits=2, seed=9,
    )
    assert result.selected_kind == "scalar"
    assert result.calibrator == ScalarGaussianCalibrator(2.0)
    assert result.to_dict() == repeat.to_dict()
    for group in set(groups):
        assert len({result.fold_assignment[i] for i, value in enumerate(groups)
                    if value == group}) == 1


def test_grouped_selection_uses_deterministic_simpler_model_tie_rule():
    mean = torch.zeros(6, dtype=torch.float64)
    deviation = torch.ones(6, dtype=torch.float64)
    labels = torch.tensor([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0], dtype=torch.float64)
    result = select_grouped_gaussian_calibrator(
        mean, deviation, labels, ["a", "a", "b", "b", "c", "c"], n_splits=3,
        nll_tolerance=1e-6, crps_tolerance=1e-6,
    )
    assert result.selected_kind == "identity"


def test_artifact_round_trip_integrity_and_identity(tmp_path):
    path = tmp_path / "calibration.json"
    identities = {
        "source": {"sha256": "abc", "row_count": 3},
        "split": {"sha256": "def"},
        "checkpoint": "run-7",
    }
    calibrator = ScalarGaussianCalibrator(1.25)
    saved = save_calibration_artifact(path, calibrator, identities=identities, validation_count=3)
    loaded, payload = load_calibration_artifact(path, expected_identities=identities)
    assert loaded == calibrator
    assert payload == saved
    assert payload["fit_partition"] == "validation"

    with pytest.raises(ValueError, match="identity mismatch"):
        load_calibration_artifact(path, expected_identities={"source": "other"})

    tampered = json.loads(path.read_text(encoding="utf-8"))
    tampered["calibrator"]["scale"] = 9.0
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(ValueError, match="integrity"):
        load_calibration_artifact(path)


def test_v2_artifact_round_trip_includes_deployable_conformal_intervals(tmp_path):
    path = tmp_path / "calibration-v2.json"
    intervals = tuple(
        SplitConformalCalibrator(radius=radius, alpha=alpha, normalized=normalized)
        for alpha, radius in ((0.1, 1.5), (0.05, 2.0))
        for normalized in (False, True)
    )
    saved = save_calibration_artifact(
        path,
        ScalarGaussianCalibrator(1.25),
        identities={"checkpoint": "abc", "task": "log10_kcat"},
        validation_count=20,
        conformal_calibrators=intervals,
        selection={"selected_kind": "scalar"},
    )
    calibrator, loaded = load_calibration_artifact(path)
    assert calibrator == ScalarGaussianCalibrator(1.25)
    assert loaded == saved
    assert loaded["format_version"] == 2
    assert set(loaded["conformal_intervals"]) == {
        "90_normalized", "90_unnormalized", "95_normalized", "95_unnormalized",
    }
