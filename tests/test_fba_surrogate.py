"""Tests for the FBA surrogate distillation.

The surrogate is the answer to "the original model should have the GEMs'
capabilities": ranker checkpoints trained on scenario grids so they predict
growth without cobra. Tests cover the scenario contract, the inference-side
carrier that must match training records, and a miniature end-to-end training
run (small grids keep it fast).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("cobra")

from biocandidate import fba_surrogate  # noqa: E402

IML1515 = ROOT / "data/iML1515.xml"
IJO1366 = ROOT / "data/iJO1366.json"
pytestmark = pytest.mark.skipif(
    not (IML1515.is_file() and IJO1366.is_file()),
    reason="model files under the gitignored data/ tree are required",
)


# --- scenario contract -----------------------------------------------------------


def test_scenario_features_match_the_checkpoint_contract():
    assert len(fba_surrogate.SCENARIO_FEATURE_IDS) == 8
    assert len(set(fba_surrogate.SCENARIO_FEATURE_IDS)) == 8


def test_paper_scenario_vector_encodes_knockout_and_overexpression():
    vector = fba_surrogate.scenario_input_vector(
        "iML1515", "paper_stack", glucose=15.0, oxygen=20.0
    )
    assert vector["features"] == [15.0, 20.0, 0.0, 0.0, 0.0, 1.0, 1.5, 0.0]
    assert vector["metadata"]["model_id"] == "iML1515"
    assert vector["metadata"]["feature_ids"] == list(fba_surrogate.SCENARIO_FEATURE_IDS)


def test_scenario_vector_is_deterministic():
    first = fba_surrogate.scenario_input_vector("iML1515", "baseline", 12.0, 9.0)
    second = fba_surrogate.scenario_input_vector("iML1515", "baseline", 12.0, 9.0)
    assert first["features"] == second["features"]
    assert first["metadata"]["condition_id"] == second["metadata"]["condition_id"]


def test_fusion_scenario_requires_iJO1366():
    with pytest.raises(ValueError, match="iJO1366"):
        fba_surrogate.scenario_input_vector("iML1515", "fusion_M5_NOG_PtxD", 10.0, 20.0)
    ok = fba_surrogate.scenario_input_vector("iJO1366", "fusion_M5_NOG_PtxD", 10.0, 20.0)
    assert ok["features"][2] == 3.0  # NOG cap
    assert ok["features"][3] == 5.0  # phosphite cap
    assert ok["features"][4] == 5.0  # CO2 cap


def test_unknown_scenario_is_rejected():
    with pytest.raises(ValueError, match="unknown scenario"):
        fba_surrogate.scenario_input_vector("iML1515", "not_a_preset", 10.0, 20.0)


# --- grid sweep ------------------------------------------------------------------


def test_grid_sweep_covers_the_paper_operating_point():
    scenarios = fba_surrogate.build_scenarios("iML1515")
    combos = {(row["glucose"], row["oxygen"]) for row in scenarios}
    assert (15.0, 20.0) in combos, "the paper's medium must be inside the training grid"


def test_grid_rows_are_unique_and_solved(iml1515_model):
    grid = fba_surrogate.run_fba_grid("iML1515", model=iml1515_model)
    conditions = [row["condition_id"] for row in grid["rows"]]
    assert len(conditions) == len(set(conditions))
    assert all(row["growth"] > 0 for row in grid["rows"])
    assert len(grid["rows"]) + len(grid["dropped"]) == len(fba_surrogate.build_scenarios("iML1515"))


@pytest.fixture(scope="module")
def iml1515_model():
    return fba_context_module().load_model("iML1515", IML1515)


def fba_context_module():
    from biocandidate import fba_context

    return fba_context


# --- miniature end-to-end training ------------------------------------------------


def test_training_produces_a_working_surrogate(tmp_path, monkeypatch, iml1515_model):
    """A tiny grid still must train, evaluate, and beat the mean baseline."""
    monkeypatch.setattr(fba_surrogate, "GLUCOSE_GRID", (6.0, 15.0))
    monkeypatch.setattr(fba_surrogate, "OXYGEN_GRID", (5.0, 20.0))
    monkeypatch.setattr(fba_surrogate, "IML1515_SCENARIOS",
                        ({"name": "baseline"}, {"name": "ptsG_knockout"}))

    metrics = fba_surrogate.train_surrogate(
        "iML1515", tmp_path, epochs=30, seed=7,
    )
    assert metrics["rows"]["train"] >= 2
    assert metrics["rows"]["test"] >= 1
    assert metrics["test_rmse_log10"] < metrics["test_mean_baseline_rmse_log10"]
    assert (tmp_path / "best.pt").is_file()
    assert (tmp_path / "metrics.json").is_file()
    # the checkpoint must expose a trained log10_flux task for run_prediction
    checkpoint_metrics = __import__("json").loads((tmp_path / "metrics.json").read_text())
    assert checkpoint_metrics["validation_tasks"]["log10_flux"]["count"] >= 1
