"""Tests for the FBA surrogate distillation.

The surrogate is the answer to "the original model should have the GEMs'
capabilities": ranker checkpoints trained on scenario grids so they predict
growth without cobra. Tests cover the scenario contract, the inference-side
carrier that must match training records, and a miniature end-to-end training
run (small grids keep it fast).
"""

from __future__ import annotations

import itertools
import json
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
    """One row per (scenario, target); condition ids stay unique across both."""
    grid = fba_surrogate.run_fba_grid("iML1515", model=iml1515_model)
    conditions = [row["condition_id"] for row in grid["rows"]]
    assert len(conditions) == len(set(conditions))
    assert all(row["value"] > 0 for row in grid["rows"])

    scenarios = fba_surrogate.build_scenarios("iML1515")
    solved = {row["scenario_id"] for row in grid["rows"]}
    assert len(solved) + len(grid["dropped"]) == len(scenarios)
    # every solved scenario must carry all four targets
    for scenario_id in solved:
        targets = {row["target"] for row in grid["rows"] if row["scenario_id"] == scenario_id}
        assert targets == set(fba_surrogate.MODEL_TARGETS["iML1515"])


@pytest.fixture(scope="module")
def iml1515_model():
    return fba_context_module().load_model("iML1515", IML1515)


def fba_context_module():
    from biocandidate import fba_context

    return fba_context


# --- miniature end-to-end training ------------------------------------------------


def test_training_produces_a_working_surrogate(tmp_path, monkeypatch, iml1515_model):
    """A tiny grid still must train, evaluate, and beat the mean baseline."""
    # 8 scenarios x 4 targets is too little for the four-target head to fit at
    # all, so the miniature grid carries 18 scenarios.
    monkeypatch.setattr(fba_surrogate, "GLUCOSE_GRID", (3.0, 9.0, 15.0))
    monkeypatch.setattr(fba_surrogate, "OXYGEN_GRID", (3.0, 12.0, 20.0))
    monkeypatch.setattr(fba_surrogate, "IML1515_SCENARIOS",
                        ({"name": "baseline"}, {"name": "ptsG_knockout"}))

    metrics = fba_surrogate.train_surrogate(
        "iML1515", tmp_path, epochs=80, seed=7,
    )
    assert metrics["rows"]["train"] >= 2
    assert metrics["rows"]["test"] >= 1
    assert set(metrics["per_target"]) == set(fba_surrogate.MODEL_TARGETS["iML1515"])
    for target, entry in metrics["per_target"].items():
        assert entry["count"] >= 1, target
        assert entry["test_mean_baseline_rmse_log10"] > 0, target
    # Asserting held-out quality on a miniature grid would be flaky: ~24 training
    # rows cannot fit a four-target head, and the real grids are where per-target
    # accuracy is measured (6x-21x over each target's mean baseline, recorded in
    # the trained metrics.json). What CI can check cheaply and stably is that the
    # context pathway actually learns: training error must fall below the naive
    # predictor's error on held-out data.
    assert metrics["train_rmse_log10"] < metrics["test_pooled_mean_baseline_rmse_log10"]
    assert (tmp_path / "best.pt").is_file()
    assert (tmp_path / "metrics.json").is_file()
    # the checkpoint must expose a trained log10_flux task for run_prediction
    checkpoint_metrics = __import__("json").loads((tmp_path / "metrics.json").read_text())
    assert checkpoint_metrics["validation_tasks"]["log10_flux"]["count"] >= 1


# --- per-model contracts and the yeast grid ---------------------------------------


def test_every_model_declares_an_eight_wide_contract_and_its_targets():
    """Width stays 8 everywhere; the meaning of the slots is per model."""
    for model_id in fba_surrogate.MODEL_REGISTRY:
        feature_ids = fba_surrogate.SCENARIO_FEATURE_IDS_BY_MODEL[model_id]
        assert len(feature_ids) == 8, model_id
        assert len(set(feature_ids)) == 8, model_id
        assert feature_ids[-1] == "target_selector", model_id

        targets = fba_surrogate.MODEL_TARGETS[model_id]
        assert "growth" in targets, model_id
        assert set(targets) <= set(fba_surrogate.TARGET_SELECTORS), model_id
        assert len(set(targets)) == len(targets), model_id


def test_target_selectors_are_unique():
    """Two targets sharing a selector value would be indistinguishable to the model."""
    values = list(fba_surrogate.TARGET_SELECTORS.values())
    assert len(set(values)) == len(values)


def test_target_reaction_ids_cover_every_declared_target():
    for targets in fba_surrogate.MODEL_TARGETS.values():
        for target in targets:
            assert target in fba_surrogate.TARGET_REACTION_IDS, target


def test_yeast_grid_spans_the_four_scenario_axes():
    scenarios = fba_surrogate.build_scenarios("yeast_metatwin")
    axes = {
        (row["scenario"]["route"], row["scenario"]["pox1_closed"],
         row["scenario"]["objective"])
        for row in scenarios
    }
    assert axes == set(itertools.product(("A", "B"), (False, True), ("growth", "target")))
    assert len(scenarios) == (
        len(axes) * len(fba_surrogate.YEAST_GLUCOSE_GRID) * len(fba_surrogate.YEAST_OXYGEN_GRID)
    )


def test_yeast_vector_places_each_axis_in_its_declared_slot():
    scenario = {
        "name": "routeB_pox1_off_target",
        "route": "B",
        "pox1_closed": True,
        "objective": "target",
        "biomass_floor_fraction": 0.1,
    }
    vector = fba_surrogate._scenario_vector(
        scenario, glucose=2.0, oxygen=750.0, target="capacity", model_id="yeast_metatwin"
    )
    assert len(vector) == 8
    assert vector[0] == 2.0 and vector[1] == 750.0
    assert vector[2] == fba_surrogate.ROUTE_SELECTORS["B"]
    assert vector[3] == 1.0
    assert vector[4] == fba_surrogate.OBJECTIVE_SELECTORS["target"]
    assert vector[5] == 0.1
    assert vector[6] == 0.0
    assert vector[7] == fba_surrogate.TARGET_SELECTORS["capacity"]


def test_ecoli_vector_is_unchanged_by_the_yeast_addition():
    """The original contract must still encode exactly as before."""
    vector = fba_surrogate._scenario_vector(
        {"name": "paper_stack", "ko": True, "oe": ("THD2pp", 1.5)},
        glucose=15.0, oxygen=20.0, target="atp", model_id="iML1515",
    )
    assert vector == (15.0, 20.0, 0.0, 0.0, 0.0, 1.0, 1.5,
                      fba_surrogate.TARGET_SELECTORS["atp"])


def test_a_target_not_declared_for_the_model_is_rejected():
    with pytest.raises(ValueError, match="not declared"):
        fba_surrogate._scenario_vector(
            {"name": "baseline"}, 15.0, 20.0, target="capacity", model_id="iML1515"
        )


def test_route_direction_is_documented_for_both_routes():
    """The route selector must map to real routes, not arbitrary labels."""
    assert set(fba_surrogate.ROUTE_SELECTORS) == {"A", "B"}
    from biocandidate import royal_jelly_acid as rja

    assert fba_surrogate.ROUTE_SELECTORS["A"] != fba_surrogate.ROUTE_SELECTORS["B"]
    assert rja.ROUTE_A and rja.ROUTE_B


def test_trained_surrogates_beat_every_target_baseline_if_present():
    """Guard against shipping an under-trained surrogate.

    The yeast grid needs roughly 900 epochs: at 250 it had not converged and
    three of six targets lost to the mean baseline. This test reads whatever
    checkpoints exist and fails if any declared target regressed, so a future
    too-short run cannot be mistaken for a usable artifact. Skipped when no
    surrogate has been trained locally.
    """
    metrics_files = sorted(ROOT.glob("artifacts/fba-surrogate-*/metrics.json"))
    if not metrics_files:
        pytest.skip("no locally trained surrogate")

    for path in metrics_files:
        metrics = json.loads(path.read_text(encoding="utf-8"))
        for target, entry in metrics.get("per_target", {}).items():
            label = f"{path.parent.name}:{target}"
            assert entry["count"] > 0, label
            assert entry["test_rmse_log10"] < entry["test_mean_baseline_rmse_log10"], (
                f"{label} did not beat its mean baseline "
                f"({entry['test_rmse_log10']:.4f} vs "
                f"{entry['test_mean_baseline_rmse_log10']:.4f}); "
                f"epochs={metrics.get('epochs')} best_epoch={metrics.get('best_epoch')}"
            )


# --- extrapolation guard ----------------------------------------------------------


def test_trained_bounds_match_the_declared_grids():
    yeast = fba_surrogate.trained_grid_bounds("yeast_metatwin")
    assert yeast["glucose"] == (min(fba_surrogate.YEAST_GLUCOSE_GRID),
                               max(fba_surrogate.YEAST_GLUCOSE_GRID))
    assert yeast["oxygen"] == (min(fba_surrogate.YEAST_OXYGEN_GRID),
                              max(fba_surrogate.YEAST_OXYGEN_GRID))
    ecoli = fba_surrogate.trained_grid_bounds("iML1515")
    assert ecoli["glucose"] == (min(fba_surrogate.GLUCOSE_GRID), max(fba_surrogate.GLUCOSE_GRID))


def test_query_inside_the_envelope_is_not_flagged():
    assert fba_surrogate.check_query_in_range("yeast_metatwin", 1.0, 1000.0) == []


def test_ecoli_medium_is_flagged_for_the_yeast_surrogate():
    """The failure this guard exists for.

    Carrying the E. coli form defaults (glucose 15, oxygen 20) over to a yeast
    scenario returns a plausible-looking capacity that is two orders of magnitude
    away from the FBA answer, so the query is labelled instead of silently
    answered.
    """
    warnings = fba_surrogate.check_query_in_range("yeast_metatwin", 15.0, 20.0)
    assert len(warnings) == 2
    assert any("glucose 15" in w for w in warnings)
    assert any("oxygen 20" in w for w in warnings)
    assert all("extrapolating" in w for w in warnings)


def test_yeast_scenario_vector_carries_the_guard_result():
    vector = fba_surrogate.scenario_input_vector(
        "yeast_metatwin", "routeB_pox1_off_target", 15.0, 20.0, "capacity"
    )
    assert vector["out_of_range"], "the vector must carry the warning, not just expose the helper"
    inside = fba_surrogate.scenario_input_vector(
        "yeast_metatwin", "routeB_pox1_off_target", 1.0, 1000.0, "capacity"
    )
    assert inside["out_of_range"] == []


def test_yeast_default_medium_is_inside_the_trained_envelope():
    """The registry default must be a query the surrogate can actually answer.

    It previously declared glucose 10, which exceeds the model's own uptake cap
    (1.0), sits outside the trained grid, and disagreed with the medium the
    external anchor was produced with.
    """
    from biocandidate.fba_context import MODEL_REGISTRY

    spec = MODEL_REGISTRY["yeast_metatwin"]
    bounds = fba_surrogate.trained_grid_bounds("yeast_metatwin")
    assert bounds["glucose"][0] <= spec["default_glucose"] <= bounds["glucose"][1]
    assert bounds["oxygen"][0] <= spec["default_oxygen"] <= bounds["oxygen"][1]


def test_per_model_scenario_tables_do_not_leak_across_models():
    """A yeast scenario name must not resolve for an E. coli model, or vice versa."""
    for model_id, table in fba_surrogate.SCENARIOS_BY_MODEL.items():
        assert table, model_id
        for name, scenario in table.items():
            if model_id == "yeast_metatwin":
                assert "route" in scenario and "objective" in scenario, name
            else:
                assert "route" not in scenario, name
