"""Tests for the genome-scale FBA context fusion.

Two external stacks are fused here — the pantothenate paper's iML1515 workflow
and the arginine project's iJO1366 fusion model — so the tests cover the seams:
the registry contract, the paper's modification semantics, the synthetic fusion
layer, and (when a checkpoint is present) the shape the ranker consumes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

pytest.importorskip("cobra")

from biocandidate import fba_context  # noqa: E402

IML1515 = ROOT / "data/iML1515.xml"
IJO1366 = ROOT / "data/iJO1366.json"
pytestmark = pytest.mark.skipif(
    not (IML1515.is_file() and IJO1366.is_file()),
    reason="model files under the gitignored data/ tree are required",
)


@pytest.fixture(scope="module")
def iml1515():
    return fba_context.load_model("iML1515", IML1515)


@pytest.fixture(scope="module")
def ijo1366():
    return fba_context.load_model("iJO1366", IJO1366)


# --- pure logic ------------------------------------------------------------------


@pytest.mark.parametrize(
    "value,expected",
    [(0.0, 0.0), (1.0, pytest.approx(0.3010, abs=1e-4)),
     (-1.0, pytest.approx(-0.3010, abs=1e-4)), (9.0, pytest.approx(1.0, abs=1e-4))],
)
def test_signed_log10_is_symmetric_and_zero_preserving(value, expected):
    assert fba_context.signed_log10(value) == expected


def test_signed_log10_rejects_non_finite():
    with pytest.raises(ValueError, match="finite"):
        fba_context.signed_log10(float("inf"))


def test_registry_contract():
    """Both registered models emit exactly 8 feature ids with a per-model objective."""
    assert set(fba_context.model_ids()) == {"iML1515", "iJO1366"}
    for model_id, spec in fba_context.MODEL_REGISTRY.items():
        assert len(spec["feature_ids"]) == 8
        assert spec["feature_ids"][0] == f"flux:{spec['objective_id']}"
        assert len(set(spec["feature_ids"])) == 8
        assert (ROOT / fba_context.MODEL_CONFIGS[model_id]).is_file()


def test_fusion_presets_are_iJO1366_only():
    assert all("fusion" not in p for p in fba_context.presets_for("iML1515"))
    assert "fusion_M5_NOG_PtxD" in fba_context.presets_for("iJO1366")


def test_unknown_model_and_preset_are_rejected(iml1515):
    with pytest.raises(ValueError, match="unknown model"):
        fba_context.run_fba_context(model_id="ecoli_core", model=iml1515)
    with pytest.raises(ValueError, match="unknown preset"):
        fba_context.run_fba_context(model_id="iML1515", preset="fusion_M2_NOG", model=iml1515)


def test_model_file_identity_is_enforced(tmp_path):
    bogus = tmp_path / "iML1515.xml"
    bogus.write_text("<not-the-model/>", encoding="utf-8")
    with pytest.raises(ValueError, match="identity mismatch"):
        fba_context.load_model("iML1515", bogus)


# --- paper semantics on iML1515 --------------------------------------------------


def test_iml1515_baseline_is_optimal_with_eight_features(iml1515):
    result = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    assert result["status"] == "optimal"
    assert result["growth"] == pytest.approx(0.9945, abs=0.01)
    assert len(result["features"]) == 8
    assert result["metadata"]["schema_version"] == 1
    assert result["metadata"]["model_id"] == "iML1515"
    # the paper's medium: glucose 15, oxygen 20
    assert result["fluxes"]["EX_glc__D_e"] == pytest.approx(-15.0)


def test_paper_medium_matches_the_manuscript(iml1515):
    """The paper sets glucose 15 and oxygen 20 mmol/gdw/h; defaults mirror it."""
    assert fba_context.DEFAULT_GLUCOSE == 15.0
    assert fba_context.DEFAULT_OXYGEN == 20.0


def test_knockout_changes_the_condition_id(iml1515):
    """Delta-ptsG is a real scenario from the paper's strain series.

    FBA fact worth documenting: at 15 mmol/gDW/h glucose the PTS knockout is
    fully compensated by the alternative uptake route, so the optimal flux
    distribution over the eight feature reactions is unchanged - which is why
    the paper's own DPAW strains pair the knockout with glk/glf overexpression.
    The scenario identity still differs, so the condition id must too.
    """
    baseline = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    knockout = fba_context.run_fba_context(
        model_id="iML1515", preset="ptsG_knockout", model=iml1515
    )
    assert knockout["condition_id"] != baseline["condition_id"]
    assert knockout["applied_modifications"] == ["knockout:b1101"]
    assert knockout["features"] == baseline["features"]
    assert knockout["growth"] == pytest.approx(baseline["growth"], abs=1e-6)


def test_overexpression_tightens_the_target_bound(iml1515):
    result = fba_context.run_fba_context(
        model_id="iML1515", preset="thd2pp_overexpression", model=iml1515
    )
    assert result["status"] == "optimal"
    assert result["applied_modifications"] == ["overexpression:THD2pp@1.5"]


def test_same_condition_is_deterministic(iml1515):
    first = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    second = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    assert first["condition_id"] == second["condition_id"]
    assert first["features"] == second["features"]


# --- the arginine fusion layer on iJO1366 ----------------------------------------


def test_fusion_layer_adds_eight_synthetic_reactions(ijo1366):

    before = len(ijo1366.reactions)
    with ijo1366 as scoped:
        added = fba_context.add_fusion_pathways(scoped)
        assert len(scoped.reactions) == before + 8
        assert set(added) == {
            "CCR_EMA", "YciA_EMA", "EMAtex", "EX_emaacid_e",
            "NOG_F6P", "PHITEtex", "EX_phite_e", "PTXD",
        }
    # the context manager reverts the layer: the base model is untouched
    assert len(ijo1366.reactions) == before


def test_fusion_scenario_runs_and_caps_are_applied(ijo1366):
    result = fba_context.run_fba_context(
        model_id="iJO1366", preset="fusion_M5_NOG_PtxD", model=ijo1366
    )
    assert result["status"] == "optimal"
    assert result["fusion_pathways_applied"] is True
    assert len(result["features"]) == 8
    assert "融合FBA" in result["provenance"]


def test_fusion_scenario_bound_scan_is_monotone(ijo1366):
    """A growth floor constrains the optimum from below without pinning it.

    This also pins the objective contract: the registry declares the WT biomass
    reaction (following the 05_FBA script), and run_fba_context must set it
    explicitly instead of trusting the model file's default "core" objective.
    """
    free = fba_context.run_fba_context(
        model_id="iJO1366", preset="fusion_M5_NOG_PtxD", model=ijo1366
    )
    floored = fba_context.run_fba_context(
        model_id="iJO1366", preset="fusion_M5_NOG_PtxD", growth_min=1.29, model=ijo1366,
    )
    assert floored["status"] == "optimal"
    assert floored["growth"] >= 1.29 - 1e-6
    # a floor cannot raise the optimum; compare with solver-noise tolerance
    assert free["growth"] >= floored["growth"] - 1e-6
    # objective identity: the WT biomass is the maximized reaction
    assert floored["metadata"]["objective_id"] == "BIOMASS_Ec_iJO1366_WT_53p95M"


# --- the ranker seam -------------------------------------------------------------


def test_fba_features_match_checkpoint_width(iml1515):
    """8 features must line up with the frozen checkpoints' fba_context_dim."""
    import torch

    from biocandidate.training import load_checkpoint

    reference = ROOT / "artifacts/esm2-t6-opt-seed7/best.pt"
    if not reference.is_file():
        pytest.skip("local checkpoint required")
    with _safe_globals():
        model, _ = load_checkpoint(reference, torch.device("cpu"))
    result = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    assert len(result["features"]) == model.config.fba_context_dim
    assert len(result["metadata"]["feature_ids"]) == model.config.fba_context_dim


def _safe_globals():
    from serve_workbench import _safe_load_globals

    return _safe_load_globals()
