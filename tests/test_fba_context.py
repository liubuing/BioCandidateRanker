"""Tests for the genome-scale FBA context fusion.

External stacks are fused here — the pantothenate paper's iML1515 workflow,
the arginine project's iJO1366 fusion model, and the yeast pathway-prediction
project's Yeast-MetaTwin model — so the tests cover the seams: the registry
contract, the paper's modification semantics, the synthetic fusion layer, the
yeast registration, and (when a checkpoint is present) the shape the ranker
consumes.
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
YEAST = ROOT / "data/Yeast-MetaTwin.yml"
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


@pytest.fixture(scope="module")
def yeast():
    if not YEAST.is_file():
        pytest.skip("yeast model file required")
    return fba_context.load_model("yeast_metatwin", YEAST)


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
    """All registered models emit exactly 8 feature ids with a per-model objective."""
    assert set(fba_context.model_ids()) == {"iML1515", "iJO1366", "yeast_metatwin"}
    for model_id, spec in fba_context.MODEL_REGISTRY.items():
        assert len(spec["feature_ids"]) == 8
        assert spec["feature_ids"][0] == f"flux:{spec['objective_id']}"
        assert len(set(spec["feature_ids"])) == 8
        assert (ROOT / fba_context.MODEL_CONFIGS[model_id]).is_file()


def test_fusion_presets_are_iJO1366_only():
    assert all("fusion" not in p for p in fba_context.presets_for("iML1515"))
    assert "fusion_M5_NOG_PtxD" in fba_context.presets_for("iJO1366")
    assert fba_context.presets_for("yeast_metatwin") == ["baseline"]


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


# --- the yeast registration ------------------------------------------------------


def test_yeast_baseline_is_optimal_with_eight_features(yeast):
    result = fba_context.run_fba_context(
        model_id="yeast_metatwin", preset="baseline", model=yeast
    )
    assert result["status"] == "optimal"
    assert len(result["features"]) == 8
    assert result["metadata"]["model_id"] == "yeast_metatwin"
    assert result["metadata"]["objective_id"] == "r_2111"
    # Declared registration defaults: glucose at the model file's own uptake cap,
    # oxygen at the file's cap. This was 10.0, a generic batch rate the model
    # cannot represent, which also put the default outside the yeast surrogate's
    # trained envelope and disagreed with the medium the external 10H2DA anchor
    # (0.274011) was produced with.
    assert result["fluxes"]["r_1714"] == pytest.approx(-1.0)
    assert result["growth"] == pytest.approx(0.0895, abs=1e-3)
    # the POX1 precursor gate that feeds the 10H2DA reference chain is in the vector
    assert "flux:r_0120" in result["metadata"]["feature_ids"]


def test_yeast_condition_id_is_deterministic(yeast):
    first = fba_context.run_fba_context(
        model_id="yeast_metatwin", preset="baseline", model=yeast
    )
    second = fba_context.run_fba_context(
        model_id="yeast_metatwin", preset="baseline", model=yeast
    )
    assert first["condition_id"] == second["condition_id"]
    assert first["features"] == second["features"]


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


# --- cofactor accounting ---------------------------------------------------------


def test_registry_declares_cofactor_targets_for_every_model():
    """Each model must declare its own cofactor metabolites, core trio included."""
    for model_id, spec in fba_context.MODEL_REGISTRY.items():
        targets = spec.get("cofactor_metabolites")
        assert targets, f"{model_id} declares no cofactor targets"
        assert set(fba_context.CORE_COFACTORS) <= set(targets), model_id
        # the extended set is what makes the registry reusable beyond the paper's
        # three cofactors, so require at least the flavin and thioester families
        assert {"fad", "fadh2", "coa", "accoa"} <= set(targets), model_id


@pytest.mark.parametrize("model_id", ["iML1515", "iJO1366", "yeast_metatwin"])
def test_cofactor_ids_exist_in_their_model(model_id):
    model = fba_context.load_model(model_id)
    for label, metabolite_id in fba_context.MODEL_REGISTRY[model_id]["cofactor_metabolites"].items():
        assert metabolite_id in model.metabolites, f"{model_id} lacks {label} ({metabolite_id})"


def test_cofactor_turnover_is_self_consistent_and_steady(iml1515):
    """Gross turnover must equal the sum of its reported producers, and net must vanish.

    Net internal production is zero by construction at steady state, so it
    doubles as a solution sanity check; a non-zero value means the solution is
    not balanced for that cofactor.
    """
    result = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    turnover = result["cofactor_turnover"]
    declared = set(fba_context.MODEL_REGISTRY["iML1515"]["cofactor_metabolites"])
    assert set(turnover) == declared
    assert set(fba_context.CORE_COFACTORS) <= set(turnover)

    for label, info in turnover.items():
        assert info["status"] == "reported"
        assert info["unit"] == "mmol gDW^-1 h^-1"
        # steady state must hold for every declared cofactor
        assert abs(info["net_production"]) < 1e-9, f"{label} is not at steady state"
        assert info["gross_turnover"] >= 0
        if info["gross_turnover"] > 0:
            assert info["top_producers"]
            # the reporting window must not exceed the total it summarises
            reported = sum(item["flux"] for item in info["top_producers"])
            assert reported <= info["gross_turnover"] + 1e-9
            assert all(item["flux"] > 0 for item in info["top_producers"])
        else:
            assert info["top_producers"] == []

    # the core trio must carry flux at any feasible baseline; the extended set
    # need not (GTP turnover is legitimately zero at this condition)
    for label in fba_context.CORE_COFACTORS:
        assert turnover[label]["gross_turnover"] > 0, label
    assert turnover["gtp"]["gross_turnover"] == 0.0


def test_nadph_and_atp_supply_match_the_paper_levers(iml1515):
    """The dominant producers must be the pathways the paper engineers."""
    turnover = fba_context.run_fba_context(
        model_id="iML1515", preset="baseline", model=iml1515
    )["cofactor_turnover"]

    nadph_producers = {item["reaction_id"] for item in turnover["nadph"]["top_producers"]}
    # pentose phosphate shunt: the paper's NADPH regeneration target
    assert "G6PDH2r" in nadph_producers or "GND" in nadph_producers
    atp_producers = {item["reaction_id"] for item in turnover["atp"]["top_producers"]}
    # oxidative phosphorylation: the paper's electron-transport-chain target
    assert "ATPS4rpp" in atp_producers


def test_cofactor_reporting_does_not_widen_the_feature_vector(iml1515):
    """The 8-wide contract every shipped checkpoint depends on must not change."""
    result = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    assert len(result["features"]) == 8
    assert len(result["metadata"]["feature_ids"]) == 8
    assert "cofactor_turnover" not in result["metadata"]
    assert result["cofactor_scope"].startswith("gross turnover")


def test_paper_stack_does_not_move_the_fba_optimum(iml1515):
    """Documented model behaviour, not a bug.

    At the paper's 15/20 medium the ptsG knockout is fully compensated by the
    alternative glucose uptake route, and the two overexpression targets sit at
    zero or near-zero baseline flux, so the LP optimum - and therefore the
    cofactor turnover - is unchanged. This is precisely why the paper had to
    measure intracellular pools rather than rely on FBA: the model-side
    quantity cannot see the effect the engineered strain exhibits.
    """
    baseline = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    stack = fba_context.run_fba_context(model_id="iML1515", preset="paper_stack", model=iml1515)

    assert stack["condition_id"] != baseline["condition_id"]
    assert stack["applied_modifications"] == [
        "knockout:b1101",
        "overexpression:THD2pp@1.5",
        "overexpression:NADK@1.5",
    ]
    # Not bit-identical: the NADK bound (1.5x a 0.0004 baseline flux) does move
    # the optimum, but by less than 2e-6, i.e. far below any reported precision.
    assert stack["growth"] == pytest.approx(baseline["growth"], abs=1e-4)
    assert abs(stack["growth"] - baseline["growth"]) < 1e-5
    for label in ("nadph", "nadh", "atp"):
        assert stack["cofactor_turnover"][label]["gross_turnover"] == pytest.approx(
            baseline["cofactor_turnover"][label]["gross_turnover"], abs=1e-3
        )


def test_extended_cofactors_carry_real_flux(iml1515):
    """A declared cofactor nobody uses would be a silently useless column.

    FAD and acetyl-CoA are the two the royal-jelly-acid route analysis depends
    on, so require that a real reaction carries them at baseline and that the
    dominant producer is identifiable.
    """
    turnover = fba_context.run_fba_context(
        model_id="iML1515", preset="baseline", model=iml1515
    )["cofactor_turnover"]

    for label in ("fad", "accoa", "malcoa", "mlthf"):
        info = turnover[label]
        assert info["status"] == "reported", label
        assert info["gross_turnover"] > 0, f"{label} carries no flux at baseline"
        assert info["top_producers"], label
    # the dominant FAD producer should be a flavin-dependent oxidation step
    assert turnover["fad"]["top_producers"][0]["reaction_id"] == "ACOAD1fr"


def test_yeast_flavin_flux_is_zero_at_baseline(yeast):
    """Documented model behaviour: no flavin-dependent reaction carries flux here.

    Kept as a test because a zero is a real result, not a reporting failure, and
    it is the reason POX1 activity in this model must be read from the r_0120
    reaction flux rather than from FAD turnover.
    """
    result = fba_context.run_fba_context(
        model_id="yeast_metatwin", preset="baseline", model=yeast
    )
    turnover = result["cofactor_turnover"]
    assert turnover["fad"]["status"] == "reported"
    # Not exactly zero: a residual 9.24e-06 flows through r_0440, which is 1e-5 of
    # the growth rate and negligible against any flux in the route analysis.
    assert turnover["fad"]["gross_turnover"] < 1e-4
    assert turnover["fad"]["gross_turnover"] / result["growth"] < 1e-4
    # while the thioester supply for fatty-acid synthesis does carry flux
    assert turnover["malcoa"]["gross_turnover"] > 0


def test_every_model_declares_oxidant_handles():
    """A model missing these silently drops the degradation-side reporting.

    This happened during development: a duplicated key left one model without an
    oxidant map, and the route analysis quietly lost its oxygen/peroxide section
    for that model instead of failing.
    """
    for model_id, spec in fba_context.MODEL_REGISTRY.items():
        oxidants = spec.get("oxidant_metabolites")
        assert oxidants, f"{model_id} declares no oxidant handles"
        assert set(oxidants) == {"h2o2", "o2"}, model_id


@pytest.mark.parametrize("model_id", ["iML1515", "iJO1366", "yeast_metatwin"])
def test_oxidant_ids_exist_in_their_model(model_id):
    model = fba_context.load_model(model_id)
    for label, metabolite_id in fba_context.MODEL_REGISTRY[model_id]["oxidant_metabolites"].items():
        assert metabolite_id in model.metabolites, f"{model_id} lacks {label} ({metabolite_id})"


def test_oxidative_burden_is_reported_with_gross_directions(iml1515):
    """Oxygen is both taken up and consumed internally, so both directions matter."""
    result = fba_context.run_fba_context(
        model_id="iML1515", preset="baseline", model=iml1515
    )
    burden = result["oxidative_burden"]
    assert set(burden) == {"h2o2", "o2"}
    oxygen = burden["o2"]
    assert oxygen["status"] == "reported"
    assert oxygen["consumed"] > 0
    assert oxygen["top_consumers"], "the respiratory consumer must be identifiable"
    assert burden["h2o2"]["status"] == "reported"


# --- marginal cost -----------------------------------------------------------------


def test_least_squares_recovers_a_known_line():
    slope, intercept, r2 = fba_context._least_squares([0.0, 1.0, 2.0], [1.0, 3.0, 5.0])
    assert slope == pytest.approx(2.0)
    assert intercept == pytest.approx(1.0)
    assert r2 == pytest.approx(1.0)


def test_least_squares_reports_a_poor_fit_honestly():
    """A curved response must not be summarised as a clean slope."""
    _, _, r2 = fba_context._least_squares([0.0, 1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 10.0])
    assert r2 < 0.9


def test_marginal_cost_needs_a_maximisable_target(iml1515):
    result = fba_context.marginal_cofactor_cost(
        iml1515,
        "EX_glc__D_e",  # a lower-bounded exchange, not a producible flux
        cofactor_metabolites={"nadph": "nadph_c"},
        level_fractions=(0.5, 1.0),
    )
    assert result["status"] != "reported"


def test_marginal_nadph_cost_is_positive_and_recorded(iml1515):
    """The core claim of the analysis: making more product costs more NADPH."""
    result = fba_context.marginal_cofactor_cost(
        iml1515,
        "BIOMASS_Ec_iML1515_core_75p37M",
        cofactor_metabolites={"nadph": "nadph_c", "atp": "atp_c"},
        objective_reaction="ATPM",
        level_fractions=(0.25, 0.5, 0.75, 1.0),
    )
    assert result["status"] == "reported"
    assert len(result["points"]) == 4
    nadph = result["marginal"]["nadph"]
    assert nadph["slope_per_unit_product"] is not None
    # levels are ascending and every point carries the tracked quantities
    levels = [point["level"] for point in result["points"]]
    assert levels == sorted(levels)
    for point in result["points"]:
        assert set(point["cofactors"]) == {"nadph", "atp"}
        assert all(value is not None for value in point["cofactors"].values())


def test_marginal_ranking_is_sorted_by_cost(iml1515):
    result = fba_context.marginal_cofactor_cost(
        iml1515,
        "BIOMASS_Ec_iML1515_core_75p37M",
        cofactor_metabolites={"nadph": "nadph_c", "atp": "atp_c"},
        objective_reaction="ATPM",
        level_fractions=(0.25, 1.0),
    )
    slopes = [item["slope_per_unit_product"] for item in result["ranking_by_cost"]]
    assert slopes == sorted(slopes, reverse=True)


def test_marginal_objective_is_reported_when_named(iml1515):
    result = fba_context.marginal_cofactor_cost(
        iml1515,
        "BIOMASS_Ec_iML1515_core_75p37M",
        cofactor_metabolites={"nadph": "nadph_c"},
        objective_reaction="ATPM",
        level_fractions=(0.25, 1.0),
    )
    assert result["objective"] is not None
    assert result["objective"]["slope_per_unit_product"] is not None
