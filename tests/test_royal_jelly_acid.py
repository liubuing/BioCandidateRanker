"""Tests for the royal jelly acid route analysis.

The anchor is the whole point of this module: the external project's Sep 2
report states a maximum target flux of 0.274011 on the model-file default
medium, and this module must reproduce it before any of its own numbers are
trusted. Skipped when the model file or the external target configuration is
unavailable.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("cobra")

from biocandidate import royal_jelly_acid as rja  # noqa: E402

MODEL = ROOT / "data/Yeast-MetaTwin.yml"
TARGET_CONFIG = Path(
    "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"
    "/10_generic_target_workflow/examples/target_workflow_10h2da_reference.json"
)

requires_model = pytest.mark.skipif(
    not MODEL.is_file(), reason="the yeast model under gitignored data/ is required"
)
requires_target = pytest.mark.skipif(
    not TARGET_CONFIG.is_file(), reason="the external target configuration is required"
)


@pytest.fixture(scope="module")
def report():
    return rja.analyze(TARGET_CONFIG)


# --- the anchor ---------------------------------------------------------------


@requires_model
@requires_target
def test_reproduces_the_documented_anchor(report):
    """The single most important guard: the external anchor must come back."""
    check = report["anchor_check"]
    assert check["matches"] is True, check
    assert check["measured"] == pytest.approx(rja.DOCUMENTED_ANCHOR, abs=1e-5)


# --- route structure ----------------------------------------------------------


@requires_model
@requires_target
def test_both_documented_routes_are_scored(report):
    assert set(report["routes"]) == {"route_A_free_acid", "route_B_coa_thioester"}
    for name, entry in report["routes"].items():
        assert entry["status"] == "optimal", name
        assert entry["target_capacity"] > 0
        assert entry["route_fluxes"], name


@requires_model
@requires_target
def test_routes_are_equivalent_before_kcat_constraints(report):
    """At the LP level neither route is capacity-limited, so both agree.

    This is asserted rather than assumed because it is the reason the two routes
    only separate once enzyme capacity binds - which is what the pilot measures.
    """
    capacities = [entry["target_capacity"] for entry in report["routes"].values()]
    assert capacities[0] == pytest.approx(capacities[1], rel=1e-9)


@requires_model
@requires_target
def test_capacity_ratio_is_unit_for_both_routes(report):
    for entry in report["routes"].values():
        assert entry["capacity_ratio_to_route_a"] == pytest.approx(1.0, rel=1e-9)


# --- the POX1 prediction ------------------------------------------------------


@requires_model
@requires_target
def test_closing_pox1_lowers_rather_than_raises_capacity(report):
    """A falsifiable, counter-intuitive model prediction, kept explicit.

    The naive expectation is that blocking the degradation entrance raises target
    accumulation. In this model at this medium it does the opposite, because the
    reaction carries redox balance the biomass floor requires.
    """
    effect = report["pox1_effect"]
    assert effect
    for name, entry in effect.items():
        assert entry["closed"] < entry["open"], name
        assert -0.05 < entry["change_fraction"] < 0, name
    assert "opposite of the naive" in report["pox1_interpretation"]


@requires_model
@requires_target
def test_pox1_is_open_in_the_default_configuration(report):
    """The anchor was produced with POX1 open, so the default must not close it."""
    assert report["pox1_closed"] is False


# --- the cofactor view --------------------------------------------------------


@requires_model
@requires_target
def test_feasibility_cofactors_are_reported_per_route(report):
    for name, entry in report["routes"].items():
        turnover = entry["cofactor_turnover"]
        assert set(turnover) == set(rja.FEASIBILITY_COFACTORS), name
        for label in ("nadph", "atp", "accoa"):
            info = turnover[label]
            assert info["status"] == "reported", (name, label)
            # the P450 step needs NADPH and the precursor chain needs acetyl-CoA,
            # so both must carry flux on any feasible route
            assert info["gross_turnover"] > 0, (name, label)
        # the core trio agrees with the FBA context reporting
        assert set(entry["core_cofactor_turnover"]) == {"nadph", "nadh", "atp"}


# --- units --------------------------------------------------------------------


def test_minimum_enzyme_follows_the_gecko_bound():
    """ub = kcat * E * 3600, so E = flux / (kcat * 3600)."""
    assert rja.minimum_enzyme_for_flux(1.0, 3600.0) == pytest.approx(1.0)
    assert rja.minimum_enzyme_for_flux(2.0, 3600.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        rja.minimum_enzyme_for_flux(0.0, 1.0)


# --- oxidative burden: the degradation-side handle flavins cannot give ------------


def test_flavin_dependent_oxidases_carry_no_flavin_in_the_model():
    """Why FAD turnover cannot stand in for oxidase activity.

    POX1 and its five paralogues are written as acyl-CoA + O2 -> enoyl-CoA +
    H2O2. That is correct chemistry for a flavin-dependent oxidase: the FAD
    cycles inside the reaction as a catalyst, so the net reaction contains no
    flavin and FAD turnover would be zero regardless of how hard the enzyme
    works. The handles that do move with oxidase flux are O2 and H2O2.
    """
    from biocandidate.fba_context import load_model
    from biocandidate.royal_jelly_acid import POX1_REACTION

    model = load_model("yeast_metatwin")
    reaction = model.reactions.get_by_id(POX1_REACTION)
    names = {(metabolite.name or "").lower() for metabolite in reaction.metabolites}
    assert not any("fad" in name for name in names), names
    assert "oxygen" in names and "hydrogen peroxide" in names, names

    # the same is true of the whole acyl-CoA oxidase family, not just POX1
    oxidases = [r for r in model.reactions if "acyl-coa oxidase" in r.name.lower()]
    assert len(oxidases) >= 5, [r.id for r in oxidases]
    for oxidase in oxidases:
        assert not any(
            "fad" in (metabolite.name or "").lower() for metabolite in oxidase.metabolites
        ), oxidase.id


@requires_model
@requires_target
def test_route_analysis_reports_oxygen_and_peroxide(report):
    for name, entry in report["routes"].items():
        burden = entry["oxidative_burden"]
        assert set(burden) == {"o2", "h2o2"}, name
        assert burden["o2"]["consumed"] > 0, name
        assert burden["h2o2"]["produced"] > 0, name
        assert burden["h2o2"]["produced"] <= burden["o2"]["consumed"] + 1e-9, name


@requires_model
@requires_target
def test_pox1_and_the_omega_hydroxylase_dominate_the_route_burden(report):
    """The two route-specific oxygen users must be identifiable, not averaged away."""
    entry = report["routes"]["route_A_free_acid"]
    capacity = entry["target_capacity"]

    o2_consumers = {
        item["reaction_id"]: item["flux"] for item in entry["oxidative_burden"]["o2"]["top_consumers"]
    }
    # the P450 omega-hydroxylation consumes exactly one O2 per product
    assert o2_consumers["CAND_T2DEC_OMEGA_HYDROXYLASE_P"] == pytest.approx(capacity, rel=1e-6)
    # POX1 is the other route-specific oxygen user
    assert o2_consumers["r_0120"] == pytest.approx(capacity, rel=1e-6)

    producers = [item["reaction_id"] for item in entry["oxidative_burden"]["h2o2"]["top_producers"]]
    assert producers[0] == "r_0120", producers
    # peroxide is produced by the oxidase family, so the total exceeds POX1 alone
    assert entry["oxidative_burden"]["h2o2"]["produced"] > capacity
