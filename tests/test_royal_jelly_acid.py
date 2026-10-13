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
