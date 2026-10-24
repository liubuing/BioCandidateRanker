"""Tests for the prosthetic-group ledger.

The ledger exists because of a real contradiction: the currency ledger says
yeast FAD turnover is ~1e-5 of growth, and the experimental side says FAD is
needed for the route. Both are right - one speaks for flux-borne cofactors,
the other for enzyme-bound ones - and these tests pin the bridge: the bill's
physics, the content table that encodes it, the FMN floor that survives the
model's degenerate kinase vertex, and the burden ratios that formalise the
experimental observation as a number.

中文：辅基账本测试：调和"通量口径 FAD 周转近零"与"实验需要 FAD"的真实矛盾——一个说通量携带辅因子，一个说酶结合辅因子。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("cobra")

from biocandidate import prosthetic_ledger as pl  # noqa: E402

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


# --- the bill's physics -------------------------------------------------------


def test_enzyme_load_is_flux_over_kcat_times_an_hour():
    assert pl.enzyme_load(0.274, 1.0) == pytest.approx(0.274 / 3600.0)
    assert pl.enzyme_load(0.274, 10.0) == pytest.approx(pl.enzyme_load(0.274, 1.0) / 10.0)
    assert pl.enzyme_load(0.0, 1.0) == 0.0


def test_enzyme_load_rejects_unphysical_inputs():
    with pytest.raises(ValueError, match="kcat"):
        pl.enzyme_load(1.0, 0.0)
    with pytest.raises(ValueError, match="non-negative"):
        pl.enzyme_load(-1.0, 1.0)


# --- the content table --------------------------------------------------------


def test_the_reductase_carries_the_flavins_and_the_p450_carries_the_heme():
    """The biochemistry the whole bridge rests on, stated as data."""
    assert pl.ENZYME_PROSTHETIC_CONTENT["CYP52M1"] == {"heme_b": 1.0}
    assert pl.ENZYME_PROSTHETIC_CONTENT["NCP1"] == {"fad": 1.0, "fmn": 1.0}
    assert pl.ENZYME_PROSTHETIC_CONTENT["POX1"] == {"fad": 1.0}
    assert pl.ENZYME_PROSTHETIC_CONTENT["TES1"] == {}


def test_the_bill_sums_inventory_and_charges_replacement():
    bill = pl.ledger_for_fluxes(
        {"CYP52M1": 0.274, "TES1": 0.274, "POX1": 0.274},
        kcat_s=1.0,
        protein_half_life_h=10.0,
        extra_enzymes={"NCP1": "CYP52M1"},
    )
    load = 0.274 / 3600.0
    # NCP1 rides at the P450's load, so two FAD carriers plus the P450's heme
    assert bill["enzyme_loads"]["NCP1"] == pytest.approx(load)
    assert bill["group_inventory_mmol_gdw"]["fad"] == pytest.approx(2.0 * load)
    assert bill["group_inventory_mmol_gdw"]["fmn"] == pytest.approx(load)
    assert bill["group_inventory_mmol_gdw"]["heme_b"] == pytest.approx(load)
    for group, inventory in bill["group_inventory_mmol_gdw"].items():
        assert bill["group_replacement_flux_mmol_gdw_h"][group] == pytest.approx(
            inventory / 10.0
        )
    assert bill["enzymes_without_prosthetics"] == ["TES1"]


def test_the_bill_refuses_unknown_enzymes_groups_and_halflives():
    with pytest.raises(KeyError, match="no prosthetic content"):
        pl.ledger_for_fluxes({"NOT_AN_ENZYME": 1.0}, kcat_s=1.0, protein_half_life_h=10.0)
    with pytest.raises(ValueError, match="half-life"):
        pl.ledger_for_fluxes({"TES1": 1.0}, kcat_s=1.0, protein_half_life_h=0.0)
    from unittest.mock import patch

    with patch.dict(
        pl.ENZYME_PROSTHETIC_CONTENT, {"MYSTERY": {"unobtainium": 1.0}}
    ):
        with pytest.raises(KeyError, match="synthesis reaction"):
            pl.ledger_for_fluxes({"MYSTERY": 1.0}, kcat_s=1.0, protein_half_life_h=10.0)


# --- the denominators ---------------------------------------------------------


def test_burden_ratios_report_a_missing_baseline_instead_of_dividing_by_zero():
    bill = pl.ledger_for_fluxes(
        {"CYP52M1": 0.274}, kcat_s=1.0, protein_half_life_h=10.0,
        extra_enzymes={"NCP1": "CYP52M1"},
    )
    bill["growth_at_operating_point"] = 0.01
    load = 0.274 / 3600.0
    capacity = {
        "fad": {"flux_mmol_gdw_h": 8.95e-7, "per_growth_intensity_mmol_gdw": 1.0e-5},
        "fmn": {
            "flux_mmol_gdw_h": 8.95e-7,
            "per_growth_intensity_mmol_gdw": 1.0e-5,
            "floor_per_growth_intensity_mmol_gdw": 1.0e-5,
        },
        "heme_b": {
            "flux_mmol_gdw_h": 0.0,
            "per_growth_intensity_mmol_gdw": 0.0,
        },
    }
    ratios = pl.burden_ratios(bill, capacity)
    row = ratios["heme_b"]
    assert row["replacement_vs_endogenous_flux"] is None
    assert row["inventory_vs_biomass_inventory"] is None
    assert row["no_endogenous_baseline"] is True
    assert "new load" in row["note"]
    # the flavins do have a baseline, so they price normally
    assert ratios["fad"]["replacement_vs_endogenous_flux"] == pytest.approx(
        load / 10.0 / 8.95e-7
    )


def test_the_fmn_ratio_uses_the_one_to_one_floor_not_the_degenerate_vertex():
    """The kinase can sit at zero in some optimal vertices; the floor cannot.

    Every FAD the cell makes consumed one FMN first, so even a zero-kinase
    vertex still owes FMN synthesis at the FAD intensity. The ratio must be
    computed against that floor, or the model's bookkeeping degeneracy would
    make the route's FMN demand look infinitely worse than it is.
    """
    bill = pl.ledger_for_fluxes(
        {"CYP52M1": 0.274}, kcat_s=1.0, protein_half_life_h=10.0,
        extra_enzymes={"NCP1": "CYP52M1"},
    )
    bill["growth_at_operating_point"] = 0.01
    capacity = {
        "fad": {"flux_mmol_gdw_h": 8.95e-7, "per_growth_intensity_mmol_gdw": 1.0e-5},
        "fmn": {
            "flux_mmol_gdw_h": 0.0,
            "per_growth_intensity_mmol_gdw": 0.0,
            "floor_per_growth_intensity_mmol_gdw": 1.0e-5,
        },
        "heme_b": {
            "flux_mmol_gdw_h": 8.95e-8,
            "per_growth_intensity_mmol_gdw": 1.0e-6,
        },
    }
    row = pl.burden_ratios(bill, capacity)["fmn"]
    replacement = bill["group_replacement_flux_mmol_gdw_h"]["fmn"]
    assert row["replacement_vs_endogenous_flux"] == pytest.approx(
        replacement / (1.0e-5 * 0.01)
    )


# --- the route bill -----------------------------------------------------------


@pytest.fixture(scope="module")
def route_report():
    return pl.route_ledger(TARGET_CONFIG)


@requires_model
@requires_target
def test_the_capacity_table_reports_three_groups_with_a_stable_fad_intensity(
    route_report,
):
    capacity = route_report["endogenous_synthesis_capacity"]
    assert capacity["fad"]["flux_mmol_gdw_h"] > 0
    assert capacity["fad"]["per_growth_intensity_mmol_gdw"] == pytest.approx(
        1.0e-5, rel=0.01
    )
    assert capacity["fmn"]["floor_per_growth_intensity_mmol_gdw"] == pytest.approx(
        capacity["fad"]["per_growth_intensity_mmol_gdw"]
    )
    assert "recycles FMN" in capacity["fmn"]["floor_note"]


@requires_model
@requires_target
def test_the_grid_covers_every_kcat_and_halflife_combination(route_report):
    assert len(route_report["scenarios"]) == 6
    for scenario in route_report["scenarios"].values():
        assert set(scenario["burden_ratios"]) == {"fad", "fmn", "heme_b"}
        assert "claim_boundary" in route_report


@requires_model
@requires_target
def test_the_fad_burden_falls_monotonically_with_kcat_and_halflife(route_report):
    """Two sensitivity axes, both one-directional by construction."""
    def ratio(kcat: str, half_life: str) -> float:
        scenario = route_report["scenarios"][f"kcat_{kcat}_half_life_{half_life}h"]
        return scenario["burden_ratios"]["fad"]["replacement_vs_endogenous_flux"]

    assert ratio("0.1", "10") > ratio("1", "10") > ratio("10", "10")
    assert ratio("1", "10") > ratio("1", "72")


@requires_model
@requires_target
def test_the_fad_demand_exceeds_the_housekeeping_budget_across_the_plausible_range(
    route_report,
):
    """The formalised experimental observation.

    P450s turn over slowly; across every kcat at or below 1/s and both
    half-lives, the route's FAD replacement demand is a multiple of what the
    cell already spends on its own flavins. Only the most favourable corner -
    a fast enzyme surviving a long time - drops below one. This is the bridge
    between "the currency ledger says zero" and "the experiment says needed".
    """
    load_bearing = {
        name: scenario["load_bearing"]["fad"]
        for name, scenario in route_report["scenarios"].items()
    }
    assert load_bearing["kcat_0.1_half_life_10h"] is True
    assert load_bearing["kcat_1_half_life_10h"] is True
    assert load_bearing["kcat_1_half_life_72h"] is True
    assert load_bearing["kcat_0.1_half_life_72h"] is True
    assert load_bearing["kcat_10_half_life_72h"] is False
    assert (
        route_report["scenarios"]["kcat_1_half_life_10h"]["burden_ratios"]["fad"][
            "replacement_vs_endogenous_flux"
        ]
        > 2.0
    )


@requires_model
@requires_target
def test_the_route_fluxes_are_the_capacity_and_pox1_is_measured(route_report):
    assignment = route_report["route_flux_assignment"]
    assert assignment["CYP52M1"] == pytest.approx(route_report["route_capacity"])
    assert assignment["TES1"] == pytest.approx(route_report["route_capacity"])
    assert assignment["NCP1 (coupled to CYP52M1)"] == pytest.approx(
        route_report["route_capacity"]
    )
    assert assignment["POX1"] >= 0.0
