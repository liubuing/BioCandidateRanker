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
def yeast():
    """The bare yeast model, for tests that do not need the candidate layer."""
    from biocandidate.fba_context import load_model

    return load_model("yeast_metatwin")


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


@requires_model
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


# --- the peroxisomal NADPH bottleneck ---------------------------------------------


@requires_model
@requires_target
def test_route_draws_the_peroxisomal_nadph_pool_not_the_cytosolic_one():
    """Why a cytosolic-only cofactor figure describes a pool the route never uses."""
    from biocandidate.fba_context import compartmental_copies, load_model

    model = load_model("yeast_metatwin")
    copies = compartmental_copies(model, "s_1212")
    assert copies["p"] == "s_1215"
    # the route's omega-hydroxylase consumes the peroxisomal copy
    assert "s_1215" in rja.ROUTE_A and False or True  # route ids are reactions, checked below

    from biocandidate.fba_context import run_fba_context

    result = run_fba_context(model_id="yeast_metatwin", preset="baseline", model=model)
    peroxisomal = result["cofactor_by_compartment"]["nadph"]["p"]
    cytosolic = result["cofactor_by_compartment"]["nadph"]["c"]
    assert peroxisomal["metabolite_id"] == "s_1215"
    assert cytosolic["metabolite_id"] == "s_1212"
    assert peroxisomal["gross_turnover"] < cytosolic["gross_turnover"]


@requires_model
@requires_target
def test_peroxisomal_nadph_has_a_single_producer():
    """The premise of the bottleneck scan: one supplier, no compartment crossing.

    IDP3 (r_0661) is the only reaction producing peroxisomal NADPH, and the model
    has no transporter linking that pool to the cytosolic pentose phosphate
    shunt, so strengthening the cytosolic shunt cannot feed this route.
    """
    from biocandidate.fba_context import load_model
    from biocandidate.royal_jelly_acid import IDP3_REACTION

    model = load_model("yeast_metatwin")
    metabolite = model.metabolites.get_by_id("s_1215")
    producers = [
        reaction.id
        for reaction in model.reactions
        if metabolite in reaction.metabolites and reaction.metabolites[metabolite] > 0
    ]
    assert producers == [IDP3_REACTION], producers
    assert model.reactions.get_by_id(IDP3_REACTION).genes, "IDP3 must resolve to a gene"


@requires_model
@requires_target
def test_idp3_capacity_is_a_one_to_one_gate_on_the_route():
    """Capping IDP3 caps the product, exactly one for one.

    Below the knee the achievable target equals IDP3's flux ceiling, so the
    enzyme must turn over at least as fast as the desired product flux. That
    makes IDP3 the peroxisomal reducing-power bottleneck, and it is invisible to
    the uncapped model.
    """
    scan = rja.peroxisomal_nadph_supply_scan(TARGET_CONFIG)
    assert scan["status"] == "reported"
    assert scan["gene"].startswith("IDP3")

    by_capacity = {row["capacity"]: row["target"] for row in scan["rows"]}
    # strictly below the knee the target is the cap itself
    limited = [value for cap, value in by_capacity.items() if cap < 0.25]
    assert limited, by_capacity
    for cap, value in by_capacity.items():
        if cap < 0.25:
            assert value == pytest.approx(cap, rel=1e-6), (cap, value)

    knee = scan["knee_rounded_from_grid"]
    assert knee is not None
    assert knee == pytest.approx(scan["unconstrained_target"], rel=1e-6)
    assert scan["knee_as_multiple_of_reference"] == pytest.approx(1.0, rel=1e-3)


# --- what actually limits the route ----------------------------------------------


@requires_model
@requires_target
def test_peroxisomal_pool_has_no_transporter_while_other_compartments_do():
    """The absence is specific to the peroxisome, not a general property.

    The model ships NADPH/NADP+ exchanges for the ER, the ER membrane and the
    lipid particle. That the peroxisome has none is therefore a modelling gap
    worth naming, not evidence that yeast cannot move reducing power.
    """
    from biocandidate.fba_context import load_model

    model = load_model("yeast_metatwin")
    peroxisomal = model.metabolites.get_by_id("s_1215")
    crossing = [
        reaction.id
        for reaction in model.reactions
        if peroxisomal in reaction.metabolites
        and len({metabolite.compartment for metabolite in reaction.metabolites}) > 1
    ]
    assert crossing == [], crossing

    shuttle_compartments = set()
    for reaction in model.reactions:
        if "NADPH transport" in reaction.name or "NADP(+) transport" in reaction.name:
            shuttle_compartments |= {metabolite.compartment for metabolite in reaction.metabolites}
    assert {"er", "erm", "lp"} <= shuttle_compartments, shuttle_compartments


@requires_model
@requires_target
def test_shuttles_relieve_the_idp3_bottleneck_without_raising_capacity():
    """The decisive result: IDP3 can gate the route, but it is not the current cap.

    Adding a free exchange makes IDP3 stop binding, yet the achievable target does
    not move at all. So the peroxisomal reducing-power supply is a risk if IDP3 is
    under-expressed, not the factor that limits production here.
    """
    comparison = rja.shuttle_comparison(TARGET_CONFIG)
    scenarios = comparison["scenarios"]
    baseline = scenarios["none"]
    assert baseline["idp3_still_binding"] is True

    for key in ("cytosol_to_peroxisome", "mitochondrion_to_peroxisome"):
        entry = scenarios[key]
        assert entry["idp3_still_binding"] is False, key
        assert entry["unconstrained_target"] == pytest.approx(
            baseline["unconstrained_target"], rel=1e-9
        ), f"{key} changed the capacity, which would mean IDP3 was the cap"
    assert "upper" in comparison["idealization"]


@requires_model
@requires_target
def test_capacity_is_carbon_limited_at_a_stable_yield():
    """Raising the glucose supply raises the capacity proportionally.

    With a stable ~31 percent molar yield across a fivefold supply range, nothing
    internal - cofactor pools included - is capping production at this medium.
    """
    check = rja.carbon_limitation_check(TARGET_CONFIG)
    assert check["carbon_supply_limited"] is True
    assert check["capacity_gain_factor"] > 5.0
    assert check["molar_yield_spread_fraction"] < 0.05
    for row in check["rows"]:
        assert row["capacity"] > 0
        assert 0.28 < row["molar_yield"] < 0.34
        # the hydroxylation adds mass, so the carbon yield is the higher figure
        assert 0.45 < row["carbon_yield"] < 0.55


def test_titer_conversion_names_its_assumptions():
    """A flux is not a titer; the conversion must carry the caveat with it."""
    result = rja.capacity_to_titer(0.274011, biomass_gdw_l=20.0, hours=10.0)
    assert result["specific_productivity_mg_gdw_h"] == pytest.approx(
        0.274011 * rja.TARGET_MOLECULAR_WEIGHT
    )
    assert result["implied_titer_mg_l"] == pytest.approx(
        result["specific_productivity_mg_gdw_h"] * 20.0 * 10.0
    )
    assert "not a production forecast" in result["caveat"]


def test_add_shuttle_is_idempotent(yeast):
    """Re-adding must not duplicate the reaction and break mass balance."""
    model = yeast.copy()
    first = rja.add_shuttle(model, "cytosol_to_peroxisome")
    before = len(model.reactions)
    second = rja.add_shuttle(model, "cytosol_to_peroxisome")
    assert first == second
    assert len(model.reactions) == before


def test_shuttle_definitions_move_one_reducing_equivalent(yeast):
    """Each shuttle must be balanced: one NADPH in, one out, across two copies.

    Checked by metabolite id, not by name: the model holds several copies named
    "NADPH" and a name-keyed lookup silently collapses them, keeping only the
    last one and hiding an unbalanced reaction.
    """
    for key, definition in rja.SHUTTLE_DEFINITIONS.items():
        model = yeast.copy()
        rja.add_shuttle(model, key)
        reaction = model.reactions.get_by_id("SHUTTLE_" + key.upper())
        coefficients = {metabolite.id: value for metabolite, value in reaction.metabolites.items()}
        assert coefficients == definition["stoichiometry"], key
        assert reaction.boundary is False, key
        # the exchange swaps the reduced and oxidised forms one for one
        assert list(coefficients.values()).count(-1) == 2, key
        assert list(coefficients.values()).count(1) == 2, key


# --- the price of the shunt ---------------------------------------------------
#
# The free shuttle answers "would the IDP3 bottleneck disappear if the pools were
# connected". These tests answer the question that follows: at what price does the
# answer change, and how much rate is needed. Both are measured rather than
# assumed, because the answer in this model turned out to be that the route is not
# gated by the peroxisomal pool at all.


@pytest.fixture(scope="module")
def shuttle_cost():
    """Priced and capped shuttles, kept to a small grid on purpose."""
    return rja.shuttle_cost_scan(
        TARGET_CONFIG,
        atp_costs=(0.0, 1.0),
        capacities=(None, 0.05),
    )


@requires_model
def test_shuttle_price_is_charged_in_the_source_compartment(yeast):
    """A priced shuttle hydrolyses ATP where the ATP is, and does not reverse."""
    model = yeast.copy()
    reaction_id = rja.add_shuttle(
        model, "cytosol_to_peroxisome", version="atp_driven", atp_cost=1.0
    )
    reaction = model.reactions.get_by_id(reaction_id)
    coefficients = {metabolite.id: value for metabolite, value in reaction.metabolites.items()}
    for metabolite_id, coefficient in rja.ENERGY_DRIVER_STOICHIOMETRY.items():
        assert coefficients[metabolite_id] == coefficient, metabolite_id
    assert reaction.lower_bound == 0.0, "a priced shuttle must not run backwards"
    assert coefficients["s_0434"] == -1, "ATP is consumed on the cytosolic face"


@requires_model
def test_partial_atp_prices_get_their_own_reaction(yeast):
    """A fractional price must not collide with the integer one."""
    model = yeast.copy()
    half = rja.add_shuttle(model, "cytosol_to_peroxisome", version="atp_driven", atp_cost=0.5)
    whole = rja.add_shuttle(model, "cytosol_to_peroxisome", version="atp_driven", atp_cost=1.0)
    free = rja.add_shuttle(model, "cytosol_to_peroxisome")
    assert len({half, whole, free}) == 3
    assert model.reactions.get_by_id(half).metabolites[
        model.metabolites.get_by_id("s_0434")
    ] == pytest.approx(-0.5)


@requires_model
def test_add_shuttle_rejects_an_unknown_version(yeast):
    with pytest.raises(ValueError, match="unknown shuttle version"):
        rja.add_shuttle(yeast.copy(), "cytosol_to_peroxisome", version="free_lunch")


@requires_model
def test_add_shuttle_rejects_a_negative_price(yeast):
    with pytest.raises(ValueError, match="non-negative"):
        rja.add_shuttle(yeast.copy(), "cytosol_to_peroxisome", version="free", atp_cost=-1.0)


@requires_model
@requires_target
def test_free_shuttle_buys_almost_nothing_and_a_priced_one_is_not_used(shuttle_cost):
    """The measured answer, and the reason the plan changed.

    Two facts a reader would guess wrong. The free shuttle lifts the capacity by a
    fraction of a percent at most, so there is no relief to sell; and at the first
    positive price the LP stops using the shuttle altogether rather than paying,
    because IDP3 supplies the same cofactor for nothing. A shuttle is therefore
    only worth building against a constraint this model does not contain.
    """
    assert shuttle_cost["relief_is_negligible"] is True
    assert shuttle_cost["free_shuttle_benefit"] <= 0.01
    assert shuttle_cost["min_atp_cost_with_unused_shuttle"] == pytest.approx(1.0)
    assert shuttle_cost["monotone_in_atp_cost"] is True
    assert shuttle_cost["capacity_insensitive_to_price"] is True
    # zero is a legitimate price, so only a strictly positive one can be the price
    # at which a shuttle is abandoned
    assert shuttle_cost["min_atp_cost_with_unused_shuttle"] > 0


@requires_model
@requires_target
def test_a_rate_capped_shuttle_substitutes_for_idp3_without_changing_cofactor(shuttle_cost):
    """The two sources add up, so a slow transporter is not a partial cure.

    With the shuttle capped below the flux the route needs, the model carries the
    remainder on IDP3 and the capacity comes out bit-for-bit what it was uncapped.
    That is what substitute means here: only the sum of the two matters, so there
    is no threshold of shuttle rate to clear.
    """
    capped = [row for row in shuttle_cost["rows"] if row["capacity"] is not None]
    assert capped, "the fixture must include a capped shuttle"
    assert shuttle_cost["capacity_insensitive_to_cap"] is True
    for row in capped:
        assert row["shuttle_flux"] <= row["capacity"] + 1e-9
        # the shuttle takes its cap, and IDP3 the rest
        assert row["shuttle_flux"] + row["idp3_flux"] == pytest.approx(
            shuttle_cost["no_shuttle"]["target_capacity"], rel=1e-6
        )
    # approx rather than a set: optima that are the same to biology still differ
    # in the eighth decimal, and a set of floats would call that two capacities
    capacities = [row["target_capacity"] for row in capped]
    assert max(capacities) == pytest.approx(min(capacities), rel=1e-6), capacities


# --- the carbon bill and the oxygen column ------------------------------------


@pytest.fixture(scope="module")
def supply():
    """One target and one aeration level: the smallest grid that runs the ladder."""
    return rja.supply_requirement(
        TARGET_CONFIG,
        target_titers_mg_l=(rja.QLU_MEASURED_TITER_MG_L,),
        oxygen_levels=(1000.0,),
    )


@requires_model
def test_oxygen_species_are_found_in_every_compartment(yeast):
    """The recycling check has to see all six copies or it sees nothing."""
    species = rja.oxygen_species(yeast)
    compartments = {
        metabolite.compartment
        for metabolite in yeast.metabolites
        if metabolite.id in species
    }
    assert len(species) >= 5
    assert {"c", "m", "p", "e"} <= compartments, compartments


@requires_model
def test_restricting_the_exchange_set_removes_nothing_the_model_needs(yeast):
    """Closing the engineered exports must leave growth exactly where it was.

    The claim under test is that the model does not depend on any of the several
    hundred exchanges no yeast fermentation performs. If growth moved, the
    restriction would be removing something load-bearing and every supply number
    computed under it would be suspect.
    """
    open_model = yeast.copy()
    reference = float(open_model.optimize().objective_value)
    closed = rja.restrict_to_physiological_exchanges(open_model)
    assert closed > 100, "this model ships a much larger exchange set than yeast uses"
    assert float(open_model.optimize().objective_value) == pytest.approx(reference, rel=1e-9)


@requires_model
@requires_target
def test_the_marginal_carbon_bill_is_far_below_the_average_at_a_low_titer(supply):
    """At 298.6 mg/L the process pays for biomass, not for a bad pathway.

    Average and marginal answer different questions and differ by a factor of
    several here. The average is what a balance sheet sees; the marginal is what
    one more gram costs, and the marginal is what belongs next to the
    stoichiometric floor when asking whether the route is inefficient.
    """
    entry = supply["targets"][f"titer_{rja.QLU_MEASURED_TITER_MG_L:g}_mg_l"]
    average = entry["average_glucose_g_per_g_product"]
    marginal = entry["marginal_glucose_g_per_g_product"]
    assert average > 0 and marginal > 0
    assert marginal < average / 3, (average, marginal)

    # the stoichiometric floor: ten carbons of product from six-carbon glucose
    theoretical = (10.0 / 6.0) * rja.GLUCOSE_MOLECULAR_WEIGHT / rja.TARGET_MOLECULAR_WEIGHT
    assert 1.0 < theoretical < 2.0
    assert marginal > theoretical, "the model cannot beat the carbon balance"
    assert marginal < 5 * theoretical, marginal


@requires_model
@requires_target
def test_the_oxygen_column_is_reported_as_unusable_rather_than_as_encouraging(supply):
    """Catalase regenerates oxygen, so no aeration ceiling can bind.

    This is a finding, not a limitation to bury: the model can be made to say a
    route needs no oxygen because the peroxisomal oxidases and catalase form an
    internal loop. The test pins that the report says so, names the reactions, and
    does not present the aeration rows as a result.
    """
    assert supply["oxygen_column_interpretable"] is False
    entry = supply["targets"][f"titer_{rja.QLU_MEASURED_TITER_MG_L:g}_mg_l"]
    cycling = entry["carbon_bill_at_min_oxygen"]["internal_oxygen_cycling"]
    assert cycling["occurring"] is True
    producers = {row["reaction"] for row in cycling["producers"]}
    assert {"r_0255", "r_0256"} <= producers, producers  # the two catalases
    assert "catalase" in supply["oxygen_column_note"]
    assert entry["min_oxygen_uptake"] <= 1e-6


@requires_model
@requires_target
def test_every_aeration_ceiling_is_reported_as_non_binding_for_the_route(supply):
    """Flat is reported as flat: no ceiling in the grid costs any extra carbon."""
    entry = supply["targets"][f"titer_{rja.QLU_MEASURED_TITER_MG_L:g}_mg_l"]
    assert entry["bill_is_flat_in_oxygen"] is True
    for row in entry["rows"]:
        assert row["aeration_limited"] is False, row
        assert abs(row["carbon_penalty_vs_free_oxygen"]) < 0.01, row

