"""Royal jelly acid (10-HDA) route analysis on the yeast genome-scale model.

Focused application of the general FBA-context capability. Royal jelly acid,
(E)-10-hydroxy-2-decenoic acid, is the target of the user's 10H2DA programme;
this module answers the questions that programme actually faces:

  * how much target flux can each documented route sustain, and what limits it
  * how much NADPH the omega-hydroxylation step needs, and how much FAD /
    acetyl-CoA the competing oxidation and precursor-supply sides carry
  * which route is cheaper in enzyme investment at equal flux

Route chemistry comes from the external project's own target configuration
(read-only) and from the candidate-extension layer the kcat-constrained FBA
pilot reconstructed. Nothing here is experimental validation: substrate supply,
enzyme levels and fluxes are model quantities.

Two modelling facts a reader must not miss, both verified against the model
rather than assumed:

  1. POX1 (r_0120) is written in this model WITHOUT an explicit flavin, so its
     activity is read from the r_0120 flux, NOT from FAD turnover.
  2. At the default yeast medium, flavin turnover is essentially zero (a residual
     9.2e-06 through r_0440, i.e. 1e-5 of the growth rate). Flavins therefore
     cannot serve as a degradation-side handle in this model; the acyl-CoA
     thioesters and the r_0120 flux can.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .fba_context import (
    CORE_COFACTORS,
    MODEL_REGISTRY,
    cofactor_turnover,
    cofactor_turnover_by_compartment,
    enzyme_capacity_scan,
    load_model,
    marginal_cofactor_cost,
    oxidative_burden,
)
from .fba_context import (
    YEAST_METATWIN_MODEL_ID as MODEL_ID,
)

BIOMASS_ID = "r_2111"
TARGET_DEMAND = "DM_CAND_10H2DA_P"
POX1_REACTION = "r_0120"
BIOMASS_FLOOR_FRACTION = 0.1
# The external project's Sep 2 report states this maximum on the model-file
# default medium; reproducing it is this module's self-check.
DOCUMENTED_ANCHOR = 0.274011
ANCHOR_TOLERANCE = 1e-5
SECONDS_PER_HOUR = 3600.0

# The route-specific reactions are created by the candidate-extension layer; the
# labels below match the external project's target configuration.
ROUTE_A = {
    "CAND_T2DEC_THIOESTERASE_P": "TES1/PTE1 (free-acid route)",
    "CAND_T2DEC_OMEGA_HYDROXYLASE_P": "CYP52M1 (free-acid route)",
}
ROUTE_B = {
    "CAND_T2DEC_COA_OMEGA_HYDROXYLASE_P": "CYP52M1 (CoA-thioester route)",
    "CAND_10H2DA_COA_THIOESTERASE_P": "TES1/PTE1 (CoA-thioester route)",
}
ROUTE_INTERMEDIATE_DEMAND = "DM_CAND_T2DEC_P"

# Cofactors reported for the feasibility view: synthesis side (NADPH for the P450),
# energy, and the fatty-acid precursor thioesters. Flavins are listed for
# completeness but are documented above as uninformative in this model.
FEASIBILITY_COFACTORS = ("nadph", "nadh", "atp", "accoa", "malcoa")

CLAIM_BOUNDARY = (
    "Route analysis on a genome-scale model. Target flux, enzyme investment and "
    "cofactor turnover are model quantities; enzyme levels are assumed and kcat "
    "inputs are uncalibrated predictions. Nothing here is experimental "
    "validation or a production forecast."
)


def load_target_config(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_candidate_extension(model, target_config: dict[str, Any]) -> list[str]:
    """Reconstruct the external project's 10H2DA candidate layer, idempotently.

    Moved here from scripts/fusion_yeast_kcat_fba_pilot.py so the pilot and the
    route analysis share one implementation instead of two copies.
    """
    import cobra

    added: list[str] = []
    for row in target_config["compounds"]:
        metabolite_id = row["model_metabolite_id"]
        if metabolite_id in model.metabolites:
            continue
        model.add_metabolites([
            cobra.Metabolite(
                metabolite_id,
                formula=row["formula"],
                charge=int(row["charge"]),
                name=row["name"],
                compartment="p",
            )
        ])
        added.append(metabolite_id)

    for row in target_config["candidate_reactions"]:
        if row["reaction_id"] in model.reactions:
            continue
        reaction = cobra.Reaction(row["reaction_id"])
        reaction.name = row["name"]
        reaction.lower_bound = 0.0
        reaction.upper_bound = 1000.0
        model.add_reactions([reaction])
        stoichiometry = row["stoichiometry"]
        if isinstance(stoichiometry, str):
            stoichiometry = json.loads(stoichiometry)
        model.reactions.get_by_id(row["reaction_id"]).add_metabolites(
            {model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()}
        )
        added.append(row["reaction_id"])

    for demand_id, stoichiometry in (
        ("DM_s_1507", {"s_1507": -1}),
        (ROUTE_INTERMEDIATE_DEMAND, {"cand_t2dec_p": -1}),
        (TARGET_DEMAND, {"cand_10h2da_p": -1}),
    ):
        if demand_id in model.reactions:
            continue
        reaction = cobra.Reaction(demand_id)
        reaction.name = f"demand {demand_id}"
        reaction.lower_bound = 0.0
        reaction.upper_bound = 1000.0
        model.add_reactions([reaction])
        reaction.add_metabolites(
            {model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()}
        )
        added.append(demand_id)
    return added


def _route_reactions(model) -> dict[str, list[str]]:
    """Only report routes whose candidate reactions actually exist."""
    routes: dict[str, list[str]] = {}
    if all(r in model.reactions for r in ROUTE_A):
        routes["route_A_free_acid"] = list(ROUTE_A)
    if all(r in model.reactions for r in ROUTE_B):
        routes["route_B_coa_thioester"] = list(ROUTE_B)
    return routes


def route_context(
    route: str,
    *,
    pox1_closed: bool,
    glucose: float = 1.0,
    oxygen: float = 1000.0,
):
    """Build the scenario applier the marginal analysis needs.

    Kept separate from analyze() so the marginal curve is measured under the same
    medium, route and gate as the capacity it is derived from.
    """
    routes = {"A": ROUTE_A, "B": ROUTE_B}
    if route not in routes:
        raise ValueError(f"unknown route {route!r}; known: {sorted(routes)}")
    other = "B" if route == "A" else "A"

    def apply(model) -> None:
        model.reactions.get_by_id("r_1714").lower_bound = -float(glucose)
        model.reactions.get_by_id("r_1992").lower_bound = -float(oxygen)
        if pox1_closed:
            model.reactions.get_by_id(POX1_REACTION).bounds = (0.0, 0.0)
        for reaction_id in routes[other]:
            model.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)

    return apply


def marginal_route_cost(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    glucose: float = 1.0,
    oxygen: float = 1000.0,
    level_fractions: tuple[float, ...] = (0.25, 0.5, 0.75, 1.0),
) -> dict[str, Any]:
    """What one more unit of 10-HDA costs the cell, per cofactor and in growth."""
    from .fba_context import MODEL_REGISTRY as registry

    model = load_model(MODEL_ID)
    build_candidate_extension(model, load_target_config(target_config_path))
    spec = registry[MODEL_ID]

    report = marginal_cofactor_cost(
        model,
        TARGET_DEMAND,
        cofactor_metabolites={
            label: spec["cofactor_metabolites"][label]
            for label in ("nadph", "nadh", "atp", "accoa", "malcoa")
        },
        oxidant_metabolites=spec["oxidant_metabolites"],
        objective_reaction=BIOMASS_ID,
        level_fractions=level_fractions,
        context=route_context(
            route, pox1_closed=pox1_closed, glucose=glucose, oxygen=oxygen
        ),
    )
    report["route"] = route
    report["pox1_closed"] = pox1_closed
    report["medium"] = {"glucose": glucose, "oxygen": oxygen}
    report["claim_boundary"] = CLAIM_BOUNDARY
    return report


IDP3_REACTION = "r_0661"  # peroxisomal NADP-isocitrate dehydrogenase


def peroxisomal_nadph_supply_scan(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    capacities: tuple[float, ...] = (
        0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0,
    ),
) -> dict[str, Any]:
    """Is the peroxisomal NADPH supply the bottleneck for this route?

    In this model the peroxisomal NADPH pool has exactly one producer, the
    peroxisomal isocitrate dehydrogenase IDP3 (r_0661), there is no
    compartment-crossing transporter connecting it to the cytosolic pentose
    phosphate shunt, and the model leaves IDP3 uncapped. Capping it and
    re-maximising the target shows whether the route depends on a flux IDP3 must
    actually be able to deliver.
    """
    model = load_model(MODEL_ID)
    build_candidate_extension(model, load_target_config(target_config_path))
    report = enzyme_capacity_scan(
        model,
        IDP3_REACTION,
        TARGET_DEMAND,
        capacities,
        context=route_context(route, pox1_closed=pox1_closed),
    )
    report["gene"] = "IDP3 (YNL009W)"
    report["role"] = "sole producer of the peroxisomal NADPH pool the route consumes"
    report["route"] = route
    report["pox1_closed"] = pox1_closed
    report["claim_boundary"] = CLAIM_BOUNDARY
    return report


# Idealized reducing-power shuttles: NADPH and its oxidised partner swapped
# across a membrane, one for one. Real shuttles are not free - they run through
# malate/oxaloacetate or a membrane protein and consume something - so these are
# upper bounds on what a shuttle could buy, not predictions. They exist to answer
# one question: if the peroxisomal pool could be fed from elsewhere, would the
# IDP3 bottleneck disappear?
SHUTTLE_DEFINITIONS = {
    "cytosol_to_peroxisome": {
        "stoichiometry": {"s_1212": -1, "s_1211": -1, "s_1207": 1, "s_1215": 1},
        "name": "cytosol-peroxisome NADPH/NADP+ exchange (idealized)",
        "source_compartment": "c",
    },
    "mitochondrion_to_peroxisome": {
        "stoichiometry": {"s_1214": -1, "s_1211": -1, "s_1210": 1, "s_1215": 1},
        "name": "mitochondrion-peroxisome NADPH/NADP+ exchange (idealized)",
        "source_compartment": "m",
    },
}


# A real shuttle is not free, and the price has a location. An ATP-hydrolysing
# driver - the ABC-transporter pattern - hydrolyses ATP on the cytosolic face, so
# ATP and water are consumed there and ADP and phosphate released there; the
# transport is the service the ATP pays for. All four ids are cytosolic because
# this model has no peroxisomal phosphate, which is itself worth knowing: a
# phosphate-releasing step cannot be written inside the peroxisome.
ENERGY_DRIVER_STOICHIOMETRY = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}

SHUTTLE_VERSIONS = {
    "free": {
        "atp_cost": 0.0,
        "description": (
            "no enzyme, no energy, no rate limit: an upper bound on what a shuttle "
            "could buy, not a mechanism"
        ),
    },
    "atp_driven": {
        "atp_cost": 1.0,
        "description": (
            "one ATP hydrolysed on the cytosolic face per NADPH moved into the "
            "peroxisome; irreversible, because paying ATP to export reducing power "
            "has no mechanism"
        ),
    },
}


def add_shuttle(
    model,
    shuttle_key: str,
    *,
    version: str = "free",
    atp_cost: float | None = None,
    capacity: float | None = None,
) -> str:
    """Add one shuttle reaction, idempotently. Returns its id.

    ``version="free"`` with no price and no cap reproduces the original idealized
    exchange. Any other version charges the transport: ``atp_cost`` ATP (default
    taken from SHUTTLE_VERSIONS) is hydrolysed per NADPH moved, and ``capacity``
    caps the rate, which is how a finite amount of transporter is represented.
    The reaction id encodes the version and the price, so two variants cannot
    silently overwrite each other.
    """
    import cobra

    definition = SHUTTLE_DEFINITIONS[shuttle_key]
    if version not in SHUTTLE_VERSIONS:
        raise ValueError(
            f"unknown shuttle version {version!r}; known: {sorted(SHUTTLE_VERSIONS)}"
        )
    price = SHUTTLE_VERSIONS[version]["atp_cost"] if atp_cost is None else float(atp_cost)
    if price < 0:
        raise ValueError(f"atp_cost must be non-negative, got {price}")

    idealized = version == "free" and price == 0.0
    reaction_id = "SHUTTLE_" + shuttle_key.upper()
    if not idealized:
        # fractional prices are legitimate (a shuttle that charges one ATP every
        # other turnover), so the id has to survive a decimal point
        reaction_id += "_" + version.upper() + f"_{price:g}".replace(".", "P")
    if reaction_id in model.reactions:
        return reaction_id

    stoichiometry = dict(definition["stoichiometry"])
    if price:
        for metabolite_id, coefficient in ENERGY_DRIVER_STOICHIOMETRY.items():
            stoichiometry[metabolite_id] = (
                stoichiometry.get(metabolite_id, 0.0) + coefficient * price
            )

    reaction = cobra.Reaction(reaction_id)
    reaction.name = (
        definition["name"]
        if idealized
        else f"{definition['name']}, {price:g} ATP per NADPH"
    )
    # The idealized exchange stays reversible: it is meant to be an upper bound,
    # and making it one-directional would quietly weaken that bound. A priced
    # shuttle is one-directional, since the price is paid per NADPH moved in.
    reaction.lower_bound = -1000.0 if idealized else 0.0
    reaction.upper_bound = 1000.0 if capacity is None else float(capacity)
    model.add_reactions([reaction])
    reaction.add_metabolites(
        {model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()}
    )
    return reaction_id


def shuttle_comparison(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    idp3_capacities: tuple[float, ...] = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0),
) -> dict[str, Any]:
    """Does opening a reducing-power connection relieve the IDP3 bottleneck?

    The model connects the cytosolic NADPH pool to the ER, the ER membrane and
    the lipid particle, but not to the peroxisome: the peroxisomal NADPH/NADP+
    pair appears in three reactions in total and has exactly one producer. If the
    route's yield is limited by that single enzyme, connecting the pool should
    lift the limit - and that is the modelling form of the compartmentalization
    question the 2025 yeast paper addressed by moving its P450 into mitochondria.
    """
    context = route_context(route, pox1_closed=pox1_closed)
    report: dict[str, Any] = {
        "route": route,
        "pox1_closed": pox1_closed,
        "scenarios": {},
        "shuttle_definitions": {
            key: value["name"] for key, value in SHUTTLE_DEFINITIONS.items()
        },
        "claim_boundary": CLAIM_BOUNDARY,
        "idealization": (
            "Shuttles are free one-for-one exchanges with no enzyme, no energy cost "
            "and no competing demand. They therefore upper-bound the benefit of "
            "connecting the pools, and a real shuttle must be shown to pay for itself."
        ),
    }

    for key in ("none", *SHUTTLE_DEFINITIONS):
        model = load_model(MODEL_ID)
        build_candidate_extension(model, load_target_config(target_config_path))
        if key != "none":
            add_shuttle(model, key)

        def apply(scoped, _key=key):
            context(scoped)

        scan = enzyme_capacity_scan(
            model,
            IDP3_REACTION,
            TARGET_DEMAND,
            idp3_capacities,
            context=apply,
        )
        report["scenarios"][key] = {
            "shuttle": None if key == "none" else SHUTTLE_DEFINITIONS[key]["name"],
            "unconstrained_target": scan["unconstrained_target"],
            "idp3_knee": scan["knee_rounded_from_grid"] or scan["knee_capacity"],
            "idp3_still_binding": scan["knee_rounded_from_grid"] is not None,
            "capacity_vs_no_shuttle": (
                scan["unconstrained_target"] / report["scenarios"]["none"]["unconstrained_target"]
                if "none" in report["scenarios"]
                else None
            ),
        }
    return report


def shuttle_cost_scan(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    shuttle_key: str = "cytosol_to_peroxisome",
    atp_costs: tuple[float, ...] = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0),
    capacities: tuple[float | None, ...] = (None,),
    biomass_floor_fraction: float = BIOMASS_FLOOR_FRACTION,
    glucose: float = 1.0,
    oxygen: float = 1000.0,
) -> dict[str, Any]:
    """What does the IDP3 relief cost, and where does it stop being worth it?

    ``shuttle_comparison`` answers the yes/no question with a free shuttle. That
    free shuttle is not a mechanism, so this scan prices it: each NADPH moved
    into the peroxisome costs ``atp_cost`` ATP, and the rate can be capped to
    represent a finite amount of transporter. Both are swept, and the achievability
    of the target is reported against the same biomass floor throughout.

    The reading that matters is the ATP price at which the relief disappears: a
    shuttle is only worth building if it comes in under that price. The cap sweep
    answers the separate question of whether a slow shuttle is enough, which is
    what decides between a native transporter and an engineered one.
    """
    from .fba_context import MODEL_REGISTRY as registry

    target_config = load_target_config(target_config_path)
    spec = registry[MODEL_ID]
    base = load_model(MODEL_ID)
    build_candidate_extension(base, target_config)
    context = route_context(
        route, pox1_closed=pox1_closed, glucose=glucose, oxygen=oxygen
    )

    def solve(*, price: float | None, cap: float | None, key: str | None) -> dict[str, Any]:
        model = base.copy()
        reaction_id = None
        if key is not None:
            reaction_id = add_shuttle(
                model, key, version="atp_driven", atp_cost=price, capacity=cap
            )
        with model:
            context(model)
            baseline = model.optimize()
            if baseline.status != "optimal":
                return {"status": f"baseline {baseline.status}"}
            model.reactions.get_by_id(BIOMASS_ID).lower_bound = (
                biomass_floor_fraction * float(baseline.objective_value)
            )
            model.objective = TARGET_DEMAND
            solution = model.optimize()
            if solution.status != "optimal":
                return {"status": solution.status, "target_capacity": None}
            capacity = float(solution.objective_value)
            shuttle_flux = (
                float(solution.fluxes[reaction_id]) if reaction_id is not None else 0.0
            )
            glucose_uptake = -float(solution.fluxes["r_1714"])
            turnover = cofactor_turnover(model, solution, spec["cofactor_metabolites"])
            gross_atp = turnover["atp"].get("gross_turnover") or 0.0
            atp_drain = shuttle_flux * (price or 0.0)
            return {
                "status": "optimal",
                "atp_cost": price,
                "capacity": cap,
                "target_capacity": capacity,
                "shuttle_flux": shuttle_flux,
                "atp_drain": atp_drain,
                "atp_drain_fraction_of_gross_atp": (
                    atp_drain / gross_atp if gross_atp else None
                ),
                "idp3_flux": float(solution.fluxes[IDP3_REACTION]),
                "growth": float(solution.fluxes[BIOMASS_ID]),
                "glucose_uptake": glucose_uptake,
                "molar_yield": capacity / glucose_uptake if glucose_uptake else None,
                "carbon_yield": (
                    capacity * 10.0 / (glucose_uptake * 6.0)
                    if glucose_uptake
                    else None
                ),
            }

    reference = solve(price=None, cap=None, key=None)
    if reference.get("status") != "optimal":
        return {
            "status": reference.get("status"),
            "error": "the route is not achievable before any shuttle is added",
            "claim_boundary": CLAIM_BOUNDARY,
        }

    rows: list[dict[str, Any]] = []
    for cap in capacities:
        for price in atp_costs:
            row = solve(price=price, cap=cap, key=shuttle_key)
            row["benefit_over_no_shuttle"] = (
                (row["target_capacity"] - reference["target_capacity"])
                / reference["target_capacity"]
                if row.get("target_capacity") is not None
                else None
            )
            rows.append(row)

    solved = [row for row in rows if row.get("target_capacity") is not None]
    free = sorted(
        (row for row in solved if row["capacity"] is None),
        key=lambda row: row["atp_cost"],
    )
    # relief has to be monotone in the price: a pricier shuttle cannot buy more
    monotone_in_price = all(
        later["target_capacity"] <= earlier["target_capacity"] + 1e-9
        for earlier, later in zip(free, free[1:])
    )
    # At a zero price the shuttle and IDP3 are interchangeable routes to the same
    # cofactor, so the LP picks one arbitrarily and the split between them carries
    # no information. What carries information is the smallest price at which the
    # shuttle stops being used at all: that is the margin a real transporter has
    # to beat, and it is measured against IDP3, not against zero.
    free_price_benefit = free[0]["benefit_over_no_shuttle"] if free else None
    unused_prices = [row["atp_cost"] for row in free if abs(row["shuttle_flux"]) <= 1e-9]
    priced_out = [
        row["atp_cost"]
        for row in free
        if (row["benefit_over_no_shuttle"] or 0.0) <= 0.01
    ]
    relief_prices = [
        row["atp_cost"]
        for row in free
        if (row["benefit_over_no_shuttle"] or 0.0) > 0.01
    ]

    cap_rows = sorted(
        (row for row in solved if row["capacity"] is not None and row["atp_cost"] == 0.0),
        key=lambda row: row["capacity"],
    )
    uncapped = free[0]["target_capacity"] if free else None
    cap_threshold = None
    for row in cap_rows:
        if uncapped and row["target_capacity"] >= 0.99 * uncapped:
            cap_threshold = row["capacity"]
            break
    # When the capacity does not move at all across a grid, the extreme value of
    # that grid is not a threshold and must not be read as one: the honest report
    # is that the variable does not matter on the tested range. The band is
    # relative and deliberately wider than solver noise (observed vertex jitter
    # between optima that are the same to biology is a few parts in 1e8).
    def flat(values: list[float]) -> bool:
        return bool(values) and max(values) - min(values) <= 1e-6 * max(
            abs(max(values)), 1.0
        )

    capacity_insensitive_to_cap = flat([row["target_capacity"] for row in cap_rows])
    capacity_insensitive_to_price = flat([row["target_capacity"] for row in free])

    return {
        "status": "reported",
        "route": route,
        "pox1_closed": pox1_closed,
        "shuttle": SHUTTLE_DEFINITIONS[shuttle_key]["name"],
        "shuttle_key": shuttle_key,
        "medium": {"glucose": glucose, "oxygen": oxygen},
        "biomass_floor_fraction": biomass_floor_fraction,
        "no_shuttle": reference,
        "rows": rows,
        "monotone_in_atp_cost": monotone_in_price,
        "free_shuttle_benefit": free_price_benefit,
        "relief_is_negligible": (
            free_price_benefit is not None and free_price_benefit <= 0.01
        ),
        "min_atp_cost_with_unused_shuttle": min(unused_prices) if unused_prices else None,
        "atp_cost_priced_out": priced_out,
        "max_atp_cost_with_relief": max(relief_prices) if relief_prices else None,
        "min_capacity_retaining_full_relief": cap_threshold,
        "capacity_insensitive_to_cap": capacity_insensitive_to_cap,
        "capacity_insensitive_to_price": capacity_insensitive_to_price,
        "interpretation": (
            "Two different statements live in this table. min_atp_cost_with_unused_shuttle "
            "says the cell will not buy the shuttle at that price while IDP3 is "
            "available, so it is the budget a transporter must beat - but only if "
            "there is anything to buy. relief_is_negligible screens for that: when a "
            "free shuttle also fails to raise the capacity, the route is not gated by "
            "the peroxisomal NADPH pool and no transporter, however cheap, changes the "
            "ceiling; the limit is elsewhere and this scan sends you to look for it "
            "rather than to build the shuttle. A capped shuttle at zero price answers "
            "the separate question of how much rate is enough, since a finite "
            "transporter only matters if it is on a binding path. When "
            "capacity_insensitive_to_cap is true the cap range is not a threshold "
            "search: it means the two NADPH sources are substitutes and only their sum "
            "matters, so the smallest cap that holds the capacity is simply the "
            "smallest one tested."
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }


# 10-HDA, C10H18O3, g/mol. Used only to express a model flux in mass terms.
TARGET_MOLECULAR_WEIGHT = 186.25
GLUCOSE_MOLECULAR_WEIGHT = 180.16
TITER_CAVEAT = (
    "A titer is biomass times time times realised specific productivity; the model "
    "has no biomass trajectory and no enzyme, toxicity or secretion limit, so this "
    "number is an upper bound on what the flux could support, not a production "
    "forecast."
)


def capacity_to_titer(
    capacity_mmol_gdw_h: float,
    *,
    biomass_gdw_l: float,
    hours: float,
    molecular_weight: float = TARGET_MOLECULAR_WEIGHT,
) -> dict[str, Any]:
    """Express a model capacity as an implied titer, with assumptions named.

    Reported so the gap against a measured titer can be stated rather than left
    implicit: genome-scale models have no enzyme, toxicity or transport limit, so
    their capacity sits well above what a real strain reaches, and quoting the
    flux without that gap invites reading it as a forecast.
    """
    specific = capacity_mmol_gdw_h * molecular_weight
    return {
        "capacity_mmol_gdw_h": capacity_mmol_gdw_h,
        "molecular_weight_g_mol": molecular_weight,
        "specific_productivity_mg_gdw_h": specific,
        "biomass_gdw_l": biomass_gdw_l,
        "hours": hours,
        "implied_titer_mg_l": specific * biomass_gdw_l * hours,
        "caveat": TITER_CAVEAT,
    }


def carbon_limitation_check(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    glucose_levels: tuple[float, ...] = (1.0, 2.0, 3.0, 5.0),
    oxygen: float = 1000.0,
) -> dict[str, Any]:
    """Is the route limited by carbon supply or by something internal?

    At the optimum the only binding constraint this model reports is the glucose
    exchange, and capacity scales with it at a stable yield, which says the
    internal machinery - cofactors included - is not what caps production. If
    capacity instead plateaued, an internal step would be limiting.
    """
    rows: list[dict[str, Any]] = []
    for glucose in glucose_levels:
        model = load_model(MODEL_ID)
        build_candidate_extension(model, load_target_config(target_config_path))
        with model:
            # the medium must be passed INTO the applier: route_context carries its
            # own defaults and setting the bounds before calling it has no effect
            route_context(
                route, pox1_closed=pox1_closed, glucose=glucose, oxygen=oxygen
            )(model)
            model.objective = TARGET_DEMAND
            solution = model.optimize()
            if solution.status != "optimal":
                rows.append({"glucose": glucose, "status": solution.status, "capacity": None})
                continue
            capacity = float(solution.objective_value)
            rows.append(
                {
                    "glucose": glucose,
                    "status": "optimal",
                    "capacity": capacity,
                    "molar_yield": capacity / glucose,
                    "carbon_yield": capacity * 10.0 / (glucose * 6.0),
                }
            )

    solved = [row for row in rows if row["capacity"] is not None]
    yields = [row["molar_yield"] for row in solved]
    capacities = [row["capacity"] for row in solved]
    gain = (max(capacities) / min(capacities)) if capacities and min(capacities) else None
    spread = (max(yields) - min(yields)) / max(yields) if yields and max(yields) else None
    # Carbon-limited means BOTH: capacity rises with the supply AND the yield holds.
    # Capacity rising would be enough on its own; a stable yield alone proves
    # nothing, because a constant capacity also gives a stable yield.
    carbon_limited = (
        gain is not None and gain > 1.05 and spread is not None and spread < 0.05
    )
    return {
        "status": "reported",
        "route": route,
        "rows": rows,
        "capacity_gain_factor": gain,
        "molar_yield_spread_fraction": spread,
        "carbon_supply_limited": carbon_limited,
        "interpretation": (
            "Capacity rising with the supply while the yield holds means carbon is the "
            "cap and nothing internal is. Capacity rising but less than proportionally "
            "- a falling yield - means an internal step starts to bind at the level "
            "where it happens, and that level is worth investigating."
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }


# Strain-level benchmark from the 2025 yeast paper (fed-batch, ethyl decanoate as
# the substrate fed to dodge decanoate toxicity). Included in the default grid so
# the model's carbon and oxygen bill for a titer someone has actually reached is
# reported rather than left to the reader to interpolate.
QLU_MEASURED_TITER_MG_L = 298.6

# Secreted products whose appearance marks the point where respiration stops
# absorbing all the carbon it is given. Read as a signature, not as a target:
# each one is carbon that did not become product.
OVERFLOW_EXCHANGES = {
    "ethanol": "r_1761",
    "acetate": "r_1634",
    "glycerol": "r_1808",
    "pyruvate": "r_2033",
    "carbon_dioxide": "r_1672",
}

# The carbon ceiling for the minimum-glucose problem. High enough not to bind at
# any plausible demand, so a failure is a real infeasibility rather than the
# ceiling being reached; the answer is the minimum, and the comparison against a
# real feed rate is the reader's.
GLUCOSE_CEILING = 1000.0

# The aeration ceiling, used when oxygen is not the variable under test.
OXYGEN_CEILING = 1000.0

# The exchanges a S. cerevisiae fermentation actually has: the carbon and energy
# gases, the mineral medium, and the products yeast is known to secrete. The
# model ships several hundred exchanges spanning engineered chemistry, and left
# open they change the answer to a redox question - the LP balances electrons by
# exporting methanol or ethyl pyruvate, which no yeast strain makes at a
# measurable rate, and reports a minimum oxygen requirement of zero as an
# artefact of that exchange list rather than as a property of the cell. Closing
# them costs nothing: growth on the reference medium is unchanged to 1e-12.
PHYSIOLOGICAL_EXCHANGES = frozenset({
    # carbon, energy and redox
    "r_1714",  # D-glucose
    "r_1992",  # oxygen
    "r_1672",  # carbon dioxide
    "r_2100",  # water
    "r_1832",  # H+
    # mineral medium
    "r_1654",  # ammonium
    "r_2005",  # phosphate
    "r_2060",  # sulphate
    "r_1861",  # iron(2+)
    "r_4593",  # chloride
    "r_4594",  # Cu2+
    "r_4595",  # Mn(2+)
    "r_4596",  # Zn(2+)
    "r_4597",  # Mg(2+)
    "r_4600",  # Ca2+
    "r_2020",  # potassium
    "r_2049",  # sodium
    # documented fermentative products
    "r_1761",  # ethanol
    "r_1634",  # acetate
    "r_1808",  # glycerol
    "r_2033",  # pyruvate
    "r_2056",  # succinate
})

# Boundary reactions created by the candidate extension: closing them would stop
# the route carrying any flux, so they are always kept open.
CANDIDATE_EXCHANGES = ("DM_s_1507", ROUTE_INTERMEDIATE_DEMAND, TARGET_DEMAND)


def restrict_to_physiological_exchanges(model, *, keep=PHYSIOLOGICAL_EXCHANGES) -> int:
    """Close every boundary reaction outside the documented set. Returns the count.

    See PHYSIOLOGICAL_EXCHANGES for why: an open engineered exchange set lets a
    redox question be answered with a secretion no strain performs.
    """
    allowed = set(keep) | {BIOMASS_ID, *CANDIDATE_EXCHANGES}
    closed = 0
    for reaction in model.reactions:
        if reaction.boundary and reaction.id not in allowed:
            reaction.bounds = (0.0, 0.0)
            closed += 1
    return closed


def oxygen_species(model) -> list[str]:
    """Every compartmental O2 species in the model, by formula rather than by id.

    This model carries six of them (cytosol, ER, extracellular, mitochondrion,
    peroxisome, ER membrane) and the recycling that matters runs through several,
    so a check that looks at one id sees nothing.
    """
    return [
        metabolite.id
        for metabolite in model.metabolites
        if str(metabolite.formula).strip().upper() == "O2"
    ]


def _internal_oxygen_cycling(
    model, solution, *, tolerance: float = 1e-6
) -> dict[str, Any]:
    """Flux that makes oxygen inside the model instead of taking it up.

    A minimum-oxygen answer can come out at zero for a reason that has nothing to
    do with the strain. The oxidases hand their electrons to O2 and leave hydrogen
    peroxide; catalase converts peroxide back to water and O2; O2 transport moves
    that regenerated oxygen between compartments. Whether that loop is net
    positive is the LP's choice, not a constraint, so a flux distribution exists
    whose net external oxygen uptake is zero while the oxidases still run.
    Closing the exchanges cannot remove it, because the cycle is internal.
    Reporting the cycling fluxes is the difference between "this route needs no
    aeration" and "this model can be made to say so".
    """
    species = set(oxygen_species(model))
    producers = []
    for reaction in model.reactions:
        if reaction.boundary:
            continue
        flux = float(solution.fluxes[reaction.id])
        if abs(flux) <= tolerance:
            continue
        produced = sum(
            coefficient
            for metabolite, coefficient in reaction.metabolites.items()
            if metabolite.id in species and coefficient > 0
        )
        if produced * flux > tolerance:
            producers.append(
                {
                    "reaction": reaction.id,
                    "name": reaction.name,
                    "flux": flux,
                    "oxygen_produced": produced * flux,
                }
            )
    producers.sort(key=lambda entry: -entry["oxygen_produced"])
    return {
        "occurring": bool(producers),
        "producers": producers[:8],
        "consequence": (
            "oxygen is regenerated inside the model, so a zero or near-zero "
            "minimum oxygen uptake is a property of this reaction set rather than "
            "evidence that the process needs no aeration"
        ),
    }


def supply_requirement(
    target_config_path: str | Path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    target_titers_mg_l: tuple[float, ...] = (QLU_MEASURED_TITER_MG_L, 1000.0, 5000.0),
    biomass_gdw_l: float = 20.0,
    hours: float = 10.0,
    oxygen_levels: tuple[float, ...] = (100.0, 250.0, 500.0, 750.0, 1000.0, 1500.0),
    biomass_floor_fraction: float = BIOMASS_FLOOR_FRACTION,
    physiological_exports: bool = True,
) -> dict[str, Any]:
    """The carbon and oxygen bill for a titer, and which of the two is binding.

    ``capacity_to_titer`` converts a flux into a titer; this runs the conversion
    the other way, because a titer target is what a process has, not what the
    model produces. Each target becomes a required specific flux, and two bills
    are computed: the least glucose that reaches it (oxygen free) and the least
    oxygen that reaches it (carbon free). Both are LP optima, so neither reports
    a search grid's resolution as the answer.

    The glucose bill is then re-solved under each aeration ceiling. The penalty
    against the oxygen-free bill is what says whether aeration is limiting: a
    ceiling that costs no extra carbon is not a constraint, however low it looks.
    Reading only the aeration rows cannot answer that, because a ceiling can be
    slack for reasons that have nothing to do with the product.

    Two numbers describe the carbon side, and they differ by a lot at low titer:
    the average bill is dominated by the biomass floor's fixed cost, while the
    marginal bill is what an extra gram of product costs. Quoting the average
    alone makes a low-titer process look inefficient in a way that scaling up does
    not fix. Byproducts are reported as what must be secreted at the optimum
    (a second LP with the carbon fixed), because the carbon optimum is degenerate
    and the solver's arbitrary pick would otherwise be read as chemistry.
    """
    target_config = load_target_config(target_config_path)
    base = load_model(MODEL_ID)
    build_candidate_extension(base, target_config)

    # The biomass floor is taken at the module's reference medium, not at the
    # generous carbon ceiling used below: the model's growth yield is high enough
    # that a ceiling of 1000 units implies a growth rate near 50/h, and a floor
    # derived from that would demand more carbon than the product does.
    reference = base.copy()
    route_context(route, pox1_closed=pox1_closed)(reference)
    reference_growth = reference.optimize()
    if reference_growth.status != "optimal":
        return {
            "status": "reference medium not optimal",
            "error": str(reference_growth.status),
            "claim_boundary": CLAIM_BOUNDARY,
        }
    floor = biomass_floor_fraction * float(reference_growth.objective_value)

    closed_exchanges = 0

    def scenarise(oxygen: float):
        """A copy under the route scenario, with the supply bounds set for a bill LP."""
        model = base.copy()
        route_context(
            route,
            pox1_closed=pox1_closed,
            glucose=GLUCOSE_CEILING,
            oxygen=oxygen,
        )(model)
        if physiological_exports:
            nonlocal closed_exchanges
            closed_exchanges = restrict_to_physiological_exchanges(model)
        model.reactions.get_by_id(BIOMASS_ID).lower_bound = floor
        return model

    supply_ids = {"glucose": "r_1714", "oxygen": "r_1992"}
    supply_ceilings = {"glucose": GLUCOSE_CEILING, "oxygen": OXYGEN_CEILING}
    byproduct_ids = [
        reaction_id
        for name, reaction_id in OVERFLOW_EXCHANGES.items()
        if name != "carbon_dioxide"
    ]

    def ladder(required_flux: float, minimise: str, oxygen_ceiling: float) -> dict[str, Any]:
        """Minimum of one supply, the other at that minimum, then the byproducts.

        Both uptakes are capped at zero and the minimised one is pinned after the
        first solve, because an optimum is a face rather than a point: reading the
        second supply off the first LP would report the solver's arbitrary pick
        within that face, and a demand for a supply can even be met by excreting
        it if the direction is left open. The third solve exists because the
        byproduct fluxes are degenerate too, and a forced byproduct is a finding
        while an arbitrary one is noise.

        Maximising an exchange flux is the same as minimising its uptake, because
        uptake is the negative flux.
        """
        first, second = (minimise, "oxygen" if minimise == "glucose" else "glucose")
        model = scenarise(oxygen_ceiling)
        with model:
            headroom = model.optimize()
            if headroom.status != "optimal":
                return {"status": f"growth not optimal: {headroom.status}"}
            max_growth = float(headroom.objective_value)
            if max_growth < floor:
                return {
                    "status": "floor_not_attainable",
                    "max_growth_at_oxygen": max_growth,
                    "note": (
                        "this aeration cannot hold the reference biomass floor, so "
                        "the bill is undefined here rather than large"
                    ),
                }

            model.reactions.get_by_id(TARGET_DEMAND).lower_bound = required_flux
            for name, reaction_id in supply_ids.items():
                reaction = model.reactions.get_by_id(reaction_id)
                reaction.upper_bound = 0.0
                reaction.lower_bound = -supply_ceilings[name]

            primary = model.reactions.get_by_id(supply_ids[first])
            model.objective = primary
            minimum = model.optimize()
            if minimum.status != "optimal":
                return {
                    "status": minimum.status,
                    "max_growth_at_oxygen": max_growth,
                    "note": (
                        "the target flux is not reachable under this supply "
                        "combination, as opposed to being reachable and expensive"
                    ),
                }
            primary_uptake = -float(minimum.fluxes[supply_ids[first]])
            # Both bounds, not just one: the flux is negative for uptake, so
            # pinning only the lower bound would forbid taking up MORE than the
            # minimum while allowing LESS, and the next LP would quietly relax
            # away from the optimum it is supposed to be measured against.
            primary.bounds = (-primary_uptake, -primary_uptake)

            secondary = model.reactions.get_by_id(supply_ids[second])
            model.objective = secondary
            paired = model.optimize()
            secondary_uptake = None
            if paired.status == "optimal":
                secondary_uptake = -float(paired.fluxes[supply_ids[second]])
                secondary.bounds = (-secondary_uptake, -secondary_uptake)

            model.objective = {
                model.reactions.get_by_id(reaction_id): -1.0
                for reaction_id in byproduct_ids
            }
            clean = model.optimize()
            uptakes = {first: primary_uptake, second: secondary_uptake}
            recycled_oxygen = None
            if clean.status == "optimal":
                for name in uptakes:
                    uptakes[name] = -float(clean.fluxes[supply_ids[name]])
                if minimise == "oxygen":
                    recycled_oxygen = _internal_oxygen_cycling(model, clean)
            return {
                "status": "optimal",
                "max_growth_at_oxygen": max_growth,
                "glucose_uptake": uptakes["glucose"],
                "oxygen_uptake": uptakes["oxygen"],
                "ceiling_binding": (
                    uptakes[first] >= supply_ceilings[first] - 1e-9
                ),
                "growth": float(clean.fluxes[BIOMASS_ID]) if clean.status == "optimal" else None,
                "internal_oxygen_cycling": recycled_oxygen,
                "forced_byproducts": (
                    {
                        name: float(clean.fluxes[reaction_id])
                        for name, reaction_id in OVERFLOW_EXCHANGES.items()
                    }
                    if clean.status == "optimal"
                    else None
                ),
                "carbon_yield": (
                    required_flux * 10.0 / (uptakes["glucose"] * 6.0)
                    if uptakes["glucose"]
                    else None
                ),
            }

    targets: dict[str, Any] = {}
    for titer in target_titers_mg_l:
        required_flux = titer / (TARGET_MOLECULAR_WEIGHT * biomass_gdw_l * hours)
        product_g = titer / 1000.0
        carbon_bill = ladder(required_flux, "glucose", OXYGEN_CEILING)

        rows: list[dict[str, Any]] = []
        if carbon_bill.get("status") == "optimal":
            best_uptake = carbon_bill["glucose_uptake"]
            glucose_g = (
                best_uptake * biomass_gdw_l * hours * GLUCOSE_MOLECULAR_WEIGHT / 1000.0
            )
            for oxygen in oxygen_levels:
                entry = ladder(required_flux, "glucose", oxygen)
                if entry.get("status") != "optimal":
                    rows.append({"oxygen": oxygen, **entry})
                    continue
                carbon_penalty = entry["glucose_uptake"] / best_uptake - 1.0
                rows.append(
                    {
                        "oxygen": oxygen,
                        "glucose_uptake": entry["glucose_uptake"],
                        "oxygen_uptake": entry["oxygen_uptake"],
                        "glucose_g_per_g_product": (
                            entry["glucose_uptake"] * biomass_gdw_l * hours
                            * GLUCOSE_MOLECULAR_WEIGHT / 1000.0 / product_g
                        ),
                        "carbon_penalty_vs_free_oxygen": carbon_penalty,
                        # a cap that costs nothing is a cap that is not binding
                        "aeration_limited": carbon_penalty > 0.01,
                        "forced_byproducts": entry["forced_byproducts"],
                    }
                )

            # One more gram of product at this operating point costs the glucose
            # of a small increment, not the average bill: the fixed cost of the
            # biomass floor dominates the average, so the two differ sharply and
            # quoting the average alone misrepresents what scaling up costs.
            step = 0.05 * required_flux
            marginal = ladder(required_flux + step, "glucose", OXYGEN_CEILING)
            marginal_glucose_per_g = None
            if marginal.get("status") == "optimal":
                extra_product_g = (
                    step * TARGET_MOLECULAR_WEIGHT * biomass_gdw_l * hours / 1000.0
                )
                extra_glucose_g = (
                    (marginal["glucose_uptake"] - best_uptake) * biomass_gdw_l * hours
                    * GLUCOSE_MOLECULAR_WEIGHT / 1000.0
                )
                if extra_product_g:
                    marginal_glucose_per_g = extra_glucose_g / extra_product_g
            # The other corner of the trade-off: the least aeration that reaches
            # the titer, and the carbon bill at that aeration. The two corners
            # together are the answer to "which supply is binding".
            oxygen_bill = ladder(required_flux, "oxygen", OXYGEN_CEILING)
            targets[f"titer_{titer:g}_mg_l"] = {
                "titer_mg_l": titer,
                "required_specific_flux_mmol_gdw_h": required_flux,
                "average_glucose_g_per_g_product": (
                    glucose_g / product_g if product_g else None
                ),
                "marginal_glucose_g_per_g_product": marginal_glucose_per_g,
                "marginal_step_fraction": 0.05,
                "min_glucose_uptake": best_uptake,
                "carbon_bill_at_min_oxygen": oxygen_bill,
                "min_oxygen_uptake": oxygen_bill.get("oxygen_uptake"),
                "carbon_cost_of_minimum_aeration": (
                    oxygen_bill["glucose_uptake"] / best_uptake - 1.0
                    if oxygen_bill.get("status") == "optimal" and best_uptake
                    else None
                ),
                "rows": rows,
                "bill_is_flat_in_oxygen": (
                    all(not row.get("aeration_limited") for row in rows if "aeration_limited" in row)
                    if rows
                    else None
                ),
            }
        else:
            targets[f"titer_{titer:g}_mg_l"] = {
                "titer_mg_l": titer,
                "required_specific_flux_mmol_gdw_h": required_flux,
                "status": carbon_bill.get("status"),
                "note": "no carbon bill exists because the target is not reachable",
                "rows": [],
            }

    # The oxygen column is only worth reading if the model cannot make its own
    # oxygen. It can: see _internal_oxygen_cycling. Saying so once, at the top,
    # is the difference between a table that looks informative and one that is.
    cycling_entries = [
        (key, entry["carbon_bill_at_min_oxygen"]["internal_oxygen_cycling"])
        for key, entry in targets.items()
        if isinstance(entry.get("carbon_bill_at_min_oxygen"), dict)
        and isinstance(
            entry["carbon_bill_at_min_oxygen"].get("internal_oxygen_cycling"), dict
        )
        and entry["carbon_bill_at_min_oxygen"]["internal_oxygen_cycling"].get("occurring")
    ]
    oxygen_interpretable = not cycling_entries
    return {
        "status": "reported",
        "route": route,
        "pox1_closed": pox1_closed,
        "physiological_exports_only": bool(physiological_exports),
        "exchanges_closed": closed_exchanges,
        "oxygen_column_interpretable": oxygen_interpretable,
        "oxygen_column_note": (
            "aeration never binds anywhere in this table, and the reason is not that "
            "the route is insensitive to oxygen: this model regenerates oxygen "
            "internally, so no ceiling on the oxygen exchange can bind. The carbon "
            "bill is unaffected and remains the usable result; the aeration columns "
            "should be read as unavailable rather than as encouraging. The offending "
            "reactions are listed under internal_oxygen_cycling.producers for each "
            "target, and fixing them (capping catalase at its peroxide input, or "
            "making oxygen transport inward-only) is a prerequisite for asking any "
            "aeration question of this model."
            if not oxygen_interpretable
            else "no internal oxygen regeneration was detected at these optima"
        ),
        "assumptions": {
            "biomass_gdw_l": biomass_gdw_l,
            "hours": hours,
            "biomass_floor": floor,
            "biomass_floor_fraction": biomass_floor_fraction,
            "biomass_floor_medium": (
                "reference medium (the model's default glucose and oxygen), so the "
                "floor is one number across the whole table"
            ),
            "glucose_ceiling": GLUCOSE_CEILING,
            "titer_to_flux": (
                "titer / (molecular weight x biomass x hours): the whole biomass "
                "pool is assumed to produce at the model's maximum specific rate"
            ),
        },
        "targets": targets,
        "interpretation": (
            "The glucose bill is the carbon cost of the titer, not a productivity "
            "forecast: it assumes the culture holds the model's maximum specific "
            "flux at the stated biomass and time. The average and marginal bills "
            "answer different questions: the average is what the process consumes "
            "per gram at that titer, the marginal is what the next gram costs, and "
            "at low titer the average is dominated by the biomass floor, so a large "
            "gap between them means the process is paying fixed costs rather than "
            "an inefficient pathway. min_oxygen_bill is the aeration the titer "
            "needs when carbon is free: if it sits below every tested ceiling, "
            "aeration is not the constraint and the carbon bill is the whole story. "
            "A row whose carbon_penalty_vs_free_oxygen is zero is a ceiling that "
            "costs nothing, which is not the same as a ceiling that is high."
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }


def analyze(
    target_config_path: str | Path,
    *,
    pox1_closed: bool = False,
    growth_floor_fraction: float = BIOMASS_FLOOR_FRACTION,
) -> dict[str, Any]:
    """Maximise target flux per route and report the cofactor view alongside.

    POX1 is left OPEN by default because that is the configuration the external
    project's documented anchor was produced with; closing it is reported as a
    separate effect rather than folded into the baseline.

    Both routes are scored with the other one blocked, so each capacity is that
    route's own. At the LP level the two come out identical because neither is
    kcat-constrained yet; the routes separate only once enzyme capacity binds.
    """
    model = load_model(MODEL_ID)
    target_config = load_target_config(target_config_path)
    build_candidate_extension(model, target_config)

    routes = _route_reactions(model)
    if not routes:
        raise ValueError(
            "no documented route is present after extending the model; check the "
            f"target configuration at {target_config_path}"
        )

    baseline = model.optimize()
    if baseline.status != "optimal":
        raise RuntimeError(f"baseline is not optimal: {baseline.status}")
    biomass_floor = growth_floor_fraction * float(baseline.objective_value)

    report: dict[str, Any] = {
        "model_id": MODEL_ID,
        "target_demand": TARGET_DEMAND,
        "baseline_growth": float(baseline.objective_value),
        "biomass_floor": biomass_floor,
        "pox1_closed": pox1_closed,
        "routes": {},
        "claim_boundary": CLAIM_BOUNDARY,
    }

    for route_name, reaction_ids in routes.items():
        with model:
            model.reactions.get_by_id(BIOMASS_ID).lower_bound = biomass_floor
            if pox1_closed:
                model.reactions.get_by_id(POX1_REACTION).bounds = (0.0, 0.0)
            # close the other route so the two are compared in isolation
            for other_name, other_ids in routes.items():
                if other_name == route_name:
                    continue
                for reaction_id in other_ids:
                    model.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)

            model.objective = TARGET_DEMAND
            solution = model.optimize()
            entry: dict[str, Any] = {"status": solution.status}
            if solution.status == "optimal":
                capacity = float(solution.objective_value)
                entry["target_capacity"] = capacity
                entry["route_reactions"] = reaction_ids
                entry["route_fluxes"] = {
                    reaction_id: float(solution.fluxes[reaction_id])
                    for reaction_id in reaction_ids
                }
                entry["pox1_flux"] = float(solution.fluxes[POX1_REACTION])
                spec = MODEL_REGISTRY[MODEL_ID]
                turnover = cofactor_turnover(
                    model, solution, spec["cofactor_metabolites"]
                )
                # flavins are useless here (see the module docstring), so the
                # degradation side is read from oxygen and peroxide instead
                burden = oxidative_burden(model, solution, spec["oxidant_metabolites"])
                # the route is peroxisomal, so the cytosolic-only figures describe
                # a pool it never touches; report the split explicitly
                entry["cofactor_by_compartment"] = cofactor_turnover_by_compartment(
                    model, solution, spec["cofactor_metabolites"]
                )
                entry["oxidative_burden"] = {
                    label: {
                        "produced": burden[label].get("produced"),
                        "consumed": burden[label].get("consumed"),
                        "top_producers": burden[label].get("top_producers", []),
                        "top_consumers": burden[label].get("top_consumers", []),
                    }
                    for label in ("h2o2", "o2")
                }
                entry["cofactor_turnover"] = {
                    label: {
                        "status": turnover[label].get("status"),
                        "metabolite_id": turnover[label].get("metabolite_id"),
                        "gross_turnover": turnover[label].get("gross_turnover"),
                        "top_producers": turnover[label].get("top_producers", []),
                    }
                    for label in FEASIBILITY_COFACTORS
                }
                entry["core_cofactor_turnover"] = {
                    label: turnover[label].get("gross_turnover")
                    for label in CORE_COFACTORS
                }
                # enzyme investment: the flux a route must sustain divided by the
                # target capacity, per unit of route flux, is unity by
                # construction; what matters is comparing the two routes'
                # capacities at the same biomass floor.
                entry["capacity_ratio_to_route_a"] = None  # filled below
            report["routes"][route_name] = entry

    rated = {
        name: entry["target_capacity"]
        for name, entry in report["routes"].items()
        if entry.get("target_capacity")
    }
    if "route_A_free_acid" in rated:
        reference = rated["route_A_free_acid"]
        for name, entry in report["routes"].items():
            if entry.get("target_capacity"):
                entry["capacity_ratio_to_route_a"] = entry["target_capacity"] / reference
    report["capacity_summary"] = rated

    # Self-check against the external report's stated anchor, and the effect of
    # closing the degradation entrance, computed on the same protocol.
    if "route_A_free_acid" in rated:
        measured = rated["route_A_free_acid"]
        report["anchor_check"] = {
            "documented": DOCUMENTED_ANCHOR,
            "measured": measured,
            "matches": abs(measured - DOCUMENTED_ANCHOR) <= ANCHOR_TOLERANCE,
        }

    if not pox1_closed:
        closed = analyze(
            target_config_path,
            pox1_closed=True,
            growth_floor_fraction=growth_floor_fraction,
        )
        closed_rated = closed.get("capacity_summary", {})
        report["pox1_effect"] = {
            name: {
                "open": open_capacity,
                "closed": closed_rated.get(name),
                "change_fraction": (
                    (closed_rated[name] - open_capacity) / open_capacity
                    if closed_rated.get(name) else None
                ),
            }
            for name, open_capacity in rated.items()
        }
        report["pox1_interpretation"] = (
            "A negative change means closing the acyl-CoA oxidase entrance LOWERS "
            "the achievable target flux in this model at this medium, which is the "
            "opposite of the naive 'block the degradation sink' expectation. The "
            "reaction contributes to redox balance that the biomass floor requires, "
            "so removing it costs more than it saves. This is a falsifiable model "
            "prediction for a POX1 knockout experiment."
        )
    return report


def minimum_enzyme_for_flux(kcat_s: float, flux: float) -> float:
    """Enzyme load (mmol gDW-1) needed to sustain a flux, ub = kcat * E * 3600."""
    if kcat_s <= 0:
        raise ValueError(f"kcat must be positive, got {kcat_s}")
    return flux / (kcat_s * SECONDS_PER_HOUR)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-config",
        default=(
            "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"
            "/10_generic_target_workflow/examples/target_workflow_10h2da_reference.json"
        ),
    )
    parser.add_argument("--close-pox1", action="store_true",
                        help="block the acyl-CoA oxidase entrance (POX1 knockout)")
    parser.add_argument("--output", default=None)
    parser.add_argument("--marginal", action="store_true",
                        help="report the marginal cofactor cost per unit product")
    parser.add_argument("--route", default="A", choices=("A", "B"))
    parser.add_argument("--supply-scan", action="store_true",
                        help="scan the peroxisomal NADPH producer's capacity ceiling")
    parser.add_argument("--shuttle-compare", action="store_true",
                        help="compare the route with and without reducing-power shuttles")
    parser.add_argument("--carbon-check", action="store_true",
                        help="test whether capacity is carbon-limited or internally limited")
    parser.add_argument("--shuttle-cost-scan", action="store_true",
                        help="price the peroxisomal shuttle in ATP and cap its rate")
    parser.add_argument("--supply-requirement", action="store_true",
                        help="carbon and oxygen bill for a target titer, with overflow markers")
    parser.add_argument("--prosthetic-ledger", action="store_true",
                        help="the enzyme-bound cofactor bill (FAD/FMN/heme) that FBA cannot see")
    parser.add_argument("--titer", action="store_true",
                        help="convert the anchor capacity to an implied titer")
    parser.add_argument("--biomass-gdw-l", type=float, default=20.0)
    parser.add_argument("--hours", type=float, default=10.0)
    args = parser.parse_args()

    if args.titer:
        import json as _json

        anchor = analyze(args.target_config, pox1_closed=args.close_pox1)
        capacity = anchor["routes"][f"route_{args.route}_" + (
            "free_acid" if args.route == "A" else "coa_thioester")]["target_capacity"]
        print(_json.dumps(
            capacity_to_titer(capacity, biomass_gdw_l=args.biomass_gdw_l, hours=args.hours),
            indent=2,
        ))
        return 0
    if args.carbon_check:
        report = carbon_limitation_check(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    elif args.shuttle_cost_scan:
        report = shuttle_cost_scan(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    elif args.supply_requirement:
        report = supply_requirement(
            args.target_config,
            route=args.route,
            pox1_closed=args.close_pox1,
            biomass_gdw_l=args.biomass_gdw_l,
            hours=args.hours,
        )
    elif args.prosthetic_ledger:
        from .prosthetic_ledger import route_ledger

        report = route_ledger(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    elif args.shuttle_compare:
        report = shuttle_comparison(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    elif args.supply_scan:
        report = peroxisomal_nadph_supply_scan(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    elif args.marginal:
        report = marginal_route_cost(
            args.target_config, route=args.route, pox1_closed=args.close_pox1
        )
    else:
        report = analyze(args.target_config, pox1_closed=args.close_pox1)
    text = json.dumps(report, indent=2)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
