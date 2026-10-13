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


def add_shuttle(model, shuttle_key: str) -> str:
    """Add one idealized shuttle reaction, idempotently. Returns its id."""
    import cobra

    definition = SHUTTLE_DEFINITIONS[shuttle_key]
    reaction_id = "SHUTTLE_" + shuttle_key.upper()
    if reaction_id in model.reactions:
        return reaction_id
    reaction = cobra.Reaction(reaction_id)
    reaction.name = definition["name"]
    reaction.lower_bound = -1000.0
    reaction.upper_bound = 1000.0
    model.add_reactions([reaction])
    reaction.add_metabolites(
        {model.metabolites.get_by_id(k): v for k, v in definition["stoichiometry"].items()}
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
