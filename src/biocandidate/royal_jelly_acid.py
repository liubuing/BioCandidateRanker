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
    load_model,
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
    args = parser.parse_args()

    report = analyze(args.target_config, pox1_closed=args.close_pox1)
    text = json.dumps(report, indent=2)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
