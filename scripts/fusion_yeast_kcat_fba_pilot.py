"""GECKO-style pilot: kcat-constrained FBA on the yeast 10H2DA reference case.

This is the first closed loop of the BioCandidateRanker / Yeast-MetaTwin fusion:

1. reproduce the external project's 10H2DA candidate-route FBA anchor
   (biomass floor 10%, maximize DM_CAND_10H2DA_P; the Sep 2 report states
   0.274011 model flux units on the default medium);
2. re-solve with enzyme-capacity upper bounds of the GECKO form
   ub = kcat [s^-1] x E [mmol/gDW] x 3600 [s/h], where the kcat values are
   BioCandidateRanker log10(kcat) predictions for the exact enzyme-substrate
   pairs (uncalibrated), and E is an assumed enzyme level that is scanned;
3. report the minimum enzyme load each step needs to sustain the unconstrained
   optimum, so the model-flux number gains an explicit enzyme-investment cost.

Claim boundary: everything here is simulation. kcat values are uncalibrated
model predictions, not measurements; enzyme levels are assumed, not measured;
FBA fluxes are model values, not experimental rates. The output is an
engineering pilot that attaches a quantitative enzyme cost to the pathway
prediction, not a validated production forecast.

The external project is opened read-only; the candidate reaction definitions
are re-constructed here from its target workflow config, whose SHA256 is
recorded in the output for provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import cobra

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from biocandidate import fba_context  # noqa: E402

MODEL_ID = "yeast_metatwin"
BIOMASS_ID = "r_2111"
TARGET_DEMAND = "DM_CAND_10H2DA_P"
BIOMASS_FLOOR_FRACTION = 0.1
SECONDS_PER_HOUR = 3600.0
AVG_RESIDUE_MASS_DA = 110.0

EXTERNAL_TARGET_CONFIG = Path(
    "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"
    "/10_generic_target_workflow/examples/target_workflow_10h2da_reference.json"
)
PREDICTIONS = ROOT / "artifacts/10h2da-yeast-prediction/task1_10h2da_pathway-predictions.json"
SEQUENCES = ROOT / "artifacts/10h2da-yeast-prediction/sequences"

# enzyme capacity scenarios: reaction id -> (enzyme label, candidate_id of the
# BioCandidateRanker prediction that supplies the kcat)
ROUTE_A = {
    "CAND_T2DEC_THIOESTERASE_P": ("TES1/PTE1", "TES1_P41903_dec2enoylCoA"),
    "CAND_T2DEC_OMEGA_HYDROXYLASE_P": ("CYP52M1", "CYP52M1_B8QHP1_dec2enoicAcid_routeA"),
}
ROUTE_B = {
    "CAND_T2DEC_COA_OMEGA_HYDROXYLASE_P": ("CYP52M1", "CYP52M1_B8QHP1_dec2enoylCoA_routeB"),
    "CAND_10H2DA_COA_THIOESTERASE_P": ("TES1/PTE1", "TES1_P41903_10HOdec2enoylCoA"),
}
POX1_GATE = {
    "r_0120": ("POX1", "POX1_P13711_decanoylCoA"),
}

ENZYME_LEVELS_MMOL_GDW = (1e-5, 3e-5, 1e-4, 1e-3, 1e-2)

CLAIM_BOUNDARY = (
    "Simulation pilot. kcat values are uncalibrated BioCandidateRanker predictions; "
    "enzyme levels are assumed, not measured; FBA fluxes are model values. Nothing "
    "here is experimental validation or a production forecast."
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_kcat_predictions() -> dict[str, float]:
    payload = json.loads(PREDICTIONS.read_text(encoding="utf-8"))
    kcat: dict[str, float] = {}
    for row in payload["predictions"]:
        task = row["tasks"]["log10_kcat"]
        if task.get("mean") is None:
            raise SystemExit(f"candidate {row['candidate_id']} has no kcat prediction")
        kcat[row["candidate_id"]] = float(10.0 ** task["mean"])
    return kcat


def load_sequence_lengths() -> dict[str, int]:
    lengths = {}
    for fasta in sorted(SEQUENCES.glob("*.fasta")):
        lines = fasta.read_text().splitlines()
        lengths[fasta.stem] = sum(len(line.strip()) for line in lines if not line.startswith(">"))
    return lengths


def build_candidate_extension(model, target_config: dict) -> None:
    """Re-construct the external project's candidate layer (idempotent)."""
    for row in target_config["compounds"]:
        met_id = row["model_metabolite_id"]
        if met_id in model.metabolites:
            continue
        model.add_metabolites([
            cobra.Metabolite(
                met_id,
                formula=row["formula"],
                charge=int(row["charge"]),
                name=row["name"],
                compartment="p",
            )
        ])
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
    demands = {
        "DM_s_1507": {"s_1507": -1},
        "DM_CAND_T2DEC_P": {"cand_t2dec_p": -1},
        TARGET_DEMAND: {"cand_10h2da_p": -1},
    }
    for demand_id, stoichiometry in demands.items():
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


def solve_target(model) -> dict[str, float | str]:
    """Biomass-floor + target-maximization protocol of the external validation."""
    work = model.copy()
    growth_max = float(work.slim_optimize())
    if growth_max <= 0:
        return {"status": "infeasible", "target_flux": 0.0, "growth": 0.0}
    biomass = work.reactions.get_by_id(BIOMASS_ID)
    work.add_cons_vars(work.problem.Constraint(
        biomass.flux_expression, lb=growth_max * BIOMASS_FLOOR_FRACTION,
        name="10h2da_pilot_biomass_floor",
    ))
    work.objective = TARGET_DEMAND
    solution = work.optimize()
    if solution.status != "optimal":
        return {"status": solution.status, "target_flux": 0.0, "growth": 0.0}
    return {
        "status": "optimal",
        "target_flux": float(solution.fluxes[TARGET_DEMAND]),
        "growth": float(solution.fluxes[BIOMASS_ID]),
    }


def solve_kcat_constrained(model, kcat: dict[str, float], constraint_map: dict[str, dict],
                           enzyme_level: float, blocked: tuple[str, ...] = ()) -> dict:
    work = model.copy()
    for reaction_id in blocked:
        work.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
    applied = []
    for reaction_id, (enzyme, candidate_id) in constraint_map.items():
        ub = kcat[candidate_id] * enzyme_level * SECONDS_PER_HOUR
        work.reactions.get_by_id(reaction_id).upper_bound = ub
        applied.append({
            "reaction_id": reaction_id,
            "enzyme": enzyme,
            "kcat_s_inv": round(kcat[candidate_id], 4),
            "enzyme_mmol_gdw": enzyme_level,
            "upper_bound_mmol_gdw_h": round(ub, 6),
        })
    result = solve_target(work)
    result["constraints"] = applied
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/fusion-yeast-pilot")
    args = parser.parse_args()
    out_dir = ROOT / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    kcat = load_kcat_predictions()
    lengths = load_sequence_lengths()

    target_config = json.loads(EXTERNAL_TARGET_CONFIG.read_text(encoding="utf-8"))
    model = fba_context.load_model(MODEL_ID)  # SHA256-verified tracked copy
    build_candidate_extension(model, target_config)

    results: dict[str, object] = {
        "pilot": "fusion_yeast_kcat_fba_pilot",
        "date": "2026-10-04",
        "claim_boundary": CLAIM_BOUNDARY,
        "model_id": MODEL_ID,
        "model_sha256": sha256(ROOT / "data/Yeast-MetaTwin.yml"),
        "external_target_config": {
            "path": str(EXTERNAL_TARGET_CONFIG),
            "sha256": sha256(EXTERNAL_TARGET_CONFIG),
        },
        "kcat_predictions": {
            "path": str(PREDICTIONS),
            "sha256": sha256(PREDICTIONS),
            "checkpoint": json.loads(PREDICTIONS.read_text(encoding="utf-8"))["checkpoint"],
        },
        "medium": "model file defaults (glucose cap 1.0, oxygen cap 1000), matching the external anchor",
        "biomass_floor_fraction": BIOMASS_FLOOR_FRACTION,
    }

    anchor = solve_target(model)
    results["unconstrained_anchor"] = anchor
    print(f"anchor (no kcat constraints): target={anchor['target_flux']:.6f} "
          f"status={anchor['status']}")

    # minimum enzyme load that still sustains the unconstrained optimum
    flux = float(anchor["target_flux"])
    enzyme_cost = {}
    for constraint_map in (ROUTE_A, ROUTE_B, POX1_GATE):
        for reaction_id, (enzyme, candidate_id) in constraint_map.items():
            e_min = flux / (kcat[candidate_id] * SECONDS_PER_HOUR)  # mmol/gDW
            uniprot = candidate_id.split("_")[1]
            mw_g_mol = lengths[uniprot] * AVG_RESIDUE_MASS_DA  # ~110 Da per residue
            mg_gdw = e_min * 1e-3 * mw_g_mol * 1000.0  # mmol -> mol -> g
            enzyme_cost[reaction_id] = {
                "enzyme": enzyme,
                "uniprot": uniprot,
                "kcat_s_inv": round(kcat[candidate_id], 4),
                "min_enzyme_mmol_gdw": float(f"{e_min:.3e}"),
                "min_enzyme_mg_gdw": round(mg_gdw, 3),
            }
    results["enzyme_cost_to_sustain_anchor"] = enzyme_cost

    scenarios = {}
    for name, constraint_map, blocked in (
        ("route_A_free_acid", ROUTE_A, tuple(ROUTE_B)),
        ("route_B_coa_bound", ROUTE_B, tuple(ROUTE_A)),
        ("combined_with_pox1_gate", {**ROUTE_A, **ROUTE_B, **POX1_GATE}, ()),
    ):
        scans = []
        for level in ENZYME_LEVELS_MMOL_GDW:
            scan = solve_kcat_constrained(model, kcat, constraint_map, level, blocked)
            scans.append({"enzyme_level_mmol_gdw": level, **scan})
            print(f"{name} E={level:g}: target={scan['target_flux']:.6f}")
        scenarios[name] = scans
    results["kcat_constrained_scans"] = scenarios

    out_path = out_dir / "pilot_results.json"
    out_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"-> {out_path}")
    return 0 if anchor["status"] == "optimal" else 1


if __name__ == "__main__":
    raise SystemExit(main())
