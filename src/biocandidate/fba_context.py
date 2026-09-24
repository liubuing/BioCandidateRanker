"""Genome-scale FBA context features for the multimodal ranker.

This module imports the technology stack of two E. coli metabolic-engineering
efforts and exposes it as governed model inputs:

* the cofactor-engineering paper by Wang et al. (D-pantothenic acid production):
  the iML1515 model, FBA with a biomass objective, knockouts and overexpression
  simulated through reaction bounds, medium glucose 15 / oxygen 20 mmol/gDW/h,
  and optional flux variability analysis;
* the user's arginine / carbon-fixation project (05_FBA 融合FBA计算.py): the
  iJO1366 model plus eight synthetic reactions (CCR_EMA, NOG_F6P, PTXD, ...) and
  the NOG / phosphite / CO2 scenario caps.

Both produce the same contract: an 8-wide fba_context vector whose reaction
identifiers are declared in FEATURE_IDS, with hash-identified model files and a
condition id that binds the vector to the exact scenario that produced it.

Claim boundary: every value produced here is a simulation under the declared
model assumptions. Nothing here is an experimental flux measurement, and no
checkpoint shipped with this project was trained with flux labels, so attaching
these features to prediction is an engineering integration, not a scientifically
validated input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .data import FBA_FEATURE_SCHEMA_VERSION

MODEL_CONFIGS = {
    "iML1515": Path("configs/fba_context_model.json"),
    "iJO1366": Path("configs/fba_context_model_ijo1366.json"),
}

GLUCOSE_EXCHANGE = "EX_glc__D_e"
OXYGEN_EXCHANGE = "EX_o2_e"

# The paper monitors EMP, PPP, ED, and TCA fluxes plus redox and energy state;
# these eight reactions cover the same gates and match ModelConfig
# .fba_context_dim of every shipped checkpoint. Both registered models carry
# every reaction below except the biomass objective, which is per-model.
SHARED_FLUX_REACTIONS = ("PFK", "G6PDH2r", "EDD", "CS", "THD2pp", "ATPM")

MODEL_REGISTRY: dict[str, dict[str, Any]] = {
    "iML1515": {
        "objective_id": "BIOMASS_Ec_iML1515_core_75p37M",
        "citation": "Monk et al. 2017, Nature Microbiology 2:17094",
    },
    "iJO1366": {
        "objective_id": "BIOMASS_Ec_iJO1366_WT_53p95M",
        "citation": "Orth et al. 2011, Mol Syst Biol 7:530",
    },
}
for _spec in MODEL_REGISTRY.values():
    _spec["feature_reactions"] = (_spec["objective_id"], GLUCOSE_EXCHANGE) + SHARED_FLUX_REACTIONS
    _spec["feature_ids"] = tuple(f"flux:{r}" for r in _spec["feature_reactions"])

DEFAULT_GLUCOSE = 15.0
DEFAULT_OXYGEN = 20.0

# Scenario caps reproduced from the arginine fusion script (05_FBA).
FUSION_SCENARIO_CAPS = {
    "M2_NOG": {"nog": 3.0, "phosphite": 0.0},
    "M3_PtxD": {"nog": 0.0, "phosphite": 5.0},
    "M5_NOG_PtxD": {"nog": 3.0, "phosphite": 5.0},
}
FUSION_CO2_UPTAKE_CAP = 5.0

CLAIM_BOUNDARY = (
    "FBA context values are simulations under the declared model assumptions. "
    "They are not experimental flux measurements, and no shipped checkpoint was "
    "trained with flux labels, so attaching them to prediction is an engineering "
    "integration rather than a validated input."
)


def signed_log10(value: float) -> float:
    """Signed log10(1 + |v|): keeps sign and zero, bounds the dynamic range."""
    if not math.isfinite(value):
        raise ValueError(f"flux must be finite, got {value}")
    magnitude = math.log10(1.0 + abs(value))
    return math.copysign(magnitude, value)


@dataclass(frozen=True)
class Modification:
    """One genetic change, in the vocabulary of the paper.

    The paper realizes knockouts by setting reaction bounds to zero and
    overexpression by constraining a reaction to a proportion of its baseline
    flux; both are reproduced here.
    """

    kind: str  # "knockout" | "overexpression"
    target: str  # gene name/id for knockout, reaction id for overexpression
    factor: float = 1.5


# Paper-derived presets: model-agnostic, resolved against the loaded model.
BASE_PRESETS: dict[str, tuple[Modification, ...]] = {
    "baseline": (),
    "ptsG_knockout": (Modification("knockout", "ptsG"),),
    "thd2pp_overexpression": (Modification("overexpression", "THD2pp", 1.5),),
    "nadK_overexpression": (Modification("overexpression", "NADK", 1.5),),
    "gltA_overexpression": (Modification("overexpression", "CS", 1.5),),
    "paper_stack": (
        Modification("knockout", "ptsG"),
        Modification("overexpression", "THD2pp", 1.5),
        Modification("overexpression", "NADK", 1.5),
    ),
}

# Arginine-project presets: require the iJO1366 fusion pathway layer.
FUSION_PRESETS: dict[str, dict[str, Any]] = {
    "fusion_M2_NOG": {"model_id": "iJO1366", **FUSION_SCENARIO_CAPS["M2_NOG"]},
    "fusion_M3_PtxD": {"model_id": "iJO1366", **FUSION_SCENARIO_CAPS["M3_PtxD"]},
    "fusion_M5_NOG_PtxD": {"model_id": "iJO1366", **FUSION_SCENARIO_CAPS["M5_NOG_PtxD"]},
}

PROVENANCE = {
    "iML1515_presets": "Wang et al., Integrated cofactor-centric engineering and energy flux optimization enable high-efficient D-pantothenic acid production in Escherichia coli",
    "iJO1366_fusion": "user project 05_FBA 融合FBA计算.py: add_fusion_pathways, configure_model scenario caps",
}


def model_ids() -> list[str]:
    return sorted(MODEL_REGISTRY)


def presets_for(model_id: str) -> list[str]:
    _require_model(model_id)
    base = sorted(BASE_PRESETS)
    if model_id == "iJO1366":
        return base + sorted(FUSION_PRESETS)
    return base


def preset_requires_fusion(preset: str) -> bool:
    return preset in FUSION_PRESETS


def load_model_config(model_id: str, root: Path | None = None) -> dict[str, Any]:
    _require_model(model_id)
    path = (root or Path.cwd()) / MODEL_CONFIGS[model_id]
    return json.loads(path.read_text(encoding="utf-8"))


def _require_model(model_id: str) -> None:
    if model_id not in MODEL_REGISTRY:
        raise ValueError(f"unknown model {model_id!r}; known: {model_ids()}")


def verify_model_file(model_path: Path, expected: dict[str, Any]) -> None:
    digest = hashlib.sha256(model_path.read_bytes()).hexdigest()
    if digest != expected["sha256"]:
        raise ValueError(
            f"model file identity mismatch: {model_path} hashes to {digest[:16]}..., "
            f"config expects {expected['sha256'][:16]}..."
        )
    if model_path.stat().st_size != expected["size_bytes"]:
        raise ValueError(f"model file size mismatch: {model_path}")


def load_model(model_id: str, model_path: str | Path | None = None, *, verify: bool = True):
    """Load a registered model and verify it against the tracked config."""
    import cobra

    _require_model(model_id)
    config = load_model_config(model_id)
    resolved = Path(model_path) if model_path else Path(config["path"])
    if verify and resolved.is_file():
        verify_model_file(resolved, config)
    if resolved.suffix == ".json":
        model = cobra.io.load_json_model(str(resolved))
    else:
        model = cobra.io.read_sbml_model(str(resolved))
    return model


def add_fusion_pathways(model) -> list[str]:
    """Port of the arginine project's synthetic pathway layer.

    Five metabolites and eight reactions (CCR_EMA carboxylation, YciA_EMA
    thioesterase, EMAtex export, NOG_F6P non-oxidative glycolysis, phosphite
    transport and PTXD oxidation) are appended to iJO1366 exactly as defined in
    05_FBA/融合FBA计算.py. The layer is inert until a fusion preset caps its
    fluxes, so baseline behaviour is unchanged. Idempotent: identifiers that
    already exist are skipped, so a reused model object cannot be corrupted by
    duplicate entities with split mass balance.
    """
    import cobra

    existing_metabolites = {m.id for m in model.metabolites}
    existing_reactions = {r.id for r in model.reactions}

    metabolites = [
        cobra.Metabolite("ema_c", formula="C26H37N7O19P3S", charge=-5,
                         name="ethylmalonyl-CoA", compartment="c"),
        cobra.Metabolite("emaacid_c", formula="C5H6O4", charge=-2,
                         name="ethylmalonate", compartment="c"),
        cobra.Metabolite("emaacid_e", formula="C5H6O4", charge=-2,
                         name="ethylmalonate", compartment="e"),
        cobra.Metabolite("phite_c", formula="HO3P", charge=-2,
                         name="phosphite", compartment="c"),
        cobra.Metabolite("phite_e", formula="HO3P", charge=-2,
                         name="phosphite", compartment="e"),
    ]
    model.add_metabolites([m for m in metabolites if m.id not in existing_metabolites])

    definitions = {
        "CCR_EMA": "b2coa_c + co2_c + nadph_c --> ema_c + nadp_c",
        "YciA_EMA": "ema_c + h2o_c --> coa_c + emaacid_c + h_c",
        "EMAtex": "emaacid_c --> emaacid_e",
        "EX_emaacid_e": "emaacid_e -->",
        "NOG_F6P": "f6p_c + 2 pi_c --> 3 actp_c + 2 h2o_c",
        "PHITEtex": "phite_e --> phite_c",
        "EX_phite_e": "phite_e <=>",
        "PTXD": "h2o_c + nad_c + phite_c --> h_c + nadh_c + pi_c",
    }
    new_ids = [rid for rid in definitions if rid not in existing_reactions]
    reactions = []
    for reaction_id in new_ids:
        reaction = cobra.Reaction(reaction_id)
        reaction.name = reaction_id
        reaction.lower_bound = 0.0
        reaction.upper_bound = 1000.0
        reactions.append(reaction)
    model.add_reactions(reactions)
    for reaction_id in new_ids:
        model.reactions.get_by_id(reaction_id).build_reaction_from_string(
            definitions[reaction_id]
        )
    if "EX_phite_e" in model.reactions:
        model.reactions.EX_phite_e.lower_bound = -1000.0
    return list(definitions)


def _resolve_gene(model, name_or_id: str):
    gene = next(
        (g for g in model.genes if g.name == name_or_id or g.id == name_or_id), None
    )
    if gene is None:
        raise ValueError(f"gene {name_or_id!r} is not present in the model")
    return gene


def _set_medium(model, spec: dict[str, Any], glucose: float, oxygen: float) -> None:
    for exchange, bound in ((GLUCOSE_EXCHANGE, glucose), (OXYGEN_EXCHANGE, oxygen)):
        if bound < 0:
            raise ValueError(f"uptake must be non-negative, got {bound}")
        model.reactions.get_by_id(exchange).lower_bound = -float(bound)
    if spec.get("fusion_pathways"):
        model.reactions.EX_co2_e.lower_bound = -float(FUSION_CO2_UPTAKE_CAP)
        model.reactions.NOG_F6P.upper_bound = float(spec["nog"])
        model.reactions.EX_phite_e.lower_bound = -float(spec["phosphite"])
        if spec["phosphite"] == 0:
            model.reactions.PTXD.upper_bound = 0.0


def _apply_knockout(model, target: str) -> str:
    gene = _resolve_gene(model, target)
    gene.knock_out()
    return f"knockout:{gene.id}"


def _apply_overexpression(model, target: str, factor: float, baseline_fluxes) -> str:
    if factor <= 0:
        raise ValueError(f"overexpression factor must be positive, got {factor}")
    reaction = model.reactions.get_by_id(target)
    optimum = float(baseline_fluxes.get(target, 0.0))
    if optimum > 0:
        reaction.lower_bound = factor * optimum
    elif optimum < 0:
        reaction.upper_bound = factor * optimum
    # A zero-flux optimum leaves nothing to scale; record that instead of guessing.
    return f"overexpression:{target}@{factor:g}"


def run_fba_context(
    *,
    model_id: str = "iML1515",
    preset: str = "baseline",
    glucose: float = DEFAULT_GLUCOSE,
    oxygen: float = DEFAULT_OXYGEN,
    growth_min: float = 0.0,
    with_fva: bool = False,
    model=None,
    model_sha256: str | None = None,
) -> dict[str, Any]:
    """Run FBA and return the 8-wide governed context vector.

    Every run operates on its own copy of the model, so structural changes
    (the fusion layer) and bound changes can never leak between calls through
    a shared object; the caller may pass a cached base model safely.
    Fusion presets add the synthetic pathway layer and its scenario caps
    before solving.
    """
    if model_id not in MODEL_REGISTRY:
        raise ValueError(f"unknown model {model_id!r}; known: {model_ids()}")
    spec = MODEL_REGISTRY[model_id]
    known = presets_for(model_id)
    if preset not in known:
        raise ValueError(f"unknown preset {preset!r} for {model_id}; known: {known}")

    fusion_spec = FUSION_PRESETS.get(preset)
    modifications = BASE_PRESETS.get(preset, ())

    if model is None:
        model = load_model(model_id)
    solver_id = type(model.solver.interface).__name__
    scoped = model.copy()

    result: dict[str, Any] = {
        "model_id": model_id,
        "objective_id": spec["objective_id"],
        "preset": preset,
        "fusion_pathways_applied": fusion_spec is not None,
        "medium": {GLUCOSE_EXCHANGE: glucose, OXYGEN_EXCHANGE: oxygen},
        "solver_id": solver_id,
        "feature_ids": list(spec["feature_ids"]),
        "claim_boundary": CLAIM_BOUNDARY,
        "provenance": (
            PROVENANCE["iJO1366_fusion"] if fusion_spec else PROVENANCE["iML1515_presets"]
        ),
    }
    if model_sha256:
        result["model_sha256"] = model_sha256

    if fusion_spec:
        add_fusion_pathways(scoped)
    # The objective is part of the declared condition: never rely on the model
    # file's default, which for iJO1366 is the "core" biomass while this
    # contract follows the arginine project's "WT" biomass.
    scoped.objective = spec["objective_id"]
    _set_medium(scoped, fusion_spec or {}, glucose, oxygen)
    if growth_min > 0:
        scoped.reactions.get_by_id(spec["objective_id"]).lower_bound = growth_min

    baseline = scoped.optimize()
    if baseline.status != "optimal":
        result.update({"status": baseline.status, "growth": None})
        return result

    applied: list[str] = []
    for modification in modifications:
        if modification.kind == "knockout":
            applied.append(_apply_knockout(scoped, modification.target))
        elif modification.kind == "overexpression":
            applied.append(
                _apply_overexpression(
                    scoped, modification.target, modification.factor, baseline.fluxes
                )
            )
        else:
            raise ValueError(f"unknown modification kind {modification.kind!r}")
    result["applied_modifications"] = applied

    solution = scoped.optimize()
    if solution.status != "optimal":
        result.update({"status": solution.status, "growth": None})
        return result
    result["status"] = "optimal"
    result["growth"] = float(solution.objective_value)

    fluxes = {}
    for reaction in spec["feature_reactions"]:
        fluxes[reaction] = float(solution.fluxes[reaction])
    result["fluxes"] = fluxes
    result["features"] = [signed_log10(fluxes[r]) for r in spec["feature_reactions"]]

    condition = {
        "model_id": model_id,
        "preset": preset,
        "medium": result["medium"],
        "growth_min": growth_min,
        "fusion_spec": fusion_spec,
        "model_sha256": model_sha256,
    }
    result["condition_id"] = (
        "sha256:" + hashlib.sha256(
            json.dumps(condition, sort_keys=True).encode("utf-8")
        ).hexdigest()
    )
    result["metadata"] = {
        "schema_version": FBA_FEATURE_SCHEMA_VERSION,
        "feature_ids": list(spec["feature_ids"]),
        "model_id": model_id,
        "solver_id": solver_id,
        "objective_id": spec["objective_id"],
        "condition_id": result["condition_id"],
    }

    if with_fva:
        from cobra.flux_analysis import flux_variability_analysis

        fva = flux_variability_analysis(
            scoped, reaction_list=list(spec["feature_reactions"]),
            fraction_of_optimum=0.9,
        )
        result["fva"] = {
            reaction: {
                "minimum": float(fva.at[(reaction, "minimum")]),
                "maximum": float(fva.at[(reaction, "maximum")]),
            }
            for reaction in spec["feature_reactions"]
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-id", default="iML1515", choices=model_ids())
    parser.add_argument("--preset", default="baseline")
    parser.add_argument("--glucose", type=float, default=DEFAULT_GLUCOSE)
    parser.add_argument("--oxygen", type=float, default=DEFAULT_OXYGEN)
    parser.add_argument("--fva", action="store_true")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    config = load_model_config(args.model_id)
    payload = run_fba_context(
        model_id=args.model_id,
        preset=args.preset,
        glucose=args.glucose,
        oxygen=args.oxygen,
        with_fva=args.fva,
        model_sha256=config["sha256"],
    )
    text = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if payload.get("status") == "optimal" else 1


if __name__ == "__main__":
    raise SystemExit(main())
