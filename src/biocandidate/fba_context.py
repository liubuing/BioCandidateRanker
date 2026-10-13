"""Genome-scale FBA context features for the multimodal ranker.

This module imports metabolic-model technology stacks and exposes them as
governed model inputs:

* the cofactor-engineering paper by Wang et al. (D-pantothenic acid production):
  the iML1515 model, FBA with a biomass objective, knockouts and overexpression
  simulated through reaction bounds, medium glucose 15 / oxygen 20 mmol/gDW/h,
  and optional flux variability analysis;
* the user's arginine / carbon-fixation project (05_FBA 融合FBA计算.py): the
  iJO1366 model plus eight synthetic reactions (CCR_EMA, NOG_F6P, PTXD, ...) and
  the NOG / phosphite / CO2 scenario caps;
* the user's yeast pathway-prediction project (Yeast-MetaTwin): the yeast
  genome-scale model of the 10H2DA reference case, registered read-only with a
  yeast-native 8-wide feature set (see MODEL_REGISTRY).

All produce the same contract: an 8-wide fba_context vector whose reaction
identifiers are declared per model, with hash-identified model files and a
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
from typing import Any, Callable, Sequence

from .data import FBA_FEATURE_SCHEMA_VERSION

MODEL_CONFIGS = {
    "iML1515": Path("configs/fba_context_model.json"),
    "iJO1366": Path("configs/fba_context_model_ijo1366.json"),
    "yeast_metatwin": Path("configs/fba_context_model_yeast_metatwin.json"),
}

GLUCOSE_EXCHANGE = "EX_glc__D_e"
OXYGEN_EXCHANGE = "EX_o2_e"

# The paper monitors EMP, PPP, ED, and TCA fluxes plus redox and energy state;
# these six reactions cover the same gates. Together with the per-model biomass
# objective and glucose exchange they form the 8-wide vector that matches
# ModelConfig.fba_context_dim of every shipped checkpoint. Both E. coli models
# carry every reaction below except the biomass objective, which is per-model.
SHARED_FLUX_REACTIONS = ("PFK", "G6PDH2r", "EDD", "CS", "THD2pp", "ATPM")

DEFAULT_GLUCOSE = 15.0
DEFAULT_OXYGEN = 20.0

# Yeast-MetaTwin registration: feature reactions are yeast-native r_ identifiers.
# Yeast has no ED pathway and no membrane transhydrogenase, so the two slots that
# E. coli fills with EDD/THD2pp carry the oxygen exchange and the peroxisomal
# POX1 acyl-CoA oxidase step (r_0120) that gates the 10H2DA precursor chain.
# Glucose 10 mmol/gDW/h is the standard aerobic glucose batch rate; oxygen 1000
# is the model file's own exchange cap.
YEAST_METATWIN_MODEL_ID = "yeast_metatwin"
YEAST_FEATURE_REACTIONS = (
    "r_2111",  # growth (biomass objective)
    "r_1714",  # D-glucose exchange
    "r_1992",  # oxygen exchange
    "r_0886",  # phosphofructokinase (glycolysis)
    "r_0466",  # glucose 6-phosphate dehydrogenase (PPP)
    "r_0300",  # citrate synthase (TCA)
    "r_4046",  # non-growth-associated maintenance
    "r_0120",  # acyl-CoA oxidase (decanoyl-CoA), POX1 precursor gate
)

# Cofactor accounting targets, declared per model. The paper's endpoint is the
# intracellular NADPH / NADH / ATP pool; FBA cannot report a pool level, so the
# governed quantity here is gross turnover (production == consumption at steady
# state) plus the reactions that dominate it.
#
# CORE_COFACTORS is the subset the distilled surrogates train on and the set the
# workbench panel leads with. Everything declared in a model's
# cofactor_metabolites map is reported by run_fba_context, so the registry can
# grow without touching the surrogate contract.
CORE_COFACTORS = ("nadph", "nadh", "atp")

# Oxidative burden handles. Flavin-dependent OXIDASES (the acyl-CoA oxidase
# family) carry their FAD as an internal cycling catalyst, so the net reaction
# contains O2 and H2O2 and no flavin - which is why FAD turnover cannot stand in
# for oxidase activity. What such a step actually consumes is oxygen and what it
# creates is peroxide, and both are metabolites the models track.
OXIDANT_METABOLITES = {
    "iML1515": {"h2o2": "h2o2_c", "o2": "o2_c"},
    "iJO1366": {"h2o2": "h2o2_c", "o2": "o2_c"},
    "yeast_metatwin": {"h2o2": "s_0840", "o2": "s_1279"},
}

E_COLI_COFACTORS = {
    # core redox and energy
    "nadph": "nadph_c",
    "nadh": "nadh_c",
    "atp": "atp_c",
    # oxidised partners: turnover is reportable, the redox ratio is not (FBA has
    # no concentrations)
    "nad": "nad_c",
    "nadp": "nadp_c",
    # flavins: 29 FAD and 25 FADH2 reactions in iML1515
    "fad": "fad_c",
    "fadh2": "fadh2_c",
    # thioesters: CoA appears in 80 reactions, acetyl-CoA in 41, malonyl-CoA is
    # the fatty-acid synthase substrate
    "coa": "coa_c",
    "accoa": "accoa_c",
    "malcoa": "malcoa_c",
    # one-carbon carriers (the pantothenate paper's third lever)
    "thf": "thf_c",
    "mlthf": "mlthf_c",
    "gtp": "gtp_c",
}

# The yeast model uses s_ ids and keeps separate copies per compartment; the
# cytosolic one is declared unless the compartment is functionally required.
# Note POX1 (r_0120) is written without an explicit flavin in this model, so its
# activity must be read from the reaction flux itself, not from FAD turnover.
YEAST_COFACTORS = {
    "nadph": "s_1212",
    "nadh": "s_1203",
    "atp": "s_0434",
    "nad": "s_1198",
    "nadp": "s_1207",
    "fad": "s_0687",
    "fadh2": "s_0689",
    "coa": "s_0529",
    "accoa": "s_0373",
    "malcoa": "s_1101",
    "thf": "s_1487",
    "mlthf": "s_0306",
    "gtp": "s_0785",
}

MODEL_REGISTRY: dict[str, dict[str, Any]] = {
    "iML1515": {
        "objective_id": "BIOMASS_Ec_iML1515_core_75p37M",
        "citation": "Monk et al. 2017, Nature Microbiology 2:17094",
        "exchanges": {"glucose": GLUCOSE_EXCHANGE, "oxygen": OXYGEN_EXCHANGE},
        "feature_reactions": ("BIOMASS_Ec_iML1515_core_75p37M", GLUCOSE_EXCHANGE)
        + SHARED_FLUX_REACTIONS,
        "presets": None,  # None means every BASE_PRESETS entry applies
        "cofactor_metabolites": E_COLI_COFACTORS,
        "oxidant_metabolites": OXIDANT_METABOLITES["iML1515"],
        "default_glucose": DEFAULT_GLUCOSE,
        "default_oxygen": DEFAULT_OXYGEN,
        "provenance_key": "iML1515_presets",
    },
    "iJO1366": {
        "objective_id": "BIOMASS_Ec_iJO1366_WT_53p95M",
        "citation": "Orth et al. 2011, Mol Syst Biol 7:530",
        "exchanges": {"glucose": GLUCOSE_EXCHANGE, "oxygen": OXYGEN_EXCHANGE},
        "feature_reactions": ("BIOMASS_Ec_iJO1366_WT_53p95M", GLUCOSE_EXCHANGE)
        + SHARED_FLUX_REACTIONS,
        "presets": None,
        "cofactor_metabolites": E_COLI_COFACTORS,
        "oxidant_metabolites": OXIDANT_METABOLITES["iJO1366"],
        "default_glucose": DEFAULT_GLUCOSE,
        "default_oxygen": DEFAULT_OXYGEN,
        "provenance_key": "iJO1366_fusion",
    },
    YEAST_METATWIN_MODEL_ID: {
        "objective_id": "r_2111",
        "citation": (
            "Yeast-MetaTwin yeast genome-scale model (user metabolic-modeling "
            "project asset, Yeast-GEM style r_ identifiers); opened read-only, "
            "re-hashed copy under data/"
        ),
        "exchanges": {"glucose": "r_1714", "oxygen": "r_1992"},
        "feature_reactions": YEAST_FEATURE_REACTIONS,
        "presets": ("baseline",),
        "cofactor_metabolites": YEAST_COFACTORS,
        "oxidant_metabolites": OXIDANT_METABOLITES["yeast_metatwin"],
        # 1.0 is the model file's own glucose uptake cap and the medium the
        # external 10H2DA anchor (0.274011) was produced with. An earlier value
        # of 10.0 was a generic batch rate that this model cannot represent and
        # that fell outside the surrogate's trained grid, so a default-configured
        # call and the documented reference case disagreed.
        "default_glucose": 1.0,
        "default_oxygen": 1000.0,
        "provenance_key": "yeast_metatwin",
    },
}
for _spec in MODEL_REGISTRY.values():
    _spec["feature_reactions"] = tuple(_spec["feature_reactions"])
    _spec["feature_ids"] = tuple(f"flux:{r}" for r in _spec["feature_reactions"])

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
    "yeast_metatwin": (
        "user project Yeast-MetaTwin yeast model (10H2DA reference case); "
        "registered 2026-10-04 from the read-only source, SHA256 pinned in "
        "configs/fba_context_model_yeast_metatwin.json"
    ),
}


def model_ids() -> list[str]:
    return sorted(MODEL_REGISTRY)


def presets_for(model_id: str) -> list[str]:
    _require_model(model_id)
    spec = MODEL_REGISTRY[model_id]
    names = sorted(BASE_PRESETS) if spec["presets"] is None else list(spec["presets"])
    if model_id == "iJO1366":
        names += sorted(FUSION_PRESETS)
    return names


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
    elif resolved.suffix in {".yml", ".yaml"}:
        model = cobra.io.load_yaml_model(str(resolved))
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


def _set_medium(model, spec: dict[str, Any], fusion_spec: dict[str, Any] | None,
                glucose: float, oxygen: float) -> None:
    exchanges = spec["exchanges"]
    for exchange, bound in ((exchanges["glucose"], glucose), (exchanges["oxygen"], oxygen)):
        if bound < 0:
            raise ValueError(f"uptake must be non-negative, got {bound}")
        model.reactions.get_by_id(exchange).lower_bound = -float(bound)
    if fusion_spec:
        model.reactions.EX_co2_e.lower_bound = -float(FUSION_CO2_UPTAKE_CAP)
        model.reactions.NOG_F6P.upper_bound = float(fusion_spec["nog"])
        model.reactions.EX_phite_e.lower_bound = -float(fusion_spec["phosphite"])
        if fusion_spec["phosphite"] == 0:
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


def cofactor_turnover(
    model, solution, cofactor_metabolites: dict[str, str], top_n: int = 5
) -> dict[str, Any]:
    """Gross turnover per cofactor, plus the reactions that dominate it.

    Net production of an internal cofactor is exactly zero at FBA steady state,
    so the governed quantity is gross turnover: the summed positive
    contributions, which equal the summed consumption. That is the closest
    model-side analogue of the paper's "NADPH regeneration capacity", and the
    dominant producers identify which pathway supplies it (for iML1515 at
    baseline, the pentose phosphate shunt for NADPH and oxidative
    phosphorylation for ATP - the same levers the paper engineers).

    Cofactor pools themselves are not reported: FBA gives fluxes, not
    concentrations, and no metabolite-level measurement is bundled here.
    """
    report: dict[str, Any] = {}
    for label, metabolite_id in cofactor_metabolites.items():
        if metabolite_id not in model.metabolites:
            report[label] = {
                "status": "metabolite_absent",
                "metabolite_id": metabolite_id,
            }
            continue
        metabolite = model.metabolites.get_by_id(metabolite_id)
        gross = 0.0
        net = 0.0
        producers: list[tuple[str, float]] = []
        for reaction in model.reactions:
            coefficient = reaction.metabolites.get(metabolite)
            if coefficient is None:
                continue
            contribution = coefficient * float(solution.fluxes[reaction.id])
            net += contribution
            if coefficient > 0 and contribution > 0:
                gross += contribution
                producers.append((reaction.id, contribution))
        producers.sort(key=lambda item: -item[1])
        report[label] = {
            "status": "reported",
            "metabolite_id": metabolite_id,
            "unit": "mmol gDW^-1 h^-1",
            "gross_turnover": gross,
            # Kept as a self-check: a non-zero value means the solution does not
            # satisfy steady state for this metabolite.
            "net_production": net,
            "top_producers": [
                {"reaction_id": reaction_id, "flux": flux}
                for reaction_id, flux in producers[:top_n]
            ],
        }
    return report




def oxidative_burden(
    model,
    solution,
    oxidant_metabolites: dict[str, str],
    top_n: int = 5,
) -> dict[str, Any]:
    """Oxygen consumed and peroxide produced, plus the reactions responsible.

    Reported per direction rather than net: a cell consumes oxygen for both
    respiration and oxidases, and produces peroxide that its peroxiredoxins and
    catalase must clear, so the gross figures are the ones with engineering
    meaning.
    """
    report: dict[str, Any] = {}
    for label, metabolite_id in oxidant_metabolites.items():
        if metabolite_id not in model.metabolites:
            report[label] = {"status": "metabolite_absent", "metabolite_id": metabolite_id}
            continue
        metabolite = model.metabolites.get_by_id(metabolite_id)
        produced = consumed = 0.0
        producers: list[tuple[str, float]] = []
        consumers: list[tuple[str, float]] = []
        for reaction in model.reactions:
            coefficient = reaction.metabolites.get(metabolite)
            if coefficient is None:
                continue
            contribution = coefficient * float(solution.fluxes[reaction.id])
            if contribution > 0:
                produced += contribution
                producers.append((reaction.id, contribution))
            elif contribution < 0:
                consumed += -contribution
                consumers.append((reaction.id, -contribution))
        producers.sort(key=lambda item: -item[1])
        consumers.sort(key=lambda item: -item[1])
        report[label] = {
            "status": "reported",
            "metabolite_id": metabolite_id,
            "unit": "mmol gDW^-1 h^-1",
            "produced": produced,
            "consumed": consumed,
            "top_producers": [
                {"reaction_id": rid, "flux": flux} for rid, flux in producers[:top_n]
            ],
            "top_consumers": [
                {"reaction_id": rid, "flux": flux} for rid, flux in consumers[:top_n]
            ],
        }
    return report


def _least_squares(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Slope, intercept and R^2 of a simple linear fit.

    R^2 is reported because a "marginal cost" is only a meaningful single number
    when the response is close to linear in the product level; a poor fit means
    the cost is regime-dependent and the slope should not be quoted alone.
    """
    count = len(xs)
    if count < 2:
        return float("nan"), float("nan"), float("nan")
    mean_x = sum(xs) / count
    mean_y = sum(ys) / count
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0:
        return float("nan"), float("nan"), float("nan")
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    if ss_tot == 0:
        return slope, intercept, 1.0
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    return slope, intercept, 1.0 - ss_res / ss_tot


def marginal_cofactor_cost(
    model,
    target_reaction: str,
    *,
    cofactor_metabolites: dict[str, str],
    oxidant_metabolites: dict[str, str] | None = None,
    objective_reaction: str | None = None,
    level_fractions: tuple[float, ...] = (0.25, 0.5, 0.75, 1.0),
    context: Callable[[Any], None] | None = None,
) -> dict[str, Any]:
    """Cost of producing one more unit of the target, in cofactor turnover.

    Reporting total turnover at a single optimum does not answer what to
    engineer: the LP sits at the same vertex, so the numbers do not move when the
    product level changes. The informative quantity is the marginal one, obtained
    by fixing the target flux at graded levels and optimising something else -
    biomass by default - then regressing each tracked quantity on the level.

    The context callable applies the scenario (medium, route, gate) to the model
    copy, so this stays model- and scenario-agnostic. A biomass floor is NOT
    imposed here: fixing the product flux at a level already determines how much
    growth is possible, and keeping a floor as well would make the curve a
    constrained slice rather than the trade-off itself.
    """
    base = model.copy()
    if context is not None:
        context(base)

    probe = base.copy()
    probe.objective = target_reaction
    probe_solution = probe.optimize()
    if probe_solution.status != "optimal":
        return {"status": probe_solution.status, "error": "target is not maximisable"}
    capacity = float(probe_solution.objective_value)
    if capacity <= 0:
        return {"status": "no_capacity", "capacity": capacity}

    # Without a named objective the target itself would be maximised, which pins
    # the LP at one vertex and yields no trade-off curve, so the objective column
    # is left empty rather than reported as a meaningless constant.
    objective_named = objective_reaction is not None

    tracked = dict(cofactor_metabolites)
    if oxidant_metabolites:
        tracked.update({f"oxidant:{name}": mid for name, mid in oxidant_metabolites.items()})

    points: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for fraction in level_fractions:
        level = capacity * fraction
        with base:
            target = base.reactions.get_by_id(target_reaction)
            target.lower_bound = level
            target.upper_bound = level
            if objective_reaction is not None:
                base.objective = objective_reaction
            solution = base.optimize()
            if solution.status != "optimal":
                skipped.append({"fraction": fraction, "level": level,
                                "status": solution.status})
                continue
            entry: dict[str, Any] = {
                "fraction": fraction,
                "level": level,
                "objective_value": (
                    float(solution.fluxes[objective_reaction]) if objective_named else None
                ),
                "cofactors": {},
            }
            turnover = cofactor_turnover(base, solution, tracked)
            for label, info in turnover.items():
                entry["cofactors"][label] = info.get("gross_turnover")
            points.append(entry)

    if len(points) < 2:
        return {
            "status": "insufficient_points",
            "capacity": capacity,
            "points": points,
            "skipped": skipped,
        }

    levels = [p["level"] for p in points]
    marginal: dict[str, Any] = {}
    for label in points[0]["cofactors"]:
        series = [p["cofactors"][label] for p in points]
        if any(value is None for value in series):
            marginal[label] = {"status": "incomplete"}
            continue
        slope, _, r2 = _least_squares(levels, series)
        marginal[label] = {
            "slope_per_unit_product": slope,
            "r_squared": r2,
            "linear": r2 == r2 and r2 > 0.9,
        }

    objective_series = [p["objective_value"] for p in points]
    objective_summary = None
    if all(value is not None for value in objective_series):
        slope, _, r2 = _least_squares(levels, objective_series)
        objective_summary = {
            "slope_per_unit_product": slope,
            "r_squared": r2,
            "linear": r2 == r2 and r2 > 0.9,
        }

    ranked = sorted(
        (
            (label, info["slope_per_unit_product"])
            for label, info in marginal.items()
            if info.get("slope_per_unit_product") is not None
        ),
        key=lambda item: -item[1],
    )
    return {
        "status": "reported",
        "target_reaction": target_reaction,
        "capacity": capacity,
        "level_fractions": list(level_fractions),
        "points": points,
        "skipped": skipped,
        "marginal": marginal,
        "objective": objective_summary,
        "largest_positive_cost": ranked[0][0] if ranked and ranked[0][1] > 0 else None,
        "ranking_by_cost": [
            {"label": label, "slope_per_unit_product": slope} for label, slope in ranked
        ],
        "interpretation": (
            "A positive slope is cofactor turnover the network must add per unit of "
            "product; a negative slope means the product route is cheaper in that "
            "cofactor than the activity it displaces. Slopes with a low R-squared are "
            "regime-dependent and should not be quoted as a single number."
        ),
    }


def compartmental_copies(model, metabolite_id: str) -> dict[str, str]:
    """Every compartmental copy of one chemical species, keyed by compartment.

    Copies are matched by metabolite name, which both models use consistently
    ("NADPH", "NADH", "ATP" repeated per compartment). The declared id anchors
    identity; only its own name is expanded, so an unrelated metabolite that
    merely shares a substring cannot leak in.
    """
    reference = model.metabolites.get_by_id(metabolite_id)
    name = (reference.name or "").strip()
    if not name:
        return {reference.compartment: metabolite_id}
    return {
        metabolite.compartment: metabolite.id
        for metabolite in model.metabolites
        if (metabolite.name or "").strip() == name
    }


def cofactor_turnover_by_compartment(
    model,
    solution,
    cofactor_metabolites: dict[str, str],
    top_n: int = 3,
) -> dict[str, Any]:
    """Cofactor turnover split by compartment.

    A single-pool report is misleading whenever a pathway sits outside the
    cytosol: the royal-jelly-acid route is peroxisomal and consumes the
    peroxisomal NADPH copy, so a cytosolic-only figure describes a pool the route
    never touches. Each compartment is labelled active when it carries flux.
    """
    report: dict[str, Any] = {}
    for label, anchor_id in cofactor_metabolites.items():
        copies = compartmental_copies(model, anchor_id)
        per_compartment: dict[str, Any] = {}
        for compartment, metabolite_id in sorted(copies.items()):
            single = cofactor_turnover(
                model, solution, {label: metabolite_id}, top_n=top_n
            )[label]
            if single.get("status") != "reported":
                per_compartment[compartment] = {"status": single.get("status")}
                continue
            per_compartment[compartment] = {
                "status": "reported",
                "metabolite_id": metabolite_id,
                "gross_turnover": single["gross_turnover"],
                "net_production": single["net_production"],
                "active": single["gross_turnover"] > 1e-9,
                "top_producers": single["top_producers"],
            }
        report[label] = per_compartment
    return report


def enzyme_capacity_scan(
    model,
    reaction_id: str,
    target_reaction: str,
    capacities: Sequence[float],
    *,
    context: Callable[[Any], None] | None = None,
    reference_flux: float | None = None,
) -> dict[str, Any]:
    """Maximum target flux as a function of one reaction's flux ceiling.

    Genome-scale models assume every enzyme is unlimited, which hides the step
    that a real cell cannot keep up with. Capping one reaction and re-maximising
    the target exposes it: the target stays flat while the cap is slack and falls
    once the cap binds, so the knee is the capacity that reaction must reach.

    `reference_flux` (typically the reaction's flux at the unconstrained optimum)
    lets the knee be expressed as a multiple of the demand, which is the number
    an experimentalist can compare against a measured or predicted kcat.
    """
    base = model.copy()
    if context is not None:
        context(base)

    probe = base.copy()
    probe.objective = target_reaction
    unconstrained = float(probe.optimize().objective_value)
    if reference_flux is None:
        reference_flux = float(probe.optimize().fluxes[reaction_id])

    rows: list[dict[str, Any]] = []
    for capacity in capacities:
        with base:
            reaction = base.reactions.get_by_id(reaction_id)
            reaction.upper_bound = float(capacity)
            base.objective = target_reaction
            solution = base.optimize()
            rows.append(
                {
                    "capacity": float(capacity),
                    "status": solution.status,
                    "target": (
                        float(solution.objective_value) if solution.status == "optimal" else None
                    ),
                }
            )

    feasible = [row for row in rows if row["target"] is not None]
    knee = None
    # the knee is the smallest capacity already reaching the unconstrained
    # optimum: every capacity at or above it leaves the target untouched
    for row in sorted(feasible, key=lambda item: item["capacity"]):
        if unconstrained - row["target"] <= 1e-9:
            knee = row["capacity"]
            break

    # Refine by bisection between the last limited and the first unlimited grid
    # point, because a grid-resolution answer can overstate the requirement by a
    # wide margin (a coarse run reported 0.4 where the true knee is 0.303).
    refined = None
    if knee is not None:
        below = [row["capacity"] for row in feasible
                 if row["capacity"] < knee and unconstrained - row["target"] > 1e-9]
        if below:
            low, high = max(below), knee
            for _ in range(24):
                mid = (low + high) / 2.0
                with base:
                    base.reactions.get_by_id(reaction_id).upper_bound = mid
                    base.objective = target_reaction
                    trial = base.optimize()
                    reaches = (
                        trial.status == "optimal"
                        and unconstrained - float(trial.objective_value) <= 1e-9
                    )
                if reaches:
                    high = mid
                else:
                    low = mid
            refined = high

    return {
        "status": "reported",
        "reaction_id": reaction_id,
        "target_reaction": target_reaction,
        "unconstrained_target": unconstrained,
        "reference_flux": reference_flux,
        "knee_capacity": knee,
        "knee_rounded_from_grid": refined if refined is not None else None,
        "knee_as_multiple_of_reference": (
            (refined if refined is not None else knee) / reference_flux
            if reference_flux
            else None
        ),
        "rows": rows,
        "interpretation": (
            "A target that falls as the cap tightens means the capped reaction is a "
            "bottleneck at that capacity. The knee is the smallest capacity that "
            "still supports the unconstrained optimum, i.e. the minimum this enzyme "
            "must deliver."
        ),
    }


def run_fba_context(
    *,
    model_id: str = "iML1515",
    preset: str = "baseline",
    glucose: float | None = None,
    oxygen: float | None = None,
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
    before solving. Glucose/oxygen default to the registered per-model values
    when omitted.
    """
    if model_id not in MODEL_REGISTRY:
        raise ValueError(f"unknown model {model_id!r}; known: {model_ids()}")
    spec = MODEL_REGISTRY[model_id]
    glucose = spec.get("default_glucose", DEFAULT_GLUCOSE) if glucose is None else float(glucose)
    oxygen = spec.get("default_oxygen", DEFAULT_OXYGEN) if oxygen is None else float(oxygen)
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
        "medium": {
            spec["exchanges"]["glucose"]: glucose,
            spec["exchanges"]["oxygen"]: oxygen,
        },
        "solver_id": solver_id,
        "feature_ids": list(spec["feature_ids"]),
        "claim_boundary": CLAIM_BOUNDARY,
        "provenance": (
            PROVENANCE["iJO1366_fusion"] if fusion_spec else PROVENANCE[spec["provenance_key"]]
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
    _set_medium(scoped, spec, fusion_spec, glucose, oxygen)
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
    # Cofactor accounting rides alongside the 8-wide feature vector rather than
    # extending it: every shipped checkpoint declares fba_context_dim = 8, so
    # widening the vector would invalidate all of them.
    if spec.get("cofactor_metabolites"):
        result["cofactor_turnover"] = cofactor_turnover(
            scoped, solution, spec["cofactor_metabolites"]
        )
        result["oxidative_burden"] = oxidative_burden(
            scoped, solution, spec.get("oxidant_metabolites", {})
        )
        # the single-pool figure above is anchored on the declared (cytosolic)
        # copy; this split shows which compartment's pool a route actually uses
        result["cofactor_by_compartment"] = cofactor_turnover_by_compartment(
            scoped, solution, spec["cofactor_metabolites"]
        )
        result["cofactor_scope"] = (
            "gross turnover (production equals consumption at steady state); "
            "pool levels and concentrations are not modelled"
        )

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
    parser.add_argument("--glucose", type=float, default=None,
                        help="glucose uptake cap; default is per-model")
    parser.add_argument("--oxygen", type=float, default=None,
                        help="oxygen uptake cap; default is per-model")
    parser.add_argument("--fva", action="store_true")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    config = load_model_config(args.model_id)
    spec = MODEL_REGISTRY[args.model_id]
    payload = run_fba_context(
        model_id=args.model_id,
        preset=args.preset,
        glucose=args.glucose if args.glucose is not None else spec.get("default_glucose", DEFAULT_GLUCOSE),
        oxygen=args.oxygen if args.oxygen is not None else spec.get("default_oxygen", DEFAULT_OXYGEN),
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
