"""Prosthetic-group ledger: the enzyme-bound cofactor bill that FBA cannot see.

A genome-scale model tracks cofactors as currency: metabolites consumed and
regenerated at steady state, so their requirement scales with flux. The
flavins in this project's yeast model measure near zero on that ledger, which
is correct as far as it goes - and is exactly what the experimental side
contradicts, because the route's P450 needs its reductase, and the reductase
carries one FAD and one FMN per enzyme molecule regardless of how fast it
turns over. That demand is a prosthetic-group inventory: it scales with enzyme
expression level, it is spent on assembly and repaid on protein turnover, and
it never appears in any reaction's stoichiometry.

This ledger makes that demand computable from the project's own numbers:

  enzyme load  = route flux / (kcat x 3600)          mmol/gDW
  inventory    = load x prosthetics per enzyme       mmol/gDW
  replacement  = inventory / protein half-life       mmol/gDW/h
  burden ratio = replacement / endogenous synthesis  dimensionless

The burden ratio is the number the experimental observation translates into.
Below one, the cell's own flavin-synthesis capacity covers assembly and the
experiment should see nothing; above one, expressing the route enzymes asks
for more prosthetic synthesis than the housekeeping demand, and precursor
feeding or RIB-gene overexpression becomes a visible lever. Neither FBA nor
this ledger models the assembly itself; what it prices is the gap between the
two ledgers, currency and prosthetic.

All kcat values are uncalibrated ranker-style guesses exposed as a scenario
grid, not measurements; the coupling ratio of reductase to P450 is idealized
at one. This is a bridge between two ledgers, not a flux prediction.
"""

from __future__ import annotations

from typing import Any

SECONDS_PER_HOUR = 3600.0

# Prosthetic content per enzyme molecule, one entry per enzyme the route
# expresses. Textbook biochemistry, not model output: P450s carry one
# protoheme IX; the yeast NADPH-cytochrome P450 reductase NCP1 carries one FAD
# and one FMN and shuttles electrons NADPH -> FAD -> FMN -> P450 heme; POX1 is
# a flavoprotein oxidase with one tightly bound FAD (which is why its GEM
# reaction correctly shows no flavin - see the royal jelly module); TES1 is a
# serine hydrolase with none. Dissociable cofactors (NADPH, CoA, ATP) are
# deliberately absent: the currency ledger already prices them per flux.
ENZYME_PROSTHETIC_CONTENT: dict[str, dict[str, float]] = {
    "CYP52M1": {"heme_b": 1.0},
    "NCP1": {"fad": 1.0, "fmn": 1.0},
    "POX1": {"fad": 1.0},
    "TES1": {},
}

# The yeast model's dedicated biosynthesis reactions per prosthetic group.
# These give the endogenous capacity the route demand is compared against.
# FMN's kinase can sit at zero in some optimal vertices because the model
# recycles FMN through an FMNH2-dependent hydroxylase cycle; the chemically
# required floor for FMN is the FAD demand, one-to-one, and the report states
# it that way instead of trusting the vertex.
YEAST_SYNTHESIS_REACTIONS: dict[str, str] = {
    "fad": "r_0440",  # FMN adenylyltransferase
    "fmn": "r_0965",  # riboflavin kinase
    "heme_b": "r_0436",  # ferrochelatase
}

# Groups are only meaningful with these names; the ledger refuses a group it
# has no synthesis reaction for rather than silently dropping it.
KNOWN_GROUPS = frozenset(YEAST_SYNTHESIS_REACTIONS)


def enzyme_load(flux_mmol_gdw_h: float, kcat_s: float) -> float:
    """Enzyme load (mmol/gDW) that sustains a flux: ub = kcat x E x 3600."""
    if kcat_s <= 0:
        raise ValueError(f"kcat must be positive, got {kcat_s}")
    if flux_mmol_gdw_h < 0:
        raise ValueError(f"flux must be non-negative, got {flux_mmol_gdw_h}")
    return flux_mmol_gdw_h / (kcat_s * SECONDS_PER_HOUR)


def ledger_for_fluxes(
    fluxes: dict[str, float],
    *,
    kcat_s: float,
    protein_half_life_h: float,
    reductase_coupling: float = 1.0,
    extra_enzymes: dict[str, float] | None = None,
) -> dict[str, Any]:
    """The prosthetic bill for one operating point.

    ``fluxes`` maps enzyme name to the flux its reaction carries (mmol/gDW/h);
    ``extra_enzymes`` adds enzymes whose load is driven by coupling rather than
    by their own reaction - the P450 reductase, whose load is the P450's load
    times ``reductase_coupling``. The bill reports per-enzyme loads, the
    inventory and replacement flux per prosthetic group, and carries the
    assumption block it was computed under so no consumer can quote the
    numbers without the kcat and half-life they rest on.
    """
    if protein_half_life_h <= 0:
        raise ValueError(f"protein half-life must be positive, got {protein_half_life_h}")

    loads: dict[str, float] = {}
    for enzyme, flux in fluxes.items():
        if enzyme not in ENZYME_PROSTHETIC_CONTENT:
            raise KeyError(
                f"no prosthetic content recorded for {enzyme!r}; "
                f"known: {sorted(ENZYME_PROSTHETIC_CONTENT)}"
            )
        loads[enzyme] = enzyme_load(flux, kcat_s)
    if extra_enzymes:
        for enzyme, reference in extra_enzymes.items():
            if enzyme not in ENZYME_PROSTHETIC_CONTENT:
                raise KeyError(
                    f"no prosthetic content recorded for {enzyme!r}; "
                    f"known: {sorted(ENZYME_PROSTHETIC_CONTENT)}"
                )
            loads[enzyme] = loads[reference] * reductase_coupling

    inventory: dict[str, float] = {}
    per_enzyme_inventory: dict[str, dict[str, float]] = {}
    for enzyme, load in loads.items():
        content = ENZYME_PROSTHETIC_CONTENT[enzyme]
        per_enzyme_inventory[enzyme] = {
            group: load * count for group, count in content.items()
        }
        for group, amount in per_enzyme_inventory[enzyme].items():
            if group not in KNOWN_GROUPS:
                raise KeyError(
                    f"prosthetic group {group!r} has no synthesis reaction to "
                    "compare against; extend YEAST_SYNTHESIS_REACTIONS"
                )
            inventory[group] = inventory.get(group, 0.0) + amount

    replacement = {
        group: amount / protein_half_life_h for group, amount in inventory.items()
    }
    return {
        "kcat_s": kcat_s,
        "protein_half_life_h": protein_half_life_h,
        "reductase_coupling": reductase_coupling,
        "enzyme_loads": loads,
        "enzyme_prosthetic_inventory": per_enzyme_inventory,
        "group_inventory_mmol_gdw": inventory,
        "group_replacement_flux_mmol_gdw_h": replacement,
        "enzymes_without_prosthetics": [
            enzyme for enzyme, content in ENZYME_PROSTHETIC_CONTENT.items()
            if enzyme in loads and not content
        ],
    }


def reference_synthesis_capacity(model, solution, growth: float) -> dict[str, Any]:
    """What the cell's own biomass asks of each prosthetic group, per growth.

    Read from the dedicated biosynthesis reactions at the reference optimum.
    The per-growth intensity is the stable quantity: it is the prosthetic
    content of the biomass composition, independent of the medium, so the
    burden ratio compares route demand against composition rather than against
    one LP vertex. The FMN row is reported as a one-to-one floor of FAD rather
    than as the kinase vertex, which the FMN/FMNH2 cycle can park at zero.
    """
    if growth <= 0:
        raise ValueError(f"growth must be positive, got {growth}")
    rows: dict[str, Any] = {}
    for group, reaction_id in YEAST_SYNTHESIS_REACTIONS.items():
        flux = float(solution.fluxes[reaction_id])
        rows[group] = {
            "synthesis_reaction": reaction_id,
            "flux_mmol_gdw_h": flux,
            "per_growth_intensity_mmol_gdw": flux / growth,
        }
    fad_intensity = rows["fad"]["per_growth_intensity_mmol_gdw"]
    rows["fmn"]["floor_note"] = (
        "the kinase vertex can be zero because the model recycles FMN through "
        "FMNH2-dependent hydroxylases; every FAD made needs one FMN first, so "
        "the FMN floor is the FAD intensity"
    )
    rows["fmn"]["floor_per_growth_intensity_mmol_gdw"] = fad_intensity
    return rows


def burden_ratios(bill: dict[str, Any], capacity: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Route replacement demand against endogenous intensity, per group.

    The replacement flux is a steady-state demand, so it is compared against
    the endogenous *flux* row directly; the inventory is additionally reported
    against the biomass inventory (the endogenous intensity scaled to the
    route's own growth), because a tank that takes hours to fill can bind even
    when the refill rate looks comfortable.

    A group the reference biomass does not use at all - heme b in this model,
    whose biomass carries heme a and heme o instead - has a zero baseline by
    construction. Its ratios are reported as null rather than as a division by
    zero, with a note saying the route demand is new load from a standing
    start: there is no housekeeping budget to compare against, which is a
    stronger statement than any finite ratio.
    """
    growth = bill.get("growth_at_operating_point")
    out: dict[str, dict[str, Any]] = {}
    for group, replacement in bill["group_replacement_flux_mmol_gdw_h"].items():
        endogenous_flux = capacity[group]["flux_mmol_gdw_h"]
        # FMN's kinase vertex can park at zero while the chemistry still needs
        # one FMN per FAD, so the FMN comparison uses the stated floor.
        if group == "fmn":
            floor = capacity["fmn"]["floor_per_growth_intensity_mmol_gdw"]
            if growth:
                endogenous_flux = max(endogenous_flux, floor * growth)
        inventory = bill["group_inventory_mmol_gdw"][group]
        biomass_inventory = (
            capacity[group]["per_growth_intensity_mmol_gdw"] * growth if growth else 0.0
        )
        no_baseline = endogenous_flux <= 0 and biomass_inventory <= 0
        out[group] = {
            "replacement_vs_endogenous_flux": (
                replacement / endogenous_flux
                if endogenous_flux > 0
                else None
            ),
            "inventory_vs_biomass_inventory": (
                inventory / biomass_inventory if biomass_inventory > 0 else None
            ),
            "no_endogenous_baseline": no_baseline,
            "note": (
                "the reference biomass carries none of this group, so every "
                "route demand is new load and the ratios are undefined rather "
                "than favourable"
            )
            if no_baseline
            else None,
        }
    return out


def route_ledger(
    target_config_path,
    *,
    route: str = "A",
    pox1_closed: bool = False,
    kcat_grid_s: tuple[float, ...] = (0.1, 1.0, 10.0),
    half_life_grid_h: tuple[float, ...] = (10.0, 72.0),
    reductase_coupling: float = 1.0,
) -> dict[str, Any]:
    """The prosthetic bill for the royal-jelly-acid route, as a scenario grid.

    Fluxes come from the module's own route analysis at the documented anchor
    configuration, so the bill's numerator is the same capacity every other
    analysis quotes. The P450 reductase rides at ``reductase_coupling`` times
    the hydroxylase load. The grid varies the two assumptions the bill is most
    sensitive to - kcat and protein half-life - because the honest headline is
    the boundary where the burden ratio crosses one, not a single number.
    """
    from .fba_context import load_model
    from .royal_jelly_acid import analyze

    report = analyze(target_config_path, pox1_closed=pox1_closed)
    route_name = f"route_{route}_" + ("free_acid" if route == "A" else "coa_thioester")
    entry = report["routes"][route_name]
    if entry.get("status") != "optimal":
        return {"status": entry.get("status"), "error": "route capacity unavailable"}
    capacity = float(entry["target_capacity"])

    model = load_model("yeast_metatwin")
    reference = model.optimize()
    growth_ref = float(reference.objective_value)
    capacity_table = reference_synthesis_capacity(model, reference, growth_ref)

    # Every route step carries the full route flux; POX1 rides at its measured
    # flux and drops to zero when the entrance is closed. NCP1 (the P450
    # reductase) is driven by coupling to the hydroxylase, not by its own LP.
    pox1_flux = 0.0 if pox1_closed else max(float(entry["pox1_flux"]), 0.0)
    fluxes = {"CYP52M1": capacity, "TES1": capacity, "POX1": pox1_flux}

    scenarios: dict[str, Any] = {}
    for kcat in kcat_grid_s:
        for half_life in half_life_grid_h:
            bill = ledger_for_fluxes(
                fluxes,
                kcat_s=kcat,
                protein_half_life_h=half_life,
                reductase_coupling=reductase_coupling,
                extra_enzymes={"NCP1": "CYP52M1"},
            )
            # the capacity the route analysis quotes is computed at the biomass
            # floor, so that floor - not the unconstrained optimum - is the
            # growth the inventory is compared against
            bill["growth_at_operating_point"] = float(report["biomass_floor"])
            bill["burden_ratios"] = burden_ratios(bill, capacity_table)
            bill["load_bearing"] = {
                group: (
                    row["replacement_vs_endogenous_flux"] is not None
                    and row["replacement_vs_endogenous_flux"] > 1.0
                )
                for group, row in bill["burden_ratios"].items()
            }
            scenarios[f"kcat_{kcat:g}_half_life_{half_life:g}h"] = bill

    return {
        "status": "reported",
        "route": route,
        "pox1_closed": pox1_closed,
        "route_capacity": capacity,
        "route_flux_assignment": {
            "CYP52M1": capacity,
            "TES1": capacity,
            "POX1": pox1_flux,
            "NCP1 (coupled to CYP52M1)": capacity * reductase_coupling,
        },
        "growth_at_operating_point": float(report["biomass_floor"]),
        "reference_growth": growth_ref,
        "endogenous_synthesis_capacity": capacity_table,
        "scenarios": scenarios,
        "reading": (
            "burden_ratios[group].replacement_vs_endogenous_flux is the number "
            "the experimental observation translates into: the route's "
            "prosthetic replacement demand as a multiple of what the cell "
            "already spends on that group for its own biomass. Above one, "
            "assembly demand exceeds the housekeeping budget and precursor "
            "feeding or RIB-gene overexpression has something to relieve; the "
            "kcat and half-life at which the ratio crosses one are the honest "
            "headline, not any single scenario."
        ),
        "claim_boundary": (
            "A ledger bridging two accounts, not a flux prediction. kcat values "
            "are uncalibrated scenario inputs; the reductase-to-P450 coupling is "
            "idealized at one; prosthetic assembly, import into the peroxisome "
            "and apoenzyme stability are not modelled anywhere. The currency "
            "ledger (FBA) remains the only flux-consistent account."
        ),
    }
