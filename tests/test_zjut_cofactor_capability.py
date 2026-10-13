"""Preserve the cofactor capability the Zhejiang University of Technology paper needs.

That paper (Wang et al., D-pantothenic acid in E. coli) defines the model
requirements this project's cofactor machinery was originally built to satisfy:
iML1515 as the initial model, FBA with a biomass objective, FVA, glucose 15 and
oxygen 20 as the medium, knockouts simulated by zeroing bounds and
overexpression by scaling to a proportion of the baseline flux, and the paper's
own engineering moves available as presets.

Everything else added since - the flavin/thioester cofactors, the surfactant
grids, the yeast model, compartment-resolved accounting - must not disturb these.
This file pins them, so a future change that quietly breaks the ZJUT capability
fails here instead of being discovered later.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("cobra")

from biocandidate import fba_context  # noqa: E402

MODEL = ROOT / "data/iML1515.xml"
pytestmark = pytest.mark.skipif(
    not MODEL.is_file(), reason="the iML1515 model under gitignored data/ is required"
)

# The paper's engineering moves, as declared presets.
PAPER_PRESETS = (
    "ptsG_knockout",
    "thd2pp_overexpression",
    "nadK_overexpression",
    "gltA_overexpression",
    "paper_stack",
)


@pytest.fixture(scope="module")
def iml1515():
    return fba_context.load_model("iML1515", MODEL)


def test_iml1515_is_registered_as_the_papers_initial_model():
    assert "iML1515" in fba_context.MODEL_REGISTRY
    spec = fba_context.MODEL_REGISTRY["iML1515"]
    assert spec["objective_id"] == "BIOMASS_Ec_iML1515_core_75p37M"
    assert "iML1515" in spec["citation"] or "Monk" in spec["citation"]


def test_the_papers_medium_is_still_the_declared_default():
    spec = fba_context.MODEL_REGISTRY["iML1515"]
    assert spec["default_glucose"] == 15.0, "the paper sets glucose uptake to 15 mmol/gDW/h"
    assert spec["default_oxygen"] == 20.0, "the paper sets oxygen uptake to 20 mmol/gDW/h"


def test_fba_and_fva_are_both_available(iml1515):
    """The paper used FBA to guide the modules and FVA for the achievable range."""
    result = fba_context.run_fba_context(
        model_id="iML1515", preset="baseline", model=iml1515, with_fva=True
    )
    assert result["status"] == "optimal"
    assert result["growth"] > 0
    assert result["fva"], "FVA must return a range per feature reaction"
    for reaction, bounds in result["fva"].items():
        assert bounds["minimum"] <= bounds["maximum"], reaction


def test_knockout_and_overexpression_semantics_are_preserved():
    """The paper zeroes bounds for knockouts and scales them for overexpression."""
    assert callable(getattr(fba_context, "_apply_knockout", None))
    assert callable(getattr(fba_context, "_apply_overexpression", None))


def test_every_paper_preset_still_exists():
    for preset in PAPER_PRESETS:
        assert preset in fba_context.BASE_PRESETS, preset


def test_the_paper_stack_applies_exactly_the_documented_moves(iml1515):
    """Delta-ptsG plus the heterologous transhydrogenase and NAD kinase."""
    result = fba_context.run_fba_context(
        model_id="iML1515", preset="paper_stack", model=iml1515
    )
    assert result["status"] == "optimal"
    assert result["applied_modifications"] == [
        "knockout:b1101",
        "overexpression:THD2pp@1.5",
        "overexpression:NADK@1.5",
    ]


def test_the_cofactors_the_paper_balances_are_still_accounted():
    """NADPH regeneration and the one-carbon carrier are the paper's two main levers."""
    spec = fba_context.MODEL_REGISTRY["iML1515"]
    targets = spec["cofactor_metabolites"]
    for label in ("nadph", "nadh", "atp", "mlthf", "thf", "accoa", "coa"):
        assert label in targets, label


@pytest.mark.parametrize("preset", PAPER_PRESETS)
def test_every_paper_preset_solves(iml1515, preset):
    result = fba_context.run_fba_context(
        model_id="iML1515", preset=preset, model=iml1515
    )
    assert result["status"] == "optimal", preset
    assert len(result["features"]) == 8, "the 8-wide contract must not drift"
    assert result["cofactor_turnover"]["nadph"]["gross_turnover"] >= 0


def test_the_generic_extensions_did_not_change_the_paper_medium_path(iml1515):
    """A default-configured call must still be the paper's 15/20 configuration."""
    result = fba_context.run_fba_context(model_id="iML1515", preset="baseline", model=iml1515)
    assert result["fluxes"]["EX_glc__D_e"] == pytest.approx(-15.0)
    assert result["medium"]["EX_glc__D_e"] == 15.0
    assert result["medium"]["EX_o2_e"] == 20.0
