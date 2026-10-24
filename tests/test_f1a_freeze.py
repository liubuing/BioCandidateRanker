"""Tests for the F1a prospective prediction freeze.

The manifest exists so the predicted ordering cannot be adjusted after the assay.
These tests check that it resolves to real stored predictions, records their
identity, and states honestly whether the ordering is distinguishable at all.

中文：F1a 前瞻预测冻结测试：清单解析到真实存档预测、记录其身份、诚实陈述排序是否可区分，防止测定后调整预测顺序。
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/freeze_f1a_predictions.py"
PREDICTIONS = ROOT / "artifacts/10h2da-yeast-prediction/task2_10hda_terminal_oxidation-predictions.json"

pytestmark = pytest.mark.skipif(
    not PREDICTIONS.is_file(),
    reason="stored 10HDA predictions under gitignored artifacts/ are required",
)


def _load():
    spec = importlib.util.spec_from_file_location("f1a_freeze", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def freeze():
    return _load().build_freeze()


def test_freezes_exactly_the_four_f1a_enzymes(freeze):
    genes = {entry["gene"] for entry in freeze["predictions"]}
    assert genes == {"HFD1", "ALD2", "ALD3", "ALD6"}
    assert all(entry["uniprot"] for entry in freeze["predictions"])


def test_predictions_resolve_to_stored_values(freeze):
    """Values must come from the stored file, not be recomputed or rounded."""
    stored = json.loads(PREDICTIONS.read_text(encoding="utf-8"))["predictions"]
    by_id = {entry["candidate_id"]: entry["tasks"]["log10_kcat"]["mean"] for entry in stored}
    for entry in freeze["predictions"]:
        assert entry["predicted_log10_kcat"] == pytest.approx(by_id[entry["candidate_id"]])
        assert entry["predicted_kcat_s-1"] == pytest.approx(10 ** entry["predicted_log10_kcat"])


def test_ranks_are_a_strict_ordering(freeze):
    ranks = [entry["predicted_rank"] for entry in freeze["predictions"]]
    assert sorted(ranks) == [1, 2, 3, 4]
    values = [entry["predicted_log10_kcat"] for entry in freeze["predictions"]]
    assert values == sorted(values, reverse=True), "predictions must be listed best-first"


def test_records_the_identity_of_every_input(freeze):
    identity = freeze["prediction_identity"]
    assert len(identity["predictions_sha256"]) == 64
    assert identity["predictions_sha256"] == _load().sha256(PREDICTIONS)
    if identity["input_spec_sha256"]:
        assert len(identity["input_spec_sha256"]) == 64
    for name, digest in identity["sequence_files"]:
        assert name and len(digest) == 64


def test_decision_rule_is_declared_before_the_assay(freeze):
    rule = freeze["decision_rule"]
    assert rule["pre_registered"] is True
    assert "spearman" in rule["primary"].lower()
    assert rule["note"]
    assert freeze["status"] == "frozen_before_assay"


def test_discriminating_power_is_recorded_and_reported_honestly(freeze):
    """The check that bounds what this test can show, stated before measuring.

    With the current model the predicted spread across these four enzymes is
    smaller than the model's own standard deviation on each of them, so the
    ordering is not statistically distinguishable. Recording that in advance is
    the point: it stops a later disagreement from being oversold as refutation.
    """
    check = freeze["discriminating_power_check"]
    assert check["predicted_spread_log10"] > 0
    assert check["max_standard_deviation_log10"] > 0
    assert check["spread_exceeds_uncertainty"] is (
        check["predicted_spread_log10"] > check["max_standard_deviation_log10"]
    )
    assert "not statistically distinguishable" in check["note"]
