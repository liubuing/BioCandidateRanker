"""Tests for the interactive workbench server.

The workbench runs real model inference, so its tests focus on the guards that
keep a browser request equivalent to the CLI it mirrors: candidate validation,
path confinement, the governance allowlist, and — when a checkpoint is present
locally — numerical agreement with predict_command.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import serve_workbench as workbench  # noqa: E402

REFERENCE_CHECKPOINT = ROOT / "artifacts/esm2-t6-opt-seed7/best.pt"
pytestmark = pytest.mark.skipif(
    not REFERENCE_CHECKPOINT.is_file(),
    reason="local checkpoint under gitignored artifacts/ is required",
)


# --- candidate validation --------------------------------------------------------


def test_parse_candidates_normalizes_sequences():
    rows = workbench.parse_candidates(
        {"candidates": [{"sequence": "mst agk", "substrate_smiles": " CCO "}]}
    )
    assert rows[0]["sequence"] == "MSTAGK"
    assert rows[0]["substrate_smiles"] == "CCO"
    assert rows[0]["candidate_id"] == "0"


def test_parse_candidates_accepts_a_bare_list_and_preserves_metadata():
    rows = workbench.parse_candidates(
        [{"candidate_id": "x", "sequence": "ACD", "substrate_smiles": "CCO",
          "organism": "E. coli", "ec": "1.1.1.1", "reaction": "a->b"}]
    )
    assert rows[0]["candidate_id"] == "x"
    assert rows[0]["organism"] == "E. coli"
    assert rows[0]["ec"] == "1.1.1.1"
    assert rows[0]["reaction"] == "a->b"


@pytest.mark.parametrize(
    "payload,reason",
    [
        ({"candidates": []}, "empty"),
        ([], "empty"),
        ({"candidates": "nope"}, "not a list"),
        ({"candidates": [{"substrate_smiles": "CCO"}]}, "missing sequence"),
        ({"candidates": [{"sequence": "ACD"}]}, "missing smiles"),
        ({"candidates": [{"sequence": "hello1", "substrate_smiles": "CCO"}]}, "bad residue"),
        ({"candidates": ["ACD"]}, "non-object row"),
    ],
)
def test_parse_candidates_rejects_invalid_input(payload, reason):
    with pytest.raises(ValueError):
        workbench.parse_candidates(payload)


# --- path confinement ------------------------------------------------------------


def test_model_cache_rejects_paths_outside_the_repository(tmp_path):
    cache = workbench.ModelCache(tmp_path, device="cpu")
    with pytest.raises(ValueError, match="inside the repository"):
        cache.resolve("../escaped.pt")


def test_model_cache_rejects_missing_and_non_checkpoint_paths(tmp_path):
    cache = workbench.ModelCache(tmp_path, device="cpu")
    with pytest.raises(ValueError, match="not found"):
        cache.resolve_data("absent.json")
    (tmp_path / "notes.txt").write_text("not a checkpoint", encoding="utf-8")
    with pytest.raises(ValueError, match="checkpoint file"):
        cache.resolve("notes.txt")
    with pytest.raises(ValueError, match="inside the repository"):
        cache.resolve_data("../../etc/config.json")


# --- governance allowlist --------------------------------------------------------


def test_every_allowlisted_action_references_an_existing_script():
    for action, command in workbench.GOVERNANCE_ACTIONS.items():
        assert (ROOT / command[1]).is_file(), action


def test_run_governance_rejects_unknown_actions(tmp_path):
    with pytest.raises(ValueError, match="unknown governance action"):
        workbench.run_governance(ROOT, "rm -rf /")


def test_run_governance_executes_the_license_audit():
    result = workbench.run_governance(ROOT, "license_audit")
    assert result["ok"] is True
    assert result["result"]["corpus_count"] >= 1


# --- model listing ---------------------------------------------------------------


def test_list_models_reports_validation_metrics_without_errors():
    models = workbench.list_models(ROOT)
    assert models, "the local artifacts tree should contain checkpoints"
    by_path = {m["path"]: m for m in models}
    reference = by_path["artifacts/esm2-t6-opt-seed7/best.pt"]
    assert reference["tasks"] == ["log10_kcat"]
    assert reference["validation"]["log10_kcat"]["rmse"] > 0

    baselines = [m for m in models if "feature-mlp" in m["path"]]
    assert baselines, "feature-MLP baselines should be listed"
    assert all("error" not in m for m in baselines)
    assert all(m.get("feature_baseline") for m in baselines)


# --- prediction equivalence with the CLI -----------------------------------------


def test_run_prediction_matches_predict_command(tmp_path):
    """The server must return what the CLI would write for the same input."""

    example = json.loads(
        (ROOT / "configs/candidate_example.json").read_text(encoding="utf-8")
    )
    rows = workbench.parse_candidates(example)
    # The example carries a 3-wide fba_context for an older checkpoint; the
    # workbench does not expose FBA context, so drop it for a fair comparison.
    for row in rows:
        row.pop("fba_context", None)

    cache = workbench.ModelCache(ROOT, device="cpu")
    served = workbench.run_prediction(
        cache, "artifacts/esm2-t6-opt-seed7/best.pt", rows, None
    )

    from biocandidate.cli import predict_command

    input_path = tmp_path / "candidates.json"
    output_path = tmp_path / "cli-out.json"
    input_path.write_text(json.dumps({"candidates": rows}), encoding="utf-8")
    with workbench._safe_load_globals():
        predict_command(
            argparse.Namespace(
                checkpoint=str(REFERENCE_CHECKPOINT),
                input=str(input_path),
                output=str(output_path),
                calibration_artifact=None,
                fba_features=None,
                device="cpu",
            )
        )
    cli = json.loads(output_path.read_text(encoding="utf-8"))

    assert served["predictions"] == cli["predictions"]
    assert served["warning"] == cli["warning"]


def test_run_prediction_rejects_a_mismatched_calibration_artifact():
    rows = workbench.parse_candidates(
        [{"sequence": "ACDEFG", "substrate_smiles": "CCO"}]
    )
    cache = workbench.ModelCache(ROOT, device="cpu")
    with pytest.raises(ValueError, match="identity mismatch|schema"):
        workbench.run_prediction(
            cache,
            "artifacts/esm2-t6-opt-seed7/best.pt",
            rows,
            "artifacts/uncertainty-calibration-v2/grouped-cv.json",
        )


# --- frontend and server stay in sync -------------------------------------------


def test_frontend_actions_are_all_allowlisted():
    html = (ROOT / "scripts/workbench_static/index.html").read_text(encoding="utf-8")
    frontend_actions = set(re.findall(r'data-action="([a-z0-9_]+)"', html))
    assert frontend_actions <= set(workbench.GOVERNANCE_ACTIONS)


def test_frontend_is_self_contained():
    html = (ROOT / "scripts/workbench_static/index.html").read_text(encoding="utf-8")
    assert re.search(r'src\s*=\s*["\']https?://', html) is None
    assert "<link" not in html


def test_index_asset_exists_next_to_the_server():
    assert workbench.INDEX_PATH.is_file()


def test_frontend_defaults_to_chinese_with_english_toggle():
    """The UI is bilingual; Chinese is the default, switchable without reload."""
    html = (ROOT / "scripts/workbench_static/index.html").read_text(encoding="utf-8")
    assert "const LANG_DEFAULT = 'zh'" in html
    assert "zh:" in html and "en:" in html
    assert "lang-btn" in html
    # The Chinese strings are real UI text, not just a language-code marker.
    assert "预测 kcat" in html
    assert "Predict kcat" in html


def test_frontend_i18n_covers_the_same_action_ids_in_both_languages():
    html = (ROOT / "scripts/workbench_static/index.html").read_text(encoding="utf-8")
    for action in workbench.GOVERNANCE_ACTIONS:
        occurrences = html.count(f"{action}:")
        assert occurrences >= 2, f"{action} missing from one language block"
