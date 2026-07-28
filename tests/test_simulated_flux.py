import json

import pytest

from biocandidate.simulated_flux import (
    CLAIM_BOUNDARY,
    FEATURE_IDS,
    TARGET_REACTION_ID,
    read_simulated_flux_jsonl,
    run_pipeline,
)


def test_pipeline_writes_blocker_instead_of_partial_artifacts(tmp_path):
    output = tmp_path / "output"
    with pytest.raises(RuntimeError, match="model asset is not a file"):
        run_pipeline(tmp_path / "missing.yml", output)
    blocker = json.loads((output / "blocker.json").read_text(encoding="utf-8"))
    assert blocker["status"] == "blocked"
    assert blocker["task"] == "simulated_flux"
    assert "experimental flux" in blocker["claim_boundary"]
    assert sorted(path.name for path in output.iterdir()) == ["blocker.json"]


def test_generated_fixture_enforces_simulated_contract_and_no_target_leakage(tmp_path):
    row = {
        "format_version": 1,
        "source_row": 1,
        "condition_id": "condition-1",
        "split": "train",
        "features": [-0.5, -2.0],
        "fba": {
            "schema_version": 1,
            "feature_ids": list(FEATURE_IDS),
            "model_id": "sha256:" + "a" * 64,
            "solver_id": "sha256:" + "b" * 64,
            "objective_id": "sha256:" + "c" * 64,
        },
        "log10_flux": -1.0,
        "label": {
            "type": "simulated",
            "target_reaction_id": TARGET_REACTION_ID,
            "canonical_unit": "log10(mmol gDW^-1 h^-1)",
            "provenance": {"claim_boundary": CLAIM_BOUNDARY},
        },
    }
    path = tmp_path / "data.jsonl"
    path.write_text(json.dumps(row) + "\n", encoding="utf-8")
    record = read_simulated_flux_jsonl(path)[0]
    assert record.flux_label_metadata.label_type.value == "simulated"
    assert TARGET_REACTION_ID not in record.fba_feature_metadata.feature_ids
    row["fba"]["feature_ids"][0] = TARGET_REACTION_ID
    path.write_text(json.dumps(row) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="target reaction"):
        read_simulated_flux_jsonl(path)
