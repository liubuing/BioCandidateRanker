import json

import pytest
from argparse import Namespace
from pathlib import Path

from biocandidate import cli
from biocandidate.config import ModelConfig


class FakeModel:
    def __init__(self, task_names=("log10_kcat", "log10_km")):
        self.config = ModelConfig(task_names=task_names)

    def eval(self):
        return self

    def __call__(self, batch):
        import torch
        n = len(batch.get("context_ids", [])) if "context_ids" in batch else 1
        return {
            "mean": torch.zeros((n, len(self.config.task_names))),
            "log_variance": torch.zeros((n, len(self.config.task_names))),
            "standard_deviation": torch.ones((n, len(self.config.task_names))),
        }


FAKE_CHECKPOINT = {
    "data_manifest": {},
    "metrics": {"validation_tasks": {
        "log10_kcat": {"count": 1},
        "log10_km": {"count": 1},
    }},
}


def _write_candidates(path: Path, rows) -> None:
    path.write_text(json.dumps({"candidates": rows}), encoding="utf-8")


def _fba_metadata_payload(width: int) -> dict:
    return {
        "schema_version": 1,
        "feature_ids": [f"fba::{index}" for index in range(width)],
        "model_id": "yeast-meta-twin-v1",
        "solver_id": "glpk",
        "objective_id": "biomass",
        "condition_id": "glucose-limited-aerobic",
    }


def _run_predict(tmp_path, monkeypatch, *, fba_features_path=None, rows=None, task_names=("log10_kcat", "log10_km")):
    checkpoint_path = tmp_path / "model.pt"
    input_path = tmp_path / "input.json"
    output_path = tmp_path / "predictions.json"
    _write_candidates(input_path, rows or [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
    }])
    monkeypatch.setattr(cli, "load_checkpoint",
                        lambda path, device: (FakeModel(task_names=task_names), FAKE_CHECKPOINT))
    monkeypatch.setattr(cli, "EnzymeSubstrateCollator", lambda **kwargs: lambda records: {
        "mean": __import__("torch").zeros((1, 2)),
        "log_variance": __import__("torch").zeros((1, 2)),
        "standard_deviation": __import__("torch").ones((1, 2)),
    })
    cli.predict_command(Namespace(
        checkpoint=str(checkpoint_path), input=str(input_path), output=str(output_path),
        calibration_artifact=None, device="cpu",
        fba_features=str(fba_features_path) if fba_features_path is not None else None,
    ))
    return json.loads(output_path.read_text(encoding="utf-8"))


def test_predict_with_empty_context_without_fba_metadata(tmp_path, monkeypatch):
    result = _run_predict(tmp_path, monkeypatch)
    assert result["predictions"][0]["tasks"]["log10_kcat"]["status"] == "trained"


def test_predict_accepts_exact_width_fba_context(tmp_path, monkeypatch):
    fba_path = tmp_path / "fba.json"
    fba_path.write_text(json.dumps(_fba_metadata_payload(8)), encoding="utf-8")
    rows = [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
        "fba_context": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    }]
    result = _run_predict(tmp_path, monkeypatch, fba_features_path=fba_path, rows=rows)
    assert result["predictions"][0]["candidate_id"] == "candidate-1"


def test_predict_rejects_fba_context_without_metadata(tmp_path, monkeypatch):
    rows = [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
        "fba_context": [0.1],
    }]
    with pytest.raises(ValueError, match="--fba-features"):
        _run_predict(tmp_path, monkeypatch, rows=rows)


def test_predict_rejects_metadata_width_mismatch(tmp_path, monkeypatch):
    fba_path = tmp_path / "fba.json"
    fba_path.write_text(json.dumps(_fba_metadata_payload(2)), encoding="utf-8")
    rows = [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
        "fba_context": [0.1, 0.2],
    }]
    with pytest.raises(ValueError, match="feature metadata width"):
        _run_predict(tmp_path, monkeypatch, fba_features_path=fba_path, rows=rows)


def test_predict_rejects_invalid_fba_metadata(tmp_path, monkeypatch):
    fba_path = tmp_path / "fba.json"
    fba_path.write_text(json.dumps({"feature_ids": []}), encoding="utf-8")
    rows = [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
        "fba_context": [0.1],
    }]
    with pytest.raises(ValueError, match="invalid FBA feature metadata"):
        _run_predict(tmp_path, monkeypatch, fba_features_path=fba_path, rows=rows)


def test_fba_feature_width_validated_against_checkpoint_before_any_context(tmp_path, monkeypatch):
    fba_path = tmp_path / "fba.json"
    fba_path.write_text(json.dumps(_fba_metadata_payload(8)), encoding="utf-8")
    input_path = tmp_path / "input.json"
    output_path = tmp_path / "predictions.json"
    _write_candidates(input_path, [{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
        "fba_context": [0.1] * 8,
    }])
    fake_model = FakeModel(task_names=("log10_kcat", "log10_km"))
    fake_model.config = ModelConfig(task_names=("log10_kcat", "log10_km"), fba_context_dim=2)
    monkeypatch.setattr(cli, "load_checkpoint",
                        lambda path, device: (fake_model, FAKE_CHECKPOINT))
    monkeypatch.setattr(cli, "EnzymeSubstrateCollator", lambda **kwargs: lambda records: {
        "mean": __import__("torch").zeros((1, 2)),
        "log_variance": __import__("torch").zeros((1, 2)),
        "standard_deviation": __import__("torch").ones((1, 2)),
    })
    with pytest.raises(ValueError, match="feature metadata width"):
        cli.predict_command(Namespace(
            checkpoint=str(tmp_path / "model.pt"), input=str(input_path),
            output=str(output_path), calibration_artifact=None, device="cpu",
            fba_features=str(fba_path),
        ))