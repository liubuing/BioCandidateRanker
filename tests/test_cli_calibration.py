import json
from argparse import Namespace
from types import SimpleNamespace

import pytest
import torch

from biocandidate import cli
from biocandidate.calibration import (
    ScalarGaussianCalibrator,
    SplitConformalCalibrator,
    save_calibration_artifact,
)


class FakeModel:
    config = SimpleNamespace(
        task_names=("log10_kcat", "log10_km"),
        context_buckets=8,
        fba_context_dim=0,
    )

    def eval(self):
        return self

    def __call__(self, batch):
        return {
            "mean": torch.tensor([[1.0, 2.0]]),
            "standard_deviation": torch.tensor([[2.0, 3.0]]),
        }


def test_fit_calibration_uses_only_validation_and_explicit_method(tmp_path, monkeypatch):
    checkpoint_path = tmp_path / "model.pt"
    checkpoint_path.write_bytes(b"frozen checkpoint")
    data_path = tmp_path / "data.json"
    data_path.write_text("[]", encoding="utf-8")
    split_path = tmp_path / "split.json"
    split_path.write_text(json.dumps({
        "format_version": 1,
        "rows": [
            {"source_row": 0, "split": "train", "cluster": "train"},
            {"source_row": 1, "split": "validation", "cluster": "validation"},
            {"source_row": 2, "split": "test", "cluster": "test"},
        ],
    }), encoding="utf-8")
    source_identity = cli._manifest_dict(cli.build_manifest(data_path))
    source_manifest = {"identity": source_identity}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(source_manifest), encoding="utf-8")
    run_manifest = {
        "source": source_manifest,
        "split": {
            "path": str(split_path.resolve()),
            "identity": cli._manifest_dict(cli.build_manifest(split_path)),
        },
    }
    records = [SimpleNamespace(split=name) for name in ("train", "validation", "test")]
    seen = []

    monkeypatch.setattr(cli, "verify_manifest", lambda *args: None)
    selected_source_rows = []

    def read_validation(path, *, source_rows):
        selected_source_rows.append(source_rows)
        return SimpleNamespace(records=records)

    monkeypatch.setattr(cli, "read_unikp_json", read_validation)
    monkeypatch.setattr(cli, "apply_split_manifest",
                     lambda values, path, source_identity: values)
    monkeypatch.setattr(cli, "load_checkpoint", lambda path, device: (FakeModel(), {
        "data_manifest": run_manifest,
    }))

    def collect(model, values, task, batch_size, device):
        seen.extend(record.split for record in values)
        return torch.tensor([0.0]), torch.tensor([2.0]), torch.tensor([4.0])

    monkeypatch.setattr(cli, "_collect_task_outputs", collect)
    output_path = tmp_path / "calibration.json"
    cli.fit_calibration_command(Namespace(
        data=str(data_path), manifest=str(manifest_path), split_manifest=str(split_path),
        checkpoint=str(checkpoint_path), task="log10_kcat", method="scalar",
        conflict_policy="keep", batch_size=8, device="cpu", output=str(output_path),
    ))

    artifact = json.loads(output_path.read_text(encoding="utf-8"))
    assert seen == ["validation"]
    assert selected_source_rows == [{1}]
    assert artifact["fit_partition"] == "validation"
    assert artifact["validation_count"] == 1
    assert artifact["calibrator"] == {"kind": "scalar", "scale": 2.0}
    assert artifact["identities"] == {
        "checkpoint": cli._binary_identity(checkpoint_path),
        "task": "log10_kcat",
    }


def test_predict_applies_matching_task_calibration_and_emits_raw_stddev(tmp_path, monkeypatch):
    checkpoint_path = tmp_path / "model.pt"
    checkpoint_path.write_bytes(b"checkpoint")
    artifact_path = tmp_path / "calibration.json"
    save_calibration_artifact(
        artifact_path,
        ScalarGaussianCalibrator(1.5),
        identities=cli._calibration_identities(checkpoint_path, "log10_kcat"),
        validation_count=3,
        conformal_calibrators=(
            SplitConformalCalibrator(1.0, 0.1),
            SplitConformalCalibrator(2.0, 0.1, normalized=True),
            SplitConformalCalibrator(1.5, 0.05),
            SplitConformalCalibrator(2.5, 0.05, normalized=True),
        ),
    )
    input_path = tmp_path / "candidates.json"
    input_path.write_text(json.dumps([{
        "candidate_id": "candidate-1", "sequence": "ACD", "substrate_smiles": "CCO",
    }]), encoding="utf-8")
    output_path = tmp_path / "predictions.json"

    monkeypatch.setattr(cli, "load_checkpoint", lambda path, device: (FakeModel(), {
        "data_manifest": {},
        "metrics": {"validation_tasks": {
            "log10_kcat": {"count": 1}, "log10_km": {"count": 1},
        }},
    }))
    monkeypatch.setattr(cli, "EnzymeSubstrateCollator", lambda **kwargs: lambda records: {
        "placeholder": torch.tensor([0.0]),
    })
    cli.predict_command(Namespace(
        checkpoint=str(checkpoint_path), input=str(input_path), output=str(output_path),
        calibration_artifact=str(artifact_path), device="cpu",
    ))

    result = json.loads(output_path.read_text(encoding="utf-8"))
    calibrated = result["predictions"][0]["tasks"]["log10_kcat"]
    untouched = result["predictions"][0]["tasks"]["log10_km"]
    assert calibrated["standard_deviation"] == 2.0
    assert calibrated["raw_standard_deviation"] == 2.0
    assert calibrated["calibrated_standard_deviation"] == 3.0
    assert calibrated["conformal_intervals"] == {
        "90_normalized": {"lower": -5.0, "upper": 7.0},
        "90_unnormalized": {"lower": 0.0, "upper": 2.0},
        "95_normalized": {"lower": -6.5, "upper": 8.5},
        "95_unnormalized": {"lower": -0.5, "upper": 2.5},
    }
    assert untouched["standard_deviation"] == 3.0
    assert untouched["calibrated_standard_deviation"] is None
    assert untouched["conformal_intervals"] is None
    assert result["calibration"]["task"] == "log10_kcat"


def test_predict_rejects_calibration_for_another_checkpoint(tmp_path, monkeypatch):
    checkpoint_path = tmp_path / "model.pt"
    checkpoint_path.write_bytes(b"current")
    other_path = tmp_path / "other.pt"
    other_path.write_bytes(b"other")
    artifact_path = tmp_path / "calibration.json"
    save_calibration_artifact(
        artifact_path,
        ScalarGaussianCalibrator(1.5),
        identities=cli._calibration_identities(other_path, "log10_kcat"),
        validation_count=3,
    )
    monkeypatch.setattr(cli, "load_checkpoint", lambda path, device: (FakeModel(), {}))
    with pytest.raises(ValueError, match="identity mismatch"):
        cli.predict_command(Namespace(
            checkpoint=str(checkpoint_path), input="unused", output="unused",
            calibration_artifact=str(artifact_path), device="cpu",
        ))


def test_calibration_cli_arguments_are_explicit():
    parser = cli.build_parser()
    args = parser.parse_args([
        "fit-calibration", "--checkpoint", "model.pt", "--data", "data.json",
        "--manifest", "manifest.json", "--split-manifest", "split.json",
        "--task", "log10_kcat", "--method", "identity", "--output", "calibration.json",
    ])
    assert args.func is cli.fit_calibration_command
    assert args.method == "identity"
    predict = parser.parse_args([
        "predict", "--checkpoint", "model.pt", "--input", "input.json",
        "--output", "output.json", "--calibration-artifact", "calibration.json",
    ])
    assert predict.calibration_artifact == "calibration.json"


def test_grouped_cv_calibration_never_requests_test_rows(tmp_path, monkeypatch):
    checkpoint_path = tmp_path / "model.pt"
    checkpoint_path.write_bytes(b"frozen checkpoint")
    data_path = tmp_path / "data.json"
    data_path.write_text("[]", encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    source_manifest = {"identity": cli._manifest_dict(cli.build_manifest(data_path))}
    manifest_path.write_text(json.dumps(source_manifest), encoding="utf-8")
    split_path = tmp_path / "split.json"
    split_payload = {
        "format_version": 1,
        "rows": [
            {"source_row": 0, "split": "train", "cluster": "train-family"},
            {"source_row": 1, "split": "validation", "cluster": "v1"},
            {"source_row": 2, "split": "validation", "cluster": "v2"},
            {"source_row": 3, "split": "test", "cluster": "test-family"},
        ],
    }
    split_path.write_text(json.dumps(split_payload), encoding="utf-8")
    run_manifest = {
        "source": source_manifest,
        "split": {
            "path": str(split_path.resolve()),
            "identity": cli._manifest_dict(cli.build_manifest(split_path)),
        },
    }
    validation_records = [
        SimpleNamespace(source_row=1, source_rows=(), split="validation"),
        SimpleNamespace(source_row=2, source_rows=(), split="validation"),
    ]
    requested = []

    def read_rows(path, *, source_rows):
        requested.extend(sorted(source_rows))
        assert 3 not in source_rows
        return SimpleNamespace(records=validation_records)

    monkeypatch.setattr(cli, "verify_manifest", lambda *args: None)
    monkeypatch.setattr(cli, "read_unikp_json", read_rows)
    monkeypatch.setattr(cli, "apply_split_manifest",
                     lambda records, path, source_identity: records)
    monkeypatch.setattr(cli, "load_checkpoint", lambda path, device: (FakeModel(), {
        "data_manifest": run_manifest,
    }))
    monkeypatch.setattr(cli, "_collect_task_outputs", lambda *args: (
        torch.tensor([0.0, 0.0]), torch.ones(2), torch.tensor([-1.0, 1.0]),
    ))
    selected = SimpleNamespace(
        calibrator=ScalarGaussianCalibrator(1.0),
        to_dict=lambda: {"selected_kind": "scalar"},
    )
    seen_groups = []

    def select(mean, deviation, labels, groups, **kwargs):
        seen_groups.extend(groups)
        return selected

    monkeypatch.setattr(cli, "select_grouped_gaussian_calibrator", select)
    monkeypatch.setattr(
        cli.SplitConformalCalibrator, "fit",
        lambda mean, labels, *, alpha, standard_deviation=None, normalized=False:
        cli.SplitConformalCalibrator(1.0, alpha, normalized),
    )
    output_path = tmp_path / "calibration-v2.json"
    cli.fit_calibration_command(Namespace(
        data=str(data_path), manifest=str(manifest_path), split_manifest=str(split_path),
        checkpoint=str(checkpoint_path), task="log10_kcat", method="grouped-cv",
        conflict_policy="keep", batch_size=8, device="cpu", output=str(output_path),
        cv_folds=2, cv_seed=0, nll_tolerance=1e-12, crps_tolerance=1e-12,
    ))

    artifact = json.loads(output_path.read_text(encoding="utf-8"))
    assert requested == [1, 2]
    assert seen_groups == ["v1", "v2"]
    assert artifact["format_version"] == 2
    assert set(artifact["conformal_intervals"]) == {
        "90_normalized", "90_unnormalized", "95_normalized", "95_unnormalized",
    }
