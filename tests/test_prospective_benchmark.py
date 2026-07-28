import csv
import hashlib
import json
from pathlib import Path

import pytest

from biocandidate.prospective_benchmark import (
    audit_readiness,
    final_evaluation,
    freeze_prediction_deposit,
)


def identity(path: Path) -> dict:
    data = path.read_bytes()
    return {"sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def governed(tmp_path):
    protocol = tmp_path / "protocol.json"
    protocol.write_text(
        json.dumps(
            {
                "readiness_gates": {
                    "minimum_campaigns": 2,
                    "minimum_families": 2,
                    "minimum_candidates_per_campaign": 3,
                },
                "campaign_metrics": ["spearman", "pairwise_accuracy"],
            }
        ),
        encoding="utf-8",
    )
    candidates = tmp_path / "candidates.csv"
    rows = [
        {
            "candidate_id": f"{campaign}-{index}",
            "campaign_id": campaign,
            "family_id": f"family-{campaign}",
            "disposition": disposition,
        }
        for campaign in ("a", "b")
        for index, disposition in enumerate(("active", "inactive", "censored"))
    ]
    write_csv(candidates, ["candidate_id", "campaign_id", "family_id", "disposition"], rows)
    registry = tmp_path / "campaigns.json"
    registry.write_text(
        json.dumps(
            {
                "complete_candidate_roster": True,
                "labels_available_at_freeze": False,
                "selection_basis": "label_independent",
                "campaigns": {"a": 3, "b": 3},
            }
        ),
        encoding="utf-8",
    )
    family = tmp_path / "families.json"
    family.write_text(
        json.dumps(
            {
                "benchmark_family_ids": ["family-a", "family-b"],
                "all_development_and_closed_tests_covered": True,
                "references": [
                    {
                        "kind": "development",
                        "name": "development-v1",
                        "sha256": "a" * 64,
                        "overlapping_family_ids": [],
                    },
                    {
                        "kind": "closed_test",
                        "name": "all-closed-tests-v1",
                        "sha256": "b" * 64,
                        "overlapping_family_ids": [],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    checkpoint = tmp_path / "model.ckpt"
    checkpoint.write_bytes(b"frozen model")
    models = tmp_path / "models.json"
    models.write_text(
        json.dumps(
            {
                "frozen_before_labels": True,
                "labels_available_at_freeze": False,
                "models": [
                    {
                        "model_id": "primary",
                        "checkpoint": checkpoint.name,
                        "identity": identity(checkpoint),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return protocol, candidates, registry, family, models, rows


def freeze_all(tmp_path, governed):
    protocol, candidates, registry, family, models, rows = governed
    readiness = audit_readiness(protocol, candidates, registry, family, models)
    readiness_path = tmp_path / "readiness.json"
    readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
    predictions = tmp_path / "predictions.csv"
    write_csv(
        predictions,
        ["model_id", "candidate_id", "prediction"],
        [
            {"model_id": "primary", "candidate_id": row["candidate_id"], "prediction": index}
            for index, row in enumerate(rows)
        ],
    )
    deposit = freeze_prediction_deposit(
        readiness_path, protocol, candidates, registry, family, models, predictions
    )
    deposit_path = tmp_path / "deposit.json"
    deposit_path.write_text(json.dumps(deposit), encoding="utf-8")
    return readiness_path, predictions, deposit_path


def test_passing_readiness_freezes_identities_and_counts_all_dispositions(governed):
    protocol, candidates, registry, family, models, _ = governed
    report = audit_readiness(protocol, candidates, registry, family, models)
    assert report["ready"] is True
    assert report["counts"]["dispositions"] == {"active": 2, "censored": 2, "inactive": 2}
    assert report["claim_boundary"].startswith("Identity/readiness")


def test_readiness_rejects_labels_in_roster_and_family_overlap(governed):
    protocol, candidates, registry, family, models, rows = governed
    for row in rows:
        row["score"] = "1"
    write_csv(
        candidates,
        ["candidate_id", "campaign_id", "family_id", "disposition", "score"],
        rows,
    )
    family_payload = json.loads(family.read_text(encoding="utf-8"))
    family_payload["references"][1]["overlapping_family_ids"] = ["family-a"]
    family.write_text(json.dumps(family_payload), encoding="utf-8")
    report = audit_readiness(protocol, candidates, registry, family, models)
    assert report["ready"] is False
    assert any("label-like" in issue for issue in report["issues"])
    assert any("zero family overlap" in issue for issue in report["issues"])


def test_prediction_deposit_requires_complete_frozen_cross_product(tmp_path, governed):
    protocol, candidates, registry, family, models, rows = governed
    readiness = audit_readiness(protocol, candidates, registry, family, models)
    readiness_path = tmp_path / "readiness.json"
    readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
    predictions = tmp_path / "predictions.csv"
    write_csv(
        predictions,
        ["model_id", "candidate_id", "prediction"],
        [
            {"model_id": "primary", "candidate_id": row["candidate_id"], "prediction": 0}
            for row in rows[:-1]
        ],
    )
    with pytest.raises(ValueError, match="incomplete"):
        freeze_prediction_deposit(
            readiness_path, protocol, candidates, registry, family, models, predictions
        )


def test_final_join_reports_campaign_metrics_and_refuses_second_receipt(tmp_path, governed):
    _, candidates, _, _, _, rows = governed
    readiness, predictions, deposit = freeze_all(tmp_path, governed)
    labels = tmp_path / "labels.csv"
    write_csv(
        labels,
        ["candidate_id", "disposition", "label"],
        [
            {
                "candidate_id": row["candidate_id"],
                "disposition": row["disposition"],
                "label": "" if row["disposition"] == "censored" else index,
            }
            for index, row in enumerate(rows)
        ],
    )
    output = tmp_path / "final-receipt.json"
    result = final_evaluation(readiness, deposit, candidates, predictions, labels, output)
    assert set(result["models"]["primary"]["campaign_metrics"]) == {"a", "b"}
    assert result["models"]["primary"]["campaign_metrics"]["a"]["count"] == 2
    assert result["candidate_dispositions"]["censored"] == 2
    with pytest.raises(FileExistsError):
        final_evaluation(
            readiness,
            deposit,
            candidates,
            predictions,
            labels,
            tmp_path / "different-output.json",
        )


def test_final_join_requires_complete_labels_including_censored(tmp_path, governed):
    _, candidates, _, _, _, rows = governed
    readiness, predictions, deposit = freeze_all(tmp_path, governed)
    labels = tmp_path / "labels.csv"
    write_csv(
        labels,
        ["candidate_id", "disposition", "label"],
        [
            {"candidate_id": row["candidate_id"], "disposition": row["disposition"], "label": 1}
            for row in rows[:-1]
        ],
    )
    with pytest.raises(ValueError, match="every frozen candidate"):
        final_evaluation(
            readiness, deposit, candidates, predictions, labels, tmp_path / "result.json"
        )


def test_final_join_rejects_roster_changed_after_deposit(tmp_path, governed):
    _, candidates, _, _, _, rows = governed
    readiness, predictions, deposit = freeze_all(tmp_path, governed)
    rows[0]["family_id"] = "changed-family"
    write_csv(
        candidates,
        ["candidate_id", "campaign_id", "family_id", "disposition"],
        rows,
    )
    labels = tmp_path / "labels.csv"
    write_csv(
        labels,
        ["candidate_id", "disposition", "label"],
        [
            {
                "candidate_id": row["candidate_id"],
                "disposition": row["disposition"],
                "label": 1,
            }
            for row in rows
        ],
    )
    with pytest.raises(ValueError, match="candidate roster changed"):
        final_evaluation(
            readiness, deposit, candidates, predictions, labels, tmp_path / "result.json"
        )
