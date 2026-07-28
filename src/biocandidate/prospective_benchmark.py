"""Fail-closed governance for a prospective candidate-ranking benchmark.

This module freezes identities and audits protocol state. It does not train a model,
invent candidates or labels, or permit scoring before all readiness gates pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .evaluation import campaign_ranking_metrics


DISPOSITIONS = {"active", "inactive", "censored"}
LABEL_LIKE_FIELDS = {
    "activity",
    "fitness",
    "label",
    "measurement",
    "outcome",
    "response",
    "score",
    "value",
}


def _identity(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"sha256": digest.hexdigest(), "size_bytes": size}


def _load_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        return list(reader), reader.fieldnames


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value.lower())
    )


def _add_identity_issue(issues: list[str], identity: Any, path: Path, description: str) -> None:
    if not isinstance(identity, dict) or not _is_sha256(identity.get("sha256")):
        issues.append(f"{description} lacks a valid frozen SHA256 identity")
        return
    if identity != _identity(path):
        issues.append(f"{description} identity does not match {path}")


def audit_readiness(
    protocol_path: Path,
    candidates_path: Path,
    campaign_registry_path: Path,
    family_audit_path: Path,
    model_freeze_path: Path,
) -> dict[str, Any]:
    """Audit and freeze a pre-label candidate roster; never compute model metrics."""
    protocol = _load_object(protocol_path)
    registry = _load_object(campaign_registry_path)
    family_audit = _load_object(family_audit_path)
    model_freeze = _load_object(model_freeze_path)
    candidates, fields = _read_csv(candidates_path)
    issues: list[str] = []

    required_fields = {"candidate_id", "campaign_id", "family_id", "disposition"}
    missing_fields = sorted(required_fields - set(fields))
    if missing_fields:
        issues.append(f"candidate roster lacks required fields: {', '.join(missing_fields)}")
    label_fields = sorted(
        field
        for field in fields
        if field.lower() in LABEL_LIKE_FIELDS or field.lower().startswith("label_")
    )
    if label_fields:
        issues.append(
            f"candidate roster contains prohibited label-like fields: {', '.join(label_fields)}"
        )

    candidate_ids = [row.get("candidate_id", "").strip() for row in candidates]
    if any(not value for value in candidate_ids):
        issues.append("every candidate must have a nonempty candidate_id")
    duplicates = sorted(
        item for item, count in Counter(candidate_ids).items() if item and count > 1
    )
    if duplicates:
        issues.append(f"duplicate candidate_ids: {', '.join(duplicates)}")
    bad_dispositions = sorted(
        {
            row.get("disposition", "").strip()
            for row in candidates
            if row.get("disposition", "").strip() not in DISPOSITIONS
        }
    )
    if bad_dispositions:
        issues.append(f"invalid dispositions: {', '.join(bad_dispositions)}")

    campaigns: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        campaigns[row.get("campaign_id", "").strip()].append(row)
    campaigns.pop("", None)
    families = {row.get("family_id", "").strip() for row in candidates}
    families.discard("")
    if any(
        not row.get("campaign_id", "").strip() or not row.get("family_id", "").strip()
        for row in candidates
    ):
        issues.append("every candidate must have nonempty campaign_id and family_id")

    gates = protocol.get("readiness_gates", {})
    for key in ("minimum_campaigns", "minimum_families", "minimum_candidates_per_campaign"):
        if (
            not isinstance(gates.get(key), int)
            or isinstance(gates.get(key), bool)
            or gates[key] < 1
        ):
            issues.append(f"protocol readiness_gates.{key} must be a positive integer")
    if (
        isinstance(gates.get("minimum_campaigns"), int)
        and len(campaigns) < gates["minimum_campaigns"]
    ):
        issues.append(f"campaigns {len(campaigns)} < required {gates['minimum_campaigns']}")
    if isinstance(gates.get("minimum_families"), int) and len(families) < gates["minimum_families"]:
        issues.append(f"families {len(families)} < required {gates['minimum_families']}")
    minimum_candidates = gates.get("minimum_candidates_per_campaign")
    if isinstance(minimum_candidates, int) and not isinstance(minimum_candidates, bool):
        for campaign_id, rows in sorted(campaigns.items()):
            if len(rows) < minimum_candidates:
                issues.append(
                    f"campaign {campaign_id} candidates {len(rows)} < required {minimum_candidates}"
                )

    declared_campaigns = registry.get("campaigns")
    if not registry.get("complete_candidate_roster"):
        issues.append("campaign registry does not declare a complete candidate roster")
    if registry.get("labels_available_at_freeze") is not False:
        issues.append("campaign registry must declare labels_available_at_freeze false")
    if registry.get("selection_basis") != "label_independent":
        issues.append("campaign selection_basis must be label_independent")
    if not isinstance(declared_campaigns, dict):
        issues.append("campaign registry campaigns must be an object of declared counts")
    else:
        actual_counts = {name: len(rows) for name, rows in sorted(campaigns.items())}
        if declared_campaigns != actual_counts:
            issues.append("campaign registry counts do not exactly match the candidate roster")

    audited_families = family_audit.get("benchmark_family_ids")
    if not isinstance(audited_families, list) or set(audited_families) != families:
        issues.append("family audit does not exactly cover roster family_ids")
    references = family_audit.get("references")
    required_reference_kinds = {"development", "closed_test"}
    found_kinds: set[str] = set()
    if not isinstance(references, list) or not references:
        issues.append("family audit must contain development and closed-test references")
    else:
        for index, reference in enumerate(references):
            if not isinstance(reference, dict):
                issues.append(f"family reference {index} must be an object")
                continue
            kind = reference.get("kind")
            if isinstance(kind, str):
                found_kinds.add(kind)
            if not reference.get("name") or not _is_sha256(reference.get("sha256")):
                issues.append(f"family reference {index} lacks name or SHA256")
            overlap = reference.get("overlapping_family_ids")
            if overlap != []:
                issues.append(f"family reference {index} does not prove zero family overlap")
        missing_kinds = sorted(required_reference_kinds - found_kinds)
        if missing_kinds:
            issues.append(f"family audit lacks reference kinds: {', '.join(missing_kinds)}")
    if family_audit.get("all_development_and_closed_tests_covered") is not True:
        issues.append("family audit does not declare coverage of all development and closed tests")

    models = model_freeze.get("models")
    if model_freeze.get("frozen_before_labels") is not True:
        issues.append("model freeze does not declare frozen_before_labels true")
    if model_freeze.get("labels_available_at_freeze") is not False:
        issues.append("model freeze must declare labels_available_at_freeze false")
    if not isinstance(models, list) or not models:
        issues.append("model freeze must contain at least one model")
        models = []
    model_ids: list[str] = []
    for index, model in enumerate(models):
        if not isinstance(model, dict):
            issues.append(f"model {index} must be an object")
            continue
        model_id = model.get("model_id")
        checkpoint = model.get("checkpoint")
        if not isinstance(model_id, str) or not model_id:
            issues.append(f"model {index} lacks model_id")
        else:
            model_ids.append(model_id)
        if not isinstance(checkpoint, str) or not checkpoint:
            issues.append(f"model {index} lacks checkpoint path")
            continue
        checkpoint_path = (model_freeze_path.parent / checkpoint).resolve()
        if not checkpoint_path.is_file():
            issues.append(f"model {model_id or index} checkpoint is not a file: {checkpoint_path}")
            continue
        _add_identity_issue(issues, model.get("identity"), checkpoint_path, f"model {model_id}")
    if len(model_ids) != len(set(model_ids)):
        issues.append("model_ids must be unique")

    metrics = protocol.get("campaign_metrics")
    supported_metrics = {
        "spearman",
        "pairwise_accuracy",
        "top_10pct_recall",
        "top_10pct_enrichment",
    }
    if not isinstance(metrics, list) or not metrics or not set(metrics) <= supported_metrics:
        issues.append("campaign_metrics must be a nonempty list of supported ranking metrics")

    identities = {
        "protocol": _identity(protocol_path),
        "candidate_roster": {**_identity(candidates_path), "row_count": len(candidates)},
        "campaign_registry": _identity(campaign_registry_path),
        "family_audit": _identity(family_audit_path),
        "model_freeze": _identity(model_freeze_path),
    }
    return {
        "schema_version": 1,
        "stage": "pre_label_readiness_freeze",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "ready": not issues,
        "issues": issues,
        "counts": {
            "campaigns": len(campaigns),
            "families": len(families),
            "candidates": len(candidates),
            "dispositions": dict(
                sorted(Counter(row.get("disposition", "").strip() for row in candidates).items())
            ),
        },
        "campaign_counts": {name: len(rows) for name, rows in sorted(campaigns.items())},
        "model_ids": model_ids,
        "campaign_metrics": metrics if isinstance(metrics, list) else [],
        "identities": identities,
        "claim_boundary": (
            "Identity/readiness manifest only. No labels were read and no model was scored."
        ),
    }


def _verify_readiness(
    readiness_path: Path,
    protocol_path: Path,
    candidates_path: Path,
    campaign_registry_path: Path,
    family_audit_path: Path,
    model_freeze_path: Path,
) -> dict[str, Any]:
    readiness = _load_object(readiness_path)
    if readiness.get("stage") != "pre_label_readiness_freeze" or readiness.get("ready") is not True:
        raise ValueError("readiness manifest is not a passing pre-label freeze")
    paths = {
        "protocol": protocol_path,
        "candidate_roster": candidates_path,
        "campaign_registry": campaign_registry_path,
        "family_audit": family_audit_path,
        "model_freeze": model_freeze_path,
    }
    identities = readiness.get("identities", {})
    for name, path in paths.items():
        expected = identities.get(name)
        actual = _identity(path)
        if name == "candidate_roster" and isinstance(expected, dict):
            expected = {key: expected[key] for key in ("sha256", "size_bytes") if key in expected}
        if expected != actual:
            raise ValueError(f"{name} changed after readiness freeze")
    return readiness


def freeze_prediction_deposit(
    readiness_path: Path,
    protocol_path: Path,
    candidates_path: Path,
    campaign_registry_path: Path,
    family_audit_path: Path,
    model_freeze_path: Path,
    predictions_path: Path,
) -> dict[str, Any]:
    """Validate and freeze a complete blind deposit without reading labels."""
    readiness = _verify_readiness(
        readiness_path,
        protocol_path,
        candidates_path,
        campaign_registry_path,
        family_audit_path,
        model_freeze_path,
    )
    candidates, _ = _read_csv(candidates_path)
    predictions, fields = _read_csv(predictions_path)
    if set(fields) != {"model_id", "candidate_id", "prediction"}:
        raise ValueError(
            "prediction deposit must have exactly model_id,candidate_id,prediction columns"
        )
    expected = {
        (model_id, row["candidate_id"].strip())
        for model_id in readiness["model_ids"]
        for row in candidates
    }
    observed: set[tuple[str, str]] = set()
    for row in predictions:
        key = (row["model_id"].strip(), row["candidate_id"].strip())
        if key in observed:
            raise ValueError(f"duplicate prediction for model/candidate {key}")
        observed.add(key)
        try:
            prediction = float(row["prediction"])
        except ValueError as error:
            raise ValueError(f"non-numeric prediction for model/candidate {key}") from error
        if not math.isfinite(prediction):
            raise ValueError(f"non-finite prediction for model/candidate {key}")
    if observed != expected:
        raise ValueError(
            f"prediction deposit is incomplete or contains unknown rows: "
            f"expected {len(expected)}, observed {len(observed)}"
        )
    return {
        "schema_version": 1,
        "stage": "blind_prediction_deposit",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "readiness_identity": _identity(readiness_path),
        "prediction_identity": {**_identity(predictions_path), "row_count": len(predictions)},
        "models": readiness["model_ids"],
        "candidates": len(candidates),
        "complete_cross_product": True,
        "claim_boundary": "Blind prediction identity only. No labels were read and no metrics computed.",
    }


def final_evaluation(
    readiness_path: Path,
    deposit_path: Path,
    candidates_path: Path,
    predictions_path: Path,
    labels_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    """Perform the only label join and atomically claim the final receipt path."""
    readiness = _load_object(readiness_path)
    deposit = _load_object(deposit_path)
    if readiness.get("ready") is not True:
        raise ValueError("cannot evaluate a benchmark that did not pass readiness")
    if deposit.get("stage") != "blind_prediction_deposit":
        raise ValueError("invalid blind prediction deposit manifest")
    if deposit.get("readiness_identity") != _identity(readiness_path):
        raise ValueError("deposit does not bind the supplied readiness manifest")
    candidate_identity = readiness.get("identities", {}).get("candidate_roster", {})
    if {key: candidate_identity.get(key) for key in ("sha256", "size_bytes")} != _identity(
        candidates_path
    ):
        raise ValueError("candidate roster changed after readiness freeze")
    prediction_identity = deposit.get("prediction_identity", {})
    if {key: prediction_identity.get(key) for key in ("sha256", "size_bytes")} != _identity(
        predictions_path
    ):
        raise ValueError("predictions changed after blind deposit")

    candidates, _ = _read_csv(candidates_path)
    predictions, _ = _read_csv(predictions_path)
    labels, label_fields = _read_csv(labels_path)
    if set(label_fields) != {"candidate_id", "disposition", "label"}:
        raise ValueError("labels must have exactly candidate_id,disposition,label columns")
    roster = {row["candidate_id"].strip(): row for row in candidates}
    if len(labels) != len(roster) or {row["candidate_id"].strip() for row in labels} != set(roster):
        raise ValueError("label file must contain every frozen candidate exactly once")

    numeric_labels: dict[str, float] = {}
    disposition_counts: Counter[str] = Counter()
    for row in labels:
        candidate_id = row["candidate_id"].strip()
        disposition = row["disposition"].strip()
        if disposition != roster[candidate_id]["disposition"].strip():
            raise ValueError(f"disposition changed for candidate {candidate_id}")
        disposition_counts[disposition] += 1
        value = row["label"].strip()
        if value:
            try:
                parsed = float(value)
            except ValueError as error:
                raise ValueError(f"non-numeric label for candidate {candidate_id}") from error
            if not math.isfinite(parsed):
                raise ValueError(f"non-finite label for candidate {candidate_id}")
            numeric_labels[candidate_id] = parsed
        elif disposition != "censored":
            raise ValueError(f"non-censored candidate {candidate_id} lacks a numeric label")

    consumption_path = deposit_path.with_name(f"{deposit_path.name}.consumed.json")
    consumption = {
        "schema_version": 1,
        "stage": "final_evaluation_consumption",
        "deposit_identity": _identity(deposit_path),
        "labels_identity": _identity(labels_path),
        "output": str(output_path.resolve()),
    }
    with consumption_path.open("x", encoding="ascii") as handle:
        json.dump(consumption, handle, indent=2, sort_keys=True)
        handle.write("\n")

    prediction_map = {
        (row["model_id"].strip(), row["candidate_id"].strip()): float(row["prediction"])
        for row in predictions
    }
    model_results: dict[str, Any] = {}
    metric_names = readiness["campaign_metrics"]
    for model_id in readiness["model_ids"]:
        grouped: dict[str, list[str]] = defaultdict(list)
        for candidate_id, row in roster.items():
            if candidate_id in numeric_labels:
                grouped[row["campaign_id"].strip()].append(candidate_id)
        campaign_metrics: dict[str, dict[str, Any]] = {}
        for campaign_id, candidate_ids in sorted(grouped.items()):
            metrics = campaign_ranking_metrics(
                [prediction_map[(model_id, candidate_id)] for candidate_id in candidate_ids],
                [numeric_labels[candidate_id] for candidate_id in candidate_ids],
                seed=int(hashlib.sha256(campaign_id.encode("utf-8")).hexdigest()[:8], 16),
            )
            campaign_metrics[campaign_id] = {
                name: value
                for name, value in metrics.items()
                if name == "count" or name in metric_names
            }
        macro = {
            name: statistics.mean(values)
            if (
                values := [
                    result[name]
                    for result in campaign_metrics.values()
                    if result.get(name) is not None
                ]
            )
            else None
            for name in metric_names
        }
        model_results[model_id] = {
            "campaign_metrics": campaign_metrics,
            "macro_campaign_metrics": macro,
        }

    result = {
        "schema_version": 1,
        "stage": "final_single_join_evaluation",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "readiness_identity": _identity(readiness_path),
        "deposit_identity": _identity(deposit_path),
        "labels_identity": {**_identity(labels_path), "row_count": len(labels)},
        "candidate_dispositions": dict(sorted(disposition_counts.items())),
        "models": model_results,
        "claim_boundary": (
            "Final campaign-level evaluation receipt; it must not be used for model selection."
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output_path.open("x", encoding="ascii") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    except Exception:
        # Keep the consumption marker: a failed final evaluation must fail closed, not retry.
        raise
    return result


def _write_new(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="ascii") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--readiness", type=Path, required=True)
    common.add_argument("--candidates", type=Path, required=True)

    audit = commands.add_parser("audit-freeze")
    audit.add_argument("--protocol", type=Path, required=True)
    audit.add_argument("--candidates", type=Path, required=True)
    audit.add_argument("--campaign-registry", type=Path, required=True)
    audit.add_argument("--family-audit", type=Path, required=True)
    audit.add_argument("--model-freeze", type=Path, required=True)
    audit.add_argument("--output", type=Path, required=True)

    deposit = commands.add_parser("deposit-freeze", parents=[common])
    deposit.add_argument("--protocol", type=Path, required=True)
    deposit.add_argument("--campaign-registry", type=Path, required=True)
    deposit.add_argument("--family-audit", type=Path, required=True)
    deposit.add_argument("--model-freeze", type=Path, required=True)
    deposit.add_argument("--predictions", type=Path, required=True)
    deposit.add_argument("--output", type=Path, required=True)

    evaluate = commands.add_parser("evaluate-once", parents=[common])
    evaluate.add_argument("--deposit", type=Path, required=True)
    evaluate.add_argument("--predictions", type=Path, required=True)
    evaluate.add_argument("--labels", type=Path, required=True)
    evaluate.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.command == "audit-freeze":
        payload = audit_readiness(
            args.protocol,
            args.candidates,
            args.campaign_registry,
            args.family_audit,
            args.model_freeze,
        )
        _write_new(args.output, payload)
        return 0 if payload["ready"] else 1
    if args.command == "deposit-freeze":
        payload = freeze_prediction_deposit(
            args.readiness,
            args.protocol,
            args.candidates,
            args.campaign_registry,
            args.family_audit,
            args.model_freeze,
            args.predictions,
        )
        _write_new(args.output, payload)
        return 0
    final_evaluation(
        args.readiness, args.deposit, args.candidates, args.predictions, args.labels, args.output
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
