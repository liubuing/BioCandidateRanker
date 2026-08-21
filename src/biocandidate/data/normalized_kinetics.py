from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Union

from .schema import EnzymeSubstrateRecord, EvidenceTier, KineticLabelMetadata


PathLike = Union[str, Path]
NORMALIZED_KINETICS_SCHEMA_VERSION = 1
CANONICAL_KINETIC_UNITS = {
    "log10_kcat": "log10(s^-1)",
    "log10_km": "log10(M)",
    "log10_kcat_per_km": "log10(M^-1 s^-1)",
}

# Factors convert positive linear source values to the physical canonical unit.
_LINEAR_UNIT_FACTORS = {
    "log10_kcat": {
        "s^-1": 1.0,
        "min^-1": 1.0 / 60.0,
        "h^-1": 1.0 / 3600.0,
    },
    "log10_km": {
        "M": 1.0,
        "mM": 1e-3,
        "uM": 1e-6,
        "nM": 1e-9,
    },
    "log10_kcat_per_km": {
        "M^-1 s^-1": 1.0,
        "mM^-1 s^-1": 1e3,
        "uM^-1 s^-1": 1e6,
        "M^-1 min^-1": 1.0 / 60.0,
        "mM^-1 min^-1": 1e3 / 60.0,
        "uM^-1 min^-1": 1e6 / 60.0,
    },
}


@dataclass(frozen=True, slots=True)
class NormalizedKineticsReadResult:
    records: tuple[EnzymeSubstrateRecord, ...]
    schema_version: int = NORMALIZED_KINETICS_SCHEMA_VERSION


def _required_text(row: Mapping[str, Any], name: str, line_number: int) -> str:
    value = row.get(name)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"line {line_number}: {name} must be a non-empty string")
    return value.strip()


def _normalize_label(
    task: str,
    label: Any,
    line_number: int,
) -> tuple[float, KineticLabelMetadata]:
    if task not in CANONICAL_KINETIC_UNITS:
        raise ValueError(
            f"line {line_number}: unsupported kinetic task {task!r}; activity and flux "
            "measurements must not be represented as absolute kinetics"
        )
    if not isinstance(label, dict):
        raise ValueError(f"line {line_number}: label {task!r} must be an object")
    if isinstance(label.get("value"), bool):
        raise ValueError(f"line {line_number}: label {task!r} value must be numeric")
    try:
        value = float(label["value"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"line {line_number}: label {task!r} value must be numeric") from exc
    if not math.isfinite(value):
        raise ValueError(f"line {line_number}: label {task!r} value must be finite")

    unit = label.get("unit")
    if not isinstance(unit, str):
        raise ValueError(f"line {line_number}: label {task!r} unit must be a string")
    canonical_unit = CANONICAL_KINETIC_UNITS[task]
    if unit == canonical_unit:
        normalized_value = value
    else:
        factor = _LINEAR_UNIT_FACTORS[task].get(unit)
        if factor is None:
            allowed = sorted((*_LINEAR_UNIT_FACTORS[task], canonical_unit))
            raise ValueError(
                f"line {line_number}: unsupported unit {unit!r} for {task}; "
                f"expected one of {allowed}"
            )
        if value <= 0:
            raise ValueError(
                f"line {line_number}: linear label {task!r} must be greater than zero")
        normalized_value = math.log10(value * factor)

    evidence = label.get("evidence")
    if evidence != EvidenceTier.DIRECT.value:
        raise ValueError(
            f"line {line_number}: label {task!r} evidence must be 'direct', got {evidence!r}")
    provenance = label.get("provenance")
    if not isinstance(provenance, dict) or not provenance:
        raise ValueError(
            f"line {line_number}: label {task!r} provenance must be a non-empty object")
    try:
        provenance_json = json.dumps(provenance, sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"line {line_number}: label {task!r} provenance must be JSON-serializable") from exc
    return normalized_value, KineticLabelMetadata(
        task=task,
        canonical_unit=canonical_unit,
        evidence_tier=EvidenceTier.DIRECT,
        provenance_json=provenance_json,
    )


def read_normalized_kinetics_jsonl(path: PathLike) -> NormalizedKineticsReadResult:
    """Read strict v1 absolute-kinetics JSONL into sparse model records.

    Each row has ``schema_version: 1``, model input fields, and a non-empty ``labels``
    object keyed by canonical model task. Label objects contain ``value``, ``unit``,
    ``evidence: \"direct\"``, and a non-empty ``provenance`` object.
    """
    records = []
    try:
        handle = Path(path).open("r", encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read normalized kinetics JSONL {path}: {exc}") from exc

    with handle:
        for line_number, text in enumerate(handle, start=1):
            if not text.strip():
                raise ValueError(f"line {line_number}: blank lines are not allowed")
            try:
                row = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"line {line_number}: invalid JSON: {exc.msg}") from exc
            if not isinstance(row, dict):
                raise ValueError(f"line {line_number}: row must be an object")
            schema_version = row.get("schema_version")
            if (
                isinstance(schema_version, bool)
                or schema_version != NORMALIZED_KINETICS_SCHEMA_VERSION
            ):
                raise ValueError(
                    f"line {line_number}: schema_version must be "
                    f"{NORMALIZED_KINETICS_SCHEMA_VERSION}"
                )
            labels = row.get("labels")
            if not isinstance(labels, dict) or not labels:
                raise ValueError(f"line {line_number}: labels must be a non-empty object")

            normalized_labels: dict[str, float] = {}
            metadata = []
            for task, label in labels.items():
                value, item_metadata = _normalize_label(task, label, line_number)
                normalized_labels[task] = value
                metadata.append(item_metadata)

            source_row = row.get("source_row", line_number - 1)
            if isinstance(source_row, bool) or not isinstance(source_row, int):
                raise ValueError(f"line {line_number}: source_row must be an integer")
            evidence_tier = EvidenceTier.DIRECT
            records.append(EnzymeSubstrateRecord(
                sequence=_required_text(row, "sequence", line_number),
                substrate_smiles=_required_text(row, "substrate_smiles", line_number),
                organism=str(row.get("organism", "")),
                ec=str(row.get("ec", "")),
                enzyme_type=str(row.get("enzyme_type", "unknown")),
                candidate_id=str(row.get("candidate_id", "")),
                substrate_name=str(row.get("substrate_name", "")),
                reaction=str(row.get("reaction", "")),
                evidence_tier=evidence_tier,
                source_dataset=_required_text(row, "source_dataset", line_number),
                source_row=source_row,
                split=row.get("split"),
                kinetic_label_metadata=tuple(metadata),
                **normalized_labels,
            ))
    return NormalizedKineticsReadResult(records=tuple(records))
