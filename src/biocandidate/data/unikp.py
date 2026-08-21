from __future__ import annotations

import json
import math
import statistics
from dataclasses import replace
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

from ._rdkit import canonicalize_smiles
from .schema import EnzymeSubstrateRecord, EvidenceTier


PathLike = Union[str, Path]


@dataclass(frozen=True, slots=True)
class RejectedRow:
    source_row: int
    reason: str


@dataclass(frozen=True, slots=True)
class DuplicateRow:
    source_row: int
    first_source_row: int
    conflicting_value: bool


@dataclass(frozen=True, slots=True)
class UniKPAuditReport:
    total_rows: int
    accepted_rows: int
    rejected: Tuple[RejectedRow, ...]
    duplicates: Tuple[DuplicateRow, ...]

    @property
    def conflict_count(self) -> int:
        return sum(item.conflicting_value for item in self.duplicates)


@dataclass(frozen=True, slots=True)
class UniKPReadResult:
    records: Tuple[EnzymeSubstrateRecord, ...]
    audit: UniKPAuditReport


def aggregate_pair_measurements(
    records: Sequence[EnzymeSubstrateRecord],
) -> Tuple[EnzymeSubstrateRecord, ...]:
    """Aggregate repeated sequence-SMILES kcat measurements on the log scale."""
    groups: dict[tuple[str, str], list[EnzymeSubstrateRecord]] = {}
    for record in records:
        groups.setdefault((record.sequence, record.substrate_smiles), []).append(record)
    aggregated = []
    for _, group in sorted(groups.items(), key=lambda item: min(r.source_row for r in item[1])):
        values = [record.log10_kcat for record in group if record.log10_kcat is not None]
        first = min(group, key=lambda record: record.source_row)
        source_rows = tuple(sorted(record.source_row for record in group))
        aggregated.append(replace(
            first,
            log10_kcat=statistics.median(values) if values else None,
            source_rows=source_rows,
            replicate_count=len(group),
            label_stddev=statistics.stdev(values) if len(values) > 1 else 0.0,
        ))
    return tuple(aggregated)


def _first(row: Mapping[str, Any], names: Sequence[str], default: Any = "") -> Any:
    for name in names:
        if name in row and row[name] is not None:
            return row[name]
    return default


def _load_rows(path: Path) -> list[Mapping[str, Any]]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read UniKP JSON {path}: {exc}") from exc
    if isinstance(payload, dict):
        payload = next(
            (payload[key] for key in ("data", "records", "rows") if isinstance(payload.get(key), list)),
            None,
        )
    if not isinstance(payload, list) or not all(isinstance(row, dict) for row in payload):
        raise ValueError("UniKP JSON must contain an array of record objects")
    return payload


def read_unikp_json(
    path: PathLike,
    *,
    source_dataset: str = "UniKP",
    evidence_tier: EvidenceTier = EvidenceTier.UNKNOWN,
    split: Optional[str] = None,
    source_rows: set[int] | None = None,
) -> UniKPReadResult:
    rows = _load_rows(Path(path))
    records = []
    rejected = []
    duplicates = []
    first_seen: dict[tuple[str, ...], tuple[int, float]] = {}

    for source_row, row in enumerate(rows):
        # Partition filtering precedes all label access, including type conversion.
        if source_rows is not None and source_row not in source_rows:
            continue
        try:
            value = float(_first(row, ("Value", "value", "kcat")))
        except (TypeError, ValueError):
            rejected.append(RejectedRow(source_row, "Value is missing or non-numeric"))
            continue
        if not math.isfinite(value) or value <= 0:
            rejected.append(RejectedRow(source_row, "Value must be finite and greater than zero"))
            continue

        sequence = str(_first(row, ("Sequence", "sequence"))).strip()
        smiles = str(_first(row, ("Smiles", "SMILES", "smiles", "substrate_smiles"))).strip()
        if not sequence:
            rejected.append(RejectedRow(source_row, "Sequence is missing"))
            continue
        if not smiles:
            rejected.append(RejectedRow(source_row, "substrate SMILES is missing"))
            continue
        if "." in smiles:
            rejected.append(RejectedRow(source_row, "multi-component substrate SMILES is not supported"))
            continue
        try:
            smiles = canonicalize_smiles(smiles)
        except ValueError as exc:
            rejected.append(RejectedRow(source_row, str(exc)))
            continue

        organism = str(_first(row, ("Organism", "organism"))).strip()
        ec = str(_first(row, ("ECNumber", "EC", "ec"))).strip()
        enzyme_type = str(
            _first(row, ("Type", "EnzymeType", "enzyme_type"), "unknown")
        ).strip()
        substrate_name = str(_first(row, ("Substrate", "substrate_name"))).strip()
        record = EnzymeSubstrateRecord(
            sequence=sequence,
            substrate_smiles=smiles,
            organism=organism,
            ec=ec,
            enzyme_type=enzyme_type,
            substrate_name=substrate_name,
            log10_kcat=math.log10(value),
            evidence_tier=evidence_tier,
            source_dataset=source_dataset,
            source_row=source_row,
            split=split,
        )
        # Leakage and conflicting-label audits operate at the actual model input
        # pair, not at metadata-decorated rows.
        key = (record.sequence, smiles)
        if key in first_seen:
            first_row, first_value = first_seen[key]
            duplicates.append(DuplicateRow(source_row, first_row, value != first_value))
        else:
            first_seen[key] = (source_row, value)
        records.append(record)

    return UniKPReadResult(
        records=tuple(records),
        audit=UniKPAuditReport(
            total_rows=len(rows),
            accepted_rows=len(records),
            rejected=tuple(rejected),
            duplicates=tuple(duplicates),
        ),
    )
