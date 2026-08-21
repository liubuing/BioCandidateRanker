from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Tuple, Union

from ._rdkit import require_rdkit
from .schema import EnzymeSubstrateRecord, EvidenceTier


PathLike = Union[str, Path]
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
DATASET_DOI = "10.5281/zenodo.17310823"
DATASET_LICENSE = "CC-BY-4.0"
ARCHIVE_SHA256 = "8013ad81586db2187162aada0709c1cabc7e7d69f03dd5c199776aaf000dd6ea"
ARCHIVE_SIZE_BYTES = 30_050_527


@dataclass(frozen=True, slots=True)
class EnzEngDBRejectedRow:
    source_file: str
    source_row: int
    reason: str


@dataclass(frozen=True, slots=True)
class EnzEngDBObservation:
    record: EnzymeSubstrateRecord
    campaign_id: str
    fitness_value: float
    endpoint_name: str
    source_file: str


@dataclass(frozen=True, slots=True)
class EnzEngDBAuditReport:
    source_files: int
    total_rows: int
    accepted_rows: int
    campaigns: int
    rejected: Tuple[EnzEngDBRejectedRow, ...]


@dataclass(frozen=True, slots=True)
class EnzEngDBReadResult:
    observations: Tuple[EnzEngDBObservation, ...]
    audit: EnzEngDBAuditReport


def campaign_rank_records(
    observations: Iterable[EnzEngDBObservation],
) -> tuple[EnzymeSubstrateRecord, ...]:
    """Attach average-tie percentile ranks while preserving observation order."""
    items = tuple(observations)
    campaign_indices: dict[str, list[int]] = defaultdict(list)
    for index, observation in enumerate(items):
        campaign_indices[observation.campaign_id].append(index)

    ranks = [0.0] * len(items)
    for indices in campaign_indices.values():
        ordered = sorted(indices, key=lambda index: items[index].fitness_value)
        start = 0
        while start < len(ordered):
            end = start + 1
            value = items[ordered[start]].fitness_value
            while end < len(ordered) and items[ordered[end]].fitness_value == value:
                end += 1
            percentile = 0.5 if len(ordered) == 1 else ((start + end - 1) / 2) / (len(ordered) - 1)
            for position in range(start, end):
                ranks[ordered[position]] = percentile
            start = end

    return tuple(
        replace(
            observation.record,
            activity_rank=ranks[index],
            campaign_group=observation.campaign_id,
        )
        for index, observation in enumerate(items)
    )


def campaign_representative_records(
    observations: Iterable[EnzEngDBObservation],
) -> tuple[EnzymeSubstrateRecord, ...]:
    """Select one actual modal-length sequence nearest consensus per campaign."""
    campaigns: dict[str, list[EnzEngDBObservation]] = defaultdict(list)
    for observation in observations:
        campaigns[observation.campaign_id].append(observation)

    representatives = []
    for campaign_id in sorted(campaigns):
        campaign = campaigns[campaign_id]
        length_counts = Counter(len(item.record.sequence) for item in campaign)
        modal_length = min(length_counts, key=lambda length: (-length_counts[length], length))
        candidates = [item for item in campaign if len(item.record.sequence) == modal_length]
        consensus = "".join(
            min(Counter(item.record.sequence[position] for item in candidates).items(),
                key=lambda pair: (-pair[1], pair[0]))[0]
            for position in range(modal_length)
        )

        def representative_key(item: EnzEngDBObservation) -> tuple:
            record = item.record
            distance = sum(left != right for left, right in zip(record.sequence, consensus))
            return (
                distance,
                record.sequence,
                record.candidate_id,
                record.source_dataset,
                record.source_row,
                repr(record),
            )

        selected = min(candidates, key=representative_key).record
        representatives.append(replace(selected, campaign_group=campaign_id))
    return tuple(representatives)


def _first(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and value.strip():
            return value.strip()
    return ""


@lru_cache(maxsize=4096)
def select_primary_reactant(reaction_smiles: str) -> str:
    """Select a deterministic molecular input from an EnzEngDB reaction.

    EnzEngDB campaigns use one reaction for all variants. The largest concrete
    carbon-containing reactant is retained; wildcard-only assay definitions are
    rejected because they cannot form a meaningful molecular graph.
    """
    if reaction_smiles.count(">>") != 1:
        raise ValueError("reaction SMILES must contain exactly one >> separator")
    reactants, _ = reaction_smiles.split(">>", 1)
    Chem = require_rdkit()
    from rdkit import rdBase

    candidates = []
    for component in reactants.split("."):
        component = component.strip()
        if not component:
            continue
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(component)
        if molecule is None:
            continue
        atoms = list(molecule.GetAtoms())
        if any(atom.GetAtomicNum() == 0 for atom in atoms):
            continue
        carbon_count = sum(atom.GetAtomicNum() == 6 for atom in atoms)
        if carbon_count == 0:
            continue
        canonical = Chem.MolToSmiles(molecule, canonical=True)
        candidates.append((molecule.GetNumHeavyAtoms(), carbon_count, canonical))
    if not candidates:
        raise ValueError("reaction has no concrete carbon-containing reactant")
    return max(candidates)[2]


def _read_campaign(path: Path) -> tuple[list[EnzEngDBObservation], list[EnzEngDBRejectedRow], int]:
    observations = []
    rejected = []
    total = 0
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read EnzEngDB CSV {path}: {exc}") from exc
    with handle:
        try:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ValueError(f"EnzEngDB CSV has no header: {path}")
            for source_row, row in enumerate(reader, start=1):
                if not any((value or "").strip() for value in row.values()):
                    continue
                total += 1
                sequence = _first(row, "aa_sequence")
                if not sequence:
                    rejected.append(EnzEngDBRejectedRow(path.name, source_row, "aa_sequence is missing"))
                    continue
                sequence = "".join(sequence.split()).upper()
                invalid = sorted(set(sequence) - CANONICAL_AMINO_ACIDS)
                if invalid:
                    rejected.append(EnzEngDBRejectedRow(
                        path.name, source_row,
                        "aa_sequence contains non-canonical residues: " + "".join(invalid)))
                    continue
                raw_fitness = _first(row, "fitness_value")
                try:
                    fitness = float(raw_fitness)
                except (TypeError, ValueError):
                    rejected.append(EnzEngDBRejectedRow(
                        path.name, source_row, "fitness_value is missing or non-numeric"))
                    continue
                if not math.isfinite(fitness):
                    rejected.append(EnzEngDBRejectedRow(
                        path.name, source_row, "fitness_value must be finite"))
                    continue
                reaction = _first(row, "reaction_smiles", "smiles_reaction")
                try:
                    substrate = select_primary_reactant(reaction)
                except ValueError as exc:
                    rejected.append(EnzEngDBRejectedRow(path.name, source_row, str(exc)))
                    continue
                raw_id = _first(row, "id") or path.stem
                candidate_id = f"{path.stem}:{source_row}:{raw_id}"
                substitutions = _first(row, "amino_acid_substitutions")
                enzyme_type = "parent" if substitutions == "#PARENT#" else "engineered_variant"
                record = EnzymeSubstrateRecord(
                    sequence=sequence,
                    substrate_smiles=substrate,
                    organism="",
                    ec="",
                    enzyme_type=enzyme_type,
                    candidate_id=candidate_id,
                    reaction=reaction,
                    evidence_tier=EvidenceTier.CURATED,
                    source_dataset=f"EnzEngDB v1 ({DATASET_DOI})",
                    source_row=source_row,
                )
                observations.append(EnzEngDBObservation(
                    record=record,
                    campaign_id=path.stem,
                    fitness_value=fitness,
                    endpoint_name="fitness_value",
                    source_file=path.name,
                ))
        except (csv.Error, UnicodeError) as exc:
            raise ValueError(f"cannot parse EnzEngDB CSV {path}: {exc}") from exc
    return observations, rejected, total


def read_enzengdb_directory(path: PathLike) -> EnzEngDBReadResult:
    root = Path(path)
    if not root.is_dir():
        raise ValueError(f"EnzEngDB experiments path is not a directory: {root}")
    files = sorted(root.glob("*.csv"), key=lambda item: item.name.lower())
    if not files:
        raise ValueError(f"EnzEngDB experiments directory has no CSV files: {root}")
    observations = []
    rejected = []
    total_rows = 0
    for source in files:
        accepted, campaign_rejected, campaign_total = _read_campaign(source)
        observations.extend(accepted)
        rejected.extend(campaign_rejected)
        total_rows += campaign_total
    campaign_count = len({item.campaign_id for item in observations})
    return EnzEngDBReadResult(
        observations=tuple(observations),
        audit=EnzEngDBAuditReport(
            source_files=len(files),
            total_rows=total_rows,
            accepted_rows=len(observations),
            campaigns=campaign_count,
            rejected=tuple(rejected),
        ),
    )
