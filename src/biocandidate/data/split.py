from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from typing import Sequence

from ._rdkit import require_rdkit
from .schema import EnzymeSubstrateRecord


SPLIT_NAMES = ("train", "validation", "test")


def molecular_scaffold(smiles: str) -> str:
    Chem = require_rdkit()
    try:
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:
        raise RuntimeError("installed RDKit does not provide Murcko scaffold support") from exc
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return scaffold or Chem.MolToSmiles(molecule, canonical=True)


def _components(records: Sequence[EnzymeSubstrateRecord], scaffolds: Sequence[str]) -> list[list[int]]:
    adjacency: dict[str, set[str]] = {}
    rows_by_node: dict[str, list[int]] = {}
    for index, (record, scaffold) in enumerate(zip(records, scaffolds)):
        protein_node = "p:" + record.sequence
        scaffold_node = "s:" + scaffold
        adjacency.setdefault(protein_node, set()).add(scaffold_node)
        adjacency.setdefault(scaffold_node, set()).add(protein_node)
        rows_by_node.setdefault(protein_node, []).append(index)
        rows_by_node.setdefault(scaffold_node, []).append(index)

    components = []
    visited = set()
    for start in sorted(adjacency):
        if start in visited:
            continue
        stack = [start]
        indices = set()
        visited.add(start)
        while stack:
            node = stack.pop()
            indices.update(rows_by_node[node])
            for neighbor in adjacency[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        components.append(sorted(indices))
    return components


def _group_rows(
    records: Sequence[EnzymeSubstrateRecord], strategy: str
) -> list[list[int]]:
    if strategy not in {"protein_cold", "scaffold_cold", "double_cold"}:
        raise ValueError("strategy must be protein_cold, scaffold_cold, or double_cold")
    if strategy == "protein_cold":
        def key(record: EnzymeSubstrateRecord) -> str:
            return record.sequence
        groups: dict[str, list[int]] = {}
        for index, record in enumerate(records):
            groups.setdefault(key(record), []).append(index)
        return list(groups.values())

    scaffolds = [molecular_scaffold(record.substrate_smiles) for record in records]
    if strategy == "scaffold_cold":
        groups = {}
        for index, scaffold in enumerate(scaffolds):
            groups.setdefault(scaffold, []).append(index)
        return list(groups.values())
    if strategy == "double_cold":
        return _components(records, scaffolds)
    raise AssertionError("unreachable")


def split_records(
    records: Sequence[EnzymeSubstrateRecord],
    strategy: str,
    *,
    ratios: tuple[float, float, float] = (0.8, 0.1, 0.1),
    seed: int = 0,
) -> tuple[EnzymeSubstrateRecord, ...]:
    if len(ratios) != 3 or any(ratio < 0 for ratio in ratios) or abs(sum(ratios) - 1.0) > 1e-9:
        raise ValueError("ratios must contain three non-negative values summing to one")
    if not records:
        return ()

    groups = _group_rows(records, strategy)
    groups.sort(
        key=lambda group: (
            -len(group),
            hashlib.sha256(f"{seed}:{','.join(map(str, group))}".encode("ascii")).hexdigest(),
        )
    )
    targets = [ratio * len(records) for ratio in ratios]
    counts = [0, 0, 0]
    assignment = [""] * len(records)
    for group in groups:
        candidates = [index for index, ratio in enumerate(ratios) if ratio > 0]
        split_index = min(
            candidates,
            key=lambda index: (
                counts[index] / targets[index] if targets[index] else float("inf"),
                counts[index],
                index,
            ),
        )
        for row in group:
            assignment[row] = SPLIT_NAMES[split_index]
        counts[split_index] += len(group)
    return tuple(replace(record, split=assignment[index]) for index, record in enumerate(records))


def write_cold_split_manifest(records: Sequence[EnzymeSubstrateRecord], strategy: str,
                              path: str, *, seed: int, source_identity: dict) -> dict:
    scaffolds = [molecular_scaffold(record.substrate_smiles) for record in records]
    rows = []
    split_counts: dict[str, int] = {}
    proteins_by_split: dict[str, set[str]] = {}
    scaffolds_by_split: dict[str, set[str]] = {}
    for record, scaffold in zip(records, scaffolds):
        split = record.split or "unset"
        split_counts[split] = split_counts.get(split, 0) + 1
        protein_hash = hashlib.sha256(record.sequence.encode("ascii")).hexdigest()
        proteins_by_split.setdefault(split, set()).add(protein_hash)
        scaffolds_by_split.setdefault(split, set()).add(scaffold)
        rows.append({
            "source_row": record.source_row,
            "split": record.split,
            "protein_sha256": protein_hash,
            "scaffold": scaffold,
        })

    names = sorted(split_counts)
    protein_crossings = sum(
        len(proteins_by_split[names[i]] & proteins_by_split[names[j]])
        for i in range(len(names)) for j in range(i + 1, len(names)))
    scaffold_crossings = sum(
        len(scaffolds_by_split[names[i]] & scaffolds_by_split[names[j]])
        for i in range(len(names)) for j in range(i + 1, len(names)))
    payload = {
        "format_version": 1,
        "parameters": {
            "method": strategy,
            "seed": seed,
            "source_identity": source_identity,
        },
        "split_counts": split_counts,
        "audit": {
            "protein_cross_split_count": protein_crossings,
            "scaffold_cross_split_count": scaffold_crossings,
        },
        "rows": rows,
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return payload
