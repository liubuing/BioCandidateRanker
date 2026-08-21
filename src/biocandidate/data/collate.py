from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

import torch

from ._rdkit import require_rdkit
from .schema import EnzymeSubstrateRecord, TASK_NAMES


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class AminoAcidTokenizer:
    pad_id = 0
    unknown_id = 1

    def __init__(self) -> None:
        self._ids = {amino_acid: index + 2 for index, amino_acid in enumerate(AMINO_ACIDS)}

    @property
    def vocabulary_size(self) -> int:
        return len(self._ids) + 2

    @lru_cache(maxsize=16384)
    def encode(self, sequence: str) -> tuple[int, ...]:
        return tuple(self._ids.get(amino_acid, self.unknown_id) for amino_acid in sequence.upper())

    def pad(self, sequences: Sequence[str]) -> tuple[torch.Tensor, torch.Tensor]:
        encoded = [self.encode(sequence) for sequence in sequences]
        width = max((len(item) for item in encoded), default=0)
        tokens = torch.full((len(encoded), width), self.pad_id, dtype=torch.long)
        mask = torch.zeros((len(encoded), width), dtype=torch.bool)
        for row, item in enumerate(encoded):
            if item:
                tokens[row, : len(item)] = torch.tensor(item, dtype=torch.long)
                mask[row, : len(item)] = True
        return tokens, mask


@dataclass(frozen=True, slots=True)
class MolecularGraph:
    atom_features: torch.Tensor
    edge_index: torch.Tensor
    edge_features: torch.Tensor


@lru_cache(maxsize=8192)
def smiles_to_sparse_graph(smiles: str) -> MolecularGraph:
    Chem = require_rdkit()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    features = [
        [
            atom.GetAtomicNum(),
            atom.GetDegree(),
            atom.GetFormalCharge(),
            int(atom.GetIsAromatic()),
            int(atom.GetHybridization()),
            atom.GetTotalNumHs(),
        ]
        for atom in molecule.GetAtoms()
    ]
    edges = []
    edge_features = []
    for bond in molecule.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edges.extend(((start, end), (end, start)))
        bond_type = {
            "SINGLE": 1,
            "DOUBLE": 2,
            "TRIPLE": 3,
            "AROMATIC": 4,
        }.get(str(bond.GetBondType()), 0)
        bond_features = (bond_type, int(bond.GetIsConjugated()), int(bond.IsInRing()))
        edge_features.extend((bond_features, bond_features))
    atom_features = torch.tensor(features, dtype=torch.long).reshape(-1, 6)
    edge_index = torch.tensor(edges, dtype=torch.long).reshape(-1, 2).t().contiguous()
    edge_features_tensor = torch.tensor(edge_features, dtype=torch.long).reshape(-1, 3)
    return MolecularGraph(atom_features, edge_index, edge_features_tensor)


def _hash_id(value: str, buckets: int) -> int:
    if not value:
        return 0
    digest = hashlib.sha256(value.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (buckets - 1) + 1


class EnzymeSubstrateCollator:
    def __init__(
        self,
        *,
        context_buckets: int = 4096,
        fba_context_dim: int = 8,
        task_names: Sequence[str] = TASK_NAMES,
    ) -> None:
        if context_buckets < 2:
            raise ValueError("context_buckets must be at least 2")
        if not task_names:
            raise ValueError("task_names must not be empty")
        if fba_context_dim < 1:
            raise ValueError("fba_context_dim must be positive")
        self.context_buckets = context_buckets
        self.fba_context_dim = fba_context_dim
        self.task_names = tuple(task_names)
        self.tokenizer = AminoAcidTokenizer()

    def __call__(self, records: Sequence[EnzymeSubstrateRecord]) -> dict[str, torch.Tensor]:
        if not records:
            raise ValueError("cannot collate an empty batch")
        sequence_tokens, sequence_mask = self.tokenizer.pad([item.sequence for item in records])
        graphs = [smiles_to_sparse_graph(item.substrate_smiles) for item in records]
        atom_features = []
        edge_indices = []
        edge_features = []
        graph_batch = []
        offset = 0
        for graph_id, graph in enumerate(graphs):
            atom_features.append(graph.atom_features)
            edge_indices.append(graph.edge_index + offset)
            edge_features.append(graph.edge_features)
            graph_batch.append(torch.full((len(graph.atom_features),), graph_id, dtype=torch.long))
            offset += len(graph.atom_features)

        labels = torch.zeros((len(records), len(self.task_names)), dtype=torch.float32)
        label_mask = torch.zeros_like(labels, dtype=torch.bool)
        for row, record in enumerate(records):
            for column, task in enumerate(self.task_names):
                value = getattr(record, task)
                if value is not None:
                    labels[row, column] = value
                    label_mask[row, column] = True

        context_ids = torch.tensor(
            [
                [
                    _hash_id(record.organism, self.context_buckets),
                    _hash_id(record.ec, self.context_buckets),
                    _hash_id(record.enzyme_type, self.context_buckets),
                    _hash_id(record.reaction, self.context_buckets),
                ]
                for record in records
            ],
            dtype=torch.long,
        )
        fba_context = torch.zeros(len(records), self.fba_context_dim, dtype=torch.float32)
        fba_context_mask = torch.zeros(len(records), dtype=torch.bool)
        for row, record in enumerate(records):
            width = len(record.fba_context)
            if width:
                if width != self.fba_context_dim:
                    raise ValueError(
                        f"record {row} FBA context width {width} does not match configured "
                        f"width {self.fba_context_dim}"
                    )
                fba_context[row] = torch.tensor(record.fba_context, dtype=torch.float32)
                fba_context_mask[row] = True
        evidence_weights = {
            "direct": 1.0,
            "curated": 0.75,
            "inferred": 0.25,
            "unknown": 0.1,
        }
        return {
            "sequence_tokens": sequence_tokens,
            "sequence_mask": sequence_mask,
            "atom_features": torch.cat(atom_features),
            "edge_index": torch.cat(edge_indices, dim=1),
            "edge_features": torch.cat(edge_features, dim=0),
            "graph_batch": torch.cat(graph_batch),
            "context_ids": context_ids,
            "fba_context": fba_context,
            "fba_context_mask": fba_context_mask,
            "labels": labels,
            "label_mask": label_mask,
            "campaign_ids": torch.tensor(
                [_stable_campaign_id(item.campaign_group) for item in records],
                dtype=torch.int64,
            ),
            "evidence_weight": torch.tensor(
                [evidence_weights[item.evidence_tier.value] for item in records],
                dtype=torch.float32,
            ),
        }


def _stable_campaign_id(value: str) -> int:
    if not value:
        return 0
    digest = hashlib.sha256(value.encode("utf-8")).digest()
    return (int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)) or 1
