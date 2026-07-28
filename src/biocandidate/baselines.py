from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from .data.schema import EnzymeSubstrateRecord


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def amino_acid_composition(records: list[EnzymeSubstrateRecord]) -> np.ndarray:
    features = np.zeros((len(records), len(AMINO_ACIDS)), dtype=np.float32)
    indices = {residue: index for index, residue in enumerate(AMINO_ACIDS)}
    for row, record in enumerate(records):
        for residue in record.sequence:
            if residue in indices:
                features[row, indices[residue]] += 1.0
        features[row] /= len(record.sequence)
    return features


def morgan_fingerprints(
    records: list[EnzymeSubstrateRecord], *, radius: int = 2, bits: int = 2048,
) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    features = np.empty((len(records), bits), dtype=np.uint8)
    for row, record in enumerate(records):
        molecule = Chem.MolFromSmiles(record.substrate_smiles)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse substrate at source row {record.source_row}")
        features[row] = generator.GetFingerprintAsNumPy(molecule)
    return features
