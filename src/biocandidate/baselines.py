"""Classical feature builders for the non-neural ridge baselines.

Deterministic, label-free featurization shared by every baseline run so the
shallow models and the deep ranker are scored on byte-identical rows: relative
amino-acid composition of the protein sequence, and binary Morgan fingerprints
of the substrate SMILES. Unparseable substrate structures are a hard error, not
a skipped row.

中文：面向非神经 ridge 基线的经典特征构建器。提供确定性的、不依赖标签的特征化
逻辑，由所有基线运行共享，从而保证浅层模型与深度排序器在字节级相同的行上打分：
蛋白质序列的相对氨基酸组成，以及底物 SMILES 的二值 Morgan 指纹。无法解析的
底物结构被视为硬错误，而不是跳过该行。
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from .data.schema import EnzymeSubstrateRecord

# The 20 standard residues in a fixed column order. Any other character in a
# sequence adds to no column but still enters the len(record.sequence)
# denominator, so composition rows sum to slightly below one when present.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def amino_acid_composition(records: list[EnzymeSubstrateRecord]) -> np.ndarray:
    """Relative amino-acid frequencies over the 20 standard residues, one row per record.

    中文：计算 20 种标准氨基酸的相对频率，每条记录一行。
    """
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
    """Binary Morgan fingerprints of each substrate SMILES, one row per record.

    中文：为每条底物 SMILES 生成二值 Morgan 指纹，每条记录一行。
    """
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    features = np.empty((len(records), bits), dtype=np.uint8)
    for row, record in enumerate(records):
        molecule = Chem.MolFromSmiles(record.substrate_smiles)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse substrate at source row {record.source_row}")
        features[row] = generator.GetFingerprintAsNumPy(molecule)
    return features
