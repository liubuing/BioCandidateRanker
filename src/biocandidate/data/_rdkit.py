from __future__ import annotations


class RDKitUnavailableError(RuntimeError):
    pass


def require_rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RDKitUnavailableError(
            "RDKit is required for molecular graphs and scaffold-based splits; "
            "install RDKit in this environment before using this operation"
        ) from exc
    return Chem


def canonicalize_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem
    except ImportError:
        return smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse SMILES")
    return Chem.MolToSmiles(molecule, canonical=True)
