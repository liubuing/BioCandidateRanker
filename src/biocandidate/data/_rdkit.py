"""Lazily guarded RDKit access for molecular graphs and scaffold operations.

RDKit is an optional dependency; every entry point funnels through
:func:`require_rdkit` so a missing installation surfaces as one consistent,
fail-closed error instead of an import-time crash.
"""
from __future__ import annotations


class RDKitUnavailableError(RuntimeError):
    """Raised when an operation needs RDKit but it is not installed."""


def require_rdkit():
    """Return the ``rdkit.Chem`` module, or raise :class:`RDKitUnavailableError`."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RDKitUnavailableError(
            "RDKit is required for molecular graphs and scaffold-based splits; "
            "install RDKit in this environment before using this operation"
        ) from exc
    return Chem


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize SMILES; returns the input unchanged when RDKit is absent."""
    try:
        from rdkit import Chem
    except ImportError:
        return smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse SMILES")
    return Chem.MolToSmiles(molecule, canonical=True)
