"""Lazily guarded RDKit access for molecular graphs and scaffold operations.

RDKit is an optional dependency; every entry point funnels through
:func:`require_rdkit` so a missing installation surfaces as one consistent,
fail-closed error instead of an import-time crash.

中文：分子图与骨架操作的惰性守卫式 RDKit 访问层。RDKit 是可选依赖；所有入口
统一经过 require_rdkit，缺失安装时呈现为同一个一致的 fail-closed（拒绝式）
错误，而不是在 import 时崩溃。
"""
from __future__ import annotations


class RDKitUnavailableError(RuntimeError):
    """Raised when an operation needs RDKit but it is not installed.
    中文：当某项操作需要 RDKit 而环境中未安装时抛出。
    """


def require_rdkit():
    """Return the ``rdkit.Chem`` module, or raise :class:`RDKitUnavailableError`.
    中文：返回 rdkit.Chem 模块；RDKit 不可用时抛出 RDKitUnavailableError。
    """
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RDKitUnavailableError(
            "RDKit is required for molecular graphs and scaffold-based splits; "
            "install RDKit in this environment before using this operation"
        ) from exc
    return Chem


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize SMILES; returns the input unchanged when RDKit is absent.
    中文：规范化 SMILES（canonical SMILES）；RDKit 缺失时原样返回输入。
    """
    try:
        from rdkit import Chem
    except ImportError:
        return smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse SMILES")
    return Chem.MolToSmiles(molecule, canonical=True)
