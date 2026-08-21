from __future__ import annotations

import re


_SUBSTITUTION = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])")


def apply_substitutions(reference_sequence: str, notation: str) -> str:
    sequence = "".join(reference_sequence.split()).upper()
    if not sequence or any(residue not in "ACDEFGHIKLMNPQRSTVWY" for residue in sequence):
        raise ValueError("reference sequence must contain only canonical amino acids")
    variant = notation.strip().upper()
    if variant == "WT":
        return sequence
    if not variant:
        raise ValueError("variant notation must not be empty")

    residues = list(sequence)
    changed_positions: set[int] = set()
    for substitution in variant.split("-"):
        match = _SUBSTITUTION.fullmatch(substitution)
        if not match:
            raise ValueError(f"invalid substitution notation: {substitution!r}")
        reference, position_text, alternate = match.groups()
        position = int(position_text)
        if not 1 <= position <= len(residues):
            raise ValueError(f"substitution position is outside the reference: {substitution!r}")
        if position in changed_positions:
            raise ValueError(f"duplicate substitution position: {position}")
        if residues[position - 1] != reference:
            raise ValueError(
                f"reference residue mismatch at {position}: expected {residues[position - 1]}, "
                f"found {reference}"
            )
        residues[position - 1] = alternate
        changed_positions.add(position)
    return "".join(residues)
