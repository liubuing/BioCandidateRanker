from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ._rdkit import canonicalize_smiles
from .schema import EnzymeSubstrateRecord, EvidenceTier


DATASET_DOI = "10.5061/dryad.7nd70"
ARTICLE_DOI = "10.1126/science.1115649"
DATASET_LICENSE = "CC0-1.0"
SOURCE_SHA256 = "7bfd71d78235ade06a9e117b8229ad99b1fe7adceb0a5ed72d708c3e925c2a69"
SOURCE_SIZE_BYTES = 38_761
UNIPROT_ACCESSION = "P30125"
MUTATION_POSITIONS = (236, 289, 290, 296, 337, 341)
REFERENCE_SEQUENCE = (
    "MSKNYHIAVLPGDGIGPEVMTQALKVLDAVRNRFAMRITTSHYDVGGAAIDNHGQPLPPATVEGCEQADAV"
    "LFGSVGGPKWEHLPPDQQPERGALLPLRKHFKLFSNLRPAKLYQGLEAFCPLRADIAANGFDILCVRELTGGI"
    "YFGQPKGREGSGQYEKAFDTEVYHRFEIERIARIAFESARKRRHKVTSIDKANVLQSSILWREIVNEIATEYP"
    "DVELAHMYIDNATMQLIKDPSQFDVLLCSNLFGDILSDECAMITGSMGMLPSASLNEQGFGLYEPAGGSAPDI"
    "AGKNIANPIAQILSLALLLRYSLDADDAACAIERAINRALEEGIRTGDLARGAAAVSTDEMGDIIARYVAEGV"
)
COFACTORS = {
    "NAD": (
        "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)"
        "OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N"
    ),
    "NADP": (
        "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)"
        "OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)"
        "C(=O)N"
    ),
}


@dataclass(frozen=True, slots=True)
class LunzerObservation:
    record: EnzymeSubstrateRecord
    genotype: str
    cofactor: str
    ln_km: float
    ln_efficiency: float
    derived_ln_kcat: float


def _variant_sequence(genotype: str) -> str:
    if len(genotype) != len(MUTATION_POSITIONS):
        raise ValueError(f"genotype must contain six residues: {genotype!r}")
    residues = list(REFERENCE_SEQUENCE)
    for position, residue in zip(MUTATION_POSITIONS, genotype):
        if residue not in "ACDEFGHIKLMNPQRSTVWY":
            raise ValueError(f"genotype contains a non-canonical residue: {genotype!r}")
        residues[position - 1] = residue
    return "".join(residues)


def read_lunzer_tsv(path: str | Path) -> tuple[LunzerObservation, ...]:
    source = Path(path)
    try:
        text = source.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read Lunzer TSV {source}: {exc}") from exc
    reader = csv.DictReader(io.StringIO(text, newline=None), delimiter="\t")
    expected = {
        "236", "289", "290", "296", "337", "341",
        "lnKmNAD", "lnKmNADP", "lnNAD", "lnNADP",
    }
    if not reader.fieldnames or not expected.issubset(reader.fieldnames):
        raise ValueError("Lunzer TSV does not contain the frozen kinetic columns")

    observations = []
    seen = set()
    for source_row, row in enumerate(reader, start=1):
        genotype = "".join(row[str(position)].strip().upper() for position in MUTATION_POSITIONS)
        if genotype in seen:
            raise ValueError(f"duplicate genotype at source row {source_row}: {genotype}")
        seen.add(genotype)
        sequence = _variant_sequence(genotype)
        for cofactor in COFACTORS:
            try:
                ln_km = float(row[f"lnKm{cofactor}"])
                ln_efficiency = float(row[f"ln{cofactor}"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"invalid {cofactor} kinetic value at source row {source_row}"
                ) from exc
            if not math.isfinite(ln_km) or not math.isfinite(ln_efficiency):
                raise ValueError(f"non-finite kinetic value at source row {source_row}")
            derived_ln_kcat = ln_km + ln_efficiency
            record = EnzymeSubstrateRecord(
                sequence=sequence,
                substrate_smiles=canonicalize_smiles(COFACTORS[cofactor]),
                substrate_name=cofactor,
                organism="Escherichia coli K-12",
                ec="1.1.1.85",
                enzyme_type="wildtype" if genotype == "DDIAGR" else "mutant",
                candidate_id=f"lunzer-{genotype}-{cofactor.lower()}",
                log10_kcat=derived_ln_kcat / math.log(10.0),
                evidence_tier=EvidenceTier.DIRECT,
                source_dataset="Lunzer et al. 2005 / Dryad 7nd70",
                source_row=source_row,
                campaign_group=cofactor,
            )
            observations.append(LunzerObservation(
                record=record,
                genotype=genotype,
                cofactor=cofactor,
                ln_km=ln_km,
                ln_efficiency=ln_efficiency,
                derived_ln_kcat=derived_ln_kcat,
            ))
    if len(seen) != 512 or len(observations) != 1024:
        raise ValueError(
            f"frozen Lunzer landscape must contain 512 genotypes and 1024 observations; "
            f"found {len(seen)} and {len(observations)}"
        )
    return tuple(observations)


def cofactor_counts(observations: Sequence[LunzerObservation]) -> dict[str, int]:
    return {
        cofactor: sum(item.cofactor == cofactor for item in observations)
        for cofactor in COFACTORS
    }
