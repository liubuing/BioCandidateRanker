from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "biostudies-S-EPMC12412097"
)
RAW = SOURCE / "raw"
DOI = "10.1021/jacs.5c05573"
SOURCE_ID = "biostudies-S-EPMC12412097"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

RAW_FILES = (
    "PMC12412097-fullText.xml",
    "ja5c05573_si_001.pdf",
    "ja5c05573_0005.jpg",
    "ja5c05573_0006.jpg",
    "9J47.cif",
    "Q2FK11.fasta",
    "ABD21118.1.fasta",
    "pubchem-nitrite.json",
)

# Figure, variant, sequence ID, mutations, oligomeric state, Km, Vmax, kcat, efficiency.
ROWS = (
    ("Figure 5C", "WT", "scda-wt", (), "partially oligomeric", 681, 32.2, 32.3, 4.75e4),
    ("Figure 5C", "S77E", "scda-s77e", ((77, "S", "E"),), "partially oligomeric", 1035, 62.9, 61.8, 5.98e4),
    ("Figure 6A", "CF", "scda-cf", ((30, "C", "A"), (31, "C", "A"), (191, "C", "A")), "predominantly dimeric", 822, 94.4, 94.3, 11.47e4),
    ("Figure 6A", "CF-S77E", "scda-cf-s77e", ((30, "C", "A"), (31, "C", "A"), (77, "S", "E"), (191, "C", "A")), "predominantly monomeric", 1468, 118, 117.9, 8.03e4),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def substitute(sequence: str, mutations: tuple[tuple[int, str, str], ...]) -> str:
    result = sequence
    for position, expected, replacement in mutations:
        if result[position - 1] != expected:
            raise ValueError(
                f"Expected {expected}{position}, found {result[position - 1]}{position}"
            )
        result = result[: position - 1] + replacement + result[position:]
    return result


def build_sequences() -> dict[str, tuple[str, str]]:
    wt = read_fasta(RAW / "Q2FK11.fasta")
    if wt != read_fasta(RAW / "ABD21118.1.fasta"):
        raise ValueError("UniProt Q2FK11 and GenBank ABD21118.1 sequences differ")
    if len(wt) != 224:
        raise ValueError(f"Expected 224-aa ScdA, found {len(wt)} aa")

    sequences = {
        sequence_id: (variant, substitute(wt, mutations))
        for _, variant, sequence_id, mutations, *_ in ROWS
    }
    cif = (RAW / "9J47.cif").read_text(encoding="utf-8")
    match = re.search(
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s+;([^;]+);", cif, re.DOTALL
    )
    if not match:
        raise ValueError("Could not recover canonical polymer sequence from 9J47.cif")
    pdb_sequence = re.sub(r"\s+", "", match.group(1))
    if pdb_sequence != sequences["scda-cf"][1]:
        raise ValueError("PDB 9J47 sequence does not match the C30A/C31A/C191A CF variant")
    return sequences


def read_substrate() -> tuple[int, str]:
    payload = json.loads((RAW / "pubchem-nitrite.json").read_text(encoding="utf-8"))
    compound = payload["PropertyTable"]["Properties"][0]
    if compound["CID"] != 946 or compound["SMILES"] != "N(=O)[O-]":
        raise ValueError("Unexpected PubChem nitrite record")
    return compound["CID"], compound["SMILES"]


def build_records() -> list[dict[str, object]]:
    cid, smiles = read_substrate()
    records = []
    for index, row in enumerate(ROWS, 1):
        figure, variant, sequence_id, mutations, oligomeric_state, km, vmax, kcat, efficiency = row
        mutation_text = "/".join(f"{old}{position}{new}" for position, old, new in mutations)
        records.append(
            {
                "candidate_id": f"scda-{index:03d}",
                "article_doi": DOI,
                "source_figure": figure,
                "source_row": f"ScdA {variant}; nitrite",
                "organism": "Staphylococcus aureus subsp. aureus USA300_FPR3757",
                "sequence_accession": "Q2FK11; ABD21118.1",
                "sequence_id": sequence_id,
                "enzyme_variant": variant,
                "amino_acid_changes": mutation_text or "WT",
                "construct": "full-length ScdA with N-terminal His6 tag; vector-derived tag sequence not reported",
                "oligomeric_state_in_mv_assay": oligomeric_state,
                "redox_ligand_state": "as-isolated di-iron protein with reduced methyl viologen electron donor; no DTT",
                "variable_substrate": "nitrite",
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "endpoint": "kcat_s-1",
                "source_kcat_min-1": kcat,
                "source_kcat_uncertainty_percent": 5,
                "kcat_s-1": kcat / 60,
                "kcat_error_s-1": kcat * 0.05 / 60,
                "km_uM": km,
                "km_error_uM": km * 0.05,
                "vmax_uM_min-1": vmax,
                "vmax_error_uM_min-1": vmax * 0.05,
                "kcat_per_km_M-1_min-1": efficiency,
                "kcat_per_km_M-1_s-1": efficiency / 60,
                "assay_pH": 7.6,
                "assay_temperature_C": "not reported",
                "enzyme_concentration_uM": 1,
                "reduced_methyl_viologen_uM": 116,
                "nitrite_range_uM": "1-1000" if variant in {"WT", "CF"} else "1-1500",
                "monitoring_wavelength_nm": 600,
                "methyl_viologen_extinction_mM-1_cm-1": 13.70,
                "maximum_monitoring_time_s": 1200,
                "replicates": ">=3 independent measurements",
                "error_type": "estimated uncertainty from propagated independent measurements",
                "fit_method": "Michaelis-Menten nonlinear regression; confirmed by Lineweaver-Burk transformation",
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def apply_homology_status(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if not hit_path.exists():
        return
    hit_ids = {
        line.split("\t", 1)[0]
        for line in hit_path.read_text(encoding="utf-8").splitlines()
        if line
    }
    for record in records:
        record["status_at_normalization"] = (
            "excluded_homology"
            if record["sequence_id"] in hit_ids
            else "accepted_homology_cold_pool"
        )


def write_outputs() -> None:
    missing = [name for name in RAW_FILES if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw source files: {', '.join(missing)}")

    sequences = build_sequences()
    records = build_records()
    apply_homology_status(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (variant, sequence) in sequences.items():
            handle.write(
                f">{sequence_id} Q2FK11 ABD21118.1 | ScdA {variant} catalytic chain; "
                "unreported N-terminal His6-tag context excluded\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    raw_hashes = {name: sha256(RAW / name) for name in RAW_FILES}
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12412097",
        "article_doi": DOI,
        "article_published": "2025-08-22",
        "license": "CC-BY-4.0",
        "reported_finite_steady_state_kcat_records": len(records),
        "accepted_records": len(accepted),
        "unique_construct_sequences": len(sequences),
        "unique_substrates": 1,
        "kinetics_sources": [
            "raw/PMC12412097-fullText.xml, Figures 5-6 and Nitrite Reductase Assay methods",
            "raw/ja5c05573_0005.jpg, Figure 5 numerical table",
            "raw/ja5c05573_0006.jpg, Figure 6 numerical table",
            "raw/ja5c05573_si_001.pdf, Figures S5-S7 controls and state validation",
        ],
        "sequence_sources": [
            "raw/Q2FK11.fasta",
            "raw/ABD21118.1.fasta",
            "raw/9J47.cif (exact CF variant sequence)",
        ],
        "sequence_mapping": "Q2FK11 and ABD21118.1 are identical 224-aa USA300 ScdA sequences. PDB 9J47 is tagless C30A/C31A/C191A CF ScdA and exactly matches the generated CF catalytic chain.",
        "construct_caveat": "Kinetic proteins were full-length N-terminal His6-tagged constructs, but the complete vector-derived tag/TEV sequence was not reported. FASTA therefore contains exact accession-mapped catalytic chains and does not invent tag residues.",
        "variant_state_policy": "WT, S77E, CF (C30A/C31A/C191A), and CF-S77E are enzyme variants. Partially oligomeric, predominantly dimeric, predominantly monomeric, as-isolated di-iron, DTT-reduced, and iron-nitrosyl are physical/redox/ligand states and are not represented as extra sequences.",
        "assay": {
            "method": "Anaerobic reduced-methyl-viologen consumption; early A600 initial rates fitted to Michaelis-Menten and checked by Lineweaver-Burk",
            "pH": 7.6,
            "buffer_composition": "not explicitly reported for the reaction buffer",
            "temperature_C": None,
            "enzyme_concentration_uM": 1,
            "reduced_methyl_viologen_uM": 116,
            "nitrite_concentration_uM": "1-1500 across figure panels",
            "monitoring_wavelength_nm": 600,
            "methyl_viologen_extinction_mM-1_cm-1": 13.70,
            "maximum_monitoring_time_s": 1200,
            "replicates": ">=3 independent measurements",
            "uncertainty": "+/-5% for reported kinetic parameters, estimated by error propagation",
        },
        "normalization": "kcat_s-1 = author-reported kcat_min-1 / 60; uncertainty converted by the same factor",
        "excluded_data": [
            "Figure 6B fixed-100-uM-nitrite specific activities: not steady-state absolute kcat",
            "H87L, H132L, E136L, H167L, H211L, and E215L: activity described as abolished with no direct finite kcat",
            "NTD-only and CTD-only fragments: relative/specific activity only, no direct finite kcat",
            "UV-vis/EPR redox and nitrosyl-state observations: mechanistic state data, not steady-state kcat",
            "In vivo growth and ATP viability measurements: non-kinetic endpoints",
        ],
        "raw_sha256": raw_hashes,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        audit = {
            "audited_on": "2026-07-22",
            "source_id": SOURCE_ID,
            "article_doi": DOI,
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": len(records),
            "unique_sequences": len(sequences),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}),
            "accepted_unique_substrates": len({record["variable_substrate"] for record in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )


if __name__ == "__main__":
    write_outputs()
