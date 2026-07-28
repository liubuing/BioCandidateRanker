from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "europepmc-PMC12817470"
)
RAW = SOURCE / "raw"
TARGET = ROOT / "artifacts" / "external" / "enzengdb-v1" / "homology" / "unikp_reference.fasta"
DOI = "10.1002/pro.70436"
SOURCE_ID = "europepmc-PMC12817470"
TARGET_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

NRPINZ = (
    "MTRIVITGASGNYGRGVADALVAMGRAADLILITRKPEKLAERADQGCTVRQGDFDHPATLPQAMAGGDVLLLISGTRV"
    "GARVVQHKAAIDAAVAAGLRHIVYTSFIGIDDPANPAEVRHDHIETERLIRASGLAFTMLRDAHYADAMLLMAGPQVMQS"
    "GKWFANAGQGREAMVWRDDCIASAVAVLTTPGHENRIYNITGPELQTFAEVAAIMAEVTGCPVDYVDLDDDAQYALFDGL"
    "GIPRRPVDDQTVAGVPWNSDDMVTFGRAIREGFLEICTDDVERLTGRPARSTRAMVEANVAMLRAAAGR"
)

# Enzyme, organism, substrate, Km, Km SD, Vmax, Vmax SD, kcat, efficiency.
ROWS = (
    ("NrPinZ", "Novosphingobium rhizosphaerae sp. LY", "racemic medioresinol", 30.9, 8.7, 383.8, 31.5, 9.3, 3.00e5),
    ("NrPinZ", "Novosphingobium rhizosphaerae sp. LY", "racemic syringaresinol", 28.8, 6.4, 794.4, 53.0, 27.8, 9.65e5),
    ("NaPinZ", "Novosphingobium aromaticivorans F199", "racemic pinoresinol", 97.6, 15.8, 1842.0, 125.1, 64.5, 9.76e5),
    ("NaPinZ", "Novosphingobium aromaticivorans F199", "racemic medioresinol", 24.5, 7.9, 71.3, 7.6, 2.5, 1.02e5),
    ("NaPinZ", "Novosphingobium aromaticivorans F199", "racemic syringaresinol", 3.1, 0.9, 596.7, 28.8, 20.9, 6.80e6),
    ("SlPinZ", "Sphingobium lignivorans SYK-6", "racemic pinoresinol", 22.4, 2.5, 135.9, 4.3, 4.8, 2.24e5),
    ("SlPinZ", "Sphingobium lignivorans SYK-6", "racemic medioresinol", 13.3, 2.2, 237.2, 7.5, 8.3, 6.26e5),
    ("SlPinZ", "Sphingobium lignivorans SYK-6", "racemic syringaresinol", 11.2, 3.2, 173.95, 9.6, 6.1, 5.46e5),
)

CANDIDATE_FIELDS = [
    "candidate_id",
    "article_doi",
    "source_table",
    "source_row",
    "organism",
    "sequence_id",
    "variable_substrate",
    "substrate_pubchem_cid",
    "substrate_isomeric_smiles",
    "endpoint",
    "kcat_s-1",
    "status_at_normalization",
]

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def pdb_sequence(path: Path) -> str:
    residues: dict[int, str] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
            residues[int(line[22:26])] = AA3[line[17:20]]
    return "".join(residues[index] for index in sorted(residues))


def supplement_text() -> str:
    with zipfile.ZipFile(RAW / "Supplementary File 1.docx") as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    return " ".join(text.strip() for text in root.itertext() if text.strip())


def pubchem_pair(stem: str) -> tuple[str, str, str]:
    plus = json.loads((RAW / f"pubchem-{stem}.json").read_text(encoding="utf-8"))
    minus = json.loads((RAW / f"pubchem-minus-{stem}.json").read_text(encoding="utf-8"))
    plus_rows = plus["PropertyTable"]["Properties"]
    minus_rows = minus["PropertyTable"]["Properties"]
    plus_row = next(row for row in plus_rows if "(+" in row["Title"] or ", (+" in row["Title"])
    minus_row = next(row for row in minus_rows if "(-" in row["Title"])
    if plus_row["ConnectivitySMILES"] != minus_row["ConnectivitySMILES"]:
        raise ValueError(f"PubChem enantiomers have different connectivity for {stem}")
    return (
        f"{plus_row['CID']};{minus_row['CID']}",
        f"{plus_row['SMILES']};{minus_row['SMILES']}",
        plus_row["ConnectivitySMILES"],
    )


def validate_sources() -> None:
    required = (
        "PMC12817470-fullText.xml",
        "PRO-35-e70436-s001.zip",
        "Supplementary File 1.docx",
        "pinz_nadph_model_0.pdb",
        "NrPinZ-source-sequence.fasta",
        "biostudies-metadata.json",
        "pubchem-pinoresinol.json",
        "pubchem-minus-pinoresinol.json",
        "pubchem-medioresinol.json",
        "pubchem-minus-medioresinol.json",
        "pubchem-syringaresinol.json",
        "pubchem-minus-syringaresinol.json",
    )
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {', '.join(missing)}")
    if sha256(TARGET) != TARGET_SHA256:
        raise ValueError("Frozen UniKP development target identity mismatch")

    article = ET.parse(RAW / "PMC12817470-fullText.xml").getroot()
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected article DOI")
    article_text = " ".join(article.itertext())
    for phrase in (
        "This is an open access article under the terms",
        "creativecommons.org/licenses/by/4.0",
        "10 different substrate concentrations",
        "Reactions were carried out in triplicate",
        "NADPH (10 mM, 40",
    ):
        if phrase not in re.sub(r"\s+", " ", article_text):
            raise ValueError(f"Missing article evidence: {phrase}")

    text = re.sub(r"\s+", "", supplement_text())
    if NRPINZ not in text:
        raise ValueError("Exact NrPinZ DATA S1 sequence is absent")
    if read_fasta(RAW / "NrPinZ-source-sequence.fasta") != NRPINZ:
        raise ValueError("Extracted NrPinZ FASTA differs from DATA S1")
    if pdb_sequence(RAW / "pinz_nadph_model_0.pdb") != NRPINZ:
        raise ValueError("NrPinZ Boltz coordinate sequence differs from DATA S1")


def build_exclusions() -> list[dict[str, object]]:
    structures = {
        name: pubchem_pair(name)
        for name in ("pinoresinol", "medioresinol", "syringaresinol")
    }
    rows = []
    for index, row in enumerate(ROWS, 1):
        enzyme, organism, substrate, km, km_sd, vmax, vmax_sd, kcat, efficiency = row
        stem = substrate.removeprefix("racemic ")
        cids, smiles, connectivity = structures[stem]
        sequence_status = (
            "exact_DATA_S1_and_Boltz_coordinate_sequence"
            if enzyme == "NrPinZ"
            else "blocked_no_sequence_or_accession_in_article_or_supplement"
        )
        rows.append(
            {
                "exclusion_id": f"pinz-hold-{index:03d}",
                "article_doi": DOI,
                "source_table": "Table 1",
                "source_row": f"{enzyme}; {substrate}",
                "primary_measurement": True,
                "organism": organism,
                "enzyme": enzyme,
                "sequence_id": "pinz-nr" if enzyme == "NrPinZ" else "",
                "sequence_mapping_status": sequence_status,
                "variable_substrate": substrate,
                "substrate_pubchem_cids_enantiomer_pair": cids,
                "substrate_isomeric_smiles_enantiomer_pair": smiles,
                "substrate_connectivity_smiles": connectivity,
                "source_endpoint": "apparent_kcat_s-1",
                "source_km_uM": km,
                "source_km_sd_uM": km_sd,
                "source_vmax_pkat_per_ug": vmax,
                "source_vmax_sd_pkat_per_ug": vmax_sd,
                "source_kcat_s-1": kcat,
                "source_kcat_per_km_M-1_s-1": efficiency,
                "efficiency_transcription_note": (
                    "Table XML displays '1.02 +/- 10^5'; narrative states 1.02 x 10^5, retained."
                    if enzyme == "NaPinZ" and stem == "medioresinol"
                    else ""
                ),
                "assay_buffer": "20 mM Tris-HCl",
                "assay_pH": 7.0,
                "assay_temperature_C": 30,
                "substrate_range_uM": "0.5-400; 10 concentrations",
                "nadph_concentration_mM": 1.6,
                "nadph_matching_km": "not reported",
                "nadph_saturation_demonstrated": False,
                "replicates": 3,
                "fit_method": "Origin Hyperbola Michaelis-Menten fit",
                "exclusion_reason": "fixed_NADPH_saturation_not_demonstrated",
                "candidate_label_created": False,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def raw_hashes() -> dict[str, dict[str, object]]:
    return {
        path.relative_to(RAW).as_posix(): {
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(RAW.rglob("*"))
        if path.is_file()
    }


def write_outputs() -> None:
    validate_sources()
    exclusions = build_exclusions()

    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(">pinz-nr NrPinZ | exact DATA S1 catalytic sequence; N-terminal His6 context unreported\n")
        for start in range(0, len(NRPINZ), 80):
            handle.write(NRPINZ[start : start + 80] + "\n")
    write_csv(SOURCE / "candidate_records.csv", [], CANDIDATE_FIELDS)
    write_csv(SOURCE / "exclusions.csv", exclusions)
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes(), indent=2) + "\n", encoding="ascii"
    )
    license_path = SOURCE / "SOURCE-LICENSE.txt"
    license_path.write_text(
        "Smith CA et al. Molecular basis and biological relevance of bacterial and plant "
        "pinoresinol/lariciresinol reductase specificities. Protein Science (2026).\n"
        "DOI: 10.1002/pro.70436\n"
        "Copyright (c) 2026 The Authors. CC BY 4.0.\n"
        "License: https://creativecommons.org/licenses/by/4.0/\n"
        "Source: https://europepmc.org/articles/PMC12817470\n",
        encoding="ascii",
    )
    commands_path = SOURCE / "homology-commands.txt"
    commands_path.write_text(
        "WSL distribution: Ubuntu-24.04\n"
        "mmseqs easy-search construct_sequences.fasta "
        "unikp_reference.fasta homology/homology_hits.tsv TMP "
        "--min-seq-id 0.3 -c 0.8 --cov-mode 0 "
        "--format-output query,target,pident,qcov,tcov,evalue,bits\n"
        "mmseqs easy-linclust construct_sequences.fasta "
        "homology/family-cluster/proteins TMP --min-seq-id 0.3 -c 0.8 --cov-mode 0\n",
        encoding="ascii",
    )

    blocker = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "status": "scientifically_ineligible_absolute_kcat",
        "blocking_gates": [
            {
                "gate": "absolute_turnover_cosubstrate_saturation",
                "affected_rows": 8,
                "missing_information": "A matching NADPH Km, NADPH titration showing saturation, or another quantitative saturation multiple for the 1.6 mM fixed NADPH condition.",
                "resolution_required": "Primary assay evidence demonstrating that fixed NADPH was saturating for NrPinZ, NaPinZ, and SlPinZ.",
            },
            {
                "gate": "unambiguous_sequence_mapping",
                "affected_enzymes": ["NaPinZ", "SlPinZ"],
                "affected_rows": 6,
                "missing_information": "Exact assayed amino-acid sequences or unambiguous accession/locus identifiers for the synthesized full-length genes.",
                "resolution_required": "Author-deposited construct sequences or accession/locus mappings; nearest-proteome homologs are not accepted as exact constructs.",
            },
        ],
        "candidate_records_emitted": 0,
        "model_predictions_run": False,
    }
    (SOURCE / "blocker.json").write_text(
        json.dumps(blocker, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.exists() or not cluster_path.exists():
        raise FileNotFoundError("Run the frozen MMseqs homology commands before curation")
    hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    clusters = [line for line in cluster_path.read_text(encoding="utf-8").splitlines() if line]

    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC12817470",
        "biostudies_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12817470",
        "article_doi": DOI,
        "article_published": "2026-01-20",
        "license": "CC-BY-4.0",
        "license_scope": "Article, Table 1 measurement data, and supporting files",
        "reported_new_direct_finite_apparent_kcat_rows": len(exclusions),
        "saturation_qualified_records": 0,
        "candidate_records": 0,
        "kinetics_source": "raw/PMC12817470-fullText.xml, Table 1 rows marked This study and Sections 4.2-4.4",
        "sequence_source": "raw/Supplementary File 1.docx DATA S1; independently identical to raw/pinz_nadph_model_0.pdb",
        "sequence_mapping": "NrPinZ is exact at 308 aa. The article and DATA S1 omit NaPinZ and SlPinZ sequences/accessions; strain-proteome nearest homologs Q2G4H9 (79.3% identity) and WP_014075058.1 (74.6% identity) are audit leads only and were not substituted for assayed constructs.",
        "construct_caveat": "All proteins had in-frame N-terminal His6 tags, but complete vector-derived residues were not reported; FASTA does not invent them.",
        "assay": {
            "method": "Reverse-phase HPLC initial-rate substrate series; Origin Hyperbola Michaelis-Menten fit",
            "buffer": "20 mM Tris-HCl pH 7.0",
            "temperature_C": 30,
            "variable_substrate_range_uM": "0.5-400; 10 concentrations",
            "fixed_nadph_mM": 1.6,
            "enzyme_ng": "6-10",
            "replicates": 3,
        },
        "saturation_policy": "A fixed reaction cosubstrate requires a matching Km or quantitative saturation evidence. No NADPH Km/titration is reported, so all apparent kcat rows are exclusions.",
        "excluded_prior_measurements": "The two NrPinZ pinoresinol rows, plant homolog rows, and native Forsythia rows cite earlier articles and are not primary measurements of this 2026 source.",
        "raw_file_hashes": "raw-file-hashes.json",
        "source_license_artifact": "SOURCE-LICENSE.txt",
        "source_license_sha256": sha256(license_path),
        "construct_sequences_sha256": sha256(fasta),
        "blocker_artifact": "blocker.json",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    audit = {
        "audited_on": "2026-07-27",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust",
        "wsl_distribution": "Ubuntu-24.04",
        "mmseqs_binary": "mmseqs",
        "mmseqs_version": MMSEQS_VERSION,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": sha256(TARGET),
        "reported_new_direct_finite_apparent_kcat_rows": len(exclusions),
        "saturation_exclusions": len(exclusions),
        "candidate_records": 0,
        "unique_exact_source_sequences_audited": 1,
        "construct_sequences_sha256": sha256(fasta),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
        "homology_hits_sha256": sha256(hit_path),
        "commands_artifact": "homology-commands.txt",
        "commands_sha256": sha256(commands_path),
        "candidate_mmseqs_families": len({line.split("\t", 1)[0] for line in clusters}),
        "family_cluster_sha256": sha256(cluster_path),
        "accepted_records": 0,
        "readiness_gate_passes": False,
        "claim_boundary": "Saturation-exclusion and sequence-availability audit only; no candidate labels or model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )
    verify_outputs(exclusions)


def verify_outputs(exclusions: list[dict[str, object]]) -> None:
    if len(exclusions) != 8 or len(NRPINZ) != 308:
        raise ValueError("Unexpected PinZ curation cardinality")
    for row in exclusions:
        for field in ("source_km_uM", "source_vmax_pkat_per_ug", "source_kcat_s-1"):
            if not math.isfinite(float(row[field])) or float(row[field]) <= 0:
                raise ValueError(f"Invalid {field} in {row['exclusion_id']}")
        if row["nadph_saturation_demonstrated"] or row["candidate_label_created"]:
            raise ValueError("A saturation-failing row was promoted to a label")
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        if list(csv.DictReader(handle)):
            raise ValueError("Candidate CSV must contain no ineligible labels")
    required = (
        "candidate_records.csv", "exclusions.csv", "construct_sequences.fasta",
        "raw-file-hashes.json", "SOURCE-LICENSE.txt", "homology-commands.txt",
        "provenance.json", "blocker.json", "homology-audit.json",
    )
    if not all((SOURCE / name).is_file() and (SOURCE / name).stat().st_size for name in required):
        raise ValueError("A required curation artifact is missing or empty")


if __name__ == "__main__":
    write_outputs()
