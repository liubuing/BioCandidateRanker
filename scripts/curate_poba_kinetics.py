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
    / "biostudies-S-EPMC10770758"
)
RAW = SOURCE / "raw"
DOI = "10.1016/j.jbc.2023.105508"
SOURCE_ID = "biostudies-S-EPMC10770758"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"

RAW_FILES = (
    "PMC10770758-fullText.xml",
    "jbc.2023.105508.pdf",
    "mmc1.docx",
    "biostudies-metadata.json",
    "OR631772.1.gb",
    "OR631773.1.gb",
    "OR631774.1.gb",
    "WOL22124.1.fasta",
    "WP_028601822.1.fasta",
    "WP_130361895.1.fasta",
    "pubchem-reaction-compounds.json",
)

STRUCTURES = {
    "protocatechuate": (72, "C1=CC(=C(C=C1C(=O)O)O)O"),
    "4-hydroxybenzoate": (135, "C1=CC(=CC=C1C(=O)O)O"),
    "NADPH": (
        5884,
        "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O",
    ),
    "oxygen": (977, "O=O"),
    "gallic acid": (370, "C1=C(C=C(C(=C1O)O)O)C(=O)O"),
}

# Table, enzyme, sequence ID, substrate, kcat, SD, Km, SD, efficiency, SD, coupling.
FINITE_ROWS = (
    ("Table 1", "XaPobA WT", "poba-xa-wt", "protocatechuate", 1.8, 0.1, 0.37, 0.09, 4.9, 1.1, 53),
    ("Table 1", "XaPobA WT", "poba-xa-wt", "4-hydroxybenzoate", 0.047, 0.004, 0.43, 0.03, 0.11, 0.03, 18),
    ("Table 1", "OtPobA WT", "poba-ot-wt", "protocatechuate", 3.5, 0.3, 0.53, 0.19, 7.6, 2.8, 36),
    ("Table 1", "OtPobA WT", "poba-ot-wt", "4-hydroxybenzoate", 0.63, 0.19, 0.077, 0.026, 8.2, 0.4, 22),
    ("Table 1", "PkPobA WT", "poba-pk-wt", "protocatechuate", 1.1, 0.1, 0.32, 0.07, 3.7, 0.6, 36),
    ("Table 1", "PkPobA WT", "poba-pk-wt", "4-hydroxybenzoate", 0.11, 0.03, 0.37, 0.18, 0.34, 0.09, 9),
    ("Table 2", "XaPobA W201A", "poba-xa-w201a", "protocatechuate", 0.093, 0.011, 0.65, 0.27, 0.14, 0.04, 35),
    ("Table 2", "XaPobA W201Y", "poba-xa-w201y", "protocatechuate", 8.1, 0.6, 0.50, 0.11, 16.0, 3.0, 46),
    ("Table 2", "XaPobA W201Y", "poba-xa-w201y", "4-hydroxybenzoate", 0.79, 0.07, 0.072, 0.022, 11.0, 3.0, 42),
    ("Table 2", "XaPobA M210A", "poba-xa-m210a", "protocatechuate", 0.16, 0.01, 0.96, 0.03, 0.16, 0.04, 2),
    ("Table 2", "XaPobA C347T", "poba-xa-c347t", "protocatechuate", 0.083, 0.004, 0.32, 0.10, 0.29, 0.09, 19),
    ("Table 2", "XaPobA C347T", "poba-xa-c347t", "4-hydroxybenzoate", 0.086, 0.016, 0.91, 0.34, 0.10, 0.03, 0.4),
)

ACCEPTED_FIELDS = [
    "candidate_id",
    "article_doi",
    "source_table",
    "source_row",
    "organism",
    "sequence_id",
    "variable_substrate",
    "substrate_pubchem_cid",
    "substrate_isomeric_smiles",
    "kcat_s-1",
    "status_at_normalization",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def genbank_translation(path: Path) -> str:
    text = path.read_text(encoding="ascii")
    match = re.search(r'/translation="([A-Z\s]+)"', text)
    if not match:
        raise ValueError(f"No CDS translation in {path.name}")
    return re.sub(r"\s+", "", match.group(1))


def substitute(sequence: str, position: int, expected: str, replacement: str) -> str:
    if sequence[position - 1] != expected:
        raise ValueError(f"Expected {expected}{position}, found {sequence[position - 1]}{position}")
    return sequence[: position - 1] + replacement + sequence[position:]


def load_sequences() -> dict[str, tuple[str, str, str]]:
    xa = genbank_translation(RAW / "OR631772.1.gb")
    ot = genbank_translation(RAW / "OR631773.1.gb")
    pk = genbank_translation(RAW / "OR631774.1.gb")
    if (len(xa), len(ot), len(pk)) != (391, 391, 388):
        raise ValueError("Unexpected deposited OR631772-OR631774 translation lengths")
    if xa != read_fasta(RAW / "WOL22124.1.fasta"):
        raise ValueError("OR631772.1 and WOL22124.1 translations differ")
    if ot != read_fasta(RAW / "WP_028601822.1.fasta"):
        raise ValueError("OR631773.1 construct and WP_028601822.1 differ")
    native_pk = read_fasta(RAW / "WP_130361895.1.fasta")
    if pk[3:] != native_pk[6:] or len(native_pk) != 391:
        raise ValueError("OR631774.1 does not map to the WP_130361895.1 catalytic chain")
    return {
        "poba-xa-wt": ("XaPobA WT", "OR631772.1/WOL22124.1; A0A978C2P2 mapping reported", xa),
        "poba-ot-wt": ("OtPobA WT", "OR631773.1/WOL22125.1; WP_028601822.1 exact", ot),
        "poba-pk-wt": ("PkPobA WT", "OR631774.1/WOL22126.1; WP_130361895.1 source", pk),
        "poba-xa-w201a": ("XaPobA W201A", "OR631772.1 + W201A", substitute(xa, 201, "W", "A")),
        "poba-xa-w201y": ("XaPobA W201Y", "OR631772.1 + W201Y", substitute(xa, 201, "W", "Y")),
        "poba-xa-m210a": ("XaPobA M210A", "OR631772.1 + M210A", substitute(xa, 210, "M", "A")),
        "poba-xa-c347t": ("XaPobA C347T", "OR631772.1 + C347T", substitute(xa, 347, "C", "T")),
    }


def verify_raw_evidence() -> None:
    missing = [name for name in RAW_FILES if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw files: {', '.join(missing)}")
    if (RAW / "jbc.2023.105508.pdf").read_bytes()[:5] != b"%PDF-":
        raise ValueError("Article original is not a PDF")

    article = ET.parse(RAW / "PMC10770758-fullText.xml").getroot()
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected article DOI")
    article_text = " ".join(article.itertext())
    required_text = (
        "Initial velocity of NADPH consumption",
        "0.25 mM NADPH",
        "did not saturate up to 0.25 mM NAD(P)H",
        "Apparent steady-state kinetic parameters",
    )
    for text in required_text:
        if text not in article_text:
            raise ValueError(f"Missing assay evidence in source XML: {text}")

    metadata = json.loads((RAW / "biostudies-metadata.json").read_text(encoding="utf-8"))
    if metadata["accno"] != "S-EPMC10770758" or metadata["section"]["files"][0]["size"] != 29237:
        raise ValueError("Unexpected BioStudies metadata")
    with zipfile.ZipFile(RAW / "mmc1.docx") as archive:
        supplement_text = ET.fromstring(archive.read("word/document.xml"))
        if "Table S1" not in " ".join(supplement_text.itertext()):
            raise ValueError("Supporting information lacks Table S1")

    payload = json.loads((RAW / "pubchem-reaction-compounds.json").read_text(encoding="utf-8"))
    observed = {row["CID"]: row["SMILES"] for row in payload["PropertyTable"]["Properties"]}
    expected = {cid: smiles for cid, smiles in STRUCTURES.values()}
    if observed != expected:
        raise ValueError(f"PubChem structure mismatch: {observed}")


def build_exclusions() -> list[dict[str, object]]:
    organisms = {
        "XaPobA": "Xylophilus ampelinus CCH5-B3",
        "OtPobA": "Ottowia thiooxydans DSM 14619",
        "PkPobA": "Pigmentiphaga kullae DSM 13708",
    }
    rows = []
    for index, row in enumerate(FINITE_ROWS, 1):
        table, enzyme, sequence_id, substrate, kcat, kcat_sd, km, km_sd, efficiency, efficiency_sd, coupling = row
        cid, smiles = STRUCTURES[substrate]
        rows.append(
            {
                "exclusion_id": f"poba-hold-{index:03d}",
                "article_doi": DOI,
                "source_table": table,
                "source_row": f"{enzyme}; {substrate}",
                "organism": organisms[enzyme.split()[0]],
                "sequence_id": sequence_id,
                "enzyme_variant": enzyme,
                "variable_substrate": substrate,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "source_endpoint": "apparent_kcat_s-1",
                "source_kcat_s-1": kcat,
                "source_kcat_sd_s-1": kcat_sd,
                "source_km_mM": km,
                "source_km_sd_mM": km_sd,
                "source_kcat_per_km_mM-1_s-1": efficiency,
                "source_kcat_per_km_sd_mM-1_s-1": efficiency_sd,
                "coupling_ratio_percent": coupling,
                "assay_nadph_mM": 0.25,
                "nadph_saturation_demonstrated": False,
                "nadph_audit": "Rates did not saturate up to 0.25 mM NAD(P)H; no NADPH Km was determined.",
                "oxygen_pubchem_cid": STRUCTURES["oxygen"][0],
                "oxygen_isomeric_smiles": STRUCTURES["oxygen"][1],
                "oxygen_concentration_reported": False,
                "oxygen_saturation_demonstrated": False,
                "oxygen_audit": "Steady-state method reports no O2 concentration, equilibration, titration, or O2 Km.",
                "assay_buffer": "20 mM Tris-HCl",
                "assay_pH": 7.9,
                "assay_temperature_C": 25,
                "reaction_volume_uL": 100,
                "variable_substrate_range_mM": "0-2",
                "monitoring_wavelength_nm": 340,
                "nadph_extinction_coefficient_M-1_cm-1": 6300,
                "replicates": 3,
                "error_type": "SD",
                "fit_method": "Michaelis-Menten nonlinear regression",
                "exclusion_reason": "nadph_not_saturated_and_oxygen_saturation_not_demonstrated",
                "candidate_label_created": False,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs() -> None:
    verify_raw_evidence()
    sequences = load_sequences()
    exclusions = build_exclusions()

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (variant, accession, sequence) in sequences.items():
            handle.write(f">{sequence_id} {variant} | {accession} | deposited catalytic chain\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    write_csv(SOURCE / "candidate_records.csv", [], ACCEPTED_FIELDS)
    write_csv(SOURCE / "exclusions.csv", exclusions)
    raw_hashes = {
        name: {"size_bytes": (RAW / name).stat().st_size, "sha256": sha256(RAW / name)}
        for name in RAW_FILES
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )

    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC10770758",
        "article_doi": DOI,
        "article_published": "2023-11-28",
        "license": "CC-BY-4.0",
        "scope": "Novel clade-F XaPobA, OtPobA, PkPobA and finite XaPobA mutant steady-state rows",
        "reported_direct_finite_apparent_steady_state_rows": len(exclusions),
        "accepted_records": 0,
        "held_saturation_records": len(exclusions),
        "unique_deposited_wild_type_constructs": 3,
        "unique_audited_sequences_including_mutants": len(sequences),
        "kinetics_sources": [
            "raw/PMC10770758-fullText.xml, Tables 1-2 and Steady state kinetics of NAD(P)H oxidation",
            "raw/jbc.2023.105508.pdf, CC-BY article original",
            "raw/mmc1.docx, supporting Tables S1-S3 and Figures S1-S10",
        ],
        "sequence_sources": [
            "raw/OR631772.1.gb (XaPobA assayed codon-optimized CDS; WOL22124.1)",
            "raw/OR631773.1.gb (OtPobA assayed codon-optimized CDS; WOL22125.1)",
            "raw/OR631774.1.gb (PkPobA assayed codon-optimized CDS; WOL22126.1)",
            "raw/WP_028601822.1.fasta and raw/WP_130361895.1.fasta (native accession audit)",
        ],
        "construct_audit": "OR631772-OR631774 are the exact deposited codon-optimized assayed CDS translations. The paper reports modified-pET28a proteins with His6, but does not disclose enough vector-derived residues to reconstruct a complete tag context; FASTA does not invent those residues.",
        "assay": {
            "method": "A340 initial NADPH-consumption rates fitted by Michaelis-Menten nonlinear regression",
            "buffer": "20 mM Tris-HCl pH 7.9",
            "temperature_C": 25,
            "volume_uL": 100,
            "nadph_mM": 0.25,
            "substrate_range_mM": "0-2",
            "replicates": 3,
            "error_type": "SD",
        },
        "nadph_saturation_audit": "Failed. The article explicitly says initial velocities did not saturate through 0.25 mM NAD(P)H and no NAD(P)H Km or kcat was determined.",
        "oxygen_saturation_audit": "Failed. Molecular oxygen is a reaction cosubstrate, but the steady-state method reports no dissolved-O2 concentration, air-equilibration protocol, O2 titration, O2 Km, or saturation multiple. Ambient O2 was not inferred.",
        "selection_policy": "Accept direct finite steady-state kcat only when NADPH and O2 saturation are demonstrated. All finite apparent values are retained solely as exclusion evidence; none is a candidate label.",
        "excluded_censored_rows": [
            "Table 2 XaPobA V199A with PCA and 4-hydroxybenzoate: kcat <1e-3 s-1",
            "Table 2 XaPobA W201A with 4-hydroxybenzoate: kcat <1e-3 s-1",
            "Table 2 XaPobA M210A with 4-hydroxybenzoate: kcat <1e-3 s-1",
        ],
        "other_exclusions": "PaPobA clade-A control rows; pre-steady-state kred; specific activities; coupling-only, product, structural, phylogenetic, and detection-limit results.",
        "pubchem_structure_source": "raw/pubchem-reaction-compounds.json; PubChem PUG REST records for protocatechuate, 4-hydroxybenzoate, NADPH, oxygen, and gallic acid; retrieved 2026-07-22",
        "construct_sequences_sha256": sha256(fasta_path),
        "raw_file_hashes": "raw-file-hashes.json",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hit_rows = [line.split("\t") for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        cluster_rows = [
            line.split("\t")
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        ]
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
            "reported_finite_apparent_rows": len(exclusions),
            "saturation_exclusions": len(exclusions),
            "candidate_records": 0,
            "unique_sequences_audited": len(sequences),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({row[0] for row in hit_rows}),
            "exact_sequence_overlap": sum(float(row[2]) >= 0.999 for row in hit_rows),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len({row[0] for row in cluster_rows}),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": 0,
            "accepted_unique_sequences": 0,
            "accepted_unique_substrates": 0,
            "readiness_gate_passes": False,
            "claim_boundary": "Saturation-exclusion audit only; no candidate labels or model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )

    verify_outputs(exclusions, sequences)


def verify_outputs(
    exclusions: list[dict[str, object]], sequences: dict[str, tuple[str, str, str]]
) -> None:
    if len(exclusions) != 12 or len(sequences) != 7:
        raise ValueError("Unexpected PobA curation cardinality")
    if any(row["candidate_label_created"] for row in exclusions):
        raise ValueError("A saturation-failing value was promoted to a label")
    for row in exclusions:
        for field in ("source_kcat_s-1", "source_kcat_sd_s-1", "source_km_mM"):
            if not math.isfinite(float(row[field])) or float(row[field]) <= 0:
                raise ValueError(f"Invalid {field} in {row['exclusion_id']}")
        if row["nadph_saturation_demonstrated"] or row["oxygen_saturation_demonstrated"]:
            raise ValueError("Unexpected saturation pass")
    if (SOURCE / "candidate_records.csv").read_text(encoding="utf-8").count("\n") != 1:
        raise ValueError("Candidate records must contain a header and zero labels")
    required = [
        SOURCE / "candidate_records.csv",
        SOURCE / "exclusions.csv",
        SOURCE / "construct_sequences.fasta",
        SOURCE / "raw-file-hashes.json",
        SOURCE / "provenance.json",
    ]
    if not all(path.is_file() and path.stat().st_size for path in required):
        raise ValueError("A standard curation artifact is missing or empty")
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if hit_path.exists() and not (SOURCE / "homology-audit.json").is_file():
        raise ValueError("Homology output exists without its audit")


if __name__ == "__main__":
    write_outputs()
