from __future__ import annotations

import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POOL = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics"
EDMOND = POOL / "edmond-S0HJ48"
FRACTAL = POOL / "europepmc-PMC11041685"
SEQUENCE_DOCX = EDMOND / "raw" / "41467_2024_54408_MOESM3_ESM.docx"

SUBSTRATES = {
    "acetyl-CoA": {
        "cid": 444493,
        "smiles": "CC(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O",
    },
    "oxaloacetate": {"cid": 970, "smiles": "C(C(=O)C(=O)O)C(=O)O"},
}

CODONS = {
    codon: amino_acid
    for amino_acid, codons in {
        "F": "TTT TTC",
        "L": "TTA TTG CTT CTC CTA CTG",
        "I": "ATT ATC ATA",
        "M": "ATG",
        "V": "GTT GTC GTA GTG",
        "S": "TCT TCC TCA TCG AGT AGC",
        "P": "CCT CCC CCA CCG",
        "T": "ACT ACC ACA ACG",
        "A": "GCT GCC GCA GCG",
        "Y": "TAT TAC",
        "*": "TAA TAG TGA",
        "H": "CAT CAC",
        "Q": "CAA CAG",
        "N": "AAT AAC",
        "K": "AAA AAG",
        "D": "GAT GAC",
        "E": "GAA GAG",
        "C": "TGT TGC",
        "W": "TGG",
        "R": "CGT CGC CGA CGG AGA AGG",
        "G": "GGT GGC GGA GGG",
    }.items()
    for codon in codons.split()
}


def docx_paragraphs(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    return [
        "".join(node.text or "" for node in paragraph.iter(namespace + "t")).strip()
        for paragraph in root.iter(namespace + "p")
    ]


def extract_dna(paragraphs: list[str], label: str) -> str:
    index = paragraphs.index(label) + 1
    while index < len(paragraphs):
        dna = paragraphs[index].replace(" ", "").upper()
        if dna and re.fullmatch(r"[ACGT]+", dna):
            if len(dna) % 3:
                raise ValueError(f"{label}: coding sequence length is not divisible by three")
            return dna
        index += 1
    raise ValueError(f"{label}: DNA sequence not found")


def translate(dna: str) -> str:
    protein = "".join(CODONS[dna[index : index + 3]] for index in range(0, len(dna), 3))
    if "*" in protein[:-1]:
        raise ValueError("Internal stop codon in construct sequence")
    return protein.removesuffix("*")


def substitute(sequence: str, changes: list[tuple[int, str, str]]) -> str:
    residues = list(sequence)
    for position, expected, replacement in changes:
        if residues[position - 1] != expected:
            raise ValueError(
                f"Expected {expected} at position {position}, found {residues[position - 1]}"
            )
        residues[position - 1] = replacement
    return "".join(residues)


def build_sequences() -> dict[str, dict[str, str]]:
    paragraphs = docx_paragraphs(SEQUENCE_DOCX)
    source_labels = {
        "cs-anc3a-wt": ("Anc 3a", "ancestral reconstruction"),
        "cs-anc3b-wt": ("Anc 3b", "ancestral reconstruction"),
        "cs-dp-wt": ("Deinococcus pimensis", "WP_027481744.1"),
        "cs-ac-wt": ("Ananas comosus", "XP_020089322.1"),
        "cs-ms-wt": ("Methylophaga sulfidovorans", "WP_091715369.1"),
        "cs-secs-wt": ("Synechococcus elongatus", "ABB56644.1"),
    }
    sequences = {
        sequence_id: {
            "sequence": translate(extract_dna(paragraphs, label)),
            "source_label": label,
            "accession": accession,
        }
        for sequence_id, (label, accession) in source_labels.items()
    }
    sequences["cs-anc3a-plus5"] = {
        "sequence": substitute(
            sequences["cs-anc3a-wt"]["sequence"],
            [(76, "A", "R"), (78, "P", "K"), (79, "E", "Y"), (90, "K", "E"), (145, "Q", "D")],
        ),
        "source_label": "Anc 3a + a76R/p78K/e79Y/k90E/q145D",
        "accession": "ancestral reconstruction plus five substitutions",
    }
    sequences["cs-ms-w150a"] = {
        "sequence": substitute(sequences["cs-ms-wt"]["sequence"], [(150, "W", "A")]),
        "source_label": "Methylophaga sulfidovorans W150A",
        "accession": "WP_091715369.1 W150A",
    }
    sequences["cs-ac-delta487-513"] = {
        "sequence": sequences["cs-ac-wt"]["sequence"][:486],
        "source_label": "Ananas comosus delta487-513",
        "accession": "XP_020089322.1 delta487-513",
    }
    sequences["cs-secs-l18q"] = {
        "sequence": substitute(sequences["cs-secs-wt"]["sequence"], [(18, "L", "Q")]),
        "source_label": "Synechococcus elongatus L18Q",
        "accession": "ABB56644.1 L18Q",
    }
    return sequences


def kinetic_record(
    candidate_id: str,
    article_doi: str,
    source_table: str,
    source_row: str,
    organism: str,
    sequence_id: str,
    construct: str,
    variable_substrate: str,
    kcat: float,
    kcat_sem: float,
    km_uM: float,
    km_sem_uM: float,
    note: str = "",
) -> dict[str, object]:
    fixed_substrate = "oxaloacetate" if variable_substrate == "acetyl-CoA" else "acetyl-CoA"
    variable = SUBSTRATES[variable_substrate]
    fixed = SUBSTRATES[fixed_substrate]
    return {
        "candidate_id": candidate_id,
        "article_doi": article_doi,
        "source_table": source_table,
        "source_row": source_row,
        "organism": organism,
        "sequence_id": sequence_id,
        "construct": construct,
        "variable_substrate": variable_substrate,
        "variable_substrate_pubchem_cid": variable["cid"],
        "variable_substrate_isomeric_smiles": variable["smiles"],
        "fixed_substrate": fixed_substrate,
        "fixed_substrate_pubchem_cid": fixed["cid"],
        "fixed_substrate_concentration_uM": 1000
        if fixed_substrate == "oxaloacetate"
        else 500,
        "endpoint": "kcat_s-1",
        "kcat_s-1": kcat,
        "kcat_sem_s-1": kcat_sem,
        "km_uM": km_uM,
        "km_sem_uM": km_sem_uM,
        "assay_pH": 7.5,
        "assay_temperature_C": 25,
        "enzyme_concentration_nM": 25,
        "error_type": "SEM",
        "source_anomaly": note,
        "status_at_normalization": "pending_homology",
    }


def edmond_records() -> list[dict[str, object]]:
    records = []
    dp_rows = [
        ("non-fractionated", 36.6, 1.1, 25.4, 3.1),
        ("F1", 39.9, 1.3, 29.6, 3.9),
        ("F2", 23.2, 0.6, 15.0, 1.8),
        ("F3", 16.6, 0.5, 17.7, 2.5),
        ("F4", 18.4, 0.6, 18.4, 2.5),
        ("F5", 14.2, 0.7, 25.9, 5.0),
    ]
    for index, (fraction, kcat, kcat_sem, km, km_sem) in enumerate(dp_rows, 1):
        records.append(
            kinetic_record(
                f"cs-edmond-{index:03d}",
                "10.1038/s41467-024-54408-6",
                "Supplementary Table S3",
                f"D. pimensis - {fraction}",
                "Deinococcus pimensis",
                "cs-dp-wt",
                f"WT; SEC fraction {fraction}",
                "oxaloacetate",
                kcat,
                kcat_sem,
                km,
                km_sem,
            )
        )

    rows = [
        ("anc3a", "ancestral reconstruction", "cs-anc3a-wt", (9.1, 0.3, 22.4, 2.6), (10.3, 0.1, 23.9, 1.1), ""),
        ("anc3b", "ancestral reconstruction", "cs-anc3b-wt", (5.6, 0.2, 8.6, 1.9), (7.5, 0.2, 13.2, 2.2), ""),
        ("anc3a+5", "ancestral reconstruction", "cs-anc3a-plus5", (21.9, 0.8, 56.5, 6.8), (23.3, 61.1, 61.1, 10.2), "The source PDF reports oxaloacetate kcat as 23.3 +/- 61.1 s-1; the unusually large SEM is retained without correction."),
        ("M. sulfidovorans WT", "Methylophaga sulfidovorans", "cs-ms-wt", (42.4, 1.0, 127.4, 7.7), (33.3, 1.8, 50.5, 5.2), ""),
        ("M. sulfidovorans W150A", "Methylophaga sulfidovorans", "cs-ms-w150a", (34.1, 1.0, 104.7, 8.9), (26.7, 1.1, 36.7, 5.7), ""),
        ("A. comosus WT", "Ananas comosus", "cs-ac-wt", (41.1, 1.5, 13.2, 2.3), (49.2, 2.2, 72.1, 10.1), ""),
        ("A. comosus delta487-513", "Ananas comosus", "cs-ac-delta487-513", (41.8, 1.5, 11.3, 2.0), (58.9, 2.6, 87.2, 11.7), ""),
    ]
    for label, organism, sequence_id, acetyl, oxaloacetate, anomaly in rows:
        for substrate, values in (("acetyl-CoA", acetyl), ("oxaloacetate", oxaloacetate)):
            records.append(
                kinetic_record(
                    f"cs-edmond-{len(records) + 1:03d}",
                    "10.1038/s41467-024-54408-6",
                    "Supplementary Table 4",
                    f"{label}; {substrate} columns",
                    organism,
                    sequence_id,
                    label,
                    substrate,
                    *values,
                    note=anomaly if substrate == "oxaloacetate" else "",
                )
            )
    return records


def fractal_records() -> list[dict[str, object]]:
    rows = [
        ("SeCS WT", "cs-secs-wt", (24.2, 1.3, 91.8, 9.2), (19.5, 1.7, 65.6, 10.6)),
        ("SeCS L18Q", "cs-secs-l18q", (21.7, 1.0, 78.5, 6.2), (20.5, 1.5, 57.2, 7.8)),
    ]
    records = []
    for label, sequence_id, acetyl, oxaloacetate in rows:
        for substrate, values in (("acetyl-CoA", acetyl), ("oxaloacetate", oxaloacetate)):
            records.append(
                kinetic_record(
                    f"cs-fractal-{len(records) + 1:03d}",
                    "10.1038/s41586-024-07287-2",
                    "Supplementary Table 1",
                    f"{label}; {substrate} columns",
                    "Synechococcus elongatus PCC 7942",
                    sequence_id,
                    label,
                    substrate,
                    *values,
                )
            )
    return records


def write_records(path: Path, records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_fasta(path: Path, sequence_ids: list[str], sequences: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id in sequence_ids:
            entry = sequences[sequence_id]
            handle.write(f">{sequence_id} {entry['accession']} | {entry['source_label']}\n")
            sequence = entry["sequence"]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def write_provenance(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def apply_homology_status(records: list[dict[str, object]], hit_path: Path) -> set[str]:
    if not hit_path.exists():
        return set()
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
    return hit_ids


def write_homology_audit(
    output_dir: Path,
    records: list[dict[str, object]],
    sequence_path: Path,
    article_doi: str,
    dataset_reference: str,
) -> None:
    hit_path = output_dir / "homology" / "homology_hits.tsv"
    cluster_path = output_dir / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.exists() or not cluster_path.exists():
        return
    hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    clusters = {
        line.split("\t", 1)[0]
        for line in cluster_path.read_text(encoding="utf-8").splitlines()
        if line
    }
    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    payload = {
        "audited_on": "2026-07-21",
        "source_id": output_dir.name,
        "article_doi": article_doi,
        "dataset_reference": dataset_reference,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target": {
            "path": "artifacts/external/absolute-kinetics-screen/dryad-4964723/homology/unikp_reference.fasta",
            "sha256": "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9",
        },
        "candidate_records": len(records),
        "unique_sequences": len({record["sequence_id"] for record in records}),
        "construct_sequences_sha256": sha256(sequence_path),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
        "homology_hits_sha256": sha256(hit_path),
        "candidate_mmseqs_families": len(clusters),
        "family_cluster_sha256": sha256(cluster_path),
        "exact_sequence_overlap": 0,
        "accepted_records": len(accepted),
        "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}),
        "accepted_unique_substrates": len(
            {record["variable_substrate"] for record in accepted}
        ),
        "readiness_gate_passes": False,
        "claim_boundary": "Curation pool only; no model predictions were generated.",
    }
    write_provenance(output_dir / "homology-audit.json", payload)


def main() -> None:
    sequences = build_sequences()
    edmond = edmond_records()
    fractal = fractal_records()

    edmond_ids = sorted({str(record["sequence_id"]) for record in edmond})
    fractal_ids = sorted({str(record["sequence_id"]) for record in fractal})
    apply_homology_status(edmond, EDMOND / "homology" / "homology_hits.tsv")
    apply_homology_status(fractal, FRACTAL / "homology" / "homology_hits.tsv")
    write_records(EDMOND / "candidate_records.csv", edmond)
    write_fasta(EDMOND / "construct_sequences.fasta", edmond_ids, sequences)
    write_provenance(
        EDMOND / "provenance.json",
        {
            "source_id": "edmond-S0HJ48",
            "dataset_doi": "10.17617/3.S0HJ48",
            "article_doi": "10.1038/s41467-024-54408-6",
            "article_published": "2024-12-03",
            "dataset_license": "CC0-1.0",
            "article_license": "CC-BY-4.0",
            "record_count": len(edmond),
            "unique_construct_sequences": len(edmond_ids),
            "variable_substrates": sorted(SUBSTRATES),
            "kinetics_source": "raw/41467_2024_54408_MOESM1_ESM.pdf, Supplementary Tables S3 and 4",
            "sequence_source": "raw/41467_2024_54408_MOESM3_ESM.docx",
            "raw_trace_sources": [
                "raw/Kinetic_Acomosus_and_variant.xlsx",
                "raw/Kinetic_Anc3a_Anc3b_Anc3aplus5.xlsx",
                "raw/Kinetic_Dpimensis_SEC_fractions.xlsx",
                "raw/Kinetic_Msulfidovorans_and_variant.xlsx",
            ],
            "assay": {
                "method": "DTNB colorimetric citrate synthase assay at 412 nm",
                "buffer": "50 mM Tris pH 7.5, 10 mM KCl, 0.1 mg/mL DTNB",
                "temperature_C": 25,
                "enzyme_concentration_nM": 25,
                "fit_software": "GraphPad Prism 8.4.3",
            },
            "normalization": "Author-reported kcat values in s-1; no refitting or unit conversion",
            "known_source_anomaly": "Supplementary Table 4 reports anc3a+5 oxaloacetate kcat as 23.3 +/- 61.1 s-1. The value and SEM are retained verbatim and flagged per record.",
            "model_predictions_run": False,
        },
    )

    write_records(FRACTAL / "candidate_records.csv", fractal)
    write_fasta(FRACTAL / "construct_sequences.fasta", fractal_ids, sequences)
    write_provenance(
        FRACTAL / "provenance.json",
        {
            "source_id": "europepmc-PMC11041685",
            "stable_record_url": "https://europepmc.org/articles/PMC11041685",
            "article_doi": "10.1038/s41586-024-07287-2",
            "article_published": "2024-04-10",
            "license": "CC-BY-4.0",
            "record_count": len(fractal),
            "unique_construct_sequences": len(fractal_ids),
            "variable_substrates": sorted(SUBSTRATES),
            "kinetics_source": "../edmond-S0HJ48/raw/41586_2024_7287_MOESM1_ESM.pdf, Supplementary Table 1",
            "sequence_source": "../edmond-S0HJ48/raw/41586_2024_7287_MOESM1_ESM.pdf, Supplementary Table 2 (ABB56644.1); exact 1,158-nt WT DNA independently matches the later sequence DOCX",
            "assay": {
                "method": "DTNB colorimetric citrate synthase assay at 412 nm",
                "buffer": "50 mM Tris pH 7.5, 10 mM KCl, 0.1 mg/mL DTNB",
                "temperature_C": 25,
                "enzyme_concentration_nM": 25,
                "fit_software": "GraphPad Prism 8.4.3",
            },
            "normalization": "Author-reported kcat values in s-1; no refitting or unit conversion",
            "model_predictions_run": False,
        },
    )
    write_homology_audit(
        EDMOND,
        edmond,
        EDMOND / "construct_sequences.fasta",
        "10.1038/s41467-024-54408-6",
        "doi:10.17617/3.S0HJ48",
    )
    write_homology_audit(
        FRACTAL,
        fractal,
        FRACTAL / "construct_sequences.fasta",
        "10.1038/s41586-024-07287-2",
        "https://europepmc.org/articles/PMC11041685",
    )


if __name__ == "__main__":
    main()
