from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12284513"


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )


def build_sequences() -> dict[str, tuple[str, str]]:
    oletla = read_fasta(SOURCE / "raw" / "A0A1M7CBV6.fasta")
    oletje = read_fasta(SOURCE / "raw" / "E9NSU2.fasta")
    if oletla[177] != "I":
        raise ValueError(f"Expected OleTLA I178, found {oletla[177]}178")
    return {
        "oletla-wt": ("A0A1M7CBV6", oletla),
        "oletla-i178l": ("A0A1M7CBV6 I178L", oletla[:177] + "L" + oletla[178:]),
        "oletje-wt": ("E9NSU2", oletje),
    }


def build_records() -> list[dict[str, object]]:
    rows = [
        ("oletla-wt", "OleTLA WT", 2.5, 20.6, 0.5),
        ("oletla-wt", "OleTLA WT", 10.0, 88.4, 1.6),
        ("oletje-wt", "OleTJE WT", 2.5, 38.8, 0.9),
        ("oletje-wt", "OleTJE WT", 10.0, 55.1, 0.2),
        ("oletla-i178l", "OleTLA I178L", 10.0, 96.0, 2.5),
    ]
    records = []
    for index, (sequence_id, construct, ethanol, kcat_min, error_min) in enumerate(rows, 1):
        records.append(
            {
                "candidate_id": f"cyp152-{index:03d}",
                "article_doi": "10.1016/j.jbc.2025.110397",
                "source_table": "Table 1",
                "source_row": f"{construct}; C16 fatty acid; {ethanol:g}% ethanol",
                "organism": "Labilithrix luteola" if sequence_id.startswith("oletla") else "Jeotgalicoccus sp. ATCC 8456",
                "sequence_id": sequence_id,
                "construct": construct,
                "variable_substrate": "hexadecanoic acid",
                "substrate_pubchem_cid": 985,
                "substrate_isomeric_smiles": "CCCCCCCCCCCCCCCC(=O)O",
                "endpoint": "kcat_s-1",
                "source_kcat_min-1": kcat_min,
                "source_kcat_error_min-1": error_min,
                "kcat_s-1": kcat_min / 60,
                "kcat_error_s-1": error_min / 60,
                "assay_pH": 7.5,
                "assay_temperature_C": 25,
                "enzyme_concentration_uM": 0.3,
                "hydrogen_peroxide_uM": 400,
                "substrate_range_uM": "5-250",
                "ethanol_percent_vv": ethanol,
                "error_type": "SD",
                "replicates": 3,
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def apply_homology_status(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "all-constructs" / "homology_hits.tsv"
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_outputs() -> None:
    sequences = build_sequences()
    records = build_records()
    apply_homology_status(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (accession, sequence) in sequences.items():
            handle.write(f">{sequence_id} {accession} | native catalytic sequence\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [record for record in records if record["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "europepmc-PMC12284513",
        "stable_record_url": "https://europepmc.org/articles/PMC12284513",
        "article_doi": "10.1016/j.jbc.2025.110397",
        "article_published": "2025-06-19",
        "license": "CC-BY-4.0",
        "record_count": len(records),
        "unique_construct_sequences": len(sequences),
        "unique_substrates": 1,
        "kinetics_source": "raw/PMC12284513-fullText.xml, Table 1",
        "sequence_sources": ["raw/A0A1M7CBV6.fasta", "raw/E9NSU2.fasta"],
        "construct_caveat": "Native catalytic sequences are used. OleTLA PDB 9JQM includes vector-derived LEHHHHHH, but the article does not explicitly print the complete tagged chain for every kinetic construct.",
        "assay": {
            "method": "Initial-rate C15 alkene formation measured by GC-MS and Michaelis-Menten fitting",
            "buffer": "50 mM NaH2PO4, 300 mM NaCl, pH 7.5",
            "temperature_C": 25,
            "enzyme_concentration_uM": 0.3,
            "hydrogen_peroxide_uM": 400,
            "substrate_range_uM": "5-250",
            "replicates": 3,
            "error_type": "SD",
        },
        "unit_conversion": "kcat_s-1 = author-reported kcat_min-1 / 60; errors converted by the same factor",
        "excluded_data": "Endpoint fatty-acid scope yields, OleTJE L176I endpoint data, and the imported pre-2023 OleTJE/C14 kcat row",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    hit_path = SOURCE / "homology" / "all-constructs" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        audit = {
            "audited_on": "2026-07-21",
            "source_id": "europepmc-PMC12284513",
            "article_doi": "10.1016/j.jbc.2025.110397",
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9",
            "candidate_records": len(records),
            "unique_sequences": len(sequences),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}),
            "accepted_unique_substrates": len({record["variable_substrate"] for record in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    write_outputs()
