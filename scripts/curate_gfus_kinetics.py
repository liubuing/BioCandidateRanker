from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12362431"
LEADER = "MHHHHHHSSGVDLGTENLYFQS"


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )


def substitute(sequence: str, position: int, expected: str, replacement: str) -> str:
    if sequence[position - 1] != expected:
        raise ValueError(f"Expected {expected}{position}, found {sequence[position - 1]}{position}")
    return sequence[: position - 1] + replacement + sequence[position:]


def build_sequences() -> dict[str, tuple[str, str]]:
    canonical = read_fasta(SOURCE / "raw" / "Q13630.fasta")
    variants = {
        "gfus-wt": ("WT", canonical),
        "gfus-h186k": ("H186K", substitute(canonical, 186, "H", "K")),
        "gfus-c116s": ("C116S", substitute(canonical, 116, "C", "S")),
        "gfus-y143f": ("Y143F", substitute(canonical, 143, "Y", "F")),
    }
    return {
        sequence_id: (variant, LEADER + sequence[7:])
        for sequence_id, (variant, sequence) in variants.items()
    }


def build_records() -> list[dict[str, object]]:
    rows = [
        ("gfus-wt", "WT", 2.11, 0.03, 1.7, 0.2, 0.75, 3.18, ""),
        ("gfus-h186k", "H186K", 0.36, 0.03, 0.8, 0.2, "", 31.8, ""),
        ("gfus-c116s", "C116S", 0.0029, 0.0003, 2.6, 0.7, "", 609, "Approximately 30% GDP-L-fucose and 70% GDP-D-altrose are formed; NADPH-consumption kcat is not product-specific."),
        ("gfus-y143f", "Y143F", 0.025, 0.002, 3.3, 0.6, "", 609, ""),
    ]
    records = []
    for index, (sequence_id, variant, kcat, kcat_error, km, km_error, km_nadph, enzyme_nm, note) in enumerate(rows, 1):
        records.append(
            {
                "candidate_id": f"gfus-{index:03d}",
                "article_doi": "10.1021/acscatal.5c02722",
                "source_table": "Table 1",
                "source_row": f"human GFUS {variant}",
                "organism": "Homo sapiens",
                "sequence_id": sequence_id,
                "construct": f"pNIC28_hgfs; His6/TEV leader + Q13630 residues 8-321; {variant}",
                "variable_substrate": "GDP-4-keto-6-deoxy-D-mannose",
                "substrate_pubchem_cid": 135398621,
                "substrate_isomeric_smiles": "C[C@@H]1C(=O)[C@@H]([C@@H]([C@H](O1)OP(=O)(O)OP(=O)(O)OC[C@@H]2[C@H]([C@H]([C@@H](O2)N3C=NC4=C3N=C(NC4=O)N)O)O)O)O",
                "cosubstrate": "NADPH",
                "cosubstrate_pubchem_cid": 5884,
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_substrate_uM": km,
                "km_substrate_error_uM": km_error,
                "km_nadph_uM": km_nadph,
                "assay_pH": 8.0,
                "assay_temperature_C": 37,
                "enzyme_concentration_nM": enzyme_nm,
                "error_type": "SD from three independent fits",
                "replicates": 3,
                "source_note": note,
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_outputs() -> None:
    sequences = build_sequences()
    records = build_records()
    apply_homology_status(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (variant, sequence) in sequences.items():
            handle.write(f">{sequence_id} Q13630_{variant} | exact pNIC28_hgfs construct\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [record for record in records if record["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "europepmc-PMC12362431",
        "stable_record_url": "https://europepmc.org/articles/PMC12362431",
        "article_doi": "10.1021/acscatal.5c02722",
        "article_published": "2025-07-29",
        "license": "CC-BY-4.0",
        "record_count": len(records),
        "unique_construct_sequences": len(sequences),
        "kinetics_source": "raw/PMC12362431-fullText.xml, Table 1",
        "sequence_source": "raw/Q13630.fasta plus article-defined pNIC28_hgfs leader and residues 8-321",
        "assay": {
            "method": "Continuous NADPH A340 Michaelis-Menten assay",
            "buffer": "10 mM Tris, 25 mM NaCl, pH 8.0",
            "temperature_C": 37,
            "volume_uL": 500,
            "substrate_range_uM": "1-50",
            "nadph_range_uM": "2-30",
            "replicates": 3,
            "error_type": "SD from three independent fits",
        },
        "excluded_data": "Stopped-flow/global-simulation duplicates, isotope-effect ratios, and non-finite C116A/H186A rows",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
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
            "source_id": "europepmc-PMC12362431",
            "article_doi": "10.1021/acscatal.5c02722",
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
