from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "dryad-ncjsxkt6b"

GLUCOSE_SMILES = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"
ATP_SMILES = (
    "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)"
    "COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
)


def read_fasta(path: Path) -> str:
    lines = path.read_text(encoding="ascii").splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if not sequence or not sequence.isalpha() or not sequence.isupper():
        raise ValueError(f"Invalid FASTA sequence: {path}")
    return sequence


def substitute(sequence: str, changes: list[tuple[int, str, str]]) -> str:
    residues = list(sequence)
    for position, expected, replacement in changes:
        observed = residues[position - 1]
        if observed != expected:
            raise ValueError(f"Expected {expected}{position}, observed {observed}{position}")
        residues[position - 1] = replacement
    return "".join(residues)


def build_sequences() -> dict[str, tuple[str, str]]:
    ps_wt = read_fasta(SOURCE / "raw" / "H7CHS4.fasta")
    ec_wt = read_fasta(SOURCE / "raw" / "P0A6V8.fasta")
    return {
        "gk-ps-wt": ("H7CHS4", ps_wt),
        "gk-ps-c325s": ("H7CHS4 C325S", substitute(ps_wt, [(325, "C", "S")])),
        "gk-ec-wt": ("P0A6V8", ec_wt),
        "gk-ec-ds-s": (
            "P0A6V8 C20S/C65S/S309C",
            substitute(ec_wt, [(20, "C", "S"), (65, "C", "S"), (309, "S", "C")]),
        ),
    }


def build_records() -> list[dict[str, object]]:
    rows = [
        ("PsGK WT", "Pseudoalteromonas sp. AS-131", "gk-ps-wt", "Table 1", 1, 0.12, 0.01, 0.83, 0.14, 11.25, 0.28, 92.15, 5.17, "SD", 3),
        ("PsGK WT", "Pseudoalteromonas sp. AS-131", "gk-ps-wt", "Table 1", 25, 0.36, 0.02, 1.55, 0.17, 33.49, 0.77, 93.99, 3.98, "SD", 3),
        ("PsGK WT", "Pseudoalteromonas sp. AS-131", "gk-ps-wt", "Table 1", 40, 0.60, 0.11, 1.54, 0.20, 40.01, 5.23, 67.93, 4.98, "SD", 3),
        ("PsGK C325S", "Pseudoalteromonas sp. AS-131", "gk-ps-c325s", "Supplementary Table S2", 1, 0.10, 0.03, 0.52, 0.04, 15.74, 0.46, 159.88, 33.79, "not explicitly stated", ""),
        ("PsGK C325S", "Pseudoalteromonas sp. AS-131", "gk-ps-c325s", "Supplementary Table S2", 25, 0.21, 0.03, 0.96, 0.09, 48.57, 3.23, 237.85, 18.73, "not explicitly stated", ""),
        ("PsGK C325S", "Pseudoalteromonas sp. AS-131", "gk-ps-c325s", "Supplementary Table S2", 40, 0.29, 0.17, 1.86, 0.25, 100.42, 1.85, 247.99, 1.74, "not explicitly stated", ""),
        ("EcGK WT", "Escherichia coli K-12", "gk-ec-wt", "Table 1", 1, 1.24, 0.08, 0.78, 0.03, 16.75, 0.84, 13.57, 0.21, "SD", 3),
        ("EcGK WT", "Escherichia coli K-12", "gk-ec-wt", "Table 1", 25, 0.73, 0.06, 1.11, 0.08, 51.36, 2.00, 70.41, 3.07, "SD", 3),
        ("EcGK WT", "Escherichia coli K-12", "gk-ec-wt", "Table 1", 40, 0.81, 0.15, 2.25, 0.28, 70.42, 6.88, 88.18, 8.00, "SD", 3),
        ("EcGK DS-S", "Escherichia coli K-12", "gk-ec-ds-s", "Supplementary Table S2", 1, 0.83, 0.12, 0.42, 0.03, 23.45, 2.87, 28.29, 0.74, "not explicitly stated", ""),
        ("EcGK DS-S", "Escherichia coli K-12", "gk-ec-ds-s", "Supplementary Table S2", 25, 0.37, 0.03, 0.78, 0.15, 42.94, 2.28, 116.65, 3.92, "not explicitly stated", ""),
        ("EcGK DS-S", "Escherichia coli K-12", "gk-ec-ds-s", "Supplementary Table S2", 40, 0.91, 0.13, 1.46, 0.19, 63.08, 6.89, 69.12, 2.66, "not explicitly stated", ""),
    ]
    records = []
    for index, row in enumerate(rows, 1):
        (
            construct,
            organism,
            sequence_id,
            source_table,
            temperature,
            km_glucose,
            km_glucose_error,
            km_atp,
            km_atp_error,
            kcat,
            kcat_error,
            efficiency,
            efficiency_error,
            error_type,
            replicates,
        ) = row
        records.append(
            {
                "candidate_id": f"gk-{index:03d}",
                "article_doi": "10.1111/febs.70367",
                "source_table": source_table,
                "source_row": f"{construct}; {temperature} C",
                "organism": organism,
                "sequence_id": sequence_id,
                "construct": construct,
                "reaction_substrate_1": "D-glucose",
                "substrate_1_pubchem_cid": 5793,
                "substrate_1_isomeric_smiles": GLUCOSE_SMILES,
                "reaction_substrate_2": "ATP",
                "substrate_2_pubchem_cid": 5957,
                "substrate_2_isomeric_smiles": ATP_SMILES,
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_glucose_mM": km_glucose,
                "km_glucose_error_mM": km_glucose_error,
                "km_atp_mM": km_atp,
                "km_atp_error_mM": km_atp_error,
                "kcat_per_km_glucose_mM-1_s-1": efficiency,
                "kcat_per_km_glucose_error": efficiency_error,
                "assay_pH": 7.6,
                "assay_temperature_C": temperature,
                "error_type": error_type,
                "biological_replicates": replicates,
                "assay_concentration_note": "Methods state glucose and ATP varied from 0.17-1.0 mM, while also listing 5 mM glucose and 2 mM ATP; fixed-cosubstrate assignment is not explicit.",
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def homology_hit_ids() -> set[str] | None:
    path = SOURCE / "homology" / "homology_hits.tsv"
    if not path.exists():
        return None
    return {
        line.split("\t", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_homology_audit(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.exists() or not cluster_path.exists():
        return
    hit_rows = [
        line.split("\t")
        for line in hit_path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    families = {
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
        "source_id": "dryad-ncjsxkt6b",
        "article_doi": "10.1111/febs.70367",
        "dataset_doi": "10.5061/dryad.ncjsxkt6b",
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
        "homology_hit_sequences": len({row[0] for row in hit_rows}),
        "exact_sequence_overlap": len({row[0] for row in hit_rows if float(row[2]) == 1.0}),
        "candidate_mmseqs_families": len(families),
        "homology_hits_sha256": sha256(hit_path),
        "family_cluster_sha256": sha256(cluster_path),
        "accepted_records": len(accepted),
        "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}),
        "accepted_unique_substrates": 0 if not accepted else 2,
        "readiness_gate_passes": False,
        "claim_boundary": "Excluded from the temporal pool because every construct has a development-corpus homology hit; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="ascii"
    )


def write_outputs() -> None:
    sequences = build_sequences()
    records = build_records()
    hits = homology_hit_ids()
    if hits is not None:
        for record in records:
            record["status_at_normalization"] = (
                "excluded_homology"
                if record["sequence_id"] in hits
                else "accepted_homology_cold_pool"
            )

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    with (SOURCE / "construct_sequences.fasta").open(
        "w", encoding="ascii", newline="\n"
    ) as handle:
        for sequence_id, (accession, sequence) in sequences.items():
            handle.write(f">{sequence_id} {accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    provenance = {
        "source_id": "dryad-ncjsxkt6b",
        "dataset_doi": "10.5061/dryad.ncjsxkt6b",
        "article_doi": "10.1111/febs.70367",
        "article_published": "2025-12-15",
        "dataset_license": "CC0-1.0",
        "article_license": "CC-BY-4.0",
        "record_count": len(records),
        "unique_construct_sequences": len(sequences),
        "kinetics_sources": [
            "raw/PMC13147314-fullText.xml, Table 1",
            "raw/FEBS-293-2617-s001.pdf, Supplementary Table S2",
        ],
        "sequence_sources": ["raw/H7CHS4.fasta", "raw/P0A6V8.fasta"],
        "assay": {
            "method": "Two-step glucokinase/G6PD assay; NADPH measured at 340 nm; Lineweaver-Burk analysis",
            "buffer": "20 mM Tris-HCl pH 7.6, 10 mM MgCl2, 1 mM DTT, 1.2 mM NADP+, 50 uM G6PD",
            "temperatures_C": [1, 25, 40],
            "concentration_caveat": "Methods state glucose and ATP varied from 0.17-1.0 mM, while also listing 5 mM glucose and 2 mM ATP.",
        },
        "normalization": "Author-reported reaction-level kcat in s-1; one record per construct-temperature, not duplicated across the two reported Km values",
        "error_note": "Table 1 explicitly reports SD from three independent experiments. Supplementary Table S2 reports +/- values without independently naming the error type.",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )
    write_homology_audit(records)


if __name__ == "__main__":
    write_outputs()
