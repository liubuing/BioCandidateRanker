from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "zenodo-14055918"


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
    wt = read_fasta(SOURCE / "raw" / "Q5S007.fasta")
    return {
        "lrrk2-wt": ("WT", wt),
        "lrrk2-r1441g": ("R1441G", substitute(wt, 1441, "R", "G")),
        "lrrk2-g2019s": ("G2019S", substitute(wt, 2019, "G", "S")),
        "lrrk2-k1906m": ("K1906M", substitute(wt, 1906, "K", "M")),
        "lrrk2-t1343a": ("T1343A", substitute(wt, 1343, "T", "A")),
    }


def build_records() -> list[dict[str, object]]:
    rows = [
        ("lrrk2-wt", "WT", False, 554, 62, 0.36, 0.02, 5),
        ("lrrk2-wt", "WT", True, 1036, 169, 0.14, 0.01, 2),
        ("lrrk2-r1441g", "R1441G", False, 272, 28, 0.43, 0.02, 4),
        ("lrrk2-g2019s", "G2019S", False, 867, 110, 0.37, 0.03, 4),
        ("lrrk2-k1906m", "K1906M", False, 181, 58, 0.28, 0.03, 3),
        ("lrrk2-t1343a", "T1343A", False, 265, 25, 0.10, 0.01, 4),
        ("lrrk2-t1343a", "T1343A", True, 328, 34, 0.10, 0.01, 4),
    ]
    records = []
    for index, (sequence_id, variant, atp, km, km_error, kcat_min, error_min, replicates) in enumerate(rows, 1):
        condition = "2 h ATP preincubation at 30 C" if atp else "no ATP preincubation"
        records.append(
            {
                "candidate_id": f"lrrk2-{index:03d}",
                "article_doi": "10.7554/eLife.91083",
                "dataset_doi": "10.5281/zenodo.14055918",
                "source_table": "Table 1",
                "source_row": f"full-length LRRK2 {variant}; {condition}",
                "organism": "Homo sapiens",
                "sequence_id": sequence_id,
                "construct": f"full-length Q5S007 {variant}; N-terminal FLAG/tandem-STREP-II fusion with complete linker sequence not recovered",
                "variable_substrate": "GTP",
                "substrate_pubchem_cid": 135398633,
                "substrate_isomeric_smiles": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
                "endpoint": "kcat_s-1",
                "source_kcat_min-1": kcat_min,
                "source_kcat_error_min-1": error_min,
                "kcat_s-1": kcat_min / 60,
                "kcat_error_s-1": error_min / 60,
                "km_gtp_uM": km,
                "km_gtp_error_uM": km_error,
                "enzyme_concentration_uM": 0.1,
                "gtp_range_uM": "0,25,75,150,250,500,1000,2000,3000,5000",
                "assay_pH": "",
                "assay_temperature_C": "",
                "atp_preincubation": atp,
                "atp_preincubation_temperature_C": 30 if atp else "",
                "atp_preincubation_duration_h": 2 if atp else "",
                "error_type": "fit standard error",
                "replicates": replicates,
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
            handle.write(f">{sequence_id} Q5S007_{variant} | full-length native sequence\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [record for record in records if record["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "zenodo-14055918",
        "dataset_doi": "10.5281/zenodo.14055918",
        "article_doi": "10.7554/eLife.91083",
        "article_published": "2024-12-19",
        "dataset_license": "CC-BY-4.0",
        "article_license": "CC-BY-4.0",
        "record_count": len(records),
        "unique_native_sequences": len(sequences),
        "kinetics_source": "raw/PMC11658767-fullText.xml, Table 1",
        "raw_curve_source": "raw/Michaelis_Menten_Kinetics_RAW_data.zip",
        "sequence_source": "raw/Q5S007.fasta plus explicit variants",
        "construct_caveat": "The paper identifies an N-terminal FLAG/tandem-STREP-II fusion, but the complete fusion linker sequence was not recovered; native full-length Q5S007 sequences are used.",
        "assay": {
            "method": "GDP production monitored by reversed-phase C18 HPLC and Michaelis-Menten fitting",
            "enzyme_concentration_uM": 0.1,
            "gtp_range_uM": [0, 25, 75, 150, 250, 500, 1000, 2000, 3000, 5000],
            "reaction_pH": None,
            "reaction_temperature_C": None,
            "fit_software": "GraFit 5.0.13",
            "error_type": "fit standard error",
        },
        "unit_conversion": "kcat_s-1 = author-reported kcat_min-1 / 60; errors converted by the same factor",
        "excluded_data": "MBP-RocCOR construct rows and supplemental screening curves not represented in main Table 1",
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
            "source_id": "zenodo-14055918",
            "article_doi": "10.7554/eLife.91083",
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
