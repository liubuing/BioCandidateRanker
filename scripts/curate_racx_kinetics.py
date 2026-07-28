from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12838360"
FAMILY_CAP = 20

SUBSTRATES = {
    "L-Lys": (5962, "C(CCN)C[C@@H](C(=O)O)N"),
    "D-Lys": (57449, "C(CCN)C[C@H](C(=O)O)N"),
    "L-Trp": (6305, "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N"),
    "L-Thr": (6288, "C[C@H]([C@@H](C(=O)O)N)O"),
    "L-Val": (6287, "CC(C)[C@@H](C(=O)O)N"),
    "L-Leu": (6106, "CC(C)C[C@@H](C(=O)O)N"),
    "L-Ile": (6306, "CC[C@H](C)[C@@H](C(=O)O)N"),
    "L-Ala": (5950, "C[C@@H](C(=O)O)N"),
    "L-Tyr": (6057, "C1=CC(=CC=C1C[C@@H](C(=O)O)N)O"),
    "L-Arg": (6322, "C(C[C@@H](C(=O)O)N)CN=C(N)N"),
    "L-Asn": (6267, "C([C@@H](C(=O)O)N)C(=O)N"),
    "L-Ser": (5951, "C([C@@H](C(=O)O)N)O"),
    "L-Orn": (6262, "C(C[C@@H](C(=O)O)N)CN"),
    "L-Phe": (6140, "C1=CC=C(C=C1)C[C@@H](C(=O)O)N"),
    "L-Gln": (5961, "C(CC(=O)N)[C@@H](C(=O)O)N"),
    "L-Pro": (145742, "C1C[C@H](NC1)C(=O)O"),
}


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
    wt = read_fasta(SOURCE / "reference_sequence.fasta")
    variants = {
        "racx-wt": ("WT", wt),
        "racx-a79c": ("A79C", substitute(wt, 79, "A", "C")),
        "racx-n80a": ("N80A", substitute(wt, 80, "N", "A")),
        "racx-t81a": ("T81A", substitute(wt, 81, "T", "A")),
        "racx-n121a": ("N121A", substitute(wt, 121, "N", "A")),
        "racx-t124a": ("T124A", substitute(wt, 124, "T", "A")),
        "racx-c193s": ("C193S", substitute(wt, 193, "C", "S")),
    }
    return variants


def build_records() -> list[dict[str, object]]:
    table_1 = [
        ("L-Lys", 1.6, 0.5, 241.8, 12.4, 151.2),
        ("L-Trp", 0.3, 0.2, 44.2, 1.0, 140.6),
        ("L-Thr", 1.1, 0.7, 40.9, 1.1, 36.9),
        ("L-Val", 23.9, 3.1, 323.0, 20.3, 13.5),
        ("L-Leu", 0.6, 0.2, 57.8, 5.8, 96.9),
        ("L-Ile", 17.1, 5.8, 140.0, 17.4, 8.2),
        ("L-Ala", 1.0, 0.5, 113.3, 15.6, 113.7),
        ("L-Tyr", 0.6, 0.2, 56.1, 3.1, 94.5),
        ("L-Arg", 0.7, 0.1, 42.1, 6.7, 60.0),
        ("L-Asn", 0.6, 0.5, 84.1, 13.7, 141.9),
        ("L-Ser", 3.4, 0.5, 71.1, 15.7, 20.7),
        ("L-Orn", 0.9, 0.5, 54.6, 6.1, 59.7),
        ("L-Phe", 2.5, 1.4, 63.3, 14.6, 25.1),
        ("L-Gln", 0.5, 0.1, 56.7, 3.6, 121.6),
        ("L-Pro", 0.6, 0.2, 51.1, 9.5, 82.6),
    ]
    table_3 = [
        ("WT", "racx-wt", "D-Lys", 25.3, 10.6, 291.2, 50.5, 11.5),
        ("A79C", "racx-a79c", "L-Lys", 25.1, 6.1, 403.7, 40.6, 16.1),
        ("A79C", "racx-a79c", "D-Lys", 40.0, 13.7, 662.7, 109.3, 16.6),
        ("C193S", "racx-c193s", "L-Lys", 14.6, 2.7, 596.7, 39.2, 40.9),
        ("N80A", "racx-n80a", "L-Lys", 3.7, 2.1, 128.5, 17.6, 34.7),
        ("T81A", "racx-t81a", "L-Lys", 7.1, 1.8, 205.9, 15.1, 29.0),
        ("N121A", "racx-n121a", "L-Lys", 8.7, 1.6, 198.0, 11.0, 22.8),
        ("T124A", "racx-t124a", "L-Lys", 6.2, 0.8, 591.3, 20.9, 95.4),
    ]
    rows = [
        ("Table 1", "WT", "racx-wt", substrate, *values)
        for substrate, *values in table_1
    ]
    rows.extend(("Table 3", construct, sequence_id, substrate, *values) for construct, sequence_id, substrate, *values in table_3)

    records = []
    for index, row in enumerate(rows, 1):
        source_table, construct, sequence_id, substrate, km, km_error, kcat, kcat_error, efficiency = row
        cid, smiles = SUBSTRATES[substrate]
        anomaly = ""
        if source_table == "Table 3" and construct == "WT" and substrate == "D-Lys":
            anomaly = "Narrative text reports kcat/Km 17.8, while Table 3 reports 11.5; 291.2/25.3 = 11.5, so the table value is retained."
        elif source_table == "Table 3" and construct == "A79C" and substrate == "D-Lys":
            anomaly = "Narrative text reports kcat/Km 25.6, while Table 3 reports 16.6; 662.7/40.0 = 16.6, so the table value is retained."
        records.append(
            {
                "candidate_id": f"racx-{index:03d}",
                "article_doi": "10.1128/aem.02015-25",
                "source_table": source_table,
                "source_row": f"RacX {construct}; {substrate}",
                "organism": "Halocola ammonii DA487T",
                "sequence_id": sequence_id,
                "construct": construct,
                "variable_substrate": substrate,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_mM": km,
                "km_error_mM": km_error,
                "kcat_per_km_mM-1_s-1": efficiency,
                "assay_pH": 7.5,
                "assay_temperature_C": 37,
                "enzyme_mass_ug": 2.5,
                "reaction_volume_uL": 50,
                "replicates": 3,
                "error_type": "not explicitly stated for kinetic tables",
                "source_anomaly": anomaly,
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def apply_selection_status(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if not hit_path.exists():
        return
    hit_ids = {
        line.split("\t", 1)[0]
        for line in hit_path.read_text(encoding="utf-8").splitlines()
        if line
    }
    accepted = 0
    for record in records:
        if record["sequence_id"] in hit_ids:
            record["status_at_normalization"] = "excluded_homology"
        elif accepted < FAMILY_CAP:
            record["status_at_normalization"] = "accepted_homology_cold_pool"
            accepted += 1
        else:
            record["status_at_normalization"] = "excluded_family_cap"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_outputs() -> None:
    sequences = build_sequences()
    records = build_records()
    apply_selection_status(records)

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (construct, sequence) in sequences.items():
            handle.write(f">{sequence_id} RacX_{construct}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    provenance = {
        "source_id": "europepmc-PMC12838360",
        "stable_record_url": "https://europepmc.org/articles/PMC12838360",
        "article_doi": "10.1128/aem.02015-25",
        "article_published": "2025-12-22",
        "license": "CC-BY-4.0",
        "raw_record_count": len(records),
        "family_cap": FAMILY_CAP,
        "accepted_record_count": len(accepted),
        "accepted_unique_substrates": len({record["variable_substrate"] for record in accepted}),
        "unique_construct_sequences": len(sequences),
        "kinetics_source": "raw/PMC12838360-fullText.xml, Tables 1 and 3",
        "sequence_source": "raw/aem.02015-25-s0001.pdf, RacX gene sequence; translated from 684 nt",
        "assay": {
            "method": "Endpoint racemase assay verified linear through 60 min; Michaelis-Menten fit",
            "buffer": "50 mM Tris-HCl pH 7.5",
            "temperature_C": 37,
            "enzyme_mass_ug": 2.5,
            "reaction_volume_uL": 50,
            "duration_min": 30,
            "replicates": 3,
        },
        "selection_rule": "After homology exclusion, retain the first 20 records in source-table and source-row order. The rule does not inspect kinetic values.",
        "normalization": "Author-reported kcat values in s-1; no refitting or unit conversion",
        "known_source_anomalies": [
            "For WT D-Lys, narrative kcat/Km is 17.8 but Table 3 and kcat/Km arithmetic give 11.5.",
            "For A79C D-Lys, narrative kcat/Km is 25.6 but Table 3 and kcat/Km arithmetic give 16.6."
        ],
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if hit_path.exists():
        hit_rows = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        audit = {
            "audited_on": "2026-07-21",
            "source_id": "europepmc-PMC12838360",
            "article_doi": "10.1128/aem.02015-25",
            "method": "MMseqs2 easy-search",
            "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9",
            "candidate_records": len(records),
            "unique_sequences": len(sequences),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_rows}),
            "homology_hits_sha256": sha256(hit_path),
            "accepted_records_after_family_cap": len(accepted),
            "accepted_unique_substrates": len({record["variable_substrate"] for record in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Accepted curation subset only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )


if __name__ == "__main__":
    write_outputs()
