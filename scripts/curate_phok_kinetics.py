from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12444516"
PDF = SOURCE / "raw" / "ao5c05616_si_001.pdf"

SUBSTRATES = {
    "4-nitrophenyl phosphate": (
        378,
        "C1=CC(=CC=C1[N+](=O)[O-])OP(=O)(O)O",
    ),
    "phenyl phosphate": (12793, "C1=CC=C(C=C1)OP(=O)(O)O"),
}

ALIASES = {
    "Sm-PhoK": "sm-phok",
    "Sb-PhoK": "sb-phok",
    "SmTDK1-PhoK": "smtdk1-phok",
    "Sp-PhoK": "st-phok",
    "Sy-PhoK": "sy-phok",
    "Na-PhoK": "na-phok",
    "No-PhoK": "no-phok",
    "SmSRS2-PhoK": "smsrs2-phok",
    "Ng-Phok": "ng-phok",
}

ORGANISMS = {
    "sb-phok": "Sphingobium sp. TCM1",
    "smtdk1-phok": "Sphingomonas sp. TDK1",
    "st-phok": "Sphingopyxis terrae",
    "sy-phok": "Sphingobium yanoikuyae",
    "na-phok": "Novosphingobium aromaticivorans",
    "no-phok": "Novosphingobium sp. EMRT2",
    "sm-phok": "Sphingomonas sp. BSAR-1",
    "ng-phok": "Novosphingobium guangzhouense",
    "smsrs2-phok": "Sphingomonas sp. SRS2",
}


def extract_sequences() -> dict[str, str]:
    text = fitz.open(PDF)[3].get_text()
    sequences = {sequence_id: "" for sequence_id in ALIASES.values()}
    names = "|".join(re.escape(name) for name in ALIASES)
    pattern = re.compile(rf"^({names})\s+(?:\d+)?([A-Z-]+)")
    for line in text.splitlines():
        match = pattern.match(line)
        if match:
            sequences[ALIASES[match.group(1)]] += match.group(2).replace("-", "")
    expected_lengths = {
        "sm-phok": 540,
        "sb-phok": 555,
        "smtdk1-phok": 526,
        "st-phok": 534,
        "sy-phok": 540,
        "na-phok": 539,
        "no-phok": 541,
        "smsrs2-phok": 529,
        "ng-phok": 542,
    }
    observed = {name: len(sequence) for name, sequence in sequences.items()}
    if observed != expected_lengths:
        raise ValueError(f"Unexpected Figure S2 sequence lengths: {observed}")
    return sequences


def build_records() -> list[dict[str, object]]:
    rows = [
        ("sb-phok", "4-nitrophenyl phosphate", 416, 3, 69, 2, 6.1e6, 0.2e6, 0.165, 2.75, 8.6, ""),
        ("sb-phok", "phenyl phosphate", 460, 10, 650, 30, 7.1e5, 0.4e5, 3.3, 2.32, 7.25, ""),
        ("smtdk1-phok", "4-nitrophenyl phosphate", 460, 5, 450, 10, 1.02e6, 0.03e6, 0.20, 2.75, 8.6, ""),
        ("smtdk1-phok", "phenyl phosphate", 365, 7, 600, 20, 6.1e5, 0.3e5, 2.0, 2.0, 37.5, ""),
        ("st-phok", "4-nitrophenyl phosphate", 440, 5, 78, 3, 5.7e6, 0.2e6, 0.13, 2.75, 8.6, "Table S2 prints enzyme concentration as 0.13 Nm; normalized to 0.13 nM from context."),
        ("st-phok", "phenyl phosphate", 314, 6, 81, 5, 3.9e6, 0.2e6, 1.8, 1.5, 4.8, ""),
        ("sy-phok", "4-nitrophenyl phosphate", 1230, 10, 176, 4, 7.0e6, 0.2e6, 0.066, 2.75, 12.9, ""),
        ("sy-phok", "phenyl phosphate", 930, 30, 670, 40, 1.4e6, 0.1e6, 0.63, 2.32, 7.25, ""),
        ("na-phok", "4-nitrophenyl phosphate", 303, 3, 89, 3, 3.4e6, 0.1e6, 0.31, 2.75, 8.6, ""),
        ("na-phok", "phenyl phosphate", 336, 6, 390, 20, 8.7e5, 0.4e5, 3.0, 2.32, 7.25, ""),
        ("no-phok", "4-nitrophenyl phosphate", 621, 6, 153, 5, 4.1e6, 1.1e6, 0.16, 2.75, 8.6, "The reported kcat/Km error (1.1e6) is unusually large; retained without correction."),
        ("no-phok", "phenyl phosphate", 259, 6, 300, 20, 8.7e5, 0.6e5, 3.2, 2.32, 7.25, ""),
        ("sm-phok", "4-nitrophenyl phosphate", 1370, 10, 174, 4, 7.9e6, 0.2e6, 0.094, 2.75, 12.9, ""),
        ("sm-phok", "phenyl phosphate", 1100, 10, 600, 20, 1.84e6, 0.06e6, 0.77, 2.32, 7.25, ""),
        ("ng-phok", "4-nitrophenyl phosphate", 297, 4, 90, 4, 3.3e6, 0.1e6, 0.24, 2.75, 8.6, ""),
        ("ng-phok", "phenyl phosphate", 190, 3, 250, 9, 7.5e5, 0.3e5, 2.3, 2.32, 7.25, ""),
        ("smsrs2-phok", "4-nitrophenyl phosphate", 12.0, 0.1, 88, 4, 1.36e5, 0.07e5, 11.3, 2.75, 12.9, ""),
        ("smsrs2-phok", "phenyl phosphate", 16.6, 0.2, 104, 5, 1.59e5, 0.07e5, 22.5, 2.0, 6.25, ""),
    ]
    records = []
    for index, row in enumerate(rows, 1):
        sequence_id, substrate, kcat, kcat_error, km, km_error, efficiency, efficiency_error, enzyme_nm, maximum_mm, minimum_um, anomaly = row
        cid, smiles = SUBSTRATES[substrate]
        construct = "native/core sequence; untagged"
        if sequence_id != "smtdk1-phok":
            construct = "native/core sequence; pET29a C-terminal 6xHis fusion with unreported linker/tag sequence"
        records.append(
            {
                "candidate_id": f"phok-{index:03d}",
                "article_doi": "10.1021/acsomega.5c05616",
                "source_table": "Supplementary Table S2",
                "source_row": f"{sequence_id}; {substrate}",
                "organism": ORGANISMS[sequence_id],
                "sequence_id": sequence_id,
                "construct": construct,
                "variable_substrate": substrate,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_uM": km,
                "km_error_uM": km_error,
                "kcat_per_km_M-1_s-1": efficiency,
                "kcat_per_km_error": efficiency_error,
                "assay_pH": 9.0,
                "assay_temperature_C": "",
                "enzyme_concentration_nM": enzyme_nm,
                "max_substrate_concentration_mM": maximum_mm,
                "min_substrate_concentration_uM": minimum_um,
                "error_type": "not reported",
                "replicates": "",
                "source_anomaly": anomaly,
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
    sequences = extract_sequences()
    records = build_records()
    apply_homology_status(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in sequences.items():
            handle.write(f">{sequence_id} {ORGANISMS[sequence_id]} | Figure S2 native/core sequence\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [record for record in records if record["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "europepmc-PMC12444516",
        "stable_record_url": "https://europepmc.org/articles/PMC12444516",
        "article_doi": "10.1021/acsomega.5c05616",
        "article_published": "2025-08-29",
        "license": "CC-BY-4.0",
        "record_count": len(records),
        "unique_core_sequences": len(sequences),
        "unique_absolute_kcat_substrates": len(SUBSTRATES),
        "kinetics_source": "raw/ao5c05616_si_001.pdf, Supplementary Table S2",
        "sequence_source": "raw/ao5c05616_si_001.pdf, Figure S2 ungapped sequences",
        "construct_caveat": "Eight proteins were pET29a C-terminal 6xHis fusions, but linker/tag residues were not reported. SmTDK1 retained its native stop codon and was untagged.",
        "assay": {
            "method": "UV/Vis initial-rate Michaelis-Menten fits over 10 min",
            "buffer": "50 mM CHES pH 9.0",
            "reaction_volume_uL": 250,
            "temperature_C": None,
            "error_type": "not reported",
            "replicates": None,
        },
        "normalization": "Author-reported kcat in s-1; no refitting or unit conversion",
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
            "source_id": "europepmc-PMC12444516",
            "article_doi": "10.1021/acsomega.5c05616",
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
