from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "biostudies-S-EPMC13033379"

SUBSTRATES = {
    "caffeoyl-CoA": (11966126, "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCSC(=O)/C=C/C4=CC(=C(C=C4)O)O)O"),
    "p-coumaroyl-CoA": (6440013, "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCSC(=O)/C=C/C4=CC=C(C=C4)O)O"),
    "shikimic acid": (8742, "C1[C@H]([C@@H]([C@@H](C=C1C(=O)O)O)O)O"),
    "quinic acid": (6508, "C1[C@H](C([C@@H](CC1(C(=O)O)O)O)O)O"),
    "3-hydroxyanthranilic acid": (86, "C1=CC(=C(C(=C1)O)N)C(=O)O"),
}

SEQUENCES = {
    "NbHCT1": ("PV344582", "XQH47715.1"),
    "NbHCT3": ("PV344584", "XQH47717.1"),
    "NbHCT8": ("PV344589", "XQH47722.1"),
}

# enzyme, table, variable, constant, Km, Km SEM, kcat, kcat SEM, saturation status
ROWS = [
    ("NbHCT1", "Table 2", "caffeoyl-CoA", "shikimic acid", 27.64, 0.73, 0.0085, 0.0001, "accepted"),
    ("NbHCT1", "Table 2", "shikimic acid", "caffeoyl-CoA", 1134.67, 153.56, 0.0114, 0.0009, "accepted"),
    ("NbHCT1", "Table 2", "p-coumaroyl-CoA", "shikimic acid", 60.13, 2.92, 0.1967, 0.0174, "accepted"),
    ("NbHCT1", "Table 2", "shikimic acid", "p-coumaroyl-CoA", 1556.00, 31.34, 0.1777, 0.0063, "excluded_cosubstrate_not_saturated"),
    ("NbHCT1", "Table 2", "p-coumaroyl-CoA", "quinic acid", 29.36, 1.73, 0.0240, 0.0030, "accepted"),
    ("NbHCT1", "Table 2", "quinic acid", "p-coumaroyl-CoA", 30146.46, 2954.84, 0.0157, 0.0003, "accepted"),
    ("NbHCT1", "Table 2", "p-coumaroyl-CoA", "3-hydroxyanthranilic acid", 4.92, 0.42, 0.0178, 0.0011, "accepted"),
    ("NbHCT1", "Table 2", "3-hydroxyanthranilic acid", "p-coumaroyl-CoA", 226.29, 138.37, 0.0140, 0.0015, "accepted"),
    ("NbHCT3", "Table 3", "caffeoyl-CoA", "shikimic acid", 168.97, "", 0.0379, "", "excluded_variable_substrate_not_saturated"),
    ("NbHCT3", "Table 3", "shikimic acid", "caffeoyl-CoA", 233975.01, "", 0.0365, "", "excluded_variable_substrate_not_saturated"),
    ("NbHCT3", "Table 3", "p-coumaroyl-CoA", "shikimic acid", 176.87, "", 0.0345, "", "excluded_variable_substrate_not_saturated"),
    ("NbHCT3", "Table 3", "shikimic acid", "p-coumaroyl-CoA", 146798.21, "", 0.0235, "", "excluded_variable_substrate_not_saturated"),
    ("NbHCT3", "Table 3", "caffeoyl-CoA", "quinic acid", 11.53, 1.31, 0.4338, 0.0299, "accepted"),
    ("NbHCT3", "Table 3", "quinic acid", "caffeoyl-CoA", 1751.67, 109.48, 0.2778, 0.0125, "accepted"),
    ("NbHCT3", "Table 3", "p-coumaroyl-CoA", "quinic acid", 34.06, 4.19, 0.6080, 0.0808, "accepted"),
    ("NbHCT3", "Table 3", "quinic acid", "p-coumaroyl-CoA", 4550.67, 2166.63, 1.1759, 0.0870, "accepted"),
    ("NbHCT8", "Table 4", "caffeoyl-CoA", "shikimic acid", 6.00, 1.03, 0.0448, 0.0012, "accepted"),
    ("NbHCT8", "Table 4", "shikimic acid", "caffeoyl-CoA", 2329.33, 153.17, 0.0855, 0.0021, "accepted"),
    ("NbHCT8", "Table 4", "p-coumaroyl-CoA", "shikimic acid", 19.07, 1.09, 0.1096, 0.0062, "accepted"),
    ("NbHCT8", "Table 4", "shikimic acid", "p-coumaroyl-CoA", 1027.33, 67.14, 0.1170, 0.0059, "accepted"),
    ("NbHCT8", "Table 4", "p-coumaroyl-CoA", "quinic acid", 10.53, 0.61, 0.0384, 0.0009, "accepted"),
    ("NbHCT8", "Table 4", "quinic acid", "p-coumaroyl-CoA", 54753.33, 13037.17, 0.0348, 0.0045, "accepted"),
    ("NbHCT8", "Table 4", "p-coumaroyl-CoA", "3-hydroxyanthranilic acid", 3.51, 0.24, 0.1085, 0.0050, "accepted"),
    ("NbHCT8", "Table 4", "3-hydroxyanthranilic acid", "p-coumaroyl-CoA", 307.60, 10.39, 0.1073, 0.0045, "accepted"),
]

CONDITIONS = {
    ("NbHCT1", "caffeoyl-CoA", "shikimic acid"): (30, 15, "1-200 uM", "10 mM", 8.8, 2.0, 360),
    ("NbHCT1", "shikimic acid", "caffeoyl-CoA"): (30, 5, "0.5-16 mM", "200 uM", 7.2, 2.0, 360),
    ("NbHCT1", "p-coumaroyl-CoA", "shikimic acid"): (30, 20, "10-400 uM", "10 mM", 6.4, 0.26, 360),
    ("NbHCT1", "shikimic acid", "p-coumaroyl-CoA"): (30, 20, "0.5-16 mM", "200 uM", 3.3, 0.26, 360),
    ("NbHCT1", "p-coumaroyl-CoA", "quinic acid"): (30, 10, "10-400 uM", "450 mM", 15, 2.0, 130),
    ("NbHCT1", "quinic acid", "p-coumaroyl-CoA"): (30, 10, "10-450 mM", "200 uM", 6.7, 2.0, 130),
    ("NbHCT1", "p-coumaroyl-CoA", "3-hydroxyanthranilic acid"): (30, 20, "5-200 uM", "12 mM", 53, 6.0, 32),
    ("NbHCT1", "3-hydroxyanthranilic acid", "p-coumaroyl-CoA"): (30, 20, "0.3-12 mM", "200 uM", 41, 6.0, 32),
    ("NbHCT3", "caffeoyl-CoA", "shikimic acid"): (35, 30, "5-200 uM", "600 mM", 2.6, 2.5, 32),
    ("NbHCT3", "shikimic acid", "caffeoyl-CoA"): (35, 30, "20-600 mM", "200 uM", 12, 2.5, 32),
    ("NbHCT3", "p-coumaroyl-CoA", "shikimic acid"): (35, 30, "5-200 uM", "600 mM", 4.1, 2.5, 32),
    ("NbHCT3", "shikimic acid", "p-coumaroyl-CoA"): (35, 60, "20-600 mM", "200 uM", 11, 2.5, 32),
    ("NbHCT3", "caffeoyl-CoA", "quinic acid"): (35, 5, "1-100 uM", "100 mM", 57, 0.2, 32),
    ("NbHCT3", "quinic acid", "caffeoyl-CoA"): (35, 5, "0.1-100 mM", "100 uM", 8.7, 0.2, 32),
    ("NbHCT3", "p-coumaroyl-CoA", "quinic acid"): (35, 5, "1-200 uM", "200 mM", 44, 0.2, 32),
    ("NbHCT3", "quinic acid", "p-coumaroyl-CoA"): (35, 5, "1-200 mM", "200 uM", 6.0, 0.2, 32),
    ("NbHCT8", "caffeoyl-CoA", "shikimic acid"): (40, 10, "1-100 uM", "20 mM", 8.6, 1.0, 32),
    ("NbHCT8", "shikimic acid", "caffeoyl-CoA"): (40, 10, "0.5-30 mM", "200 uM", 18, 1.0, 32),
    ("NbHCT8", "p-coumaroyl-CoA", "shikimic acid"): (40, 20, "2-200 uM", "8 mM", 8.0, 0.3, 32),
    ("NbHCT8", "shikimic acid", "p-coumaroyl-CoA"): (40, 10, "0.2-20 mM", "200 uM", 10, 0.3, 32),
    ("NbHCT8", "p-coumaroyl-CoA", "quinic acid"): (40, 8, "10-300 uM", "300 mM", 63, 1.0, 32),
    ("NbHCT8", "quinic acid", "p-coumaroyl-CoA"): (40, 8, "5-200 mM", "200 uM", 19, 1.0, 32),
    ("NbHCT8", "p-coumaroyl-CoA", "3-hydroxyanthranilic acid"): (40, 10, "2-100 uM", "3000 uM", 9.8, 2.0, 32),
    ("NbHCT8", "3-hydroxyanthranilic acid", "p-coumaroyl-CoA"): (40, 10, "75-3000 uM", "100 uM", 29, 2.0, 32),
}


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    key = ""
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            key = next(protein for protein in SEQUENCES if f"gene={protein[2:]}" in line)
            records[key] = ""
        else:
            records[key] += line.strip()
    return records


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_outputs() -> None:
    sequences = parse_fasta(SOURCE / "raw" / "fern-HCT-cds-aa.fasta")
    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        for enzyme, sequence in sequences.items():
            accession, protein_id = SEQUENCES[enzyme]
            handle.write(f">{enzyme.lower()} {accession} {protein_id} | native coding sequence\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    records = []
    for index, (enzyme, table, variable, constant, km, km_error, kcat, kcat_error, initial_status) in enumerate(ROWS, 1):
        temp, duration, variable_range, constant_conc, ratio, protein_ug, buffer_mm = CONDITIONS[(enzyme, variable, constant)]
        cid, smiles = SUBSTRATES[variable]
        records.append({
            "candidate_id": f"fern-hct-{index:03d}",
            "article_doi": "10.1111/tpj.70837",
            "source_table": table,
            "source_conditions": "Supporting Table S16",
            "organism": "Neoblechnum brasiliense",
            "sequence_id": enzyme.lower(),
            "sequence_accession": SEQUENCES[enzyme][0],
            "construct": f"full-length {enzyme} with 6xHis tag; native coding sequence mapped, complete vector-derived tag context not recovered",
            "variable_substrate": variable,
            "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles,
            "constant_cosubstrate": constant,
            "constant_cosubstrate_concentration": constant_conc,
            "constant_cosubstrate_multiple_of_km": ratio,
            "variable_substrate_range": variable_range,
            "km_uM": km,
            "km_error_uM": km_error,
            "kcat_s-1": kcat,
            "kcat_error_s-1": kcat_error,
            "error_type": "SEM" if kcat_error != "" else "not reported",
            "assay_pH": 7.0,
            "assay_temperature_C": temp,
            "assay_duration_min": duration,
            "protein_amount_ug": protein_ug,
            "assay_volume_uL": 125,
            "potassium_phosphate_mM": buffer_mm,
            "replicates": "3x3" if kcat_error != "" else 3,
            "status_at_normalization": initial_status,
        })

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if hit_path.exists():
        hit_ids = {line.split("\t", 1)[0] for line in hit_path.read_text(encoding="utf-8").splitlines() if line}
        for record in records:
            if record["sequence_id"] in hit_ids:
                record["status_at_normalization"] = "excluded_homology"
            elif record["status_at_normalization"] == "accepted":
                record["status_at_normalization"] = "accepted_homology_cold_pool"

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "biostudies-S-EPMC13033379",
        "article_doi": "10.1111/tpj.70837",
        "article_published": "2026-03-29",
        "article_license": "CC-BY-4.0",
        "reported_kcat_rows": len(records),
        "accepted_records": len(accepted),
        "excluded_saturation_records": len(records) - len(accepted),
        "sequence_source": "GenBank PV344582, PV344584, and PV344589",
        "kinetics_source": "Main Tables 2-4 and Supporting Table S16",
        "selection_policy": "All rows passing the frozen saturation gate; no kinetic-value-based selection and family cap not reached.",
        "saturation_rule": "Constant cosubstrate must be at least 5x its fitted Km and the variable substrate must provide an interpretable saturation fit.",
        "construct_caveat": "Author kcat uses molecular mass including a 6xHis tag; native coding sequences are mapped, but complete vector-derived tag context was not recovered.",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    cluster = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster.exists():
        hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        families = {line.split("\t", 1)[0] for line in cluster.read_text(encoding="utf-8").splitlines() if line}
        audit = {
            "audited_on": "2026-07-22",
            "source_id": "biostudies-S-EPMC13033379",
            "article_doi": "10.1111/tpj.70837",
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9",
            "candidate_records": len(records),
            "unique_sequences": len(sequences),
            "construct_sequences_sha256": sha256(fasta),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster),
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
            "accepted_unique_substrates": len({row["variable_substrate"] for row in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    write_outputs()
