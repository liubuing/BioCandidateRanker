from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "biostudies-S-EPMC12364357"
)
RAW = SOURCE / "raw"
MOLECULAR_MASS_KDA = 46.5
ENZYME_AMOUNT_UG = 1.11
ASSAY_VOLUME_ML = 0.5
SATURATION_MULTIPLE = 5.0

STRUCTURES = {
    "NADH": (
        439153,
        "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O",
    ),
    "NADPH": (
        5884,
        "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O",
    ),
    "2-oxoglutarate": (51, "C(CC(=O)O)C(=O)C(=O)O"),
    "ammonium chloride": (25517, "[NH4+].[Cl-]"),
    "NAD+": (
        5892,
        "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
    ),
    "NADP+": (
        5885,
        "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)N",
    ),
    "L-glutamate": (33032, "C(CC(=O)O)[C@@H](C(=O)O)N"),
}

# direction, cofactor, variable, Km, Km SEM, Vmax, Vmax SEM, reported kcat
ROWS = [
    ("reductive amination", "NADH", "NADH", 0.070, 0.014, 5.43, 0.34, 4.21),
    ("reductive amination", "NADPH", "NADPH", 0.078, 0.017, 6.29, 0.46, 4.87),
    ("reductive amination", "NADH", "2-oxoglutarate", 1.30, 0.21, 5.93, 0.26, 4.60),
    ("reductive amination", "NADPH", "2-oxoglutarate", 0.98, 0.24, 4.16, 0.33, 3.22),
    ("reductive amination", "NADH", "ammonium chloride", 5.87, 0.20, 5.01, 0.04, 3.88),
    ("reductive amination", "NADPH", "ammonium chloride", 12.88, 3.23, 6.05, 0.48, 4.69),
    ("oxidative deamination", "NAD+", "NAD+", 0.90, 0.06, 11.57, 0.30, 8.97),
    ("oxidative deamination", "NADP+", "NADP+", 0.44, 0.05, 0.67, 0.02, 0.52),
    ("oxidative deamination", "NAD+", "L-glutamate", 99.1, 15.9, 10.48, 1.24, 8.12),
    ("oxidative deamination", "NADP+", "L-glutamate", 27.6, 5.4, 0.50, 0.08, 0.39),
]

KM = {
    (direction, cofactor, variable): km
    for direction, cofactor, variable, km, *_ in ROWS
}


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fixed_cosubstrates(direction: str, cofactor: str, variable: str) -> list[dict[str, object]]:
    if direction == "reductive amination":
        fixed = [("2-oxoglutarate", 5.0), ("ammonium chloride", 50.0), (cofactor, 0.2)]
    else:
        fixed = [("L-glutamate", 150.0), (cofactor, 2.0)]
    evidence = []
    for name, concentration_mm in fixed:
        if name == variable:
            continue
        km_mm = KM[(direction, cofactor, name)]
        evidence.append(
            {
                "name": name,
                "concentration_mM": concentration_mm,
                "measured_km_mM": km_mm,
                "multiple_of_km": round(concentration_mm / km_mm, 6),
            }
        )
    return evidence


def source_rows() -> list[dict[str, object]]:
    rows = []
    for index, row in enumerate(ROWS, 1):
        direction, cofactor, variable, km, km_sem, vmax, vmax_sem, reported_kcat = row
        converted_kcat = vmax * MOLECULAR_MASS_KDA / 60
        if round(converted_kcat, 2) != reported_kcat:
            raise ValueError(
                f"Table 2 row {index} kcat mismatch: {converted_kcat} versus {reported_kcat}"
            )
        evidence = fixed_cosubstrates(direction, cofactor, variable)
        saturated = all(item["multiple_of_km"] >= SATURATION_MULTIPLE for item in evidence)
        cid, smiles = STRUCTURES[variable]
        variable_range = {
            "2-oxoglutarate": "0-15",
            "NADH": "0-0.5",
            "NADPH": "0-0.5",
            "ammonium chloride": "0-100",
            "NAD+": "0-5",
            "NADP+": "0-5",
            "L-glutamate": "0-250",
        }[variable]
        if direction == "reductive amination":
            buffer, ph = "80 mM Tris-HCl", 8.0
            reaction = "2-oxoglutarate + NH4+ + NAD(P)H + H+ -> L-glutamate + NAD(P)+ + H2O"
        else:
            buffer = "80 mM CHES" if cofactor == "NAD+" else "80 mM imidazole"
            ph = 9.5 if cofactor == "NAD+" else 6.0
            reaction = "L-glutamate + NAD(P)+ + H2O -> 2-oxoglutarate + NH4+ + NAD(P)H + H+"
        rows.append(
            {
                "source_row": index,
                "direction": direction,
                "reaction": reaction,
                "cofactor_system": cofactor,
                "variable_substrate": variable,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "km_mM": km,
                "km_error_mM": km_sem,
                "source_vmax_U_mg-1": vmax,
                "source_vmax_error_U_mg-1": vmax_sem,
                "source_reported_kcat_s-1": reported_kcat,
                "derived_kcat_error_s-1": round(vmax_sem * MOLECULAR_MASS_KDA / 60, 6),
                "cosubstrate_saturation_evidence": json.dumps(evidence, separators=(",", ":")),
                "cosubstrates_demonstrated_saturated": saturated,
                "variable_substrate_range_mM": variable_range,
                "assay_buffer": buffer,
                "assay_pH": ph,
                "assay_temperature_C": 25,
                "assay_volume_mL": ASSAY_VOLUME_ML,
                "enzyme_amount_ug": ENZYME_AMOUNT_UG,
                "enzyme_concentration_nM": round(
                    ENZYME_AMOUNT_UG * 1e3 / (MOLECULAR_MASS_KDA * ASSAY_VOLUME_ML), 6
                ),
                "concentrations_tested": "at least 8",
                "error_type": "SEM",
                "replicates": 3,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    names = fieldnames or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def homology_hit_ids() -> set[str]:
    path = SOURCE / "homology" / "homology_hits.tsv"
    if not path.exists():
        return set()
    return {
        line.split("\t", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    }


def write_outputs() -> None:
    metadata = json.loads((RAW / "metadata.json").read_text(encoding="utf-8"))
    deposited_files = metadata["section"]["files"]
    for deposited in deposited_files:
        path = RAW / deposited["path"]
        if not path.is_file() or path.stat().st_size != deposited["size"]:
            raise ValueError(f"BioStudies acquisition mismatch for {path.name}")

    bak = read_fasta(RAW / "BAK13910.1.fasta")
    wp = read_fasta(RAW / "WP_013027799.1.fasta")
    if bak != wp or len(bak) != 424:
        raise ValueError("BAK13910.1 and WP_013027799.1 must be identical 424-aa sequences")
    if not bak.startswith("MDKLSYASDSSTSAWSTYL") or not bak.endswith("MARKDRGIYPG"):
        raise ValueError("Unexpected GdhPa accession sequence termini")

    # TEV cleaves ENLYFQ|G. The S1 primer junction is ENLYFQG followed by native Met.
    construct = "G" + bak
    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(
            ">gdhpa-tev BAK13910.1=WP_013027799.1 | pET-15-TEV, cleaved G+native construct\n"
        )
        for start in range(0, len(construct), 80):
            handle.write(construct[start : start + 80] + "\n")

    extracted = source_rows()
    eligible = [row for row in extracted if row["cosubstrates_demonstrated_saturated"]]
    excluded = [row for row in extracted if not row["cosubstrates_demonstrated_saturated"]]
    if len(eligible) != 1 or eligible[0]["variable_substrate"] != "NADP+":
        raise ValueError("Expected only the NADP+ oxidative-deamination row to pass saturation")

    hits = homology_hit_ids()
    candidate = eligible[0]
    records = [
        {
            "candidate_id": "pantoea-gdh-001",
            "article_doi": "10.1371/journal.pone.0328289",
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12364357",
            "source_table": "Table 2",
            "source_row": "oxidative deamination; NADP+-dependent; NADP+ varied",
            "organism": "Pantoea ananatis AJ13355",
            "sequence_id": "gdhpa-tev",
            "sequence_accessions": "BAK13910.1;WP_013027799.1",
            "construct": "pET-15-TEV-GdhA(Pa), affinity tag removed; G+BAK13910.1",
            "direction": candidate["direction"],
            "reaction": "L-glutamate + NADP+ + H2O -> 2-oxoglutarate + NH4+ + NADPH + H+",
            "variable_substrate": candidate["variable_substrate"],
            "substrate_pubchem_cid": candidate["substrate_pubchem_cid"],
            "substrate_isomeric_smiles": candidate["substrate_isomeric_smiles"],
            "fixed_cosubstrates": candidate["cosubstrate_saturation_evidence"],
            "cosubstrate_saturation_rule": ">=5x each measured Km",
            "cosubstrates_demonstrated_saturated": True,
            "variable_substrate_range_mM": "0-5",
            "km_uM": float(candidate["km_mM"]) * 1000,
            "km_error_uM": float(candidate["km_error_mM"]) * 1000,
            "source_vmax_U_mg-1": candidate["source_vmax_U_mg-1"],
            "source_vmax_error_U_mg-1": candidate["source_vmax_error_U_mg-1"],
            "kcat_s-1": candidate["source_reported_kcat_s-1"],
            "kcat_error_s-1": candidate["derived_kcat_error_s-1"],
            "kcat_error_derivation": "Vmax_SEM*46.5_kDa/60",
            "molecular_mass_kDa": MOLECULAR_MASS_KDA,
            "enzyme_amount_ug": ENZYME_AMOUNT_UG,
            "enzyme_concentration_nM": round(
                ENZYME_AMOUNT_UG * 1e3 / (MOLECULAR_MASS_KDA * ASSAY_VOLUME_ML), 6
            ),
            "assay_volume_mL": ASSAY_VOLUME_ML,
            "assay_buffer": "80 mM imidazole",
            "assay_pH": 6.0,
            "assay_temperature_C": 25,
            "assay_method": "A340 initial-rate assay; epsilon=6.22e3 M-1 cm-1",
            "concentrations_tested": "at least 8",
            "error_type": "SEM",
            "replicates": 3,
            "status_at_normalization": (
                "excluded_homology" if "gdhpa-tev" in hits else "accepted_homology_cold_pool"
            )
            if (SOURCE / "homology" / "homology_hits.tsv").exists()
            else "pending_homology",
        }
    ]
    write_csv(SOURCE / "candidate_records.csv", records)

    exclusion_rows = []
    for row in excluded:
        failed = [
            item
            for item in json.loads(str(row["cosubstrate_saturation_evidence"]))
            if item["multiple_of_km"] < SATURATION_MULTIPLE
        ]
        exclusion_rows.append(
            {
                **row,
                "exclusion_reason": "cosubstrate_not_demonstrated_saturated",
                "failed_cosubstrates": json.dumps(failed, separators=(",", ":")),
                "candidate_label_created": False,
            }
        )
    write_csv(SOURCE / "exclusions.csv", exclusion_rows)

    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir())
        if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )

    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": "biostudies-S-EPMC12364357",
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12364357",
        "article_doi": "10.1371/journal.pone.0328289",
        "article_published": "2025-08-19",
        "license": "CC-BY-4.0",
        "kinetics_source": "raw/PMC12364357-fullText.xml, Table 2 and Enzymatic assays",
        "construct_sources": [
            "raw/PMC12364357-fullText.xml, expression/purification method",
            "raw/pone.0328289.s006.docx, S1 Table primers",
            "raw/BAK13910.1.fasta",
            "raw/WP_013027799.1.fasta",
        ],
        "sequence_recovery": "BAK13910.1 and WP_013027799.1 are identical (424 aa). The S1 primer junction encodes ENLYFQG immediately before native Met; TEV cleavage at Q|G yields G+native (425 aa).",
        "reported_table_rows": len(extracted),
        "saturation_eligible_rows": len(records),
        "saturation_exclusions": len(excluded),
        "accepted_records": len(accepted),
        "conversion": {
            "formula": "kcat_s-1 = Vmax_U_mg-1 * subunit_molecular_mass_kDa / 60",
            "subunit_molecular_mass_kDa": MOLECULAR_MASS_KDA,
            "enzyme_amount_ug": ENZYME_AMOUNT_UG,
            "assay_volume_mL": ASSAY_VOLUME_ML,
            "enzyme_concentration_nM": records[0]["enzyme_concentration_nM"],
            "central_values": "Author-reported kcat retained after every row reproduced at two decimals.",
            "errors": "kcat SEM derived as Vmax SEM * 46.5 / 60.",
        },
        "saturation_policy": "Every fixed cosubstrate must be at least 5x its measured Km in the matching cofactor/direction assay. Only oxidative NADP+ variation passes (150 mM L-glutamate / 27.6 mM Km = 5.434783).",
        "excluded_data": "Nine Table 2 rows failing demonstrated cosubstrate saturation; specific-activity comparisons, pH profiles, effectors, crude-extract activities, and in-silico results.",
        "pubchem_structure_source": "PubChem PUG REST, retrieved 2026-07-22",
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
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        audit = {
            "audited_on": "2026-07-22",
            "source_id": "biostudies-S-EPMC12364357",
            "article_doi": "10.1371/journal.pone.0328289",
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9",
            "candidate_records": len(records),
            "source_table_rows": len(extracted),
            "saturation_exclusions": len(excluded),
            "unique_sequences": 1,
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({row[0] for row in hit_rows}),
            "exact_sequence_overlap": sum(float(row[2]) >= 0.999 for row in hit_rows),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
            "accepted_unique_substrates": len({row["variable_substrate"] for row in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )


if __name__ == "__main__":
    if not math.isfinite(MOLECULAR_MASS_KDA) or MOLECULAR_MASS_KDA <= 0:
        raise ValueError("A finite positive molecular mass is required")
    write_outputs()
