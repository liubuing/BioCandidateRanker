from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "biostudies-S-EPMC11705230"
)
RAW = SOURCE / "raw"
TARGET = (
    ROOT
    / "artifacts"
    / "external"
    / "absolute-kinetics-screen"
    / "dryad-4964723"
    / "homology"
    / "unikp_reference.fasta"
)
DOI = "10.1021/acscatal.4c04935"
SOURCE_ID = "biostudies-S-EPMC11705230"
SATURATION_MULTIPLE = 5.0
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
TARGET_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"

STRUCTURES = {
    "hexanal": (6184, "CCCCCC=O", "Hexanal"),
    "allylamine": (7853, "C=CCN", "Allylamine"),
    "NADPH": (
        5884,
        "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O",
        "NADPH",
    ),
    "NADH": (
        439153,
        "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O",
        "Nadh",
    ),
}

# Main Table 2: cofactor, Km, Km SD, kcat, kcat SD, efficiency, SI range.
KINETICS = (
    ("NADPH", 15.0, 4.0, 3.6, 0.2, 241.0, "0.025-0.4 mM"),
    ("NADH", 247.0, 24.0, 9.0, 0.3, 36.0, "0.025-0.5 mM"),
)

FIELDS = [
    "source_record_id",
    "article_doi",
    "stable_record_url",
    "source_table",
    "source_figure",
    "organism",
    "enzyme",
    "sequence_id",
    "sequence_accessions",
    "kinetic_construct",
    "exact_sequence_basis",
    "reaction",
    "carbonyl_substrate",
    "carbonyl_pubchem_cid",
    "carbonyl_isomeric_smiles",
    "amine_substrate",
    "amine_pubchem_cid",
    "amine_isomeric_smiles",
    "variable_substrate",
    "substrate_pubchem_cid",
    "substrate_isomeric_smiles",
    "variable_substrate_range",
    "km_uM",
    "km_error_uM",
    "kcat_s-1",
    "kcat_error_s-1",
    "kcat_per_km_s-1_mM-1",
    "error_type",
    "assay_method",
    "assay_buffer",
    "assay_pH",
    "assay_temperature_C",
    "enzyme_amount_ug",
    "hexanal_concentration_mM",
    "hexanal_matching_km",
    "allylamine_concentration_mM",
    "allylamine_matching_km",
    "cosubstrate_saturation_rule",
    "cosubstrates_demonstrated_saturated",
    "saturation_provenance",
    "status_at_normalization",
    "exclusion_reason",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdb_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s+;([^;]+);", text, re.DOTALL
    )
    if not match:
        raise ValueError(f"Canonical polymer sequence missing from {path.name}")
    sequence = re.sub(r"\s+", "", match.group(1))
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(f"Invalid protein polymer in {path.name}")
    for marker in (
        "_struct_ref.pdbx_db_accession          A0A3G9A570",
        "_entity_src_gen.pdbx_gene_src_gene                 EHW12_26885",
        "_pdbx_struct_assembly.oligomeric_details   dimeric",
    ):
        if marker not in text:
            raise ValueError(f"Expected construct metadata absent from {path.name}: {marker}")
    return sequence


def validate_sources() -> str:
    metadata = json.loads((RAW / "metadata.json").read_text(encoding="utf-8"))
    if metadata["accno"] != "S-EPMC11705230":
        raise ValueError("Unexpected BioStudies accession")
    files = metadata["section"]["files"]
    if files != [{"path": "cs4c04935_si_001.pdf", "size": 2804916, "type": "file"}]:
        raise ValueError("Unexpected BioStudies deposited-file manifest")
    if (RAW / files[0]["path"]).stat().st_size != files[0]["size"]:
        raise ValueError("BioStudies SI size mismatch")

    with zipfile.ZipFile(RAW / "PMC11705230-supplementaryFiles.zip") as archive:
        info = archive.getinfo("cs4c04935_si_001.pdf")
        if info.file_size != files[0]["size"]:
            raise ValueError("Europe PMC supplementary archive SI size mismatch")

    xml = (RAW / "PMC11705230-fullText.xml").read_text(encoding="utf-8")
    required_xml = (
        "10.1021/acscatal.4c04935",
        "https://creativecommons.org/licenses/by/4.0/",
        "<td style=\"border:none;\" align=\"left\" colspan=\"1\" rowspan=\"1\">NADPH</td>",
        "<td style=\"border:none;\" align=\"left\" colspan=\"1\" rowspan=\"1\">NADH</td>",
        "10 mM hexanal, 100 mM allylamine",
    )
    if any(marker not in xml for marker in required_xml):
        raise ValueError("Article XML does not contain expected license/kinetics evidence")

    sequence_9fm8 = pdb_sequence(RAW / "9FM8.cif")
    sequence_9fm7 = pdb_sequence(RAW / "9FM7.cif")
    if sequence_9fm8 != sequence_9fm7 or len(sequence_9fm8) != 291:
        raise ValueError("9FM8 and 9FM7 must contain the same exact 291-aa polymer")
    if not sequence_9fm8.startswith("MDVSILGTGLMGTALAQAL") or not sequence_9fm8.endswith(
        "KGIFAQIETLSANPQSAI"
    ):
        raise ValueError("Unexpected RytRedAm sequence termini")

    for name, (expected_cid, expected_smiles, expected_title) in STRUCTURES.items():
        payload = json.loads((RAW / f"pubchem-{name.lower()}.json").read_text(encoding="utf-8"))
        item = payload["PropertyTable"]["Properties"][0]
        if (item["CID"], item["SMILES"], item["Title"]) != (
            expected_cid,
            expected_smiles,
            expected_title,
        ):
            raise ValueError(f"Unexpected PubChem structure for {name}")
    return sequence_9fm8


def build_rows() -> list[dict[str, object]]:
    hexanal_cid, hexanal_smiles, _ = STRUCTURES["hexanal"]
    allylamine_cid, allylamine_smiles, _ = STRUCTURES["allylamine"]
    rows = []
    for index, (cofactor, km, km_error, kcat, kcat_error, efficiency, range_text) in enumerate(
        KINETICS, 1
    ):
        cid, smiles, _ = STRUCTURES[cofactor]
        row = {
            "source_record_id": f"rytredam-source-{index:03d}",
            "article_doi": DOI,
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11705230",
            "source_table": "Main Table 2",
            "source_figure": "SI Figure S12" if cofactor == "NADPH" else "SI Figure S13",
            "organism": "Rhodococcus erythropolis",
            "enzyme": "RytRedAm",
            "sequence_id": "rytredam-pdb-exact",
            "sequence_accessions": "WP_020971038.1;A0A3G9A570;EHW12_26885;PDB 9FM8;PDB 9FM7",
            "kinetic_construct": "pET-28a(+) N-terminal His-tagged RytRedAm; complete vector-derived tag sequence not reported",
            "exact_sequence_basis": "291-aa native coding sequence exactly shared by SI, 9FM8, and 9FM7; crystallographic HRV3C-cleaved protein",
            "reaction": "hexanal + allylamine + NAD(P)H + H+ -> N-allylhexan-1-amine + NAD(P)+ + H2O",
            "carbonyl_substrate": "hexanal",
            "carbonyl_pubchem_cid": hexanal_cid,
            "carbonyl_isomeric_smiles": hexanal_smiles,
            "amine_substrate": "allylamine",
            "amine_pubchem_cid": allylamine_cid,
            "amine_isomeric_smiles": allylamine_smiles,
            "variable_substrate": cofactor,
            "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles,
            "variable_substrate_range": range_text,
            "km_uM": km,
            "km_error_uM": km_error,
            "kcat_s-1": kcat,
            "kcat_error_s-1": kcat_error,
            "kcat_per_km_s-1_mM-1": efficiency,
            "error_type": "one standard deviation",
            "assay_method": "direct NAD(P)H-consumption initial-rate Michaelis-Menten fit; A340, with A370 used at high NADH",
            "assay_buffer": "100 mM potassium phosphate",
            "assay_pH": 7.0,
            "assay_temperature_C": 30,
            "enzyme_amount_ug": 7.5,
            "hexanal_concentration_mM": 10,
            "hexanal_matching_km": "not reported",
            "allylamine_concentration_mM": 100,
            "allylamine_matching_km": "not reported",
            "cosubstrate_saturation_rule": ">=5x each matching fitted Km",
            "cosubstrates_demonstrated_saturated": False,
            "saturation_provenance": "SI Section 3.2 and Figures S12-S13 fix hexanal at 10 mM and allylamine at 100 mM; the article reports no fitted Km for either fixed substrate, so neither multiple can be established.",
            "status_at_normalization": "excluded_cosubstrate_not_demonstrated_saturated",
            "exclusion_reason": "Matching Km values are absent for both fixed organic cosubstrates; the frozen >=5x gate cannot be demonstrated.",
        }
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs() -> None:
    if sha256(TARGET) != TARGET_SHA256:
        raise ValueError("Frozen UniKP development target identity mismatch")
    sequence = validate_sources()
    rows = build_rows()

    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(
            ">rytredam-pdb-exact WP_020971038.1 A0A3G9A570 EHW12_26885 | "
            "exact 9FM8/9FM7 HRV3C-cleaved polymer; homodimer subunit\n"
        )
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")

    # No row passes the frozen cosubstrate gate; candidate_records remains schema-only.
    write_csv(SOURCE / "candidate_records.csv", [])
    write_csv(SOURCE / "exclusions.csv", rows)

    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir())
        if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )

    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11705230",
        "article_doi": DOI,
        "article_published": "2024-12-16",
        "license": "CC-BY-4.0",
        "reported_direct_finite_steady_state_kcat_rows": len(rows),
        "saturation_eligible_rows": 0,
        "saturation_exclusions": len(rows),
        "accepted_records": 0,
        "kinetics_sources": [
            "raw/PMC11705230-fullText.xml, Main Table 2",
            "raw/cs4c04935_si_001.pdf, Section 3.2 and Figures S12-S13",
        ],
        "construct_sources": [
            "raw/cs4c04935_si_001.pdf, Sections 1.1, 1.3, and 1.5",
            "raw/9FM8.cif",
            "raw/9FM7.cif",
        ],
        "sequence_mapping": "The SI WP_020971038.1 native sequence and the canonical polymer entities in 9FM8/9FM7 are identical 291-aa RytRedAm. Both PDB entries map to A0A3G9A570 and gene EHW12_26885 and describe a dimer.",
        "construct_boundary": "Kinetics used the pET-28a(+) N-terminal His-tagged protein without a reported full tag sequence. The exact 9FM8/9FM7 construct was produced separately in pETYSBLIC-3C and HRV3C-cleaved; its deposited polymer is the exact native 291-aa sequence. No unreported tag residues were invented.",
        "assay": {
            "method": "Direct NAD(P)H-consumption initial-rate Michaelis-Menten fit",
            "buffer": "100 mM potassium phosphate, pH 7.0",
            "temperature_C": 30,
            "enzyme_amount_ug": 7.5,
            "hexanal_mM": 10,
            "allylamine_mM": 100,
            "NADPH_range_mM": "0.025-0.4",
            "NADH_range_mM": "0.025-0.5 as stated in the SI Figure S13 caption; plotted points extend higher",
            "error_type": "one standard deviation",
        },
        "saturation_policy": "Every fixed cosubstrate must be >=5x its matching fitted Km. Concentrations alone do not establish saturation. Neither hexanal nor allylamine has a reported matching Km, so both finite cofactor-kcat rows are excluded.",
        "excluded_data": [
            "Main Table 2 NADPH row: fixed hexanal and allylamine lack matching fitted Km values.",
            "Main Table 2 NADH row: fixed hexanal and allylamine lack matching fitted Km values.",
            "Specific activities and pH/temperature profiles: not saturation-derived absolute kcat rows.",
            "Biotransformation conversions, time courses, isolated yields, and substrate-scope screens: endpoint/progress data, not direct steady-state kcat.",
        ],
        "pubchem_structure_source": "PubChem PUG REST, retrieved 2026-07-22",
        "raw_file_hashes": "raw-file-hashes.json",
        "construct_sequences_sha256": sha256(fasta),
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        cluster_lines = [
            line for line in cluster_path.read_text(encoding="utf-8").splitlines() if line
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
            "development_target_sha256": sha256(TARGET),
            "reported_direct_finite_steady_state_kcat_rows": len(rows),
            "saturation_exclusions": len(rows),
            "candidate_records": 0,
            "unique_exact_pdb_construct_sequences": 1,
            "construct_sequences_sha256": sha256(fasta),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(
                {line.split("\t", 1)[0] for line in cluster_lines}
            ),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": 0,
            "readiness_gate_passes": False,
            "claim_boundary": "Saturation audit only; no accepted labels or model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )

    verify_outputs(rows, sequence)


def verify_outputs(rows: list[dict[str, object]], sequence: str) -> None:
    if len(rows) != 2 or len(sequence) != 291:
        raise ValueError("Unexpected RytRedAm curation cardinality")
    if {row["variable_substrate"] for row in rows} != {"NADPH", "NADH"}:
        raise ValueError("Unexpected kinetic cofactors")
    for row in rows:
        for field in ("km_uM", "km_error_uM", "kcat_s-1", "kcat_error_s-1"):
            value = float(row[field])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"Invalid {field} in {row['source_record_id']}")
        if row["cosubstrates_demonstrated_saturated"] is not False:
            raise ValueError("Frozen saturation gate was not applied strictly")
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        if list(csv.DictReader(handle)):
            raise ValueError("No accepted candidate record should have been emitted")
    with (SOURCE / "exclusions.csv").open(newline="", encoding="utf-8") as handle:
        if len(list(csv.DictReader(handle))) != 2:
            raise ValueError("Expected two deterministic exclusion rows")
    required = [
        SOURCE / "construct_sequences.fasta",
        SOURCE / "candidate_records.csv",
        SOURCE / "exclusions.csv",
        SOURCE / "raw-file-hashes.json",
        SOURCE / "provenance.json",
    ]
    if not all(path.is_file() and path.stat().st_size for path in required):
        raise ValueError("A required standard output is missing or empty")


if __name__ == "__main__":
    write_outputs()
