"""Curate Acinetobacter baumannii ATP phosphoribosyltransferase kinetics from DOI
10.1021/acs.biochem.3c00551.

Evidence written under artifacts/external/temporal-absolute-kinetics/biostudies-S-EPMC10795190:
candidate_records.csv (nine tabulated kcat rows for HisGS and the HisGS4:HisZ4 ATPPRT
holoenzyme), provenance.json, construct_sequences.fasta, and homology-audit.json.

A row is admitted only when at least one saturation curve held its fixed cosubstrate at
>=5x its fitted Km; surviving rows must then be homology-cold against the frozen
development reference and stay within the per-family cap of 20 records.

中文：整理 DOI 10.1021/acs.biochem.3c00551 的 Acinetobacter baumannii ATP 磷酸核糖转移酶
（ATPPRT）动力学。接受/排除逻辑要点：仅当至少一条饱和曲线把固定共底物维持在 >=5 倍其
拟合 Km 时该行被接受；幸存行还须对冻结开发参考库同源冷（零命中），并遵守每家族 20 条
的上限。证据写入 artifacts/external/temporal-absolute-kinetics/biostudies-S-EPMC10795190：
candidate_records.csv（HisGS 与 HisGS4:HisZ4 ATPPRT 全酶的 9 条表列 kcat 行）、
provenance.json、construct_sequences.fasta 与 homology-audit.json。
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "biostudies-S-EPMC10795190"
FAMILY_CAP = 20
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"

HISGS_ID = "ab-hisgs-wp001000724"
HISZ_ID = "ab-hisz-wp000155680"
SEQUENCE_SOURCES = {
    HISGS_ID: ("WP_001000724.1", "raw/WP_001000724.1.fasta", 227),
    HISZ_ID: ("WP_000155680.1", "raw/WP_000155680.1.fasta", 388),
}

SUBSTRATES = {
    "ATP": (
        5957,
        "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
    ),
    "ADP": (
        6022,
        "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
    ),
    "PRPP": (
        7339,
        "C([C@@H]1[C@H]([C@H]([C@H](O1)OP(=O)(O)OP(=O)(O)O)O)O)OP(=O)(O)O",
    ),
}

# table, temperature, state, nucleotide, metal, kcat, error, Km nucleotide,
# error, Km PRPP, error, nucleotide efficiency, error, PRPP efficiency, error
ROWS = [
    ("Table 1", 25, "HisGS", "ATP", "Mg2+", 0.384, 0.006, 0.83, 0.06, 0.60, 0.06, 460, 30, 640, 60),
    ("Table 1", 25, "HisGS", "ADP", "Mg2+", 0.48, 0.02, 1.5, 0.3, 1.2, 0.1, 320, 70, 400, 40),
    ("Table 1", 25, "HisGS", "ATP", "Mn2+", 0.94, 0.03, 0.39, 0.07, 0.60, 0.08, 2400, 400, 1600, 200),
    ("Table 1", 25, "HisGS", "ADP", "Mn2+", 3.3, 0.1, 2.2, 0.3, 0.44, 0.06, 1500, 200, 8000, 1000),
    ("Table 1", 25, "ATPPRT", "ATP", "Mg2+", 10.8, 0.3, 0.19, 0.02, 0.14, 0.01, 57000, 6000, 77000, 6000),
    ("Table 1", 25, "ATPPRT", "ADP", "Mg2+", 16.6, 0.3, 0.36, 0.03, 0.096, 0.007, 46000, 4000, 170000, 10000),
    ("Supporting Table S3", 5, "HisGS", "ATP", "Mg2+", 0.0413, 0.0006, 1.07, 0.06, 0.34, 0.03, 39, 2, 120, 10),
    ("Supporting Table S3", 5, "HisGS", "ATP", "Mn2+", 0.18, 0.01, 0.49, 0.06, 1.4, 0.2, 370, 50, 130, 20),
    ("Supporting Table S3", 5, "ATPPRT", "ATP", "Mg2+", 2.95, 0.05, 0.23, 0.02, 0.111, 0.008, 13000, 1000, 27000, 2000),
]


def sha256(path: Path) -> str:
    """Return the SHA-256 hex digest of the file at ``path``.

    中文：返回 ``path`` 文件的 SHA-256 十六进制摘要。
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path, expected_accession: str, expected_length: int) -> str:
    """Read a FASTA, failing unless the header carries ``expected_accession`` and the sequence has
    the exact length.

    中文：读取 FASTA，除非头部含 ``expected_accession`` 且序列长度精确一致，否则失败。
    """
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or expected_accession not in lines[0]:
        raise ValueError(f"Expected {expected_accession} in {path.name} header")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if len(sequence) != expected_length or not sequence.isalpha() or not sequence.isupper():
        raise ValueError(f"Invalid {expected_accession} sequence in {path}")
    return sequence


def load_sequences() -> dict[str, tuple[str, str]]:
    """Load and validate the exact native RefSeq chains; return sequence_id -> (accession,
    sequence).

    中文：加载并校验精确的天然 RefSeq 链；返回 sequence_id -> (accession, sequence)。
    """
    sequences = {}
    for sequence_id, (accession, relative_path, length) in SEQUENCE_SOURCES.items():
        sequences[sequence_id] = (
            accession,
            read_fasta(SOURCE / relative_path, accession, length),
        )
    return sequences


def verify_raw_evidence() -> None:
    """Fail closed unless the article XML, supplement PDF, archive, and PubChem records carry the
    tabulated values.

    中文：除非文章 XML、补充 PDF、补充材料压缩包与 PubChem 记录都带有表列数值，否则
    fail closed。
    """
    raw = SOURCE / "raw"
    article_root = ET.parse(raw / "PMC10795190-fullText.xml").getroot()
    doi = article_root.findtext(".//article-id[@pub-id-type='doi']")
    if doi != "10.1021/acs.biochem.3c00551":
        raise ValueError(f"Unexpected article DOI: {doi}")
    table_1 = article_root.find(".//table-wrap[@id='tbl1']")
    table_1_text = " ".join(table_1.itertext()) if table_1 is not None else ""
    for value in ("0.384", "0.48", "0.94", "3.3", "10.8", "16.6"):
        if value not in table_1_text:
            raise ValueError(f"Table 1 value {value} not found in source XML")

    supplement = raw / "bi3c00551_si_001.pdf"
    if supplement.read_bytes()[:5] != b"%PDF-":
        raise ValueError("Supporting information is not a PDF")
    with fitz.open(supplement) as document:
        supplement_text = "\n".join(page.get_text() for page in document)
    for value in ("0.0413", "0.18", "2.95"):
        if value not in supplement_text:
            raise ValueError(f"Supporting Table S3 value {value} not found in source PDF")

    with zipfile.ZipFile(raw / "PMC10795190-supplementaryFiles.zip") as archive:
        if not any(name.endswith("bi3c00551_si_001.pdf") for name in archive.namelist()):
            raise ValueError("Europe PMC supplementary archive lacks the source PDF")

    pubchem = json.loads((raw / "pubchem-substrates.json").read_text(encoding="utf-8"))
    observed = {
        item["CID"]: item.get("SMILES", item.get("IsomericSMILES"))
        for item in pubchem["PropertyTable"]["Properties"]
    }
    expected = {cid: smiles for cid, smiles in SUBSTRATES.values()}
    if observed != expected:
        raise ValueError(f"PubChem structure mismatch: {observed}")


def assay_limits(state: str) -> tuple[float, float]:
    """Return the fixed cosubstrate concentrations (mM) used by the paper for (PRPP, nucleotide)
    curves.

    中文：返回论文中（PRPP 曲线、核苷酸曲线）使用的固定共底物浓度（mM）。
    """
    if state == "ATPPRT":
        return 1.6, 1.6
    return 3.2, 6.4


def build_records() -> list[dict[str, object]]:
    """Turn the frozen ROWS table into candidate records, flagging saturation-qualifying curves per
    row.

    中文：把冻结的 ROWS 表转换为候选记录，并逐行标记符合饱和判定的曲线。
    """
    records = []
    for index, row in enumerate(ROWS, 1):
        (
            table,
            temperature,
            state,
            nucleotide,
            metal,
            kcat,
            kcat_error,
            km_nucleotide,
            km_nucleotide_error,
            km_prpp,
            km_prpp_error,
            efficiency_nucleotide,
            efficiency_nucleotide_error,
            efficiency_prpp,
            efficiency_prpp_error,
        ) = row
        fixed_prpp, fixed_nucleotide = assay_limits(state)
        nucleotide_curve_fixed_multiple = fixed_prpp / km_prpp
        prpp_curve_fixed_multiple = fixed_nucleotide / km_nucleotide
        qualifying_curves = []
        if nucleotide_curve_fixed_multiple >= 5:
            qualifying_curves.append(f"{nucleotide} varied")
        if prpp_curve_fixed_multiple >= 5:
            qualifying_curves.append("PRPP varied")

        holoenzyme = state == "ATPPRT"
        components = f"{HISGS_ID};{HISZ_ID}" if holoenzyme else HISGS_ID
        component_accessions = (
            "WP_001000724.1;WP_000155680.1" if holoenzyme else "WP_001000724.1"
        )
        product = "PRATP" if nucleotide == "ATP" else "PRADP"
        replicates = 3 if temperature == 5 and (holoenzyme or metal == "Mn2+") else 2
        nucleotide_cid, nucleotide_smiles = SUBSTRATES[nucleotide]
        prpp_cid, prpp_smiles = SUBSTRATES["PRPP"]
        records.append(
            {
                "candidate_id": f"atpprt-{index:03d}",
                "article_doi": "10.1021/acs.biochem.3c00551",
                "source_table": table,
                "source_row": f"Ab{state}; {nucleotide}/{metal}; {temperature} C",
                "organism": "Acinetobacter baumannii ATCC 17978",
                "sequence_id": HISGS_ID,
                "component_sequence_ids": components,
                "component_accessions": component_accessions,
                "catalytic_state": "HisGS homodimer" if not holoenzyme else "ATPPRT hetero-octamer",
                "component_stoichiometry": "HisGS2" if not holoenzyme else "HisGS4:HisZ4",
                "construct": "tag-cleaved recombinant protein; one residual N-terminal Gly per assayed chain; FASTA stores exact native RefSeq chains",
                "reaction_direction": "forward phosphoribosyl transfer",
                "reaction_substrate_1": nucleotide,
                "substrate_1_pubchem_cid": nucleotide_cid,
                "substrate_1_isomeric_smiles": nucleotide_smiles,
                "reaction_substrate_2": "PRPP",
                "substrate_2_pubchem_cid": prpp_cid,
                "substrate_2_isomeric_smiles": prpp_smiles,
                "reaction_products": f"{product} + pyrophosphate",
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_nucleotide_mM": km_nucleotide,
                "km_nucleotide_error_mM": km_nucleotide_error,
                "km_prpp_mM": km_prpp,
                "km_prpp_error_mM": km_prpp_error,
                "kcat_per_km_nucleotide_M-1_s-1": efficiency_nucleotide,
                "kcat_per_km_nucleotide_error": efficiency_nucleotide_error,
                "kcat_per_km_prpp_M-1_s-1": efficiency_prpp,
                "kcat_per_km_prpp_error": efficiency_prpp_error,
                "assay_pH": 8.5,
                "assay_temperature_C": temperature,
                "metal": metal,
                "metal_concentration_mM": 15,
                "hisgs_concentration_uM": 0.04 if holoenzyme else (3 if temperature == 5 else 1),
                "hisz_concentration_uM": 2 if holoenzyme else "",
                "calculated_atpprt_concentration_uM": 0.039 if holoenzyme else "",
                "fixed_prpp_mM_for_nucleotide_curve": fixed_prpp,
                "fixed_prpp_multiple_of_km_prpp": round(nucleotide_curve_fixed_multiple, 3),
                "fixed_nucleotide_mM_for_prpp_curve": fixed_nucleotide,
                "fixed_nucleotide_multiple_of_km_nucleotide": round(prpp_curve_fixed_multiple, 3),
                "qualifying_saturation_curves": ";".join(qualifying_curves),
                "independent_measurements": replicates,
                "error_type": "nonlinear-regression fitting error",
                "status_at_normalization": "accepted" if qualifying_curves else "excluded_cosubstrate_not_saturated",
            }
        )
    return records


def component_hits(path: Path) -> set[str]:
    """Return the distinct query names (first TSV column) with a frozen-development homology hit.

    中文：返回在冻结开发库中有同源命中的不同查询名（TSV 第一列）。
    """
    if not path.exists():
        return set()
    return {
        line.split("\t", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    }


def apply_selection(records: list[dict[str, object]]) -> None:
    """Re-status accepted rows fail-closed: homology hits exclude, then a per-catalytic-family cap
    of FAMILY_CAP applies.

    中文：以 fail-closed 方式重定已接受行的状态：同源命中即排除，随后应用每催化家族
    FAMILY_CAP 上限。
    """
    hits = component_hits(SOURCE / "homology" / "homology_hits.tsv")
    family_counts: dict[str, int] = {}
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    component_families = {}
    if cluster_path.exists():
        for line in cluster_path.read_text(encoding="utf-8").splitlines():
            if line:
                representative, member = line.split("\t")[:2]
                component_families[member] = representative

    for record in records:
        if record["status_at_normalization"] != "accepted":
            continue
        components = str(record["component_sequence_ids"]).split(";")
        if any(component in hits for component in components):
            record["status_at_normalization"] = "excluded_homology"
            continue
        catalytic_family = component_families.get(str(record["sequence_id"]), str(record["sequence_id"]))
        used = family_counts.get(catalytic_family, 0)
        if used >= FAMILY_CAP:
            record["status_at_normalization"] = "excluded_family_cap"
        else:
            record["status_at_normalization"] = "accepted_homology_cold_pool"
            family_counts[catalytic_family] = used + 1


def write_outputs() -> None:
    """Verify raw evidence, build and select records, then write the CSV/FASTA/provenance/audit
    evidence bundle.

    中文：校验原始证据、构建并筛选记录，然后写出 CSV/FASTA/provenance/audit 证据包。
    """
    verify_raw_evidence()
    sequences = load_sequences()
    records = build_records()
    apply_selection(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (accession, sequence) in sequences.items():
            role = "catalytic HisGS" if sequence_id == HISGS_ID else "regulatory HisZ"
            handle.write(f">{sequence_id} {accession} | exact native RefSeq chain | {role}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    raw_files = sorted(path for path in (SOURCE / "raw").iterdir() if path.is_file())
    provenance = {
        "source_id": "biostudies-S-EPMC10795190",
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC10795190",
        "article_doi": "10.1021/acs.biochem.3c00551",
        "article_published": "2023-12-27",
        "article_license": "CC-BY-4.0",
        "reported_direct_finite_steady_state_kcat_rows": len(records),
        "accepted_records": len(accepted),
        "excluded_saturation_records": sum(
            record["status_at_normalization"] == "excluded_cosubstrate_not_saturated"
            for record in records
        ),
        "family_cap": FAMILY_CAP,
        "accepted_substrates": sorted(
            {str(record["reaction_substrate_1"]) for record in accepted} | {"PRPP"}
        ) if accepted else [],
        "kinetics_sources": [
            "raw/PMC10795190-fullText.xml, Table 1",
            "raw/bi3c00551_si_001.pdf, Supporting Table S3",
        ],
        "sequence_sources": [
            "raw/WP_001000724.1.fasta",
            "raw/WP_000155680.1.fasta",
            "raw/8OY0.cif (assayed-construct residual-G evidence)",
        ],
        "structure_source": "raw/pubchem-substrates.json; PubChem PUG REST isomeric SMILES for CID 5957 (ATP), 6022 (ADP), and 7339 (PRPP), retrieved 2026-07-22",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "selection_policy": "All direct finite steady-state absolute kcat rows with at least one saturation curve whose fixed cosubstrate was >=5x its fitted Km; no kinetic-value-based selection.",
        "saturation_interpretation": "A table row is reaction-level and is not duplicated across nucleotide- and PRPP-varied curves. It passes when either independently fitted curve has fixed cosubstrate >=5x Km; all tabulated curves have finite fitted parameters and errors.",
        "direction_audit": "All retained rows are forward ATP/ADP + PRPP phosphoribosyl transfer. ATP produces PRATP; alternative substrate ADP produces PRADP. No reverse-direction or equilibrium values are included.",
        "complex_audit": "HisGS is the catalytic homodimer. ATPPRT is the HisGS4:HisZ4 hetero-octamer assayed with 2 uM HisZ, reported as saturated; kcat is normalized to bound HisGS/ATPPRT active sites.",
        "construct_audit": "PDB 8OY0 reports an N-terminal Gly remaining after purification-tag cleavage on both recombinant chains. FASTA deliberately stores the exact requested native RefSeq sequences and records the assay residual Gly as construct metadata.",
        "exclusions": {
            "HisGS_ADP_Mg2+_25_C": "Neither saturation experiment had fixed cosubstrate >=5x its fitted Km (PRPP 3.2/1.2=2.667; ADP 6.4/1.5=4.267).",
            "pre_steady_state": "Burst, lag, single-turnover, and apparent approach-to-steady-state constants are not direct steady-state kcat rows.",
            "pH_profile": "Figure-only fitted/graphical kcat values are not directly tabulated exact rows.",
            "viscosity": "Apparent rates under glycerol/PEG conditions are not tabulated absolute kcat rows.",
        },
        "normalization": "Author-reported kcat in s-1; no refitting, graph digitization, or unit conversion.",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        accepted_catalytic_families = {
            record["sequence_id"] for record in accepted
        }
        audit = {
            "audited_on": "2026-07-22",
            "source_id": "biostudies-S-EPMC10795190",
            "article_doi": "10.1021/acs.biochem.3c00551",
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": len(records),
            "unique_component_sequences": len(sequences),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_component_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
            "exact_sequence_overlap": sum(float(line.split("\t")[2]) == 1.0 for line in hit_lines),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_component_mmseqs_families": len(families),
            "accepted_catalytic_families": len(accepted_catalytic_families),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": len(accepted),
            "accepted_unique_catalytic_sequences": len({record["sequence_id"] for record in accepted}),
            "accepted_unique_substrates": len(
                {str(record["reaction_substrate_1"]) for record in accepted} | ({"PRPP"} if accepted else set())
            ),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )

    verify_outputs(records, sequences)


def verify_outputs(
    records: list[dict[str, object]], sequences: dict[str, tuple[str, str]]
) -> None:
    """Fail closed on cardinality, finite parameters, the 8-accepted/1-excluded split, and required
    output files.

    中文：对记录基数、参数有限性、8 条接受/1 条排除的拆分以及必需输出文件做
    fail-closed 校验。
    """
    if len(records) != 9 or len(sequences) != 2:
        raise ValueError("Unexpected curation cardinality")
    for record in records:
        for field in ("kcat_s-1", "kcat_error_s-1", "km_nucleotide_mM", "km_prpp_mM"):
            value = float(record[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Non-finite or negative {field} in {record['candidate_id']}")
    statuses = [str(record["status_at_normalization"]) for record in records]
    if statuses.count("accepted_homology_cold_pool") != 8:
        raise ValueError(f"Expected eight accepted records, observed {statuses}")
    if statuses.count("excluded_cosubstrate_not_saturated") != 1:
        raise ValueError("Expected one saturation exclusion")
    if component_hits(SOURCE / "homology" / "homology_hits.tsv"):
        raise ValueError("Expected zero frozen-development homology hits")
    required = [
        SOURCE / "candidate_records.csv",
        SOURCE / "construct_sequences.fasta",
        SOURCE / "provenance.json",
        SOURCE / "homology-audit.json",
    ]
    if not all(path.is_file() and path.stat().st_size for path in required):
        raise ValueError("A required standard output is missing or empty")


if __name__ == "__main__":
    write_outputs()
