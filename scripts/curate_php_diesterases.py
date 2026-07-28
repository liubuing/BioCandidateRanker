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
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "europepmc-PMC11656708"
)
RAW = SOURCE / "raw"
DOI = "10.1021/acs.biochem.4c00350"
SOURCE_ID = "europepmc-PMC11656708"
FAMILY_CAP = 20
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

ACCESSIONS = {
    "BAX25588.1": (441, "Sphingobium sp. TCM1", "Sb-PDE"),
    "QCI95559.1": (441, "Novosphingobium sp. EMRT-2", "No-PDE"),
    "CCW15874.1": (442, "Sphingobium indicum BiD32", "Si-PDE"),
    "ABR41129.1": (363, "Phocaeicola vulgatus ATCC 8482", "3E38-PHP"),
}

SUBSTRATES = {
    1: (
        "bis(4-nitrophenyl) phosphate",
        255,
        "C1=CC(=CC=C1[N+](=O)[O-])OP(=O)(O)OC2=CC=C(C=C2)[N+](=O)[O-]",
        "diester",
        32,
    ),
    2: (
        "diphenyl phosphate",
        13282,
        "C1=CC=C(C=C1)OP(=O)(O)OC2=CC=CC=C2",
        "diester",
        18,
    ),
    11: (
        "4-nitrophenyl phosphate",
        378,
        "C1=CC(=CC=C1[N+](=O)[O-])OP(=O)(O)O",
        "monoester",
        32,
    ),
    12: (
        "phenyl phosphate",
        12793,
        "C1=CC=C(C=C1)OP(=O)(O)O",
        "monoester",
        18,
    ),
}

MUTATIONS = {
    "WT": (),
    "I96Y": ((96, "I", "Y"),),
    "I96A": ((96, "I", "A"),),
    "T137A": ((137, "T", "A"),),
    "I206W": ((206, "I", "W"),),
    "I206A": ((206, "I", "A"),),
    "S209Q": ((209, "S", "Q"),),
    "S209A": ((209, "S", "A"),),
    "F239L": ((239, "F", "L"),),
    "F239A": ((239, "F", "A"),),
    "L263Q": ((263, "L", "Q"),),
    "L263A": ((263, "L", "A"),),
    "E268D": ((268, "E", "D"),),
    "E268A": ((268, "E", "A"),),
}

ZINC_PER_ENZYME = {
    "WT": 2.7,
    "Sb-WT": 2.9,
    "I96Y": 2.3,
    "I96A": 2.6,
    "T137A": 3.0,
    "I206W": 2.8,
    "I206A": 3.0,
    "S209Q": 3.0,
    "S209A": 2.7,
    "F239L": 2.9,
    "F239A": 3.0,
    "L263Q": 2.9,
    "L263A": 3.0,
    "E268D": 3.0,
    "E268A": 3.0,
}

# Enzyme, variant, compound, kcat, kcat error, Km, Km error, enzyme nM,
# maximum substrate mM, minimum substrate uM. Values are transcribed from SI
# Tables S2-S4; only finite Michaelis-Menten fits are represented.
ROWS = (
    ("No-PDE", "WT", 1, 6.7, 0.2, 0.62, 0.04, 47, 2.5, 47),
    ("No-PDE", "WT", 2, 49, 3, 2.3, 0.2, 154, 5.5, 130),
    ("No-PDE", "WT", 11, 0.70, 0.04, 1.4, 0.2, 47, 5.0, 94),
    ("Sb-PDE", "WT", 1, 19.8, 0.6, 0.14, 0.02, 15, 2.5, 48),
    ("Sb-PDE", "WT", 2, 49, 4, 0.9, 0.1, 149, 1.4, 26),
    ("Sb-PDE", "WT", 11, 0.63, 0.02, 0.14, 0.01, 149, 2.6, 48),
    ("Sb-PDE", "WT", 12, 0.213, 0.008, 1.36, 0.09, 1500, 2.0, 38),
    ("No-PDE", "I96Y", 11, 0.21, 0.01, 2.11, 0.08, 84, 5.32, 67),
    ("No-PDE", "I96A", 1, 19, 2, 0.49, 0.06, 23, 0.2, 1),
    ("No-PDE", "I96A", 11, 1.0, 0.1, 5.0, 0.5, 231, 5.32, 17),
    ("No-PDE", "I96A", 2, 56, 8, 5.0, 1.0, 96, 3.4, 64),
    ("No-PDE", "I96A", 12, 0.20, 0.04, 13, 3, 1900, 5.0, 94),
    ("No-PDE", "T137A", 1, 1.89, 0.06, 0.77, 0.06, 27, 2.4, 8),
    ("No-PDE", "T137A", 11, 0.09, 0.01, 3.7, 0.4, 2700, 5.2, 16),
    ("No-PDE", "T137A", 2, 12, 2, 3.3, 0.7, 225, 3.4, 64),
    ("No-PDE", "I206W", 1, 0.60, 0.01, 0.05, 0.004, 10.4, 2.4, 8),
    ("No-PDE", "I206W", 11, 0.30, 0.01, 0.22, 0.02, 104, 5.2, 16),
    ("No-PDE", "I206W", 2, 5.0, 0.2, 0.88, 0.09, 433, 2.72, 64),
    ("No-PDE", "I206W", 12, 0.21, 0.01, 0.37, 0.04, 4300, 4.0, 188),
    ("No-PDE", "I206A", 11, 1.83, 0.06, 1.4, 0.1, 23, 5.2, 16),
    ("No-PDE", "I206A", 2, 85, 7, 6.2, 0.7, 94, 3.4, 64),
    ("No-PDE", "S209Q", 1, 6.2, 0.3, 0.30, 0.03, 7.3, 2.4, 8),
    ("No-PDE", "S209Q", 11, 0.60, 0.02, 0.60, 0.07, 73, 5.2, 16),
    ("No-PDE", "S209Q", 2, 77, 3, 2.1, 0.2, 30, 3.4, 64),
    ("No-PDE", "S209A", 1, 10.9, 0.5, 0.90, 0.08, 10.4, 2.4, 8),
    ("No-PDE", "S209A", 11, 1.09, 0.09, 2.32, 0.09, 52, 5.2, 16),
    ("No-PDE", "S209A", 2, 65, 7, 1.7, 0.4, 22, 3.4, 64),
    ("No-PDE", "F239L", 11, 1.34, 0.03, 2.7, 0.1, 69, 5.2, 16),
    ("No-PDE", "F239L", 2, 68, 8, 4.0, 0.8, 29, 3.4, 64),
    ("No-PDE", "F239A", 1, 6.8, 0.4, 4.3, 0.4, 39, 2.4, 8),
    ("No-PDE", "F239A", 11, 1.50, 0.04, 2.2, 0.1, 79, 5.2, 16),
    ("No-PDE", "F239A", 12, 1.8, 0.1, 17, 1, 1600, 5.0, 188),
    ("No-PDE", "L263Q", 1, 8.7, 0.35, 4.9, 0.4, 53, 2.4, 8),
    ("No-PDE", "L263Q", 11, 1.7, 0.1, 13, 1, 533, 5.2, 16),
    ("No-PDE", "L263A", 11, 5.8, 0.3, 6.3, 0.4, 32, 5.2, 16),
    ("No-PDE", "E268D", 11, 0.83, 0.03, 4.8, 0.2, 405, 5.2, 16),
    ("No-PDE", "E268A", 1, 16.9, 0.8, 2.5, 0.2, 65, 2.54, 8),
    ("No-PDE", "E268A", 11, 1.87, 0.07, 1.8, 0.01, 65, 5.32, 17),
    ("No-PDE", "E268A", 2, 70, 6, 8.1, 0.9, 269, 3.4, 64),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(accession: str) -> str:
    path = RAW / f"{accession}.fasta"
    lines = path.read_text(encoding="ascii").splitlines()
    expected_length = ACCESSIONS[accession][0]
    if not lines or accession not in lines[0]:
        raise ValueError(f"Expected {accession} in {path.name} header")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if len(sequence) != expected_length or not sequence.isalpha() or not sequence.isupper():
        raise ValueError(f"Invalid {accession} sequence")
    return sequence


def substitute(sequence: str, mutations: tuple[tuple[int, str, str], ...]) -> str:
    result = sequence
    for position, expected, replacement in mutations:
        if result[position - 1] != expected:
            raise ValueError(
                f"Expected {expected}{position}, found {result[position - 1]}{position}"
            )
        result = result[: position - 1] + replacement + result[position:]
    return result


def sequence_id(enzyme: str, variant: str) -> str:
    prefix = "no-pde" if enzyme == "No-PDE" else "sb-pde"
    return f"{prefix}-{variant.lower()}"


def build_sequences() -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    native = {accession: read_fasta(accession) for accession in ACCESSIONS}
    no_pde = native["QCI95559.1"]
    constructs = {
        "no-pde-wt": ("QCI95559.1", no_pde),
        "sb-pde-wt": ("BAX25588.1", native["BAX25588.1"]),
    }
    for variant, mutations in MUTATIONS.items():
        if variant != "WT":
            constructs[sequence_id("No-PDE", variant)] = (
                "QCI95559.1",
                substitute(no_pde, mutations),
            )
    return constructs, native


def verify_raw_evidence() -> dict[int, tuple[int, str]]:
    article = ET.parse(RAW / "PMC11656708-fullText.xml").getroot()
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected article DOI")
    license_text = " ".join(article.find(".//permissions").itertext())
    if "creativecommons.org/licenses/by/4.0" not in license_text:
        raise ValueError("Article is not verified CC-BY-4.0")
    table_text = " ".join(article.find(".//table-wrap[@id='tbl1']").itertext())
    for value in ("6.7", "19.8", "0.63", "0.21"):
        if value not in table_text:
            raise ValueError(f"Main Table 1 value {value} not found")

    supplement = RAW / "bi4c00350_si_001.pdf"
    if supplement.read_bytes()[:5] != b"%PDF-":
        raise ValueError("Supporting information is not a PDF")
    with fitz.open(supplement) as document:
        text = "\n".join(page.get_text() for page in document)
    for phrase in (
        "Table S2. Experimental conditions and rate constants",
        "6.7 ± 0.2",
        "0.213 ± 0.008",
        "Table S3: Kinetic constants for variants",
        "16.9 ± 0.8",
        "All reactions were conducted at 22 °C",
    ):
        if phrase not in text:
            raise ValueError(f"Supporting evidence not found: {phrase}")
    with zipfile.ZipFile(RAW / "PMC11656708-supplementaryFiles.zip") as archive:
        if "bi4c00350_si_001.pdf" not in archive.namelist():
            raise ValueError("Supplementary archive lacks the source PDF")

    payload = json.loads((RAW / "pubchem-substrates.json").read_text(encoding="utf-8"))
    observed = {
        item["CID"]: item.get("SMILES", item.get("IsomericSMILES"))
        for item in payload["PropertyTable"]["Properties"]
    }
    expected = {values[1]: values[2] for values in SUBSTRATES.values()}
    if observed != expected:
        raise ValueError(f"PubChem structure mismatch: {observed}")
    return {number: (values[1], values[2]) for number, values in SUBSTRATES.items()}


def build_records() -> list[dict[str, object]]:
    records = []
    for index, row in enumerate(ROWS, 1):
        enzyme, variant, compound, kcat, kcat_error, km, km_error, enzyme_nm, high, low = row
        name, cid, smiles, ester_class, points = SUBSTRATES[compound]
        accession = "QCI95559.1" if enzyme == "No-PDE" else "BAX25588.1"
        mutations = MUTATIONS[variant]
        records.append(
            {
                "candidate_id": f"php-pde-{index:03d}",
                "article_doi": DOI,
                "source_table": "Supporting Table S2" if variant == "WT" else "Supporting Tables S3-S4",
                "source_row": f"{enzyme} {variant}; compound {compound}",
                "organism": ACCESSIONS[accession][1],
                "sequence_accession": accession,
                "sequence_id": sequence_id(enzyme, variant),
                "enzyme_variant": variant,
                "amino_acid_changes": "/".join(
                    f"{old}{position}{new}" for position, old, new in mutations
                ) or "WT",
                "construct": "pET-28a synthetic full precursor; native Sec/SPI signal peptide residues 1-28 efficiently processed in E. coli; assayed catalytic chain residues 29-441",
                "expressed_precursor_length_aa": 441,
                "assayed_catalytic_chain_coordinates": "29-441",
                "assayed_catalytic_chain_length_aa": 413,
                "variable_substrate": name,
                "compound_number": compound,
                "phosphate_ester_class": ester_class,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "reaction_product": "4-nitrophenol" if compound in {1, 11} else "phenol",
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "km_mM": km,
                "km_error_mM": km_error,
                "assay_method": "direct initial-rate UV/visible spectrophotometric Michaelis-Menten fit",
                "assay_buffer": "50 mM CHES",
                "assay_pH": 9.0,
                "assay_temperature_C": 22,
                "enzyme_concentration_nM": enzyme_nm,
                "substrate_min_uM": low,
                "substrate_max_mM": high,
                "substrate_concentrations": points,
                "monitoring_wavelength_nm": 400 if compound in {1, 11} else 275,
                "metal": "trinuclear Zn site as purified; no metal added to kinetic assay",
                "measured_zinc_per_enzyme": ZINC_PER_ENZYME["Sb-WT" if enzyme == "Sb-PDE" else variant],
                "error_type": "author-reported nonlinear-regression fitting error",
                "replicates": "not reported for spectrophotometric titrations",
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def read_hits(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {
        line.split("\t", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    }


def apply_selection(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.exists() or not cluster_path.exists():
        return
    hits = read_hits(hit_path)
    clusters = {}
    for line in cluster_path.read_text(encoding="utf-8").splitlines():
        if line:
            representative, member = line.split("\t")[:2]
            clusters[member] = representative

    family_counts: dict[str, int] = {}
    for record in records:
        sequence = str(record["sequence_id"])
        if sequence in hits:
            record["status_at_normalization"] = "excluded_homology"
            continue
        family = clusters.get(sequence, sequence)
        if family_counts.get(family, 0) >= FAMILY_CAP:
            record["status_at_normalization"] = "excluded_family_cap"
        else:
            record["status_at_normalization"] = "accepted_homology_cold_pool"
            family_counts[family] = family_counts.get(family, 0) + 1


def write_fasta(path: Path, entries: dict[str, tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, (description, sequence) in entries.items():
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_outputs() -> None:
    verify_raw_evidence()
    constructs, native = build_sequences()
    records = build_records()
    apply_selection(records)

    construct_path = SOURCE / "construct_sequences.fasta"
    write_fasta(
        construct_path,
        {
            identifier: (
                f"{accession} | exact synthetic full precursor; residues 1-28 processed; variant {identifier.rsplit('-', 1)[-1].upper()}",
                sequence,
            )
            for identifier, (accession, sequence) in constructs.items()
        },
    )
    reference_path = SOURCE / "family_reference_sequences.fasta"
    write_fasta(
        reference_path,
        {
            accession: (
                f"{name} | exact NCBI accession sequence | {organism}",
                native[accession],
            )
            for accession, (_, organism, name) in ACCESSIONS.items()
        },
    )

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    raw_files = sorted(
        path for path in RAW.iterdir() if path.is_file() and path.name != "PMC11656708.tar.gz"
    )
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC11656708",
        "article_doi": DOI,
        "article_published": "2024-12-02",
        "article_license": "CC-BY-4.0",
        "reported_finite_direct_steady_state_kcat_records": len(records),
        "accepted_records": len(accepted),
        "unique_assayed_construct_sequences": len(constructs),
        "family_cap": FAMILY_CAP,
        "kinetics_sources": [
            "raw/PMC11656708-fullText.xml, Table 1",
            "raw/bi4c00350_si_001.pdf, Tables S1-S4 and Spectrophotometric Kinetic Assays",
        ],
        "sequence_sources": [f"raw/{accession}.fasta" for accession in ACCESSIONS],
        "structure_source": "raw/pubchem-substrates.json; PubChem PUG REST structures for CIDs 255, 13282, 378, and 12793, retrieved 2026-07-22",
        "construct_audit": "No-PDE and Sb-PDE were synthetic full-length precursor genes in pET-28a. The native residues 1-28 Sec/SPI signal peptide was efficiently processed in E. coli, giving assayed residues 29-441. FASTA stores exact translated precursors so article-numbered mutations and accession identity remain exact.",
        "reference_sequence_audit": "BAX25588.1 (Sb-PDE), QCI95559.1 (No-PDE), CCW15874.1 (Si-PDE), and ABR41129.1 (3E38-PHP) were downloaded exactly from NCBI. Si-PDE and 3E38-PHP are family context only and have no candidate kinetic rows.",
        "assay": {
            "method": "Direct initial-rate UV/visible spectrophotometry fitted to Michaelis-Menten kinetics",
            "buffer": "50 mM CHES, pH 9.0",
            "temperature_C": 22,
            "substrate_points": "32 for compounds 1/11; 18 for compounds 2/12",
            "metal": "No metal added; as-purified proteins contained 2.3-3.0 Zn/enzyme by atomic absorption",
            "error_type": "author-reported nonlinear-regression fitting error",
            "replicates": "not reported for spectrophotometric titrations",
        },
        "selection_policy": "All and only finite directly fitted steady-state Michaelis-Menten kcat rows; source-table order retained before homology and family-cap selection; no kinetic-value-based selection.",
        "exclusions": {
            "efficiency_only_or_non_saturating_rows": 33,
            "censored_or_no_reaction_rows": 12,
            "efficiency_only_detail": "Compounds 3-10 were pseudo-first-order 31P-NMR kcat/Km measurements. No-PDE compound 12 and mutant nd rows used low-substrate linear slopes because saturation was insufficient. None supplies finite direct kcat.",
            "relative_selectivity": "Figure 4 diesterase/phosphatase ratios are derived relative endpoints and are not independent absolute kinetics.",
            "inactive_Si_PDE": "Si-PDE did not express active protein; no kinetic parameter was reported.",
            "metal_tests": "Metal supplementation, chelator resistance, and Zn stoichiometry are mechanistic/constitution measurements, not kcat rows.",
        },
        "normalization": "Author-reported kcat in s-1 and Km in mM; no refitting, graph digitization, or unit conversion.",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        audit = {
            "audited_on": "2026-07-22",
            "source_id": SOURCE_ID,
            "article_doi": DOI,
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": len(records),
            "unique_assayed_construct_sequences": len(constructs),
            "construct_sequences_sha256": sha256(construct_path),
            "family_reference_sequences_sha256": sha256(reference_path),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
            "exact_sequence_overlap": sum(float(line.split("\t")[2]) == 1.0 for line in hits),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster_path),
            "family_cap": FAMILY_CAP,
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}),
            "accepted_unique_substrates": len({record["variable_substrate"] for record in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )

    verify_outputs(records, constructs, native)


def verify_outputs(
    records: list[dict[str, object]],
    constructs: dict[str, tuple[str, str]],
    native: dict[str, str],
) -> None:
    if len(records) != 39 or len(constructs) != 15 or len(native) != 4:
        raise ValueError("Unexpected curation cardinality")
    if set(native) != set(ACCESSIONS):
        raise ValueError("A requested native accession is missing")
    for record in records:
        for field in ("kcat_s-1", "kcat_error_s-1", "km_mM", "km_error_mM"):
            value = float(record[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid {field} in {record['candidate_id']}")
    statuses = [str(record["status_at_normalization"]) for record in records]
    if "pending_homology" in statuses:
        raise ValueError("Frozen homology artifacts have not been generated")
    if read_hits(SOURCE / "homology" / "homology_hits.tsv"):
        expected_accepted = 0
    else:
        expected_accepted = FAMILY_CAP
    if statuses.count("accepted_homology_cold_pool") != expected_accepted:
        raise ValueError(f"Unexpected accepted count: {statuses}")
    required_nonempty = (
        SOURCE / "candidate_records.csv",
        SOURCE / "construct_sequences.fasta",
        SOURCE / "family_reference_sequences.fasta",
        SOURCE / "provenance.json",
        SOURCE / "homology-audit.json",
        SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv",
    )
    if not all(path.is_file() and path.stat().st_size for path in required_nonempty):
        raise ValueError("A required standard output is missing or empty")
    if not (SOURCE / "homology" / "homology_hits.tsv").is_file():
        raise ValueError("The homology-hit artifact is missing")


if __name__ == "__main__":
    write_outputs()
