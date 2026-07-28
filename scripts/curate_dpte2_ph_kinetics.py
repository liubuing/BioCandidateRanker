from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "biostudies-S-EPMC12961738"
RAW = SOURCE / "raw"
REFERENCE = ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723" / "homology" / "unikp_reference.fasta"
DOI = "10.1021/acs.biochem.5c00768"
SOURCE_ID = "biostudies-S-EPMC12961738"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
FAMILY_CAP = 20

# Exact post-SuperTEV product from SI protein sequences 1/2. The first ten residues
# are the retained SNIGSGITNS cleavage scar; native P0A434 residue 34 begins at G11.
CLEAVED_D_PTE2 = (
    "SNIGSGITNSGDRINTVRGPITISEAGFTLMHEHICGSSAGFLRAWPEFFGSRDALAEKAVRGLRRARAAGVRTIVDVSTFDIGRDVELLAEVSEAADVHIVAATGLWFDPPLSMRLRSVEELTQFFLREIQYGIEDTGIRAGIIKVATTGKATPFQERVLRAAARASLATGVPVTTHTDASQRDGEQQADIFESEGLDPSRVCIGHSDDTDDLDYLTALAARGYLIGLDHIPHSAIGLEDNASAAALLGLRSWQTRALLIKALIDQGYADQILVSNDWLFGFSSYVTNIMDVMDRVNPDGMAFIPLRVIPFLREKGVPDETLETIMVDNPARFLSPTLRAS"
)
FUSION_PREFIX = "MGSSWSHPQFEKGSSMASMSDSEVNQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGGIEENLYFQ"

SUBSTRATES = {
    "II": {
        "name": "3-fluoro-4-nitrophenyl diethyl phosphate",
        "cid": 85613363,
        "smiles": "CCOP(=O)(OCC)OC1=CC(=C(C=C1)[N+](=O)[O-])F",
        "figure": "Supporting Figure 6",
        "table": "Supporting Table 2",
    },
    "III": {
        "name": "paraoxon",
        "cid": 9395,
        "smiles": "CCOP(=O)(OCC)OC1=CC=C(C=C1)[N+](=O)[O-]",
        "figure": "Supporting Figure 7",
        "table": "Supporting Table 3",
    },
    "IV": {
        "name": "4-acetylphenyl diethyl phosphate",
        "cid": 10880230,
        "smiles": "CCOP(=O)(OCC)OC1=CC=C(C=C1)C(=O)C",
        "figure": "Supporting Figure 8",
        "table": "Supporting Table 4",
    },
    "V": {
        "name": "methyl 4-((diethoxyphosphoryl)oxy)benzoate",
        "cid": 12741195,
        "smiles": "CCOP(=O)(OCC)OC1=CC=C(C=C1)C(=O)OC",
        "figure": "Supporting Figure 9",
        "table": "Supporting Table 5",
    },
}
SUBSTRATE_I = {
    "name": "2,6-difluoro-4-nitrophenyl diethyl phosphate",
    "cid": "",
    "smiles": "CCOP(=O)(OCC)OC1=C(F)C=C([N+](=O)[O-])C=C1F",
    "figure": "Supporting Figure 5",
    "table": "Supporting Table 1",
}

VARIANTS = ("H55(Zn)", "H55piMH(Zn)", "H55(Co)", "H55piMH(Co)")

# SI Tables 2-5, columns in VARIANTS order. Values are author-reported
# log10(kcat / s-1) +/- SD over three biological replicates.
TABLES = {
    "II": """
8.5 3.34 0.05 3.41 0.03 3.82 0.03 3.86 0.02
8.0 3.26 0.04 3.35 0.05 3.72 0.02 3.80 0.03
7.5 3.20 0.06 3.34 0.02 3.61 0.03 3.78 0.02
7.0 3.05 0.05 3.20 0.04 3.48 0.01 3.67 0.03
6.5 2.82 0.06 3.05 0.04 3.34 0.07 3.57 0.03
6.0 2.59 0.05 2.84 0.03 3.10 0.09 3.39 0.07
5.5 2.37 0.08 2.56 0.03 2.90 0.08 3.15 0.10
5.0 1.95 0.04 2.31 0.03 2.72 0.12 3.00 0.01
4.5 - - - - 2.53 0.04 2.89 0.06
""",
    "III": """
8.5 3.36 0.04 3.38 0.01 3.89 0.02 3.87 0.06
8.0 3.30 0.05 3.34 0.01 3.84 0.01 3.84 0.06
7.5 3.23 0.05 3.30 0.01 3.73 0.01 3.79 0.06
7.0 3.10 0.04 3.22 0.02 3.59 0.04 3.69 0.05
6.5 2.90 0.04 3.09 0.02 3.47 0.01 3.60 0.06
6.0 2.62 0.05 2.88 0.01 3.29 0.01 3.43 0.05
5.5 2.34 0.02 2.57 0.03 3.09 0.02 3.25 0.04
5.0 2.00 0.05 2.24 0.02 2.88 0.04 3.02 0.03
4.5 - - - - 2.66 0.02 2.78 0.03
""",
    "IV": """
8.5 2.88 0.06 2.67 0.04 3.45 0.02 3.15 0.01
8.0 2.88 0.06 2.69 0.03 3.43 0.01 3.16 0.01
7.5 2.76 0.05 2.60 0.03 3.31 0.01 3.07 0.01
7.0 2.65 0.04 2.53 0.02 3.17 0.00 2.98 0.03
6.5 2.43 0.03 2.41 0.08 2.95 0.01 2.82 0.04
6.0 2.10 0.09 2.20 0.05 2.65 0.05 2.56 0.01
5.5 1.63 0.02 1.89 0.05 2.27 0.05 2.27 0.07
""",
    "V": """
8.5 2.54 0.05 2.23 0.06 2.86 0.02 2.53 0.03
8.0 2.54 0.06 2.40 0.05 2.79 0.04 2.50 0.02
7.5 2.62 0.05 2.48 0.05 2.79 0.07 2.54 0.00
""",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_table(text: str) -> list[tuple[float, str, float, float]]:
    rows = []
    for line in text.strip().splitlines():
        cells = line.split()
        pH = float(cells[0])
        for index, variant in enumerate(VARIANTS):
            value, error = cells[1 + 2 * index : 3 + 2 * index]
            if value != "-":
                rows.append((pH, variant, float(value), float(error)))
    return rows


def sequence_id(variant: str) -> str:
    return "dpte2-h55-pimh" if "piMH" in variant else "dpte2-h55"


def write_fasta(path: Path) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, description in (
            ("dpte2-h55", "SI-exact post-SuperTEV chain; native H55; 342 aa"),
            ("dpte2-h55-pimh", "SI-exact post-SuperTEV chain; H55piMH represented canonically as H for MMseqs2; 342 residues"),
        ):
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(CLEAVED_D_PTE2), 80):
                handle.write(CLEAVED_D_PTE2[start : start + 80] + "\n")


def write_fusion_fasta(path: Path) -> None:
    sequence = FUSION_PREFIX + CLEAVED_D_PTE2
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, description in (
            ("strep-sumo-tev-dpte2-h55", "SI protein sequence 1; full expressed fusion; 466 aa"),
            ("strep-sumo-tev-dpte2-h55-pimh", "SI protein sequence 2; full expressed fusion; piMH represented canonically as H; 466 residues"),
        ):
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def verify_raw() -> None:
    required = (
        "PMC12961738-fullText.xml",
        "PMC12961738-supplementaryFiles.zip",
        "bi5c00768_si_001.pdf",
        "P0A434.fasta",
        "9TI1-polymer-entity-1.json",
        "9TI2-polymer-entity-1.json",
        "pubchem-substrate-I-no-hit.json",
        "pubchem-substrate-II.json",
        "pubchem-substrate-III.json",
        "pubchem-substrate-IV.json",
        "pubchem-substrate-V.json",
    )
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {', '.join(missing)}")
    article = ET.parse(RAW / "PMC12961738-fullText.xml").getroot()
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected article DOI")
    permissions = " ".join(article.find(".//permissions").itertext())
    if "creativecommons.org/licenses/by/4.0" not in permissions:
        raise ValueError("Article and SI are not verified CC-BY-4.0")
    with zipfile.ZipFile(RAW / "PMC12961738-supplementaryFiles.zip") as archive:
        if "bi5c00768_si_001.pdf" not in archive.namelist():
            raise ValueError("Europe PMC supplement archive lacks SI PDF")
    with fitz.open(RAW / "bi5c00768_si_001.pdf") as document:
        pages = [page.get_text() for page in document]
    for page, phrase in (
        (5, "Data points represent a single measurement"),
        (17, "Supplementary Table 1. log(kcat, s-1)"),
        (18, "Supplementary Table 4. log(kcat, s-1)"),
        (35, "Protein sequence 1. strepTagII-SUMO-(tev)-dPTE2-H55"),
    ):
        if phrase not in pages[page]:
            raise ValueError(f"Missing SI evidence on PDF page {page + 1}: {phrase}")
    compact = re.sub(r"[^A-Z]", "", pages[35].split("Protein sequence 1.", 1)[1].split("Protein sequence 2.", 1)[0])
    if not compact.endswith(CLEAVED_D_PTE2):
        raise ValueError("Exact cleaved dPTE2 sequence was not recovered from SI")
    if len(CLEAVED_D_PTE2) != 342 or CLEAVED_D_PTE2[31] != "H":
        raise ValueError("Unexpected cleaved construct coordinates")
    for roman, data in SUBSTRATES.items():
        payload = json.loads((RAW / f"pubchem-substrate-{roman}.json").read_text(encoding="utf-8"))
        item = payload["PropertyTable"]["Properties"][0]
        if item["CID"] != data["cid"] or item["SMILES"] != data["smiles"]:
            raise ValueError(f"PubChem mismatch for substrate {roman}")
    no_hit = json.loads((RAW / "pubchem-substrate-I-no-hit.json").read_text(encoding="utf-8"))
    if no_hit != {"IdentifierList": {"CID": [0]}}:
        raise ValueError("Substrate I PubChem exact-identity result changed")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference hash changed")


def build_records() -> list[dict[str, object]]:
    records = []
    for roman, table in TABLES.items():
        substrate = SUBSTRATES[roman]
        for pH, variant, log_kcat, log_error in parse_table(table):
            metal = "Co2+" if "(Co)" in variant else "Zn2+"
            ncaa = "piMH" in variant
            records.append(
                {
                    "candidate_id": f"dpte2-ph-{len(records) + 1:03d}",
                    "article_doi": DOI,
                    "source_file": "raw/bi5c00768_si_001.pdf",
                    "source_table": substrate["table"],
                    "source_figure": substrate["figure"],
                    "source_row": f"pH {pH:.1f}; {variant}",
                    "organism": "Brevundimonas diminuta (reported as Pseudomonas diminuta)",
                    "sequence_accession": "P0A434",
                    "sequence_id": sequence_id(variant),
                    "enzyme_variant": "dPTE2-H55piMH" if ncaa else "dPTE2-H55",
                    "position_55_residue": "N-pi-methyl-L-histidine" if ncaa else "L-histidine",
                    "position_55_construct_coordinate": 32,
                    "background_mutations_vs_P0A434": "T54M/K77D/S111E/R118E/L182R/K185R/A203D/A214D/S222D/S238D/S269A/I274L/M293A/K294D/Q343D/A347E/G348T/T350M/T352D",
                    "construct": "SI-exact 342-residue post-SuperTEV product: SNIGSGITNS scar + P0A434 residues 34-365 carrying dPTE2 background",
                    "construct_length_aa": 342,
                    "variable_substrate": substrate["name"],
                    "substrate_roman": roman,
                    "substrate_pubchem_cid": substrate["cid"],
                    "substrate_isomeric_smiles": substrate["smiles"],
                    "endpoint": "log10_kcat_s-1",
                    "log10_kcat_s-1": log_kcat,
                    "log10_kcat_sd": log_error,
                    "kcat_s-1": 10**log_kcat,
                    "assay_method": "direct continuous phenolate initial-rate substrate titration fitted to Michaelis-Menten equation",
                    "saturation_evidence": f"{substrate['figure']}; individual curves over 0-1500 uM",
                    "assay_buffer": "40 mM each acetate/MES/HEPES/Tris, 100 mM NaCl, 1 mg/mL BSA, 5% acetonitrile final",
                    "assay_pH": pH,
                    "assay_temperature_C": 25,
                    "enzyme_concentration_nM": 1,
                    "substrate_min_uM": 0,
                    "substrate_max_uM": 1500,
                    "monitoring_wavelength_nm": "",
                    "active_site_metal": metal,
                    "biological_replicates": 3,
                    "error_type": "SD of author-reported log10(kcat / s-1)",
                    "status_at_normalization": "pending_homology",
                }
            )
    return records


def exclusions() -> list[dict[str, object]]:
    return [
        {
            "source_rows": "Supporting Table 1, all 34 finite substrate-I H2O kcat cells",
            "count": 34,
            "reason": "defined_structure_has_no_exact_PubChem_CID; PubChem identity query returned CID 0; no guessed mapping",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Supporting Tables 11-12, all 51 finite D2O kcat cells",
            "count": 51,
            "reason": "D2O_solvent_isotope_experiment_outside_requested_H2O_pH_resolved_labels",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Supporting Tables 6-10 and 13-14",
            "count": 3,
            "reason": "kcat_per_Km_only_cells_without_finite_kcat_not_promoted_to_kcat_labels",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Supporting Figure 1 ncAA screen",
            "count": 7,
            "reason": "uncleaved_full_fusion_single_biological_replicate_screen; not the exact cleaved assay construct",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Supporting Tables 15-18 and 32-34",
            "count": 124,
            "reason": "pH_profile_or_Bronsted_model_derived_parameters_not_direct_row_specific_Michaelis_Menten_kcat",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Supporting Tables 19-31",
            "count": 13,
            "reason": "D2O_pKa_KSIE_delta_pKa_and_model_derived_summary_tables_not_direct_H2O_Michaelis_Menten_kcat",
            "candidate_label_created": False,
        },
    ]


def read_hit_ids(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    return {line.split("\t", 1)[0] for line in path.read_text(encoding="utf-8").splitlines() if line}


def apply_homology(records: list[dict[str, object]]) -> None:
    hits = read_hit_ids(SOURCE / "homology" / "homology_hits.tsv")
    if not (SOURCE / "homology" / "homology_hits.tsv").is_file():
        return
    for record in records:
        record["status_at_normalization"] = "excluded_homology" if record["sequence_id"] in hits else "accepted_homology_cold_pool"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def substrate_mapping_rows() -> list[dict[str, object]]:
    rows = []
    for roman, data in (("I", SUBSTRATE_I), *SUBSTRATES.items()):
        rows.append(
            {
                "substrate_roman": roman,
                "si_name": data["name"],
                "si_structure_smiles": data["smiles"],
                "pubchem_cid": data["cid"],
                "pubchem_mapping_status": "exact_structure_no_hit_CID_0" if roman == "I" else "exact",
                "kinetics_table": data["table"],
                "saturation_figure": data["figure"],
                "candidate_rows_created": roman != "I",
            }
        )
    return rows


def write_outputs(*, prepare_only: bool) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    verify_raw()
    records = build_records()
    construct_path = SOURCE / "construct_sequences.fasta"
    write_fasta(construct_path)
    fusion_path = SOURCE / "expressed_fusion_sequences.fasta"
    write_fusion_fasta(fusion_path)
    write_csv(SOURCE / "exclusions.csv", exclusions())
    write_csv(SOURCE / "substrate-mapping-audit.csv", substrate_mapping_rows())
    if prepare_only:
        write_csv(SOURCE / "candidate_records.csv", records)
        return
    apply_homology(records)
    if any(record["status_at_normalization"] == "pending_homology" for record in records):
        raise ValueError("Pinned MMseqs2 output must exist before finalization")
    write_csv(SOURCE / "candidate_records.csv", records)
    write_audits(records, construct_path, fusion_path)
    verify_outputs(records)


def write_audits(records: list[dict[str, object]], construct_path: Path, fusion_path: Path) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    raw_files = sorted(path for path in RAW.iterdir() if path.is_file() and path.stat().st_size)
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12961738",
        "article_doi": DOI,
        "pmc_id": "PMC12961738",
        "article_published": "2026-02-18",
        "article_and_supplement_license": "CC-BY-4.0",
        "finite_direct_H2O_kcat_cells_in_SI_tables_1_to_5": 142,
        "structure_compliant_candidate_records_before_homology": len(records),
        "accepted_records": len(accepted),
        "construct_audit": {
            "fusion": "SI sequence is StrepTagII-SUMO-(tev)-dPTE2; main Methods report overnight SuperTEV cleavage followed by removal of protease and cleaved N terminus. Both full expressed fusion variants are retained in expressed_fusion_sequences.fasta.",
            "cleavage": "The fusion junction is ...EENLYFQ|SNIGSGITNSGDR...; cleavage after Q leaves the exact 342-residue product beginning SNIGSGITNSGDR.",
            "background": "The cleaved product is a 10-residue SNIGSGITNS scar plus mature P0A434 residues 34-365 with the SI-exact dPTE2 substitutions. The full sequence is stored verbatim in construct_sequences.fasta.",
            "H55": "Native precursor H55 maps to cleaved-construct position 32. SI protein sequence 2 marks this residue ncAA; LC-MS/MS and PDB 9TI2 identify it as piMH/MHS.",
            "canonical_homology_representation": "MMseqs2 cannot encode piMH. Both exact backbones are queried with H at construct position 32, matching the RCSB canonical sequence representation; chemical residue identity remains explicit per row.",
        },
        "kinetics_evidence": "raw/bi5c00768_si_001.pdf, Supporting Figures 6-9 and Tables 2-5; direct initial-rate Michaelis-Menten fits at each pH, mean of three biological replicates",
        "selection_policy": "All and only finite direct H2O pH-resolved kcat cells for substrates II-V with exact PubChem identity, before pinned homology exclusion; no value-based selection, graph digitization, refitting, prediction, or imputation.",
        "substrate_I_exclusion": "2,6-difluoro-4-nitrophenyl diethyl phosphate is structurally defined and synthesis-confirmed in the SI, but exact PubChem structure lookup returned CID 0. All 34 rows fail closed rather than borrowing an analogue CID.",
        "normalization": "Primary label retained exactly as author-reported log10(kcat / s-1); linear kcat is deterministic 10**label and is not treated as an independently reported value.",
        "acquisition_note": "NCBI OA API advertised PMC12961738.tar.gz but its HTTPS endpoint returned 404 on 2026-07-27. Europe PMC CC-BY full-text XML and supplementary archive/PDF were retained instead. Main PDF acquisition returned an empty response and is not an evidentiary original.",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "artifact_sha256": {
            "candidate_records.csv": sha256(SOURCE / "candidate_records.csv"),
            "exclusions.csv": sha256(SOURCE / "exclusions.csv"),
            "construct_sequences.fasta": sha256(construct_path),
            "expressed_fusion_sequences.fasta": sha256(fusion_path),
            "substrate-mapping-audit.csv": sha256(SOURCE / "substrate-mapping-audit.csv"),
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")
    audit = {
        "audited_on": "2026-07-27",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search",
        "mmseqs_version": MMSEQS_VERSION,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": len(records),
        "unique_assayed_construct_sequences": 2,
        "unique_canonical_homology_sequences": 1,
        "construct_sequences_sha256": sha256(construct_path),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
        "exact_sequence_overlap": sum(len(parts) > 2 and float(parts[2]) == 1.0 for parts in (line.split("\t") for line in hits)),
        "homology_hits_sha256": sha256(hit_path),
        "candidate_mmseqs_families": 0 if hits else 1,
        "family_cluster_sha256": sha256(SOURCE / "homology" / "family-cluster.tsv"),
        "family_cap": FAMILY_CAP,
        "accepted_records": len(accepted),
        "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
        "accepted_unique_substrates": len({row["variable_substrate"] for row in accepted}),
        "readiness_gate_passes": False,
        "claim_boundary": "Curation/exclusion evidence only; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


def verify_outputs(records: list[dict[str, object]]) -> None:
    if len(records) != 108 or len(parse_table(TABLES["II"])) != 34 or len(parse_table(TABLES["III"])) != 34:
        raise ValueError("Unexpected candidate cardinality")
    for row in records:
        if not math.isclose(math.log10(float(row["kcat_s-1"])), float(row["log10_kcat_s-1"]), abs_tol=1e-12):
            raise ValueError(f"kcat conversion mismatch in {row['candidate_id']}")
    with (SOURCE / "exclusions.csv").open(newline="", encoding="utf-8") as handle:
        excluded = list(csv.DictReader(handle))
    if any(row["candidate_label_created"] != "False" for row in excluded):
        raise ValueError("Excluded evidence leaked into candidate labels")
    required = (
        SOURCE / "candidate_records.csv",
        SOURCE / "exclusions.csv",
        SOURCE / "construct_sequences.fasta",
        SOURCE / "expressed_fusion_sequences.fasta",
        SOURCE / "substrate-mapping-audit.csv",
        SOURCE / "provenance.json",
        SOURCE / "homology-audit.json",
        SOURCE / "homology" / "homology_hits.tsv",
        SOURCE / "homology" / "family-cluster.tsv",
    )
    if not all(path.is_file() and path.stat().st_size for path in required):
        raise ValueError("A required curator artifact is missing or empty")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    write_outputs(prepare_only=args.prepare_only)
