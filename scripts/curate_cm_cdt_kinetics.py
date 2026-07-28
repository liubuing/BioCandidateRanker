from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
import urllib.request
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
    / "europepmc-PMC10520331"
)
RAW = SOURCE / "raw"
DOI = "10.1016/j.jbc.2023.105161"
SOURCE_ID = "europepmc-PMC10520331"
REFERENCE = (
    ROOT
    / "artifacts"
    / "external"
    / "absolute-kinetics-screen"
    / "dryad-4964723"
    / "homology"
    / "unikp_reference.fasta"
)
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
MIN_IDENTITY = 0.30
COVERAGE = 0.80
COVERAGE_MODE = 0
MAX_SUBSTRATE_UM = 100.0
SATURATION_MULTIPLE = 5.0

DOWNLOADS = {
    "PMC10520331-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10520331/fullTextXML",
    "PMC10520331-article.pdf": "https://europepmc.org/articles/PMC10520331?pdf=render",
    "PMC10520331-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10520331/supplementaryFiles?includeInlineImage=false",
    "8CQ3.cif": "https://files.rcsb.org/download/8CQ3.cif",
    "8CQ4.cif": "https://files.rcsb.org/download/8CQ4.cif",
    "8CQ6.cif": "https://files.rcsb.org/download/8CQ6.cif",
    "pubchem-chorismate.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/12039/property/IsomericSMILES,Title/JSON",
    "pubchem-prephenate.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1028/property/IsomericSMILES,Title/JSON",
}

COMPOUNDS = {
    "chorismate": {
        "cid": 12039,
        "smiles": "C=C(C(=O)O)O[C@@H]1C=C(C=C[C@H]1O)C(=O)O",
    },
    "prephenate": {
        "cid": 1028,
        "smiles": "C1=CC(C=CC1O)(CC(=O)C(=O)O)C(=O)O",
    },
}

PARENTS = {
    "Af": {
        "sequence_id": "cmcdt-af-wt",
        "accession": "WP_083814300.1",
        "organism": "Aequoribacter fuscus",
        "topology": "CM-CDT",
        "mature_marker": "DNHSEQ",
        "pdb_id": "8CQ3",
        "cm_nM": 10.0,
        "cdt_nM": 5.0,
    },
    "Sc": {
        "sequence_id": "cmcdt-sc-wt",
        "accession": "WP_116808336.1",
        "organism": "Steroidobacter cummioxidans",
        "topology": "CM-CDT",
        "mature_marker": "ASFTGP",
        "pdb_id": "",
        "cm_nM": 2.5,
        "cdt_nM": 2.5,
    },
    "Sb": {
        "sequence_id": "cmcdt-sb-wt",
        "accession": "ABN63218.1",
        "organism": "Shewanella baltica",
        "topology": "CM-CDT",
        "mature_marker": "DTDENA",
        "pdb_id": "",
        "cm_nM": 5.0,
        "cdt_nM": 2.5,
    },
    "Jb": {
        "sequence_id": "cdtcm-jb-wt",
        "accession": "ELX09769.1",
        "organism": "Janthinobacterium sp. HH01",
        "topology": "CDT-CM",
        "mature_marker": "DHRLDD",
        "pdb_id": "8CQ4",
        "cm_nM": 2.5,
        "cdt_nM": 2.5,
    },
    "Ds": {
        "sequence_id": "cdtcm-ds-wt",
        "accession": "WP_072786685.1",
        "organism": "Duganella sacchari",
        "topology": "CDT-CM",
        "mature_marker": "GRLEEI",
        "pdb_id": "8CQ6",
        "cm_nM": 2.5,
        "cdt_nM": 2.5,
    },
    "Mp": {
        "sequence_id": "cdtcm-mp-wt",
        "accession": "PQP01982.1",
        "organism": "Massilia phosphatilytica",
        "topology": "CDT-CM",
        "mature_marker": "GHLDDI",
        "pdb_id": "",
        "cm_nM": 2.5,
        "cdt_nM": 2.5,
    },
}

# Source cells from article Tables 1 and 3: code, activity, kcat, SD, Km, SD.
DIRECT_ROWS = (
    ("Af", "CM", 2.2, 0.3, 1.8, 0.1),
    ("Af", "CDT", 29.5, 2.4, 39.9, 0.5),
    ("Sc", "CM", 6.8, 0.7, 6.3, 3.9),
    ("Sc", "CDT", 32.0, 8.4, 16.6, 2.5),
    ("Sb", "CM", 6.8, 0.1, 3.2, 2.7),
    ("Sb", "CDT", 12.2, 2.8, 3.6, 1.1),
    ("Jb", "CM", 17.1, 0.9, 3.0, 0.1),
    ("Jb", "CDT", 13.0, 1.3, 7.3, 2.0),
    ("Ds", "CM", 19.5, 0.2, 3.0, 0.8),
    ("Ds", "CDT", 24.9, 2.0, 14.6, 3.1),
    ("Mp", "CM", 26.3, 1.4, 3.2, 0.0),
    ("Mp", "CDT", 41.8, 3.5, 18.1, 1.7),
    ("Af-E353Q", "CM", 1.9, 0.1, 1.7, 0.3),
    ("Af-K48A", "CDT", 33.3, 2.9, 17.9, 0.7),
    ("Af-CM", "CM", 2.9, 0.5, 4.5, 0.4),
    ("Af-CDT", "CDT", 20.2, 4.1, 28.2, 1.6),
    ("Jb-E200Q", "CM", 15.6, 0.4, 2.4, 1.3),
    ("Jb-K287A", "CDT", 13.9, 0.3, 4.2, 0.4),
    ("Jb-CM", "CM", 16.4, 1.5, 2.2, 0.3),
)

VARIANT_META = {
    "Af-E353Q": ("cmcdt-af-e353q", "Af", "full fusion; CDT active-site KO E353Q", 10.0),
    "Af-K48A": ("cmcdt-af-k48a", "Af", "full fusion; CM active-site KO K48A", 2.5),
    "Af-CM": ("cm-af-split", "Af", "split CM domain; N-terminal MHHHHHHSSG tag", 10.0),
    "Af-CDT": ("cdt-af-split", "Af", "split CDT domain; C-terminal LEHHHHHH tag", 2.5),
    "Jb-E200Q": ("cdtcm-jb-e200q", "Jb", "full fusion; CDT active-site KO E200Q", 2.5),
    "Jb-K287A": ("cdtcm-jb-k287a", "Jb", "full fusion; CM active-site KO K287A", 2.5),
    "Jb-CM": ("cm-jb-split", "Jb", "split CM domain; C-terminal LEHHHHHH tag", 2.5),
}

CODON_TABLE = {
    codon: amino_acid
    for amino_acid, codons in {
        "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG",
        "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
        "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG",
        "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
        "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC",
        "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
        "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC",
        "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
        "G": "GGT GGC GGA GGG",
    }.items()
    for codon in codons.split()
}


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in DOWNLOADS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size > 100:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    payload = response.read()
                if len(payload) <= 100:
                    raise ValueError(f"Empty or truncated response for {url}")
                path.write_bytes(payload)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    with zipfile.ZipFile(RAW / "PMC10520331-supplementaryFiles.zip") as archive:
        payload = archive.read("mmc1.pdf")
    (RAW / "mmc1.pdf").write_bytes(payload)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def translate(dna: str) -> str:
    start = dna.find("ATG")
    if start < 0:
        raise ValueError("Expression DNA lacks a start codon")
    protein = []
    for index in range(start, len(dna) - 2, 3):
        residue = CODON_TABLE[dna[index : index + 3]]
        if residue == "*":
            return "".join(protein)
        protein.append(residue)
    raise ValueError("Expression DNA lacks an in-frame stop codon")


def si_parent_sequences() -> dict[str, str]:
    expected_lengths = {
        "WP_083814300.1": 424,
        "WP_116808336.1": 435,
        "ABN63218.1": 429,
        "ELX09769.1": 425,
        "WP_072786685.1": 416,
        "PQP01982.1": 429,
    }
    with fitz.open(RAW / "mmc1.pdf") as document:
        text = "\n".join(document[index].get_text() for index in range(10, 15))
    entries = re.split(r"^>\s*", text, flags=re.MULTILINE)[1:]
    by_accession = {}
    for entry in entries:
        lines = entry.splitlines()
        header = lines[0]
        match = re.search(r"accession no\. ([A-Z0-9_.]+)", header)
        if not match or match.group(1) not in expected_lengths:
            continue
        dna = "".join(
            line.strip().upper()
            for line in lines[1:]
            if re.fullmatch(r"[ACGTacgt]{20,}", line.strip())
        )
        by_accession[match.group(1)] = translate(dna)
    observed = {accession: len(by_accession.get(accession, "")) for accession in expected_lengths}
    if observed != expected_lengths:
        raise ValueError(f"Unexpected SI expression sequences: {observed}")
    return by_accession


def cif_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s*\n;([^;]+);", text, re.DOTALL
    )
    if not match:
        raise ValueError(f"Canonical polymer sequence missing from {path.name}")
    return re.sub(r"\s+", "", match.group(1))


def mutate_precursor(sequence: str, position: int, expected: str, replacement: str) -> str:
    residues = list(sequence)
    if residues[position - 1] != expected:
        raise ValueError(f"Expected precursor residue {expected}{position}, found {residues[position - 1]}")
    residues[position - 1] = replacement
    return "".join(residues)


def leaderless(sequence: str, marker: str) -> str:
    if not sequence.endswith("LEHHHHHH") or sequence.count(marker) != 1:
        raise ValueError(f"Cannot locate exact mature tagged construct using marker {marker}")
    return "M" + sequence[sequence.index(marker) :]


def build_sequences() -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    parent_precursors = si_parent_sequences()
    sequences: dict[str, dict[str, str]] = {}
    full_by_code = {}
    for code, meta in PARENTS.items():
        full = parent_precursors[str(meta["accession"])]
        full_by_code[code] = full
        if meta["pdb_id"] and cif_sequence(RAW / f"{meta['pdb_id']}.cif") != full:
            raise ValueError(f"PDB {meta['pdb_id']} does not match SI sequence for {code}")
        sequences[str(meta["sequence_id"])] = {
            "sequence": leaderless(full, str(meta["mature_marker"])),
            "construct": "leaderless cytoplasmic full fusion; C-terminal LEHHHHHH tag",
            "source": f"SI Figure S2 and cloning primer for {code}; {meta['accession']}",
        }

    mutations = {
        "Af-E353Q": (353, "E", "Q"),
        "Af-K48A": (48, "K", "A"),
        "Jb-E200Q": (200, "E", "Q"),
        "Jb-K287A": (287, "K", "A"),
    }
    for key, change in mutations.items():
        sequence_id, parent, construct, _ = VARIANT_META[key]
        mutated = mutate_precursor(full_by_code[parent], *change)
        sequences[sequence_id] = {
            "sequence": leaderless(mutated, str(PARENTS[parent]["mature_marker"])),
            "construct": construct,
            "source": f"SI active-site KO primers; mutation numbered on full precursor {PARENTS[parent]['accession']}",
        }

    af = full_by_code["Af"]
    if af[178:180] != "LD":
        raise ValueError("Af split boundary L179/D180 changed")
    sequences[VARIANT_META["Af-CM"][0]] = {
        "sequence": "MHHHHHHSSG" + af[af.index("DNHSEQ") : 179],
        "construct": VARIANT_META["Af-CM"][2],
        "source": "SI AfCM split primers; full-precursor residues through L179",
    }
    sequences[VARIANT_META["Af-CDT"][0]] = {
        "sequence": "M" + af[179:],
        "construct": VARIANT_META["Af-CDT"][2],
        "source": "SI AfCDT split primers; full-precursor residues D180-end",
    }
    jb = full_by_code["Jb"]
    boundary = jb.index("WLDFPWG") + len("WLDFPW")
    sequences[VARIANT_META["Jb-CM"][0]] = {
        "sequence": "M" + jb[boundary:],
        "construct": VARIANT_META["Jb-CM"][2],
        "source": "SI JbCM split primer; Gly immediately after WLDFPW linker through C-terminal tag",
    }
    return sequences, full_by_code


def direct_record(
    candidate_id: str,
    code: str,
    activity: str,
    kcat: float,
    kcat_sd: float,
    km: float,
    km_sd: float,
    sequences: dict[str, dict[str, str]],
) -> dict[str, object]:
    if code in PARENTS:
        parent = code
        sequence_id = str(PARENTS[code]["sequence_id"])
        construct = str(sequences[sequence_id]["construct"])
        enzyme_nM = float(PARENTS[code]["cm_nM" if activity == "CM" else "cdt_nM"])
        source_table = "Table 1"
        source_row = f"*{code}{PARENTS[code]['topology'].replace('-', '')}; {activity} row"
    else:
        sequence_id, parent, construct, enzyme_nM = VARIANT_META[code]
        source_table = "Table 3"
        source_row = f"*{code}; {activity} assay cell"
    parent_meta = PARENTS[parent]
    substrate = "chorismate" if activity == "CM" else "prephenate"
    compound = COMPOUNDS[substrate]
    saturation = MAX_SUBSTRATE_UM / km
    return {
        "candidate_id": candidate_id,
        "article_doi": DOI,
        "source_table": source_table,
        "source_cell": source_row,
        "organism": parent_meta["organism"],
        "sequence_id": sequence_id,
        "accession": parent_meta["accession"],
        "construct": construct,
        "fusion_topology": parent_meta["topology"],
        "catalytic_activity": activity,
        "ec_number": "5.4.99.5" if activity == "CM" else "4.2.1.51",
        "reaction": "chorismate -> prephenate" if activity == "CM" else "prephenate -> phenylpyruvate + H2O",
        "variable_substrate": substrate,
        "substrate_pubchem_cid": compound["cid"],
        "substrate_isomeric_smiles": compound["smiles"],
        "endpoint": "kcat_s-1",
        "kcat_s-1": kcat,
        "kcat_sd_s-1": kcat_sd,
        "km_uM": km,
        "km_sd_uM": km_sd,
        "error_type": "SD of independently fitted parameters",
        "biological_replicates": 2,
        "fit": "individual Michaelis-Menten fit per biological replicate; mean and SD of fitted parameters",
        "assay_type": "direct discontinuous single-activity initial-rate assay",
        "assay_buffer": "50 mM potassium phosphate, 0.1 mg/mL BSA",
        "assay_pH": 7.5,
        "assay_temperature_C": 30,
        "reaction_volume_uL": 200,
        "enzyme_concentration_nM": enzyme_nM,
        "substrate_range_uM": "2.5,5,10,25,50,100",
        "maximum_substrate_uM": MAX_SUBSTRATE_UM,
        "maximum_substrate_multiple_of_km": round(saturation, 3),
        "saturation_threshold_multiple_of_km": SATURATION_MULTIPLE,
        "saturation_passes": saturation >= SATURATION_MULTIPLE,
        "lowest_substrate_brackets_km_2_5_fold": 2.5 <= km / 2.5,
        "source_italicized_low_km_reliability_warning": not (2.5 <= km / 2.5),
        "initial_rate_constraints": "four time points over 0-4 min; substrate turnover kept below 25%; spontaneous background corrected",
        "readout": "phenylpyruvate enolate A320; extinction coefficient 17500 M-1 cm-1",
        "status_at_normalization": "accepted_saturation_pending_homology",
    }


def build_records(sequences: dict[str, dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    direct = [
        direct_record(f"cm-cdt-{index:03d}", *row, sequences)
        for index, row in enumerate(DIRECT_ROWS, 1)
    ]
    accepted = [row for row in direct if row["saturation_passes"]]
    exclusions = []
    for row in direct:
        if row["saturation_passes"]:
            continue
        exclusions.append(
            {
                "article_doi": DOI,
                "source_table": row["source_table"],
                "source_cell": row["source_cell"],
                "sequence_id": row["sequence_id"],
                "catalytic_activity": row["catalytic_activity"],
                "source_kcat_s-1": row["kcat_s-1"],
                "source_kcat_sd_s-1": row["kcat_sd_s-1"],
                "source_km_uM": row["km_uM"],
                "maximum_substrate_multiple_of_km": row["maximum_substrate_multiple_of_km"],
                "exclusion_reason": "substrate_maximum_below_5x_km",
                "candidate_label_created": False,
            }
        )

    coupled_table_1 = {
        "Af": (2.9, 0.3), "Sc": (8.8, 3.0), "Sb": (7.2, 0.8),
        "Jb": (7.7, 0.4), "Ds": (11.7, 0.2), "Mp": (19.0, 3.1),
    }
    for code, (value, error) in coupled_table_1.items():
        exclusions.append(exclusion("Table 1", f"*{code}; CM+CDT row", value, error, "sequential_coupled_apparent_kcat_not_single_activity"))

    nonfinite = (
        ("*AfCMCDT E353Q; CDT assay", "no_activity_detected"),
        ("*AfCMCDT K48A; CM assay", "no_activity_detected"),
        ("*AfCM; CDT assay", "not_measured"),
        ("*AfCDT; CM assay", "not_measured"),
        ("*JbCDTCM E200Q; CDT assay", "no_activity_detected"),
        ("*JbCDTCM K287A; CM assay", "no_activity_detected"),
        ("*JbCM; CDT assay", "not_measured"),
    )
    for cell, reason in nonfinite:
        exclusions.append(exclusion("Table 3", cell, "", "", reason))

    table_4 = (
        ("Af control", 2.9, 0.3, "duplicate_of_table_1_sequential_coupled"),
        ("Jb control", 7.7, 0.4, "duplicate_of_table_1_sequential_coupled"),
        ("Af + Phe", 2.9, 0.1, "sequential_coupled_apparent_kcat_not_single_activity"),
        ("Jb + Phe", 7.8, 1.2, "sequential_coupled_apparent_kcat_not_single_activity"),
        ("Af + Tyr", 3.0, 0.3, "sequential_coupled_apparent_kcat_not_single_activity"),
        ("Jb + Tyr", 8.0, 1.5, "sequential_coupled_apparent_kcat_not_single_activity"),
    )
    for cell, value, error, reason in table_4:
        exclusions.append(exclusion("Table 4", cell, value, error, reason))
    exclusions.append(exclusion("Results text / Figure S12", "mixed Af K48A + E353Q", 2.8, 0.2, "mixed_construct_sequential_coupled_apparent_kcat"))
    exclusions.append(exclusion("Results text / Figure S12", "mixed AfCM + AfCDT", 1.3, 0.2, "mixed_construct_sequential_coupled_apparent_kcat"))
    return accepted, exclusions


def exclusion(table: str, cell: str, value: object, error: object, reason: str) -> dict[str, object]:
    return {
        "article_doi": DOI,
        "source_table": table,
        "source_cell": cell,
        "sequence_id": "",
        "catalytic_activity": "CM+CDT" if "coupled" in reason or "mixed" in reason else "",
        "source_kcat_s-1": value,
        "source_kcat_sd_s-1": error,
        "source_km_uM": "",
        "maximum_substrate_multiple_of_km": "",
        "exclusion_reason": reason,
        "candidate_label_created": False,
    }


def verify_raw() -> None:
    missing = [name for name in DOWNLOADS if not (RAW / name).is_file()]
    if missing or not (RAW / "mmc1.pdf").is_file():
        raise FileNotFoundError(f"Missing raw files: {missing}; run with --download")
    root = ET.parse(RAW / "PMC10520331-fullText.xml").getroot()
    if root.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Article DOI mismatch")
    table_1 = root.find(".//table-wrap[@id='tbl1']")
    table_3 = root.find(".//table-wrap[@id='tbl3']")
    evidence = " ".join("".join(item.itertext()) for item in (table_1, table_3) if item is not None)
    for token in ("29.5", "41.8", "33.3", "20.2", "16.4"):
        if token not in evidence:
            raise ValueError(f"Required kinetic source value missing: {token}")
    for name, expected in COMPOUNDS.items():
        payload = json.loads((RAW / f"pubchem-{name}.json").read_text(encoding="utf-8"))
        item = payload["PropertyTable"]["Properties"][0]
        if item["CID"] != expected["cid"] or item.get("SMILES", item.get("IsomericSMILES")) != expected["smiles"]:
            raise ValueError(f"Unexpected PubChem record for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference hash changed")
    if (RAW / "PMC10520331-article.pdf").read_bytes()[:5] != b"%PDF-":
        raise ValueError("Downloaded article is not a PDF")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty artifact {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, sequences: dict[str, dict[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id in sorted(sequences):
            item = sequences[sequence_id]
            handle.write(f">{sequence_id} | {item['construct']} | {item['source']}\n")
            sequence = item["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def homology_hits(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    return {line.split("\t", 1)[0] for line in path.read_text(encoding="utf-8").splitlines() if line}


def write_outputs() -> None:
    verify_raw()
    sequences, full_by_code = build_sequences()
    records, exclusions = build_records(sequences)
    SOURCE.mkdir(parents=True, exist_ok=True)
    fasta_path = SOURCE / "construct_sequences.fasta"
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    hits = homology_hits(hit_path)
    for record in records:
        record["status_at_normalization"] = (
            "excluded_homology" if record["sequence_id"] in hits else "accepted_homology_cold_pool"
        )
    write_csv(SOURCE / "candidate_records.csv", records)
    write_csv(SOURCE / "exclusions.csv", exclusions)
    write_fasta(fasta_path, sequences)

    pdb_rows = []
    for code, meta in PARENTS.items():
        if not meta["pdb_id"]:
            continue
        full = full_by_code[code]
        pdb_rows.append(
            {
                "pdb_id": meta["pdb_id"],
                "parent_sequence_id": meta["sequence_id"],
                "accession": meta["accession"],
                "pdb_full_precursor_tagged_length": len(full),
                "pdb_full_precursor_tagged_sha256": hashlib.sha256(full.encode("ascii")).hexdigest(),
                "assay_construct_difference": "kinetics construct removes native signal peptide and prepends initiator Met",
                "raw_cif_sha256": sha256(RAW / f"{meta['pdb_id']}.cif"),
            }
        )
    write_csv(SOURCE / "pdb_construct_mapping.csv", pdb_rows)

    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    raw_files = sorted(
        [RAW / name for name in DOWNLOADS] + [RAW / "mmc1.pdf"], key=lambda path: path.name
    )
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC10520331",
        "article_doi": DOI,
        "pmc_id": "PMC10520331",
        "article_published": "2023-08-14",
        "article_and_supporting_information_license": "CC-BY-4.0",
        "pdb_coordinate_license": "CC0-1.0",
        "reported_finite_direct_single_activity_kcat_rows": len(DIRECT_ROWS),
        "saturation_qualified_rows": len(records),
        "accepted_after_frozen_homology": len(accepted),
        "saturation_exclusions": sum(row["exclusion_reason"] == "substrate_maximum_below_5x_km" for row in exclusions),
        "unique_exact_assayed_constructs": len(sequences),
        "selection_policy": "Retain finite directly fitted steady-state kcat for an isolated CM or CDT activity only when the 100 uM maximum varied substrate is at least 5x the fitted Km. Exclude coupled sequential CM+CDT apparent kcat, nonnumeric no-activity/not-measured cells, and duplicate controls.",
        "activity_boundary": {
            "CM": "EC 5.4.99.5; direct chorismate to prephenate assay; PubChem CID 12039",
            "CDT": "EC 4.2.1.51; direct prephenate to phenylpyruvate + water assay; PubChem CID 1028",
            "CM+CDT": "sequential chorismate to phenylpyruvate readout; excluded as an apparent two-site turnover endpoint",
        },
        "kinetics_sources": [
            "raw/PMC10520331-fullText.xml, Tables 1 and 3 and Experimental procedures",
            "raw/mmc1.pdf, Supporting Experimental Procedures and Table S4",
        ],
        "sequence_sources": [
            "raw/mmc1.pdf, Figure S2 expression DNA and construct-specific cloning primers",
            "raw/8CQ3.cif, raw/8CQ4.cif, raw/8CQ6.cif, exact full precursor/tagged parent checks",
        ],
        "construct_audit": "The kinetic proteins were cytoplasmic leaderless constructs. SI PCR primers define each replacement initiator Met and mature start; Figure S2 defines the codon-derived parent sequence and C-terminal LE-His6 tag. Active-site substitutions retain full-precursor numbering. Split constructs follow the exact SI primer boundaries and tag orientations.",
        "assay": {
            "type": "direct discontinuous initial-rate assay for each isolated activity",
            "buffer": "50 mM potassium phosphate pH 7.5, 0.1 mg/mL BSA",
            "temperature_C": 30,
            "substrate_concentrations_uM": [2.5, 5, 10, 25, 50, 100],
            "time_points": 4,
            "maximum_time_min": 4,
            "maximum_substrate_turnover_percent": 25,
            "background_correction_s-1": {"chorismate": 1.15e-5, "prephenate": 2.5e-5},
            "fit": "Michaelis-Menten in GraphPad Prism, separately for two biological replicates",
            "error": "mean and sample SD of the two independently fitted parameters",
        },
        "saturation_audit": "The SI states an ideal range extends 5-fold above Km. Af WT CDT (100/39.9=2.506) and split AfCDT (100/28.2=3.546) fail; all other direct rows reach at least 5x Km.",
        "low_concentration_audit": "The source italicizes parameters when 2.5 uM is not at least 2.5-fold below Km. This reliability flag is retained per record but is distinct from the high-substrate saturation gate.",
        "pubchem_structure_source": "PubChem PUG REST isomeric SMILES for chorismate CID 12039 and prephenate CID 1028, retrieved 2026-07-22",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "artifact_sha256": {
            "candidate_records.csv": sha256(SOURCE / "candidate_records.csv"),
            "exclusions.csv": sha256(SOURCE / "exclusions.csv"),
            "construct_sequences.fasta": sha256(fasta_path),
            "pdb_construct_mapping.csv": sha256(SOURCE / "pdb_construct_mapping.csv"),
        },
        "invalid_preexisting_raw_artifact": {
            "path": "raw/main.pdf",
            "reason": "HHS access-denial placeholder; not used",
            "sha256": sha256(RAW / "main.pdf") if (RAW / "main.pdf").is_file() else "",
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )
    write_homology_audit(records, sequences, fasta_path)
    verify_outputs(records, exclusions, sequences)


def write_homology_audit(
    records: list[dict[str, object]],
    sequences: dict[str, dict[str, str]],
    fasta_path: Path,
) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        return
    hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    cluster_rows = [line.split("\t")[:2] for line in cluster_path.read_text(encoding="utf-8").splitlines() if line]
    member_to_family = {member: representative for representative, member in cluster_rows}
    if set(member_to_family) != set(sequences):
        raise ValueError("MMseqs cluster output does not cover every exact construct")
    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    audit = {
        "audited_on": "2026-07-22",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": MMSEQS_VERSION,
        "min_identity": MIN_IDENTITY,
        "coverage": COVERAGE,
        "coverage_mode": COVERAGE_MODE,
        "development_target_sha256": REFERENCE_SHA256,
        "saturation_qualified_candidate_records": len(records),
        "unique_exact_assayed_construct_sequences": len(sequences),
        "construct_sequences_sha256": sha256(fasta_path),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
        "exact_sequence_overlap": sum(float(line.split("\t")[2]) == 1.0 for line in hit_lines),
        "homology_hits_sha256": sha256(hit_path),
        "candidate_mmseqs_clusters": len(set(member_to_family.values())),
        "cluster_counting_policy": "Count only distinct representatives emitted by frozen linclust across all exact assayed constructs; do not infer domain families or merge clusters manually.",
        "family_cluster_sha256": sha256(cluster_path),
        "accepted_records": len(accepted),
        "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
        "accepted_cm_records": sum(row["catalytic_activity"] == "CM" for row in accepted),
        "accepted_cdt_records": sum(row["catalytic_activity"] == "CDT" for row in accepted),
        "readiness_gate_passes": False,
        "claim_boundary": "Isolated curation artifacts only; no registry integration or model predictions.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )


def verify_outputs(
    records: list[dict[str, object]],
    exclusions: list[dict[str, object]],
    sequences: dict[str, dict[str, str]],
) -> None:
    if len(records) != 17 or len(sequences) != 13 or len(exclusions) != 23:
        raise ValueError(
            f"Unexpected cardinality: {len(records)} records, {len(sequences)} sequences, {len(exclusions)} exclusions"
        )
    if len({item["sequence"] for item in sequences.values()}) != 13:
        raise ValueError("Expected 13 distinct exact assay construct sequences")
    expected_activities = {"CM": 10, "CDT": 7}
    observed_activities = {
        activity: sum(row["catalytic_activity"] == activity for row in records)
        for activity in expected_activities
    }
    if observed_activities != expected_activities:
        raise ValueError(f"Activity boundary failed: {observed_activities}")
    for row in records:
        for field in ("kcat_s-1", "kcat_sd_s-1", "km_uM", "km_sd_uM"):
            value = float(row[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid {field} in {row['candidate_id']}")
        if not row["saturation_passes"] or float(row["maximum_substrate_multiple_of_km"]) < 5:
            raise ValueError(f"Unsaturated row leaked into candidates: {row['candidate_id']}")
        substrate = "chorismate" if row["catalytic_activity"] == "CM" else "prephenate"
        if int(row["substrate_pubchem_cid"]) != COMPOUNDS[substrate]["cid"]:
            raise ValueError(f"Activity/substrate mismatch in {row['candidate_id']}")
    if sum(row["exclusion_reason"] == "substrate_maximum_below_5x_km" for row in exclusions) != 2:
        raise ValueError("Expected exactly two direct saturation exclusions")
    if any(str(row["candidate_label_created"]).lower() not in ("false", "0") for row in exclusions):
        raise ValueError("An exclusion was marked as a candidate label")
    required = (
        SOURCE / "candidate_records.csv",
        SOURCE / "exclusions.csv",
        SOURCE / "construct_sequences.fasta",
        SOURCE / "pdb_construct_mapping.csv",
        SOURCE / "provenance.json",
    )
    if not all(path.is_file() and path.stat().st_size for path in required):
        raise ValueError("Required standard artifact missing")
    audit_path = SOURCE / "homology-audit.json"
    if audit_path.is_file():
        audit = json.loads(audit_path.read_text(encoding="ascii"))
        if audit["coverage_mode"] != 0 or audit["min_identity"] != 0.3 or audit["coverage"] != 0.8:
            raise ValueError("Frozen MMseqs parameters changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        download_raw()
    write_outputs()


if __name__ == "__main__":
    main()
