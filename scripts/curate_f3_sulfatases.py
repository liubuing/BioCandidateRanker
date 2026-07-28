from __future__ import annotations

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
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "europepmc-PMC11659886"
)
RAW = SOURCE / "raw"
REFERENCE = (
    ROOT
    / "artifacts"
    / "external"
    / "absolute-kinetics-screen"
    / "dryad-4964723"
    / "homology"
    / "unikp_reference.fasta"
)
DOI = "10.1021/acsbiomedchemau.4c00088"
SOURCE_ID = "europepmc-PMC11659886"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
FAMILY_CAP = 20

SUBSTRATES = {
    "4MU-6S-GlcNAc": (
        18740504,
        "CC1=CC(=O)OC2=C1C=CC(=C2)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)COS(=O)(=O)O)O)O)NC(=O)C",
        "pubchem-4MU-6S-GlcNAc.json",
    ),
    "4MU-GlcNAc": (
        2733787,
        "CC1=CC(=O)OC2=C1C=CC(=C2)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)NC(=O)C",
        "pubchem-4MU-GlcNAc.json",
    ),
}

ENZYMES = {
    "WG": {
        "uniprot": "R6ARV4",
        "genbank": "CDA43927.1",
        "organism": "Prevotella sp. CAG:5226",
        "pH": 6.0,
        "variants": {
            "WT": (),
            "W437F": ((437, "W", "F"),),
            "W437A": ((437, "W", "A"),),
            "W437Q": ((437, "W", "Q"),),
            "G438I": ((438, "G", "I"),),
            "N443D": ((443, "N", "D"),),
            "R444A": ((444, "R", "A"),),
        },
    },
    "YG": {
        "uniprot": "L1MPC2",
        "genbank": "EKX92880.1",
        "organism": "Alloprevotella sp. oral taxon 473 str. F0040",
        "pH": 6.0,
        "variants": {
            "WT": (),
            "Y439F": ((439, "Y", "F"),),
            "Y439A": ((439, "Y", "A"),),
            "Y439Q": ((439, "Y", "Q"),),
            "G440I": ((440, "G", "I"),),
            "N445D": ((445, "N", "D"),),
            "R446A": ((446, "R", "A"),),
        },
    },
    "F3-ORF26": {
        "uniprot": "A0A4R4I8J5",
        "genbank": "GKH82764.1",
        "organism": "Phocaeicola dorei",
        "pH": 6.0,
        "variants": {
            "WT": (),
            "Q433E": ((433, "Q", "E"),),
            "S438A": ((438, "S", "A"),),
            "N439D": ((439, "N", "D"),),
            "Y442F": ((442, "Y", "F"),),
            "N443D": ((443, "N", "D"),),
            "R444A": ((444, "R", "A"),),
        },
    },
}

# Protein, variant, table label, substrate, Km +/- uncertainty (uM), kcat +/-
# uncertainty (s-1), enzyme nM, substrate minimum/maximum (uM), SI figure.
ROWS = (
    ("WG", "WT", "wild-type (4MU-6S-GlcNAc)", "4MU-6S-GlcNAc", 5.1, 0.3, 39.9, 0.6, 0.5, 1, 150, "Figure 4c"),
    ("WG", "WT", "wild-type (4MU-GlcNAc)", "4MU-GlcNAc", 2706, 344, 11.6, 0.5, 40, 100, 10000, "Figure 4d"),
    ("WG", "W437F", "W437F", "4MU-6S-GlcNAc", 8.0, 0.8, 63, 2, 0.5, 1.25, 200, "Figure S7a"),
    ("WG", "W437A", "W437A", "4MU-6S-GlcNAc", 17, 2, 45, 1, 2, 2.5, 300, "Figure S7b"),
    ("WG", "W437Q", "W437Q", "4MU-6S-GlcNAc", 20, 3, 48, 2, 2, 2.5, 300, "Figure S7c"),
    ("WG", "G438I", "G438I", "4MU-6S-GlcNAc", 32, 3, 39, 1, 2, 2.5, 500, "Figure S7d"),
    ("WG", "N443D", "N443D", "4MU-6S-GlcNAc", 56, 6, 69, 4, 1, 5, 750, "Figure S7e"),
    ("WG", "R444A", "R444A", "4MU-6S-GlcNAc", 47, 6, 62, 2, 2, 5, 800, "Figure S7f"),
    ("YG", "WT", "wild-type (4MU-6S-GlcNAc)", "4MU-6S-GlcNAc", 12, 1, 76, 2, 2, 2.5, 300, "Figure 4h"),
    ("YG", "WT", "wild-type (4MU-GlcNAc)", "4MU-GlcNAc", 2209, 355, 11.7, 0.7, 31.8, 100, 10000, "Figure 4i"),
    ("YG", "Y439F", "Y439F", "4MU-6S-GlcNAc", 9.6, 0.7, 69, 1, 0.46, 1.25, 200, "Figure S8a"),
    ("YG", "Y439A", "Y439A", "4MU-6S-GlcNAc", 10.7, 0.8, 30.2, 0.5, 2, 2.5, 300, "Figure S8b"),
    ("YG", "Y439Q", "Y439Q", "4MU-6S-GlcNAc", 24, 2, 44, 1, 2, 2.5, 300, "Figure S8c"),
    ("YG", "G440I", "G440I", "4MU-6S-GlcNAc", 22, 1, 42.8, 0.7, 2, 2.5, 300, "Figure S8d"),
    ("YG", "N445D", "N445D", "4MU-6S-GlcNAc", 130, 14, 98, 3, 2, 2.5, 300, "Figure S8e"),
    ("YG", "R446A", "R446A", "4MU-6S-GlcNAc", 90, 9, 97, 3, 2, 5, 800, "Figure S8f"),
    ("F3-ORF26", "WT", "wild-type (4MU-6S-GlcNAc)", "4MU-6S-GlcNAc", 31, 3, 61, 2, 1, 5, 1000, "Figure 1d"),
    ("F3-ORF26", "WT", "wild-type (4MU-GlcNAc)", "4MU-GlcNAc", 985, 247, 6.9, 0.7, 100, 100, 8000, "Figure 1e"),
    ("F3-ORF26", "Q433E", "Q443E", "4MU-6S-GlcNAc", 33, 2, 64.1, 0.8, 2, 10, 1000, "Figure S2a"),
    ("F3-ORF26", "S438A", "S438A", "4MU-6S-GlcNAc", 944, 94, 98, 4, 2, 5, 3000, "Figure S2b"),
    ("F3-ORF26", "N439D", "N439D", "4MU-6S-GlcNAc", 1987, 195, 40, 2, 4, 10, 5000, "Figure S2c"),
    ("F3-ORF26", "Y442F", "Y442F", "4MU-6S-GlcNAc", 33, 2, 69, 1, 2, 10, 1000, "Figure S2d"),
    ("F3-ORF26", "N443D", "N443D", "4MU-6S-GlcNAc", 139, 16, 64, 2, 2, 10, 1500, "Figure S2e"),
    ("F3-ORF26", "R444A", "R444A", "4MU-6S-GlcNAc", 50, 6, 82, 3, 2, 5, 750, "Figure S2f"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def read_latest_unisave_sequence(path: Path, accession: str) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or f"AC   {accession};" not in lines:
        raise ValueError(f"UniSave record does not identify {accession}")
    start = next(i for i, line in enumerate(lines) if line.startswith("SQ   SEQUENCE"))
    sequence_lines = []
    for line in lines[start + 1 :]:
        if line == "//":
            break
        if line.startswith("     "):
            sequence_lines.append(re.sub(r"[^A-Z]", "", line))
    sequence = "".join(sequence_lines)
    if not sequence or not sequence.isalpha() or not sequence.isupper():
        raise ValueError(f"Could not recover first archived {accession} sequence")
    return sequence


def mutate(sequence: str, changes: tuple[tuple[int, str, str], ...], accession: str) -> str:
    residues = list(sequence)
    for position, expected, replacement in changes:
        if residues[position - 1] != expected:
            raise ValueError(
                f"{accession} expected {expected}{position}, found {residues[position - 1]}"
            )
        residues[position - 1] = replacement
    return "".join(residues)


def make_construct(enzyme: str, native: str) -> str:
    if enzyme == "WG":
        return "M" + native[21:] + "LEHHHHHH"
    if enzyme == "YG":
        # SI sequence uses native 23-1284, then vector-derived TLKVL-E-His6.
        return "M" + native[22:1284] + "TLKVLEHHHHHH"
    if enzyme == "F3-ORF26":
        return "M" + native[22:] + "HHHHHH"
    raise ValueError(f"Unknown enzyme {enzyme}")


def sequence_id(enzyme: str, variant: str) -> str:
    prefix = {"WG": "wg", "YG": "yg", "F3-ORF26": "f3-orf26"}[enzyme]
    return f"{prefix}-{variant.lower()}"


def build_sequences() -> tuple[dict[str, tuple[str, str, str]], dict[str, str]]:
    native = {
        str(data["uniprot"]): read_latest_unisave_sequence(
            RAW / f"{data['uniprot']}.unisave.txt", str(data["uniprot"])
        )
        for data in ENZYMES.values()
    }
    if len(native["A0A4R4I8J5"]) != 773:
        raise ValueError("Expected 773-aa A0A4R4I8J5")
    if native["R6ARV4"] != read_fasta(RAW / "CDA43927.1.fasta"):
        raise ValueError("R6ARV4 does not exactly match CDA43927.1")
    if native["L1MPC2"] != read_fasta(RAW / "EKX92880.1.fasta"):
        raise ValueError("L1MPC2 does not exactly match EKX92880.1")

    constructs = {}
    for enzyme, data in ENZYMES.items():
        accession = str(data["uniprot"])
        for variant, changes in data["variants"].items():
            variant_native = mutate(native[accession], changes, accession)
            constructs[sequence_id(enzyme, variant)] = (
                enzyme,
                variant,
                make_construct(enzyme, variant_native),
            )
    return constructs, native


def parse_measurement(text: str) -> tuple[float, float]:
    parts = re.split(r"\s*±\s*", text.replace("\xa0", " ").strip())
    if len(parts) != 2:
        raise ValueError(f"Expected value +/- uncertainty, found {text!r}")
    return float(parts[0]), float(parts[1])


def table_measurements(article: ET.Element) -> dict[tuple[str, str], tuple[float, ...]]:
    table = article.find(".//table-wrap[@id='tbl1']")
    if table is None:
        raise ValueError("Article Table 1 is missing")
    observed = {}
    protein = ""
    for row in table.findall(".//tbody/tr"):
        cells = [" ".join("".join(cell.itertext()).split()) for cell in row]
        if cells[0].startswith("Bt4394"):
            protein = "Bt4394"
            cells.pop(0)
        elif cells[0] in {"WG", "YG", "F3-ORF26"}:
            protein = cells.pop(0)
        label = cells[0]
        km, km_error = parse_measurement(cells[1])
        kcat, kcat_error = parse_measurement(cells[2])
        observed[(protein, label)] = (km, km_error, kcat, kcat_error)
    return observed


def verify_raw(constructs: dict[str, tuple[str, str, str]]) -> None:
    required = (
        "PMC11659886-fullText.xml",
        "PMC11659886-supplementaryFiles.zip",
        "PMC11659886-main.pdf",
        "bg4c00088_si_001.pdf",
        "A0A4R4I8J5.unisave.txt",
        "R6ARV4.unisave.txt",
        "CDA43927.1.fasta",
        "L1MPC2.unisave.txt",
        "EKX92880.1.fasta",
        "pubchem-4MU-6S-GlcNAc.json",
        "pubchem-4MU-GlcNAc.json",
    )
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw originals: {', '.join(missing)}")
    for name in ("PMC11659886-main.pdf", "bg4c00088_si_001.pdf"):
        if (RAW / name).read_bytes()[:5] != b"%PDF-":
            raise ValueError(f"{name} is not a PDF")

    article = ET.parse(RAW / "PMC11659886-fullText.xml").getroot()
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected article DOI")
    permissions = " ".join(article.find(".//permissions").itertext())
    if "creativecommons.org/licenses/by/4.0" not in permissions:
        raise ValueError("Article is not verified CC-BY-4.0")
    observed = table_measurements(article)
    for enzyme, _, label, _, km, km_error, kcat, kcat_error, *_ in ROWS:
        expected = (float(km), float(km_error), float(kcat), float(kcat_error))
        if observed.get((enzyme, label)) != expected:
            raise ValueError(f"Table 1 mismatch for {enzyme} {label}")

    with fitz.open(RAW / "bg4c00088_si_001.pdf") as document:
        si_text = "".join(page.get_text() for page in document)
    compact_si = re.sub(r"[^A-Z]", "", si_text)
    for sequence_id_value in ("wg-wt", "yg-wt", "f3-orf26-wt"):
        sequence = constructs[sequence_id_value][2]
        if sequence not in compact_si:
            raise ValueError(f"Exact SI construct not recovered for {sequence_id_value}")
    for phrase in (
        "Figure S2. Michaelis-Menten plots",
        "Figure S7. Michaelis-Menten plots",
        "Figure S8. Michaelis-Menten plots",
        "All experiments were performed in triplicate",
    ):
        if phrase not in si_text and phrase not in " ".join(article.itertext()):
            raise ValueError(f"Required kinetic evidence is missing: {phrase}")
    with zipfile.ZipFile(RAW / "PMC11659886-supplementaryFiles.zip") as archive:
        if "bg4c00088_si_001.pdf" not in archive.namelist():
            raise ValueError("PMC supplementary archive lacks the SI PDF")

    for name, (cid, smiles, filename) in SUBSTRATES.items():
        payload = json.loads((RAW / filename).read_text(encoding="utf-8"))
        compound = payload["PropertyTable"]["Properties"][0]
        if compound["CID"] != cid or compound["SMILES"] != smiles:
            raise ValueError(f"Unexpected PubChem identity for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP development reference hash changed")


def build_records() -> list[dict[str, object]]:
    records = []
    for index, row in enumerate(ROWS, 1):
        (
            enzyme, variant, source_label, substrate, km, km_error, kcat,
            kcat_error, enzyme_nm, substrate_min, substrate_max, source_figure,
        ) = row
        data = ENZYMES[enzyme]
        cid, smiles, _ = SUBSTRATES[substrate]
        changes = data["variants"][variant]
        records.append(
            {
                "candidate_id": f"f3-sulfoglycosidase-{index:03d}",
                "article_doi": DOI,
                "source_table": "Table 1",
                "source_figure": source_figure,
                "source_row": f"{enzyme} {source_label}",
                "organism": data["organism"],
                "sequence_accession": data["uniprot"],
                "cross_reference_accession": data["genbank"],
                "sequence_id": sequence_id(enzyme, variant),
                "enzyme_name": enzyme,
                "enzyme_variant": variant,
                "amino_acid_changes": "/".join(
                    f"{old}{position}{new}" for position, old, new in changes
                ) or "WT",
                "construct": {
                    "WG": "SI-exact M + R6ARV4/CDA43927.1 residues 22-1284 + LE-His6",
                    "YG": "SI-exact M + L1MPC2/EKX92880.1 residues 23-1284 + TLKVL-E-His6",
                    "F3-ORF26": "SI-exact M + A0A4R4I8J5 residues 23-773 + His6",
                }[enzyme],
                "variable_substrate": substrate,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "reaction_product": "4-methylumbelliferone + 6-sulfo-GlcNAc" if substrate == "4MU-6S-GlcNAc" else "4-methylumbelliferone + GlcNAc",
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_uncertainty_s-1": kcat_error,
                "km_uM": km,
                "km_uncertainty_uM": km_error,
                "assay_method": "direct continuous 4MU fluorescence initial-rate Michaelis-Menten fit",
                "assay_buffer": "25 mM Bis-tris propane, 25 mM citrate, 300 mM NaCl",
                "assay_pH": data["pH"],
                "assay_temperature_C": 25,
                "assay_volume_uL": 100,
                "enzyme_concentration_nM": enzyme_nm,
                "substrate_min_uM": substrate_min,
                "substrate_max_uM": substrate_max,
                "excitation_nm": 360,
                "emission_nm": 450,
                "replicates": 3,
                "uncertainty_type": "author-reported +/-; statistic not defined",
                "fit": "Michaelis-Menten equation in GraphPad Prism 6.01",
                "source_nomenclature_note": "Table/SI label Q443E corrected to sequence- and text-consistent Q433E" if enzyme == "F3-ORF26" and variant == "Q433E" else "",
                "status_at_normalization": "pending_homology",
            }
        )
    return records


def read_hits(path: Path) -> set[str]:
    if not path.is_file():
        return set()
    return {
        line.split("\t", 1)[0]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    }


def apply_selection(records: list[dict[str, object]]) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        return
    hits = read_hits(hit_path)
    clusters = {}
    for line in cluster_path.read_text(encoding="utf-8").splitlines():
        if line:
            representative, member = line.split("\t")[:2]
            clusters[member] = representative
    family_counts: dict[str, int] = {}
    for record in records:
        identifier = str(record["sequence_id"])
        if identifier in hits:
            record["status_at_normalization"] = "excluded_homology"
            continue
        family = clusters.get(identifier, identifier)
        if family_counts.get(family, 0) >= FAMILY_CAP:
            record["status_at_normalization"] = "excluded_family_cap"
        else:
            record["status_at_normalization"] = "accepted_homology_cold_pool"
            family_counts[family] = family_counts.get(family, 0) + 1


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, entries: dict[str, tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, (description, sequence) in entries.items():
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def exclusion_records() -> list[dict[str, object]]:
    return [
        {
            "source_rows": "Table 1 Bt4394 WT (two substrates), Q431W/I432G (two substrates), and Q431Y/I432G (two substrates)",
            "count": 6,
            "reason": "comparator_outside_requested_A0A4R4I8J5_R6ARV4_CDA43927.1_L1MPC2_EKX92880.1_scope",
            "candidate_label_created": False,
        },
        {
            "source_rows": "Figures 1b, 4b, and 4g pH-rate profiles",
            "count": 3,
            "reason": "fixed_substrate_activity_profiles_not_saturation_kcat",
            "candidate_label_created": False,
        },
        {
            "source_rows": "SI Figure S14 truncated YG constructs",
            "count": 3,
            "reason": "insoluble_unproductive_constructs_no_finite_kinetics",
            "candidate_label_created": False,
        },
    ]


def write_outputs() -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    constructs, native = build_sequences()
    verify_raw(constructs)
    records = build_records()
    apply_selection(records)

    construct_path = SOURCE / "construct_sequences.fasta"
    write_fasta(
        construct_path,
        {
            identifier: (
                f"{ENZYMES[enzyme]['uniprot']} {ENZYMES[enzyme]['genbank']} | SI-exact tagged construct | {variant}",
                sequence,
            )
            for identifier, (enzyme, variant, sequence) in constructs.items()
        },
    )
    reference_path = SOURCE / "family_reference_sequences.fasta"
    write_fasta(
        reference_path,
        {
            accession: ("archived UniProt sequence version 1", sequence)
            for accession, sequence in native.items()
        }
        | {
            "CDA43927.1": ("exact GenBank sequence; identical to R6ARV4", read_fasta(RAW / "CDA43927.1.fasta")),
            "EKX92880.1": ("exact GenBank sequence; identical to L1MPC2", read_fasta(RAW / "EKX92880.1.fasta")),
        },
    )
    write_csv(SOURCE / "candidate_records.csv", records)
    write_csv(SOURCE / "exclusions.csv", exclusion_records())

    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    ]
    raw_files = sorted(path for path in RAW.iterdir() if path.is_file())
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC11659886",
        "article_doi": DOI,
        "pmc_id": "PMC11659886",
        "article_published": "2024-11-19",
        "article_and_supplement_license": "CC-BY-4.0",
        "uniprot_unisave_license": "CC-BY-4.0",
        "genbank_notice": "NCBI sequence downloads retained verbatim",
        "pubchem_notice": "PubChem PUG REST records retained verbatim",
        "reported_requested_direct_finite_saturation_kcat_records": len(records),
        "accepted_records": len(accepted),
        "unique_exact_assayed_construct_sequences": len(constructs),
        "unique_substrates": len(SUBSTRATES),
        "kinetics_sources": [
            "raw/PMC11659886-fullText.xml, Table 1, Figures 1 and 4, and Michaelis-Menten Kinetics methods",
            "raw/bg4c00088_si_001.pdf, exact constructs and Figures S2, S7, and S8",
            "raw/PMC11659886-main.pdf, CC-BY Europe PMC rendered article PDF",
        ],
        "sequence_sources": [
            "raw/A0A4R4I8J5.unisave.txt (archived sequence version 1)",
            "raw/R6ARV4.unisave.txt and raw/CDA43927.1.fasta (identical sequence)",
            "raw/L1MPC2.unisave.txt and raw/EKX92880.1.fasta (identical sequence)",
            "raw/bg4c00088_si_001.pdf (complete tagged assayed constructs)",
        ],
        "construct_audit": {
            "F3-ORF26": "SI exact: initiator M + A0A4R4I8J5 residues 23-773 + C-terminal His6; all six variants use accession coordinates.",
            "WG": "SI exact: initiator M + R6ARV4/CDA43927.1 residues 22-1284 + LE-His6; all six variants use accession coordinates.",
            "YG": "SI exact: initiator M + L1MPC2/EKX92880.1 residues 23-1284 + TLKVL-E-His6. This differs from the methods label YG(1-1232) and from the accession terminal KP; the explicit SI protein sequence is retained without inference.",
            "F3_Q433E_nomenclature": "Table 1 and SI primer name say Q443E, but A0A4R4I8J5 residue 443 is N. Main text identifies Q433E, residue 433 is Q, and the primer sequence maps to Q433. Curated variant is therefore Q433E while source label is retained.",
        },
        "structure_source": "raw/pubchem-4MU-6S-GlcNAc.json and raw/pubchem-4MU-GlcNAc.json; PubChem CIDs 18740504 and 2733787 retrieved 2026-07-22",
        "assay": {
            "method": "continuous fluorescent 4MU initial rates fitted directly to Michaelis-Menten",
            "buffer": "25 mM Bis-tris propane, 25 mM citrate, 300 mM NaCl, pH 6.0",
            "temperature_C": 25,
            "volume_uL": 100,
            "excitation_nm": 360,
            "emission_nm": 450,
            "replicates": 3,
            "uncertainty": "author-reported +/- values; statistic is not defined in article or SI",
        },
        "selection_policy": "All and only finite direct substrate-saturation kcat rows in Table 1 for the explicitly requested A0A4R4I8J5, R6ARV4/CDA43927.1, and L1MPC2/EKX92880.1 enzymes and variants; no graph digitization, refitting, or inferred labels.",
        "excluded_or_skipped": [
            "Six finite Bt4394 comparator rows were outside the explicitly requested accession scope and are recorded in exclusions.csv.",
            "Fixed-substrate pH profiles and relative activities are not saturation kcat endpoints.",
            "Truncated YG constructs were insoluble and yielded no finite kinetics.",
            "Current UniProt FASTA endpoints for retired A0A4R4I8J5, R6ARV4, and L1MPC2 returned empty files; exact archived version-1 sequences were recovered from CC-BY UniSave records instead.",
            "The NCBI OA API advertised ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/4d/fc/PMC11659886.tar.gz, but the HTTPS endpoint returned 404 on 2026-07-22; the Europe PMC CC-BY XML, rendered main PDF, and supplementary archive/PDF were retained instead.",
            "raw/bg4c00088.pdf is a 1,817-byte failed ACS response, not a PDF, and is excluded from evidentiary originals.",
        ],
        "normalization": "Author-reported kcat already in s-1 and Km in uM; values and uncertainties transcribed without conversion or refitting.",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "artifact_sha256": {
            "candidate_records.csv": sha256(SOURCE / "candidate_records.csv"),
            "exclusions.csv": sha256(SOURCE / "exclusions.csv"),
            "construct_sequences.fasta": sha256(construct_path),
            "family_reference_sequences.fasta": sha256(reference_path),
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )
    write_homology_audit(records, constructs, construct_path, reference_path)
    verify_outputs(records, constructs)


def write_homology_audit(
    records: list[dict[str, object]],
    constructs: dict[str, tuple[str, str, str]],
    construct_path: Path,
    reference_path: Path,
) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        return
    hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    families = {
        line.split("\t", 1)[0]
        for line in cluster_path.read_text(encoding="utf-8").splitlines()
        if line
    }
    accepted = [
        record
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
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
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": len(records),
        "unique_assayed_construct_sequences": len(constructs),
        "construct_sequences_sha256": sha256(construct_path),
        "family_reference_sequences_sha256": sha256(reference_path),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hits}),
        "exact_sequence_overlap": sum(
            len(parts) > 2 and float(parts[2]) == 1.0
            for parts in (line.split("\t") for line in hits)
        ),
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


def verify_outputs(
    records: list[dict[str, object]], constructs: dict[str, tuple[str, str, str]]
) -> None:
    if len(records) != 24 or len(constructs) != 21:
        raise ValueError("Unexpected curation cardinality")
    if len({sequence for _, _, sequence in constructs.values()}) != 21:
        raise ValueError("Expected 21 distinct exact assayed constructs")
    for record in records:
        for field in ("kcat_s-1", "kcat_uncertainty_s-1", "km_uM", "km_uncertainty_uM"):
            value = float(record[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid {field} in {record['candidate_id']}")
    with (SOURCE / "exclusions.csv").open(newline="", encoding="utf-8") as handle:
        exclusions = list(csv.DictReader(handle))
    if len(exclusions) != 3 or any(row["candidate_label_created"] != "False" for row in exclusions):
        raise ValueError("Excluded or unproductive rows leaked into candidate labels")

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.is_file() or cluster_path.is_file():
        if not hit_path.is_file() or not cluster_path.is_file():
            raise ValueError("Homology output is incomplete")
        if any(record["status_at_normalization"] == "pending_homology" for record in records):
            raise ValueError("Frozen homology artifacts were not applied")
        required = (
            SOURCE / "candidate_records.csv",
            SOURCE / "construct_sequences.fasta",
            SOURCE / "family_reference_sequences.fasta",
            SOURCE / "exclusions.csv",
            SOURCE / "provenance.json",
            SOURCE / "homology-audit.json",
            cluster_path,
        )
        if not all(path.is_file() and path.stat().st_size for path in required):
            raise ValueError("A required standard artifact is missing or empty")


if __name__ == "__main__":
    write_outputs()
