from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC13210114"
DOI = "10.3390/microorganisms14051070"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
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
DEFAULT_MMSEQS = "mmseqs"
NATIVE_ACCESSION = "AYM85572.1"
NATIVE_SEQUENCE = (
    "MYSLKAPVHKSARKGIILAGGSGTRLYPLTKVVSKQLMPVYDKPMIFYPVSTLMMAGITEILIISTPAELPRFKELLGDG"
    "SAWGISFEYVEQPSPDGLAQAFLLAEDFLQGQSAALVLGDNLFYGHDLSVSLQNATVCEYGATVFGYHVANPKSYGVVEF"
    "DENGKAISIEEKPDKPKSHYAVPGLYFFDNRVVEFAKNVKPSERGELEITDVIEQYLNNKELNVEIMGRGTAWLDTGTLD"
    "DLLDAANFIRAIEKRQGLKINCPEEIAYRMGYINAEELKKLAKPLKKSGYGKYLLSLLDQTVF"
)
SUBSTRATES = {
    "Glc-1-P": (
        65533,
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OP(=O)(O)O)O)O)O)O",
    ),
    "dTTP": (
        64968,
        "CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O",
    ),
}
KINETIC_ROWS = (
    ("rmla-glc1p", "Glc-1-P", 67.375, 3.305, 0.299, 0.054, 33.685, 1.665, "dTTP"),
    ("rmla-dttp", "dTTP", 42.845, 1.525, 0.025, 0.004, 171.4, 6.3, "Glc-1-P"),
)
URLS = {
    "PMC13210114-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13210114/fullTextXML"
    ),
    "PMC13210114-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13210114/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "GCA_003668795.1.zip": (
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        "GCA_003668795.1/download?include_annotation_type=PROT_FASTA"
        "&include_annotation_type=GENOME_GFF&include_annotation_type=GENOME_GBFF"
        "&include_annotation_type=SEQUENCE_REPORT"
    ),
    "pubchem-dTTP-64968.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/64968/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "pubchem-alpha-D-glucose-1-phosphate-65533.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/65533/"
        "property/IsomericSMILES,Title/JSON"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)

    supplement = RAW / "microorganisms-4240565-supplementary.pdf"
    if not supplement.is_file():
        with zipfile.ZipFile(RAW / "PMC13210114-supplementaryFiles.zip") as outer:
            nested_name = next(name for name in outer.namelist() if name.endswith("-s001.zip"))
            nested_path = RAW / "microorganisms-14-01070-s001.zip"
            nested_path.write_bytes(outer.read(nested_name))
        with zipfile.ZipFile(nested_path) as nested:
            pdf_name = next(name for name in nested.namelist() if name.endswith("supplementary.pdf"))
            supplement.write_bytes(nested.read(pdf_name))

    assembly_files = {
        "genomic.gbff": "ncbi_dataset/data/GCA_003668795.1/genomic.gbff",
        "genomic.gff": "ncbi_dataset/data/GCA_003668795.1/genomic.gff",
        "protein.faa": "ncbi_dataset/data/GCA_003668795.1/protein.faa",
        "sequence_report.jsonl": "ncbi_dataset/data/GCA_003668795.1/sequence_report.jsonl",
        "assembly_data_report.jsonl": "ncbi_dataset/data/assembly_data_report.jsonl",
    }
    with zipfile.ZipFile(RAW / "GCA_003668795.1.zip") as assembly:
        for output_name, member in assembly_files.items():
            output = RAW / output_name
            if not output.is_file():
                output.write_bytes(assembly.read(member))


def read_fasta_record(path: Path, accession: str) -> str:
    sequence: list[str] = []
    selected = False
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if selected:
                break
            selected = line[1:].split()[0] == accession
        elif selected:
            sequence.append(line.strip())
    return "".join(sequence)


def validate_raw() -> None:
    required = [*URLS, "microorganisms-14-01070-s001.zip", "microorganisms-4240565-supplementary.pdf", "genomic.gbff", "genomic.gff", "protein.faa", "sequence_report.jsonl", "assembly_data_report.jsonl"]
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {missing}; rerun with --acquire")

    xml = (RAW / "PMC13210114-fullText.xml").read_text(encoding="utf-8")
    for marker in (
        DOI,
        "creativecommons.org/licenses/by/4.0/",
        "GCA_003668795.1",
        "pET-22b",
        "pET16b-",
        "67.375",
        "33.685",
        "171.4",
        "excess Glc-1-P",
        "dTTP was provided in excess",
    ):
        if marker not in xml:
            raise ValueError(f"Article XML lacks expected evidence marker: {marker}")

    supplement_text = "\n".join(
        page.get_text() for page in fitz.open(RAW / "microorganisms-4240565-supplementary.pdf")
    )
    for marker in (
        "pET16b-PaRmlA",
        "CGCATATGATGTATTCGTTAAAAGCCCCAG",
        "CGCGGATCCTTAAAAAACAGTTTGATCTAAAA",
        "Table S2. Enzymatic Assay System of Pa-RmlA",
    ):
        if marker not in supplement_text.replace("\n", ""):
            raise ValueError(f"Supplement lacks expected evidence marker: {marker}")

    if read_fasta_record(RAW / "protein.faa", NATIVE_ACCESSION) != NATIVE_SEQUENCE:
        raise ValueError("GCA_003668795.1 AYM85572.1 sequence changed")
    gff = (RAW / "genomic.gff").read_text(encoding="utf-8")
    for marker in ("CP033065.1", "473759", "474670", "D9T18_02105", "AYM85572.1"):
        if marker not in gff:
            raise ValueError(f"Assembly mapping evidence missing: {marker}")
    for name, raw_name in (
        ("Glc-1-P", "pubchem-alpha-D-glucose-1-phosphate-65533.json"),
        ("dTTP", "pubchem-dTTP-64968.json"),
    ):
        compound = json.loads((RAW / raw_name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if compound["CID"] != SUBSTRATES[name][0] or compound["SMILES"] != SUBSTRATES[name][1]:
            raise ValueError(f"PubChem mapping changed for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development target is absent or changed")


def command(executable: str, *args: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *args]
    _, distribution, binary = executable.split(":", 2)
    return ["wsl", "-d", distribution, binary, *args]


def tool_path(path: Path, executable: str) -> str:
    if not executable.startswith("wsl:"):
        return str(path)
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def run_checked(args: list[str]) -> str:
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(args)}\n{result.stderr}")
    return (result.stdout or result.stderr).strip()


def run_mmseqs(query: Path, executable: str, threads: int) -> str:
    version = run_checked(command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs version {version!r} is not frozen version {MMSEQS_VERSION!r}")
    if HOMOLOGY.exists():
        if executable.startswith("wsl:"):
            run_checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(HOMOLOGY, executable)])
        else:
            shutil.rmtree(HOMOLOGY)
    HOMOLOGY.mkdir(parents=True)
    run_checked(
        command(
            executable,
            "easy-search",
            tool_path(query, executable),
            tool_path(REFERENCE, executable),
            tool_path(HOMOLOGY / "homology_hits.tsv", executable),
            tool_path(HOMOLOGY / "search-tmp", executable),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.8",
            "--cov-mode",
            "0",
            "--format-output",
            "query,target,fident,alnlen,qcov,tcov,evalue,bits",
            "--threads",
            str(threads),
        )
    )
    return version


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(executable: str, threads: int, run_homology: bool) -> None:
    validate_raw()
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "native_prescreen_sequence.fasta"
    query.write_text(
        ">AYM85572.1 | GCA_003668795.1 native CDS translation; NOT exact assay construct\n"
        + "\n".join(NATIVE_SEQUENCE[start : start + 80] for start in range(0, len(NATIVE_SEQUENCE), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    version = run_mmseqs(query, executable, threads) if run_homology else MMSEQS_VERSION
    hits_path = HOMOLOGY / "homology_hits.tsv"
    if not hits_path.is_file():
        raise RuntimeError("Pinned MMseqs evidence is required before curation output")
    hits = [line for line in hits_path.read_text(encoding="utf-8").splitlines() if line]

    candidate_fields = [
        "candidate_id", "article_doi", "source_file", "organism", "enzyme_identity",
        "sequence_id", "variable_substrate", "substrate_pubchem_cid", "endpoint",
        "value", "unit", "status_at_normalization",
    ]
    write_csv(SOURCE / "candidate_records.csv", [], candidate_fields)

    molecular_weight_kda = 33.6
    stated_mass_ug_ml = 1.67
    stated_enzyme_um = stated_mass_ug_ml / molecular_weight_kda
    exclusions: list[dict[str, object]] = []
    for row_id, substrate, vmax, vmax_sd, km, km_sd, kcat, kcat_sd, fixed in KINETIC_ROWS:
        implied_um = vmax / kcat / 60
        exclusions.append(
            {
                "evidence_id": row_id,
                "article_doi": DOI,
                "source_file": "raw/PMC13210114-fullText.xml",
                "source_section": "Table 1 and Results 3.5",
                "source_row": f"Table 1: {substrate}",
                "organism": "Pseudoalteromonas agarivorans Hao 2018",
                "enzyme_identity": "Pa-RmlA glucose-1-phosphate thymidylyltransferase",
                "native_sequence_accession": NATIVE_ACCESSION,
                "assay_construct_status": "unresolved",
                "variable_substrate": substrate,
                "substrate_pubchem_cid": SUBSTRATES[substrate][0],
                "substrate_isomeric_smiles": SUBSTRATES[substrate][1],
                "reported_vmax": vmax,
                "reported_vmax_sd": vmax_sd,
                "reported_vmax_unit": "uM/min",
                "reported_km": km,
                "reported_km_sd": km_sd,
                "reported_km_unit": "mM",
                "reported_kcat": kcat,
                "reported_kcat_sd": kcat_sd,
                "reported_kcat_unit": "s-1",
                "fixed_cosubstrate": fixed,
                "fixed_cosubstrate_concentration": "",
                "fixed_cosubstrate_saturation_evidence": "author says excess; concentration and saturation series absent",
                "assay_temperature_C": 50,
                "assay_pH": 8,
                "assay_Mg2_mM": 10,
                "assay_time_min": 3,
                "stated_enzyme_mass_concentration_ug_mL": stated_mass_ug_ml,
                "reported_protein_molecular_weight_kDa": molecular_weight_kda,
                "derived_stated_enzyme_concentration_uM": round(stated_enzyme_um, 9),
                "diagnostic_kcat_from_vmax_and_stated_enzyme_s-1": round(vmax / stated_enzyme_um / 60, 9),
                "diagnostic_implied_enzyme_concentration_uM": round(implied_um, 9),
                "diagnostic_implied_enzyme_mass_ug_mL": round(implied_um * molecular_weight_kda, 9),
                "status_at_normalization": "excluded_fail_closed",
                "exclusion_codes": "DEVELOPMENT_HOMOLOGY_HIT;EXACT_ASSAY_CONSTRUCT_UNRESOLVED;FIXED_COSUBSTRATE_SATURATION_UNRESOLVED;REPORTED_KCAT_VMAX_ENZYME_CONCENTRATION_INCONSISTENT",
                "candidate_label_created": False,
            }
        )
    exclusion_fields = list(exclusions[0])
    write_csv(SOURCE / "exclusions.csv", exclusions, exclusion_fields)

    sequence_mapping = [{
        "assembly": "GCA_003668795.1",
        "replicon": "CP033065.1",
        "coordinates_1_based_inclusive": "473759-474670",
        "strand": "+",
        "gene": "rfbA",
        "locus_tag": "D9T18_02105",
        "protein_accession": NATIVE_ACCESSION,
        "native_length_aa": len(NATIVE_SEQUENCE),
        "native_sequence_sha256": hashlib.sha256(NATIVE_SEQUENCE.encode("ascii")).hexdigest(),
        "primer_forward": "CGCATATGATGTATTCGTTAAAAGCCCCAG",
        "primer_reverse": "CGCGGATCCTTAAAAAACAGTTTGATCTAAAAG",
        "amplicon_mapping": "full native ORF; NdeI/BamHI primers; forward primer has NdeI ATG followed by native ATG; reverse primer encodes a stop codon",
        "assay_construct_mapping": "unresolved; native sequence is prescreen evidence only",
    }]
    write_csv(SOURCE / "sequence_mapping.csv", sequence_mapping, list(sequence_mapping[0]))

    construct = {
        "status": "unresolved_fail_closed",
        "native_sequence": {"accession": NATIVE_ACCESSION, "length_aa": 303, "sequence_file": "native_prescreen_sequence.fasta"},
        "consistent_evidence": [
            "Table S1 primers amplify the complete native AYM85572.1 ORF through NdeI and BamHI.",
            "The forward primer is CGC-CATATG-ATGTAT..., retaining the native ATG after the NdeI start codon.",
            "The RmlA reverse primer contains a stop codon before BamHI.",
            "Protein was purified by Ni-NTA and the article reports an apparent mass of 33.6 kDa.",
        ],
        "contradictions": [
            "Methods 2.3 says pET-22b and a C-terminal 6xHis tag.",
            "Methods 2.3 later says pET16b-PaRmlA; Supplement Figure S1 and its caption also say pET16b-PaRmlA.",
            "The RmlA reverse primer's stop codon prevents the claimed vector-encoded C-terminal His fusion.",
            "The duplicated start codon context and unspecified junction sequence prevent an exact N-terminal fusion reconstruction.",
            "No deposited plasmid or junction sequencing is available; the data statement says underlying data are available only on request.",
        ],
        "decision": "No exact assayed fusion sequence can be reconstructed without choosing among contradictory source statements. The native 303-aa sequence is not substituted for the assay construct.",
        "model_predictions_run": False,
    }
    (SOURCE / "construct-resolution.json").write_text(json.dumps(construct, indent=2) + "\n", encoding="ascii")

    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path), "url": URLS.get(path.name, "derived from a downloaded archive")}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    write_csv(SOURCE / "homology_hits.csv", [dict(zip(("query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"), line.split("\t"))) for line in hits], ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"])

    audit = {
        "audited_on": "2026-07-27",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search native-sequence prescreen; exact assay construct remained unresolved",
        "mmseqs_version": version,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "search_command": "mmseqs easy-search native_prescreen_sequence.fasta unikp_reference.fasta homology_hits.tsv search-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": {"count": 0, "sha256": sha256(SOURCE / "candidate_records.csv")},
        "reported_numerical_rows_retained_as_exclusions": len(exclusions),
        "homology_queries_sha256": sha256(query),
        "homology_query_scope": "AYM85572.1 native sequence only; sufficient to exclude, not asserted as exact assay construct",
        "homology_hit_sequences": 1 if hits else 0,
        "homology_hit_alignments": len(hits),
        "homology_hits_sha256": sha256(hits_path),
        "accepted_records": 0,
        "status": "excluded_fail_closed",
        "exclusion_codes": ["DEVELOPMENT_HOMOLOGY_HIT", "EXACT_ASSAY_CONSTRUCT_UNRESOLVED", "FIXED_COSUBSTRATE_SATURATION_UNRESOLVED", "REPORTED_KCAT_VMAX_ENZYME_CONCENTRATION_INCONSISTENT"],
        "exact_sequence_substrate_overlap": "not evaluated after mandatory homology exclusion",
        "readiness_gate_passes": False,
        "claim_boundary": "Curation/exclusion evidence only; no recalculated labels and no model predictions.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")

    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC13210114",
        "article_doi": DOI,
        "article_published": "2026-05-09",
        "license": "CC-BY-4.0",
        "license_scope": "Article and embedded supporting information",
        "reported_direct_kcat_rows": 2,
        "candidate_records": 0,
        "accepted_records": 0,
        "kinetics_source": "raw/PMC13210114-fullText.xml, Methods 2.7, Results 3.5, Figure 7, Table 1",
        "supplement_source": "raw/microorganisms-4240565-supplementary.pdf, Figure S1 and Tables S1-S2",
        "sequence_source": "raw GCA_003668795.1 package: genomic.gbff, genomic.gff, protein.faa",
        "construct_resolution": "construct-resolution.json",
        "sequence_mapping": "sequence_mapping.csv",
        "substrate_structures": {name: {"pubchem_cid": value[0], "isomeric_smiles": value[1]} for name, value in SUBSTRATES.items()},
        "saturation_audit": [
            {"variable": "Glc-1-P", "fixed": "dTTP", "fixed_concentration": None, "status": "unresolved", "source_wording": "dTTP was provided in excess"},
            {"variable": "dTTP", "fixed": "Glc-1-P", "fixed_concentration": None, "status": "unresolved", "source_wording": "excess Glc-1-P ensured"},
        ],
        "arithmetic_audit": {
            "stated_enzyme_mass_concentration_ug_mL": stated_mass_ug_ml,
            "reported_molecular_weight_kDa": molecular_weight_kda,
            "derived_stated_enzyme_concentration_uM": stated_enzyme_um,
            "formula": "kcat = Vmax_uM_per_min / enzyme_uM / 60",
            "decision": "Reported kcat values are retained verbatim only as excluded evidence. Diagnostic recalculations are not labels.",
        },
        "homology_status": "excluded; three frozen-threshold development hits to the native sequence",
        "final_disposition": "excluded_fail_closed",
        "raw_file_hashes": "raw-file-hashes.json",
        "model_predictions_run": False,
        "recalculated_labels_created": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    parser.add_argument("--run-mmseqs", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.acquire:
        download_raw()
    write_outputs(args.mmseqs, args.threads, args.run_mmseqs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
