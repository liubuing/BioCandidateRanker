from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC13223937"
DOI = "10.1016/j.jbc.2026.113078"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"
DEVELOPMENT_CORPUS_SHA256 = "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719"
DEVELOPMENT_CORPUS_SIZE = 12132719
DEFAULT_DEVELOPMENT_CORPUS = Path(
    r"D:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment"
    r"\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json"
)

URLS = {
    "PMC13223937-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13223937/fullTextXML",
    "PMC13223937-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13223937/supplementaryFiles?includeInlineImage=false",
    "PMC7877847-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7877847/fullTextXML",
    "PMC7877847-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7877847/supplementaryFiles?includeInlineImage=false",
    "P22413.fasta": "https://rest.uniprot.org/uniprotkb/P22413.fasta",
    "P22413.txt": "https://rest.uniprot.org/uniprotkb/P22413.txt",
    "Q6UWV6.fasta": "https://rest.uniprot.org/uniprotkb/Q6UWV6.fasta",
    "pubchem-cgamp-135564529.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/135564529/property/Title,IsomericSMILES,InChI,InChIKey/JSON",
    "pubchem-papg-137069451.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/137069451/property/Title,IsomericSMILES,InChI,InChIKey/JSON",
    "ENP-HM102-product.html": "https://www.kactusbio.com/products/human-enpp-1-protein-enp-hm102",
}

STRUCTURES = {
    "2-prime-3-prime-cGAMP": {
        "cid": 135564529,
        "smiles": "C1[C@@H]2[C@H]([C@H]([C@@H](O2)N3C=NC4=C3N=C(NC4=O)N)OP(=O)(OC[C@@H]5[C@H]([C@H]([C@@H](O5)N6C=NC7=C(N=CN=C76)N)O)OP(=O)(O1)O)O)O",
        "inchi_key": "XRILCFTWUCUKJR-INFSMZHSSA-N",
    },
    "pApG": {
        "cid": 137069451,
        "smiles": "C1=NC(=C2C(=N1)N(C=N2)[C@H]3C(C([C@H](O3)COP(=O)([O-])[O-])OP(=O)([O-])OC[C@@H]4C(C([C@@H](O4)N5C=NC6=C5N=C(NC6=O)N)O)O)O)N",
        "inchi_key": "SXJAFMAZBFVAJD-KNKPYYCGSA-K",
    },
}

KINETICS = (
    ("enpp1-cgamp-001", "2-prime-3-prime-cGAMP", 1.2, 0.1, 6.5, 1.2),
    ("enpp1-papg-002", "pApG", 5.2, 0.3, 0.4, 0.3),
)


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
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)
    for archive_name in ("PMC13223937-supplementaryFiles.zip", "PMC7877847-supplementaryFiles.zip"):
        with zipfile.ZipFile(RAW / archive_name) as archive:
            for member in archive.namelist():
                if member.lower().endswith((".docx", ".pdf")):
                    target = RAW / Path(member).name
                    if not target.is_file():
                        with archive.open(member) as source, target.open("wb") as output:
                            shutil.copyfileobj(source, output)


def fasta_sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines() if line and not line.startswith(">"))


def docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    return " ".join(text.strip() for text in root.itertext() if text.strip())


def exact_construct() -> tuple[str, str, str]:
    enpp1 = fasta_sequence(RAW / "P22413.fasta")
    enpp7 = fasta_sequence(RAW / "Q6UWV6.fasta")
    if len(enpp1) != 925 or len(enpp7) != 458:
        raise ValueError("Unexpected UniProt sequence length")
    text = docx_text(RAW / "CTS-14-362-s001.docx")
    compact = re.sub(r"\s+", "", text)
    match = re.search(r"MRGPAVLLTVALATLLAPGAGAPSCAKEVKSCKGRCFERTF.*?NYKTTPPVLD", compact)
    if not match:
        raise ValueError("Construct 770 sequence block is absent")
    # Word stores the final six sequence lines before the main sequence text box.
    fc_tail = "SDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    if not all(peptide in compact for peptide in ("SDGSFFLYSK", "LTVDKSRWQQ", "GNVFSCSVMH", "EALHYTQKSL", "SLSPGK")):
        raise ValueError("Construct 770 Fc tail evidence is incomplete")
    precursor = match.group() + fc_tail
    signal = enpp7[:21]
    component = enpp1[98:]
    if precursor[:21] != signal or precursor[22:849] != component or len(precursor) != 1078:
        raise ValueError("Construct 770 does not map to ENPP7(1-21)/ENPP1(99-925)")
    mature = precursor[21:]
    if hashlib.sha256(mature.encode("ascii")).hexdigest() != "c671dc50b0daa0c1f6511d7603d2cc19969e3d6c04e2b5156e840e14db9a0b0e":
        raise ValueError("Construct 770 mature sequence identity changed")
    return precursor, mature, component


def validate_evidence() -> tuple[str, str, str]:
    xml = (RAW / "PMC13223937-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI, "creativecommons.org/licenses/by/4.0", "construct 770", "Glycosylated human ENPP1",
        "Lys98-Asp925", "1.2", "5.2", "6.5", "0.4", "Briggs-Haldane",
        "final mononucleotide products", "at least one step before pApG hydrolysis",
    )
    missing = [marker for marker in required if marker not in xml]
    if missing:
        raise ValueError(f"Primary article evidence is absent: {missing}")
    product = (RAW / "ENP-HM102-product.html").read_text(encoding="utf-8")
    if not all(marker in product for marker in ("Lys98-Asp925", "HEK293", "C-His")):
        raise ValueError("Commercial comparator boundaries are not corroborated")
    for substrate, raw_name in (
        ("2-prime-3-prime-cGAMP", "pubchem-cgamp-135564529.json"),
        ("pApG", "pubchem-papg-137069451.json"),
    ):
        compound = json.loads((RAW / raw_name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        expected = STRUCTURES[substrate]
        if (compound["CID"], compound["SMILES"], compound["InChIKey"]) != (expected["cid"], expected["smiles"], expected["inchi_key"]):
            raise ValueError(f"PubChem mapping changed for {substrate}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference changed")
    return exact_construct()


def exact_overlap(mature: str, development_corpus: Path) -> dict[str, object]:
    if development_corpus.stat().st_size != DEVELOPMENT_CORPUS_SIZE or sha256(development_corpus) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("Frozen development corpus is absent or changed")
    rows = json.loads(development_corpus.read_text(encoding="utf-8"))
    if len(rows) != 17010:
        raise ValueError("Frozen development row count changed")
    sequence_rows = [row for row in rows if row.get("Sequence") == mature]
    pairs = {
        substrate: sum(row.get("Smiles") == structure["smiles"] for row in sequence_rows)
        for substrate, structure in STRUCTURES.items()
    }
    return {
        "method": "Exact string comparison of the mature assayed fusion sequence and PubChem isomeric SMILES against every frozen development row",
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_corpus_rows": len(rows),
        "exact_sequence_rows": len(sequence_rows),
        "exact_sequence_substrate_pairs": pairs,
        "exact_sequence_substrate_overlap": sum(pairs.values()),
    }


def mmseqs_command(executable: str, *arguments: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *arguments]
    _, distribution, binary = executable.split(":", 2)
    return ["wsl", "-d", distribution, binary, *arguments]


def tool_path(path: Path, executable: str) -> str:
    if not executable.startswith("wsl:"):
        return str(path)
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def run_checked(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(command)}\n{result.stderr}")
    return (result.stdout or result.stderr).strip()


def run_homology(query: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(mmseqs_command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 {version!r} differs from frozen version")
    homology = SOURCE / "homology"
    if homology.exists():
        run_checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(homology, executable)]) if executable.startswith("wsl:") else shutil.rmtree(homology)
    homology.mkdir(parents=True)
    hits = homology / "homology_hits.tsv"
    prefix = homology / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    common = ("--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads))
    run_checked(mmseqs_command(executable, "easy-search", tool_path(query, executable), tool_path(REFERENCE, executable), tool_path(hits, executable), tool_path(homology / "search-tmp", executable), *common, "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"))
    run_checked(mmseqs_command(executable, "easy-linclust", tool_path(query, executable), tool_path(prefix, executable), tool_path(homology / "cluster-tmp", executable), *common))
    return hits, prefix.with_name("proteins_cluster.tsv"), version


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(acquire: bool, mmseqs: str, threads: int, development_corpus: Path) -> None:
    if acquire:
        download_raw()
    precursor, mature, component = validate_evidence()
    overlap = exact_overlap(mature, development_corpus)
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "construct_sequences.fasta"
    sequences = {"enpp1-fc-770-mature": mature, "enpp1-fc-770-catalytic-component": component}
    with query.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            handle.write("\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
    hits, cluster, version = run_homology(query, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="utf-8").splitlines() if line]
    hit_ids = {row[0] for row in hit_rows}
    excluded = bool(hit_ids & set(sequences))
    status = "excluded_homology" if excluded else "accepted_homology_cold_pool"

    rows = []
    for candidate_id, substrate, kcat, kcat_error, km_um, km_error_um in KINETICS:
        structure = STRUCTURES[substrate]
        rows.append({
            "candidate_id": candidate_id, "article_doi": DOI,
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC13223937",
            "source_file": "raw/PMC13223937-fullText.xml", "source_row": f"Table 2; {substrate}",
            "organism": "Homo sapiens", "enzyme_identity": "glycosylated human ENPP1-Fc construct 770",
            "sequence_id": "enpp1-fc-770-mature", "construct": "CHO-produced mature ENPP7 signal-cleaved ENPP1(99-925)-human IgG1 Fc fusion; 1057 aa",
            "variable_substrate": substrate, "substrate_pubchem_cid": structure["cid"],
            "substrate_isomeric_smiles": structure["smiles"], "endpoint": "kcat_s-1",
            "kcat_s-1": kcat, "kcat_standard_error_s-1": kcat_error,
            "km_uM": km_um, "km_standard_error_uM": km_error_um,
            "assay_temperature_C": 25, "assay_pH": 7.4,
            "assay_buffer": "20 mM Tris-HCl, 154 mM NaCl, 0.014 mM ZnCl2, 1 mM MgCl2, 1 mM CaCl2, 4.5 mM KCl",
            "fit": "initial GMP liberation rate per active site versus substrate fit directly to Briggs-Haldane rectangular hyperbola",
            "endpoint_semantics": "overall cGAMP-to-AMP+GMP two-cleavage cycle" if substrate.startswith("2-prime") else "pApG-to-AMP+GMP cleavage",
            "substrate_series_uM": "2-50", "saturation_evidence": f"direct hyperbolic fit; maximum observed concentration is {50 / km_um:.1f}x fitted KM",
            "status_at_normalization": status,
        })
    write_csv(SOURCE / "excluded_records.csv" if excluded else SOURCE / "candidate_records.csv", rows)
    empty_path = SOURCE / ("candidate_records.csv" if excluded else "excluded_records.csv")
    empty_path.write_text("", encoding="ascii")

    raw_hashes = {path.name: {"url": URLS.get(path.name), "size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    sequence_rows = [
        {"sequence_id": "enpp1-fc-770-precursor", "length_aa": len(precursor), "role": "expressed precursor", "mapping": "ENPP7 1-21 signal + A junction + ENPP1 P22413 99-925 + RS junction + human IgG1 Fc", "sha256": hashlib.sha256(precursor.encode()).hexdigest()},
        {"sequence_id": "enpp1-fc-770-mature", "length_aa": len(mature), "role": "exact mature glycoprotein polypeptide", "mapping": "signal-cleaved A junction + ENPP1 P22413 99-925 + RS junction + human IgG1 Fc", "sha256": hashlib.sha256(mature.encode()).hexdigest()},
        {"sequence_id": "enpp1-fc-770-catalytic-component", "length_aa": len(component), "role": "homology component", "mapping": "P22413 residues 99-925", "sha256": hashlib.sha256(component.encode()).hexdigest()},
    ]
    write_csv(SOURCE / "sequence_mapping.csv", sequence_rows)
    accepted = 0 if excluded else len(rows)
    provenance = {
        "source_id": SOURCE_ID, "article_doi": DOI, "article_published": "2026-04-27", "license": "CC-BY-4.0",
        "biostudies_status": "Stable accession page resolves; API record returned 404 during acquisition. Europe PMC's BioStudies-backed mirror supplied the hashed article supplement package.",
        "direct_table_2_kcat_rows": 2, "accepted_records": accepted,
        "selection_policy": "All direct finite Table 2 Briggs-Haldane kcat rows; no graph digitization or value-based selection.",
        "construct_resolution": "Construct 770 Supplementary Figure 1 gives a 1078-aa precursor: ENPP7 signal 1-21, A junction at 22, ENPP1 P22413 residues 99-925 at construct 23-849, RS junction at 850-851, and IgG1 Fc at 852-1078. Signal cleavage yields the exact 1057-aa assayed polypeptide. CHO production and source glycan analyses establish glycosylation; glycans are PTMs and are not encoded as amino acids.",
        "commercial_comparator_boundary": "ENP-HM102 is a separate HEK293-produced P22413 Lys98-Asp925 C-His glycoprotein used only in Figure 2C validation; it yielded an identical time course but was not used for the Table 2 fits.",
        "endpoint_audit": "At 260 nm cGAMP and pApG have nearly equal extinction coefficients, so the real-time signal excludes pApG accumulation and measures final AMP/GMP liberation. The cGAMP kcat is therefore the complete two-step cycle, not ring opening. The 4-fold faster pApG kcat shows at least one pre-pApG step limits overall cGAMP turnover.",
        "saturation_audit": "Both rows are direct substrate-series hyperbolic fits over 2-50 uM. Max/KM is 7.69 for cGAMP and 125 for pApG; fitted kcat is the explicitly reported saturating-substrate asymptote.",
        "structure_mapping": STRUCTURES, "exact_overlap_audit": overlap,
        "final_disposition": "excluded_by_frozen_homology_gate" if excluded else status,
        "model_predictions_run": False, "raw_file_hashes": "raw-file-hashes.json",
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")
    blocker = {
        "schema_version": 1, "source_id": SOURCE_ID, "article_doi": DOI,
        "status": "excluded_homology", "accepted_records": accepted,
        "blocker_code": "FROZEN_DEVELOPMENT_HOMOLOGY",
        "evidence": [
            "Pinned MMseqs2 found a qualifying hit for the exact 827-aa ENPP1 catalytic component at 80.0% identity, 100% query coverage, and 91.2% target coverage.",
            "The full 1057-aa mature ENPP1-Fc fusion has no threshold hit because the Fc fusion dilutes query coverage; the frozen fail-closed component rule prevents this artificial escape.",
            "The source establishes a glycosylated CHO product and an exact 1057-aa mature polypeptide, but does not define one homogeneous glycan composition. Glycan microheterogeneity is recorded as a PTM limitation and is not represented as invented amino-acid sequence.",
        ],
        "otherwise_valid_rows": 2,
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"],
        "resolution_required": "No source clarification can reverse the frozen development-homology exclusion.",
        "model_predictions_run": False,
    }
    (SOURCE / "blocker-evidence.json").write_text(json.dumps(blocker, indent=2) + "\n", encoding="ascii")
    cluster_rows = [line.split("\t") for line in cluster.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust of exact mature fusion plus catalytic ENPP1 component",
        "fusion_homology_rule": "Exclude when either the exact mature chain or its mapped catalytic component has a qualifying development hit; prevents Fc coverage dilution.",
        "mmseqs_version": version, "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "query_sequences": len(sequences),
        "construct_sequences_sha256": sha256(query), "homology_hit_query_sequences": len(hit_ids),
        "homology_hit_alignments": len(hit_rows), "homology_hits_sha256": sha256(hits),
        "candidate_mmseqs_families": len({row[0] for row in cluster_rows}), "family_cluster_sha256": sha256(cluster),
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"],
        "candidate_records": len(rows), "accepted_records": accepted, "status": status,
        "readiness_gate_passes": False, "claim_boundary": "Curation/exclusion evidence only; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-acquire", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--development-corpus", type=Path, default=DEFAULT_DEVELOPMENT_CORPUS)
    args = parser.parse_args()
    curate(not args.no_acquire, args.mmseqs, args.threads, args.development_corpus)
