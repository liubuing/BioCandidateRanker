from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import date
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC13151041"
DOI = "10.1021/acs.biochem.6c00183"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
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
    "PMC13151041-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13151041/fullTextXML",
    "PMC13151041-article.pdf": "https://europepmc.org/articles/PMC13151041?pdf=render",
    "PMC13151041-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13151041/supplementaryFiles?includeInlineImage=false",
    "B3PHM6.fasta": "https://rest.uniprot.org/uniprotkb/B3PHM6.fasta",
    "B3PHM6.json": "https://rest.uniprot.org/uniprotkb/B3PHM6.json",
    "pubchem-ABTS-35688.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/35688/property/IsomericSMILES,Title/JSON",
    "pubchem-DMP-7041.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/7041/property/IsomericSMILES,Title/JSON",
    "pubchem-Fe2-27284.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/27284/property/IsomericSMILES,Title/JSON",
}

SUBSTRATES = {
    "ABTS": (35688, "CCN1C2=C(C=C(C=C2)S(=O)(=O)O)SC1=NN=C3N(C4=C(S3)C=C(C=C4)S(=O)(=O)O)CC"),
    "Fe(II)": (27284, "[Fe+2]"),
    "2,6-DMP": (7041, "COC1=C(C(=CC=C1)OC)O"),
}

# Table S2, in source order: pH, buffer, kcat, kcat SD, Km, Km SD, kcat/Km, kcat/Km SD.
DMP_ROWS = (
    (4.0, "50 mM sodium acetate", 0.074, 0.004, 2.4, 0.3, 30, 4),
    (4.5, "50 mM sodium acetate", 0.18, 0.02, 2.1, 0.4, 90, 20),
    (5.0, "50 mM sodium acetate", 0.30, 0.02, 2.9, 0.4, 100, 20),
    (5.5, "50 mM sodium acetate", 0.214, 0.006, 2.4, 0.1, 89, 6),
    (5.5, "50 mM MES", 0.12, 0.02, 6.0, 1.0, 18, 5),
    (6.0, "50 mM MES", 0.07, 0.01, 3.0, 1.0, 22, 8),
    (6.5, "50 mM MES", 0.076, 0.008, 1.9, 0.5, 40, 10),
    (6.5, "50 mM MOPS", 0.074, 0.004, 3.7, 0.4, 20, 2),
    (7.0, "50 mM MOPS", 0.051, 0.007, 1.5, 0.5, 30, 10),
    (7.5, "50 mM MOPS", 0.065, 0.007, 1.8, 0.5, 40, 10),
    (7.5, "50 mM TRIS", 0.26, 0.02, 1.9, 0.3, 140, 20),
    (8.0, "50 mM TRIS", 0.21, 0.01, 1.4, 0.2, 160, 20),
    (8.5, "50 mM TRIS", 0.137, 0.007, 1.1, 0.2, 130, 20),
    (9.0, "50 mM TRIS", 0.08, 0.01, 0.7, 0.3, 110, 40),
    (5.0, "50 mM Citrate", 0.100, 0.008, 2.5, 0.3, 40, 6),
    (5.0, "50 mM Citrate + 100 mM KCl", 0.032, 0.001, 0.89, 0.09, 35, 4),
    (5.0, "50 mM Citrate + 500 mM KCl", 0.008, 0.000, 0.750, 0.000, 10.76, 0.00),
    (5.0, "500 mM Citrate", 0.09, 0.01, 3.4, 0.8, 25, 7),
    (6.0, "50 mM Citrate", 0.100, 0.008, 3.3, 0.5, 30, 5),
    (6.0, "50 mM Phosphate", 0.165, 0.005, 1.28, 0.09, 130, 10),
    (7.0, "50 mM Phosphate", 0.101, 0.001, 0.91, 0.03, 111, 4),
    (7.0, "50 mM Phosphate + 100 mM KCl", 0.116, 0.002, 0.66, 0.02, 176, 6),
    (7.0, "50 mM Phosphate + 500 mM KCl", 0.095, 0.001, 0.29, 0.01, 330, 10),
    (7.0, "500 mM Phosphate", 0.143, 0.004, 0.92, 0.04, 155, 8),
    (8.0, "50 mM Phosphate", 0.111, 0.002, 0.57, 0.03, 196, 9),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            (RAW / name).write_bytes(response.read())
    with zipfile.ZipFile(RAW / "PMC13151041-supplementaryFiles.zip") as archive:
        member = next(name for name in archive.namelist() if name.endswith("_si_001.pdf"))
        (RAW / "bi6c00183_si_001.pdf").write_bytes(archive.read(member))

    # Preserve negative exact-structure lookups as machine-readable HTTP evidence.
    lookups = {}
    for database, url in {
        "rcsb_uniprot": "https://data.rcsb.org/rest/v1/core/uniprot/B3PHM6",
        "pdbe_uniprot_mapping": "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/B3PHM6",
    }.items():
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                lookups[database] = {"url": url, "http_status": response.status, "body": response.read().decode("utf-8")}
        except urllib.error.HTTPError as exc:
            lookups[database] = {"url": url, "http_status": exc.code, "body": exc.read().decode("utf-8", errors="replace")}
    (RAW / "experimental-structure-lookups.json").write_text(
        json.dumps(lookups, indent=2) + "\n", encoding="ascii"
    )


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines() if not line.startswith(">"))


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
    completed = subprocess.run(args, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"Command failed: {' '.join(args)}\n{completed.stderr}")
    return (completed.stdout or completed.stderr).strip()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def verify_raw() -> tuple[str, str]:
    missing = [name for name in (*URLS, "bi6c00183_si_001.pdf", "experimental-structure-lookups.json") if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {missing}; run with --download")
    root = ET.parse(RAW / "PMC13151041-fullText.xml").getroot()
    if root.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Article DOI mismatch")
    sequence = read_fasta(RAW / "B3PHM6.fasta")
    if len(sequence) != 468:
        raise ValueError("Unexpected B3PHM6 length")
    mature = sequence[24:468]
    if len(mature) != 444:
        raise ValueError("Mature B3PHM6(25-468) boundary failed")
    uniprot = json.loads((RAW / "B3PHM6.json").read_text(encoding="utf-8"))
    pdb_refs = [item for item in uniprot["uniProtKBCrossReferences"] if item["database"] == "PDB"]
    if pdb_refs:
        raise ValueError(f"Unexpected PDB cross-reference requires manual exact-structure audit: {pdb_refs}")
    si_text = "\n".join(page.get_text() for page in fitz.open(RAW / "bi6c00183_si_001.pdf"))
    for token in ("Table S2", "0.074 ± 0.004", "0.111 ± 0.002", "0.57 ± 0.03"):
        if token not in si_text:
            raise ValueError(f"Required Table S2 evidence missing: {token}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference changed")
    for name, (expected_cid, expected_smiles) in SUBSTRATES.items():
        raw_name = {"ABTS": "pubchem-ABTS-35688.json", "Fe(II)": "pubchem-Fe2-27284.json", "2,6-DMP": "pubchem-DMP-7041.json"}[name]
        item = json.loads((RAW / raw_name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if item["CID"] != expected_cid or item.get("SMILES", item.get("IsomericSMILES")) != expected_smiles:
            raise ValueError(f"Unexpected PubChem structure for {name}")
    return sequence, mature


def exact_overlap(mature: str, development_corpus: Path) -> dict[str, object]:
    if not development_corpus.is_file() or development_corpus.stat().st_size != DEVELOPMENT_CORPUS_SIZE:
        raise FileNotFoundError("Frozen 17,010-row development corpus is absent or has changed size")
    if sha256(development_corpus) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("Frozen development corpus hash changed")
    rows = json.loads(development_corpus.read_text(encoding="utf-8"))
    if not isinstance(rows, list) or len(rows) != 17010:
        raise ValueError("Frozen development corpus row count changed")
    exact_sequence_rows = [row for row in rows if isinstance(row, dict) and row.get("Sequence") == mature]
    substrate_smiles = {item[1] for item in SUBSTRATES.values()}
    exact_pairs = [row for row in exact_sequence_rows if row.get("Smiles") in substrate_smiles]
    return {
        "method": "Exact string comparison of mature B3PHM6(25-468) and three PubChem isomeric SMILES against all frozen development rows",
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_corpus_size_bytes": DEVELOPMENT_CORPUS_SIZE,
        "development_corpus_rows": len(rows),
        "exact_sequence_rows": len(exact_sequence_rows),
        "exact_sequence_substrate_overlap": len(exact_pairs),
    }


def run_homology(fasta: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise RuntimeError(f"MMseqs2 version {version!r} != frozen version {MMSEQS_VERSION!r}")
    workspace = SOURCE / "homology"
    if workspace.exists():
        if executable.startswith("wsl:"):
            run_checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(workspace, executable)])
        else:
            shutil.rmtree(workspace)
    workspace.mkdir(parents=True)
    hits = workspace / "homology_hits.tsv"
    prefix = workspace / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    run_checked(command(executable, "easy-search", tool_path(fasta, executable), tool_path(REFERENCE, executable), tool_path(hits, executable), tool_path(workspace / "search-tmp", executable), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits", "--threads", str(threads)))
    run_checked(command(executable, "easy-linclust", tool_path(fasta, executable), tool_path(prefix, executable), tool_path(workspace / "cluster-tmp", executable), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads)))
    return hits, prefix.with_name("proteins_cluster.tsv"), version


def excluded_records() -> list[dict[str, object]]:
    reason = "oxygen_concentration_and_saturation_not_established_for_obligate_O2_cosubstrate"
    common = {
        "article_doi": DOI,
        "organism": "Cellvibrio japonicus Ueda107",
        "sequence_id": "B3PHM6-mature-25-468",
        "uniprot_accession": "B3PHM6",
        "construct": "untagged mature B3PHM6 residues 25-468",
        "endpoint": "reported_kcat_s-1",
        "assay_temperature_C": 25,
        "oxygen_saturation_status": "not_demonstrated",
        "status_at_normalization": "excluded_cosubstrate_saturation_unresolved",
        "exclusion_reason": reason,
        "candidate_label_created": False,
    }
    rows = [
        {**common, "evidence_id": "cjmco-abts-main", "source_table": "Main text and Figure 3A", "source_row": "ABTS substrate-inhibition fit", "substrate": "ABTS", "substrate_pubchem_cid": 35688, "substrate_isomeric_smiles": SUBSTRATES["ABTS"][1], "reported_kcat_s-1": 0.35, "reported_kcat_sd_s-1": "", "assay_pH": 5.5, "assay_buffer": "100 mM sodium acetate", "donor_substrate_fit": "substrate inhibition"},
        {**common, "evidence_id": "cjmco-feii-main", "source_table": "Main text and Figure 3B", "source_row": "Fe(II) Michaelis-Menten fit", "substrate": "Fe(II)", "substrate_pubchem_cid": 27284, "substrate_isomeric_smiles": SUBSTRATES["Fe(II)"][1], "reported_kcat_s-1": 0.43, "reported_kcat_sd_s-1": "", "assay_pH": 6.5, "assay_buffer": "100 mM MES", "donor_substrate_fit": "Michaelis-Menten"},
    ]
    for index, (ph, buffer, kcat, kcat_sd, km, km_sd, efficiency, efficiency_sd) in enumerate(DMP_ROWS, 1):
        rows.append({**common, "evidence_id": f"cjmco-dmp-s2-{index:02d}", "source_table": "Supporting Information Table S2", "source_row": f"Table S2 row {index}: pH {ph}; {buffer}", "substrate": "2,6-DMP", "substrate_pubchem_cid": 7041, "substrate_isomeric_smiles": SUBSTRATES["2,6-DMP"][1], "reported_kcat_s-1": kcat, "reported_kcat_sd_s-1": kcat_sd, "assay_pH": ph, "assay_buffer": buffer, "donor_substrate_fit": "Michaelis-Menten", "reported_km_mM": km, "reported_km_sd_mM": km_sd, "reported_kcat_per_km_M-1_s-1": efficiency, "reported_kcat_per_km_sd_M-1_s-1": efficiency_sd})
    return rows


def write_outputs(mmseqs: str, threads: int, development_corpus: Path) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    sequence, mature = verify_raw()
    overlap = exact_overlap(mature, development_corpus)
    fasta = SOURCE / "construct_sequences.fasta"
    fasta.write_text(">B3PHM6-mature-25-468 | exact untagged kinetic construct\n" + "\n".join(mature[i:i + 80] for i in range(0, len(mature), 80)) + "\n", encoding="ascii")
    hits, cluster, version = run_homology(fasta, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="utf-8").splitlines() if line]
    rows = excluded_records()
    write_csv(SOURCE / "excluded_records.csv", rows)
    (SOURCE / "candidate_records.csv").write_text("", encoding="ascii")
    write_csv(SOURCE / "sequence_mapping.csv", [{"sequence_id": "B3PHM6-mature-25-468", "uniprot_accession": "B3PHM6", "native_length": len(sequence), "assayed_product_length": len(mature), "native_mapping": "exact B3PHM6 residues 25-468 inclusive; signal peptide residues 1-24 removed", "expression_tag": "none", "sequence_sha256": hashlib.sha256(mature.encode("ascii")).hexdigest()}])
    write_csv(SOURCE / "structure_mapping.csv", [{"uniprot_accession": "B3PHM6", "kinetic_construct": "mature residues 25-468", "exact_experimental_structure_status": "none_identified", "uniprot_pdb_cross_references": 0, "rcsb_uniprot_http_status": 404, "pdbe_uniprot_mapping_http_status": 404, "predicted_structure_status": "excluded_not_acquired", "boundary": "Article AlphaFold 3 trimer is a model prediction and is not admissible as an exact experimental structure."}])
    derived_sources = {
        "bi6c00183_si_001.pdf": "extracted byte-for-byte from raw/PMC13151041-supplementaryFiles.zip",
        "experimental-structure-lookups.json": "RCSB and PDBe API responses generated by this script",
    }
    raw_hashes = {path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path), "url_or_derivation": URLS[path.name] if path.name in URLS else derived_sources[path.name]} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    blocker = {
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "status": "blocked_fail_closed",
        "blocker_code": "OBLIGATE_O2_COSUBSTRATE_SATURATION_UNRESOLVED",
        "candidate_records_released": False,
        "reported_direct_kcat_evidence_rows_retained": len(rows),
        "affected_substrates": ["ABTS", "Fe(II)", "2,6-DMP"],
        "decision": "Reported values are apparent donor-substrate turnover parameters, not defensible oxygen-saturated direct kcat labels under the frozen protocol.",
        "evidence": [
            "Kinetic reactions used 1 uM enzyme, 0-5 mM donor substrate, 25 C, and initial rates through 30 s; the method reports no dissolved O2 concentration, gas equilibration, vessel headspace protocol, or O2 saturation series.",
            "The separate stopped-flow oxidation experiment used measured 500 uM O2, but this condition is not stated for the ABTS, Fe(II), or DMP steady-state assays.",
            "The article attributes the approximately 0.4 s-1 turnover to slow enzyme reoxidation by O2, making unverified O2 independence mechanistically material.",
        ],
        "resolution_required": "Source evidence showing steady-state kcat is independent of O2 concentration under each assay geometry/condition, or donor-substrate fits repeated at demonstrably saturating measured O2.",
        "model_predictions_run": False,
        "recorded_on": str(date.today()),
    }
    (SOURCE / "blocker-evidence.json").write_text(json.dumps(blocker, indent=2) + "\n", encoding="ascii")
    structure = json.loads((RAW / "experimental-structure-lookups.json").read_text(encoding="ascii"))
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC13151041",
        "article_doi": DOI,
        "pmc_id": "PMC13151041",
        "article_published": "2026-04-27",
        "article_and_supporting_information_license": "CC-BY-4.0",
        "kinetics_sources": ["raw/PMC13151041-fullText.xml, Methods 2.11, Results and Figure 3", "raw/bi6c00183_si_001.pdf, Table S2"],
        "sequence_source": "raw/B3PHM6.fasta and raw/B3PHM6.json; article explicitly defines mature untagged residues 25-468",
        "substrate_structure_source": "raw PubChem PUG REST snapshots for CIDs 35688, 27284, and 7041",
        "exact_structure_audit": {"status": "no_exact_experimental_structure_identified", "uniprot_pdb_cross_references": 0, "api_lookups": structure, "prediction_boundary": "AlphaFoldDB and article AlphaFold 3 models excluded; no prediction downloaded or used."},
        "exact_overlap_audit": overlap,
        "oxygen_saturation_audit": blocker,
        "reported_rows": len(rows),
        "accepted_records": 0,
        "selection_policy": "Audit all reported direct ABTS and Fe(II) kcat values and every Table S2 DMP kcat row without value-based selection; fail closed on the shared unresolved oxygen-cosubstrate saturation condition.",
        "raw_file_hashes": "raw-file-hashes.json",
        "artifact_sha256": {name: sha256(SOURCE / name) for name in ("excluded_records.csv", "construct_sequences.fasta", "sequence_mapping.csv", "structure_mapping.csv", "blocker-evidence.json")},
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")
    cluster_rows = [line.split("\t") for line in cluster.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust", "mmseqs_version": version,
        "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "construct_sequences_sha256": sha256(fasta),
        "homology_hit_sequences": len({row[0] for row in hit_rows}), "homology_hits": len(hit_rows),
        "exact_sequence_overlap": sum(float(row[2]) >= 0.999 for row in hit_rows),
        "homology_hits_sha256": sha256(hits), "candidate_mmseqs_families": len({row[0] for row in cluster_rows}),
        "family_cluster_sha256": sha256(cluster), "accepted_records": 0,
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"],
        "exact_overlap_development_corpus_sha256": overlap["development_corpus_sha256"],
        "homology_disposition": "Reported kinetics remain oxygen-blocked regardless of homology result.",
        "model_predictions_run": False,
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--development-corpus", type=Path, default=DEFAULT_DEVELOPMENT_CORPUS)
    args = parser.parse_args()
    if args.download:
        download()
    write_outputs(args.mmseqs, args.threads, args.development_corpus)


if __name__ == "__main__":
    main()
