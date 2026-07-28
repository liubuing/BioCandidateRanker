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


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC12781114"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1021/acs.biochem.5c00559"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
DEVELOPMENT_CORPUS = Path(
    r"C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment"
    r"\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json"
)
DEVELOPMENT_CORPUS_SHA256 = "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719"
DEVELOPMENT_CORPUS_SIZE = 12132719
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

URLS = {
    "biostudies-metadata.json": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12781114",
    "PMC12781114-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12781114/fullTextXML",
    "PMC12781114-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12781114/supplementaryFiles",
    "Q9I4V0.fasta": "https://rest.uniprot.org/uniprotkb/Q9I4V0.fasta",
    "Q9I4V0.json": "https://rest.uniprot.org/uniprotkb/Q9I4V0.json",
    "6E2A-polymer-entity.json": "https://data.rcsb.org/rest/v1/core/polymer_entity/6E2A/1",
    "pubchem-NADH.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/439153/property/IsomericSMILES,Title/JSON",
    "pubchem-CoQ0.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/69068/property/IsomericSMILES,Title/JSON",
}
STRUCTURES = {
    "NADH": (439153, "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"),
    "coenzyme Q0": (69068, "CC1=CC(=O)C(=C(C1=O)OC)OC"),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)
    supplement = RAW / "bi5c00559_si_001.pdf"
    if not supplement.is_file():
        with zipfile.ZipFile(RAW / "PMC12781114-supplementaryFiles.zip") as archive:
            member = next(item for item in archive.namelist() if item.endswith("bi5c00559_si_001.pdf"))
            with archive.open(member) as source, supplement.open("wb") as output:
                shutil.copyfileobj(source, output)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def validate_raw() -> tuple[dict[str, str], dict[str, tuple[int, str]]]:
    metadata = json.loads((RAW / "biostudies-metadata.json").read_text(encoding="utf-8"))
    if metadata.get("accno") != "S-EPMC12781114":
        raise ValueError("Unexpected BioStudies accession")
    files = metadata["section"]["files"]
    if not any(item.get("path") == "bi5c00559_si_001.pdf" for item in files):
        raise ValueError("BioStudies metadata lacks the expected supporting PDF")

    xml = (RAW / "PMC12781114-fullText.xml").read_text(encoding="utf-8")
    markers = (
        DOI, "creativecommons.org/licenses/by/4.0/", "pET20b", "NQO-P78G",
        "11 ± 1", "5.4", "130 ±", "35 ± 6", "10 ± 1", "30–100",
        "30–150", "10 to 100", "2 to 10", "saturated concentrations for both substrates",
    )
    missing = [marker for marker in markers if marker not in xml]
    if missing:
        raise ValueError(f"Primary XML lacks expected construct/kinetic evidence: {missing}")
    if (RAW / "bi5c00559_si_001.pdf").read_bytes()[:5] != b"%PDF-":
        raise ValueError("BioStudies supporting information is not a PDF")

    native = read_fasta(RAW / "Q9I4V0.fasta")
    uniprot = json.loads((RAW / "Q9I4V0.json").read_text(encoding="utf-8"))
    if uniprot["primaryAccession"] != "Q9I4V0" or uniprot["sequence"]["value"] != native:
        raise ValueError("UniProt JSON/FASTA sequence mismatch")
    if len(native) != 328 or native[77] != "P":
        raise ValueError("Q9I4V0 is not the expected 328-aa P78 parent")

    pdb = json.loads((RAW / "6E2A-polymer-entity.json").read_text(encoding="utf-8"))
    pdb_product = pdb["entity_poly"]["pdbx_seq_one_letter_code_can"]
    wt = native + "HHHHHH"
    if pdb_product != wt:
        raise ValueError("6E2A does not resolve Q9I4V0 followed directly by His6")
    p78g = list(wt)
    p78g[77] = "G"
    constructs = {"nqo-wt-his6": wt, "nqo-p78g-his6": "".join(p78g)}

    structures = {}
    for name, raw_name in (("NADH", "pubchem-NADH.json"), ("coenzyme Q0", "pubchem-CoQ0.json")):
        item = json.loads((RAW / raw_name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        expected_cid, expected_smiles = STRUCTURES[name]
        if item["CID"] != expected_cid or item["SMILES"] != expected_smiles:
            raise ValueError(f"Unexpected PubChem mapping for {name}")
        structures[name] = (item["CID"], item["SMILES"])
    return constructs, structures


def exact_overlap(constructs: dict[str, str], structures: dict[str, tuple[int, str]]) -> dict[str, object]:
    if not DEVELOPMENT_CORPUS.is_file() or DEVELOPMENT_CORPUS.stat().st_size != DEVELOPMENT_CORPUS_SIZE:
        raise FileNotFoundError("Frozen 17,010-row development corpus is absent or changed size")
    if sha256(DEVELOPMENT_CORPUS) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("Frozen development corpus SHA256 changed")
    rows = json.loads(DEVELOPMENT_CORPUS.read_text(encoding="utf-8"))
    if not isinstance(rows, list) or len(rows) != 17010:
        raise ValueError("Frozen development corpus row count changed")
    candidate_sequences = set(constructs.values())
    substrate_smiles = {smiles for _, smiles in structures.values()}
    sequence_rows = [row for row in rows if isinstance(row, dict) and row.get("Sequence") in candidate_sequences]
    pairs = [row for row in sequence_rows if row.get("Smiles") in substrate_smiles]
    return {
        "method": "Exact string comparison of both complete His6 assay products and both PubChem isomeric SMILES against all frozen development rows",
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_corpus_size_bytes": DEVELOPMENT_CORPUS_SIZE,
        "development_corpus_rows": len(rows),
        "exact_sequence_rows": len(sequence_rows),
        "exact_sequence_substrate_overlap": len(pairs),
    }


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive.rstrip(':').lower()}/{resolved.as_posix()[3:]}"


def mmseqs_command(*arguments: str) -> list[str]:
    _, distribution, executable = MMSEQS.split(":", 2)
    return ["wsl", "-d", distribution, executable, *arguments]


def reset_wsl(path: Path) -> None:
    if path.exists():
        subprocess.run(["wsl", "-d", MMSEQS.split(":", 2)[1], "rm", "-rf", "--", wsl_path(path)], check=True)


def run_mmseqs(query: Path) -> str:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development FASTA is absent or changed")
    version = subprocess.run(mmseqs_command("version"), check=True, capture_output=True, text=True).stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 version {version!r} is not frozen version {MMSEQS_VERSION!r}")
    reset_wsl(HOMOLOGY)
    HOMOLOGY.mkdir(parents=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    subprocess.run(mmseqs_command(
        "easy-search", wsl_path(query), wsl_path(REFERENCE), wsl_path(hits), wsl_path(HOMOLOGY / "search-tmp"),
        "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0",
        "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
    ), check=True)
    prefix = HOMOLOGY / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    subprocess.run(mmseqs_command(
        "easy-linclust", wsl_path(query), wsl_path(prefix), wsl_path(HOMOLOGY / "cluster-tmp"),
        "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0",
    ), check=True)
    return version


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(*, acquire: bool, run_homology: bool) -> None:
    if acquire:
        download_raw()
    constructs, structures = validate_raw()
    overlap = exact_overlap(constructs, structures)
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "construct_sequences.fasta"
    query.write_text("".join(
        f">{sequence_id} exact assayed pET20b(+) product\n"
        + "\n".join(sequence[start:start + 80] for start in range(0, len(sequence), 80)) + "\n"
        for sequence_id, sequence in constructs.items()
    ), encoding="ascii", newline="\n")
    version = run_mmseqs(query) if run_homology else MMSEQS_VERSION
    hit_path = HOMOLOGY / "homology_hits.tsv"
    cluster_path = HOMOLOGY / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        raise RuntimeError("Pinned MMseqs evidence is required before label disposition")
    hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    hit_ids = {line.split("\t", 1)[0] for line in hits}

    nadh_cid, nadh_smiles = structures["NADH"]
    coq_cid, coq_smiles = structures["coenzyme Q0"]
    specifications = (
        ("nqo-001", "NQO-WT", "nqo-wt-his6", "WT", 11.0, 1.0, 0.4, 0.3, 35.0, 6.0, "30-100", "30-150"),
        ("nqo-002", "NQO-P78G", "nqo-p78g-his6", "P78G", 5.4, 0.2, 130.0, 10.0, 10.0, 1.0, "10-100", "2-10"),
    )
    rows = []
    for candidate_id, enzyme, sequence_id, variant, kcat, error, knadh, knadh_error, kcoq, kcoq_error, nadh_range, coq_range in specifications:
        status = "excluded_homology" if sequence_id in hit_ids else "accepted_homology_cold_pool"
        rows.append({
            "candidate_id": candidate_id, "article_doi": DOI,
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12781114",
            "source_file": "raw/PMC12781114-fullText.xml", "source_table": "Table 1",
            "source_row": enzyme, "organism": "Pseudomonas aeruginosa PAO1",
            "enzyme_identity": "NADH:quinone reductase (EC 1.6.5.9)", "sequence_accession": "Q9I4V0",
            "sequence_id": sequence_id, "variant": variant,
            "construct": "Q9I4V0 residues 1-328 + direct C-terminal His6; no pelB signal/linker; P78G at native residue 78 where applicable",
            "reaction": "NADH + H+ + coenzyme Q0 -> NAD+ + reduced coenzyme Q0",
            "reaction_substrate_1": "NADH", "substrate_1_pubchem_cid": nadh_cid,
            "substrate_1_isomeric_smiles": nadh_smiles,
            "reaction_substrate_2": "coenzyme Q0", "substrate_2_pubchem_cid": coq_cid,
            "substrate_2_isomeric_smiles": coq_smiles,
            "endpoint": "kcat_s-1", "kcat_s-1": kcat, "kcat_error_s-1": error,
            "kcat_error_type": "reported plus/minus; statistical type not specified",
            "km_nadh_uM": knadh, "km_nadh_error_uM": knadh_error,
            "km_coq0_uM": kcoq, "km_coq0_error_uM": kcoq_error,
            "nadh_range_uM": nadh_range, "coq0_range_uM": coq_range,
            "fit_model": "global ping-pong bi-bi; WT includes NADH substrate inhibition",
            "saturation_semantics": "global-fit kcat at saturated concentrations of both substrates",
            "assay_pH": 6.0, "assay_temperature_C": 25, "enzyme_concentration_nM": 100,
            "status_at_normalization": status,
        })
    write_csv(SOURCE / "candidate_records.csv", rows)

    exclusions = [
        {"source_row": "Table 1; NQO-WT K_NADH and kcat/K_NADH", "reported_endpoint": "Km and kcat/Km", "value": "0.4 +/- 0.3; 2.8e7 +/- 2.1e7", "unit": "uM; M-1 s-1", "exclusion_reason": "authors explicitly call values estimated and highly incorrect due to NADH substrate inhibition", "candidate_label_created": False},
        {"source_row": "Reductive half-reaction; NQO-WT and P78G", "reported_endpoint": "kred and Kd", "value": "12.9 +/- 0.3; 4.8 +/- 0.3; 450 +/- 40", "unit": "s-1; uM", "exclusion_reason": "pre-steady-state half-reaction endpoints are not direct steady-state kcat/Km labels", "candidate_label_created": False},
        {"source_row": "Molecular-dynamics analyses", "reported_endpoint": "simulated conformational metrics", "value": "", "unit": "", "exclusion_reason": "predictions/simulations are prohibited as labels", "candidate_label_created": False},
    ]
    write_csv(SOURCE / "exclusions.csv", exclusions)

    raw_hashes = {path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    accepted = [row for row in rows if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12781114",
        "article_doi": DOI, "article_published": "2025-12-13", "license": "CC-BY-4.0",
        "license_scope": "Primary article, embedded Table 1, and supporting information",
        "reported_direct_finite_global_fit_kcat_rows": 2, "curated_reaction_level_kcat_rows": 2,
        "accepted_records": len(accepted),
        "kinetics_source": "raw/PMC12781114-fullText.xml, Experimental Procedures 2.7/2.9 and Table 1",
        "construct_resolution": {
            "vector": "pET20b(+)", "pelB_signal_retained": False,
            "wt_formula": "Q9I4V0(1-328) + HHHHHH", "p78g_formula": "Q9I4V0 P78G + HHHHHH",
            "length_aa_each": 334,
            "evidence": "Current article names pET20b(+)/NQO and sequence-confirmed P78G; same-protein experimental PDB 6E2A polymer entity independently records Q9I4V0(1-328) followed directly by residues 329-334 His6 and marks only those six residues as expression tag.",
            "boundary": "6E2A is construct evidence, not evidence that the crystallized aliquot was used in this assay; P78G was purified by the cited WT method.",
        },
        "reaction_level_label_semantics": "Table 1 gives one global two-substrate kcat per construct. It is not duplicated as separate NADH-varied and CoQ0-varied labels.",
        "structure_mapping": {name: {"pubchem_cid": cid, "isomeric_smiles": smiles} for name, (cid, smiles) in structures.items()},
        "exact_overlap_audit": overlap, "raw_file_hashes": "raw-file-hashes.json",
        "homology_status": "excluded_homology" if hits else "accepted_homology_cold_pool",
        "final_disposition": "excluded_by_frozen_homology_gate" if hits else "accepted_homology_cold_pool",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    families = {line.split("\t", 1)[0] for line in cluster_path.read_text(encoding="utf-8").splitlines() if line}
    audit = {
        "audited_on": "2026-07-27", "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust", "mmseqs_version": version,
        "wsl_distribution": "Ubuntu-24.04", "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "search_command": "mmseqs easy-search construct_sequences.fasta unique_proteins.fasta homology_hits.tsv search-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "cluster_command": "mmseqs easy-linclust construct_sequences.fasta proteins cluster-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0",
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": {"count": len(rows), "sha256": sha256(SOURCE / "candidate_records.csv")},
        "unique_sequences": len(constructs), "construct_sequences_sha256": sha256(query),
        "homology_hit_sequences": len(hit_ids), "homology_hit_alignments": len(hits),
        "homology_hits_sha256": sha256(hit_path), "candidate_mmseqs_families": len(families),
        "family_cluster_sha256": sha256(cluster_path), "exact_overlap_audit": overlap,
        "accepted_records": len(accepted), "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
        "accepted_unique_substrates": 2 if accepted else 0,
        "status": "excluded_homology" if hits else "accepted_homology_cold_pool",
        "exclusion_reason": "At least one frozen-threshold UniKP development hit" if hits else None,
        "readiness_gate_passes": False,
        "claim_boundary": "Direct-label curation/exclusion evidence only; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-acquire", action="store_true")
    parser.add_argument("--skip-mmseqs", action="store_true")
    args = parser.parse_args()
    curate(acquire=not args.no_acquire, run_homology=not args.skip_mmseqs)
