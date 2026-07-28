from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC12594388"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1371/journal.pone.0336433"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

URLS = {
    "PMC12594388-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12594388/fullTextXML"
    ),
    "OPC79493.1.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=protein&id=OPC79493.1&rettype=fasta&retmode=text"
    ),
    "pubchem-L-asparagine.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6267/"
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
        with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def validate_raw() -> tuple[str, int, str]:
    xml = (RAW / "PMC12594388-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI,
        "creativecommons.org/licenses/by/4.0/",
        "OPC79493.1",
        "pET28a-AnsSs",
        "40.832",
        "7.361",
        "Substrate (asparagine) concentrations ranging from 0.5 to 20 mM",
    )
    missing = [marker for marker in required if marker not in xml]
    if missing:
        raise ValueError(f"Article XML lacks expected evidence: {missing}")

    sequence = read_fasta(RAW / "OPC79493.1.fasta")
    if len(sequence) != 327 or not sequence.startswith("MRPEDTVPHESTT"):
        raise ValueError("Unexpected OPC79493.1 sequence")
    if set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("OPC79493.1 contains unsupported residues")

    payload = json.loads((RAW / "pubchem-L-asparagine.json").read_text(encoding="utf-8"))
    compound = payload["PropertyTable"]["Properties"][0]
    if compound["CID"] != 6267 or compound["Title"] != "(-)-Asparagine":
        raise ValueError("Unexpected PubChem L-asparagine record")
    return sequence, compound["CID"], compound["SMILES"]


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"Cannot map path to WSL: {resolved}")
    return f"/mnt/{drive}/{resolved.as_posix()[3:]}"


def mmseqs_command(*arguments: str) -> list[str]:
    _, distribution, executable = MMSEQS.split(":", 2)
    return ["wsl", "-d", distribution, executable, *arguments]


def reset_wsl_workspace(path: Path) -> None:
    if not path.exists():
        return
    _, distribution, _ = MMSEQS.split(":", 2)
    subprocess.run(["wsl", "-d", distribution, "rm", "-rf", "--", wsl_path(path)], check=True)
    if path.exists():
        raise RuntimeError(f"Could not reset managed MMseqs workspace: {path}")


def run_mmseqs(query: Path) -> None:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development FASTA is missing or has changed")
    version = subprocess.run(
        mmseqs_command("version"), check=True, capture_output=True, text=True
    ).stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 version {version!r} is not the frozen version")

    HOMOLOGY.mkdir(parents=True, exist_ok=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    search_tmp = HOMOLOGY / "search-tmp"
    family = HOMOLOGY / "family-cluster"
    reset_wsl_workspace(search_tmp)
    reset_wsl_workspace(family)
    family.mkdir()
    subprocess.run(
        mmseqs_command(
            "easy-search", wsl_path(query), wsl_path(REFERENCE), wsl_path(hits),
            wsl_path(search_tmp), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0",
            "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        ),
        check=True,
    )
    subprocess.run(
        mmseqs_command("createdb", wsl_path(query), wsl_path(family / "input_db")), check=True
    )
    subprocess.run(
        mmseqs_command(
            "linclust", wsl_path(family / "input_db"), wsl_path(family / "clusters_db"),
            wsl_path(family / "cluster_tmp"), "--min-seq-id", "0.3", "-c", "0.8",
            "--cov-mode", "0",
        ),
        check=True,
    )
    subprocess.run(
        mmseqs_command(
            "createtsv", wsl_path(family / "input_db"), wsl_path(family / "input_db"),
            wsl_path(family / "clusters_db"), wsl_path(family / "proteins_cluster.tsv"),
        ),
        check=True,
    )


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(*, acquire: bool, run_homology: bool) -> None:
    if acquire:
        download_raw()
    sequence, cid, smiles = validate_raw()
    SOURCE.mkdir(parents=True, exist_ok=True)

    query = SOURCE / "native_sequences.fasta"
    query.write_text(
        ">OPC79493.1 native database sequence; not asserted as the exact tagged assay construct\n"
        + "\n".join(sequence[start : start + 80] for start in range(0, len(sequence), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    if run_homology:
        run_mmseqs(query)

    fields = [
        "candidate_id", "article_doi", "source_row", "organism", "sequence_accession",
        "variable_substrate", "substrate_pubchem_cid", "substrate_isomeric_smiles",
        "kcat_s-1", "kcat_error_s-1", "km_mM", "km_error_mM", "assay_pH",
        "assay_temperature_C", "status_at_normalization",
    ]
    # The reported His-tag was not removed and its complete vector-derived residues are absent.
    # Keep the otherwise eligible measurement out of candidate_records rather than guessing them.
    write_csv(SOURCE / "candidate_records.csv", [], fields)
    excluded = {
        "candidate_id": "scabrispora-asparaginase-001",
        "article_doi": DOI,
        "source_row": "Figure 2B and Table 1; optimum-pH kinetics",
        "organism": "Streptomyces scabrisporus (current NCBI name: Embleya scabrispora)",
        "sequence_accession": "OPC79493.1",
        "variable_substrate": "L-asparagine",
        "substrate_pubchem_cid": cid,
        "substrate_isomeric_smiles": smiles,
        "kcat_s-1": 40.832,
        "kcat_error_s-1": 2.0,
        "km_mM": 7.361,
        "km_error_mM": 1.78,
        "assay_pH": 10,
        "assay_temperature_C": 37,
        "status_at_normalization": "excluded_exact_construct_unresolved",
    }
    write_csv(SOURCE / "exclusions.csv", [excluded], fields)

    hit_path = HOMOLOGY / "homology_hits.tsv"
    cluster_path = HOMOLOGY / "family-cluster" / "proteins_cluster.tsv"
    hit_lines = (
        [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        if hit_path.is_file() else []
    )
    raw_hashes = {
        name: {"sha256": sha256(RAW / name), "size_bytes": (RAW / name).stat().st_size}
        for name in URLS
    }
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC12594388",
        "article_doi": DOI,
        "article_published": "2025-11-07",
        "license": "CC-BY-4.0",
        "reported_direct_finite_kcat_rows": 1,
        "accepted_records": 0,
        "construct_audit": (
            "The paper identifies synthesized OPC79493.1 in pET-28a(+) and HisTrap purification, "
            "but does not report complete vector-derived tag/linker residues or tag removal. "
            "OPC79493.1 is retained only as a native homology query, not an exact assay mapping."
        ),
        "kinetics_source": "raw/PMC12594388-fullText.xml, Figure 2B and Table 1",
        "selection_policy": "Fail closed when the complete assayed construct sequence is unavailable.",
        "raw_files": raw_hashes,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    if hit_path.is_file() and cluster_path.is_file():
        audit = {
            "audited_on": "2026-07-27",
            "source_id": SOURCE_ID,
            "article_doi": DOI,
            "status": "excluded_exact_construct_unresolved",
            "method": "MMseqs2 easy-search and linclust of explicitly non-construct native query",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": {"count": 0, "sha256": sha256(SOURCE / "candidate_records.csv")},
            "reported_direct_finite_kcat_rows": 1,
            "construct_resolved": False,
            "homology_queries_sha256": sha256(query),
            "homology_query_scope": "OPC79493.1 native sequence only; not an assay-construct claim",
            "homology_hit_sequences": 1 if hit_lines else 0,
            "homology_hit_alignments": len(hit_lines),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": 1,
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": 0,
            "readiness_gate_passes": False,
            "claim_boundary": "Construct-resolution and homology prescreen only; no labels or model scores.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-acquire", action="store_true")
    parser.add_argument("--skip-mmseqs", action="store_true")
    args = parser.parse_args()
    curate(acquire=not args.no_acquire, run_homology=not args.skip_mmseqs)
