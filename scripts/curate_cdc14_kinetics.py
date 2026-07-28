from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import subprocess
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "purdue-429k-qe94"
RAW = SOURCE / "raw"
DOI = "10.1016/j.jbc.2024.107644"
DATASET_DOI = "10.4231/429k-qe94"
SOURCE_ID = "purdue-429k-qe94"
PMC_ID = "PMC11407943"
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

DOWNLOADS = {
    "datacite-metadata.json": "https://api.datacite.org/dois/10.4231/429k-qe94",
    "PMC11407943-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11407943/fullTextXML"
    ),
    "PMC11407943-article.pdf": "https://europepmc.org/articles/PMC11407943?pdf=render",
    "PMC11407943-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11407943/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "pubchem-difmup.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/219176/"
        "property/IsomericSMILES,Title/JSON"
    ),
}

# PURR currently redirects the DOI here. Alternate forms are retained because they produce
# independently useful machine-readable evidence when the canonical endpoint is unavailable.
PURR_ENDPOINTS = (
    "https://purr.purdue.edu/publications/4552/1",
    "https://purr.purdue.edu/publications/4552/1/files",
    "https://purr.purdue.edu/api/publications/4552",
    "https://purr.purdue.edu/publications/429k-qe94",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def request(url: str, timeout: int = 180) -> tuple[bytes | None, dict[str, object]]:
    started = utc_now()
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json, application/xml, text/html, */*",
            "User-Agent": "BioCandidateRanker/0.1 temporal-curation (contact: local-audit)",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context()) as response:
            payload = response.read()
            attempt = {
                "requested_url": url,
                "started_at_utc": started,
                "completed_at_utc": utc_now(),
                "outcome": "success",
                "http_status": response.status,
                "final_url": response.url,
                "content_type": response.headers.get_content_type(),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            return payload, attempt
    except urllib.error.HTTPError as error:
        payload = error.read()
        return None, {
            "requested_url": url,
            "started_at_utc": started,
            "completed_at_utc": utc_now(),
            "outcome": "http_error",
            "http_status": error.code,
            "final_url": error.url,
            "error_type": type(error).__name__,
            "error": str(error),
            "response_bytes": len(payload),
            "response_sha256": hashlib.sha256(payload).hexdigest(),
        }
    except Exception as error:
        return None, {
            "requested_url": url,
            "started_at_utc": started,
            "completed_at_utc": utc_now(),
            "outcome": "transport_error",
            "http_status": None,
            "error_type": type(error).__name__,
            "error": str(error),
        }


def acquire() -> list[dict[str, object]]:
    RAW.mkdir(parents=True, exist_ok=True)
    attempts: list[dict[str, object]] = []
    for name, url in DOWNLOADS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size > 100:
            attempts.append(
                {
                    "requested_url": url,
                    "completed_at_utc": utc_now(),
                    "outcome": "reused_local",
                    "path": str(path.relative_to(SOURCE)).replace("\\", "/"),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
            continue
        payload, attempt = request(url)
        attempt["intended_path"] = str(path.relative_to(SOURCE)).replace("\\", "/")
        attempts.append(attempt)
        if payload is not None and len(payload) > 100:
            path.write_bytes(payload)

    for url in PURR_ENDPOINTS:
        payload, attempt = request(url, timeout=60)
        attempts.append(attempt)
        if payload is not None and len(payload) > 100:
            # Preserve responses for later parser development without claiming they are data files.
            index = PURR_ENDPOINTS.index(url) + 1
            path = RAW / f"purr-response-{index}.bin"
            path.write_bytes(payload)
            attempt["saved_path"] = str(path.relative_to(SOURCE)).replace("\\", "/")
        time.sleep(1)

    (SOURCE / "acquisition-attempts.json").write_text(
        json.dumps(
            {
                "source_id": SOURCE_ID,
                "dataset_doi": DATASET_DOI,
                "generated_at_utc": utc_now(),
                "attempts": attempts,
            },
            indent=2,
        )
        + "\n",
        encoding="ascii",
    )
    return attempts


def verify_open_metadata() -> dict[str, object]:
    metadata_path = RAW / "datacite-metadata.json"
    xml_path = RAW / "PMC11407943-fullText.xml"
    supplement_path = RAW / "PMC11407943-supplementaryFiles.zip"
    missing = [path.name for path in (metadata_path, xml_path, supplement_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Open metadata/article acquisition incomplete: {missing}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))["data"]["attributes"]
    if metadata["doi"].lower() != DATASET_DOI or metadata["publicationYear"] != 2024:
        raise ValueError("DataCite dataset identity or year mismatch")
    licenses = {
        item.get("rightsIdentifier", "").lower() for item in metadata.get("rightsList", [])
    }
    if "cc0-1.0" not in licenses:
        raise ValueError("Dataset is not explicitly CC0-1.0 in DataCite")

    xml = xml_path.read_text(encoding="utf-8")
    if f'<article-id pub-id-type="doi">{DOI}</article-id>' not in xml:
        raise ValueError("Europe PMC article DOI mismatch")
    if "creativecommons.org/licenses/by/4.0" not in xml:
        raise ValueError("Europe PMC article does not state CC BY 4.0")
    with zipfile.ZipFile(supplement_path) as archive:
        supplement_members = sorted(archive.namelist())
    return {
        "dataset_title": metadata["titles"][0]["title"],
        "dataset_license": "CC0-1.0",
        "dataset_files_reported": metadata.get("sizes", []),
        "dataset_formats_reported": metadata.get("formats", []),
        "article_license": "CC-BY-4.0",
        "europe_pmc_supplement_members": supplement_members,
    }


def repository_payload_available(attempts: list[dict[str, object]]) -> bool:
    return any(
        attempt.get("requested_url") in PURR_ENDPOINTS
        and attempt.get("outcome") == "success"
        and int(attempt.get("bytes", 0)) > 100
        for attempt in attempts
    )


def write_blocker(
    attempts: list[dict[str, object]], metadata: dict[str, object]
) -> None:
    failed = [
        {
            "url": attempt["requested_url"],
            "outcome": attempt["outcome"],
            "http_status": attempt.get("http_status"),
            "error_type": attempt.get("error_type"),
            "error": attempt.get("error"),
            "completed_at_utc": attempt.get("completed_at_utc"),
        }
        for attempt in attempts
        if attempt.get("requested_url") in PURR_ENDPOINTS
        and attempt.get("outcome") != "success"
    ]
    blocker = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "dataset_doi": DATASET_DOI,
        "article_doi": DOI,
        "generated_at_utc": utc_now(),
        "status": "acquisition_blocked",
        "blocked_resource": "Purdue repository file inventory and deposited data files",
        "reason_code": "repository_transport_failure",
        "reason": (
            "DataCite confirms a 43-file CC0 dataset, but none of the canonical Purdue "
            "record/file endpoints returned retrievable content. Exact deposited rows and their "
            "repository file/row provenance therefore cannot be audited."
        ),
        "failed_attempts": failed,
        "independent_metadata_evidence": metadata,
        "curation_effect": {
            "candidate_records_written": 0,
            "candidate_labels_created": False,
            "human_rows_status": "not_curated",
            "yeast_rows_status": "excluded_known_development_overlap",
            "homology_run_status": "not_run_no_exact_eligible_candidate_sequences",
        },
        "resolution_requirement": (
            "Retrieve the Purdue file manifest and deposited kinetic workbooks, then map each "
            "human measurement to an exact construct sequence and source cell before MMseqs2."
        ),
        "no_guessing_statement": (
            "No kinetic value, variant, sequence, substrate structure, or candidate count was "
            "inferred from figures, AlphaFold predictions, or discovery estimates."
        ),
        "model_predictions_run": False,
    }
    (SOURCE / "blocker-evidence.json").write_text(
        json.dumps(blocker, indent=2) + "\n", encoding="ascii"
    )


def write_provenance(metadata: dict[str, object], blocked: bool) -> None:
    raw_files = sorted(path for path in RAW.iterdir() if path.is_file())
    provenance = {
        "source_id": SOURCE_ID,
        "dataset_doi": DATASET_DOI,
        "stable_record_url": "https://purr.purdue.edu/publications/4552/1",
        "article_doi": DOI,
        "pmc_id": PMC_ID,
        "article_published": "2024-08-08",
        "dataset_published": "2024-07-29",
        "dataset_license": metadata["dataset_license"],
        "article_license": metadata["article_license"],
        "status": "acquisition_blocked" if blocked else "repository_payload_acquired_pending_curation",
        "scope": "homology-cold Homo sapiens Cdc14 rows only",
        "known_exclusion": (
            "All Saccharomyces cerevisiae rows are excluded before normalization because prior "
            "discovery established exact UniKP development overlap."
        ),
        "article_package_use": (
            "Identity, methods, license, and supporting context only while the deposited Purdue "
            "file inventory is inaccessible; not a substitute for repository row provenance."
        ),
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "candidate_records_written": 0,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )


def write_candidate_artifact(blocked: bool) -> None:
    artifact = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "dataset_doi": DATASET_DOI,
        "article_doi": DOI,
        "generated_at_utc": utc_now(),
        "status": "acquisition_blocked" if blocked else "pending_exact_repository_curation",
        "scope": "Homo sapiens direct experimental absolute-kinetics rows only",
        "candidates": [],
        "candidate_count": 0,
        "labels_created": False,
        "excluded_known_groups": [
            {
                "organism": "Saccharomyces cerevisiae",
                "reason": "known_exact_development_sequence_overlap",
                "normalized_rows_created": 0,
            }
        ],
        "human_rows": {
            "status": "not_curated_repository_unavailable",
            "estimated_discovery_count_used": False,
        },
        "model_predictions_run": False,
    }
    (SOURCE / "candidate-artifact.json").write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="ascii"
    )


def run_frozen_wsl_mmseqs(blocked: bool) -> None:
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference hash changed")
    query = SOURCE / "construct_sequences.fasta"
    command = [
        "wsl.exe",
        "-d",
        "Debian",
        "--",
        "mmseqs",
        "easy-search",
        "<WSL_QUERY_FASTA>",
        "<WSL_FROZEN_UNIKP_FASTA>",
        "<WSL_HOMOLOGY_HITS_TSV>",
        "<WSL_TMP_DIR>",
        "--min-seq-id",
        str(MIN_IDENTITY),
        "-c",
        str(COVERAGE),
        "--cov-mode",
        str(COVERAGE_MODE),
        "--format-output",
        "query,target,fident,qcov,tcov,evalue,bits",
    ]
    version_command = ["wsl.exe", "-d", "Debian", "--", "mmseqs", "version"]
    attempted_at = utc_now()
    result = subprocess.run(version_command, capture_output=True, text=True, check=False)
    installed_version = result.stdout.strip() if result.returncode == 0 else None
    if installed_version and installed_version != MMSEQS_VERSION:
        status = "blocked_version_mismatch"
    elif not query.is_file():
        status = "not_run_no_exact_candidate_sequences"
    elif blocked:
        status = "not_run_acquisition_blocked"
    elif result.returncode != 0:
        status = "blocked_mmseqs_unavailable"
    else:
        status = "ready_not_run_parser_pending"
    audit = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "attempted_at_utc": attempted_at,
        "status": status,
        "wsl_distribution": "Debian",
        "version_probe_command": version_command,
        "version_probe_returncode": result.returncode,
        "version_probe_stdout": result.stdout.strip(),
        "version_probe_stderr": result.stderr.strip(),
        "frozen_search_command": command,
        "required_mmseqs_version": MMSEQS_VERSION,
        "observed_mmseqs_version": installed_version,
        "min_identity": MIN_IDENTITY,
        "coverage": COVERAGE,
        "coverage_mode": COVERAGE_MODE,
        "development_target": str(REFERENCE.relative_to(ROOT)).replace("\\", "/"),
        "development_target_sha256": REFERENCE_SHA256,
        "query_fasta_present": query.is_file(),
        "query_sequences": 0,
        "homology_hits_created": False,
        "homology_cold_claimed": False,
        "reason": (
            "Exact human candidate constructs cannot be created before repository acquisition; "
            "the WSL version probe additionally reports that mmseqs is not installed."
        ),
        "model_predictions_run": False,
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Acquire and gate the Cdc14 kinetics source")
    parser.add_argument("--no-download", action="store_true", help="Only validate local raw files")
    args = parser.parse_args()
    SOURCE.mkdir(parents=True, exist_ok=True)
    attempts = [] if args.no_download else acquire()
    if args.no_download:
        attempt_path = SOURCE / "acquisition-attempts.json"
        if attempt_path.is_file():
            attempts = json.loads(attempt_path.read_text(encoding="utf-8"))["attempts"]
    metadata = verify_open_metadata()
    blocked = not repository_payload_available(attempts)
    if blocked:
        write_blocker(attempts, metadata)
    write_provenance(metadata, blocked)
    write_candidate_artifact(blocked)
    run_frozen_wsl_mmseqs(blocked)
    if blocked:
        print("Purdue acquisition blocked; wrote machine-readable evidence and no candidates.")
    else:
        print("Purdue response acquired; exact file inventory still requires parser verification.")


if __name__ == "__main__":
    main()
