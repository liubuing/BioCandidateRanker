from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC12995904"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1016/j.jbc.2026.111301"
PMC_ID = "PMC12995904"
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
MMSEQS = "mmseqs"
MIN_IDENTITY = 0.30
COVERAGE = 0.80
COVERAGE_MODE = 0

DOWNLOADS = {
    f"{PMC_ID}-fullText.xml": (
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMC_ID}/fullTextXML"
    ),
    f"{PMC_ID}-article.pdf": f"https://europepmc.org/articles/{PMC_ID}?pdf=render",
    f"{PMC_ID}-supplementaryFiles.zip": (
        f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMC_ID}/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "crossref.json": f"https://api.crossref.org/works/{DOI}",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    attempts = []
    for name, url in DOWNLOADS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size > 100:
            outcome = "reused_local"
        else:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "BioCandidateRanker/0.1 temporal-curation"},
            )
            with urllib.request.urlopen(request, timeout=300) as response:
                path.write_bytes(response.read())
            outcome = "downloaded"
        attempts.append(
            {
                "url": url,
                "path": f"raw/{name}",
                "outcome": outcome,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    (SOURCE / "acquisition-attempts.json").write_text(
        json.dumps({"acquired_at_utc": utc_now(), "files": attempts}, indent=2) + "\n",
        encoding="ascii",
    )


def docx_paragraphs(payload: bytes) -> list[str]:
    temp = RAW / "_docx-inspection.zip"
    temp.write_bytes(payload)
    try:
        with zipfile.ZipFile(temp) as archive:
            root = ET.fromstring(archive.read("word/document.xml"))
    finally:
        temp.unlink()
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    paragraphs = []
    for paragraph in root.findall(".//w:p", namespace):
        text = "".join(node.text or "" for node in paragraph.findall(".//w:t", namespace))
        if text.strip():
            paragraphs.append(text.strip())
    return paragraphs


def extract_supplements() -> list[dict[str, object]]:
    extracted = RAW / "supplements"
    extracted.mkdir(exist_ok=True)
    inventory = []
    with zipfile.ZipFile(RAW / f"{PMC_ID}-supplementaryFiles.zip") as archive:
        for member in archive.infolist():
            if member.is_dir():
                continue
            name = Path(member.filename).name
            payload = archive.read(member)
            path = extracted / name
            path.write_bytes(payload)
            item: dict[str, object] = {
                "archive_member": member.filename,
                "path": f"raw/supplements/{name}",
                "size_bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            if name.lower().endswith(".docx"):
                paragraphs = docx_paragraphs(payload)
                text_path = extracted / f"{name}.txt"
                text_path.write_text("\n".join(paragraphs) + "\n", encoding="utf-8")
                item["extracted_text_path"] = f"raw/supplements/{text_path.name}"
                item["extracted_paragraphs"] = len(paragraphs)
                item["extracted_text_sha256"] = sha256(text_path)
            inventory.append(item)
    (SOURCE / "supplement-inventory.json").write_text(
        json.dumps({"files": inventory}, indent=2) + "\n", encoding="ascii"
    )
    return inventory


def validate_article() -> str:
    xml = (RAW / f"{PMC_ID}-fullText.xml").read_text(encoding="utf-8")
    required = (
        f'<article-id pub-id-type="doi">{DOI}</article-id>',
        "creativecommons.org/licenses/by/4.0",
        "normalized to the corresponding active enzyme concentrations",
        "kcat = 16.9",
        "kcat = 2.3",
        "3.5",
    )
    missing = [token for token in required if token not in xml]
    if missing:
        raise ValueError(f"Article validation failed; missing tokens: {missing}")
    return xml


def inspect_authoritative_text() -> dict[str, object]:
    text_files = sorted((RAW / "supplements").glob("*.docx.txt"))
    joined = "\n".join(path.read_text(encoding="utf-8") for path in text_files)
    evidence = {}
    for term in ("VaPAL2", "OaAEP1b", "GN10-SL", "GYSGSDAL", "GIPNSL"):
        lines = [line for line in joined.splitlines() if term.lower() in line.lower()]
        evidence[term] = lines
    return {"text_files": [path.name for path in text_files], "matches": evidence}


def exact_source_sequences() -> tuple[str, str, str]:
    text = (RAW / "supplements" / "mmc1.docx.txt").read_text(encoding="utf-8")
    match = re.search(r"^>VaPAL2\s*\n([A-Z]+)$", text, re.MULTILINE)
    if not match:
        raise ValueError("Data S1 VaPAL2 FASTA record was not recovered")
    native = match.group(1)
    if len(native) != 483 or native[243] != "I":
        raise ValueError(f"Unexpected VaPAL2 translation length/site: {len(native)}, {native[243]}")

    # Table S4's forward primer translates to His6-Gly-Ser followed by native residue 2;
    # the native initiator Met/signal peptide start is therefore replaced by the tag.
    if "GGCCACCACCACCACCACCACGGAAGTAAGCTATTCGCCGCC" not in text:
        raise ValueError("VaPAL2 His6 construct primer changed or is absent")
    expressed_wt = "HHHHHHGS" + native[1:]
    expressed_i244a = expressed_wt[:250] + "A" + expressed_wt[251:]
    if expressed_wt[250] != "I":
        raise ValueError("Expressed-construct I244 coordinate reconciliation failed")
    return native, expressed_wt, expressed_i244a


def windows_to_wsl(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def write_sequence_evidence(native: str, wt: str, variant: str) -> Path:
    source_fasta = SOURCE / "authoritative_translated_sequences.fasta"
    source_fasta.write_text(
        ">VaPAL2_native_Data_S1 length=483; supplement_Figure_S2_caption_claims_1551bp\n"
        + "\n".join(native[i : i + 80] for i in range(0, len(native), 80))
        + "\n>VaPAL2_His6_expressed_zymogen exact_from_Table_S4 length=490\n"
        + "\n".join(wt[i : i + 80] for i in range(0, len(wt), 80))
        + "\n>VaPAL2_I244A_His6_expressed_zymogen exact_from_Table_S4 length=490\n"
        + "\n".join(variant[i : i + 80] for i in range(0, len(variant), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    query = SOURCE / "construct_sequences.fasta"
    query.write_text(
        ">vapal2-wt exact_His6_expressed_zymogen; active_autocleavage_termini_unreported\n"
        + "\n".join(wt[i : i + 80] for i in range(0, len(wt), 80))
        + "\n>vapal2-i244a exact_His6_expressed_zymogen; active_autocleavage_termini_unreported\n"
        + "\n".join(variant[i : i + 80] for i in range(0, len(variant), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    return query


def write_substrate_evidence() -> None:
    evidence = {
        "name": "GN10-SL",
        "linear_sequence": "GLRRGYSGSNSL",
        "termini": {"n_terminus": "free amine", "c_terminus": "amide"},
        "biln": "H-GLRRGYSGSNSL-NH2",
        "monoisotopic_neutral_mass_Da": 1265.6476,
        "monoisotopic_M_plus_H_Da": 1266.6549,
        "reported_calculated_M_plus_H_Da": 1266,
        "sequence_basis": (
            "The cited source DOI 10.1016/j.jbc.2021.101325 defines GN10-AL as "
            "GLRRGYSGSNAL; the measured source substitutes the required SL exit motif."
        ),
        "terminus_basis": "Rink amide resin in article Experimental procedures",
        "mass_check": "Exact sequence and amidation agree with Supplementary Figure S4 mass",
        "registry_structure_identifier": None,
        "protocol_effect": (
            "Sequence/BILN and mass define the peptide, but the frozen registry only recognizes "
            "PubChem CID structure identifiers; no PubChem record was found for GN10-SL."
        ),
    }
    (SOURCE / "substrate-structure-evidence.json").write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="ascii"
    )


def write_measurement_audit() -> None:
    rows = [
        ("VaPAL2", 100, 3.5, 0.2, 140, 30, 2.5e4),
        ("VaPAL2(I244A)", 50, 16.9, 0.8, 87, 20, 1.9e5),
        ("OaAEP1b(C247A)", 200, 2.3, 0.3, None, None, 1.8e4),
    ]
    audit = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "source_location": "Article Figure 3C and Results paragraph p0085",
        "measurement_type": "direct initial-rate Michaelis-Menten fit",
        "normalization": "kcat = Vmax/[E], normalized to corresponding active enzyme concentration",
        "assay": {
            "variable_substrate": "GN10-SL",
            "substrate_concentration_uM": "50-1000",
            "buffer": "0.1 M sodium acetate, 50 mM NaCl, 1 mM EDTA",
            "pH": "6.5-7.0",
            "temperature_C": 37,
            "sampling_interval_s": "20-30",
            "readout": "UPLC product formation",
            "fit_software": "GraphPad Prism",
            "uncertainty": "fitting error",
        },
        "reported_rows": [
            {
                "enzyme": enzyme,
                "active_enzyme_concentration_nM": concentration,
                "kcat_s-1": kcat,
                "kcat_fitting_error_s-1": kcat_error,
                "km_uM": km,
                "km_fitting_error_uM": km_error,
                "kcat_per_km_M-1_s-1": efficiency,
                "curated_candidate": False,
            }
            for enzyme, concentration, kcat, kcat_error, km, km_error, efficiency in rows
        ],
        "candidate_records_written": 0,
        "reason_not_curated": (
            "Exact processed active-enzyme termini are unreported and GN10-SL lacks a structure "
            "identifier accepted by the frozen registry."
        ),
        "model_predictions_run": False,
    }
    (SOURCE / "measurement-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )


def write_blocker(reason_code: str, reason: str, inspection: dict[str, object]) -> None:
    blocker = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "stable_record_url": f"https://europepmc.org/articles/{PMC_ID}",
        "generated_at_utc": utc_now(),
        "status": "curation_blocked",
        "reason_code": reason_code,
        "reason": reason,
        "authoritative_docx_inspection": inspection,
        "candidate_records_written": 0,
        "labels_created": False,
        "homology_cold_claimed": False,
        "model_predictions_run": False,
    }
    (SOURCE / "blocker-evidence.json").write_text(
        json.dumps(blocker, indent=2) + "\n", encoding="ascii"
    )


def write_provenance(status: str, inspection: dict[str, object]) -> None:
    raw_files = sorted(path for path in RAW.rglob("*") if path.is_file())
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": f"https://europepmc.org/articles/{PMC_ID}",
        "article_doi": DOI,
        "pmc_id": PMC_ID,
        "article_published": "2026-02-17",
        "license": "CC-BY-4.0",
        "status": status,
        "measurement": (
            "Direct initial-rate Michaelis-Menten kcat normalized by reported active-enzyme "
            "concentrations; no values inferred from structural models."
        ),
        "authoritative_assets": inspection,
        "raw_file_sha256": {
            str(path.relative_to(SOURCE)).replace("\\", "/"): sha256(path) for path in raw_files
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )


def run_frozen_homology(query: Path) -> dict[str, object]:
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference hash changed")
    probe = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu-24.04", "--", MMSEQS, "version"],
        capture_output=True,
        text=True,
        check=True,
    )
    version = probe.stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 version {version!r} is not the frozen version")
    HOMOLOGY.mkdir(parents=True, exist_ok=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    temp = HOMOLOGY / "search-tmp"
    command = [
        MMSEQS, "easy-search", windows_to_wsl(query), windows_to_wsl(REFERENCE),
        windows_to_wsl(hits), windows_to_wsl(temp), "--min-seq-id", str(MIN_IDENTITY),
        "-c", str(COVERAGE), "--cov-mode", str(COVERAGE_MODE), "--format-output",
        "query,target,fident,alnlen,qcov,tcov,evalue,bits",
    ]
    subprocess.run(
        ["wsl.exe", "-d", "Ubuntu-24.04", "--", *command], check=True
    )
    hit_lines = [line for line in hits.read_text(encoding="utf-8").splitlines() if line]
    hit_queries = sorted({line.split("\t", 1)[0] for line in hit_lines})
    commands = SOURCE / "homology-commands.txt"
    commands.write_text(
        f"{MMSEQS} version\n{' '.join(command)}\n", encoding="ascii", newline="\n"
    )
    audit = {
        "schema_version": 1,
        "source_id": SOURCE_ID,
        "attempted_at_utc": utc_now(),
        "status": "completed_exclusion_audit",
        "method": "MMseqs2 easy-search",
        "wsl_distribution": "Ubuntu-24.04",
        "mmseqs_binary": MMSEQS,
        "mmseqs_version": version,
        "min_identity": MIN_IDENTITY,
        "coverage": COVERAGE,
        "coverage_mode": COVERAGE_MODE,
        "development_target_sha256": REFERENCE_SHA256,
        "query_sequences": 2,
        "construct_sequences_sha256": sha256(query),
        "homology_hit_sequences": len(hit_queries),
        "homology_hit_query_ids": hit_queries,
        "homology_hits_sha256": sha256(hits),
        "commands_artifact": commands.name,
        "commands_sha256": sha256(commands),
        "candidate_records": 0,
        "accepted_records": 0,
        "homology_cold_claimed": False,
        "reason": (
            "Expressed zymogens were searched conservatively because exact active-core "
            "autocleavage termini are unreported. Any hit excludes the corresponding enzyme."
        ),
        "model_predictions_run": False,
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description="Acquire and gate VaPAL2 absolute kinetics")
    parser.add_argument("--no-download", action="store_true")
    args = parser.parse_args()
    SOURCE.mkdir(parents=True, exist_ok=True)
    if not args.no_download:
        download()
    validate_article()
    extract_supplements()
    inspection = inspect_authoritative_text()
    native, wt, variant = exact_source_sequences()
    query = write_sequence_evidence(native, wt, variant)
    write_substrate_evidence()
    write_measurement_audit()
    audit = run_frozen_homology(query)
    if audit["homology_hit_sequences"]:
        reason_code = "development_homology_hit"
        reason = (
            f"Frozen MMseqs2 found development-corpus hits for {audit['homology_hit_sequences']} "
            "of 2 exact expressed zymogen sequences. The protocol excludes every sequence with "
            "a hit, so no kcat row can enter the homology-cold pool."
        )
    else:
        reason_code = "exact_active_construct_and_registry_structure_id_unresolved"
        reason = (
            "No homology hit was found, but autocatalytic processing produced an active core "
            "whose exact termini were not reported, and GN10-SL has no structure identifier "
            "recognized by the frozen PubChem-based registry."
        )
    write_blocker(reason_code, reason, inspection)
    write_provenance("curation_blocked", inspection)
    print("Acquired VaPAL2 assets; completed pinned homology audit and wrote blocker evidence.")


if __name__ == "__main__":
    main()
