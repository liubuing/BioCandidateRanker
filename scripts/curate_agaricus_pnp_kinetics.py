from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC12028538"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.3390/jof11040268"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

URLS = {
    "PMC12028538-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12028538/fullTextXML"
    ),
    "ON012779.1.gb": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nuccore&id=ON012779&rettype=gb&retmode=text"
    ),
    "WEP24514.1.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=protein&id=WEP24514.1&rettype=fasta&retmode=text"
    ),
    "XP_007331146.1.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=protein&id=XP_007331146.1&rettype=fasta&retmode=text"
    ),
    "pubchem-inosine.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/135398641/"
        "property/IsomericSMILES,Title/JSON"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if not path.exists():
            request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as output:
                shutil.copyfileobj(response, output)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def validate_sources() -> tuple[str, str, int, str]:
    xml = (RAW / "PMC12028538-fullText.xml").read_text(encoding="utf-8")
    required_article_text = (
        DOI,
        "ON012779",
        "007331146",
        "14.38",
        "7.1",
        "2.02",
        "Creative Commons Attribution",
    )
    missing = [text for text in required_article_text if text not in xml]
    if missing:
        raise ValueError(f"Article XML lacks expected evidence: {missing}")

    construct = read_fasta(RAW / "WEP24514.1.fasta")
    native = read_fasta(RAW / "XP_007331146.1.fasta")
    if len(construct) != 311 or construct != "EF" + native[1:] + "HHHHHH":
        raise ValueError("Deposited WEP24514.1 construct does not match EF + native residues 2-304 + His6")

    substrate = json.loads((RAW / "pubchem-inosine.json").read_text(encoding="utf-8"))
    compound = substrate["PropertyTable"]["Properties"][0]
    if compound["CID"] != 135398641 or compound["Title"] != "Inosine":
        raise ValueError("Unexpected PubChem inosine record")
    return construct, native, compound["CID"], compound["SMILES"]


def run_mmseqs(query_path: Path) -> None:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development FASTA is missing or has changed")
    executable = shutil.which("mmseqs")
    if not executable:
        return
    version = subprocess.run(
        [executable, "version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs version is not frozen version: {version}")

    HOMOLOGY.mkdir(parents=True, exist_ok=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    search_tmp = HOMOLOGY / "search-tmp"
    subprocess.run(
        [
            executable,
            "easy-search",
            str(query_path),
            str(REFERENCE),
            str(hits),
            str(search_tmp),
            "--min-seq-id", "0.3",
            "-c", "0.8",
            "--cov-mode", "0",
            "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        ],
        check=True,
    )
    family = HOMOLOGY / "family-cluster"
    family.mkdir(exist_ok=True)
    subprocess.run([executable, "createdb", str(query_path), str(family / "input_db")], check=True)
    subprocess.run(
        [
            executable,
            "linclust",
            str(family / "input_db"),
            str(family / "clusters_db"),
            str(family / "cluster_tmp"),
            "--min-seq-id", "0.3",
            "-c", "0.8",
            "--cov-mode", "0",
        ],
        check=True,
    )
    subprocess.run(
        [
            executable,
            "createtsv",
            str(family / "input_db"),
            str(family / "input_db"),
            str(family / "clusters_db"),
            str(family / "proteins_cluster.tsv"),
        ],
        check=True,
    )


def write_outputs(*, acquire: bool = True) -> None:
    if acquire:
        download_raw()
    construct, native, cid, smiles = validate_sources()
    SOURCE.mkdir(parents=True, exist_ok=True)

    fasta_path = SOURCE / "construct_sequences.fasta"
    fasta_path.write_text(
        ">abpnp-wep24514.1 WEP24514.1 ON012779.1 | exact deposited codon-optimized construct\n"
        + "\n".join(construct[start : start + 80] for start in range(0, len(construct), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    run_mmseqs(fasta_path)

    hit_path = HOMOLOGY / "homology_hits.tsv"
    cluster_path = HOMOLOGY / "family-cluster" / "proteins_cluster.tsv"
    homology_complete = hit_path.is_file() and cluster_path.is_file()
    hit_lines = (
        [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        if hit_path.is_file()
        else []
    )
    status = "pending_homology"
    if homology_complete:
        status = "excluded_homology" if hit_lines else "accepted_homology_cold_pool"

    record = {
        "candidate_id": "abpnp-001",
        "article_doi": DOI,
        "stable_record_url": "https://europepmc.org/articles/PMC12028538",
        "source_file": "raw/PMC12028538-fullText.xml",
        "source_section": "Sections 2.5 and 3.4; Conclusions",
        "source_row": "AbPNP; inosine; pH 7.0; 60 C",
        "organism": "Agaricus bisporus var. burnettii JB137-S8",
        "enzyme_identity": "purine nucleoside phosphorylase (EC 2.4.2.1)",
        "sequence_accession": "ON012779.1; WEP24514.1",
        "sequence_id": "abpnp-wep24514.1",
        "construct": "deposited codon-optimized CDS: EF + XP_007331146.1 residues 2-304 + C-terminal His6",
        "variable_substrate": "inosine",
        "substrate_pubchem_cid": cid,
        "substrate_isomeric_smiles": smiles,
        "endpoint": "kcat_s-1",
        "kcat_s-1": 14.38,
        "km_uM": 7.10,
        "kcat_per_km_M-1_s-1": 2.02e6,
        "assay_pH": 7.0,
        "assay_temperature_C": 60,
        "replicates": 3,
        "error_type": "not reported for kinetic parameters",
        "fit_method": "classical Michaelis-Menten equation",
        "status_at_normalization": status,
    }
    csv_path = SOURCE / "candidate_records.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record))
        writer.writeheader()
        writer.writerow(record)

    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC12028538",
        "article_doi": DOI,
        "article_published": "2025-03-31",
        "license": "CC-BY-4.0",
        "license_scope": "Article and embedded measurement text/table content",
        "record_count": 1,
        "unique_construct_sequences": 1,
        "unique_substrates": 1,
        "kinetics_source": "raw/PMC12028538-fullText.xml, Sections 2.5 and 3.4; Conclusions",
        "sequence_sources": [
            "raw/ON012779.1.gb (author-deposited optimized CDS)",
            "raw/WEP24514.1.fasta (author-deposited construct translation)",
            "raw/XP_007331146.1.fasta (native reference translation)",
        ],
        "sequence_mapping": {
            "assayed_enzyme": "Article states codon-optimized AbPNP based on XM_007331084, deposited as ON012779, and confirmed by DNA sequencing.",
            "exact_audit_sequence": "WEP24514.1, the ON012779.1 CDS translation (311 aa)",
            "native_reference": "XP_007331146.1 (304 aa)",
            "relationship": "WEP24514.1 = EF + XP_007331146.1 residues 2-304 + HHHHHH",
            "article_deposit_discrepancy": "Article describes an N-terminal His6 tag; ON012779.1 encodes a C-terminal His6 tag and an EF N-terminus. The deposited exact translation is used without correction.",
        },
        "substrate_mapping": {
            "article_name": "inosine (hypoxanthine nucleoside)",
            "pubchem_cid": cid,
            "source": "raw/pubchem-inosine.json",
        },
        "assay": {
            "pH": 7.0,
            "temperature_C": 60,
            "inosine_range_mM": "0-16",
            "replicates": 3,
            "kinetic_parameter_uncertainty": None,
        },
        "normalization": "Author-reported kcat in s-1 and Km interpreted as uM from Conclusions and comparison text; no numerical conversion or refitting.",
        "quality_caveat": "Section 3.4 prints Km as 7.1 micromol, while Conclusions prints 7.10 uM; the concentration unit is retained from the unambiguous Conclusions statement.",
        "excluded_data": [
            "Guanosine and adenosine substrate-specificity values are relative activity, not absolute kinetics.",
            "Beer/wort purine concentrations are bulk process outcomes, not enzyme kinetic constants.",
        ],
        "raw_sha256": {name: sha256(RAW / name) for name in URLS},
        "homology_status": status,
        "final_disposition": "excluded_by_frozen_homology_gate" if hit_lines else status,
        "exact_sequence_substrate_overlap": "not evaluated after mandatory homology exclusion" if hit_lines else "pending global overlap audit",
        "blockers": [] if homology_complete else ["Frozen MMseqs2 executable/evidence unavailable; record remains pending_homology."],
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    if homology_complete:
        (SOURCE / "prescreen-audit.json").unlink(missing_ok=True)
        families = {
            line.split("\t", 1)[0]
            for line in cluster_path.read_text(encoding="utf-8").splitlines()
            if line
        }
        accepted = 1 if status == "accepted_homology_cold_pool" else 0
        audit = {
            "audited_on": "2026-07-27",
            "source_id": SOURCE_ID,
            "article_doi": DOI,
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "search_command": "mmseqs easy-search construct_sequences.fasta unique_proteins.fasta homology_hits.tsv search-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits",
            "cluster_command": "mmseqs linclust input_db clusters_db cluster_tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0",
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": {"count": 1, "sha256": sha256(csv_path)},
            "unique_sequences": 1,
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": 1 if hit_lines else 0,
            "homology_hit_alignments": len(hit_lines),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len(families),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": accepted,
            "accepted_unique_sequences": accepted,
            "accepted_unique_substrates": accepted,
            "status": "excluded_homology" if hit_lines else "accepted_homology_cold_pool",
            "exclusion_reason": "At least one frozen-threshold UniKP development hit" if hit_lines else None,
            "exact_sequence_substrate_overlap": "not evaluated after mandatory homology exclusion" if hit_lines else "pending global overlap audit",
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(
            json.dumps(audit, indent=2) + "\n", encoding="ascii"
        )
    else:
        blocked_audit = {
            "schema_version": 1,
            "audited_on": "2026-07-27",
            "source_id": SOURCE_ID,
            "article_doi": DOI,
            "status": "blocked_pending_frozen_mmseqs",
            "method": "MMseqs2 easy-search and linclust",
            "required_mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": {"count": 1, "sha256": sha256(csv_path)},
            "construct_sequences_sha256": sha256(fasta_path),
            "accepted_records": 0,
            "blockers": ["Frozen MMseqs2 executable was not available; no homology conclusion was inferred."],
            "claim_boundary": "Pending curation only; no model predictions were generated.",
        }
        (SOURCE / "prescreen-audit.json").write_text(
            json.dumps(blocked_audit, indent=2) + "\n", encoding="ascii"
        )


if __name__ == "__main__":
    write_outputs()
