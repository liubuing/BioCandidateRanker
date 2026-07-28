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


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC12959137"
DOI = "10.1128/aac.01069-25"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"

URLS = {
    "PMC12959137-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12959137/fullTextXML",
    "PMC12959137-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12959137/supplementaryFiles",
    "CDO50616.1.fasta": "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=CDO50616.1&db=protein&report=fasta&retmode=text",
    "5UL8.fasta": "https://www.rcsb.org/fasta/entry/5UL8/display",
    "pubchem-cephalothin.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6024/property/Title,IsomericSMILES,InChIKey/JSON",
    "pubchem-imipenem.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/104838/property/Title,IsomericSMILES,InChIKey/JSON",
    "pubchem-ampicillin.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6249/property/Title,IsomericSMILES,InChIKey/JSON",
    "pubchem-nitrocefin.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6436140/property/Title,IsomericSMILES,InChIKey/JSON",
}

STRUCTURES = {
    "cephalothin": (6024, "CC(=O)OCC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)CC3=CC=CS3)SC1)C(=O)O", "XIURVHNZVLADCM-IUODEOHRSA-N"),
    "imipenem": (104838, "C[C@H]([C@@H]1[C@H]2CC(=C(N2C1=O)C(=O)O)SCCN=CN)O", "ZSKVGTPCRGIANV-ZXFLCMHBSA-N"),
    "ampicillin": (6249, "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=CC=C3)N)C(=O)O)C", "AVKUERGKIZMTKX-NJBDSQKTSA-N"),
    "nitrocefin": (6436140, "C1C(=C(N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CS3)C(=O)O)/C=C/C4=C(C=C(C=C4)[N+](=O)[O-])[N+](=O)[O-]", "LHNIIDJCEODSHA-OQRUQETBSA-N"),
}

# Direct finite wild-type KPC-2 Michaelis-Menten fits. Table 3 has fixed 136 mM
# ionic strength; Table 1 and the supplement intentionally retain their own conditions.
KINETICS = (
    ("T1", "cephalothin", "10 mM sodium phosphate", None, 120, 3, 25, 3),
    ("T1", "cephalothin", "50 mM sodium phosphate", None, 150, 4, 86, 7),
    ("T1", "cephalothin", "50 mM sodium phosphate", 150, 300, 10, 110, 10),
    ("T1", "cephalothin", "10 mM HEPES", 150, 390, 10, 27, 4),
    ("T1", "cephalothin", "50 mM Tris-HCl", None, 160, 4, 5.5, 1),
    ("T3", "cephalothin", "10 mM sodium phosphate", 123, 187, 4, 69, 5),
    ("T3", "cephalothin", "25 mM sodium phosphate", 102, 172, 4, 145, 9),
    ("T3", "cephalothin", "50 mM sodium phosphate", 68, 159, 5, 198, 12),
    ("T3", "cephalothin", "75 mM sodium phosphate", 34, 125, 5, 189, 18),
    ("T3", "cephalothin", "100 mM sodium phosphate", 0, 167, 7, 241, 23),
    ("T3", "cephalothin", "10 mM HEPES sodium salt", 127, 186, 7, 26, 5),
    ("T3", "cephalothin", "25 mM HEPES sodium salt", 112, 197, 9, 36, 7),
    ("T3", "cephalothin", "50 mM HEPES sodium salt", 87, 208, 6, 58, 5),
    ("T3", "cephalothin", "75 mM HEPES sodium salt", 62, 245, 16, 79, 15),
    ("T3", "cephalothin", "100 mM HEPES sodium salt", 37, 272, 12, 81, 10),
    ("S3", "imipenem", "10 mM sodium phosphate", None, 36, 1, 110, 6),
    ("S3", "imipenem", "50 mM sodium phosphate", None, 39, 1, 290, 20),
    ("S3", "imipenem", "50 mM sodium phosphate", 150, 45, 2, 150, 20),
    ("S3", "imipenem", "10 mM HEPES", 150, 24, 1, 48, 6),
    ("S3", "imipenem", "50 mM Tris-HCl", None, 53, 2, 22, 4),
    ("S5", "ampicillin", "10 mM sodium phosphate", None, 160, 6, 140, 10),
    ("S5", "ampicillin", "50 mM sodium phosphate", None, 160, 10, 240, 30),
    ("S5", "ampicillin", "50 mM sodium phosphate", 150, 170, 20, 290, 80),
    ("S5", "ampicillin", "10 mM HEPES", 150, 210, 8, 72, 9),
    ("S5", "ampicillin", "50 mM Tris-HCl", None, 170, 8, 56, 10),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_raw() -> list[dict[str, object]]:
    RAW.mkdir(parents=True, exist_ok=True)
    attempts = []
    for name, url in URLS.items():
        path = RAW / name
        if not path.is_file() or not path.stat().st_size:
            request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
            with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
                shutil.copyfileobj(response, output)
        attempts.append({"url": url, "status": "acquired", "path": f"raw/{name}"})
    api_url = "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12959137"
    try:
        urllib.request.urlopen(api_url, timeout=60)
        api_status = "resolved"
    except urllib.error.HTTPError as error:
        api_status = f"HTTP {error.code}"
    attempts.append({"url": api_url, "status": api_status, "note": "Stable accession page exists; API absence does not replace the Europe PMC deposit mirror."})
    return attempts


def fasta_sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines() if line and not line.startswith(">"))


def validate_raw() -> tuple[str, str]:
    xml = (RAW / "PMC12959137-fullText.xml").read_text(encoding="utf-8")
    markers = (DOI, "creativecommons.org/licenses/by/4.0", "pET28a", "TEV protease", "removal of the His-tag", "nonlinear Michaelis", "at least two independent replicates")
    missing = [marker for marker in markers if marker not in xml]
    if missing:
        raise ValueError(f"Primary article evidence changed or is incomplete: {missing}")
    with zipfile.ZipFile(RAW / "PMC12959137-supplementaryFiles.zip") as archive:
        member = next(name for name in archive.namelist() if name.endswith("aac.01069-25-s0001.docx"))
        supplement = RAW / "aac.01069-25-s0001.docx"
        supplement.write_bytes(archive.read(member))
    with zipfile.ZipFile(supplement) as archive:
        document = " ".join(ET.fromstring(archive.read("word/document.xml")).itertext())
    for marker in ("Table S3", "Table S 5", "Michaelis", "imipenem", "ampicillin"):
        if marker not in document:
            raise ValueError(f"Supplement marker absent: {marker}")
    precursor = fasta_sequence(RAW / "CDO50616.1.fasta")
    if len(precursor) != 293 or not precursor.startswith("MSLYRRLVLLSCLSWPLAGFSATA"):
        raise ValueError("CDO50616.1 precursor identity changed")
    mature_component = precursor[24:]
    tagged_5ul8 = fasta_sequence(RAW / "5UL8.fasta")
    if not tagged_5ul8.endswith(mature_component) or len(tagged_5ul8) != 290:
        raise ValueError("5UL8 tagged KPC-2 does not map to the 269-aa mature component")
    for name, expected in STRUCTURES.items():
        observed = json.loads((RAW / f"pubchem-{name}.json").read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if (observed["CID"], observed["SMILES"], observed["InChIKey"]) != expected:
            raise ValueError(f"PubChem identity changed for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development FASTA is absent or changed")
    return mature_component, tagged_5ul8


def command(executable: str, *arguments: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *arguments]
    _, distribution, binary = executable.split(":", 2)
    return ["wsl", "-d", distribution, binary, *arguments]


def tool_path(path: Path, executable: str) -> str:
    if not executable.startswith("wsl:"):
        return str(path)
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def run_checked(arguments: list[str]) -> str:
    result = subprocess.run(arguments, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(arguments)}\n{result.stderr}")
    return (result.stdout or result.stderr).strip()


def run_homology(query: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 {version!r} differs from frozen version")
    if HOMOLOGY.exists():
        if executable.startswith("wsl:"):
            run_checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(HOMOLOGY, executable)])
        else:
            shutil.rmtree(HOMOLOGY)
    HOMOLOGY.mkdir(parents=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    cluster_prefix = HOMOLOGY / "family-cluster" / "proteins"
    cluster_prefix.parent.mkdir()
    common = ("--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads))
    run_checked(command(executable, "easy-search", tool_path(query, executable), tool_path(REFERENCE, executable), tool_path(hits, executable), tool_path(HOMOLOGY / "search-tmp", executable), *common, "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"))
    run_checked(command(executable, "easy-linclust", tool_path(query, executable), tool_path(cluster_prefix, executable), tool_path(HOMOLOGY / "cluster-tmp", executable), *common))
    return hits, cluster_prefix.with_name("proteins_cluster.tsv"), version


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(acquire: bool, mmseqs: str, threads: int) -> None:
    attempts = download_raw() if acquire else []
    mature, tagged = validate_raw()
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "construct_sequences.fasta"
    sequences = {"kpc2-mature-catalytic-component": mature, "kpc2-5ul8-his-thrombin-tagged-comparator": tagged}
    query.write_text("".join(f">{name}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n" for name, sequence in sequences.items()), encoding="ascii", newline="\n")

    # Homology is deliberately completed before any kinetic row is normalized.
    hits, cluster, version = run_homology(query, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="utf-8").splitlines() if line]
    hit_ids = {row[0] for row in hit_rows}
    excluded_homology = "kpc2-mature-catalytic-component" in hit_ids
    exact_assay_construct_resolved = False
    status = "excluded_homology_and_exact_construct" if excluded_homology else "excluded_exact_construct"

    rows = []
    for index, (table, substrate, buffer_name, nacl_mm, kcat, kcat_error, km, km_error) in enumerate(KINETICS, 1):
        cid, smiles, inchikey = STRUCTURES[substrate]
        table_label = {"T1": "1", "T3": "3", "S3": "S3", "S5": "S5"}[table]
        rows.append({
            "candidate_id": f"kpc2-buffer-{index:03d}", "article_doi": DOI,
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12959137",
            "source_file": "raw/PMC12959137-fullText.xml" if table in {"T1", "T3"} else "raw/aac.01069-25-s0001.docx",
            "source_row": f"Table {table_label}; wild-type KPC-2; {substrate}; {buffer_name}" + (f" + {nacl_mm} mM NaCl" if nacl_mm is not None else ""),
            "organism": "Klebsiella pneumoniae", "enzyme_identity": "KPC-2 beta-lactamase",
            "sequence_id": "kpc2-mature-catalytic-component", "construct": "pET28a mature-protein insert; N-terminal His tag removed by TEV before assay; exact insert boundary and residual TEV scar not reported",
            "variable_substrate": substrate, "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles, "substrate_inchikey": inchikey,
            "endpoint": "kcat_s-1", "kcat_s-1": kcat, "kcat_fit_error_s-1": kcat_error,
            "km_uM": km, "km_fit_error_uM": km_error,
            "fit_error_type": "reported nonlinear-fit error; statistical type not further specified",
            "assay_temperature_C": 25, "assay_pH": 7.0, "assay_buffer": buffer_name,
            "added_NaCl_mM": "" if nacl_mm is None else nacl_mm,
            "ionic_strength_mM": 136 if table == "T3" else "not held constant/not reported as one value",
            "BSA_mg_per_mL": 0.1, "fit": "direct nonlinear Michaelis-Menten fit of initial rates from at least two independent replicates" if table != "T3" else "direct nonlinear Michaelis-Menten fit of initial rates from at least four independent replicates",
            "saturation_evidence": "finite kcat and KM from a substrate-series nonlinear Michaelis-Menten fit",
            "status_at_normalization": status,
        })
    write_csv(SOURCE / "excluded_records.csv", rows)
    (SOURCE / "candidate_records.csv").write_text("", encoding="ascii")

    exclusions = [
        {"source_scope": "KPC-2 T237G, CTX-M-14, and TEM-1 table rows", "endpoint": "kcat/KM/KM", "reason": "outside requested wild-type KPC-2 scope", "candidate_label_created": False},
        {"source_scope": "Table 4 and Figure 3", "endpoint": "apparent Ki from Morrison inhibition fits", "reason": "inhibition endpoint, not kcat", "candidate_label_created": False},
        {"source_scope": "Figure 4", "endpoint": "IC50 for avibactam/clavulanate", "reason": "inhibition endpoint, not kcat", "candidate_label_created": False},
        {"source_scope": "global phosphate/HEPES fits", "endpoint": "competitive/mixed inhibition Ki", "reason": "global inhibition endpoint, not row-level kcat", "candidate_label_created": False},
        {"source_scope": "nitrocefin fixed 50 uM reporter experiments", "endpoint": "initial velocity versus buffer concentration", "reason": "no substrate saturation fit or kcat", "candidate_label_created": False},
    ]
    write_csv(SOURCE / "exclusions.csv", exclusions)
    sequence_rows = [
        {"sequence_id": "kpc2-precursor-CDO50616.1", "length_aa": 293, "role": "native precursor reference", "mapping": "NCBI CDO50616.1; residues 1-24 signal peptide", "sha256": hashlib.sha256(("MSLYRRLVLLSCLSWPLAGFSATA" + mature).encode()).hexdigest()},
        {"sequence_id": "kpc2-mature-catalytic-component", "length_aa": len(mature), "role": "homology component", "mapping": "CDO50616.1 residues 25-293; mature core encoded by the reported mature-protein insert", "sha256": hashlib.sha256(mature.encode()).hexdigest()},
        {"sequence_id": "kpc2-5ul8-his-thrombin-tagged-comparator", "length_aa": len(tagged), "role": "tagged structural comparator only", "mapping": "PDB 5UL8 chain A: MGSSHHHHHHSSGLVPRGSHM + CDO50616.1 residues 25-293; not asserted to be the TEV-cleaved assay chain", "sha256": hashlib.sha256(tagged.encode()).hexdigest()},
    ]
    write_csv(SOURCE / "sequence_mapping.csv", sequence_rows)
    raw_hashes = {path.name: {"url": URLS.get(path.name), "size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    (SOURCE / "acquisition-attempts.json").write_text(json.dumps(attempts, indent=2) + "\n", encoding="ascii")

    best_hits = sorted(hit_rows, key=lambda row: float(row[7]), reverse=True)
    blocker = {
        "schema_version": 1, "source_id": SOURCE_ID, "article_doi": DOI,
        "status": status, "accepted_records": 0,
        "blocker_codes": (["FROZEN_DEVELOPMENT_HOMOLOGY"] if excluded_homology else []) + ["EXACT_ASSAY_CONSTRUCT_UNRESOLVED"],
        "otherwise_direct_finite_kcat_rows": len(rows),
        "evidence": [
            f"Pinned MMseqs2 found {len(hit_rows)} qualifying alignments across {len(hit_ids)} query sequences; complete alignments are retained in homology/homology_hits.tsv.",
            "The mature 269-aa KPC-2 catalytic component is sequence-resolved from CDO50616.1 and the article's mature-protein insert statement, and was searched independently so short affinity-tag differences cannot dilute coverage.",
            "The article states that the pET28a N-terminal His tag was removed with TEV, but does not give insert boundaries, cloning junction sequence, or the residual post-cleavage scar. The exact assayed polypeptide therefore cannot be asserted.",
        ],
        "top_homology_hits": [dict(zip(("query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"), row)) for row in best_hits[:10]],
        "resolution_required": "No construct clarification can reverse a qualifying frozen-development hit; exact plasmid/insert sequence would additionally be required to assert the assayed chain.",
        "model_predictions_run": False,
    }
    (SOURCE / "blocker-evidence.json").write_text(json.dumps(blocker, indent=2) + "\n", encoding="ascii")
    cluster_rows = [line.split("\t") for line in cluster.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust of the mature KPC-2 catalytic component and an explicitly labeled tagged structural comparator",
        "mmseqs_version": version, "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "query_sequences": len(sequences),
        "construct_sequences_sha256": sha256(query), "homology_hit_query_sequences": len(hit_ids),
        "homology_hit_alignments": len(hit_rows), "homology_hits_sha256": sha256(hits),
        "candidate_mmseqs_families": len({row[0] for row in cluster_rows}), "family_cluster_sha256": sha256(cluster),
        "exact_assay_construct_resolved": exact_assay_construct_resolved,
        "direct_finite_wildtype_kcat_rows": len(rows), "accepted_records": 0, "status": status,
        "readiness_gate_passes": False, "claim_boundary": "Curation/exclusion evidence only; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12959137",
        "article_doi": DOI, "article_published": "2026-01-30", "license": "CC-BY-4.0",
        "reported_direct_finite_wildtype_kcat_rows": len(rows), "accepted_records": 0,
        "selection_policy": "All direct finite wild-type KPC-2 kcat values from substrate-series nonlinear Michaelis-Menten fits in main Tables 1 and 3 and Supplementary Tables S3 and S5; no digitization, prediction, or value-based selection.",
        "construct_resolution": "The 269-aa mature catalytic component maps exactly to CDO50616.1 residues 25-293. pET28a provided an N-terminal His tag and TEV site, and the tag was removed. The article omits exact cloning junctions and residual TEV scar, so the exact assay chain fails closed. PDB 5UL8 is retained only as a distinct His/thrombin-tagged comparator.",
        "saturation_audit": "Every retained row reports finite kcat and KM from direct nonlinear Michaelis-Menten substrate-series fits. The source does not publish concentration grids for every fit, so saturation is evidenced by the direct finite asymptotic fit rather than an independently calculated maximum-substrate/KM ratio.",
        "structure_mapping": {name: {"pubchem_cid": values[0], "isomeric_smiles": values[1], "inchikey": values[2]} for name, values in STRUCTURES.items()},
        "final_disposition": status, "raw_file_hashes": "raw-file-hashes.json",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-acquire", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    curate(not args.no_acquire, args.mmseqs, args.threads)
