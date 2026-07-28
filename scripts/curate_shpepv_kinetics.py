from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC11742315"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
DOI = "10.1016/j.jbc.2024.107928"
TAG = "ENLYFQGHHHHHHHHHH"
REFERENCE = (
    ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723"
    / "homology" / "unikp_reference.fasta"
)
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"

SUBSTRATES = {
    "Cys-Gly-3M3SH": (
        59475291,
        "CCCC(C)(CCO)SC[C@@H](C(=O)NCC(=O)O)N",
    ),
    "Cys-Gly": (439498, "C([C@@H](C(=O)NCC(=O)O)N)S"),
}

# Unique finite ShPepV rows from main Tables 1 and 2. The WT malodor-substrate row
# appears in both tables and is retained once.
ROWS = (
    ("WT", "Cys-Gly-3M3SH", 3.41, 0.07, 1.29, 0.04, "Tables 1 and 2"),
    ("WT", "Cys-Gly", 405.40, 23.56, 2.19, 0.22, "Table 1"),
    ("A174V", "Cys-Gly-3M3SH", 2.30, 0.05, 1.58, 0.14, "Table 2"),
    ("I413W", "Cys-Gly-3M3SH", 3.01, 0.11, 7.79, 0.67, "Table 2"),
    ("M439A", "Cys-Gly-3M3SH", 4.71, 0.05, 0.70, 0.01, "Table 2"),
    ("D437A", "Cys-Gly-3M3SH", 5.79, 0.06, 0.37, 0.03, "Table 2"),
    ("E175A", "Cys-Gly-3M3SH", 4.73, 0.21, 0.85, 0.11, "Table 2"),
    ("M271R", "Cys-Gly-3M3SH", 2.96, 0.03, 2.31, 0.08, "Table 2"),
)
MUTATIONS = {
    "WT": (),
    "A174V": ((174, "A", "V"),),
    "I413W": ((413, "I", "W"),),
    "M439A": ((439, "M", "A"),),
    "D437A": ((437, "D", "A"),),
    "E175A": ((175, "E", "A"),),
    "M271R": ((271, "M", "R"),),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def mutate(sequence: str, changes: tuple[tuple[int, str, str], ...]) -> str:
    result = sequence
    for position, expected, replacement in changes:
        if result[position - 1] != expected:
            raise ValueError(f"Expected {expected}{position}, found {result[position - 1]}")
        result = result[: position - 1] + replacement + result[position:]
    return result


def verify_raw() -> str:
    required = (
        "PMC11742315-fullText.xml",
        "PMC11742315-supplementaryFiles.zip",
        "jbc.2024.107928.pdf",
        "mmc1.docx",
        "biostudies-mmc1.docx",
        "biostudies-metadata.json",
        "QKQ29470.1.fasta",
        "A0A657M1L9.unisave.txt",
        "pubchem-substrates.json",
    )
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {missing}")
    xml = (RAW / "PMC11742315-fullText.xml").read_text(encoding="utf-8")
    for evidence in (
        DOI,
        "creativecommons.org/licenses/by/4.0",
        "C-terminal tag (-ENLYFQGHHHHHHHHHH)",
        "Menten plots were generated to derive the",
        "5.79",
        "405.40",
    ):
        if evidence not in xml:
            raise ValueError(f"Required article evidence is absent: {evidence}")
    metadata = json.loads((RAW / "biostudies-metadata.json").read_text(encoding="utf-8"))
    if metadata["accno"] != "S-EPMC11742315" or metadata["section"]["files"] != [
        {"path": "mmc1.docx", "size": 1114856, "type": "file"}
    ]:
        raise ValueError("Unexpected BioStudies identity or file manifest")
    if (RAW / "biostudies-mmc1.docx").stat().st_size != 1114856:
        raise ValueError("Downloaded BioStudies file does not match the deposited size")
    native = read_fasta(RAW / "QKQ29470.1.fasta")
    if len(native) != 469:
        raise ValueError("Unexpected QKQ29470.1 length")
    unisave = (RAW / "A0A657M1L9.unisave.txt").read_text(encoding="utf-8")
    sequence_record = unisave.split("SQ   SEQUENCE", 1)[1].split("//", 1)[0]
    sequence_lines = sequence_record.splitlines()[1:]
    archived = "".join("".join(character for character in line if character.isalpha()) for line in sequence_lines)
    if "DR   EMBL; CP054550; QKQ29470.1" not in unisave or native != archived:
        raise ValueError("A0A657M1L9 archive does not corroborate QKQ29470.1")
    compounds = json.loads((RAW / "pubchem-substrates.json").read_text(encoding="utf-8"))
    observed = {item["CID"]: item["SMILES"] for item in compounds["PropertyTable"]["Properties"]}
    if observed != {cid: smiles for cid, smiles in SUBSTRATES.values()}:
        raise ValueError("PubChem substrate mapping changed")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference changed")
    return native


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


def run_homology(fasta: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(mmseqs_command(executable, "version")).splitlines()[-1]
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
    common = ("--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads))
    run_checked(mmseqs_command(executable, "easy-search", tool_path(fasta, executable), tool_path(REFERENCE, executable), tool_path(hits, executable), tool_path(workspace / "search-tmp", executable), *common, "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits"))
    run_checked(mmseqs_command(executable, "easy-linclust", tool_path(fasta, executable), tool_path(prefix, executable), tool_path(workspace / "cluster-tmp", executable), *common))
    return hits, prefix.with_name("proteins_cluster.tsv"), version


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(mmseqs: str, threads: int) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    native = verify_raw()
    sequences = {variant: mutate(native, changes) + TAG for variant, changes in MUTATIONS.items()}
    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        for variant, sequence in sequences.items():
            handle.write(f">shpepv-{variant.lower()} QKQ29470.1/A0A657M1L9 + C-terminal {TAG} | {variant}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    hits, cluster, version = run_homology(fasta, mmseqs, threads)
    hit_lines = [line for line in hits.read_text(encoding="utf-8").splitlines() if line]
    rows = []
    for index, (variant, substrate, kcat, kcat_sd, km, km_sd, table) in enumerate(ROWS, 1):
        cid, smiles = SUBSTRATES[substrate]
        coupled = substrate == "Cys-Gly-3M3SH"
        rows.append({
            "exclusion_id": f"shpepv-{index:03d}", "article_doi": DOI,
            "source_table": table, "source_row": f"ShPepV {variant}; {substrate}",
            "organism": "Staphylococcus hominis", "sequence_id": f"shpepv-{variant.lower()}",
            "construct": f"QKQ29470.1/A0A657M1L9 {variant} + C-terminal {TAG}",
            "variable_substrate": substrate, "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles, "reported_kcat_s-1": kcat,
            "reported_kcat_sd_s-1": kcat_sd, "reported_km_mM": km,
            "reported_km_sd_mM": km_sd, "variable_substrate_saturation_fit": True,
            "status_at_normalization": "excluded_coupled_assay_capacity_unresolved" if coupled else "excluded_enzyme_concentration_provenance_unresolved",
            "exclusion_reason": "ShPatB-DTNB coupling capacity was not demonstrated non-rate-limiting across the substrate series" if coupled else "The Cd-ninhydrin method reports an unresolved 20 mM enzyme concentration",
            "candidate_label_created": False,
        })
    fields = list(rows[0])
    write_csv(SOURCE / "excluded_records.csv", rows, fields)
    write_csv(SOURCE / "candidate_records.csv", [], [
        "candidate_id", "article_doi", "sequence_id", "substrate_pubchem_cid",
        "kcat_s-1", "status_at_normalization",
    ])
    blocker = {
        "schema_version": 1, "source_id": SOURCE_ID, "article_doi": DOI,
        "status": "blocked_fail_closed", "blocker_code": "COUPLED_ASSAY_CAPACITY_UNRESOLVED",
        "reported_unique_finite_saturation_fit_rows": len(rows), "accepted_records": 0,
        "evidence": [
            "The steady-state assay varies Cys-Gly-3M3SH and detects product through added ShPatB and DTNB.",
            "The article states the 1:1 signal stoichiometry but reports no ShPatB concentration-series, coupling-enzyme capacity control, or proof that coupling remains faster than PepV at every fitted substrate concentration.",
            "The independent Cys-Gly Cd-ninhydrin fit is retained but the source reports an implausible methods concentration of 20 mM enzyme; no candidate label is released without resolution.",
        ],
        "resolution_required": "Primary evidence that product detection is non-rate-limiting over every fitted rate, plus clarification of the Cys-Gly enzyme concentration.",
        "model_predictions_run": False, "recorded_on": str(date.today()),
    }
    (SOURCE / "blocker-evidence.json").write_text(json.dumps(blocker, indent=2) + "\n", encoding="ascii")
    raw_hashes = {path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(RAW.iterdir()) if path.is_file()}
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11742315",
        "article_doi": DOI, "article_published": "2024-12-01", "article_and_supplement_license": "CC-BY-4.0",
        "biostudies_file_manifest": [{
            "path": "raw/biostudies-mmc1.docx", "deposited_name": "mmc1.docx",
            "size_bytes": 1114856, "sha256": sha256(RAW / "biostudies-mmc1.docx"),
        }],
        "europe_pmc_supplement_copy": {
            "path": "raw/mmc1.docx", "size_bytes": (RAW / "mmc1.docx").stat().st_size,
            "sha256": sha256(RAW / "mmc1.docx"),
            "note": "Extracted Europe PMC package copy differs in bytes from the BioStudies deposit; both are retained and identified separately.",
        },
        "reported_unique_finite_shpepv_kcat_rows": len(rows), "accepted_records": 0,
        "exact_constructs": len(sequences), "substrate_pubchem_cids": sorted(cid for cid, _ in SUBSTRATES.values()),
        "selection_policy": "All unique finite ShPepV saturation-fit kcat rows, retained without graph digitization or value-based selection; duplicate WT table row counted once.",
        "construct_mapping": f"QKQ29470.1 is cross-referenced by archived A0A657M1L9; article pBADcLIC expression retains C-terminal {TAG}; source-numbered mutations applied before tag addition.",
        "raw_file_hashes": raw_hashes, "blocker_evidence": "blocker-evidence.json", "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")
    cluster_lines = [line for line in cluster.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust", "mmseqs_version": version,
        "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "candidate_records": 0,
        "reported_exclusion_records": len(rows), "unique_assayed_construct_sequences": len(sequences),
        "construct_sequences_sha256": sha256(fasta), "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
        "homology_hits_sha256": sha256(hits), "candidate_mmseqs_families": len({line.split("\t", 1)[0] for line in cluster_lines}),
        "family_cluster_sha256": sha256(cluster), "accepted_records": 0,
        "claim_boundary": "Blocker evidence only; no candidate labels or model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    write_outputs(args.mmseqs, args.threads)


if __name__ == "__main__":
    main()
