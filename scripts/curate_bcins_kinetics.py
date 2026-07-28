from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12958101"
RAW = SOURCE / "raw"
REFERENCE = (
    ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723"
    / "homology" / "unikp_reference.fasta"
)
DOI = "10.1111/febs.70308"
PMC_ID = "PMC12958101"
UNIPROT_ACCESSION = "B5GMG2"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
DEFAULT_MMSEQS = "mmseqs"

DOWNLOADS = {
    "PMC12958101-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12958101/fullTextXML"
    ),
    "PMC12958101-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12958101/"
        "supplementaryFiles?includeInlineImage=false"
    ),
    "B5GMG2.fasta": "https://rest.uniprot.org/uniprotkb/B5GMG2.fasta",
    "B5GMG2.json": "https://rest.uniprot.org/uniprotkb/B5GMG2.json",
    "pubchem-GPP-445995.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/445995/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "pubchem-FPP-445713.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/445713/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "pubchem-GGPP-447277.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/447277/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "rcsb-5NX6.json": "https://data.rcsb.org/rest/v1/core/entry/5NX6",
    "rcsb-5NX7.json": "https://data.rcsb.org/rest/v1/core/entry/5NX7",
}

SUBSTRATES = {
    "GPP": ("geranyl diphosphate", 445995),
    "FPP": ("farnesyl diphosphate", 445713),
    "GGPP": ("geranylgeranyl diphosphate", 447277),
}

# Main text, "In vitro kinetic assays reveal substrate-specific activities". These are
# all finite apparent values reported for the purified-enzyme experiment. Product-branch
# values are kept distinct from total detected-product turnover.
REPORTED_VALUES = [
    ("WT", "GPP", "cineole", 4.14, 0.68),
    ("WT", "GPP", "total_detected_products", 4.38, 0.70),
    ("F74A", "GPP", "cineole", 1.24, 0.21),
    ("F74A", "GPP", "total_detected_products", 1.94, 0.30),
    ("F74A", "FPP", "sesquicineole", 2.68, 0.37),
    ("F74A", "FPP", "total_detected_products", 3.79, 0.63),
    ("F74A-F179A", "GPP", "cineole", 3.03, 0.48),
    ("F74A-F179A", "GPP", "total_detected_products", 5.07, 0.86),
    ("F74A-F179A", "FPP", "sesquicineole", 5.04, 0.07),
    ("F74A-F179A", "FPP", "total_detected_products", 5.49, 0.10),
    ("W58A-F74A-F179A", "GPP", "linalool", 0.24, 0.04),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def acquire() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in DOWNLOADS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size > 100:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/0.1"})
        with urllib.request.urlopen(request, timeout=180) as response:
            path.write_bytes(response.read())
    supplement = RAW / "FEBS-293-1341-s001.pdf"
    if not supplement.is_file():
        shutil.unpack_archive(RAW / "PMC12958101-supplementaryFiles.zip", RAW)


def read_fasta(path: Path) -> str:
    sequence = "".join(
        line.strip() for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(f"Invalid protein FASTA: {path}")
    return sequence


def mutate(sequence: str, changes: tuple[tuple[int, str, str], ...]) -> str:
    result = sequence
    for position, expected, replacement in changes:
        if result[position - 1] != expected:
            raise ValueError(f"Expected {expected}{position}, found {result[position - 1]}{position}")
        result = result[:position - 1] + replacement + result[position:]
    return result


def verify_and_build_sequences() -> dict[str, str]:
    required = [*DOWNLOADS, "FEBS-293-1341-s001.pdf"]
    missing = [name for name in required if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {missing}; run with --download")

    root = ET.parse(RAW / "PMC12958101-fullText.xml").getroot()
    if root.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Article DOI mismatch")
    license_url = root.find(".//license/{http://www.niso.org/schemas/ali/1.0/}license_ref")
    if license_url is None or "creativecommons.org/licenses/by/4.0" not in (license_url.text or ""):
        raise ValueError("Article package is not explicitly CC BY 4.0")

    native = read_fasta(RAW / "B5GMG2.fasta")
    if len(native) != 330:
        raise ValueError(f"Unexpected B5GMG2 length: {len(native)}")
    uniprot = json.loads((RAW / "B5GMG2.json").read_text(encoding="utf-8"))
    if uniprot["primaryAccession"] != UNIPROT_ACCESSION:
        raise ValueError("UniProt accession mismatch")

    sequences = {
        "bcins-native-wt": native,
        "bcins-native-f74a": mutate(native, ((74, "F", "A"),)),
        "bcins-native-f74a-f179a": mutate(native, ((74, "F", "A"), (179, "F", "A"))),
        "bcins-native-w58a-f74a-f179a": mutate(
            native, ((58, "W", "A"), (74, "F", "A"), (179, "F", "A"))
        ),
    }
    for code, (_, cid) in SUBSTRATES.items():
        payload = json.loads((RAW / f"pubchem-{code}-{cid}.json").read_text(encoding="utf-8"))
        if payload["PropertyTable"]["Properties"][0]["CID"] != cid:
            raise ValueError(f"PubChem CID mismatch for {code}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference changed")
    return sequences


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in sequences.items():
            handle.write(
                f">{sequence_id} B5GMG2 native-coordinate sequence; not the unresolved His6-TEV assay construct\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


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


def run_homology(fasta: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise RuntimeError(f"MMseqs2 version {version!r} != frozen version {MMSEQS_VERSION!r}")
    workspace = SOURCE / "homology"
    if workspace.exists():
        if executable.startswith("wsl:"):
            run_checked([
                "wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--",
                tool_path(workspace, executable),
            ])
        else:
            shutil.rmtree(workspace)
    workspace.mkdir(parents=True)
    hits = workspace / "homology_hits.tsv"
    prefix = workspace / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    run_checked(command(
        executable, "easy-search", tool_path(fasta, executable), tool_path(REFERENCE, executable),
        tool_path(hits, executable), tool_path(workspace / "search-tmp", executable),
        "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--format-output",
        "query,target,fident,alnlen,qcov,tcov,evalue,bits", "--threads", str(threads),
    ))
    run_checked(command(
        executable, "easy-linclust", tool_path(fasta, executable), tool_path(prefix, executable),
        tool_path(workspace / "cluster-tmp", executable), "--min-seq-id", "0.3", "-c", "0.8",
        "--cov-mode", "0", "--threads", str(threads),
    ))
    return hits, prefix.with_name("proteins_cluster.tsv"), version


def sequence_id(variant: str) -> str:
    return "bcins-native-" + variant.lower().replace("-", "-")


def substrate_structure(code: str) -> tuple[int, str]:
    cid = SUBSTRATES[code][1]
    item = json.loads((RAW / f"pubchem-{code}-{cid}.json").read_text(encoding="utf-8"))[
        "PropertyTable"
    ]["Properties"][0]
    return cid, item["SMILES"]


def excluded_records(hit_queries: set[str]) -> list[dict[str, object]]:
    rows = []
    for index, (variant, substrate, product_scope, value, error) in enumerate(REPORTED_VALUES, 1):
        cid, smiles = substrate_structure(substrate)
        seq_id = sequence_id(variant)
        product_specific = product_scope != "total_detected_products"
        construct_name = "pET-24d-His6-TEV-bCinS"
        if variant != "WT":
            construct_name += f"-{variant}"
        rows.append({
            "evidence_id": f"bcins-apparent-{index:02d}",
            "article_doi": DOI,
            "stable_record_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12958101/",
            "source_file": "raw/PMC12958101-fullText.xml",
            "source_section": "In vitro kinetic assays reveal substrate-specific activities",
            "source_table_crosscheck": "Supporting Information Tables S9-S10, page 26",
            "organism": "Streptomyces clavuligerus",
            "enzyme": "bacterial 1,8-cineole synthase (bCinS/CnsA)",
            "uniprot_accession": UNIPROT_ACCESSION,
            "variant": variant,
            "native_coordinate_sequence_id": seq_id,
            "assay_construct": construct_name,
            "exact_assayed_sequence_status": "unresolved_full_fusion_and_cleavage_state",
            "substrate": SUBSTRATES[substrate][0],
            "substrate_code": substrate,
            "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles,
            "reported_endpoint": "apparent_kcat_product_branch" if product_specific else "apparent_kcat_total_detected_products",
            "product_scope": product_scope,
            "reported_apparent_kcat_min-1": value,
            "reported_apparent_kcat_sd_min-1": error,
            "reported_apparent_kcat_s-1": round(value / 60, 9),
            "reported_apparent_kcat_sd_s-1": round(error / 60, 9),
            "replicates": 2,
            "assay_time_min": 5,
            "assay_temperature_C": 37,
            "assay_pH": 8.0,
            "assay_buffer": "25 mM Tris, 150 mM NaCl, 1 mM DTT, 4 mM MgCl2, 5% glycerol",
            "enzyme_concentration_uM": "1.175" if variant == "WT" else "1.0",
            "substrate_concentration_reported_range_mM": "0.25-1; row-specific concentration undisclosed",
            "substrate_saturation_status": "called excess/Vmax but no curve, Km, or row-specific concentration disclosed",
            "metal_saturation_status": "4 mM MgCl2 fixed; no metal titration or saturation demonstration",
            "rate_basis": "5-minute captured product endpoint within a separately established 120-minute linear phase",
            "homology_query_basis": "native-coordinate sequence only; conservative prescreen, not exact assay construct",
            "frozen_mmseqs_hit": seq_id in hit_queries,
            "status_at_normalization": "excluded_multiple_protocol_gates",
            "exclusion_reasons": (
                "apparent_endpoint_substrate_saturation_unverified;"
                "exact_his6_tev_assay_construct_sequence_unresolved;"
                + ("development_homology_hit" if seq_id in hit_queries else "")
            ).rstrip(";"),
            "candidate_label_created": False,
        })
    return rows


def write_outputs(sequences: dict[str, str], hits: Path, clusters: Path, version: str) -> None:
    fasta = SOURCE / "native_coordinate_homology_queries.fasta"
    hit_lines = [line for line in hits.read_text(encoding="utf-8").splitlines() if line]
    hit_queries = {line.split("\t", 1)[0] for line in hit_lines}
    rows = excluded_records(hit_queries)
    candidate_path = SOURCE / "candidate_records.csv"
    with candidate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    shutil.copyfile(candidate_path, SOURCE / "excluded_records.csv")

    qualitative = {
        "GGPP": "No diterpenoid was detected under kinetic conditions; GGPP solubility was problematic; no kinetic value exists.",
        "WT_FPP_GGPP": "No in vitro activity detected; absence observations are not finite kcat records.",
        "triple_FPP": "Trace (E)-beta-farnesene was unquantifiable; no kinetic value exists.",
        "in_vivo_tables": "Titres are culture/organic-overlay outcomes with unknown substrate supply, not direct absolute kcat.",
        "docking_and_CASTp": "Computational predictions are structure context only and never candidate labels.",
    }
    raw_hashes = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    provenance = {
        "schema_version": 1,
        "source_id": "europepmc-PMC12958101",
        "article_doi": DOI,
        "pmc_id": PMC_ID,
        "article_published_online": "2025-10-28",
        "article_license": "CC-BY-4.0",
        "status": "finalized_exclusion_only",
        "accepted_records": 0,
        "numeric_exclusion_records": len(rows),
        "native_orf": {"accession": UNIPROT_ACCESSION, "length": 330, "variants": list(sequences)},
        "assay_construct_evidence": "Supporting Information Table S2: pET-24d-His6-TEV-bCinS variants",
        "construct_blocker": (
            "Methods names a pETM-11 backbone while Supporting Information Table S2 names "
            "pET-24d-His6-TEV constructs. Neither TEV cleavage nor the complete retained "
            "His6-TEV fusion/scar sequence is disclosed."
        ),
        "endpoint_interpretation": (
            "Reported apparent kcat values are 5-minute product accumulation divided by enzyme and time. "
            "Table S9/S10 amounts numerically reproduce kcat*[E]*5; no saturation curve or Km is supplied."
        ),
        "substrate_saturation": "Unverified: only a global 0.25-1 mM range and the word excess are reported.",
        "metal_saturation": "Unverified: buffer contains fixed 4 mM MgCl2 without a metal titration.",
        "qualitative_exclusions": qualitative,
        "structures": {
            "5NX6": "1.63 A WT bCinS homodimer with 2-fluoroneryl diphosphate and Mg2+; prior DOI 10.1021/acscatal.7b01924",
            "5NX7": "1.51 A WT bCinS homodimer with 2-fluoroneryl/2-fluorogeranyl diphosphate and Mg2+; prior DOI 10.1021/acscatal.7b01924",
            "variant_structures": "None reported; variants used in silico mutation/docking and CASTp only.",
        },
        "selection_policy": "Retain every finite apparent kcat stated in the article as an exclusion; no value-based selection.",
        "raw_files": raw_hashes,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    cluster_lines = [line for line in clusters.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "source_id": "europepmc-PMC12958101",
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust of native-coordinate ORF variants",
        "prescreen_limitation": "Exact His6-TEV assay constructs are unresolved; native-coordinate searches cannot cure the construct gate.",
        "mmseqs_version": version,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target": str(REFERENCE.relative_to(ROOT)).replace("\\", "/"),
        "development_target_sha256": sha256(REFERENCE),
        "homology_queries_sha256": sha256(fasta),
        "candidate_records": {
            "count": len(rows),
            "sha256": sha256(candidate_path),
        },
        "query_sequences": len(sequences),
        "hit_queries": sorted(hit_queries),
        "development_hits": len(hit_lines),
        "homology_hits_sha256": sha256(hits),
        "candidate_families": len({line.split("\t", 1)[0] for line in cluster_lines}),
        "family_cluster_sha256": sha256(clusters),
        "accepted_records": 0,
        "claim_boundary": "Exclusion evidence only; no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.download:
        acquire()
    sequences = verify_and_build_sequences()
    fasta = SOURCE / "native_coordinate_homology_queries.fasta"
    write_fasta(fasta, sequences)
    hits, clusters, version = run_homology(fasta, args.mmseqs, args.threads)
    write_outputs(sequences, hits, clusters, version)


if __name__ == "__main__":
    main()
