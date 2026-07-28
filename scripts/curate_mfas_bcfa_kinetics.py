from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
import zipfile
from datetime import date
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC12267109"
DOI = "10.1002/pro.70229"
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
    r"C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment"
    r"\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json"
)

URLS = {
    "PMC12267109-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12267109/fullTextXML",
    "biostudies-metadata.json": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12267109",
    "PRO-34-e70229-s001.zip": "https://www.ebi.ac.uk/biostudies/files/S-EPMC12267109/PRO-34-e70229-s001.zip",
    "P19096.fasta": "https://rest.uniprot.org/uniprotkb/P19096.fasta",
    "P19096.json": "https://rest.uniprot.org/uniprotkb/P19096.json",
    "6ROP.fasta": "https://www.rcsb.org/fasta/entry/6ROP/display",
    "pubchem-acetyl-CoA-444493.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/444493/property/IsomericSMILES,Title/JSON",
    "pubchem-methylmalonyl-CoA-123909.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/123909/property/IsomericSMILES,Title/JSON",
    "pubchem-NADPH-5884.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5884/property/IsomericSMILES,Title/JSON",
}

COMPOUNDS = {
    "acetyl-CoA": (444493, "CC(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O"),
    "methylmalonyl-CoA": (123909, "CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O"),
    "NADPH": (5884, "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    with zipfile.ZipFile(RAW / "PRO-34-e70229-s001.zip") as archive:
        member = next(name for name in archive.namelist() if name.lower().endswith(".pdf"))
        (RAW / "si_gusenda-ochs.pdf").write_bytes(archive.read(member))


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


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


def verify_raw() -> tuple[str, str, str]:
    required_files = (*URLS, "si_gusenda-ochs.pdf")
    missing = [name for name in required_files if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw evidence: {missing}; run with --download")
    xml = (RAW / "PMC12267109-fullText.xml").read_text(encoding="utf-8")
    article_text = re.sub(r"<[^>]+>", "", xml)
    article_text = re.sub(r"\s+", " ", article_text).replace("‐", "-").replace("–", "-")
    required_text = (
        DOI,
        "10-1000 μM acetyl-CoA",
        "100 μM metmal-CoA",
        "50 μM NADPH",
        "270 μM metmal-ACP",
        "5 μM MabA",
        "KS-MAT",
        "S581A",
    )
    for token in required_text:
        if token not in article_text:
            raise ValueError(f"Required article evidence missing: {token}")
    supplement_text = "\n".join(
        page.get_text() for page in fitz.open(RAW / "si_gusenda-ochs.pdf")
    )
    for token in ("substrate inhibition with equation 3", "cooperative kinetics and substrate", "inhibition with equation 4"):
        if token not in supplement_text:
            raise ValueError(f"Required supplement evidence missing: {token}")

    native = read_fasta(RAW / "P19096.fasta")
    if len(native) != 2504 or native[580] != "S":
        raise ValueError("Unexpected mouse FASN P19096 sequence")
    pdb_sequence = read_fasta(RAW / "6ROP.fasta")
    if len(pdb_sequence) != 852 or pdb_sequence[1:] != native[1:852]:
        raise ValueError("6ROP does not verify the native-coordinate P19096(1-852) KS-MAT boundary")
    ks_mat0 = native[:580] + "A" + native[581:852]

    for name, (expected_cid, expected_smiles) in COMPOUNDS.items():
        raw_name = {
            "acetyl-CoA": "pubchem-acetyl-CoA-444493.json",
            "methylmalonyl-CoA": "pubchem-methylmalonyl-CoA-123909.json",
            "NADPH": "pubchem-NADPH-5884.json",
        }[name]
        item = json.loads((RAW / raw_name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if item["CID"] != expected_cid or item.get("SMILES", item.get("IsomericSMILES")) != expected_smiles:
            raise ValueError(f"Unexpected PubChem mapping for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference changed")
    return native, ks_mat0, supplement_text


def exact_overlap(sequences: dict[str, str], development_corpus: Path) -> dict[str, object]:
    if not development_corpus.is_file() or development_corpus.stat().st_size != DEVELOPMENT_CORPUS_SIZE:
        raise FileNotFoundError("Frozen 17,010-row development corpus is absent or changed")
    if sha256(development_corpus) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("Frozen development corpus hash changed")
    rows = json.loads(development_corpus.read_text(encoding="utf-8"))
    if not isinstance(rows, list) or len(rows) != 17010:
        raise ValueError("Frozen development row count changed")
    sequence_set = set(sequences.values())
    smiles = {value[1] for value in COMPOUNDS.values()}
    exact_sequence_rows = [row for row in rows if row.get("Sequence") in sequence_set]
    exact_pairs = [row for row in exact_sequence_rows if row.get("Smiles") in smiles]
    return {
        "method": "Exact string comparison of both native-coordinate catalytic queries and three PubChem isomeric SMILES against all frozen development rows",
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_corpus_size_bytes": DEVELOPMENT_CORPUS_SIZE,
        "development_corpus_rows": len(rows),
        "exact_sequence_rows": len(exact_sequence_rows),
        "exact_sequence_substrate_overlap": len(exact_pairs),
    }


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in sequences.items():
            handle.write(f">{sequence_id}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


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
    common = {
        "article_doi": DOI,
        "organism": "Mus musculus",
        "uniprot_accession": "P19096",
        "status_at_normalization": "excluded_fixed_component_saturation_unresolved",
        "candidate_label_created": False,
    }
    return [
        {
            **common,
            "evidence_id": "mfas-bcfa-table1",
            "endpoint_family": "integrated_metazoan_type_I_FAS_fatty_acid_cycle",
            "source_table": "Table 1, fatty acid cycle (total mFAS)",
            "sequence_id": "P19096-native-full-length-query",
            "construct": "full-length mFAS; Ni-NTA and Strep-Tactin purification imply tags whose sequences and positions are not disclosed",
            "varied_substrate": "acetyl-CoA",
            "varied_substrate_pubchem_cid": 444493,
            "varied_substrate_concentration_range_uM": "10-1000",
            "fit_semantics": "apparent substrate-inhibition fit; acetyl-CoA and methylmalonyl-CoA likely compete at MAT",
            "reported_apparent_kcat_s-1": 0.0176,
            "reported_apparent_kcat_se_s-1": 0.0014,
            "reported_apparent_half_saturation_uM": 69,
            "reported_apparent_half_saturation_se_uM": 12,
            "reported_inhibition_constant_mM": 2.9,
            "reported_inhibition_constant_se_mM": 1.2,
            "fixed_methylmalonyl_CoA_uM": 100,
            "fixed_methylmalonyl_CoA_pubchem_cid": 123909,
            "fixed_methylmalonyl_CoA_saturation": "not established; source states self-priming overlap makes accurate extender kinetic characterization unfeasible",
            "fixed_NADPH_uM": 50,
            "fixed_NADPH_pubchem_cid": 5884,
            "fixed_NADPH_saturation": "not demonstrated",
            "exclusion_reason": "fixed methylmalonyl-CoA and NADPH saturation not demonstrated; exact expressed tag sequence unresolved",
        },
        {
            **common,
            "evidence_id": "ks-mat0-bcfa-table1",
            "endpoint_family": "isolated_KS_domain_activity_in_KS_MAT0_didomain_coupled_assay",
            "source_table": "Table 1, chain elongation (KS single domain)",
            "sequence_id": "P19096-KS-MAT-residues-1-852-S581A-query",
            "construct": "KS-MAT S581A didomain used because isolated KS cannot be expressed; exact termini and affinity-tag sequence are not disclosed",
            "varied_substrate": "decanoyl-ACP",
            "varied_substrate_pubchem_cid": "",
            "varied_substrate_concentration_range_uM": "0-210",
            "fit_semantics": "apparent combined cooperativity-and-substrate-inhibition fit",
            "reported_apparent_kcat_s-1": 0.0116,
            "reported_apparent_kcat_se_s-1": 0.0006,
            "reported_apparent_half_saturation_uM": 23.7,
            "reported_apparent_half_saturation_se_uM": 1.4,
            "reported_hill_coefficient": 3.3,
            "reported_hill_coefficient_se": 0.5,
            "reported_inhibition_constant_uM": 270,
            "reported_inhibition_constant_se_uM": 20,
            "fixed_methylmalonyl_ACP_uM": 270,
            "fixed_methylmalonyl_ACP_saturation": "not demonstrated",
            "fixed_NADPH_uM": 50,
            "fixed_NADPH_pubchem_cid": 5884,
            "fixed_NADPH_saturation": "not demonstrated",
            "coupling_MabA_uM": 5,
            "coupling_capacity": "not demonstrated non-rate-limiting across 5 uM KS and the fitted decanoyl-ACP range",
            "exclusion_reason": "fixed metmal-ACP/NADPH saturation and MabA coupling capacity not demonstrated; exact ACP, MabA, termini, and tag sequences unresolved; varied substrate is an ACP conjugate without a single PubChem mapping",
        },
    ]


def write_outputs(mmseqs: str, threads: int, development_corpus: Path) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    native, ks_mat0, _ = verify_raw()
    sequences = {
        "P19096-native-full-length-query": native,
        "P19096-KS-MAT-residues-1-852-S581A-query": ks_mat0,
    }
    overlap = exact_overlap(sequences, development_corpus)
    fasta = SOURCE / "construct_sequences.fasta"
    write_fasta(fasta, sequences)
    hits, clusters, version = run_homology(fasta, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="utf-8").splitlines() if line]
    cluster_rows = [line.split("\t") for line in clusters.read_text(encoding="utf-8").splitlines() if line]

    rows = excluded_records()
    write_csv(SOURCE / "excluded_records.csv", rows)
    (SOURCE / "candidate_records.csv").write_text("", encoding="ascii")
    write_csv(SOURCE / "sequence-mapping.csv", [
        {"sequence_id": "P19096-native-full-length-query", "audited_sequence": "exact UniProt P19096 full native sequence", "length": len(native), "expressed_product_exact": False, "blocker": "article does not disclose affinity-tag sequences or positions"},
        {"sequence_id": "P19096-KS-MAT-residues-1-852-S581A-query", "audited_sequence": "P19096 residues 1-852 with source-stated S581A", "length": len(ks_mat0), "expressed_product_exact": False, "blocker": "article does not disclose exact termini/tag; 6ROP independently verifies native-coordinate residues 2-852"},
    ])
    write_csv(SOURCE / "saturation-audit.csv", [
        {"evidence_id": rows[0]["evidence_id"], "component": "acetyl-CoA", "role": "varied", "concentration": "10-1000 uM", "saturation_status": "substrate-inhibition model; apparent kcat is model asymptote, not rate at highest concentration"},
        {"evidence_id": rows[0]["evidence_id"], "component": "racemic methylmalonyl-CoA", "role": "fixed extender", "concentration": "100 uM", "saturation_status": "not established; self-priming prevents accurate extender characterization"},
        {"evidence_id": rows[0]["evidence_id"], "component": "NADPH", "role": "fixed reductant/reporting cofactor", "concentration": "50 uM", "saturation_status": "not demonstrated"},
        {"evidence_id": rows[1]["evidence_id"], "component": "decanoyl-ACP", "role": "varied starter", "concentration": "0-210 uM", "saturation_status": "combined cooperativity/substrate-inhibition model; apparent kcat is model asymptote"},
        {"evidence_id": rows[1]["evidence_id"], "component": "methylmalonyl-ACP", "role": "fixed extender", "concentration": "270 uM", "saturation_status": "not demonstrated"},
        {"evidence_id": rows[1]["evidence_id"], "component": "NADPH", "role": "fixed coupling cofactor", "concentration": "50 uM", "saturation_status": "not demonstrated"},
        {"evidence_id": rows[1]["evidence_id"], "component": "MabA", "role": "coupling enzyme", "concentration": "5 uM", "saturation_status": "non-rate-limiting coupling capacity not demonstrated for 5 uM KS"},
    ])
    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path), "url_or_derivation": URLS.get(path.name, "extracted from raw/PRO-34-e70229-s001.zip")}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    blocker = {
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "status": "excluded_familywise_fail_closed",
        "candidate_records_released": False,
        "direct_apparent_kcat_evidence_rows_retained": 2,
        "integrated_mfas_disposition": rows[0]["exclusion_reason"],
        "isolated_domain_disposition": rows[1]["exclusion_reason"],
        "no_label_inference": "The two reported apparent kcat values are retained only in excluded_records.csv and are not duplicated by each reaction component.",
        "model_predictions_run": False,
        "recorded_on": str(date.today()),
    }
    (SOURCE / "blocker-evidence.json").write_text(json.dumps(blocker, indent=2) + "\n", encoding="ascii")
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "MMseqs2 easy-search and easy-linclust of both native-coordinate catalytic queries before admission", "mmseqs_version": version,
        "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "construct_sequences_sha256": sha256(fasta),
        "candidate_records": 0, "accepted_records": 0,
        "queried_catalytic_sequences": len(sequences), "query_scope_limitation": "Exact expressed affinity-tagged products were not recoverable; no exact-construct admission was attempted.",
        "homology_hit_sequences": len({row[0] for row in hit_rows}), "homology_hits": len(hit_rows), "homology_hits_sha256": sha256(hits),
        "query_disposition": {
            sequence_id: {
                "qualifying_hits": sum(row[0] == sequence_id for row in hit_rows),
                "frozen_homology_gate": "failed" if any(row[0] == sequence_id for row in hit_rows) else "passed",
            }
            for sequence_id in sequences
        },
        "candidate_mmseqs_families": len({row[0] for row in cluster_rows}), "family_cluster_sha256": sha256(clusters),
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"], "exact_overlap_development_corpus_sha256": overlap["development_corpus_sha256"],
        "homology_disposition": "Saturation, defined-substrate, coupling-capacity, and exact-construct gates already exclude both families; MMseqs evidence is retained and does not override those blockers.",
        "model_predictions_run": False,
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12267109", "article_doi": DOI,
        "article_published": "2025-07-16", "article_and_supplement_license": "CC-BY-4.0",
        "kinetics_sources": ["raw/PMC12267109-fullText.xml, Table 1 and Methods 4.5-4.6", "raw/si_gusenda-ochs.pdf, Equations S3-S4 and Figures S1/S5"],
        "sequence_sources": ["raw/P19096.fasta and P19096.json", "raw/6ROP.fasta native-coordinate KS-MAT evidence"],
        "compound_mappings": {name: {"pubchem_cid": cid, "isomeric_smiles": smiles} for name, (cid, smiles) in COMPOUNDS.items()},
        "endpoint_separation": {"integrated_metazoan_type_I_FAS": "one apparent reaction-level kcat from acetyl-CoA substrate-inhibition fit", "isolated_domain": "one apparent KS elongation kcat measured with KS-MAT S581A and MabA coupling; not an integrated FAS endpoint"},
        "reported_direct_apparent_kcat_rows": 2, "accepted_records": 0,
        "selection_policy": "Audit both and only the two new Table 1 apparent kcat endpoints, family by family, without value selection; exclude each independently when any frozen gate fails.",
        "exact_overlap_audit": overlap, "saturation_audit": "saturation-audit.csv", "blocker_evidence": blocker,
        "raw_file_hashes": "raw-file-hashes.json", "artifact_sha256": {name: sha256(SOURCE / name) for name in ("excluded_records.csv", "candidate_records.csv", "construct_sequences.fasta", "sequence-mapping.csv", "saturation-audit.csv", "blocker-evidence.json")},
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")


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
