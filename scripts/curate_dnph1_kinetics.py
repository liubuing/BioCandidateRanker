from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC12751057"
DOI = "10.1042/BCJ20253400"
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
    "PMC12751057-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12751057/fullTextXML",
    "PMC12751057-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12751057/supplementaryFiles?includeInlineImage=false",
    "biostudies-metadata.json": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12751057",
    "PMC10483697-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10483697/fullTextXML",
    "PMC10483697-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10483697/supplementaryFiles?includeInlineImage=false",
    "O43598.fasta": "https://rest.uniprot.org/uniprotkb/O43598.fasta",
    "8QHQ.cif": "https://files.rcsb.org/download/8QHQ.cif",
    "8QHQ-entry.json": "https://data.rcsb.org/rest/v1/core/entry/8QHQ",
    "8QHQ-polymer-entity-1.json": "https://data.rcsb.org/rest/v1/core/polymer_entity/8QHQ/1",
    "8OS9.cif": "https://files.rcsb.org/download/8OS9.cif",
    "8OSC.cif": "https://files.rcsb.org/download/8OSC.cif",
    "pdb-ccd-5HU.json": "https://data.rcsb.org/rest/v1/core/chemcomp/5HU",
    "pubchem-447206.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/447206/property/IsomericSMILES,CanonicalSMILES,InChI,InChIKey,Title/JSON",
}

SUBSTRATE = {
    "name": "5-hydroxymethyl-2'-deoxyuridine 5'-monophosphate",
    "pubchem_cid": 447206,
    "isomeric_smiles": "C1[C@@H]([C@H](O[C@H]1N2C=C(C(=O)NC2=O)CO)COP(=O)(O)O)O",
    "inchi_key": "WEBVWKFGRVLCNS-XLPZGREQSA-N",
    "pdb_ccd": "5HU",
}

# Supplementary Table S2. Errors are best-fit errors; all rows are direct steady-state kcat.
KINETICS = (
    ("dnph1-wt-ph7", "WT", 7.0, 0.210, 0.001, 4),
    ("dnph1-wt-ph8p5", "WT", 8.5, 0.054, 0.002, 2),
    ("dnph1-h56a-ph7", "H56A", 7.0, 0.015, 0.001, 2),
    ("dnph1-e55q-ph7", "E55Q", 7.0, 0.0015, 0.0001, 2),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    payload = response.read()
                if not payload:
                    raise ValueError(f"Empty response for {url}")
                path.write_bytes(payload)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def mutate_native(native: str, position: int, expected: str, replacement: str) -> str:
    if native[position - 1] != expected:
        raise ValueError(f"O43598 expected {expected}{position}, found {native[position - 1]}")
    return native[: position - 1] + replacement + native[position:]


def build_sequences() -> dict[str, str]:
    native = read_fasta(RAW / "O43598.fasta")
    if len(native) != 174 or not native.startswith("MAAAMVPGRSESWERGEPGR"):
        raise ValueError("Unexpected UniProt O43598 sequence")
    variants = {
        "dnph1-wt": native,
        "dnph1-h56a": mutate_native(native, 56, "H", "A"),
        "dnph1-e55q": mutate_native(native, 55, "E", "Q"),
    }
    # MHHHHHHENLYFQG-native was purified and TEV-cleaved. Cleavage leaves one vector G.
    return {identifier: "G" + sequence for identifier, sequence in variants.items()}


def xml_text(path: Path) -> tuple[ET.Element, str]:
    root = ET.parse(path).getroot()
    return root, re.sub(r"\s+", " ", " ".join(root.itertext())).strip()


def cif_polymer_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s+\n?;([^;]+);", text, re.DOTALL
    )
    if not match:
        raise ValueError(f"Canonical polymer sequence absent from {path.name}")
    return re.sub(r"\s+", "", match.group(1))


def verify_raw(sequences: dict[str, str]) -> None:
    missing = [name for name in URLS if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw files: {', '.join(missing)}; run with --download")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP development reference hash changed")

    article, article_text = xml_text(RAW / "PMC12751057-fullText.xml")
    if article.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Unexpected measurement-article DOI")
    required_article = (
        "Creative Commons Attribution License 4.0",
        "activity under steady-state conditions",
        "1.125–40 μM 5hmdUMP",
        "Two independent measurements were performed, except for WT-",
        "Production of WT-, H56A- and E104A-",
        "variants were produced as previously described",
        "Supplementary Table S2",
    )
    for evidence in required_article:
        if evidence not in article_text:
            raise ValueError(f"Measurement evidence absent from source XML: {evidence}")

    prior, prior_text = xml_text(RAW / "PMC10483697-fullText.xml")
    if prior.findtext(".//article-id[@pub-id-type='doi']") != "10.1021/acs.biochem.3c00369":
        raise ValueError("Unexpected construct-provenance article DOI")
    for evidence in (
        "MHHHHHHENLYFQG-protein sequence",
        "TEVP-cleavable N-terminal His-tag",
        "DNA-encoding",
        "O43598",
        "19 164.9",
    ):
        if evidence not in prior_text:
            raise ValueError(f"Construct evidence absent from provenance XML: {evidence}")

    metadata = json.loads((RAW / "biostudies-metadata.json").read_text(encoding="utf-8"))
    files = metadata["section"]["files"]
    if metadata["accno"] != "S-EPMC12751057" or files != [
        {"path": "bcj-482-24-BCJ20253400-s001.pdf", "size": 1502659, "type": "file"}
    ]:
        raise ValueError("Unexpected BioStudies accession or file manifest")
    with zipfile.ZipFile(RAW / "PMC12751057-supplementaryFiles.zip") as archive:
        info = archive.getinfo("bcj-482-24-BCJ20253400-s001.pdf")
        if info.file_size != 1502659:
            raise ValueError("Unexpected DNPH1 supplementary PDF size")

    pubchem = json.loads((RAW / "pubchem-447206.json").read_text(encoding="utf-8"))
    compound = pubchem["PropertyTable"]["Properties"][0]
    if (
        compound["CID"] != SUBSTRATE["pubchem_cid"]
        or compound["SMILES"] != SUBSTRATE["isomeric_smiles"]
        or compound["InChIKey"] != SUBSTRATE["inchi_key"]
    ):
        raise ValueError("PubChem 5hmdUMP mapping changed")
    ccd = json.loads((RAW / "pdb-ccd-5HU.json").read_text(encoding="utf-8"))
    if ccd["chem_comp"]["id"] != "5HU" or not any(
        row.get("resource_name") == "PubChem" and row.get("resource_accession_code") == "447206"
        for row in ccd["rcsb_chem_comp_related"]
    ):
        raise ValueError("PDB CCD 5HU to PubChem mapping changed")

    entity = json.loads((RAW / "8QHQ-polymer-entity-1.json").read_text(encoding="utf-8"))
    pdb_sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
    if pdb_sequence != cif_polymer_sequence(RAW / "8QHQ.cif"):
        raise ValueError("8QHQ API and coordinate polymer sequences differ")
    if entity["rcsb_polymer_entity_align"][0]["reference_database_accession"] != "O43598":
        raise ValueError("8QHQ no longer maps to O43598")
    if "5-HYDROXYMETHYLURIDINE-2'-DEOXY-5'-MONOPHOSPHATE" not in (
        RAW / "8QHQ.cif"
    ).read_text(encoding="utf-8"):
        raise ValueError("8QHQ lacks bound 5hmdUMP")
    truncated = "GM" + read_fasta(RAW / "O43598.fasta")[19:162]
    for pdb_id in ("8OS9", "8OSC"):
        if cif_polymer_sequence(RAW / f"{pdb_id}.cif") != truncated:
            raise ValueError(f"{pdb_id} does not corroborate the TEV-cleaved truncated construct")
    if any(len(sequence) != 175 for sequence in sequences.values()):
        raise ValueError("Expected three 175-residue exact assay products")


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(sequences: dict[str, str]) -> Path:
    path = SOURCE / "construct_sequences.fasta"
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sequences.items():
            variant = identifier.removeprefix("dnph1-").upper()
            handle.write(
                f">{identifier} O43598 | TEV-cleaved G + full-length native sequence | {variant}\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    return path


def build_records(hit_ids: set[str]) -> list[dict[str, object]]:
    records = []
    for index, (condition_id, variant, ph, kcat, error, replicates) in enumerate(KINETICS, 1):
        sequence_id = f"dnph1-{variant.lower()}"
        status = "excluded_homology" if sequence_id in hit_ids else "accepted_homology_cold_pool"
        records.append(
            {
                "candidate_id": f"human-dnph1-{index:03d}",
                "article_doi": DOI,
                "source_table": "Supplementary Table S2",
                "source_row": f"HsDNPH1 {variant}; pH {ph:.1f}; kcat row",
                "organism": "Homo sapiens",
                "sequence_id": sequence_id,
                "uniprot_accession": "O43598",
                "variant": variant,
                "construct": "TEV-cleaved G + full-length O43598 (native residues 1-174)",
                "variable_substrate": SUBSTRATE["name"],
                "substrate_pubchem_cid": SUBSTRATE["pubchem_cid"],
                "substrate_isomeric_smiles": SUBSTRATE["isomeric_smiles"],
                "substrate_pdb_ccd": SUBSTRATE["pdb_ccd"],
                "endpoint": "kcat_s-1",
                "kcat_s-1": kcat,
                "kcat_best_fit_error_s-1": error,
                "assay_pH": ph,
                "assay_temperature_C": 25,
                "assay_buffer": "100 mM HEPES" if ph == 7.0 else "100 mM TAPS",
                "replicates": replicates,
                "status_at_normalization": status,
                "condition_id": condition_id,
            }
        )
    return records


def exact_overlap_audit(
    sequences: dict[str, str], development_corpus: Path
) -> dict[str, object]:
    if not development_corpus.is_file():
        raise FileNotFoundError(f"Frozen development corpus is absent: {development_corpus}")
    if development_corpus.stat().st_size != DEVELOPMENT_CORPUS_SIZE:
        raise ValueError("Frozen development corpus size changed")
    if sha256(development_corpus) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("Frozen development corpus SHA256 changed")
    rows = json.loads(development_corpus.read_text(encoding="utf-8"))
    if not isinstance(rows, list) or len(rows) != 17010:
        raise ValueError("Frozen development corpus row count changed")
    candidate_sequences = set(sequences.values())
    exact_sequence_rows = [
        row for row in rows if isinstance(row, dict) and row.get("Sequence") in candidate_sequences
    ]
    exact_pairs = [
        row
        for row in exact_sequence_rows
        if row.get("Smiles") == SUBSTRATE["isomeric_smiles"]
    ]
    return {
        "method": "Exact string comparison of complete assayed product sequence and isomeric SMILES against every frozen development row",
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_corpus_size_bytes": DEVELOPMENT_CORPUS_SIZE,
        "development_corpus_rows": len(rows),
        "exact_sequence_rows": len(exact_sequence_rows),
        "exact_sequence_substrate_overlap": len(exact_pairs),
    }


def mmseqs_command(executable: str, *arguments: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *arguments]
    _, distribution, command = executable.split(":", 2)
    result = ["wsl"]
    if distribution:
        result += ["-d", distribution]
    return [*result, command, *arguments]


def tool_path(path: Path, executable: str) -> str:
    if not executable.startswith("wsl:"):
        return str(path)
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}/{resolved.as_posix()[3:]}"


def run_checked(command: list[str]) -> str:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}\n{completed.stderr}")
    return (completed.stdout or completed.stderr).strip()


def run_homology(fasta: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(mmseqs_command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise RuntimeError(f"MMseqs2 version {version!r} != frozen version {MMSEQS_VERSION!r}")
    homology = SOURCE / "homology"
    if homology.exists():
        if executable.startswith("wsl:"):
            _, distribution, _ = executable.split(":", 2)
            command = ["wsl"]
            if distribution:
                command += ["-d", distribution]
            run_checked([*command, "rm", "-rf", "--", tool_path(homology, executable)])
        else:
            shutil.rmtree(homology)
    homology.mkdir(parents=True)
    hits = homology / "homology_hits.tsv"
    cluster_prefix = homology / "family-cluster" / "proteins"
    cluster_prefix.parent.mkdir()
    run_checked(
        mmseqs_command(
            executable,
            "easy-search",
            tool_path(fasta, executable),
            tool_path(REFERENCE, executable),
            tool_path(hits, executable),
            tool_path(homology / "search-tmp", executable),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.8",
            "--cov-mode",
            "0",
            "--format-output",
            "query,target,fident,alnlen,qcov,tcov,evalue,bits",
            "--threads",
            str(threads),
        )
    )
    run_checked(
        mmseqs_command(
            executable,
            "easy-linclust",
            tool_path(fasta, executable),
            tool_path(cluster_prefix, executable),
            tool_path(homology / "cluster-tmp", executable),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.8",
            "--cov-mode",
            "0",
            "--threads",
            str(threads),
        )
    )
    cluster = cluster_prefix.with_name("proteins_cluster.tsv")
    if not hits.is_file() or not cluster.is_file():
        raise RuntimeError("MMseqs2 did not emit required evidence")
    return hits, cluster, version


def write_outputs(mmseqs: str, threads: int, development_corpus: Path) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    sequences = build_sequences()
    verify_raw(sequences)
    overlap = exact_overlap_audit(sequences, development_corpus)
    fasta = write_fasta(sequences)
    hits, cluster, version = run_homology(fasta, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="utf-8").splitlines() if line]
    hit_ids = {row[0] for row in hit_rows}
    records = build_records(hit_ids)
    write_csv(SOURCE / "candidate_records.csv", records)

    mapping_rows = []
    native = read_fasta(RAW / "O43598.fasta")
    for sequence_id, sequence in sequences.items():
        variant = sequence_id.removeprefix("dnph1-").upper()
        mutation = {"WT": "none", "H56A": "O43598 H56A", "E55Q": "O43598 E55Q"}[variant]
        mapping_rows.append(
            {
                "sequence_id": sequence_id,
                "variant": variant,
                "uniprot_accession": "O43598",
                "native_length": len(native),
                "assayed_product_length": len(sequence),
                "native_mapping": "assay residues 2-175 = O43598 residues 1-174",
                "vector_residues_retained_after_tev": "G at assay position 1",
                "mutation_mapping": mutation,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
        )
    write_csv(SOURCE / "sequence_mapping.csv", mapping_rows)

    structure_rows = [
        {
            "pdb_id": "8QHQ",
            "role": "2025 article substrate-bound QM/MM starting structure; experimental E104Q construct",
            "ligand": "5HU (5hmdUMP)",
            "kinetic_construct_exact_match": False,
            "mapping": "SIFTS O43598 residues 19-162; GSM tag; E104Q",
            "source_article_doi": "10.1038/s41467-023-42544-4",
            "coordinate_sha256": sha256(RAW / "8QHQ.cif"),
        },
        {
            "pdb_id": "8OS9",
            "role": "construct-provenance unliganded truncated WT structure",
            "ligand": "none",
            "kinetic_construct_exact_match": False,
            "mapping": "GM + O43598 residues 20-162",
            "source_article_doi": "10.1021/acs.biochem.3c00369",
            "coordinate_sha256": sha256(RAW / "8OS9.cif"),
        },
        {
            "pdb_id": "8OSC",
            "role": "construct-provenance dUMP-bound truncated WT structure",
            "ligand": "UMP (dUMP)",
            "kinetic_construct_exact_match": False,
            "mapping": "GM + O43598 residues 20-162",
            "source_article_doi": "10.1021/acs.biochem.3c00369",
            "coordinate_sha256": sha256(RAW / "8OSC.cif"),
        },
    ]
    write_csv(SOURCE / "structure_mapping.csv", structure_rows)

    raw_hashes = {
        name: {
            "url": URLS[name],
            "size_bytes": (RAW / name).stat().st_size,
            "sha256": sha256(RAW / name),
        }
        for name in URLS
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )

    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12751057",
        "article_doi": DOI,
        "pmc_id": "PMC12751057",
        "article_published": "2025-12-17",
        "article_and_supplement_license": "CC-BY-4.0",
        "construct_provenance_article": "PMC10483697; DOI 10.1021/acs.biochem.3c00369; CC-BY-4.0",
        "reported_direct_steady_state_kcat_rows_curated": len(records),
        "accepted_records": len(accepted),
        "excluded_homology_records": len(records) - len(accepted),
        "kinetics_source": "raw/PMC12751057-supplementaryFiles.zip, bcj-482-24-BCJ20253400-s001.pdf, Supplementary Table S2; assay conditions in raw/PMC12751057-fullText.xml",
        "selection_policy": "All direct finite steady-state kcat rows in Supplementary Table S2 for 5hmdUMP at 25 C; no value-based selection. Binding, burst, single-turnover, computed QM/MM, inhibitor, and kcat/KM-only results are not candidate kcat labels.",
        "sequence_source": "raw/O43598.fasta plus the explicit MHHHHHHENLYFQG-protein and TEV-cleavage provenance in raw/PMC10483697-fullText.xml; source-native mutation numbering",
        "construct_audit": "The kinetic proteins are full-length O43598. TEV cleavage of MHHHHHHENLYFQG-native leaves G-native, yielding 175 residues. WT, H56A and E55Q products are reconstructed exactly; mutations are applied at O43598 native coordinates before adding G.",
        "substrate_mapping": "5hmdUMP = PubChem CID 447206 = PDB CCD 5HU, InChIKey WEBVWKFGRVLCNS-XLPZGREQSA-N; raw PubChem and RCSB CCD records retained.",
        "structure_boundary": "No deposited structure is claimed as an exact kinetic construct. 8QHQ is a GSM-tagged E104Q O43598(19-162) substrate complex from another article and is only corroborating ligand/structure evidence. 8OS9/8OSC are GM-tagged O43598(20-162) constructs corroborating the TEV residue convention.",
        "raw_file_hashes": "raw-file-hashes.json",
        "exact_overlap_audit": overlap,
        "artifact_sha256": {
            name: sha256(SOURCE / name)
            for name in (
                "candidate_records.csv",
                "construct_sequences.fasta",
                "sequence_mapping.csv",
                "structure_mapping.csv",
            )
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    cluster_rows = [line.split("\t") for line in cluster.read_text(encoding="utf-8").splitlines() if line]
    audit = {
        "audited_on": str(date.today()),
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": version,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "easy_search_format": "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "search_command": "mmseqs easy-search construct_sequences.fasta unikp_reference.fasta homology_hits.tsv search-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "cluster_command": "mmseqs easy-linclust construct_sequences.fasta proteins cluster-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0",
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": len(records),
        "unique_assayed_construct_sequences": len(sequences),
        "construct_sequences_sha256": sha256(fasta),
        "homology_hit_sequences": len(hit_ids),
        "exact_sequence_overlap": sum(float(row[2]) >= 0.999 for row in hit_rows),
        "homology_hits_sha256": sha256(hits),
        "candidate_mmseqs_families": len({row[0] for row in cluster_rows}),
        "family_cluster_sha256": sha256(cluster),
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"],
        "exact_overlap_development_corpus_sha256": overlap["development_corpus_sha256"],
        "accepted_records": len(accepted),
        "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
        "accepted_unique_substrates": len({row["substrate_pubchem_cid"] for row in accepted}),
        "readiness_gate_passes": False,
        "claim_boundary": "Source-level curation pool only; global readiness gates are not evaluated here and no model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )
    verify_outputs(records, sequences)


def verify_outputs(records: list[dict[str, object]], sequences: dict[str, str]) -> None:
    if len(records) != 4 or len(sequences) != 3 or len(set(sequences.values())) != 3:
        raise ValueError("Unexpected DNPH1 curation cardinality")
    if any(row["endpoint"] != "kcat_s-1" or float(row["kcat_s-1"]) <= 0 for row in records):
        raise ValueError("Invalid candidate endpoint")
    if any(row["substrate_pubchem_cid"] != 447206 for row in records):
        raise ValueError("Unexpected candidate substrate mapping")
    audit = json.loads((SOURCE / "homology-audit.json").read_text(encoding="ascii"))
    if audit["mmseqs_version"] != MMSEQS_VERSION or audit["coverage_mode"] != 0:
        raise ValueError("Homology audit violates frozen protocol")
    provenance = json.loads((SOURCE / "provenance.json").read_text(encoding="ascii"))
    if provenance["model_predictions_run"]:
        raise ValueError("No-peeking boundary violated")


def write_blocker(exc: Exception) -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    evidence = {
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "status": "blocked_fail_closed",
        "error_type": type(exc).__name__,
        "error": str(exc),
        "recorded_on": str(date.today()),
        "candidate_records_released": False,
        "model_predictions_run": False,
        "raw_files_present": sorted(path.name for path in RAW.glob("*") if path.is_file()),
    }
    (SOURCE / "blocker-evidence.json").write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="ascii"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--development-corpus", type=Path, default=DEFAULT_DEVELOPMENT_CORPUS)
    args = parser.parse_args()
    try:
        if args.download:
            download_raw()
        write_outputs(args.mmseqs, args.threads, args.development_corpus)
        blocker = SOURCE / "blocker-evidence.json"
        if blocker.exists():
            blocker.unlink()
    except Exception as exc:
        write_blocker(exc)
        raise


if __name__ == "__main__":
    main()
