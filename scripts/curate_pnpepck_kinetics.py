from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC12529880"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1002/pro.70326"
ACCESSION = "WP_011799537.1"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
PDB_IDS = ("9E32", "9E33", "9E34", "9E35", "9E36", "9E37", "9E38")

URLS = {
    "PMC12529880-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12529880/fullTextXML"
    ),
    "PMC12529880-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12529880/supplementaryFiles"
    ),
    f"{ACCESSION}.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=protein&id={ACCESSION}&rettype=fasta&retmode=text"
    ),
    "pubchem-phosphoenolpyruvate.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1005/"
        "property/IsomericSMILES,Title/JSON"
    ),
    **{f"{pdb_id}.cif": f"https://files.rcsb.org/download/{pdb_id}.cif" for pdb_id in PDB_IDS},
}

# These are direct author-reported PnPEPCK Table 2 rows. They remain unnormalized
# audit evidence because fixed cosubstrate/metal saturation was not established.
DIRECT_ROWS = (
    (15, 810, 87, 5.8, 0.19, 7.2e3),
    (25, 2500, 650, 23, 1.4, 9.2e3),
    (37, 3400, 290, 92, 3.2, 2.7e4),
    (45, 5900, 1000, 47, 3.0, 8.0e3),
)

STRUCTURES = (
    ("9E32", "holo", "open", 1.53),
    ("9E33", "phosphoenolpyruvate", "open", 1.74),
    ("9E34", "phosphoenolpyruvate", "open", 1.57),
    ("9E35", "beta-sulfopyruvate", "closed", 1.71),
    ("9E36", "beta-sulfopyruvate + GTP", "closed", 1.80),
    ("9E37", "oxalate + GTP", "closed", 2.00),
    ("9E38", "phosphoglycolate + GDP", "closed", 1.52),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.exists():
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)

    supplement = RAW / "PRO-34-e70326-s001.docx"
    if not supplement.exists():
        with zipfile.ZipFile(RAW / "PMC12529880-supplementaryFiles.zip") as archive:
            supplement.write_bytes(archive.read(supplement.name))


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def validate_sources() -> tuple[str, int, str]:
    xml = (RAW / "PMC12529880-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI,
        ACCESSION,
        "SUMO protease results in the full",
        "4\u2009mM MgCl",
        "100\u2009\u03bcM MnCl",
        "1\u2009mM GDP",
        "50\u2009mM KHCO",
        "TABLE 2",
        "https://creativecommons.org/licenses/by/4.0/",
    )
    missing = [text for text in required if text not in xml]
    if missing:
        raise ValueError(f"Article XML lacks expected construct/assay evidence: {missing}")

    sequence = read_fasta(RAW / f"{ACCESSION}.fasta")
    if len(sequence) != 622:
        raise ValueError(f"Expected a 622-aa {ACCESSION} sequence, found {len(sequence)} aa")

    pubchem = json.loads((RAW / "pubchem-phosphoenolpyruvate.json").read_text(encoding="utf-8"))
    compound = pubchem["PropertyTable"]["Properties"][0]
    if compound["CID"] != 1005:
        raise ValueError("Unexpected phosphoenolpyruvate PubChem record")

    with zipfile.ZipFile(RAW / "PRO-34-e70326-s001.docx") as archive:
        document = archive.read("word/document.xml").decode("utf-8")
    for pdb_id in PDB_IDS:
        cif = (RAW / f"{pdb_id}.cif").read_text()
        match = re.search(
            r"_entity_poly\.pdbx_seq_one_letter_code_can\s+;([^;]+);", cif, re.DOTALL
        )
        pdb_sequence = re.sub(r"\s+", "", match.group(1)) if match else ""
        if pdb_id not in document or f"data_{pdb_id}" not in cif:
            raise ValueError(f"Supporting-table/PDB mapping failed for {pdb_id}")
        expected_sequence = sequence if pdb_id != "9E34" else sequence[:8] + "N" + sequence[9:]
        if "_struct_ref.pdbx_db_accession          A1VIE9" not in cif or pdb_sequence != expected_sequence:
            raise ValueError(f"Unexpected PDB sequence mapping for {pdb_id}")
    return sequence, compound["CID"], compound["SMILES"]


def windows_to_wsl(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def run_frozen_mmseqs(query_path: Path) -> None:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP reference FASTA is missing or has changed")
    version = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu-24.04", "--", MMSEQS, "version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 version {version!r} is not the frozen version")

    HOMOLOGY.mkdir(parents=True, exist_ok=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    search_tmp = HOMOLOGY / "search-tmp"
    subprocess.run(
        [
            "wsl.exe", "-d", "Ubuntu-24.04", "--", MMSEQS, "easy-search",
            windows_to_wsl(query_path), windows_to_wsl(REFERENCE), windows_to_wsl(hits),
            windows_to_wsl(search_tmp), "--min-seq-id", "0.3", "-c", "0.8",
            "--cov-mode", "0", "--format-output",
            "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        ],
        check=True,
    )


def write_outputs(*, acquire: bool = True) -> None:
    if acquire:
        download_raw()
    sequence, cid, smiles = validate_sources()
    SOURCE.mkdir(parents=True, exist_ok=True)
    fasta_path = SOURCE / "construct_sequences.fasta"
    fasta_path.write_text(
        f">pnpepck-{ACCESSION.lower()} | exact full-length tag-free catalytic construct\n"
        + "\n".join(sequence[start : start + 80] for start in range(0, len(sequence), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )

    # Frozen homology is intentionally completed before any normalized record is emitted.
    run_frozen_mmseqs(fasta_path)
    hit_path = HOMOLOGY / "homology_hits.tsv"
    hits = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]

    audit_fields = (
        "source_table", "source_row", "organism", "enzyme_identity", "sequence_accession",
        "construct", "reaction", "variable_substrate", "substrate_pubchem_cid",
        "substrate_isomeric_smiles", "temperature_C", "km_uM", "km_sem_uM",
        "kcat_s-1", "kcat_sem_s-1", "kcat_per_km_M-1_s-1", "assay_pH",
        "fixed_GDP_mM", "fixed_MgCl2_mM", "fixed_MnCl2_uM", "fixed_KHCO3_mM",
        "replicates", "saturation_audit", "homology_audit", "curation_status",
    )
    with (SOURCE / "direct-kcat-audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields)
        writer.writeheader()
        for temperature, km, km_sem, kcat, kcat_sem, efficiency in DIRECT_ROWS:
            writer.writerow(
                {
                    "source_table": "Table 2",
                    "source_row": f"{temperature} C; PnPEPCK",
                    "organism": "Polaromonas naphthalenivorans CJ2",
                    "enzyme_identity": "GTP-dependent phosphoenolpyruvate carboxykinase (EC 4.1.1.32)",
                    "sequence_accession": ACCESSION,
                    "construct": "full-length 622-aa native sequence; N-terminal SUMO-His6 removed without residual amino acids",
                    "reaction": "PEP + GDP + CO2 -> oxaloacetate + GTP",
                    "variable_substrate": "phosphoenolpyruvate",
                    "substrate_pubchem_cid": cid,
                    "substrate_isomeric_smiles": smiles,
                    "temperature_C": temperature,
                    "km_uM": km,
                    "km_sem_uM": km_sem,
                    "kcat_s-1": kcat,
                    "kcat_sem_s-1": kcat_sem,
                    "kcat_per_km_M-1_s-1": efficiency,
                    "assay_pH": 7.5,
                    "fixed_GDP_mM": 1,
                    "fixed_MgCl2_mM": 4,
                    "fixed_MnCl2_uM": 100,
                    "fixed_KHCO3_mM": 50,
                    "replicates": "duplicate; error is SEM",
                    "saturation_audit": "fail: fixed GDP, Mg2+, Mn2+, and bicarbonate concentrations reported but not established as saturating for Table 2",
                    "homology_audit": "fail: frozen MMseqs2 development-corpus hit" if hits else "pass",
                    "curation_status": "excluded_unsaturated_fixed_components_and_homology" if hits else "excluded_unsaturated_fixed_components",
                }
            )

    # Keep a schema-bearing accepted file, but do not normalize excluded measurements.
    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerow(("candidate_id", "article_doi", "endpoint", "value", "unit"))

    raw_hashes = {path.name: sha256(path) for path in sorted(RAW.iterdir()) if path.is_file()}
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC12529880",
        "article_doi": DOI,
        "article_published": "2025-10-16",
        "license": "CC-BY-4.0",
        "construct_mapping": {
            "source_accession": ACCESSION,
            "pdb_uniprot_alias": "A1VIE9",
            "native_length_aa": len(sequence),
            "expression_fusion": "N-terminal SUMO 6-His",
            "mature_assayed_protein": "full-length accession sequence after SUMO-protease cleavage; no additional amino acids",
            "construct_sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            "structure_sequence_caveat": "Six PDB entries exactly match WP_011799537.1. PDB 9E34 has an L9N canonical-sequence discrepancy while declaring no mutation; it is not treated as evidence for a kinetic variant.",
        },
        "direct_kcat_rows_reported_for_construct": len(DIRECT_ROWS),
        "accepted_records": 0,
        "normalization_performed": False,
        "exclusion": "Table 2 does not establish saturation of every fixed substrate/metal; frozen MMseqs2 also finds development homology.",
        "raw_sha256": raw_hashes,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    structure_audit = {
        "source": "raw/PRO-34-e70326-s001.docx, Tables S3-S4; RCSB coordinate files",
        "construct_accession": ACCESSION,
        "pdb_uniprot_alias": "A1VIE9",
        "sequence_mapping": "Six PDB canonical entity sequences exactly match the 622-aa WP_011799537.1 sequence. PDB 9E34 differs at L9N despite cross-referencing A1VIE9 and declaring no mutation.",
        "structures": [
            {"pdb_id": pdb_id, "complex": complex_name, "conformation": state, "resolution_A": resolution,
             "sequence_relation": "L9N discrepancy versus WP_011799537.1" if pdb_id == "9E34" else "exact WP_011799537.1 match",
             "raw_file": f"raw/{pdb_id}.cif", "sha256": raw_hashes[f"{pdb_id}.cif"]}
            for pdb_id, complex_name, state, resolution in STRUCTURES
        ],
        "model_predictions_run": False,
    }
    (SOURCE / "structure-audit.json").write_text(
        json.dumps(structure_audit, indent=2) + "\n", encoding="ascii"
    )

    saturation_audit = {
        "source": "raw/PMC12529880-fullText.xml, Section 2.2 and Table 2",
        "standard_reaction_mixture": {
            "buffer": "100 mM HEPES-OH pH 7.5",
            "DTT_mM": 10,
            "NADH_uM": 300,
            "MgCl2_mM": 4,
            "GDP_mM": 1,
            "Mg_to_GDP_ratio": "4:1",
            "MnCl2_uM": 100,
            "PEP_mM": 10,
            "KHCO3_mM": 50,
            "KHCO3_preparation": "bubbled with dry ice",
            "MDH_units": 10,
            "PnPEPCK_ug": 2.5,
            "reaction_volume_mL": 1,
        },
        "table_2_pep_varied": {
            "temperatures_C": [15, 25, 37, 45],
            "variable_component": "PEP; tested concentration range not reported",
            "fixed_components": {
                "GDP_mM": 1,
                "MgCl2_mM": 4,
                "MnCl2_uM": 100,
                "KHCO3_mM": 50,
            },
            "author_saturation_statement_applies": False,
            "decision": "fail",
            "reason": "The methods say other substrates were held constant, but reserve the explicit all-substrates-saturating statement for Arrhenius/Eyring kinetics. Saturation of fixed GDP, bicarbonate, Mg2+, and Mn2+ for Table 2 was not demonstrated.",
        },
        "arrhenius_eyring": {
            "author_saturation_statement_applies": True,
            "statement": "saturating concentrations of all substrates were used",
            "interpreted_standard_concentrations": {
                "PEP_mM": 10,
                "GDP_mM": 1,
                "KHCO3_mM": 50,
                "MgCl2_mM": 4,
                "MnCl2_uM": 100,
            },
            "curation_decision": "not_curated_as_direct_kcat",
            "reason": "Table S1 reports transformed ln(kcat/T) values rather than direct kcat rows; no values were back-calculated or inferred.",
        },
        "model_predictions_run": False,
    }
    (SOURCE / "saturation-audit.json").write_text(
        json.dumps(saturation_audit, indent=2) + "\n", encoding="ascii"
    )

    homology_audit = {
        "audited_on": "2026-07-27",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search",
        "mmseqs_version": MMSEQS_VERSION,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256,
        "construct_sequences_sha256": sha256(fasta_path),
        "homology_hit_sequences": 1 if hits else 0,
        "qualifying_alignments": len(hits),
        "homology_hits": hits,
        "homology_hits_sha256": sha256(hit_path),
        "accepted_records": 0,
        "homology_cold_claimed": False,
        "claim_boundary": "Excluded curation evidence only; no normalization or model predictions were run.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(homology_audit, indent=2) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    write_outputs()
