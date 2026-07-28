from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "europepmc-PMC12676663"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1021/acs.biochem.5c00561"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

PARENT_SEQUENCE = (
    "MLPKLVITHRVHDEILQLLAPHCELMTNQTDSTLTREEILRRCRDAQAMMAFMPDRVDADFLQACPELRV"
    "VGCALKGFDNFDVDACTARGVWLTFVPDLLTVPTAELAIGLAVGLGRHLRAADAFVRSGEFQGWQPQFYG"
    "TGLDNATVGILGMGAIGLAMADRLQGWGATLQYHEAKALDTQTEQRLGLRQVACSELFASSDFILLALPL"
    "NADTQHLVNAELLALVRPGALLVNPCRGSVVDEAAVLAALERGQLGGYAADVFEMEDWARADRPRLIDPA"
    "LLAHPNTLFTPHIGSAVRAVRLEIERCAAQNIIQVLAGARPINAANRLPKAEPAAC"
)
MUTATIONS = {
    13: "E", 26: "I", 71: "I", 130: "K", 132: "R", 137: "R", 150: "F",
    175: "A", 215: "L", 275: "Q", 276: "Q", 313: "L", 315: "A", 319: "E",
    325: "V", 332: "N", 336: "D",
}
RESIDUE_MASSES = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
STRUCTURES = {
    "phosphite dianion": (177605, "OP([O-])[O-]"),
    "NAD+": (
        5892,
        "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])"
        "OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
    ),
}
URLS = {
    "PMC12676663-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12676663/fullTextXML"
    ),
    "O69054.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=protein&id=21362789&rettype=fasta&retmode=text"
    ),
    "pubchem-phosphite.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/177605/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "pubchem-NAD.json": (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5892/"
        "property/IsomericSMILES,Title/JSON"
    ),
    "4E5K-polymer-entity.json": "https://data.rcsb.org/rest/v1/core/polymer_entity/4E5K/1",
    "4E5K.pdb": "https://files.rcsb.org/download/4E5K.pdb",
    "PMC12676663-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12676663/supplementaryFiles"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)
    supplement = RAW / "bi5c00561_si_001.pdf"
    if not supplement.is_file():
        with zipfile.ZipFile(RAW / "PMC12676663-supplementaryFiles.zip") as archive:
            member = next(
                item for item in archive.namelist() if item.endswith("bi5c00561_si_001.pdf")
            )
            with archive.open(member) as extracted, supplement.open("wb") as output:
                shutil.copyfileobj(extracted, output)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def exact_construct() -> str:
    if len(PARENT_SEQUENCE) != 336:
        raise ValueError("Unexpected O69054 parent sequence length")
    sequence = list(PARENT_SEQUENCE)
    expected = "DMVEQQIEQRLIVAAEC"
    observed = "".join(sequence[position - 1] for position in MUTATIONS)
    if observed != expected:
        raise ValueError(f"Mutation parent residues differ: {observed}")
    for position, residue in MUTATIONS.items():
        sequence[position - 1] = residue
    construct = "GSH" + "".join(sequence)
    average_mass = sum(RESIDUE_MASSES[residue] for residue in construct) + 18.01528
    if len(construct) != 339 or abs(average_mass - 36764.4) > 0.1:
        raise ValueError(f"Construct mass validation failed: {len(construct)} aa, {average_mass} Da")
    return construct


def validate_raw() -> tuple[str, dict[str, tuple[int, str]]]:
    xml = (RAW / "PMC12676663-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI, "creativecommons.org/licenses/by/4.0/", "D13E", "C336D", "E175A",
        "Gly-Ser-His tripeptide", "36764.1", "36764.4", "2.3", "2.1", "0.082",
        "4.1", "1 mM NAD", "20 mM HPO",
    )
    missing = [marker for marker in required if marker not in xml]
    if missing:
        raise ValueError(f"Article XML lacks expected evidence: {missing}")
    if read_fasta(RAW / "O69054.fasta") != PARENT_SEQUENCE:
        raise ValueError("NCBI O69054 sequence differs from frozen parent sequence")

    observed_structures = {}
    for name, raw_name in (
        ("phosphite dianion", "pubchem-phosphite.json"), ("NAD+", "pubchem-NAD.json")
    ):
        compound = json.loads((RAW / raw_name).read_text(encoding="utf-8"))[
            "PropertyTable"
        ]["Properties"][0]
        expected_cid, expected_smiles = STRUCTURES[name]
        if compound["CID"] != expected_cid or compound["SMILES"] != expected_smiles:
            raise ValueError(f"Unexpected PubChem mapping for {name}")
        observed_structures[name] = (compound["CID"], compound["SMILES"])

    pdb = json.loads((RAW / "4E5K-polymer-entity.json").read_text(encoding="utf-8"))
    pdb_sequence = pdb["entity_poly"]["pdbx_seq_one_letter_code_can"]
    construct = exact_construct()
    # 4E5K is the 16X parent of the assayed 17X protein, lacking only E175A.
    expected_16x = list(construct[3:-7])
    expected_16x[174] = "E"
    if pdb_sequence != "".join(expected_16x):
        raise ValueError("PDB 4E5K does not match residues 1-329 of the expected 16X construct")
    return construct, observed_structures


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
    construct, structures = validate_raw()
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "construct_sequences.fasta"
    query.write_text(
        ">17x-psptdh-gsh exact assayed 17X-PsPTDH; retained thrombin GSH scar\n"
        + "\n".join(construct[start : start + 80] for start in range(0, len(construct), 80))
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    if run_homology:
        run_mmseqs(query)

    hit_path = HOMOLOGY / "homology_hits.tsv"
    cluster_path = HOMOLOGY / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        raise RuntimeError("Pinned MMseqs evidence is required before kinetics normalization")
    hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    status = "excluded_homology" if hit_lines else "accepted_homology_cold_pool"

    common = {
        "article_doi": DOI,
        "stable_record_url": "https://europepmc.org/articles/PMC12676663",
        "source_file": "raw/PMC12676663-fullText.xml",
        "organism": "Pseudomonas stutzeri WM88 (current genus: Stutzerimonas)",
        "enzyme_identity": "17X-PsPTDH phosphite dehydrogenase (EC 1.20.1.1)",
        "sequence_accession": "O69054",
        "sequence_id": "17x-psptdh-gsh",
        "construct": "GSH thrombin scar + O69054 with 17 substitutions; 339 aa",
        "reaction": "phosphite + NAD+ + H2O -> phosphate + NADH + H+",
        "endpoint": "kcat_s-1",
        "assay_buffer": "36 mM triethanolamine",
        "assay_pH": 7.5,
        "assay_temperature_C": 25,
        "ionic_strength_M": 0.3,
        "ionic_strength_salt": "NaCl",
        "enzyme_concentration_uM": "approximately 0.3",
        "assay_volume_mL": 1.0,
        "assay_method": "A340 initial-rate assay; NADH epsilon=6220 M-1 cm-1",
        "fit_method": "nonlinear least-squares Michaelis-Menten fit in GraphPad Prism 10",
        "status_at_normalization": status,
    }
    rows = []
    specifications = (
        ("ptdh-001", "Figure S2", "phosphite dianion", 2.3, 0.02, 4.1, 0.1,
         "NAD+", 1.0, 0.082, 12.195122, "saturating"),
        ("ptdh-002", "Figure S3", "NAD+", 2.1, 0.02, 0.082, 0.003,
         "phosphite dianion", 20.0, 4.1, 4.878049, "nearly saturating"),
    )
    for candidate_id, source_figure, variable, kcat, kcat_error, km, km_error, fixed, fixed_mm, fixed_km, multiple, wording in specifications:
        cid, smiles = structures[variable]
        rows.append(
            {
                "candidate_id": candidate_id,
                **common,
                "source_section": "Experimental Section 2.2.1 and Results",
                "source_row": f"{source_figure}; {variable} varied; reciprocal saturation experiment",
                "variable_substrate": variable,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "kcat_s-1": kcat,
                "kcat_error_s-1": kcat_error,
                "kcat_error_type": "reported fit uncertainty (+/-); statistical type not specified",
                "km_mM": km,
                "km_error_mM": km_error,
                "fixed_cosubstrate": fixed,
                "fixed_cosubstrate_concentration_mM": fixed_mm,
                "fixed_cosubstrate_measured_km_mM": fixed_km,
                "fixed_cosubstrate_multiple_of_km": multiple,
                "author_saturation_wording": wording,
                "limiting_substrate_consumption": "<=10%",
                "monitoring_time_min": 10,
            }
        )
    fields = list(rows[0])
    write_csv(SOURCE / "candidate_records.csv", rows, fields)

    exclusion_fields = [
        "source_row", "reported_endpoint", "value", "unit", "exclusion_reason",
        "candidate_label_created",
    ]
    exclusions = [
        {
            "source_row": "Table 1; all six rows",
            "reported_endpoint": "third- or fourth-order rate constant",
            "value": "0.78 to 6.8e6; 70 to 4400",
            "unit": "M-2 s-1 or M-3 s-1",
            "exclusion_reason": "not an accepted kcat, Km, or kcat/Km endpoint",
            "candidate_label_created": False,
        },
        {
            "source_row": "Fragment-activated NR/NMN experiments",
            "reported_endpoint": "initial-rate slopes and activation parameters",
            "value": "",
            "unit": "",
            "exclusion_reason": "no direct finite kcat endpoint reported",
            "candidate_label_created": False,
        },
        {
            "source_row": "Earlier-study comparison in Results",
            "reported_endpoint": "kcat",
            "value": 3.3,
            "unit": "s-1",
            "exclusion_reason": "not measured in the audited primary article",
            "candidate_label_created": False,
        },
    ]
    write_csv(SOURCE / "exclusions.csv", exclusions, exclusion_fields)

    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )
    accepted = [row for row in rows if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://europepmc.org/articles/PMC12676663",
        "article_doi": DOI,
        "article_published": "2025-11-11",
        "license": "CC-BY-4.0",
        "license_scope": "Article and embedded supporting information",
        "reported_direct_finite_kcat_rows": 2,
        "curated_reciprocal_kcat_rows": 2,
        "accepted_records": len(accepted),
        "kinetics_source": "raw/PMC12676663-fullText.xml, Experimental Section 2.2.1, Results, Figures S2-S3",
        "sequence_sources": [
            "raw/O69054.fasta (336-aa wild-type parent)",
            "raw/PMC12676663-fullText.xml (17 substitutions and retained GSH scar)",
            "raw/4E5K-polymer-entity.json (independent 16X sequence cross-check)",
            "raw/bi5c00561_si_001.pdf, Figure S1 (assayed-protein mass spectrum)",
        ],
        "exact_construct_recovery": {
            "formula": "GSH + O69054[D13E,M26I,V71I,E130K,Q132R,Q137R,I150F,E175A,Q215L,R275Q,L276Q,I313L,V315A,A319E,A325V,E332N,C336D]",
            "length_aa": len(construct),
            "calculated_average_mass_Da": 36764.37,
            "article_calculated_mass_Da": 36764.4,
            "article_observed_mass_Da": 36764.1,
            "mass_difference_observed_Da": 0.27,
        },
        "structure_mapping": {
            name: {"pubchem_cid": cid, "isomeric_smiles": smiles}
            for name, (cid, smiles) in structures.items()
        },
        "experimental_structure": {
            "pdb_id": "4E5K",
            "doi": "10.2210/pdb4e5k/pdb",
            "description": "16X thermostable PsPTDH residues 1-329 with NAD+ and sulfite; differs from assayed 17X by E175A and lacks GSH plus seven C-terminal residues",
            "use": "provenance/context only; not substituted for the exact assay sequence",
        },
        "reciprocal_saturation": [
            {"variable": "phosphite", "fixed": "1.0 mM NAD+", "fixed_over_measured_Km": 12.195122, "author_wording": "saturating"},
            {"variable": "NAD+", "fixed": "20 mM phosphite", "fixed_over_measured_Km": 4.878049, "author_wording": "nearly saturating"},
        ],
        "raw_file_hashes": "raw-file-hashes.json",
        "homology_status": status,
        "final_disposition": "excluded_by_frozen_homology_gate" if hit_lines else status,
        "exact_sequence_substrate_overlap": "not evaluated after mandatory homology exclusion" if hit_lines else "pending global overlap audit",
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    families = {
        line.split("\t", 1)[0]
        for line in cluster_path.read_text(encoding="utf-8").splitlines() if line
    }
    audit = {
        "audited_on": "2026-07-27",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": MMSEQS_VERSION,
        "wsl_distribution": "Ubuntu-24.04",
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "search_command": "mmseqs easy-search construct_sequences.fasta unique_proteins.fasta homology_hits.tsv search-tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0 --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "cluster_command": "mmseqs linclust input_db clusters_db cluster_tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0",
        "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": {"count": len(rows), "sha256": sha256(SOURCE / "candidate_records.csv")},
        "unique_sequences": 1,
        "construct_sequences_sha256": sha256(query),
        "homology_hit_sequences": 1 if hit_lines else 0,
        "homology_hit_alignments": len(hit_lines),
        "homology_hits_sha256": sha256(hit_path),
        "candidate_mmseqs_families": len(families),
        "family_cluster_sha256": sha256(cluster_path),
        "accepted_records": len(accepted),
        "accepted_unique_sequences": 1 if accepted else 0,
        "accepted_unique_substrates": len({row["variable_substrate"] for row in accepted}),
        "status": status,
        "exclusion_reason": "At least one frozen-threshold UniKP development hit" if hit_lines else None,
        "exact_sequence_substrate_overlap": "not evaluated after mandatory homology exclusion" if hit_lines else "pending global overlap audit",
        "readiness_gate_passes": False,
        "claim_boundary": "Curation/exclusion evidence only; no model predictions were generated.",
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
