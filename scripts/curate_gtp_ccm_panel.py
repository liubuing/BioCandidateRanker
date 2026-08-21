from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import urllib.request
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "biostudies-S-EPMC11955521"
DOI = "10.1038/s42003-025-07971-7"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
DEVELOPMENT_CORPUS = Path(
    r"D:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment"
    r"\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json"
)
DEVELOPMENT_CORPUS_SHA256 = "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719"
DEVELOPMENT_CORPUS_SIZE = 12132719
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

PROTEINS = {
    "xk": ("xylulokinase", "CCEL_RS17225", "Ccel_3431", "WP_015926771.1"),
    "pgk": ("phosphoglycerate kinase", "CCEL_RS11420", "Ccel_2260", "WP_015925691.1"),
    "ak": ("acetate kinase", "CCEL_RS10800", "Ccel_2136", "WP_015925578.1"),
}
PRIMERS = {
    "xk": ("TTTTCCATGGCATCTTTTCTTATTGGTATTGATCTAGG", "TTTTCTCGAGTTTCAATATTGTGCTGAGCTGGTTAA"),
    "pgk": ("TTTTCCATGGCAAGCATGATGAACAAAAAAAC", "TTTTCTCGAGTGCTTTTGCGATTATTGAGAA"),
    "ak": ("TTTTTTCCATGGCTAAGGTTTTAGTTATAAATGCGGGGAG", "TTTTTCTCGAGCTTAACCAATCTCACTGTTTCTCTTGC"),
}
URLS = {
    "biostudies-metadata.json": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC11955521",
    "PMC11955521-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11955521/fullTextXML",
    "42003_2025_7971_MOESM1_ESM.pdf": "https://www.ebi.ac.uk/biostudies/files/S-EPMC11955521/42003_2025_7971_MOESM1_ESM.pdf",
    "42003_2025_7971_MOESM2_ESM.xlsx": "https://www.ebi.ac.uk/biostudies/files/S-EPMC11955521/42003_2025_7971_MOESM2_ESM.xlsx",
    **{
        f"{protein_id}.fasta": f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={protein_id}&rettype=fasta&retmode=text"
        for _, _, _, protein_id in PROTEINS.values()
    },
}

SUBSTRATES = {
    "D-xylulose": (5289590, "C([C@H]([C@@H](C(=O)CO)O)O)O"),
    "GTP": (135398633, "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N"),
    "ATP": (5957, "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"),
    "GDP": (135398619, "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N"),
    "ADP": (6022, "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N"),
    "3-phosphoglycerate": (439183, "C([C@H](C(=O)O)O)OP(=O)(O)O"),
    "acetyl phosphate": (186, "CC(=O)OP(=O)(O)O"),
    "acetate": (176, "CC(=O)O"),
    "phosphoenolpyruvate": (1005, "C=C(C(=O)O)OP(=O)(O)O"),
}
URLS.update({
    f"pubchem-{cid}.json":
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,Title/JSON"
    for cid, _ in SUBSTRATES.values()
})


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
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def constructs() -> dict[str, str]:
    result = {}
    for key, (_, _, _, protein_id) in PROTEINS.items():
        native = read_fasta(RAW / f"{protein_id}.fasta")
        forward, reverse = PRIMERS[key]
        if "CCATGG" not in forward or "CTCGAG" not in reverse:
            raise ValueError(f"{key}: expected NcoI/XhoI cloning junctions are absent")
        # NcoI contributes MA and the reverse primer has no stop; XhoI/vector contributes LEHis6.
        result[f"{key}-pet28a-his6"] = "MA" + native[1:] + "LEHHHHHH"
    return result


def validate_raw() -> dict[str, str]:
    missing = [name for name in URLS if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing raw evidence: {missing}; run with --download")
    metadata = json.loads((RAW / "biostudies-metadata.json").read_text(encoding="utf-8"))
    if metadata.get("accno") != "S-EPMC11955521":
        raise ValueError("unexpected BioStudies accession")
    xml = (RAW / "PMC11955521-fullText.xml").read_text(encoding="utf-8")
    markers = (
        DOI, "creativecommons.org/licenses/by/4.0", "CCEL_RS17225/Ccel_3431",
        "CCEL_RS11420/Ccel_2260", "CCEL_RS12995/Ccel_2569", "CCEL_RS10800/Ccel_2136",
        "pET28a-XK", "pET28a-PGK", "pET28a-PK", "pET28a-AK",
        "fixed, saturating concentration of the co-substrate", "iu/µmol",
    )
    absent = [marker for marker in markers if marker not in xml]
    if absent:
        raise ValueError(f"primary source markers changed or are absent: {absent}")
    if (RAW / "42003_2025_7971_MOESM1_ESM.pdf").read_bytes()[:5] != b"%PDF-":
        raise ValueError("supporting primer table is not a PDF")
    for cid, expected_smiles in SUBSTRATES.values():
        item = json.loads((RAW / f"pubchem-{cid}.json").read_text(encoding="utf-8"))[
            "PropertyTable"
        ]["Properties"][0]
        if item["CID"] != cid or item["SMILES"] != expected_smiles:
            raise ValueError(f"PubChem identity changed for CID {cid}")
    return constructs()


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive.rstrip(':').lower()}/{resolved.as_posix()[3:]}"


def mmseqs_command(*arguments: str) -> list[str]:
    _, distribution, executable = MMSEQS.split(":", 2)
    return ["wsl", "-d", distribution, executable, *arguments]


def run_checked(arguments: list[str]) -> str:
    completed = subprocess.run(arguments, check=False, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f"command failed: {' '.join(arguments)}\n{completed.stderr}")
    return completed.stdout or completed.stderr


def reset_wsl(path: Path) -> None:
    if path.exists():
        run_checked(["wsl", "-d", MMSEQS.split(":", 2)[1], "rm", "-rf", "--", wsl_path(path)])


def run_homology(query: Path, threads: int) -> str:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("frozen development reference is absent or changed")
    version = run_checked(mmseqs_command("version")).strip().splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 {version!r} != frozen {MMSEQS_VERSION!r}")
    reset_wsl(HOMOLOGY)
    HOMOLOGY.mkdir(parents=True)
    combined = HOMOLOGY / "homology_hits.tsv"
    combined.write_text("", encoding="ascii")
    entries = read_fasta_entries(query)
    for family in PROTEINS:
        family_dir = HOMOLOGY / family
        family_dir.mkdir()
        family_fasta = family_dir / "query.fasta"
        sequence_id = f"{family}-pet28a-his6"
        family_fasta.write_text(f">{sequence_id}\n{entries[sequence_id]}\n", encoding="ascii")
        hits = family_dir / "homology_hits.tsv"
        run_checked(mmseqs_command(
            "easy-search", wsl_path(family_fasta), wsl_path(REFERENCE), wsl_path(hits),
            wsl_path(family_dir / "search-tmp"), "--min-seq-id", "0.3", "-c", "0.8",
            "--cov-mode", "0", "--format-output",
            "query,target,fident,alnlen,qcov,tcov,evalue,bits", "--threads", str(threads),
        ))
        with combined.open("a", encoding="ascii", newline="\n") as output:
            output.write(hits.read_text(encoding="ascii"))
    prefix = HOMOLOGY / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    run_checked(mmseqs_command(
        "easy-linclust", wsl_path(query), wsl_path(prefix), wsl_path(HOMOLOGY / "cluster-tmp"),
        "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads),
    ))
    return version


def read_fasta_entries(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    identifier = ""
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            identifier = line[1:].split()[0]
            result[identifier] = ""
        else:
            result[identifier] += line.strip()
    return result


def exact_overlap(sequences: dict[str, str]) -> dict[str, int | str]:
    if DEVELOPMENT_CORPUS.stat().st_size != DEVELOPMENT_CORPUS_SIZE:
        raise ValueError("frozen development corpus size changed")
    if sha256(DEVELOPMENT_CORPUS) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("frozen development corpus SHA256 changed")
    development = json.loads(DEVELOPMENT_CORPUS.read_text(encoding="utf-8"))
    sequence_set = set(sequences.values())
    smiles_set = {smiles for _, smiles in SUBSTRATES.values()}
    exact_sequence_rows = [row for row in development if row.get("Sequence") in sequence_set]
    return {
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_rows": len(development),
        "exact_sequence_rows": len(exact_sequence_rows),
        "exact_sequence_substrate_overlap": sum(
            row.get("Smiles") in smiles_set for row in exact_sequence_rows
        ),
    }


def specifications() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(family: str, varied: str, fixed: str, fixed_mm: float | str, kcat: float,
            error: float, km: float, disposition: str, condition: str = "") -> None:
        rows.append({
            "family_key": family, "varied": varied, "fixed": fixed, "fixed_mm": fixed_mm,
            "kcat": kcat, "error": error, "km": km, "disposition": disposition,
            "condition": condition,
        })

    for varied, fixed, fixed_mm, kcat, error, km in (
        ("GTP", "D-xylulose", 2.5, 6540, 177, 0.12),
        ("ATP", "D-xylulose", 2.5, 6653, 281, 1.18),
        ("D-xylulose", "GTP", 10, 8788, 375, 0.30),
        ("D-xylulose", "ATP", 25, 7489, 147, 0.50),
    ):
        add("xk", varied, fixed, fixed_mm, kcat, error, km, "qualified_source_declares_fixed_cosubstrate_saturating")
    add("pgk", "GDP", "1,3-bisphosphoglycerate", "generated in situ", 771, 20, 0.135,
        "excluded_real_fixed_substrate_concentration_unknown_in_coupled_assay")
    add("pgk", "ADP", "1,3-bisphosphoglycerate", "generated in situ", 733, 12, 0.08,
        "excluded_real_fixed_substrate_concentration_unknown_in_coupled_assay")
    for varied, fixed, fixed_mm, kcat, error, km in (
        ("GTP", "3-phosphoglycerate", 25, 2831, 348, 4.66),
        ("ATP", "3-phosphoglycerate", 4, 4256, 148, 0.48),
        ("3-phosphoglycerate", "GTP", 25, 3273, 365, 6.90),
        ("3-phosphoglycerate", "ATP", 10, 5767, 194, 0.65),
    ):
        add("pgk", varied, fixed, fixed_mm, kcat, error, km, "qualified_source_declares_fixed_cosubstrate_saturating")
    for fbp, values in (
        (0, (("GDP", 977, 48, 0.90), ("ADP", 9508, 328, 1.87),
             ("phosphoenolpyruvate", 3084, 312, 23.0), ("phosphoenolpyruvate", 9212, 50, 1.56))),
        (0.5, (("GDP", 2341, 102, 0.92), ("ADP", 9213, 245, 1.74),
               ("phosphoenolpyruvate", 3794, 218, 6.40), ("phosphoenolpyruvate", 9814, 152, 1.64))),
        (3, (("GDP", 5747, 218, 1.06), ("ADP", 9032, 328, 1.75),
             ("phosphoenolpyruvate", 6545, 146, 2.75), ("phosphoenolpyruvate", 9585, 143, 1.52))),
        (10, (("GDP", 7479, 287, 0.89), ("ADP", 8846, 383, 1.47),
              ("phosphoenolpyruvate", 7747, 152, 2.30), ("phosphoenolpyruvate", 9335, 90, 1.60))),
    ):
        for index, (varied, kcat, error, km) in enumerate(values):
            fixed = "phosphoenolpyruvate" if index < 2 else ("GDP" if index == 2 else "ADP")
            add("pk", varied, fixed, 10 if index < 2 else 5, kcat, error, km,
                "excluded_exact_pyruvate_kinase_construct_unresolved", f"FBP={fbp} mM")
    for varied, fixed, fixed_mm, kcat, error, km in (
        ("GDP", "acetyl phosphate", 5, 35365, 1150, 1.47),
        ("ADP", "acetyl phosphate", 5, 25437, 1668, 1.13),
        ("acetyl phosphate", "GDP", 10, 32137, 823, 1.28),
        ("acetyl phosphate", "ADP", 10, 25112, 1755, 1.49),
        ("ATP", "acetate", 10, 16.68, 0.74, 1.47),
        ("acetate", "ATP", 10, 12.85, 0.85, 2.37),
    ):
        add("ak", varied, fixed, fixed_mm, kcat, error, km, "qualified_source_declares_fixed_cosubstrate_saturating")
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(*, acquire: bool, homology: bool, threads: int) -> None:
    if acquire:
        download()
    sequences = validate_raw()
    overlap = exact_overlap(sequences)
    SOURCE.mkdir(parents=True, exist_ok=True)
    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in sequences.items():
            handle.write(f">{sequence_id} exact NcoI/XhoI pET28a product\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    version = run_homology(fasta, threads) if homology else MMSEQS_VERSION
    hits_path = HOMOLOGY / "homology_hits.tsv"
    cluster_path = HOMOLOGY / "family-cluster" / "proteins_cluster.tsv"
    if not hits_path.is_file() or not cluster_path.is_file():
        raise RuntimeError("pinned family-wise MMseqs evidence is required before normalization")
    hit_lines = [line for line in hits_path.read_text(encoding="ascii").splitlines() if line]
    hit_ids = {line.split("\t", 1)[0] for line in hit_lines}
    cluster_lines = [line for line in cluster_path.read_text(encoding="ascii").splitlines() if line]

    records = []
    for index, spec in enumerate(specifications(), 1):
        family = str(spec["family_key"])
        sequence_id = f"{family}-pet28a-his6" if family in PROTEINS else ""
        disposition = str(spec["disposition"])
        qualified = disposition.startswith("qualified_")
        status = (
            "excluded_exact_construct_unresolved" if family == "pk" else
            "excluded_assay_saturation" if not qualified else
            "excluded_homology" if sequence_id in hit_ids else
            "accepted_homology_cold_pool"
        )
        varied = str(spec["varied"])
        cid, smiles = SUBSTRATES[varied]
        protein = PROTEINS.get(family)
        records.append({
            "candidate_id": f"gtpccm-{index:03d}", "article_doi": DOI,
            "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11955521",
            "source_file": "raw/PMC11955521-fullText.xml", "source_table": "Table 1",
            "source_row": f"{family}; {varied}; {spec['condition']}",
            "organism": "Ruminiclostridium cellulolyticum H10",
            "enzyme_identity": "pyruvate kinase" if family == "pk" else protein[0],
            "locus_tag": "CCEL_RS12995" if family == "pk" else protein[1],
            "protein_accession": "unresolved assay product" if family == "pk" else protein[3],
            "sequence_id": sequence_id,
            "construct": (
                "not reconstructed: the supplementary primer table does not map the stated pyruvate kinase primers"
                if family == "pk" else
                "NcoI product MA + native residues 2-end + XhoI/vector-derived LEHHHHHH"
            ),
            "variable_substrate": varied, "substrate_pubchem_cid": cid,
            "substrate_isomeric_smiles": smiles, "fixed_reaction_substrate": spec["fixed"],
            "fixed_reaction_substrate_mM": spec["fixed_mm"], "endpoint": "kcat_s-1",
            "reported_kcat_iu_per_umol": spec["kcat"],
            "reported_kcat_sd_iu_per_umol": spec["error"],
            "kcat_s-1": float(spec["kcat"]) / 60,
            "kcat_sd_s-1": float(spec["error"]) / 60,
            "km_mM": spec["km"], "assay_temperature_C": 37,
            "saturation_disposition": disposition, "status_at_normalization": status,
        })
    write_csv(SOURCE / "candidate_records.csv", records)
    exclusions = [
        {**record, "exclusion_reason": record["status_at_normalization"].removeprefix("excluded_")}
        for record in records if record["status_at_normalization"] != "accepted_homology_cold_pool"
    ]
    write_csv(SOURCE / "exclusions.csv", exclusions)
    write_csv(SOURCE / "direct-kcat-rows.csv", [
        {**record, "candidate_label_created": record["status_at_normalization"] == "accepted_homology_cold_pool"}
        for record in records
    ])
    saturation_rows = [{
        "family": family, "source_rows": sum(row["family_key"] == family for row in specifications()),
        "disposition": (
            "two coupled forward PGK rows excluded because generated 1,3-BPG concentration is unknown; four reverse rows pass source-declared saturation"
            if family == "pgk" else
            "all finite rows excluded before saturation admission because exact construct is unresolved"
            if family == "pk" else
            "all finite rows pass the article's explicit fixed-saturating-cosubstrate statement"
        ),
    } for family in ("xk", "pgk", "pk", "ak")]
    write_csv(SOURCE / "saturation-audit.csv", saturation_rows)

    family_audit = []
    for family in ("xk", "pgk", "pk", "ak"):
        family_rows = [row for row in records if row["source_row"].startswith(f"{family};")]
        family_hits = [line.split("\t") for line in hit_lines if line.startswith(f"{family}-")]
        family_audit.append({
            "family": family, "reported_finite_kcat_rows": len(family_rows),
            "exact_construct_resolved": family != "pk", "homology_search_run": family != "pk",
            "homology_hit_sequences": len({line[0] for line in family_hits}),
            "accepted_records": sum(row["status_at_normalization"] == "accepted_homology_cold_pool" for row in family_rows),
            "disposition_counts": {
                status: sum(row["status_at_normalization"] == status for row in family_rows)
                for status in sorted({str(row["status_at_normalization"]) for row in family_rows})
            },
            "hit_file": f"homology/{family}/homology_hits.tsv" if family != "pk" else None,
            "hit_sha256": sha256(HOMOLOGY / family / "homology_hits.tsv") if family != "pk" else None,
        })
    (SOURCE / "family-admission-audit.json").write_text(
        json.dumps(family_audit, indent=2) + "\n", encoding="ascii"
    )
    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )
    provenance = {
        "source_id": SOURCE_ID, "article_doi": DOI, "article_published": "2025-03-30",
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11955521",
        "license": "CC-BY-4.0", "reported_finite_direct_kcat_rows": len(records),
        "accepted_records": len(accepted), "accepted_families": sorted({row["enzyme_identity"] for row in accepted}),
        "construct_resolution": {
            "resolved": ["xylulokinase", "phosphoglycerate kinase", "acetate kinase"],
            "formula": "M + NcoI-derived A + native residues 2-end + XhoI/vector-derived LEHHHHHH",
            "excluded": "Pyruvate kinase: no matching Ccel_2569 primer pair occurs in Table S1; no product was inferred.",
        },
        "saturation_policy": "The article states fixed cosubstrates were saturating. PGK coupled forward rows still fail because the actual generated 1,3-BPG concentration is unknown.",
        "unit_conversion": "Reported IU/umol equals umol product/min per umol enzyme, hence min-1. Normalized kcat_s-1 and SD are the reported values divided by exactly 60; reported values are retained in dedicated columns. No values were inferred or refitted.",
        "exact_overlap_audit": overlap, "model_predictions_run": False,
        "artifacts": {name: sha256(SOURCE / name) for name in (
            "candidate_records.csv", "exclusions.csv", "direct-kcat-rows.csv",
            "construct_sequences.fasta", "saturation-audit.csv", "family-admission-audit.json",
        )},
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )
    audit = {
        "audited_on": date.today().isoformat(), "source_id": SOURCE_ID, "article_doi": DOI,
        "method": "three family-independent pinned MMseqs2 easy-search runs and source-local easy-linclust",
        "mmseqs_version": version, "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "candidate_records": len(records),
        "construct_sequences_sha256": sha256(fasta), "homology_hits_sha256": sha256(hits_path),
        "homology_hit_sequences": len(hit_ids), "candidate_mmseqs_families": len({line.split("\t")[0] for line in cluster_lines}),
        "accepted_catalytic_families": len({row["enzyme_identity"] for row in accepted}),
        "family_cluster_sha256": sha256(cluster_path), "accepted_records": len(accepted),
        "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
        "accepted_unique_substrates": len({row["substrate_pubchem_cid"] for row in accepted}),
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"],
        "family_admission_audit_sha256": sha256(SOURCE / "family-admission-audit.json"),
        "readiness_gate_passes": False, "model_predictions_run": False,
        "claim_boundary": "Family-independent curation pool only; no model predictions or pseudolabels were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--run-homology", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    curate(acquire=args.download, homology=args.run_homology, threads=args.threads)


if __name__ == "__main__":
    main()
