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
SOURCE_ID = "biostudies-S-EPMC12667214"
DOI = "10.1002/2211-5463.70103"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
REFERENCE = ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723" / "homology" / "unikp_reference.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"
DEVELOPMENT_CORPUS_SHA256 = "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719"
DEVELOPMENT_CORPUS_SIZE = 12132719
DEFAULT_DEVELOPMENT_CORPUS = Path(r"D:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json")

PROTEINS = {
    "ampk": ("AMP kinase", "TTHA1671", "Q5SHQ9", "pET-3a/ttha1671", 34, "S", "A", "E", 25),
    "cmpk": ("CMP kinase", "TTHA0458", "Q5SL35", "pET-11a/ttha0458", 11, "S", "A", "E", 25),
    "umpk": ("UMP kinase", "TTHA0859", "P43891", "pET-11a/ttha0859", 11, "S", "A", "E", 25),
    "idh": ("isocitrate dehydrogenase", "TTHA1535", "P33197", "pET-11a/ttha1535", 98, "S", "A", "E", 60),
    "mdh": ("malate dehydrogenase", "TTHA0536", "Q5SKV7", "pET-11a/ttha0536", 236, "S", "A", "E", 30),
    "ppase": ("inorganic pyrophosphatase", "TTHA1965", "P38576", "pET-11a/ttha1965", 140, "Y", "F", "E", 60),
}

SUBSTRATES = {
    "AMP": (6083, "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N"),
    "CMP": (6131, "C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O"),
    "UMP": (6030, "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O"),
    "ATP": (5957, "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"),
    "isocitrate": (1198, "C(C(C(C(=O)O)O)C(=O)O)C(=O)O"),
    "NADP+": (5885, "C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)N"),
    "oxaloacetate": (970, "C(C(=O)C(=O)O)C(=O)O"),
    "NADH": (439153, "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"),
    "pyrophosphate": (1023, "OP(=O)(O)OP(=O)(O)O"),
}

URLS = {
    "biostudies-metadata.json": "https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC12667214",
    "PMC12667214-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12667214/fullTextXML",
    "FEB4-15-1987-s001.docx": "https://www.ebi.ac.uk/biostudies/files/S-EPMC12667214/FEB4-15-1987-s001.docx",
    **{f"{accession}.fasta": f"https://rest.uniprot.org/uniprotkb/{accession}.fasta" for _, _, accession, *_ in PROTEINS.values()},
    **{f"{accession}.json": f"https://rest.uniprot.org/uniprotkb/{accession}.json" for _, _, accession, *_ in PROTEINS.values()},
    **{f"pubchem-{cid}.json": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON" for cid, _ in SUBSTRATES.values()},
    **{f"rcsb-{pdb_id}-entity-1.json": f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1" for pdb_id in ("3CM0", "3AKC", "8YH1", "2D1C", "2PRD")},
}

PDB_SEQUENCE_CHECKS = {
    "3CM0": ("ampk", "exact_full_ORF"),
    "3AKC": ("cmpk", "exact_full_ORF"),
    "8YH1": ("umpk", "exact_full_ORF"),
    "2D1C": ("idh", "exact_full_ORF"),
    "2PRD": ("ppase", "native_residues_2_end_missing_initiator_Met"),
}

# family, variant, varied substrate, fixed substrate, fixed uM, Km, SD, kcat, SD, assay disposition
ROWS = (
    ("ampk", "WT", "AMP", "ATP", 1000, 23, 1, 111, 3, "qualified_author_designates_fixed_ATP_excess"),
    ("ampk", "S34A", "AMP", "ATP", 1000, 71, 2, 142, 4, "qualified_author_designates_fixed_ATP_excess"),
    ("ampk", "S34E", "AMP", "ATP", 1000, 153, 1, 3.7, 0, "qualified_author_designates_fixed_ATP_excess"),
    ("cmpk", "WT", "CMP", "ATP", 1000, 61, 19, 54, 16, "qualified_reciprocal_saturation"),
    ("cmpk", "WT", "ATP", "CMP", 1000, 53, 1, 40, 1, "qualified_reciprocal_saturation"),
    ("cmpk", "S11A", "CMP", "ATP", 1000, 53, 3, 57, 8, "qualified_reciprocal_saturation"),
    ("cmpk", "S11A", "ATP", "CMP", 1000, 73, 16, 59, 4, "qualified_reciprocal_saturation"),
    ("cmpk", "S11E", "CMP", "ATP", 1000, 79, 38, 0.11, 0.04, "qualified_reciprocal_saturation"),
    ("cmpk", "S11E", "ATP", "CMP", 1000, 60, 14, 0.05, 0.01, "qualified_reciprocal_saturation"),
    ("umpk", "WT", "UMP", "ATP", 1000, 74, 10, 21, 1, "excluded_fixed_ATP_only_4.7x_Km"),
    ("umpk", "WT", "ATP", "CMP", 1000, 214, 7, 14, 0, "excluded_fixed_component_reported_as_CMP_not_UMP"),
    ("umpk", "S11A", "UMP", "ATP", 1000, 66, 1, 0.07, 0.01, "excluded_fixed_ATP_only_1.37x_reciprocal_Km"),
    ("umpk", "S11A", "ATP", "CMP", 1000, 732, 206, 0.04, 0.02, "excluded_fixed_component_reported_as_CMP_not_UMP"),
    ("idh", "WT", "isocitrate", "NADP+", 4000, 34, 2, 26, 1, "qualified_reciprocal_saturation"),
    ("idh", "WT", "NADP+", "isocitrate", 2000, 39, 1, 30, 2, "qualified_reciprocal_saturation"),
    ("idh", "S98A", "isocitrate", "NADP+", 4000, 14, 6, 4.4, 0.1, "qualified_reciprocal_saturation"),
    ("idh", "S98A", "NADP+", "isocitrate", 2000, 150, 2, 5.0, 0.9, "qualified_reciprocal_saturation"),
    ("idh", "S98E", "isocitrate", "NADP+", 4000, 178, 28, 0.19, 0.02, "qualified_reciprocal_saturation"),
    ("idh", "S98E", "NADP+", "isocitrate", 2000, 144, 21, 0.12, 0, "qualified_reciprocal_saturation"),
    ("mdh", "WT", "oxaloacetate", "NADH", 150, 8.6, 0.4, 203, 5, "qualified_fixed_NADH_55.6x_reciprocal_Km"),
    ("mdh", "WT", "NADH", "oxaloacetate", 25, 2.7, 0.9, 126, 14, "excluded_fixed_oxaloacetate_only_2.9x_Km"),
    ("mdh", "S236A", "oxaloacetate", "NADH", 150, 52, 2, 622, 8, "qualified_fixed_NADH_19.7x_reciprocal_Km"),
    ("mdh", "S236A", "NADH", "oxaloacetate", 25, 7.6, 0.2, 468, 15, "excluded_fixed_oxaloacetate_below_Km"),
    ("mdh", "S236E", "oxaloacetate", "NADH", 150, 19800, 100, 12, 1, "excluded_no_reciprocal_NADH_parameter_for_variant"),
    ("ppase", "WT", "pyrophosphate", "", "", 17, 2, 62, 11, "qualified_single_substrate"),
    ("ppase", "Y140F", "pyrophosphate", "", "", 53, 0, 19, 3, "qualified_single_substrate"),
    ("ppase", "Y140E", "pyrophosphate", "", "", 62, 25, 3.1, 0.4, "qualified_single_substrate"),
)


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
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        if not payload:
            raise ValueError(f"empty download: {url}")
        path.write_bytes(payload)


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines() if line and not line.startswith(">"))


def mutate(sequence: str, position: int, expected: str, replacement: str) -> str:
    if sequence[position - 1] != expected:
        raise ValueError(f"expected {expected}{position}, found {sequence[position - 1]}")
    return sequence[: position - 1] + replacement + sequence[position:]


def sequences() -> dict[str, str]:
    result = {}
    for family, (_, _, accession, _, position, native, control, phosphomimetic, _) in PROTEINS.items():
        wt = read_fasta(RAW / f"{accession}.fasta")
        for variant, replacement in (("WT", native), (f"{native}{position}{control}", control), (f"{native}{position}{phosphomimetic}", phosphomimetic)):
            result[f"{family}-{variant}"] = wt if variant == "WT" else mutate(wt, position, native, replacement)
    return result


def tool_path(path: Path, executable: str) -> str:
    resolved = path.resolve()
    if not executable.startswith("wsl:"):
        return str(resolved)
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}/{resolved.as_posix()[3:]}"


def command(executable: str, *args: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *args]
    _, distribution, binary = executable.split(":", 2)
    return ["wsl", "-d", distribution, binary, *args]


def run_checked(args: list[str]) -> str:
    completed = subprocess.run(args, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {' '.join(args)}\n{completed.stderr}")
    return completed.stdout or completed.stderr


def run_homology(fasta: Path, executable: str, threads: int) -> tuple[Path, Path, str]:
    version = run_checked(command(executable, "version")).strip().splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise RuntimeError(f"MMseqs2 version {version!r} != frozen version {MMSEQS_VERSION!r}")
    homology = SOURCE / "homology"
    if homology.exists():
        if executable.startswith("wsl:"):
            run_checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(homology, executable)])
        else:
            shutil.rmtree(homology)
    homology.mkdir(parents=True)
    combined = homology / "homology_hits.tsv"
    combined.write_text("", encoding="ascii")
    for family in PROTEINS:
        family_dir = homology / family
        family_dir.mkdir()
        family_fasta = family_dir / "queries.fasta"
        entries = read_fasta_entries(fasta)
        with family_fasta.open("w", encoding="ascii", newline="\n") as handle:
            for identifier, sequence in entries.items():
                if identifier.startswith(f"{family}-"):
                    handle.write(f">{identifier}\n{sequence}\n")
        hits = family_dir / "homology_hits.tsv"
        run_checked(command(executable, "easy-search", tool_path(family_fasta, executable), tool_path(REFERENCE, executable), tool_path(hits, executable), tool_path(family_dir / "search-tmp", executable), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits", "--threads", str(threads)))
        with combined.open("a", encoding="ascii", newline="\n") as handle:
            handle.write(hits.read_text(encoding="ascii"))
    prefix = homology / "family-cluster" / "proteins"
    prefix.parent.mkdir()
    run_checked(command(executable, "easy-linclust", tool_path(fasta, executable), tool_path(prefix, executable), tool_path(homology / "cluster-tmp", executable), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--threads", str(threads)))
    return combined, prefix.with_name("proteins_cluster.tsv"), version


def read_fasta_entries(path: Path) -> dict[str, str]:
    result, identifier, parts = {}, None, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if identifier:
                result[identifier] = "".join(parts)
            identifier, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if identifier:
        result[identifier] = "".join(parts)
    return result


def exact_overlap(entries: dict[str, str], records: list[dict[str, object]], corpus: Path) -> dict[str, object]:
    if corpus.stat().st_size != DEVELOPMENT_CORPUS_SIZE or sha256(corpus) != DEVELOPMENT_CORPUS_SHA256:
        raise ValueError("frozen 17,010-row development corpus changed")
    development = json.loads(corpus.read_text(encoding="utf-8"))
    sequence_set = set(entries.values())
    smiles = {value[1] for value in SUBSTRATES.values()}
    exact_sequences = [row for row in development if row.get("Sequence") in sequence_set]
    accepted_pairs = {
        (entries[str(record["sequence_id"])], str(record["substrate_isomeric_smiles"]))
        for record in records
        if record["status_at_normalization"] == "accepted_homology_cold_pool"
    }
    return {
        "development_corpus_sha256": DEVELOPMENT_CORPUS_SHA256,
        "development_rows": len(development),
        "exact_sequence_rows": len(exact_sequences),
        "exact_sequence_substrate_overlap": sum(row.get("Smiles") in smiles for row in exact_sequences),
        "accepted_exact_sequence_substrate_overlap": sum(
            (row.get("Sequence"), row.get("Smiles")) in accepted_pairs for row in exact_sequences
        ),
    }


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(mmseqs: str, threads: int, corpus: Path) -> None:
    missing = [name for name in URLS if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing raw files: {', '.join(missing)}; run with --download")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("frozen UniKP reference changed")
    xml = (RAW / "PMC12667214-fullText.xml").read_text(encoding="utf-8")
    for token in (DOI, "pET‐3a/ttha1671", "pET‐11a/ ttha0536", "Two <italic>k</italic>", "mean ± SD", "1 m<sc>m</sc> CMP and ATP"):
        if token not in xml:
            raise ValueError(f"required primary-source evidence missing: {token}")
    for cid, smiles in SUBSTRATES.values():
        item = json.loads((RAW / f"pubchem-{cid}.json").read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if item["CID"] != cid or item.get("SMILES", item.get("IsomericSMILES")) != smiles:
            raise ValueError(f"PubChem identity changed for CID {cid}")

    pdb_audit = []
    native_entries = sequences()
    for pdb_id, (family, expected_mapping) in PDB_SEQUENCE_CHECKS.items():
        payload = json.loads((RAW / f"rcsb-{pdb_id}-entity-1.json").read_text(encoding="utf-8"))
        observed = payload["entity_poly"]["pdbx_seq_one_letter_code_can"]
        native = native_entries[f"{family}-WT"]
        expected = native if expected_mapping == "exact_full_ORF" else native[1:]
        if observed != expected:
            raise ValueError(f"{pdb_id} polymer does not match expected {family} mapping")
        pdb_audit.append({"pdb_id": pdb_id, "family": family, "mapping": expected_mapping, "polymer_length": len(observed), "polymer_sha256": hashlib.sha256(observed.encode("ascii")).hexdigest()})

    SOURCE.mkdir(parents=True, exist_ok=True)
    entries = native_entries
    fasta = SOURCE / "construct_sequences.fasta"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in entries.items():
            family = identifier.split("-", 1)[0]
            protein = PROTEINS[family]
            handle.write(f">{identifier} {protein[2]} | {protein[3]} | untagged full ORF | {identifier.split('-', 1)[1]}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    hits, cluster, version = run_homology(fasta, mmseqs, threads)
    hit_rows = [line.split("\t") for line in hits.read_text(encoding="ascii").splitlines() if line]
    hit_ids = {row[0] for row in hit_rows}
    records, exclusions = [], []
    for index, row in enumerate(ROWS, 1):
        family, variant, substrate, fixed, fixed_um, km, km_sd, kcat, kcat_sd, disposition = row
        sequence_id = f"{family}-{variant}"
        cid, smiles = SUBSTRATES[substrate]
        qualified = disposition.startswith("qualified_")
        status = "excluded_homology" if qualified and sequence_id in hit_ids else ("accepted_homology_cold_pool" if qualified else "excluded_assay_saturation_or_identity")
        record = {
            "candidate_id": f"ttphos-{index:03d}", "article_doi": DOI, "source_table": "Table 2", "source_row": f"{family} {variant}; varied {substrate}",
            "organism": "Thermus thermophilus HB8", "family": PROTEINS[family][0], "orf_id": PROTEINS[family][1], "uniprot_accession": PROTEINS[family][2],
            "sequence_id": sequence_id, "construct": f"untagged full native ORF from RIKEN {PROTEINS[family][3]}", "variant": variant,
            "variable_substrate": substrate, "substrate_pubchem_cid": cid, "substrate_isomeric_smiles": smiles, "fixed_reaction_substrate": fixed,
            "fixed_reaction_substrate_uM": fixed_um, "endpoint": "kcat_s-1", "kcat_s-1": kcat, "kcat_sd_s-1": kcat_sd, "km_uM": km, "km_sd_uM": km_sd,
            "assay_pH": 7.5, "assay_temperature_C": PROTEINS[family][8], "saturation_disposition": disposition, "status_at_normalization": status,
        }
        records.append(record)
        if status != "accepted_homology_cold_pool":
            exclusions.append({**record, "exclusion_reason": "development_corpus_homology_hit" if status == "excluded_homology" else disposition.removeprefix("excluded_")})
    write_csv(SOURCE / "candidate_records.csv", records)
    write_csv(SOURCE / "exclusions.csv", exclusions)
    write_csv(SOURCE / "direct-kcat-rows.csv", [{**record, "candidate_label_created": record["status_at_normalization"] == "accepted_homology_cold_pool"} for record in records])
    write_csv(SOURCE / "pdb-construct-audit.csv", pdb_audit)

    saturation_rows = [
        {"family": "ampk", "component": "ATP", "concentration": "1 mM", "role": "fixed target-enzyme substrate", "evidence": "authors designate ATP excess; no reciprocal ATP series", "admission_effect": "passes author-reported excess policy"},
        {"family": "ampk", "component": "phosphoenolpyruvate; NADH; pyruvate kinase; lactate dehydrogenase", "concentration": "0.3 mM; 0.15 mM; 20 U/mL; 10 U/mL", "role": "coupled reporter system", "evidence": "fixed concentrations and activities reported; not independently titrated", "admission_effect": "documented reporter capacity; not treated as target-enzyme substrates"},
        {"family": "cmpk", "component": "ATP / CMP", "concentration": "1 mM / 1 mM", "role": "reciprocal fixed target-enzyme substrate", "evidence": "both reciprocal series and variant-specific Km values reported", "admission_effect": "passes reciprocal saturation design"},
        {"family": "cmpk", "component": "phosphoenolpyruvate; NADH; pyruvate kinase; lactate dehydrogenase", "concentration": "0.3 mM; 0.15 mM; 20 U/mL; 10 U/mL", "role": "coupled reporter system", "evidence": "fixed concentrations and activities reported; not independently titrated", "admission_effect": "documented reporter capacity; not treated as target-enzyme substrates"},
        {"family": "umpk", "component": "ATP / reported CMP", "concentration": "1 mM / 1 mM", "role": "reciprocal fixed target-enzyme substrate", "evidence": "ATP is only 4.7x WT or 1.37x S11A reciprocal Km; ATP-series method explicitly says CMP rather than UMP", "admission_effect": "fails all finite rows"},
        {"family": "idh", "component": "NADP+ / isocitrate", "concentration": "4 mM / 2 mM", "role": "reciprocal fixed target-enzyme substrate", "evidence": "both reciprocal series and variant-specific Km values reported", "admission_effect": "passes reciprocal saturation design"},
        {"family": "idh", "component": "MgCl2", "concentration": "5 mM", "role": "metal cofactor", "evidence": "fixed condition reported; no metal titration", "admission_effect": "assay condition recorded, not a varied reaction substrate"},
        {"family": "mdh", "component": "NADH", "concentration": "0.15 mM", "role": "fixed target-enzyme substrate in oxaloacetate series", "evidence": "55.6x WT and 19.7x S236A reciprocal Km; S236E reciprocal Km unavailable", "admission_effect": "WT/S236A pass; S236E fails closed"},
        {"family": "mdh", "component": "oxaloacetate", "concentration": "25 uM", "role": "fixed target-enzyme substrate in NADH series", "evidence": "2.9x WT Km and below S236A Km", "admission_effect": "reciprocal NADH rows fail"},
        {"family": "mdh", "component": "MgCl2", "concentration": "5 mM", "role": "reported assay additive", "evidence": "fixed condition reported; no titration", "admission_effect": "assay condition recorded"},
        {"family": "ppase", "component": "sodium pyrophosphate", "concentration": "0-1.5 mM", "role": "only target-enzyme substrate", "evidence": "single-substrate Michaelis-Menten series", "admission_effect": "passes single-substrate design"},
        {"family": "ppase", "component": "MgCl2", "concentration": "2 mM", "role": "metal cofactor", "evidence": "fixed condition reported; no metal titration", "admission_effect": "assay condition recorded, not a second substrate"},
    ]
    write_csv(SOURCE / "saturation-audit.csv", saturation_rows)

    family_audit = []
    for family in PROTEINS:
        family_rows = [record for record in records if record["sequence_id"].startswith(f"{family}-")]
        family_hits = [row for row in hit_rows if row[0].startswith(f"{family}-")]
        top_hit = max(family_hits, key=lambda row: (float(row[4]), float(row[2]))) if family_hits else None
        family_audit.append({
            "family": family, "reported_finite_kcat_rows": len(family_rows), "assay_qualified_rows": sum(record["saturation_disposition"].startswith("qualified_") for record in family_rows),
            "homology_hit_sequences": len({row[0] for row in family_hits}), "accepted_records": sum(record["status_at_normalization"] == "accepted_homology_cold_pool" for record in family_rows),
            "top_development_hit": None if top_hit is None else {"query": top_hit[0], "target": top_hit[1], "identity": float(top_hit[2]), "query_coverage": float(top_hit[4]), "target_coverage": float(top_hit[5])},
            "family_hit_file": f"homology/{family}/homology_hits.tsv", "family_hit_sha256": sha256(SOURCE / "homology" / family / "homology_hits.tsv"),
        })
    (SOURCE / "family-admission-audit.json").write_text(json.dumps(family_audit, indent=2) + "\n", encoding="ascii")
    overlap = exact_overlap(entries, records, corpus)
    cluster_rows = [line.split("\t") for line in cluster.read_text(encoding="ascii").splitlines() if line]
    accepted = [record for record in records if record["status_at_normalization"] == "accepted_homology_cold_pool"]
    audit = {
        "audited_on": str(date.today()), "source_id": SOURCE_ID, "article_doi": DOI, "method": "six family-independent MMseqs2 easy-search runs and one source-local easy-linclust",
        "mmseqs_version": version, "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0, "development_target_sha256": REFERENCE_SHA256,
        "candidate_records": len(records), "construct_sequences_sha256": sha256(fasta), "homology_hit_sequences": len(hit_ids), "exact_sequence_overlap": sum(float(row[2]) >= 0.999 for row in hit_rows),
        "homology_hits_sha256": sha256(hits), "candidate_mmseqs_families": len({row[0] for row in cluster_rows}), "family_cluster_sha256": sha256(cluster),
        "accepted_catalytic_families": len({record["family"] for record in accepted}),
        "accepted_records": len(accepted), "accepted_unique_sequences": len({record["sequence_id"] for record in accepted}), "accepted_unique_substrates": len({record["substrate_pubchem_cid"] for record in accepted}),
        "exact_sequence_substrate_overlap": overlap["exact_sequence_substrate_overlap"], "accepted_exact_sequence_substrate_overlap": overlap["accepted_exact_sequence_substrate_overlap"], "family_admission_audit_sha256": sha256(SOURCE / "family-admission-audit.json"),
        "readiness_gate_passes": False, "model_predictions_run": False, "claim_boundary": "Family-wise curation and exclusion audit only; no model predictions, inferred labels, or pseudolabels were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    raw_hashes = {path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC12667214", "article_doi": DOI, "article_published": "2025-09-15",
        "article_and_supplement_license": "CC-BY-4.0", "reported_finite_direct_kcat_rows": len(records), "accepted_records": len(accepted), "exact_assayed_construct_sequences": len(entries),
        "construct_mapping": "RIKEN pET-3a/pET-11a plasmids express untagged full HB8 ORFs; no affinity tag or cleavage is reported. UniProt HB8 sequences are mutation-validated and four families are independently corroborated by exact deposited PDB polymers (3CM0, 3AKC, 8YH1, and 2D1C). PDB 2PRD lacks the initiator Met and was not claimed as an exact full construct; no experimental ttMDH structure is reported.",
        "source_scope_correction": "The discovery request named pyruvate kinase, phosphoglycerate mutase, and triosephosphate isomerase. The DOI does not assay those proteins; it assays CMP kinase, UMP kinase, and isocitrate dehydrogenase in addition to AMP kinase, malate dehydrogenase, and pyrophosphatase.",
        "saturation_policy": "Preserve all finite author-reported Table 2 kcat rows, but admit only rows with a defined varied substrate and a demonstrably excess fixed reaction substrate. Reciprocal curves are audited per variant; no value is inferred or refitted.",
        "umpk_source_defect": "The ATP-series method explicitly fixes 1 mM CMP rather than UMP. This was not silently corrected. UMP-series fixed ATP is also only 4.7x WT or 1.37x S11A reciprocal Km; all finite UMPK rows fail closed.",
        "mdh_row_policy": "WT and S236A oxaloacetate-series rows have fixed NADH at 55.6x and 19.7x their reciprocal Km and are assay-qualified. Reciprocal NADH rows use fixed oxaloacetate at only 2.9x or below Km; S236E has no reciprocal NADH parameter.",
        "all_fixed_component_audit": "saturation-audit.csv", "pdb_construct_audit": "pdb-construct-audit.csv", "exact_overlap_audit": overlap, "raw_file_hashes": "raw-file-hashes.json", "artifact_sha256": {name: sha256(SOURCE / name) for name in ("candidate_records.csv", "exclusions.csv", "direct-kcat-rows.csv", "construct_sequences.fasta", "family-admission-audit.json", "saturation-audit.csv", "pdb-construct-audit.csv")},
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
