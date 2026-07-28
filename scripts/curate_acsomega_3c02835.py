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
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts/external/temporal-absolute-kinetics/europepmc-PMC10620931"
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
REFERENCE = ROOT / "artifacts/external/absolute-kinetics-screen/dryad-4964723/homology/unikp_reference.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"

URLS = {
    "PMC10620931-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10620931/fullTextXML",
    "PMC10620931-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10620931/supplementaryFiles?includeInlineImage=false",
    "PMC8671467-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8671467/fullTextXML",
    "PMC8671467-supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8671467/supplementaryFiles?includeInlineImage=false",
    "WP_017726464.1.fasta": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=WP_017726464.1&rettype=fasta&retmode=text",
    "CAA67941.1.fasta": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=CAA67941&rettype=fasta&retmode=text",
    "CAA67941.1.gb": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=CAA67941&rettype=gb&retmode=text",
    "P00760.fasta": "https://rest.uniprot.org/uniprotkb/P00760.fasta",
    "3HGI.fasta": "https://www.rcsb.org/fasta/entry/3HGI/display",
    "1EB2.fasta": "https://www.rcsb.org/fasta/entry/1EB2/display",
    "pubchem-3-methylcatechol-340.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/340/property/IsomericSMILES,Title/JSON",
    "pubchem-L-BApNA-HCl-16219022.json": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/16219022/property/IsomericSMILES,Title/JSON",
}

FAMILIES = {
    "blc23o": "BLC23O-exact-construct.fasta",
    "ro12ctd": "CAA67941.1.fasta",
    "bovine_trypsin": "mature-bovine-trypsin.fasta",
}

BLC23O_PREFIX = "MGSSHHHHHHSSGLVPRGSHM"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acquire() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)
    for archive_name in ("PMC10620931-supplementaryFiles.zip", "PMC8671467-supplementaryFiles.zip"):
        with zipfile.ZipFile(RAW / archive_name) as archive:
            for member in archive.namelist():
                if member.lower().endswith((".pdf", ".docx")):
                    output = RAW / Path(member).name
                    if not output.exists():
                        output.write_bytes(archive.read(member))


def fasta_sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text(encoding="ascii").splitlines() if not line.startswith(">"))


def write_fasta(path: Path, header: str, sequence: str) -> None:
    path.write_text(">" + header + "\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n", encoding="ascii", newline="\n")


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


def checked(args: list[str]) -> str:
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Command failed: {' '.join(args)}\n{result.stderr}")
    return (result.stdout or result.stderr).strip()


def run_mmseqs(executable: str, threads: int) -> str:
    version = checked(command(executable, "version")).splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs version {version!r} does not match frozen {MMSEQS_VERSION!r}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development reference changed")
    for family, query_name in FAMILIES.items():
        family_dir = HOMOLOGY / family
        if family_dir.exists():
            if executable.startswith("wsl:"):
                checked(["wsl", "-d", executable.split(":", 2)[1], "rm", "-rf", "--", tool_path(family_dir, executable)])
            else:
                shutil.rmtree(family_dir)
        family_dir.mkdir(parents=True)
        checked(command(
            executable, "easy-search", tool_path(RAW / query_name, executable),
            tool_path(REFERENCE, executable), tool_path(family_dir / "homology_hits.tsv", executable),
            tool_path(family_dir / "tmp", executable), "--min-seq-id", "0.3", "-c", "0.8",
            "--cov-mode", "0", "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
            "--threads", str(threads),
        ))
    return version


def prepare_queries() -> None:
    native_blc23o = fasta_sequence(RAW / "WP_017726464.1.fasta")
    exact_blc23o = BLC23O_PREFIX + native_blc23o[1:]
    write_fasta(
        RAW / "BLC23O-exact-construct.fasta",
        "BLC23O | pET-28b(+) N-terminal His6-thrombin fusion; NdeI junction; uncleaved",
        exact_blc23o,
    )
    p00760 = fasta_sequence(RAW / "P00760.fasta")
    mature = p00760[p00760.index("IVGGYTCGANTV"):]
    if len(mature) != 223:
        raise ValueError("Expected 223-aa mature bovine trypsin reference")
    write_fasta(RAW / "mature-bovine-trypsin.fasta", "P00760 residues 24-246 | mature reference; commercial lot species unresolved", mature)


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(version: str) -> None:
    rows: list[dict[str, object]] = []
    definitions = [
        ("blc-hepes-76", "BLC23O", "HEPES", 7.6, 32.5, 0.54, 0.02, 0.45, 0.01, 3.0, "3-methylcatechol", 340, "OBLIGATE_O2_SATURATION_UNRESOLVED"),
        ("blc-hepes-74", "BLC23O", "HEPES", 7.4, 32.5, 0.48, 0.01, 0.32, 0.01, 3.0, "3-methylcatechol", 340, "OBLIGATE_O2_SATURATION_UNRESOLVED"),
        ("blc-tris-74", "BLC23O", "Tris-HCl", 7.4, 32.5, 0.71, 0.02, 0.33, 0.00, 3.0, "3-methylcatechol", 340, "VARIED_SUBSTRATE_MAX_BELOW_5X_KM;OBLIGATE_O2_SATURATION_UNRESOLVED"),
        ("blc-phosphate-72", "BLC23O", "Na-phosphate", 7.2, 30.0, 0.76, 0.01, 0.27, 0.00, 3.0, "3-methylcatechol", 340, "VARIED_SUBSTRATE_MAX_BELOW_5X_KM;MN2_MAX_BELOW_5X_KD;OBLIGATE_O2_SATURATION_UNRESOLVED"),
        ("blc-phosphate-74", "BLC23O", "Na-phosphate", 7.4, 32.5, 0.24, 0.01, 0.07, 0.01, 3.0, "3-methylcatechol", 340, "MN2_MAX_BELOW_5X_KD;OBLIGATE_O2_SATURATION_UNRESOLVED"),
        ("ro-hepes", "Ro1,2-CTD", "HEPES", 7.2, "room temperature", 0.00180, 0.00006, 0.64, 0.00, 0.5, "3-methylcatechol", 340, "OBLIGATE_O2_SATURATION_UNRESOLVED;EXACT_ASSAY_CONSTRUCT_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
        ("ro-tris", "Ro1,2-CTD", "Tris-HCl", 7.2, "room temperature", 0.00693, 0.00026, 1.14, 0.01, 0.5, "3-methylcatechol", 340, "OBLIGATE_O2_SATURATION_UNRESOLVED;EXACT_ASSAY_CONSTRUCT_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
        ("ro-phosphate", "Ro1,2-CTD", "Na-phosphate", 7.2, "room temperature", 0.00364, 0.00011, 1.01, 0.01, 0.5, "3-methylcatechol", 340, "OBLIGATE_O2_SATURATION_UNRESOLVED;EXACT_ASSAY_CONSTRUCT_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
        ("trypsin-hepes", "bovine pancreatic trypsin Sigma T1426", "HEPES", 8.0, "room temperature", 3.14, 0.14, 1.51, 0.02, 15.0, "L-BApNA hydrochloride", 16219022, "VARIED_SUBSTRATE_MAX_BELOW_5X_KM;COMMERCIAL_LOT_MOLECULAR_SPECIES_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
        ("trypsin-tris", "bovine pancreatic trypsin Sigma T1426", "Tris-HCl", 8.0, "room temperature", 3.07, 0.16, 1.47, 0.02, 15.0, "L-BApNA hydrochloride", 16219022, "VARIED_SUBSTRATE_MAX_BELOW_5X_KM;COMMERCIAL_LOT_MOLECULAR_SPECIES_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
        ("trypsin-phosphate", "bovine pancreatic trypsin Sigma T1426", "Na-phosphate", 8.0, "room temperature", 2.91, 0.02, 1.53, 0.01, 15.0, "L-BApNA hydrochloride", 16219022, "COMMERCIAL_LOT_MOLECULAR_SPECIES_UNRESOLVED;DEVELOPMENT_HOMOLOGY_HIT"),
    ]
    for item in definitions:
        row_id, enzyme, buffer, ph, temperature, km, km_sd, kcat, kcat_sd, maximum, substrate, cid, codes = item
        rows.append({
            "evidence_id": row_id, "article_doi": "10.1021/acsomega.3c02835",
            "source_table": "Table 1" if row_id.startswith("blc") else ("Table 2" if row_id.startswith("ro") else "Table 3"),
            "enzyme_identity": enzyme, "buffer": buffer, "assay_pH": ph, "assay_temperature_C": temperature,
            "variable_substrate": substrate, "substrate_pubchem_cid": cid, "reported_km_mM": km,
            "reported_km_sd_mM": km_sd, "reported_kcat_s-1": kcat, "reported_kcat_sd_s-1": kcat_sd,
            "maximum_substrate_mM": maximum, "maximum_over_km": round(maximum / float(km), 4),
            "status_at_normalization": "excluded_fail_closed", "exclusion_codes": codes,
            "candidate_label_created": False,
        })
    write_csv(SOURCE / "excluded_records.csv", rows, list(rows[0]))
    candidate_fields = ["candidate_id", "substrate_pubchem_cid", "endpoint", "value", "unit", "status_at_normalization"]
    write_csv(SOURCE / "candidate_records.csv", [], candidate_fields)

    construct = {
        "blc23o": {
            "status": "resolved_exact",
            "sequence_file": "raw/BLC23O-exact-construct.fasta",
            "evidence": "Supplement ORF is 852 nt and translates exactly to WP_017726464.1; cited method gives pET-28b(+), NdeI/XhoI, N-terminal His6 followed by thrombin site; no cleavage step is reported and Ni-NTA-purified protein retains the tag.",
            "junction": "MGSSHHHHHHSSGLVPRGSHM + native residues 2-283",
        },
        "ro12ctd": {
            "status": "unresolved_fail_closed",
            "native_accession": "CAA67941.1",
            "experimental_structure": "PDB 3HGI chain A exactly matches all 280 native residues",
            "blocker": "The article gives pQE-80L, BamHI/HindIII and an N-terminal His tag, but deposits neither the synthesized coding sequence nor the insert-vector junction. Whether the synthetic insert retained its initiator and exact pQE-80L junction residues cannot be established.",
        },
        "trypsin": {
            "status": "unresolved_fail_closed",
            "product": "Sigma T1426, TPCK-treated bovine pancreatic trypsin",
            "reference": "P00760 residues 24-246 and PDB 1EB2 chain A are the same 223-aa mature trypsin sequence",
            "blocker": "No lot, chain composition, autolysis state, isoenzyme proportion, or lot-specific sequence is reported for the commercial TPCK-treated powder; a canonical mature chain is not substituted for the assayed material.",
        },
        "model_predictions_run": False,
    }
    (SOURCE / "construct-resolution.json").write_text(json.dumps(construct, indent=2) + "\n", encoding="ascii")

    homology = {}
    for family, query_name in FAMILIES.items():
        hits = HOMOLOGY / family / "homology_hits.tsv"
        lines = [line for line in hits.read_text(encoding="utf-8").splitlines() if line]
        homology[family] = {
            "query": f"raw/{query_name}", "query_sha256": sha256(RAW / query_name), "hits": len(lines),
            "hits_sha256": sha256(hits), "disposition": "development_homology_hit" if lines else "homology_cold",
        }
    audit = {
        "audited_on": "2026-07-27", "source_id": "europepmc-PMC10620931",
        "article_doi": "10.1021/acsomega.3c02835", "method": "MMseqs2 easy-search, family-by-family",
        "mmseqs_version": version, "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256, "families": homology,
        "candidate_records": {"count": 0, "sha256": sha256(SOURCE / "candidate_records.csv")},
        "reported_rows_retained_as_exclusions": 11, "accepted_records": 0,
        "claim_boundary": "Prescreens on an exact BLC23O construct and source-supported native/reference sequences; unresolved sequences are sufficient only to exclude, never to admit.",
        "model_predictions_run": False,
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")

    saturation = {
        "policy": "Varied-substrate maximum must reach at least 5x fitted Km; obligate cosubstrates must have auditable saturation evidence.",
        "rows": [{"evidence_id": row["evidence_id"], "maximum_over_km": row["maximum_over_km"], "varied_substrate_decision": "pass" if row["maximum_over_km"] >= 5 else "fail", "oxygen_decision": "unresolved" if row["evidence_id"].startswith(("blc", "ro")) else "not_applicable"} for row in rows],
        "dioxygenase_blocker": "The assays were performed in microplates without reported dissolved O2, gas equilibration, sealed/open geometry control, or O2 titration. O2 is an obligate cosubstrate, so eight donor-substrate fits are apparent values only.",
        "blc23o_manganese_note": "Using the same 5x rule diagnostically, 10 uM Mn2+ is 5.59-6.71x the matching HEPES/Tris Kd and passes, whereas 100 uM is only 2.26x and 1.81x the phosphate Kd values and fails. O2 independently blocks every BLC23O row.",
        "accepted_records": 0,
    }
    (SOURCE / "saturation-audit.json").write_text(json.dumps(saturation, indent=2) + "\n", encoding="ascii")

    raw_hashes = {path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path), "url": URLS.get(path.name, "derived reproducibly by this script")} for path in sorted(RAW.iterdir()) if path.is_file()}
    (SOURCE / "raw-file-hashes.json").write_text(json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii")
    provenance = {
        "source_id": "europepmc-PMC10620931", "article_doi": "10.1021/acsomega.3c02835",
        "stable_record_url": "https://europepmc.org/articles/PMC10620931", "article_published": "2023-10-18",
        "license": "CC-BY-4.0", "reported_direct_kcat_rows": 11, "accepted_records": 0,
        "row_scope": {"BLC23O": 5, "Ro1,2-CTD": 3, "bovine_trypsin": 3},
        "pubchem": {"3-methylcatechol": {"cid": 340}, "L-BApNA hydrochloride": {"cid": 16219022}},
        "endpoint_decision": "All values are direct fitted kcat evidence under distinct buffer/pH conditions, not 11 independent protein-substrate pairs. They are retained row-by-row as excluded evidence.",
        "final_disposition": "excluded_fail_closed", "candidate_labels_created": 0,
        "blockers": ["Dioxygenase O2 saturation unresolved (8 rows)", "Varied-substrate range below 5x Km (4 rows)", "Exact Ro1,2-CTD assay junction unresolved (3 rows)", "Commercial trypsin lot molecular species unresolved (3 rows)", "Frozen development homology for Ro1,2-CTD and trypsin (6 rows)"],
        "model_predictions_run": False, "raw_file_hashes": "raw-file-hashes.json",
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")


def validate_raw() -> None:
    missing = [name for name in URLS if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw files: {missing}; run with --acquire")
    article = (RAW / "PMC10620931-fullText.xml").read_text(encoding="utf-8")
    for marker in ("10.1021/acsomega.3c02835", "CAA67941", "pQE80L", "T1426", "0.05 and 3 mM", "4 to 500", "0.1 to 15 mM"):
        if marker not in article:
            raise ValueError(f"Primary article marker absent: {marker}")
    expected = {"WP_017726464.1.fasta": 283, "CAA67941.1.fasta": 280}
    for name, length in expected.items():
        if len(fasta_sequence(RAW / name)) != length:
            raise ValueError(f"Unexpected sequence length for {name}")
    supplement_docx = next(RAW.glob("*MOESM1_ESM.docx"), None)
    if supplement_docx is None:
        raise FileNotFoundError("BLC23O supplement DOCX was not extracted")
    with zipfile.ZipFile(supplement_docx) as archive:
        document = archive.read("word/document.xml").decode("utf-8")
    document_text = re.sub(r"<[^>]+>", " ", document)
    dna_runs = [re.sub(r"\s", "", match) for match in re.findall(r"[ATGC](?:\s*[ATGC]){29,}", document_text)]
    if len(dna_runs) != 1 or len(dna_runs[0]) != 852:
        raise ValueError("Expected one deposited 852-nt BLC23O ORF in supplement")
    codons = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }
    translated = "".join(codons[dna_runs[0][index:index + 3]] for index in range(0, 852, 3))
    if translated != fasta_sequence(RAW / "WP_017726464.1.fasta") + "*":
        raise ValueError("Deposited BLC23O ORF does not translate to WP_017726464.1")
    compounds = {
        "pubchem-3-methylcatechol-340.json": (340, "CC1=C(C(=CC=C1)O)O"),
        "pubchem-L-BApNA-HCl-16219022.json": (16219022, "C1=CC=C(C=C1)C(=O)N[C@@H](CCCN=C(N)N)C(=O)NC2=CC=C(C=C2)[N+](=O)[O-].Cl"),
    }
    for name, (cid, smiles) in compounds.items():
        record = json.loads((RAW / name).read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        if record["CID"] != cid or record["SMILES"] != smiles:
            raise ValueError(f"PubChem mapping changed for {name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquire", action="store_true")
    parser.add_argument("--run-mmseqs", action="store_true")
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.acquire:
        acquire()
    validate_raw()
    prepare_queries()
    version = run_mmseqs(args.mmseqs, args.threads) if args.run_mmseqs else MMSEQS_VERSION
    missing_hits = [family for family in FAMILIES if not (HOMOLOGY / family / "homology_hits.tsv").is_file()]
    if missing_hits:
        raise FileNotFoundError(f"Pinned homology outputs missing for {missing_hits}; run with --run-mmseqs")
    write_outputs(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
