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
SOURCE_ID = "europepmc-PMC12894115"
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / SOURCE_ID
RAW = SOURCE / "raw"
HOMOLOGY = SOURCE / "homology"
DOI = "10.1007/s12010-025-05449-0"
REFERENCE = ROOT / "artifacts" / "homology-final" / "unique_proteins.fasta"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
MMSEQS = "mmseqs"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"

ORF_SEQUENCE = (
    "MLSSVVIVGGGTSGWMTAAYLRAALGTSVNVTVIESKRIKTIGVGEATFSTVRHFFDYLGLSENDWMPKCNATY"
    "KLGIRFENWRAKGHYFYHPFERLRVVDGFPLTEWWLNKKPSDRFDQDVFLMSEICDTMRSPRYLDGTLFEQDFVE"
    "HGGSMDPERSTLSEQATQFPYAYQFDASLLADFLTEYATTRGARHIEDDVVEVVRDERGWISHLKTREHGELAGDL"
    "FVDCTGFAGLLLNKTLGEPFVSYQNTLPNDSAVALRVPHDAERTRLRPCTTATAQEAGWIWTIPLFERIGTGYVYA"
    "SDYTTPEEAERTLREFVGPQAADVEANHIRMRIGRSRHSWVNNCVAIGLSSGFVEPLESTGIFFIQQGIEELVKHF"
    "PDAKWDPKLRDSYNRVVANTMDGVREFLVLHYRTAARNDNAYWRDAKTRELPDGLAARLEAWQSKLPTEETVFPHY"
    "HGFEPYSYHAMLLGLGGLDVKPAPVLAHMDDSRAAQEIQRLKDQARDIAKRLPSQYEYLAQMH"
)
STRUCTURES = {
    "L-tryptophan": (6305, "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N"),
    "6-chloro-L-tryptophan": (10062693, "C1=CC2=C(C=C1Cl)NC=C2C[C@@H](C(=O)O)N"),
    "7-chloro-L-tryptophan": (3081936, "C1=CC2=C(C(=C1)Cl)NC=C2C[C@@H](C(=O)O)N"),
    "5-chloro-L-tryptophan": (644330, "C1=CC2=C(C=C1Cl)C(=CN2)C[C@@H](C(=O)O)N"),
}
URLS = {
    "PMC12894115-fullText.xml": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12894115/fullTextXML"
    ),
    "PMC12894115-supplementaryFiles.zip": (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12894115/supplementaryFiles"
    ),
    "PNE38568.1.fasta": (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=protein&id=PNE38568.1&rettype=fasta&retmode=text"
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
    supplement = RAW / "12010_2025_5449_MOESM1_ESM.docx"
    if not supplement.is_file():
        with zipfile.ZipFile(RAW / "PMC12894115-supplementaryFiles.zip") as archive:
            with archive.open("12010_2025_5449_MOESM1_ESM.docx") as source, supplement.open("wb") as output:
                shutil.copyfileobj(source, output)
    for name, (cid, _) in STRUCTURES.items():
        path = RAW / f"pubchem-{cid}.json"
        if path.is_file() and path.stat().st_size:
            continue
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/"
            "property/IsomericSMILES,Title/JSON"
        )
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as output:
            shutil.copyfileobj(response, output)


def validate_raw() -> None:
    xml = (RAW / "PMC12894115-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI, "creativecommons.org/licenses/by/4.0/", "9.6 mM NADH", "100 µM FAD",
        "10 mM NaCl", "10 µM Fre", "5.5 µM SnFDHal", "20–1000 µM",
        "243.09", "91.273", "22.092", "Nde1 and Xho1",
    )
    missing = [marker for marker in required if marker not in xml]
    if missing:
        raise ValueError(f"Article XML lacks expected evidence: {missing}")
    with zipfile.ZipFile(RAW / "12010_2025_5449_MOESM1_ESM.docx") as archive:
        document = archive.read("word/document.xml").decode("utf-8")
    for marker in (ORF_SEQUENCE[:30], ORF_SEQUENCE[-30:], "GAATTCTCAGTGCATCTGGGCGAGG"):
        if marker not in document.replace("</w:t><w:t>", ""):
            raise ValueError(f"Supplement lacks expected sequence evidence: {marker}")
    for name, (cid, smiles) in STRUCTURES.items():
        compound = json.loads((RAW / f"pubchem-{cid}.json").read_text(encoding="utf-8"))[
            "PropertyTable"
        ]["Properties"][0]
        if compound["CID"] != cid or compound["SMILES"] != smiles:
            raise ValueError(f"Unexpected PubChem mapping for {name}")


def wsl_path(path: Path) -> str:
    resolved = path.resolve()
    return f"/mnt/{resolved.drive[0].lower()}/{resolved.as_posix()[3:]}"


def mmseqs_command(*arguments: str) -> list[str]:
    _, distribution, executable = MMSEQS.split(":", 2)
    return ["wsl", "-d", distribution, executable, *arguments]


def reset_wsl_workspace(path: Path) -> None:
    if path.exists():
        _, distribution, _ = MMSEQS.split(":", 2)
        subprocess.run(["wsl", "-d", distribution, "rm", "-rf", "--", wsl_path(path)], check=True)


def run_mmseqs(query: Path) -> str:
    if not REFERENCE.is_file() or sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen development FASTA is missing or has changed")
    version = subprocess.run(
        mmseqs_command("version"), check=True, capture_output=True, text=True
    ).stdout.strip()
    if version != MMSEQS_VERSION:
        raise ValueError(f"MMseqs2 version {version!r} is not the frozen version")
    HOMOLOGY.mkdir(parents=True, exist_ok=True)
    hits = HOMOLOGY / "homology_hits.tsv"
    workspace = HOMOLOGY / "search-tmp"
    reset_wsl_workspace(workspace)
    subprocess.run(
        mmseqs_command(
            "easy-search", wsl_path(query), wsl_path(REFERENCE), wsl_path(hits),
            wsl_path(workspace), "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0",
            "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        ),
        check=True,
    )
    return version


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def curate(*, acquire: bool, run_homology: bool) -> None:
    if acquire:
        download_raw()
    validate_raw()
    SOURCE.mkdir(parents=True, exist_ok=True)
    query = SOURCE / "disclosed_orf_sequence.fasta"
    query.write_text(
        ">SnFDHal-disclosed-orf SI_Table_1; not the unresolved His6-tagged assay construct\n"
        + "\n".join(ORF_SEQUENCE[i:i + 80] for i in range(0, len(ORF_SEQUENCE), 80)) + "\n",
        encoding="ascii", newline="\n",
    )
    version = run_mmseqs(query) if run_homology else MMSEQS_VERSION
    hits = HOMOLOGY / "homology_hits.tsv"
    if not hits.is_file():
        raise RuntimeError("Pinned MMseqs evidence is required")
    hit_lines = [line for line in hits.read_text(encoding="utf-8").splitlines() if line]

    candidate_fields = [
        "candidate_id", "article_doi", "stable_record_url", "source_file", "source_row",
        "organism", "enzyme_identity", "sequence_accession", "sequence_id", "construct",
        "variable_substrate", "substrate_pubchem_cid", "substrate_isomeric_smiles",
        "endpoint", "kcat_s-1", "assay_pH", "assay_temperature_C", "status_at_normalization",
    ]
    write_csv(SOURCE / "candidate_records.csv", [], candidate_fields)

    exclusions = []
    reported = (
        ("L-tryptophan", 0.6216, 243.09),
        ("6-chloro-L-tryptophan", 0.9451, 91.273),
        ("7-chloro-L-tryptophan", 0.1718, 22.092),
    )
    for substrate, km, kcat_min in reported:
        cid, smiles = STRUCTURES[substrate]
        vmax_mm_min = kcat_min * 0.0055
        rate_at_1mm = vmax_mm_min * 1.0 / (km + 1.0)
        exclusions.append({
            "source_row": f"Table 1; {substrate}", "substrate": substrate,
            "substrate_pubchem_cid": cid, "substrate_isomeric_smiles": smiles,
            "reported_km_mM": km, "reported_kcat_min-1": kcat_min,
            "converted_kcat_s-1": kcat_min / 60,
            "implied_vmax_at_5.5uM_enzyme_mM_min-1": round(vmax_mm_min, 6),
            "implied_product_at_1mM_after_30min_mM": round(rate_at_1mm * 30, 6),
            "exclusion_reason": (
                "Not a direct SnFDHal turnover label: 30-min endpoint coupled assay; fixed Fre/FAD/"
                "NADH/chloride and oxygen saturation not established; exact tagged construct unresolved; "
                "reported kcat is mass-balance-incompatible with the stated endpoint assay."
            ),
            "candidate_label_created": False,
        })
    write_csv(SOURCE / "exclusions.csv", exclusions, list(exclusions[0]))

    raw_hashes = {
        path.name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(RAW.iterdir()) if path.is_file()
    }
    (SOURCE / "raw-file-hashes.json").write_text(
        json.dumps(raw_hashes, indent=2) + "\n", encoding="ascii"
    )
    final_disposition = "excluded_homology_and_blocked_no_valid_direct_kcat_labels"
    blocker = {
        "audited_on": "2026-07-27", "source_id": SOURCE_ID, "article_doi": DOI,
        "disposition": final_disposition, "candidate_labels_created": 0,
        "blockers": [
            {
                "gate": "frozen_development_homology",
                "evidence": (
                    f"Pinned MMseqs2 found {len(hit_lines)} development-corpus hits for the disclosed "
                    "516-aa ORF at 30% identity, 80% coverage, coverage mode 0. The best hit is 39.0% "
                    "identity with 99.8% query coverage."
                ),
            },
            {
                "gate": "exact_assay_construct",
                "evidence": (
                    "SI Table 1 gives a 516-aa ORF, while the assay used a His6-tagged protein. The tag "
                    "sequence/terminus and retained vector residues are not disclosed. Methods say NdeI/"
                    "XhoI, but SI Table 2 provides an EcoRI reverse primer (GAATTC), preventing exact "
                    "reconstruction of the expressed pET28a product."
                ),
            },
            {
                "gate": "stable_exact_sequence_accession",
                "evidence": (
                    "No accession is assigned to the SI sequence. PNE38568.1/A0A2N8PC14 is a distinct "
                    "501-aa tryptophan halogenase from S. noursei JCM 4701 and is not the 516-aa SI ORF."
                ),
            },
            {
                "gate": "direct_turnover_and_initial_rate",
                "evidence": (
                    "Product was measured after 30 min, not over a demonstrated initial-rate interval. "
                    "At 5.5 uM enzyme, each reported kcat implies more product in 30 min at 1 mM substrate "
                    "than the reaction initially contained (see exclusions.csv), so the endpoint method "
                    "cannot support the reported values as direct turnover labels."
                ),
            },
            {
                "gate": "coupled_component_saturation",
                "evidence": (
                    "The general assay contains E. coli Fre (10 uM), FAD (100 uM), NADH (9.6 mM), and "
                    "NaCl (10 mM), but no concentration-dependence or saturation evidence is reported for "
                    "Fre, FAD/reduced-flavin supply, NADH, or chloride. Kinetic Methods only says "
                    "'halogenation buffer' and does not independently enumerate these fixed components."
                ),
            },
            {
                "gate": "oxygen_saturation",
                "evidence": (
                    "Molecular oxygen is an obligatory substrate for hydroperoxyflavin formation. Reactions "
                    "were placed in a static incubator; dissolved oxygen, headspace/mixing, oxygen dependence, "
                    "and oxygen saturation were not reported."
                ),
            },
        ],
        "claim_boundary": "Audit and homology prescreen only; no reported value was normalized as kcat and no model prediction was generated.",
    }
    (SOURCE / "blocker-evidence.json").write_text(
        json.dumps(blocker, indent=2) + "\n", encoding="ascii"
    )
    provenance = {
        "source_id": SOURCE_ID, "stable_record_url": "https://europepmc.org/articles/PMC12894115",
        "article_doi": DOI, "article_published": "2025-11-08", "license": "CC-BY-4.0",
        "reported_kinetic_rows": 3, "curated_direct_kcat_rows": 0,
        "kinetics_source": "raw/PMC12894115-fullText.xml, Kinetic Analysis and Table 1",
        "sequence_source": "raw/12010_2025_5449_MOESM1_ESM.docx, SI Table 1",
        "disclosed_orf": {"length_aa": len(ORF_SEQUENCE), "sha256": hashlib.sha256(ORF_SEQUENCE.encode()).hexdigest()},
        "sequence_accession_audit": {
            "exact_accession": None, "nearby_annotation": "PNE38568.1 / UniProt A0A2N8PC14",
            "nearby_annotation_is_exact": False,
            "reason": "Nearby record is 501 aa, belongs to strain JCM 4701, and differs from the 516-aa SI sequence.",
        },
        "assay_construct": {
            "reported": "purified His6-tagged SnFDHal from pET28a-SnFDHal in E. coli BL21(DE3)",
            "recoverable_exact_sequence": False,
            "conflict": "Methods NdeI/XhoI cloning versus SI EcoRI reverse primer; tag/vector peptide unspecified",
        },
        "assay": {
            "volume_uL": 100, "buffer": "100 mM phosphate, pH 6.0", "temperature_C": 35,
            "duration_min": 30, "substrate_range_uM": "20-1000", "SnFDHal_uM": 5.5,
            "Fre_uM": 10, "FAD_uM": 100, "NADH_mM": 9.6, "NaCl_mM": 10,
            "incubation": "static incubator", "quantification": "quenched endpoint HPLC product formation",
            "fit": "Michaelis-Menten nonlinear regression, GraphPad Prism 10.2.2",
        },
        "auxiliary_reductase": {
            "identity": "E. coli BL21(DE3) flavin reductase Fre", "exact_construct_or_accession": None,
            "purification_reference": "Zeng and Zhan 2010, DOI 10.1002/cbic.201000439",
            "saturation_established": False,
        },
        "fixed_component_saturation": {"Fre": False, "FAD": False, "NADH": False, "chloride": False, "oxygen": False},
        "structure_mapping": {name: {"pubchem_cid": cid, "isomeric_smiles": smiles} for name, (cid, smiles) in STRUCTURES.items()},
        "product_identity": {
            "L-tryptophan_product": "5-chloro-L-tryptophan confirmed by authentic-standard HPLC, MS isotope pattern, and 1H NMR",
            "chlorotryptophan_products": "putative 6,7-dichlorotryptophan based on retention time/UV comparison; not used as varied-substrate identities",
        },
        "raw_file_hashes": "raw-file-hashes.json", "homology_prescreen_sequence": "disclosed ORF only, not exact tagged construct",
        "model_predictions_run": False, "final_disposition": final_disposition,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")
    audit = {
        "audited_on": "2026-07-27", "source_id": SOURCE_ID,
        "method": "MMseqs2 easy-search", "mmseqs_version": version,
        "wsl_distribution": "Ubuntu-24.04", "min_identity": 0.3, "coverage": 0.8, "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256,
        "query_scope": "516-aa disclosed ORF; assay His6 tag/vector residues unresolved",
        "query_sha256": sha256(query), "homology_hit_alignments": len(hit_lines),
        "homology_hits_sha256": sha256(hits), "candidate_records": {"count": 0, "sha256": sha256(SOURCE / "candidate_records.csv")},
        "status": "prescreen_complete_but_label_and_construct_blocked",
        "claim_boundary": "Homology result cannot cure the direct-label or exact-construct blockers; no prediction was run.",
    }
    (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-acquire", action="store_true")
    parser.add_argument("--skip-mmseqs", action="store_true")
    args = parser.parse_args()
    curate(acquire=not args.no_acquire, run_homology=not args.skip_mmseqs)
