from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "external"
    / "temporal-absolute-kinetics"
    / "biostudies-S-EPMC10589216"
)
RAW = SOURCE / "raw"
DOI = "10.1038/s41467-023-42456-3"
SOURCE_ID = "biostudies-S-EPMC10589216"
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
TAG = "WSHPQFEK"

PDB_FORMS = {
    "8JDC": ("WT", "FAD"),
    "8JDD": ("H405A", "FAD"),
    "8JDE": ("WT", "FAD + Mn2+ + D-lactate"),
    "8JDF": ("H405A", "FAD + D-lactate"),
    "8JDG": ("H405A", "FAD + D-2-hydroxybutyrate"),
    "8JDN": ("H405A", "FAD + D-2-hydroxyvalerate"),
    "8JDO": ("H405A", "FAD + D-2-hydroxyhexanoate"),
    "8JDB": ("H405A", "FAD + D-2-hydroxyoctanoate"),
    "8JDP": ("H405A", "FAD + D-2-hydroxyisovalerate"),
    "8JDQ": ("H405A", "FAD + D-2-hydroxyisocaproate"),
    "8JDR": ("H405A", "FAD + (2R,3S)-2-hydroxy-3-methylvalerate"),
    "8JDS": ("WT", "FAD + Mn2+ + pyruvate"),
    "8JDT": ("WT", "FAD + Mn2+ + 2-ketobutyrate"),
    "8JDU": ("WT", "FAD + Mn2+ + 2-ketovalerate"),
    "8JDV": ("WT", "FAD + Mn2+ + 2-ketohexanoate"),
    "8JDX": ("WT", "FAD + Mn2+ + 2-ketoisovalerate"),
    "8JDY": ("WT", "FAD + Mn2+ + 2-ketoisocaproate"),
    "8JDZ": ("WT", "FAD + Mn2+ + 2-keto-3-methylvalerate"),
}

SUBSTRATES = {
    "D-lactate": (61503, "C[C@H](C(=O)O)O"),
    "D-2-hydroxybutyrate": (449265, "CC[C@H](C(=O)O)O"),
    "D-2-hydroxyvalerate": (6950705, "CCC[C@H](C(=O)O)O"),
    "D-2-hydroxyhexanoate": (6995511, "CCCC[C@H](C(=O)O)O"),
    "D-2-hydroxyoctanoate": (5312860, "CCCCCC[C@H](C(=O)O)O"),
    "D-2-hydroxyisovalerate": (5289545, "CC(C)[C@H](C(=O)O)O"),
    "D-2-hydroxyisocaproate": (439960, "CC(C)C[C@H](C(=O)O)O"),
    "(2R,3S)-2-hydroxy-3-methylvalerate": (
        10820562,
        "CC[C@H](C)[C@H](C(=O)O)O",
    ),
    "D-2-hydroxy-3-phenylpropionate": (
        643327,
        "C1=CC=C(C=C1)C[C@H](C(=O)O)O",
    ),
}

# Substrate, source label, product, Vmax +/- SEM, Km +/- SEM (uM), kcat +/- SEM (min-1), efficiency.
KINETICS = (
    ("D-lactate", "D-lactate (D-LAC)", "pyruvate", 1.22, 0.01, 122.0, 4.5, 61.0, 0.6, 0.50),
    ("D-2-hydroxybutyrate", "D-2-hydroxybutyrate (D-2-HB)", "2-ketobutyrate", 1.60, 0.02, 62.9, 2.8, 80.0, 0.9, 1.27),
    ("D-2-hydroxyvalerate", "DL-2-hydroxyvalerate (DL-2-HV)", "2-ketovalerate", 1.12, 0.02, 31.4, 2.7, 56.0, 1.1, 1.78),
    ("D-2-hydroxyhexanoate", "DL-2-hydroxyhexanoate (DL-2-HH)", "2-ketohexanoate", 1.48, 0.02, 34.9, 1.6, 74.0, 0.7, 2.12),
    ("D-2-hydroxyoctanoate", "DL-2-hydroxyoctanoate (DL-2-HO)", "2-ketooctanoate", 0.409, 0.010, 104.0, 9.0, 20.5, 0.5, 0.197),
    ("D-2-hydroxyisovalerate", "D-2-hydroxyisovalerate (D-2-HIV)", "2-ketoisovalerate", 1.11, 0.02, 13.9, 1.2, 55.5, 0.9, 3.99),
    ("D-2-hydroxyisocaproate", "DL-2-hydroxyisocaproate (DL-2-HIC)", "2-ketoisocaproate", 0.902, 0.017, 16.4, 0.7, 45.1, 0.8, 2.75),
    ("(2R,3S)-2-hydroxy-3-methylvalerate", "DL-2-hydroxy-3-methyl-valerate (DL-2-HMV)", "2-keto-3-methylvalerate", 1.01, 0.01, 16.5, 0.1, 50.5, 0.6, 3.06),
    ("D-2-hydroxy-3-phenylpropionate", "D-2-hydroxy-3-phenyl-propionate (D-2-HPP)", "2-keto-3-phenylpropionate", 0.788, 0.037, 337.0, 44.0, 39.4, 1.8, 0.117),
)

VARIANTS = {
    "mldhd-wt": (),
    "mldhd-r347a": ((347, "R", "A"),),
    "mldhd-h398a": ((398, "H", "A"),),
    "mldhd-h405a": ((405, "H", "A"),),
    "mldhd-e442a": ((442, "E", "A"),),
    "mldhd-h443a": ((443, "H", "A"),),
    "mldhd-t228m": ((228, "T", "M"),),
    "mldhd-r347w": ((347, "R", "W"),),
    "mldhd-w351c": ((351, "W", "C"),),
    "mldhd-t440m": ((440, "T", "M"),),
}


def downloads() -> dict[str, str]:
    result = {
        "PMC10589216-fullText.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10589216/fullTextXML",
        "Q7TNG8.fasta": "https://rest.uniprot.org/uniprotkb/Q7TNG8.fasta",
    }
    base = "https://www.ebi.ac.uk/biostudies/files/S-EPMC10589216"
    for number in range(1, 7):
        extension = "pdf" if number <= 4 else "xlsx"
        name = f"41467_2023_42456_MOESM{number}_ESM.{extension}"
        result[name] = f"{base}/{name}"
    for pdb_id in PDB_FORMS:
        result[f"{pdb_id}.cif"] = f"https://files.rcsb.org/download/{pdb_id}.cif"
    for name, (cid, _) in SUBSTRATES.items():
        slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
        result[f"pubchem-{slug}.json"] = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/property/IsomericSMILES/JSON"
        )
    return result


def download_raw() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in downloads().items():
        path = RAW / name
        if path.is_file() and path.stat().st_size:
            continue
        request = urllib.request.Request(url, headers={"User-Agent": "BioCandidateRanker/curation"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    payload = response.read()
                if not payload:
                    raise ValueError(f"Empty response for {url}")
                path.write_bytes(payload)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )


def cif_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s+\n?;([^;]+);", text, re.DOTALL
    )
    if not match:
        raise ValueError(f"Canonical polymer sequence missing from {path.name}")
    return re.sub(r"\s+", "", match.group(1))


def mutate_native(native: str, changes: tuple[tuple[int, str, str], ...]) -> str:
    residues = list(native)
    for position, expected, replacement in changes:
        if residues[position - 1] != expected:
            raise ValueError(
                f"Q7TNG8 expected {expected}{position}, found {residues[position - 1]}"
            )
        residues[position - 1] = replacement
    return TAG + "".join(residues[21:])


def build_sequences() -> dict[str, tuple[str, str]]:
    native = read_fasta(RAW / "Q7TNG8.fasta")
    if len(native) != 484 or not native.startswith("MAMLLRVATQRLSPWRSFCSR"):
        raise ValueError("Unexpected UniProt Q7TNG8 sequence")
    sequences = {
        sequence_id: (sequence_id.removeprefix("mldhd-").upper(), mutate_native(native, changes))
        for sequence_id, changes in VARIANTS.items()
    }
    for pdb_id, (variant, _) in PDB_FORMS.items():
        expected = sequences[f"mldhd-{variant.lower()}"][1]
        observed = cif_sequence(RAW / f"{pdb_id}.cif")
        if observed != expected:
            raise ValueError(f"{pdb_id} does not match exact {variant} assayed construct")
    return sequences


def verify_raw() -> None:
    missing = [name for name in downloads() if not (RAW / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing raw files: {', '.join(missing)}; run with --download")
    xml = (RAW / "PMC10589216-fullText.xml").read_text(encoding="utf-8")
    required = (
        DOI,
        "200\u2009\u03bcM PMS",
        "100\u2009\u03bcM DCIP",
        "molar ratio of FAD:mLDHD",
        "0.73\u2009\u00b1\u20090.04",
        "residues 22\u2013484",
        "Strep tag at the N-terminus",
    )
    for text in required:
        if text not in xml:
            raise ValueError(f"Required article evidence not found: {text}")
    for name, (cid, smiles) in SUBSTRATES.items():
        slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
        payload = json.loads((RAW / f"pubchem-{slug}.json").read_text(encoding="utf-8"))
        compound = payload["PropertyTable"]["Properties"][0]
        if compound["CID"] != cid or compound["SMILES"] != smiles:
            raise ValueError(f"Unexpected PubChem identity for {name}")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Frozen UniKP development reference hash changed")


FIELDS = [
    "candidate_id", "article_doi", "source_table", "source_row", "organism",
    "sequence_id", "construct", "variable_substrate", "substrate_pubchem_cid",
    "substrate_isomeric_smiles", "endpoint", "kcat_s-1", "kcat_sem_s-1",
    "km_uM", "km_sem_uM", "assay_pH", "assay_temperature_C",
    "status_at_normalization",
]


def exclusion_records() -> list[dict[str, object]]:
    records = []
    for index, row in enumerate(KINETICS, 1):
        substrate, source_label, product, vmax, vmax_sem, km, km_sem, kcat, kcat_sem, efficiency = row
        cid, smiles = SUBSTRATES[substrate]
        records.append(
            {
                "source_row": index,
                "article_doi": DOI,
                "source_table": "Table 1",
                "source_substrate_label": source_label,
                "defined_D_substrate": substrate,
                "substrate_pubchem_cid": cid,
                "substrate_isomeric_smiles": smiles,
                "product": product,
                "sequence_id": "mldhd-wt",
                "construct": "N-terminal WSHPQFEK Strep tag + Q7TNG8 residues 22-484",
                "source_vmax_umol_min-1_mg-1": vmax,
                "source_vmax_sem_umol_min-1_mg-1": vmax_sem,
                "km_uM": km,
                "km_sem_uM": km_sem,
                "source_kcat_min-1": kcat,
                "source_kcat_sem_min-1": kcat_sem,
                "normalized_kcat_s-1": kcat / 60,
                "normalized_kcat_sem_s-1": kcat_sem / 60,
                "source_kcat_per_km_min-1_uM-1": efficiency,
                "substrate_range_mM": "0-2 (D-isomer concentration; racemate doubled where applicable)",
                "prosthetic_FAD_occupancy_mol_per_mol": 0.73,
                "prosthetic_FAD_occupancy_sem": 0.04,
                "activity_corrected_for_FAD_occupancy": True,
                "exogenous_FAD_added": False,
                "primary_electron_acceptor": "phenazine methosulfate (PMS)",
                "primary_electron_acceptor_uM": 200,
                "terminal_indicator_acceptor": "2,6-dichloroindophenol (DCIP)",
                "terminal_indicator_acceptor_uM": 100,
                "electron_acceptor_saturation_tested": False,
                "assay_buffer": "50 mM Tris-HCl",
                "assay_pH": 7.4,
                "assay_temperature_C": 37,
                "assay_volume_uL": 200,
                "enzyme_amount_ug": 2,
                "manganese_uM": 50,
                "monitoring_wavelength_nm": 600,
                "fit": "Michaelis-Menten, GraphPad Prism 7.0",
                "replicates": 3,
                "error_type": "SEM",
                "exclusion_reason": "coupled_apparent_kcat_electron_acceptor_saturation_not_demonstrated",
                "candidate_label_created": False,
            }
        )
    return records


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs() -> None:
    verify_raw()
    sequences = build_sequences()
    exclusions = exclusion_records()
    SOURCE.mkdir(parents=True, exist_ok=True)

    write_csv(SOURCE / "candidate_records.csv", [], FIELDS)
    write_csv(SOURCE / "exclusions.csv", exclusions)

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, (variant, sequence) in sequences.items():
            handle.write(
                f">{sequence_id} Q7TNG8 residues 22-484 | N-terminal WSHPQFEK Strep tag | {variant}\n"
            )
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    pdb_rows = [
        {
            "pdb_id": pdb_id,
            "sequence_id": f"mldhd-{variant.lower()}",
            "variant": variant,
            "bound_state": state,
            "exact_polymer_length": len(sequences[f"mldhd-{variant.lower()}"][1]),
            "polymer_sha256": hashlib.sha256(
                sequences[f"mldhd-{variant.lower()}"][1].encode("ascii")
            ).hexdigest(),
            "raw_cif_sha256": sha256(RAW / f"{pdb_id}.cif"),
        }
        for pdb_id, (variant, state) in PDB_FORMS.items()
    ]
    write_csv(SOURCE / "pdb_construct_mapping.csv", pdb_rows)

    raw_files = sorted(path for path in RAW.iterdir() if path.is_file())
    provenance = {
        "source_id": SOURCE_ID,
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC10589216",
        "article_doi": DOI,
        "pmc_id": "PMC10589216",
        "article_published": "2023-10-20",
        "article_and_supplement_license": "CC-BY-4.0",
        "uniprot_license": "CC-BY-4.0",
        "pdb_coordinate_license": "CC0-1.0",
        "pubchem_notice": "US government work; PubChem download records retained verbatim",
        "reported_finite_steady_state_kcat_rows": len(exclusions),
        "eligible_direct_saturated_kcat_rows": 0,
        "accepted_records": 0,
        "excluded_coupled_unsaturated_acceptor_rows": len(exclusions),
        "defined_D_substrates_with_pubchem_structures": len(SUBSTRATES),
        "exact_assayed_construct_sequences": len(sequences),
        "paper_deposited_mldhd_pdb_entries": len(PDB_FORMS),
        "kinetics_sources": [
            "raw/PMC10589216-fullText.xml, Table 1 and Methods",
            "raw/41467_2023_42456_MOESM6_ESM.xlsx, Table 1 source data",
        ],
        "construct_sources": [
            "raw/Q7TNG8.fasta",
            "raw/PMC10589216-fullText.xml, Q7TNG8 residues 22-484 and N-terminal Strep-tag methods",
            "raw/8JD*.cif, exact deposited WT and H405A tagged polymer sequences",
        ],
        "construct_audit": "All proteins used Q7TNG8 residues 22-484 with the N-terminal WSHPQFEK Strep tag recovered from every deposited polymer. FASTA includes WT, active-site R347A/H398A/H405A/E442A/H443A, and disease-associated T228M/R347W/W351C/T440M assayed variants.",
        "pdb_series_audit": "The paper deposits 18 mLDHD entries: 8JDC-8JDG, 8JDN-8JDV, and 8JDX-8JDZ, plus out-of-order 8JDB. IDs 8JDH-8JDM and 8JDW are not claimed by the paper and are not downloaded.",
        "prosthetic_FAD_audit": "WT FAD:mLDHD occupancy was 0.73 +/- 0.04 mol/mol (SEM, n=3); no exogenous FAD was added. Authors state activities and kinetic parameters were corrected to active-enzyme concentration based on FAD occupancy.",
        "electron_acceptor_audit": "All nine substrate curves used the PMS-DCIP coupled assay with fixed 200 uM PMS and 100 uM DCIP. No PMS or DCIP concentration series, Km, plateau, or independent saturation evidence is reported. The physiological primary oxidant is explicitly unknown.",
        "selection_policy": "Retain only direct finite steady-state kcat for a chemically defined D-2-hydroxyacid with a PubChem isomeric structure when every coupled electron acceptor is demonstrated saturating. No row passes the electron-acceptor criterion.",
        "racemate_audit": "For five DL reagents, authors doubled reagent concentration and fitted against one-half concentration as D-isomer; structural deposits show only D-isomer bound. 8JDR identifies the bound branched substrate as (2R,3S), mapped to PubChem CID 10820562.",
        "assay": {
            "method": "PMS-DCIP coupled initial-rate readout at 600 nm",
            "buffer": "50 mM Tris-HCl pH 7.4",
            "temperature_C": 37,
            "volume_uL": 200,
            "enzyme_amount_ug": 2,
            "MnCl2_uM": 50,
            "PMS_uM": 200,
            "DCIP_uM": 100,
            "substrate_range_mM": "0-2",
            "replicates": 3,
            "error_type": "SEM",
            "fit": "Michaelis-Menten, GraphPad Prism 7.0",
        },
        "normalization": "Excluded source kcat and SEM values are divided by 60 from min-1 to s-1 without refitting; no excluded value is a candidate label.",
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "artifact_sha256": {
            "candidate_records.csv": sha256(SOURCE / "candidate_records.csv"),
            "exclusions.csv": sha256(SOURCE / "exclusions.csv"),
            "construct_sequences.fasta": sha256(fasta_path),
            "pdb_construct_mapping.csv": sha256(SOURCE / "pdb_construct_mapping.csv"),
        },
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="ascii"
    )

    write_homology_audit(sequences, exclusions, fasta_path)
    verify_outputs(sequences, exclusions)


def write_homology_audit(
    sequences: dict[str, tuple[str, str]], exclusions: list[dict[str, object]], fasta_path: Path
) -> None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if not hit_path.is_file() or not cluster_path.is_file():
        return
    hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
    families = {
        line.split("\t", 1)[0]
        for line in cluster_path.read_text(encoding="utf-8").splitlines()
        if line
    }
    audit = {
        "audited_on": "2026-07-22",
        "source_id": SOURCE_ID,
        "article_doi": DOI,
        "method": "MMseqs2 easy-search and linclust",
        "mmseqs_version": MMSEQS_VERSION,
        "min_identity": 0.3,
        "coverage": 0.8,
        "coverage_mode": 0,
        "development_target_sha256": REFERENCE_SHA256,
        "reported_finite_kcat_rows": len(exclusions),
        "candidate_records": 0,
        "excluded_acceptor_saturation_rows": len(exclusions),
        "unique_assayed_construct_sequences": len(sequences),
        "construct_sequences_sha256": sha256(fasta_path),
        "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
        "exact_sequence_overlap": sum(
            len(parts) > 2 and float(parts[2]) == 1.0
            for parts in (line.split("\t") for line in hit_lines)
        ),
        "homology_hits_sha256": sha256(hit_path),
        "candidate_mmseqs_families": len(families),
        "family_cluster_sha256": sha256(cluster_path),
        "accepted_records": 0,
        "accepted_unique_sequences": 0,
        "accepted_unique_substrates": 0,
        "readiness_gate_passes": False,
        "claim_boundary": "Exclusion-only source audit; no candidate labels or model predictions were generated.",
    }
    (SOURCE / "homology-audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii"
    )


def verify_outputs(
    sequences: dict[str, tuple[str, str]], exclusions: list[dict[str, object]]
) -> None:
    if len(sequences) != 10 or len(set(sequence for _, sequence in sequences.values())) != 10:
        raise ValueError("Expected ten distinct exact assayed constructs")
    if len(exclusions) != 9 or any(row["candidate_label_created"] for row in exclusions):
        raise ValueError("Expected nine excluded and zero accepted kinetic rows")
    with (SOURCE / "candidate_records.csv").open(newline="", encoding="utf-8") as handle:
        if list(csv.DictReader(handle)):
            raise ValueError("Excluded apparent values leaked into candidate records")
    with (SOURCE / "exclusions.csv").open(newline="", encoding="utf-8") as handle:
        written = list(csv.DictReader(handle))
    if len(written) != 9 or {int(row["substrate_pubchem_cid"]) for row in written} != {
        cid for cid, _ in SUBSTRATES.values()
    }:
        raise ValueError("Exclusion artifact failed cardinality or PubChem identity check")
    audit_path = SOURCE / "homology-audit.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text(encoding="ascii"))
        if audit["accepted_records"] != 0 or audit["coverage_mode"] != 0:
            raise ValueError("Homology audit violates frozen selection boundary")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        download_raw()
    write_outputs()


if __name__ == "__main__":
    main()
