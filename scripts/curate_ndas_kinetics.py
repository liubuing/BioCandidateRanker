from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "biostudies-S-EPMC11760959"
TARGET = ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723" / "homology" / "unikp_reference.fasta"

# Supplementary Methods pp. 4-5: Ndas1146 carries this vector-derived N terminus,
# native residues 2-161, and a C-terminal His tag. Ndas1147 is untagged.
NDAS1146_PREFIX = "MDTGSSEPDANRCPSQRSSHALQTLTTRRAVRAFADRPV"
NDAS1146_SUFFIX = "GSSHHHHHH"

SUBSTRATES = {
    "cWF": ("cyclo(L-Trp-L-Phe)", 7408486),
    "cWY": ("cyclo(L-Trp-L-Tyr)", 7408259),
    "cWS": ("cyclo(L-Trp-L-Ser)", 70694585),
    "cFG": ("cyclo(L-Phe-Gly)", 7076549),
    "cLP": ("cyclo(L-Leu-L-Pro)", 7074739),
    "cFP": ("cyclo(L-Phe-L-Pro)", 443440),
    "cHF": ("cyclo(L-His-L-Phe)", 7408279),
}

# Supplementary Table 3: complex, substrate, kcat, kcat SE, Km, Km SE.
KINETICS = [
    ("ndascdo-wt", "cWF", 0.05, 0.01, 0.08, 0.01),
    ("ndascdo-wt", "cWY", 0.11, 0.01, 0.17, 0.02),
    ("ndascdo-wt", "cWS", 0.076, 0.002, 0.57, 0.05),
    ("ndascdo-wt", "cFG", 0.87, 0.08, 4.2, 0.8),
    ("ndascdo-s58a", "cFG", 0.12, 0.02, 4.9, 1.4),
    ("ndascdo-wt", "cLP", 2.07, 0.06, 2.5, 0.2),
    ("ndascdo-wt", "cFP", 4.6, 0.2, 0.9, 0.1),
    ("ndascdo-wt", "cHF", 17.3, 0.7, 0.57, 0.08),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> str:
    sequence = "".join(
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith(">")
    )
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(f"Invalid protein FASTA: {path}")
    return sequence


def mutate(sequence: str, position: int, expected: str, replacement: str) -> str:
    if sequence[position - 1] != expected:
        raise ValueError(
            f"Expected {expected}{position}, found {sequence[position - 1]}{position}"
        )
    return sequence[: position - 1] + replacement + sequence[position:]


def build_sequences() -> tuple[dict[str, str], dict[str, str]]:
    native_a = read_fasta(SOURCE / "raw" / "D7B1W6.fasta")
    native_b = read_fasta(SOURCE / "raw" / "D7B1W7.fasta")
    if len(native_a) != 157 or len(native_b) != 105:
        raise ValueError(f"Unexpected UniProt lengths: D7B1W6={len(native_a)}, D7B1W7={len(native_b)}")

    construct_a = NDAS1146_PREFIX + native_a[1:] + NDAS1146_SUFFIX
    construct_a_s58a = mutate(construct_a, 58, "S", "A")
    native_a_s20a = mutate(native_a, 20, "S", "A")
    if construct_a_s58a != NDAS1146_PREFIX + native_a_s20a[1:] + NDAS1146_SUFFIX:
        raise ValueError("NdasCDO-A construct/native mutation mapping is inconsistent")

    constructs = {
        "ndascdo-wt-a": construct_a,
        "ndascdo-s58a-a": construct_a_s58a,
        "ndascdo-b": native_b,
    }
    homology_queries = {
        "ndascdo-wt-a": native_a,
        "ndascdo-s58a-a": native_a_s20a,
        "ndascdo-b": native_b,
    }
    return constructs, homology_queries


def pubchem_structure(code: str) -> tuple[int, str, str]:
    payload = json.loads((SOURCE / "raw" / f"pubchem-{code}.json").read_text(encoding="utf-8"))
    properties = payload["PropertyTable"]["Properties"]
    if len(properties) != 1:
        raise ValueError(f"Expected one PubChem structure for {code}")
    item = properties[0]
    expected_cid = SUBSTRATES[code][1]
    if item["CID"] != expected_cid:
        raise ValueError(f"Expected PubChem CID {expected_cid} for {code}, found {item['CID']}")
    return item["CID"], item["SMILES"], item["Title"]


def component_ids(complex_id: str) -> tuple[str, str]:
    chain_a = "ndascdo-s58a-a" if complex_id == "ndascdo-s58a" else "ndascdo-wt-a"
    return chain_a, "ndascdo-b"


def build_records() -> list[dict[str, object]]:
    records = []
    for index, (complex_id, code, kcat, kcat_error, km, km_error) in enumerate(KINETICS, 1):
        substrate_name, _ = SUBSTRATES[code]
        cid, smiles, pubchem_title = pubchem_structure(code)
        chain_a, chain_b = component_ids(complex_id)
        variant = "NdasCDO-A S58A (D7B1W6 native-coordinate S20A)" if complex_id.endswith("s58a") else "WT"
        records.append({
            "candidate_id": f"ndascdo-{index:03d}",
            "article_doi": "10.1038/s41467-025-56127-y",
            "source_table": "Supplementary Table 3",
            "source_row": f"{code}{'-S58A' if complex_id.endswith('s58a') else ''}",
            "organism": "Nocardiopsis dassonvillei",
            "sequence_id": complex_id,
            "chain_a_sequence_id": chain_a,
            "chain_a_accession": "D7B1W6",
            "chain_b_sequence_id": chain_b,
            "chain_b_accession": "D7B1W7",
            "complex_composition": "filament repeat A2B2; both chain types required",
            "construct": f"pRSFDuet-1 coexpression; Ndas1146 vector-derived N terminus/native residues 2-157/C-terminal GSSHHHHHH; untagged native Ndas1147; {variant}",
            "variable_substrate": substrate_name,
            "substrate_code": code,
            "substrate_pubchem_cid": cid,
            "substrate_pubchem_title": pubchem_title,
            "substrate_isomeric_smiles": smiles,
            "endpoint": "direct_initial_rate_steady_state_kcat_s-1",
            "kcat_s-1": kcat,
            "kcat_error_s-1": kcat_error,
            "km_mM": km,
            "km_error_mM": km_error,
            "error_type": "standard error of fit across three independent experiments",
            "assay_method": "direct absorbance initial-rate Michaelis-Menten assay; readings below 10% substrate conversion",
            "assay_buffer": "50 mM Tris, 200 mM NaCl",
            "assay_pH": 9.0,
            "assay_temperature_C": 30,
            "enzyme_concentration_uM": 0.1,
            "replicates": 3,
            "status_at_normalization": "pending_homology",
        })
    return records


def read_hit_ids() -> set[str] | None:
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    if not hit_path.exists():
        return None
    return {
        line.split("\t", 1)[0]
        for line in hit_path.read_text(encoding="utf-8").splitlines()
        if line
    }


def apply_homology_status(records: list[dict[str, object]]) -> None:
    hit_ids = read_hit_ids()
    if hit_ids is None:
        return
    for record in records:
        components = {record["chain_a_sequence_id"], record["chain_b_sequence_id"]}
        record["status_at_normalization"] = (
            "excluded_homology"
            if components & hit_ids
            else "accepted_homology_cold_pool"
        )


def write_fasta(path: Path, records: dict[str, str], descriptions: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id, sequence in records.items():
            handle.write(f">{sequence_id} {descriptions[sequence_id]}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_outputs() -> None:
    constructs, homology_queries = build_sequences()
    records = build_records()
    apply_homology_status(records)

    fasta_path = SOURCE / "construct_sequences.fasta"
    write_fasta(fasta_path, constructs, {
        "ndascdo-wt-a": "D7B1W6-derived expressed NdasCDO-A chain | WT A2 component",
        "ndascdo-s58a-a": "D7B1W6-derived expressed NdasCDO-A chain | construct S58A/native S20A A2 component",
        "ndascdo-b": "D7B1W7 exact native NdasCDO-B chain | shared B2 component",
    })
    query_path = SOURCE / "homology" / "homology_queries.fasta"
    write_fasta(query_path, homology_queries, {
        "ndascdo-wt-a": "D7B1W6 exact native sequence used for coverage-aware homology",
        "ndascdo-s58a-a": "D7B1W6 S20A native-coordinate equivalent of construct S58A",
        "ndascdo-b": "D7B1W7 exact native sequence used for coverage-aware homology",
    })

    with (SOURCE / "candidate_records.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    accepted = [row for row in records if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    raw_hashes = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted((SOURCE / "raw").iterdir())
        if path.is_file()
    }
    provenance = {
        "source_id": "biostudies-S-EPMC11760959",
        "stable_record_url": "https://www.ebi.ac.uk/biostudies/studies/S-EPMC11760959",
        "article_doi": "10.1038/s41467-025-56127-y",
        "article_published": "2025-01-24",
        "article_license": "CC-BY-4.0",
        "reported_steady_state_kcat_records": len(records),
        "accepted_records": len(accepted),
        "unique_complex_variants": 2,
        "unique_substrates": len(SUBSTRATES),
        "kinetics_source": "raw/41467_2025_56127_MOESM1_ESM.pdf, Supplementary Table 3; assay provenance in main Methods and Fig. 5",
        "sequence_sources": ["raw/D7B1W6.fasta", "raw/D7B1W7.fasta", "Supplementary Methods pp. 4-5"],
        "complex_representation": "Each candidate references separate A and B chain sequence IDs. The biological repeat is A2B2 and propagates as a filament; no artificial chain concatenation is used.",
        "assay": {
            "method": "Direct absorbance initial-rate Michaelis-Menten assay",
            "conversion_gate": "Only readings below 10% substrate conversion",
            "buffer": "50 mM Tris, 200 mM NaCl, pH 9.0",
            "temperature_C": 30,
            "enzyme_concentration_uM": 0.1,
            "replicates": 3,
            "error_type": "standard error of fit across three independent experiments",
        },
        "selection_policy": "All direct finite steady-state absolute kcat rows in Supplementary Table 3; no kinetic-value selection.",
        "excluded_data": [
            "cWG, cWW, cFF, and cHP endpoint/LC-MS activity: no saturation-derived kinetic parameters",
            "pH-rate and temperature series: condition/model analyses rather than the standard-condition substrate table",
            "solvent kinetic isotope effect and viscosity values: perturbation/ratio analyses",
            "stopped-flow multiple-turnover bursts and Kintek global fits: pre-steady-state/model-derived values",
            "LC-MS completion and time-course values: endpoint/progress-curve data, not steady-state absolute kcat",
            "docking and Eyring parameters: computational or derived model values",
        ],
        "raw_files": raw_hashes,
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.exists() and cluster_path.exists():
        hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        family_lines = [line for line in cluster_path.read_text(encoding="utf-8").splitlines() if line]
        families = {line.split("\t", 1)[0] for line in family_lines}
        excluded_complexes = {
            row["sequence_id"] for row in records if row["status_at_normalization"] == "excluded_homology"
        }
        audit = {
            "audited_on": "2026-07-22",
            "source_id": "biostudies-S-EPMC11760959",
            "article_doi": "10.1038/s41467-025-56127-y",
            "method": "MMseqs2 easy-search of native-coordinate component chains and linclust of candidate components",
            "complex_homology_rule": "Exclude a heteromeric complex when either required component chain has a qualifying hit.",
            "mmseqs_version": "5d152c612b6ad2a56f657b7a02c127eceaea2a75",
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target": str(TARGET.relative_to(ROOT)).replace("\\", "/"),
            "development_target_sha256": sha256(TARGET),
            "candidate_records": len(records),
            "candidate_complex_variants": 2,
            "component_query_sequences": len(homology_queries),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_queries_sha256": sha256(query_path),
            "homology_hit_chains": len({line.split("\t", 1)[0] for line in hit_lines}),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_component_families": len(families),
            "family_representatives": sorted(families),
            "family_cluster_sha256": sha256(cluster_path),
            "excluded_complex_variants": sorted(excluded_complexes),
            "accepted_records": len(accepted),
            "accepted_complex_variants": len({row["sequence_id"] for row in accepted}),
            "accepted_unique_substrates": len({row["variable_substrate"] for row in accepted}),
            "readiness_gate_passes": False,
            "claim_boundary": "Curation pool only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    write_outputs()
