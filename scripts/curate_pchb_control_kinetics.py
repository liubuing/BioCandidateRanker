from __future__ import annotations

import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "external" / "temporal-absolute-kinetics" / "europepmc-PMC12329711"
RAW = SOURCE / "raw"
REFERENCE = ROOT / "artifacts" / "external" / "absolute-kinetics-screen" / "dryad-4964723" / "homology" / "unikp_reference.fasta"
DOI = "10.1021/acs.biochem.5c00157"
MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
REFERENCE_SHA256 = "69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9"
CHORISMATE_SMILES = "C=C(C(=O)O)O[C@@H]1C=C(C=C[C@H]1O)C(=O)O"
UNRESOLVED_STATUS = "excluded_construct_unresolved"
UNRESOLVED_IDS = {
    "control-mtcm",
    "control-mtcm-v73e",
    "control-mtcm-r49l",
    "control-ecm",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> str:
    sequence = "".join(
        line.strip() for line in path.read_text(encoding="ascii").splitlines()
        if line and not line.startswith(">")
    )
    if not sequence or not sequence.isalpha() or not sequence.isupper():
        raise ValueError(f"Invalid FASTA: {path}")
    return sequence


def mutate(sequence: str, position: int, expected: str, replacement: str) -> str:
    if sequence[position - 1] != expected:
        raise ValueError(f"Expected {expected}{position}, found {sequence[position - 1]}{position}")
    return sequence[: position - 1] + replacement + sequence[position:]


def sequences() -> dict[str, dict[str, str]]:
    ec = read_fasta(RAW / "P0A9J8.fasta")
    mj = read_fasta(RAW / "Q57696.fasta")
    mt = read_fasta(RAW / "P9WIB9.fasta")
    sc = read_fasta(RAW / "P32178.fasta")
    mmj = read_fasta(RAW / "2GTV.fasta")

    # The current article's MtCM forward primer encodes M followed by precursor
    # residues 34-end; substitutions retain the published precursor numbering.
    def mt_core(parent: str) -> str:
        return "M" + parent[33:]
    return {
        "control-ecm": {
            "sequence": ec[:109],
            "source": "P0A9J8 residues 1-109; 1ECM exact domain sequence",
        },
        "control-mjcm-prime": {
            "sequence": mj[:93],
            "source": "Q57696 residues 1-93; cited six-residue C-terminal truncation",
        },
        "control-mmjcm": {
            "sequence": mmj,
            "source": "RCSB 2GTV polymer entity 1, engineered monomeric construct",
        },
        "control-mtcm": {
            "sequence": mt_core(mt),
            "source": "P9WIB9 precursor residues 34-199 with primer-encoded initiator Met",
        },
        "control-mtcm-v73e": {
            "sequence": mt_core(mutate(mt, 73, "V", "E")),
            "source": "P9WIB9 V73E; current-article primer-defined leaderless construct",
        },
        "control-mtcm-r49l": {
            "sequence": mt_core(mutate(mt, 49, "R", "L")),
            "source": "P9WIB9 R49L; current-article primer-defined leaderless construct",
        },
        "control-sccm": {
            "sequence": sc,
            "source": "P32178 full coding sequence; current-article ARO7 primers",
        },
    }


ROWS = (
    ("*MtCM", "Mycobacterium tuberculosis H37Rv", "control-mtcm", "P9WIB9", 49.0, 11.0, 170.0, 16.0, 30, "50 mM potassium phosphate"),
    ("*MtCM V73E", "Mycobacterium tuberculosis H37Rv", "control-mtcm-v73e", "P9WIB9 V73E", 13.6, 2.8, 560.0, 60.0, 30, "50 mM potassium phosphate"),
    ("*MtCM R49L", "Mycobacterium tuberculosis H37Rv", "control-mtcm-r49l", "P9WIB9 R49L", 0.0070, 0.0019, 220.0, 80.0, 30, "50 mM potassium phosphate"),
    ("EcCM", "Escherichia coli K-12", "control-ecm", "P0A9J8 residues 1-109", 55.0, "", 240.0, "", 30, "50 mM potassium phosphate"),
    ("MjCM'", "Methanocaldococcus jannaschii", "control-mjcm-prime", "Q57696 residues 1-93", 5.5, "", 16.0, "", 30, "50 mM potassium phosphate"),
    ("mMjCM", "Methanocaldococcus jannaschii", "control-mmjcm", "Q57696 engineered monomer; 2GTV", 3.5, "", 270.0, "", 20, "10 mM potassium phosphate, 160 mM NaCl, 0.1 mg/mL BSA"),
    ("ScCM", "Saccharomyces cerevisiae", "control-sccm", "P32178", 37.0, "", 590.0, "", 30, "50 mM potassium phosphate, 10 uM L-tryptophan"),
)


def records() -> list[dict[str, object]]:
    output = []
    for index, row in enumerate(ROWS, 1):
        enzyme, organism, sequence_id, accession, kcat, kcat_sd, km, km_sd, temperature, buffer = row
        output.append({
            "candidate_id": f"pchb-control-{index:03d}",
            "article_doi": DOI,
            "source_table": "Table 4",
            "source_cell": f"{enzyme}; CM activity kcat",
            "organism": organism,
            "enzyme_identity": enzyme,
            "sequence_id": sequence_id,
            "accession_and_variant": accession,
            "construct_sequence_scope": "authoritative catalytic polypeptide only; not the exact assayed His6-tagged chain and never admissible as a substitute",
            "reaction": "chorismate -> prephenate",
            "substrate_name": "chorismate",
            "substrate_pubchem_cid": 12039,
            "substrate_isomeric_smiles": CHORISMATE_SMILES,
            "endpoint": "kcat_s-1",
            "kcat_s-1": kcat,
            "kcat_sd_s-1": kcat_sd,
            "km_uM": km,
            "km_sd_uM": km_sd,
            "error_type": "SD of independently fitted biological-replicate parameters" if kcat_sd != "" else "not reported for current measurement",
            "biological_replicates": "at least 2",
            "fit": "direct initial rates fitted to the Michaelis-Menten equation",
            "assay_buffer": buffer,
            "assay_pH": 7.5,
            "assay_temperature_C": temperature,
            "status_at_normalization": (
                UNRESOLVED_STATUS if sequence_id in UNRESOLVED_IDS else "pending_homology"
            ),
        })
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify_sources() -> None:
    root = ET.parse(RAW / "PMC12329711-fullText.xml").getroot()
    if root.findtext(".//article-id[@pub-id-type='doi']") != DOI:
        raise ValueError("Article DOI mismatch")
    table = root.find(".//table-wrap[@id='tbl4']")
    evidence = " ".join(table.itertext()) if table is not None else ""
    for token in ("49 ± 11", "13.6 ± 2.8", "0.0070", "55/", "5.5/", "3.5/", "37/"):
        if token not in evidence:
            raise ValueError(f"Missing Table 4 evidence: {token}")
    compound = json.loads((RAW / "pubchem-chorismate.json").read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
    if compound["CID"] != 12039 or compound.get("SMILES", compound.get("IsomericSMILES")) != CHORISMATE_SMILES:
        raise ValueError("Unexpected PubChem chorismate structure")
    if sha256(REFERENCE) != REFERENCE_SHA256:
        raise ValueError("Development reference hash changed")


def write_outputs() -> None:
    verify_sources()
    sequence_map = sequences()
    candidate_rows = records()
    hit_path = SOURCE / "homology" / "homology_hits.tsv"
    hits = {
        line.split("\t", 1)[0]
        for line in hit_path.read_text(encoding="utf-8").splitlines()
        if line
    } if hit_path.is_file() else set()
    for row in candidate_rows:
        if row["sequence_id"] in hits:
            row["status_at_normalization"] = "excluded_homology"

    construct_resolution = {
        "schema_version": 1,
        "audited_on": "2026-07-27",
        "source_id": "europepmc-PMC12329711",
        "decision_rule": (
            "An accession-derived catalytic sequence cannot substitute for an unresolved "
            "full assayed construct. Only an exact, source-supported translated chain may "
            "receive a no-hit admission decision."
        ),
        "searches": [
            {
                "source": "current article and supporting information",
                "locator": "raw/PMC12329711-fullText.xml; raw/bi5c00157_si_001.pdf",
                "result": "His6 purification is stated, but complete translated tag junctions are absent",
            },
            {
                "source": "cited MtCM construct paper",
                "locator": "10.1111/j.1742-4658.2004.04478.x",
                "result": "pKTU3-HCT lineage identified; complete plasmid/translated junction unavailable",
            },
            {
                "source": "cited EcCM pATCH construct paper",
                "locator": "10.1021/bi980449t",
                "result": "pET-EcCM-pATCH identified; complete assayed fusion sequence unavailable",
            },
            {
                "source": "RCSB PDB/mmCIF polymer records",
                "locator": "1ECM entity 1; 2FP2 entity 1",
                "result": "catalytic polymers resolve domain boundaries but contain no purification tag junction",
            },
            {
                "source": "current-article cloning primers and expression vectors",
                "locator": "Experimental Procedures sections 2.2 and 2.6",
                "result": (
                    "MtCM termini are partly constrained; R49L text conflicts by naming an "
                    "XhoI-bearing reverse primer but an NdeI-SpeI ligation, so translation "
                    "through the pET21a junction is not uniquely established"
                ),
            },
            {
                "source": "public plasmid archive search",
                "locator": "Addgene catalog searches for pKTU3-HCT and pET-EcCM-pATCH",
                "result": "no archived plasmid record or sequence located",
            },
        ],
        "rows": [
            {
                "candidate_id": "pchb-control-001",
                "enzyme": "*MtCM",
                "sequence_id": "control-mtcm",
                "status": UNRESOLVED_STATUS,
                "evidence_codes": [
                    "FULL_HIS6_JUNCTION_UNDISCLOSED",
                    "UNPUBLISHED_PLASMID_SEQUENCE_UNAVAILABLE",
                    "CATALYTIC_ACCESSION_NOT_SUBSTITUTED",
                ],
            },
            {
                "candidate_id": "pchb-control-002",
                "enzyme": "*MtCM V73E",
                "sequence_id": "control-mtcm-v73e",
                "status": UNRESOLVED_STATUS,
                "evidence_codes": [
                    "FULL_HIS6_JUNCTION_UNDISCLOSED",
                    "ARCHIVED_PLASMID_SEQUENCE_NOT_FOUND",
                    "CATALYTIC_ACCESSION_NOT_SUBSTITUTED",
                ],
            },
            {
                "candidate_id": "pchb-control-003",
                "enzyme": "*MtCM R49L",
                "sequence_id": "control-mtcm-r49l",
                "status": UNRESOLVED_STATUS,
                "evidence_codes": [
                    "PRIMER_VECTOR_JUNCTION_CONFLICT",
                    "ARCHIVED_PLASMID_SEQUENCE_NOT_FOUND",
                    "CATALYTIC_ACCESSION_NOT_SUBSTITUTED",
                ],
            },
            {
                "candidate_id": "pchb-control-004",
                "enzyme": "EcCM",
                "sequence_id": "control-ecm",
                "status": UNRESOLVED_STATUS,
                "evidence_codes": [
                    "FULL_HIS6_JUNCTION_UNDISCLOSED",
                    "ARCHIVED_PLASMID_SEQUENCE_NOT_FOUND",
                    "PDB_POLYMER_LACKS_ASSAY_TAG",
                    "CATALYTIC_ACCESSION_NOT_SUBSTITUTED",
                ],
            },
        ],
        "exact_assayed_chains_reconstructed": 0,
        "exact_chain_mmseqs_runs_required": 0,
        "terminal": True,
        "model_predictions_run": False,
    }
    resolution_path = SOURCE / "construct-resolution.json"
    resolution_path.write_text(
        json.dumps(construct_resolution, indent=2) + "\n", encoding="ascii"
    )

    fasta_path = SOURCE / "construct_sequences.fasta"
    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_id in sorted(sequence_map):
            item = sequence_map[sequence_id]
            handle.write(f">{sequence_id} | {item['source']}\n")
            sequence = item["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    write_csv(SOURCE / "candidate_records.csv", candidate_rows)

    status_by_sequence = {row["sequence_id"]: row["status_at_normalization"] for row in candidate_rows}
    resolution = []
    for row in ROWS:
        status = status_by_sequence[row[2]]
        resolution.append({
            "held_row": row[0],
            "table_4_cm_kcat": row[4],
            "resolution": status,
            "sequence_id": row[2],
            "reason": (
                "frozen development-corpus MMseqs hit"
                if status == "excluded_homology"
                else "terminal exclusion: exact assayed His6-tagged chain could not be reconstructed; catalytic sequence was not substituted"
            ),
        })
    resolution.append({"held_row": "*MtCM R134L", "table_4_cm_kcat": "ND", "resolution": "excluded_nonfinite", "sequence_id": "control-mtcm-r134l", "reason": "individual kcat not determinable"})
    write_csv(SOURCE / "control-row-resolution.csv", resolution)

    exclusions = [
        {"scope": "all PchB-derived Table 2/Table 4 rows", "source_value": "19 rows", "reason": "excluded_homology: Q51507 exact development overlap; applies to WT and every derived variant", "candidate_label_created": False},
        {"scope": "*MtCM R134L CM kcat", "source_value": "ND", "reason": "nonfinite endpoint", "candidate_label_created": False},
        {"scope": "*MtCM IPL kcat", "source_value": "(0.0009 +/- 0.0002)", "reason": "parenthesized rough estimate, not eligible direct kcat", "candidate_label_created": False},
        {"scope": "*MtCM V73E/R49L/R134L and EcCM/MjCM'/mMjCM/ScCM IPL kcat", "source_value": "NA", "reason": "no finite kcat; efficiency-only value or detection limit", "candidate_label_created": False},
        {"scope": "italicized EcCM/MjCM'/mMjCM/ScCM comparison values", "source_value": "69; 5.7; 3.2; 13", "reason": "independent pre-2023 literature values, not current measurements", "candidate_label_created": False},
    ]
    write_csv(SOURCE / "exclusions.csv", exclusions)

    raw_files = sorted(RAW.glob("*.fasta")) + [
        RAW / "PMC12329711-fullText.xml",
        RAW / "PMC12329711-supplementaryFiles.zip",
        RAW / "bi5c00157_si_001.pdf",
        RAW / "pubchem-chorismate.json",
    ]
    accepted = [row for row in candidate_rows if row["status_at_normalization"] == "accepted_homology_cold_pool"]
    pending = [row for row in candidate_rows if str(row["status_at_normalization"]).startswith("pending")]
    provenance = {
        "source_id": "europepmc-PMC12329711",
        "stable_record_url": "https://europepmc.org/articles/PMC12329711",
        "article_doi": DOI,
        "article_published": "2025-07-22",
        "article_and_supporting_information_license": "CC-BY-4.0",
        "held_independent_control_rows": 8,
        "finite_direct_control_kcat_candidates": 7,
        "nonfinite_control_rows": 1,
        "pchb_derived_rows_excluded": 19,
        "selection_policy": "Only current, finite, directly fitted Table 4 CM kcat values are candidates. ND, NA, detection limits, rough parenthesized estimates, italicized historical comparison values, and every PchB-derived row are excluded.",
        "construct_boundary": "Authoritative catalytic polypeptides are retained only as conservative exclusion queries. Four no-hit rows are terminally excluded because their exact assayed His6-tagged chains could not be reconstructed; no catalytic accession was substituted.",
        "substrate_structure_source": "PubChem CID 12039 PUG REST isomeric SMILES, retrieved 2026-07-27",
        "kinetics_source": "raw/PMC12329711-fullText.xml, Table 4 and Experimental Procedures",
        "sequence_sources": [item["source"] for item in sequence_map.values()],
        "raw_file_sha256": {path.name: sha256(path) for path in raw_files},
        "artifact_sha256": {
            "candidate_records.csv": sha256(SOURCE / "candidate_records.csv"),
            "construct-resolution.json": sha256(resolution_path),
            "control-row-resolution.csv": sha256(SOURCE / "control-row-resolution.csv"),
            "construct_sequences.fasta": sha256(fasta_path),
            "exclusions.csv": sha256(SOURCE / "exclusions.csv"),
        },
        "accepted_after_frozen_homology": len(accepted),
        "model_predictions_run": False,
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="ascii")

    cluster_path = SOURCE / "homology" / "family-cluster" / "proteins_cluster.tsv"
    if hit_path.is_file() and cluster_path.is_file():
        hit_lines = [line for line in hit_path.read_text(encoding="utf-8").splitlines() if line]
        cluster_lines = [line for line in cluster_path.read_text(encoding="utf-8").splitlines() if line]
        audit = {
            "audited_on": "2026-07-27",
            "source_id": "europepmc-PMC12329711",
            "article_doi": DOI,
            "method": "MMseqs2 easy-search and linclust",
            "mmseqs_version": MMSEQS_VERSION,
            "min_identity": 0.3,
            "coverage": 0.8,
            "coverage_mode": 0,
            "development_target_sha256": REFERENCE_SHA256,
            "candidate_records": len(candidate_rows),
            "held_control_rows_resolved": 8,
            "unique_homology_query_sequences": len(sequence_map),
            "construct_sequences_sha256": sha256(fasta_path),
            "homology_hit_sequences": len({line.split("\t", 1)[0] for line in hit_lines}),
            "exact_sequence_overlap": sum(float(line.split("\t")[2]) == 1.0 for line in hit_lines),
            "homology_hits_sha256": sha256(hit_path),
            "candidate_mmseqs_families": len({line.split("\t", 1)[0] for line in cluster_lines}),
            "family_cluster_sha256": sha256(cluster_path),
            "accepted_records": len(accepted),
            "accepted_unique_sequences": len({row["sequence_id"] for row in accepted}),
            "terminal_construct_exclusions": len(construct_resolution["rows"]),
            "readiness_gate_passes": not pending,
            "blockers": [],
            "claim_boundary": "Control-row curation and homology exclusion only; no model predictions were generated.",
        }
        (SOURCE / "homology-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")

    if (
        len(candidate_rows) != 7
        or len(resolution) != 8
        or pending
        or {row["sequence_id"] for row in candidate_rows if row["status_at_normalization"] == UNRESOLVED_STATUS} != UNRESOLVED_IDS
        or any(not math.isfinite(float(row["kcat_s-1"])) for row in candidate_rows)
    ):
        raise ValueError("Control-row cardinality or finite-kcat invariant failed")


if __name__ == "__main__":
    write_outputs()
