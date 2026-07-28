"""Build the label-independent global temporal family readiness audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any

from biocandidate.data.homology import read_mmseqs_clusters, run_mmseqs_easy_cluster
from biocandidate.temporal_registry import finalized_accepted_sources


MMSEQS_VERSION = "5d152c612b6ad2a56f657b7a02c127eceaea2a75"
DEFAULT_MMSEQS = "mmseqs"

# Explicit schema decisions. Adding a source requires declaring which exact FASTA IDs
# constitute its assayed construct; no fallback to names, accessions, or similarity exists.
SOURCE_COMPONENT_COLUMNS = {
    "biostudies-S-EPMC11955521": ("sequence_id",),
    "biostudies-S-EPMC10795190": ("component_sequence_ids",),
    "biostudies-S-EPMC11760959": ("chain_a_sequence_id", "chain_b_sequence_id"),
    "biostudies-S-EPMC12412097": ("sequence_id",),
    "biostudies-S-EPMC12645435": ("uniprot_accession", "variant"),
    "biostudies-S-EPMC12667214": ("sequence_id",),
    "biostudies-S-EPMC12781114": ("sequence_id",),
    "biostudies-S-EPMC12751057": ("sequence_id",),
    "biostudies-S-EPMC13033379": ("sequence_id",),
    "edmond-S0HJ48": ("sequence_id",),
    "europepmc-PMC10520331": ("sequence_id",),
    "europepmc-PMC11656708": ("sequence_id",),
    "europepmc-PMC11659886": ("sequence_id",),
    "europepmc-PMC12284513": ("sequence_id",),
    "europepmc-PMC12329711": ("sequence_id",),
    "europepmc-PMC12362431": ("sequence_id",),
    "europepmc-PMC12444516": ("sequence_id",),
    "europepmc-PMC12838360": ("sequence_id",),
    "zenodo-14055918": ("sequence_id",),
}

SOURCE_FASTA_FILES = {
    "biostudies-S-EPMC12645435": "variant_sequences.fasta",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                entries[identifier] = "".join(parts)
            identifier = line[1:].split()[0]
            if not identifier or identifier in entries:
                raise ValueError(f"{path}: duplicate or empty FASTA ID at line {line_number}")
            parts = []
        elif identifier is None:
            raise ValueError(f"{path}: sequence before first header at line {line_number}")
        else:
            sequence = line.upper()
            if not sequence.isalpha():
                raise ValueError(f"{path}: invalid sequence characters at line {line_number}")
            parts.append(sequence)
    if identifier is not None:
        entries[identifier] = "".join(parts)
    if not entries or any(not sequence for sequence in entries.values()):
        raise ValueError(f"{path}: empty FASTA or sequence")
    return entries


def _component_ids(
    row: dict[str, str], columns: tuple[str, ...], *, source_id: str
) -> list[str]:
    if source_id == "biostudies-S-EPMC12645435":
        accession = (row.get("uniprot_accession") or "").strip()
        variant = (row.get("variant") or "").strip()
        if not accession or not variant:
            raise ValueError("CtNDT mapping requires uniprot_accession and variant")
        return [f"{accession}|{variant}"]
    identifiers: list[str] = []
    for column in columns:
        value = (row.get(column) or "").strip()
        if not value:
            raise ValueError(f"missing required mapping column {column}")
        identifiers.extend(part.strip() for part in value.split(";") if part.strip())
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ValueError("component mapping is empty or repeats an ID")
    return identifiers


def _substrate_id(row: dict[str, str]) -> str:
    for column in (
        "substrate_pubchem_cid",
        "variable_substrate_pubchem_cid",
        "substrate_1_pubchem_cid",
    ):
        value = (row.get(column) or "").strip()
        if value:
            return f"pubchem:{value}"
    raise ValueError("missing supported substrate identifier")


def _mmseqs_command(executable: str, *arguments: str) -> list[str]:
    if not executable.startswith("wsl:"):
        return [executable, *arguments]
    parts = executable.split(":", 2)
    if len(parts) != 3:
        raise ValueError("WSL MMseqs2 must use wsl:<distribution>:<command>")
    command = ["wsl"]
    if parts[1]:
        command += ["-d", parts[1]]
    return [*command, parts[2], *arguments]


def _verify_mmseqs(executable: str) -> str:
    completed = subprocess.run(
        _mmseqs_command(executable, "version"), capture_output=True, text=True, check=False
    )
    if completed.returncode:
        raise RuntimeError(f"MMseqs2 version command failed: {completed.stderr.strip()}")
    version = (completed.stdout or completed.stderr).strip().splitlines()[-1]
    if version != MMSEQS_VERSION:
        raise RuntimeError(f"MMseqs2 version {version!r} != frozen version {MMSEQS_VERSION!r}")
    return version


def _reset_mmseqs_workspace(path: Path, executable: str) -> None:
    if not path.exists():
        return
    if not executable.startswith("wsl:"):
        shutil.rmtree(path)
        return
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        raise ValueError(f"cannot map MMseqs workspace to WSL: {resolved}")
    wsl_path = f"/mnt/{drive}/{resolved.as_posix()[3:]}"
    _, distribution, _ = executable.split(":", 2)
    command = ["wsl"]
    if distribution:
        command += ["-d", distribution]
    command += ["rm", "-rf", "--", wsl_path]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode or path.exists():
        raise RuntimeError(f"could not reset managed MMseqs workspace: {completed.stderr.strip()}")


def _record_families(records: list[dict[str, Any]], clusters: dict[str, str]) -> dict[int, str]:
    parents = list(range(len(records)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    seen: dict[str, int] = {}
    for index, record in enumerate(records):
        for global_id in record["global_sequence_ids"]:
            representative = clusters[global_id]
            if representative in seen:
                union(index, seen[representative])
            else:
                seen[representative] = index
    roots = sorted({find(index) for index in range(len(records))})
    names = {root: f"family-{position:04d}" for position, root in enumerate(roots, 1)}
    return {index: names[find(index)] for index in range(len(records))}


def build_global_audit(
    root: Path,
    protocol_path: Path,
    output_dir: Path,
    *,
    mmseqs: str | None = DEFAULT_MMSEQS,
    cluster_tsv: Path | None = None,
    threads: int = 4,
) -> dict[str, Any]:
    """Map finalized rows, globally cluster components, cap families, and write artifacts."""
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    leakage = protocol["leakage_gates"]
    gates = protocol["diversity_and_size_gate"]
    frozen_clustering = {
        "mmseqs_min_identity": 0.3,
        "mmseqs_coverage": 0.8,
        "mmseqs_coverage_mode": 0,
    }
    if any(leakage.get(key) != value for key, value in frozen_clustering.items()):
        raise ValueError("protocol does not specify exact global clustering at 30%/80%/mode 0")
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "global-accepted-components.fasta"
    manifest_path = output_dir / "global-accepted-mapping.csv"
    blockers: list[str] = []
    records: list[dict[str, Any]] = []
    sequences: list[dict[str, str]] = []
    source_summaries: list[dict[str, Any]] = []

    for source in finalized_accepted_sources(root, protocol_path):
        source_id = source["source_id"]
        columns = SOURCE_COMPONENT_COLUMNS.get(source_id)
        source_blockers: list[str] = []
        fasta_file = source["source_dir"] / SOURCE_FASTA_FILES.get(
            source_id, "construct_sequences.fasta"
        )
        if columns is None:
            source_blockers.append("no explicit source mapping schema")
        if not fasta_file.is_file():
            source_blockers.append("missing construct_sequences.fasta")
        fasta: dict[str, str] = {}
        if not source_blockers:
            try:
                fasta = _read_fasta(fasta_file)
                audit = json.loads(source["audit_path"].read_text(encoding="utf-8"))
                expected_hash = audit.get("construct_sequences_sha256")
                if source_id == "biostudies-S-EPMC12645435":
                    expected_hash = audit.get("variant_sequences", {}).get("sha256")
                if not isinstance(expected_hash, str) or len(expected_hash) != 64:
                    raise ValueError("finalized audit lacks construct_sequences_sha256")
                if expected_hash != _sha256(fasta_file):
                    raise ValueError("construct FASTA SHA256 differs from finalized audit")
            except (OSError, ValueError, json.JSONDecodeError) as error:
                source_blockers.append(str(error))

        pending_records: list[dict[str, Any]] = []
        pending_sequences: list[dict[str, str]] = []
        if not source_blockers:
            for row_index, row in enumerate(source["rows"], 1):
                candidate_id = (row.get("candidate_id") or row.get("source_record_id") or "").strip()
                try:
                    if not candidate_id:
                        raise ValueError("missing candidate_id")
                    ids = _component_ids(row, columns or (), source_id=source_id)
                    missing = [identifier for identifier in ids if identifier not in fasta]
                    if missing:
                        raise ValueError(f"FASTA has no exact ID(s): {','.join(missing)}")
                    global_ids = []
                    for component_index, identifier in enumerate(ids, 1):
                        global_id = f"gseq-{len(sequences) + len(pending_sequences) + 1:06d}"
                        global_ids.append(global_id)
                        pending_sequences.append(
                            {
                                "global_sequence_id": global_id,
                                "source_id": source_id,
                                "candidate_id": candidate_id,
                                "component_index": str(component_index),
                                "source_sequence_id": identifier,
                                "sequence": fasta[identifier],
                                "sequence_sha256": hashlib.sha256(fasta[identifier].encode("ascii")).hexdigest(),
                            }
                        )
                    pending_records.append(
                        {
                            "source_id": source_id,
                            "source_row_index": row_index,
                            "candidate_id": candidate_id,
                            "substrate_id": _substrate_id(row),
                            "global_sequence_ids": global_ids,
                        }
                    )
                except ValueError as error:
                    source_blockers.append(f"{candidate_id or 'row ' + str(row_index)}: {error}")
                    break
        if source_blockers:
            blockers.extend(f"{source_id}: {reason}" for reason in source_blockers)
        else:
            records.extend(pending_records)
            sequences.extend(pending_sequences)
        source_summaries.append(
            {
                "source_id": source_id,
                "finalized_accepted_records": len(source["rows"]),
                "mapped_records": 0 if source_blockers else len(pending_records),
                "mapping_columns": list(columns or ()),
                "blockers": source_blockers,
            }
        )

    with fasta_path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence in sequences:
            handle.write(f">{sequence['global_sequence_id']}\n{sequence['sequence']}\n")

    version: str | None = None
    actual_cluster_path = cluster_tsv
    if sequences and actual_cluster_path is None and mmseqs:
        try:
            version = _verify_mmseqs(mmseqs)
            mmseqs_output = output_dir / "mmseqs"
            _reset_mmseqs_workspace(mmseqs_output, mmseqs)
            actual_cluster_path = run_mmseqs_easy_cluster(
                mmseqs, fasta_path, mmseqs_output,
                min_identity=leakage["mmseqs_min_identity"],
                coverage=leakage["mmseqs_coverage"], threads=threads,
            )
        except (OSError, RuntimeError, ValueError) as error:
            blockers.append(f"global MMseqs2 clustering unavailable: {error}")
    elif actual_cluster_path is None:
        blockers.append("global MMseqs2 clustering not run; provide --mmseqs or --cluster-tsv")

    clusters: dict[str, str] = {}
    if actual_cluster_path is not None:
        try:
            clusters = read_mmseqs_clusters(actual_cluster_path)
            expected = {sequence["global_sequence_id"] for sequence in sequences}
            if set(clusters) != expected:
                missing = expected - set(clusters)
                extra = set(clusters) - expected
                raise ValueError(
                    f"cluster membership mismatch: {len(missing)} missing, {len(extra)} extra"
                )
        except (OSError, ValueError) as error:
            blockers.append(f"invalid global cluster TSV: {error}")
            clusters = {}

    family_counts: dict[str, int] = defaultdict(int)
    retained = 0
    if clusters:
        families = _record_families(records, clusters)
        cap = gates["maximum_records_per_family"]
        for index, record in enumerate(records):
            record["global_family_id"] = families[index]
            record["retained_after_family_cap"] = family_counts[families[index]] < cap
            family_counts[families[index]] += 1
            retained += int(record["retained_after_family_cap"])

    sequence_by_id = {sequence["global_sequence_id"]: sequence for sequence in sequences}
    fields = [
        "source_id", "source_row_index", "candidate_id", "substrate_id",
        "global_family_id", "retained_after_family_cap", "component_index",
        "source_sequence_id", "global_sequence_id", "sequence_sha256",
        "mmseqs_representative",
    ]
    with manifest_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            for global_id in record["global_sequence_ids"]:
                sequence = sequence_by_id[global_id]
                writer.writerow(
                    {
                        **{field: record.get(field, sequence.get(field, "")) for field in fields},
                        "global_sequence_id": global_id,
                        "mmseqs_representative": clusters.get(global_id, ""),
                    }
                )

    substrates = {
        record["substrate_id"] for record in records if record.get("retained_after_family_cap")
    }
    counts = {
        "accepted_sources": sum(not source["blockers"] for source in source_summaries),
        "excluded_or_unresolved_sources": sum(bool(source["blockers"]) for source in source_summaries),
        "mapped_records_before_cap": len(records),
        "mapped_component_sequences": len(sequences),
        "global_families": len(family_counts) if clusters else None,
        "records_after_global_family_cap": retained if clusters else None,
        "substrates_after_global_family_cap": len(substrates) if clusters else None,
    }
    if blockers:
        blockers.append("one or more finalized accepted sources could not be audited exactly")
    if clusters:
        if retained < gates["minimum_records"]:
            blockers.append(f"capped records {retained} < required {gates['minimum_records']}")
        if len(family_counts) < gates["minimum_mmseqs_families"]:
            blockers.append(
                f"global families {len(family_counts)} < required {gates['minimum_mmseqs_families']}"
            )
        if len(substrates) < gates["minimum_unique_substrates"]:
            blockers.append(f"capped substrates {len(substrates)} < required {gates['minimum_unique_substrates']}")

    report = {
        "schema_version": 1,
        "generated_on": date.today().isoformat(),
        "claim_boundary": "Readiness audit only; no model predictions or scores were generated.",
        "parameters": {
            "min_identity": leakage["mmseqs_min_identity"],
            "coverage": leakage["mmseqs_coverage"],
            "coverage_mode": leakage["mmseqs_coverage_mode"],
            "maximum_records_per_family": gates["maximum_records_per_family"],
            "cap_order": "source_id lexical order, then accepted CSV row order",
            "cap_uses_labels": False,
            "multicomponent_rule": "records are connected when any components share a global MMseqs cluster",
        },
        "mmseqs": {
            "executable": mmseqs,
            "required_version": MMSEQS_VERSION,
            "observed_version": version,
            "cluster_tsv": str(actual_cluster_path) if actual_cluster_path else None,
            "cluster_tsv_sha256": _sha256(actual_cluster_path) if clusters and actual_cluster_path else None,
        },
        "artifacts": {
            "global_fasta": str(fasta_path),
            "global_fasta_sha256": _sha256(fasta_path),
            "mapping_manifest": str(manifest_path),
            "mapping_manifest_sha256": _sha256(manifest_path),
        },
        "counts": counts,
        "family_record_counts_before_cap": dict(sorted(family_counts.items())),
        "readiness_gate_passes": not blockers,
        "blockers": blockers,
        "sources": source_summaries,
    }
    (output_dir / "readiness-audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("artifacts/external/temporal-absolute-kinetics"))
    parser.add_argument("--protocol", type=Path, default=Path("configs/temporal_absolute_kinetics_protocol.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/external/temporal-global-family-audit"))
    parser.add_argument("--mmseqs", default=DEFAULT_MMSEQS)
    parser.add_argument("--cluster-tsv", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    report = build_global_audit(
        args.root, args.protocol, args.output_dir, mmseqs=args.mmseqs,
        cluster_tsv=args.cluster_tsv, threads=args.threads,
    )
    return 0 if report["readiness_gate_passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
