from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from .schema import EnzymeSubstrateRecord


def sequence_id(sequence: str) -> str:
    return "seq_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()[:24]


def write_unique_fasta(records: Sequence[EnzymeSubstrateRecord], path: str | Path) -> dict[str, str]:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    sequences = {sequence_id(record.sequence): record.sequence for record in records}
    with target.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sorted(sequences.items()):
            handle.write(f">{identifier}\n{sequence}\n")
    return sequences


def run_mmseqs_easy_cluster(
    executable: str | Path,
    fasta_path: str | Path,
    output_dir: str | Path,
    *,
    min_identity: float = 0.3,
    coverage: float = 0.8,
    threads: int = 4,
) -> Path:
    executable_text = str(executable)
    wsl_parts = executable_text.split(":", 2) if executable_text.startswith("wsl:") else None
    if wsl_parts:
        if len(wsl_parts) != 3:
            raise ValueError("WSL MMseqs2 must use wsl:<distribution>:<command>")
        _, distribution, command_name = wsl_parts

        def tool_path(path: Path) -> str:
            resolved = path.resolve()
            drive = resolved.drive.rstrip(":").lower()
            if not drive:
                raise ValueError(f"cannot map path to WSL: {resolved}")
            relative = resolved.as_posix()[3:]
            return f"/mnt/{drive}/{relative}"

        command_prefix = ["wsl", "-d", distribution, command_name]
    else:
        executable_path = Path(executable_text)
        if not executable_path.is_file():
            raise FileNotFoundError(f"MMseqs2 executable not found: {executable_path}")

        def tool_path(path: Path) -> str:
            return str(path.resolve())

        command_prefix = [str(executable_path)]
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    database = output / "input_db"
    clusters = output / "clusters_db"
    temporary = output / "cluster_tmp"
    cluster_path = output / "proteins_cluster.tsv"
    commands = [
        command_prefix + ["createdb", tool_path(Path(fasta_path)), tool_path(database)],
        command_prefix + [
            "linclust", tool_path(database), tool_path(clusters),
            tool_path(temporary), "--min-seq-id", str(min_identity), "-c", str(coverage),
            "--cov-mode", "0", "--threads", str(threads),
        ],
        command_prefix + [
            "createtsv", tool_path(database), tool_path(database),
            tool_path(clusters), tool_path(cluster_path), "--threads", str(threads),
        ],
    ]
    for command in commands:
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            raise RuntimeError(
                f"MMseqs2 command failed with exit code {completed.returncode}: "
                f"{completed.stderr[-2000:]}")
    if not cluster_path.is_file():
        raise RuntimeError(f"MMseqs2 did not create expected cluster file: {cluster_path}")
    return cluster_path


def run_mmseqs_easy_search(
    executable: str | Path,
    query_fasta: str | Path,
    target_fasta: str | Path,
    output_dir: str | Path,
    *,
    min_identity: float = 0.3,
    coverage: float = 0.8,
    threads: int = 4,
) -> Path:
    executable_text = str(executable)
    wsl_parts = executable_text.split(":", 2) if executable_text.startswith("wsl:") else None
    if wsl_parts:
        if len(wsl_parts) != 3:
            raise ValueError("WSL MMseqs2 must use wsl:<distribution>:<command>")
        _, distribution, command_name = wsl_parts

        def tool_path(path: Path) -> str:
            resolved = path.resolve()
            drive = resolved.drive.rstrip(":").lower()
            if not drive:
                raise ValueError(f"cannot map path to WSL: {resolved}")
            return f"/mnt/{drive}/{resolved.as_posix()[3:]}"

        command_prefix = ["wsl", "-d", distribution, command_name]
    else:
        executable_path = Path(executable_text)
        if not executable_path.is_file():
            raise FileNotFoundError(f"MMseqs2 executable not found: {executable_path}")

        def tool_path(path: Path) -> str:
            return str(path.resolve())

        command_prefix = [str(executable_path)]
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    result_path = output / "homology_hits.tsv"
    temporary = output / "search_tmp"
    command = command_prefix + [
        "easy-search", tool_path(Path(query_fasta)), tool_path(Path(target_fasta)),
        tool_path(result_path), tool_path(temporary),
        "--min-seq-id", str(min_identity), "-c", str(coverage), "--cov-mode", "0",
        "--max-seqs", "1", "--threads", str(threads),
        "--format-output", "query,target,fident,qcov,tcov,evalue,bits",
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(
            f"MMseqs2 command failed with exit code {completed.returncode}: "
            f"{completed.stderr[-2000:]}")
    if not result_path.is_file():
        raise RuntimeError(f"MMseqs2 did not create expected search file: {result_path}")
    return result_path


def read_mmseqs_clusters(path: str | Path) -> dict[str, str]:
    member_to_representative: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"invalid MMseqs2 TSV at line {line_number}")
            representative, member = parts[:2]
            previous = member_to_representative.setdefault(member, representative)
            if previous != representative:
                raise ValueError(f"member {member} appears in multiple clusters")
    return member_to_representative


def assign_homology_splits(
    records: Sequence[EnzymeSubstrateRecord],
    member_to_representative: Mapping[str, str],
    *,
    ratios: tuple[float, float, float] = (0.8, 0.1, 0.1),
    seed: int = 42,
) -> tuple[EnzymeSubstrateRecord, ...]:
    if abs(sum(ratios) - 1.0) > 1e-9 or any(value < 0 for value in ratios):
        raise ValueError("ratios must be non-negative and sum to one")
    names = ("train", "validation", "test")
    groups: dict[str, list[int]] = {}
    for index, record in enumerate(records):
        identifier = sequence_id(record.sequence)
        if identifier not in member_to_representative:
            raise ValueError(f"sequence missing from MMseqs2 output: {identifier}")
        groups.setdefault(member_to_representative[identifier], []).append(index)
    ordered = sorted(
        groups.items(),
        key=lambda item: (
            -len(item[1]),
            hashlib.sha256(f"{seed}:{item[0]}".encode("ascii")).hexdigest(),
        ),
    )
    targets = [ratio * len(records) for ratio in ratios]
    counts = [0, 0, 0]
    assignments = [""] * len(records)
    for _, indices in ordered:
        candidates = [index for index, ratio in enumerate(ratios) if ratio > 0]
        destination = min(
            candidates,
            key=lambda index: (
                counts[index] / targets[index] if targets[index] else float("inf"),
                counts[index],
                index,
            ),
        )
        for index in indices:
            assignments[index] = names[destination]
        counts[destination] += len(indices)
    return tuple(replace(record, split=assignments[index]) for index, record in enumerate(records))


def write_split_manifest(records: Sequence[EnzymeSubstrateRecord], clusters: Mapping[str, str],
                         path: str | Path, *, parameters: dict) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    split_counts: dict[str, int] = {}
    rows = []
    for record in records:
        split_counts[record.split or "unset"] = split_counts.get(record.split or "unset", 0) + 1
        identifier = sequence_id(record.sequence)
        rows.append({
            "source_row": record.source_row,
            "sequence_id": identifier,
            "cluster": clusters[identifier],
            "split": record.split,
        })
    payload = {
        "format_version": 1,
        "parameters": parameters,
        "cluster_count": len(set(clusters.values())),
        "split_counts": split_counts,
        "rows": rows,
    }
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def apply_split_manifest(records: Sequence[EnzymeSubstrateRecord], path: str | Path,
                         *, source_identity: Mapping[str, Any] | None = None):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("format_version") != 1:
        raise ValueError("unsupported split manifest format")
    declared = payload.get("parameters", {}).get("source_identity")
    if source_identity is not None:
        if declared is not None and declared != dict(source_identity):
            raise ValueError(
                "split manifest source identity does not match the supplied data source")
    elif declared is not None:
        raise ValueError(
            "split manifest declares a source_identity; pass the matching data source "
            "identity instead of relying on the raw manifest")
    assignments: dict[int, str] = {}
    allowed = {"train", "validation", "test"}
    for row in payload["rows"]:
        source_row = int(row["source_row"])
        if source_row in assignments:
            raise ValueError(f"split manifest contains duplicate source_row {source_row}")
        split = row["split"]
        if split not in allowed:
            raise ValueError(f"split manifest has invalid split {split!r} at source_row {source_row}")
        assignments[source_row] = split
    record_rows = {record.source_row for record in records}
    unknown = sorted(assignments.keys() - record_rows)
    if unknown:
        raise ValueError(
            f"split manifest assigns {len(unknown)} source rows absent from the records")
    missing = [record.source_row for record in records if record.source_row not in assignments]
    if missing:
        raise ValueError(f"split manifest is missing {len(missing)} accepted source rows")
    return tuple(replace(record, split=assignments[record.source_row]) for record in records)
