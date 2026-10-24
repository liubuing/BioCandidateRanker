"""MMseqs2 homology clustering and cold-split assignment.

Sequences are clustered at the governed 30% identity / 80% coverage thresholds
(``--cov-mode 0`` requires both query and target to be covered), and split
assignment operates on whole clusters only, so no homolog pair is ever split
across train/validation/test. Split manifests record the MMseqs2 parameters
and the data source identity so runs can replay or reject stale manifests.

中文：MMseqs2 同源聚类与同源冷切分分配。序列按受治理的 30% 序列一致性 / 80%
覆盖率阈值（--cov-mode 0 要求 query 与 target 双侧都被覆盖）聚类，切分只以整个
聚类为单位分配，因此任何同源对都不会被拆到 train/validation/test 两侧。切分
清单记录 MMseqs2 参数与数据源身份，运行方可重放或拒绝过期清单。
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from .schema import EnzymeSubstrateRecord


def sequence_id(sequence: str) -> str:
    """Stable content address ("seq_" + first 24 hex chars of the SHA256 digest).
    中文：稳定的内容寻址（"seq_" + SHA256 摘要前 24 个十六进制字符）。
    """
    return "seq_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()[:24]


def write_unique_fasta(records: Sequence[EnzymeSubstrateRecord], path: str | Path) -> dict[str, str]:
    """Write deduplicated sequences to FASTA; returns the identifier-to-sequence mapping.
    中文：把去重后的序列写入 FASTA；返回标识符到序列的映射。
    """
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
    """Cluster the FASTA with MMseqs2 linclust and return the cluster TSV path.

    ``min_identity`` and ``coverage`` are the governed homology thresholds
    (0.3 / 0.8); coverage uses ``--cov-mode 0``, meaning both cluster members
    must cover that fraction of the alignment.
    中文：用 MMseqs2 linclust 对 FASTA 聚类并返回聚类 TSV 路径。min_identity 与
    coverage 是受治理的同源阈值（0.3 / 0.8）；覆盖率使用 --cov-mode 0，即两个
    聚类成员都必须覆盖该比例的比对区域。
    """
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
        completed = subprocess.run(command, check=False, capture_output=True, text=True, errors="replace")
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
    """Search a query FASTA against a target FASTA; returns the top-hit TSV path.
    中文：用 query FASTA 对 target FASTA 做搜索；返回 top-hit TSV 路径。
    """
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
    completed = subprocess.run(command, check=False, capture_output=True, text=True, errors="replace")
    if completed.returncode:
        raise RuntimeError(
            f"MMseqs2 command failed with exit code {completed.returncode}: "
            f"{completed.stderr[-2000:]}")
    if not result_path.is_file():
        raise RuntimeError(f"MMseqs2 did not create expected search file: {result_path}")
    return result_path


def read_mmseqs_clusters(path: str | Path) -> dict[str, str]:
    """Read a createtsv cluster TSV into a member-to-representative mapping.

    Raises if any member appears under two representatives, which would make
    cold-split assignment ambiguous.
    中文：把 createtsv 聚类 TSV 读成成员到代表序列的映射。若任一成员出现在两个
    代表之下则抛错，否则同源冷切分分配会产生歧义。
    """
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
    """Assign whole MMseqs2 clusters to train/validation/test by greedy fill.

    Clusters are processed largest-first with a seed-hashed tie-break, so the
    result is deterministic and a homology cluster never straddles two splits.
    中文：按贪心填充把整个 MMseqs2 聚类分配到 train/validation/test。聚类按
    从大到小处理并以种子哈希打破并列，因此结果确定，且同源聚类绝不横跨两个
    切分。
    """
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
    """Persist per-row split assignments plus the parameters that produced them.
    中文：持久化逐行切分分配及生成这些分配所用的参数。
    """
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
    """Replay split assignments from a manifest, verifying version and coverage.

    Fails closed on unsupported formats, duplicate or invalid rows, source rows
    absent from the records, accepted rows missing from the manifest, and any
    declared ``source_identity`` that does not match the supplied data source.
    中文：从清单重放切分分配，并校验版本与覆盖情况。对不支持的格式、重复或非法
    行、记录中不存在的源行、清单缺失的已接受行，以及声明的 source_identity 与
    所供数据源不符等情况一律 fail-closed（拒绝）。
    """
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
