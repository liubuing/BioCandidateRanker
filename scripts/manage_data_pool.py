"""Inventory storage and maintain a bounded, incremental evidence index without moving files."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

POOL = "artifacts/external/temporal-absolute-kinetics"
OUTPUT = "artifacts/data-pool-management"
CACHE_PARTS = {"pip-cache", ".pytest_cache", ".ruff_cache", "__pycache__", "cluster_tmp"}
ENV_PARTS = {"venv", ".venv", "site-packages"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def encode(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n"


def source_id(path: str) -> str | None:
    parts = Path(path).as_posix().split("/")
    prefix = POOL.split("/")
    return parts[len(prefix)] if parts[:len(prefix)] == prefix and len(parts) > 4 else None


def storage_kind(path: str) -> str:
    parts = set(Path(path).parts)
    if parts & CACHE_PARTS:
        return "cache"
    if parts & ENV_PARTS:
        return "environment"
    if path.startswith("artifacts/external/") or path.startswith("data/"):
        return "evidence"
    return "research_output"


def scan(root: Path, output: Path, previous: dict, budget: int, verify: bool) -> tuple:
    files = []
    errors = []
    for folder in ("artifacts", "data"):
        if not (root / folder).exists():
            continue
        for parent, dirs, names in os.walk(root / folder, followlinks=False,
                                         onerror=lambda exc: errors.append({
                                             "path": str(exc.filename), "error": str(exc)})):
            descend = []
            for directory in sorted(dirs):
                child = Path(parent) / directory
                if child == output:
                    continue
                try:
                    info = child.lstat()
                    if child.is_symlink() or getattr(info, "st_file_attributes", 0) & 0x400:
                        names.append(directory)
                    else:
                        descend.append(directory)
                except OSError as exc:
                    errors.append({"path": str(child), "error": str(exc)})
            dirs[:] = descend
            for name in sorted(names):
                path = Path(parent) / name
                try:
                    stat = path.lstat()
                except OSError as exc:
                    errors.append({"path": str(path), "error": str(exc)})
                    continue
                relative = path.relative_to(root).as_posix()
                # Linux links on NTFS may not satisfy is_symlink() on Windows.
                # Keep their metadata without attempting to read a foreign link target.
                is_link = path.is_symlink() or bool(getattr(stat, "st_file_attributes", 0) & 0x400)
                old = previous.get(relative)
                unchanged = old and old["size"] == stat.st_size and (
                    old["mtime_ns"] == stat.st_mtime_ns
                )
                kind = "link" if is_link else storage_kind(relative)
                files.append({"path": relative, "size": stat.st_size,
                              "mtime_ns": stat.st_mtime_ns,
                              "kind": kind,
                              "source_id": source_id(relative),
                              "sha256": old["sha256"] if unchanged and not verify
                              and kind == "evidence" else None})
    # Pool evidence gets priority; large external tables can be deferred by the byte budget.
    pending = sorted((f for f in files if f["kind"] == "evidence" and not f["sha256"]),
                     key=lambda f: (not f["path"].startswith(POOL + "/"), f["size"], f["path"]))
    hashed = reused = spent = 0
    for item in pending:
        if spent + item["size"] > budget:
            continue
        path = root / item["path"]
        try:
            value = digest(path)
            after = path.stat()
            if (after.st_size, after.st_mtime_ns) != (item["size"], item["mtime_ns"]):
                raise OSError("file changed during hashing; retry next run")
            item["sha256"] = value
            hashed += 1
            spent += item["size"]
        except OSError as exc:
            errors.append({"path": item["path"], "error": str(exc)})
    reused = sum(bool(f["sha256"]) for f in files) - hashed
    return sorted(files, key=lambda f: f["path"]), errors, hashed, reused, spent


def manage(root: Path, budget: int, verify: bool = False, reviewed: tuple = ()) -> dict:
    root = root.resolve()
    output = root / OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    old_report = json.loads((output / "report.json").read_text(encoding="utf-8")) if (
        output / "report.json"
    ).exists() else {}
    db = sqlite3.connect(output / "index.sqlite")
    db.row_factory = sqlite3.Row
    try:
        db.execute("CREATE TABLE IF NOT EXISTS files (path TEXT PRIMARY KEY, size INTEGER, "
                   "mtime_ns INTEGER, kind TEXT, source_id TEXT, sha256 TEXT)")
        previous = {row["path"]: dict(row) for row in db.execute("SELECT * FROM files")}
        files, errors, hashed, reused, spent = scan(root, output, previous, budget, verify)
        by_source = defaultdict(list)
        for item in files:
            if item["source_id"]:
                by_source[item["source_id"]].append(item)

        # Reconcile the current CSVs and audit contracts through the existing read-only registry.
        sys.path.insert(0, str(root / "src"))
        from biocandidate.temporal_registry import build_registry

        registry = build_registry(root / POOL,
                                  root / "configs/temporal_absolute_kinetics_protocol.json")
        dispositions = {s["source_id"]: s for s in registry["sources"]}
        previous_sources = {s["source_id"]: s for s in old_report.get("sources", [])}
        sources = []
        for name, entries in sorted(by_source.items()):
            contract = dispositions.get(name, {})
            accepted = contract.get("accepted_evidence_final") and (
                contract.get("csv_accepted_records", 0) > 0
            )
            terminal = contract.get("disposition") == "finalized_zero_candidate_source"
            layer = "accepted_source_evidence" if accepted else (
                "archive_candidate" if terminal else "pending"
            )
            changed = any(previous.get(f["path"], {}).get("sha256") != f["sha256"]
                          or previous.get(f["path"], {}).get("size") != f["size"]
                          or previous.get(f["path"], {}).get("mtime_ns") != f["mtime_ns"]
                          for f in entries)
            removed = any(f["source_id"] == name and p not in {e["path"] for e in entries}
                          for p, f in previous.items())
            evidence_changed = any(
                previous.get(f["path"], {}).get("sha256") != f["sha256"]
                or previous.get(f["path"], {}).get("size") != f["size"]
                for f in entries if f["kind"] == "evidence"
            )
            present_paths = {e["path"] for e in entries}
            evidence_removed = any(f["source_id"] == name and f["kind"] == "evidence"
                                   and p not in present_paths for p, f in previous.items())
            sources.append({"source_id": name, "layer": layer, "files": len(entries),
                            "bytes": sum(f["size"] for f in entries),
                            "hashes_complete": all(f["sha256"] for f in entries
                                                   if f["kind"] == "evidence"),
                            "change": "new" if not any(f["source_id"] == name
                                                       for f in previous.values()) else (
                                "changed" if changed or removed else "unchanged"),
                            "accepted_records": contract.get("csv_accepted_records", 0),
                            "disposition": contract.get("disposition", "unregistered"),
                            "issues": contract.get("issues", [])})
            sources[-1]["review_required"] = bool(
                (old_report.get("schema_version") == 2
                 and previous_sources.get(name, {}).get("review_required"))
                or (name in previous_sources and (evidence_changed or evidence_removed))
                or layer == "pending"
            )
        removed_paths = sorted(set(previous) - {f["path"] for f in files})
        groups = defaultdict(lambda: {"files": 0, "bytes": 0})
        duplicates = defaultdict(list)
        for item in files:
            groups[item["kind"]]["files"] += 1
            groups[item["kind"]]["bytes"] += item["size"]
            if item["kind"] == "evidence" and item["sha256"] and item["size"]:
                duplicates[(item["sha256"], item["size"])].append(item["path"])
        duplicate_groups = [{"sha256": sha, "size": size, "paths": paths,
                             "redundant_bytes": size * (len(paths) - 1)}
                            for (sha, size), paths in duplicates.items() if len(paths) > 1]
        duplicate_groups.sort(key=lambda g: (-g["redundant_bytes"], g["sha256"]))
        # Freeze durable regular evidence files; cache and foreign links remain listed exclusions.
        pool_files = [{"path": f["path"], "size": f["size"], "sha256": f["sha256"]}
                      for f in files if f["path"].startswith(POOL + "/")
                      and f["kind"] == "evidence"]
        exclusions = [{"path": f["path"], "kind": f["kind"]}
                      for f in files if f["path"].startswith(POOL + "/")
                      and f["kind"] != "evidence"]
        controls = []
        for name in ("configs/temporal_absolute_kinetics_protocol.json",
                     "artifacts/external/temporal-global-family-audit/readiness-audit.json",
                     "configs/software_implementation_release_v6.json"):
            controls.append({"path": name, "sha256": digest(root / name)})
        controls_changed = bool(old_report) and old_report.get("controls") != controls
        removed_sources = sorted(({f["source_id"] for f in previous.values()
                                   if f["source_id"]}
                                  | set(old_report.get("removed_sources", []))) - set(by_source))
        unknown_reviews = set(reviewed) - set(by_source)
        if unknown_reviews:
            raise ValueError(f"cannot acknowledge missing sources: {sorted(unknown_reviews)}")
        for source in sources:
            if controls_changed:
                source["review_required"] = True
            if source["source_id"] in reviewed:
                if source["layer"] == "pending" or not source["hashes_complete"] or errors:
                    raise ValueError(f"cannot acknowledge unresolved source: {source['source_id']}")
                source["review_required"] = False
        snapshot = None
        if pool_files and all(f["sha256"] for f in pool_files) and not errors:
            payload = {"schema_version": 1, "files": pool_files, "controls": controls,
                       "excluded_generated_files_and_links": exclusions,
                       "scope": "Regular evidence file identity only; generated caches and "
                                "foreign links excluded. Does not authorize scoring or training."}
            identity = hashlib.sha256(encode(payload).encode()).hexdigest()
            target = output / "snapshots" / f"{identity}.json"
            target.parent.mkdir(exist_ok=True)
            if not target.exists():
                with target.open("x", encoding="utf-8") as handle:
                    handle.write(encode(payload))
            elif target.read_text(encoding="utf-8") != encode(payload):
                raise ValueError("existing snapshot content mismatch")
            snapshot = target.relative_to(root).as_posix()
        report = {"schema_version": 2,
                  "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                  "total_files": len(files), "total_bytes": sum(f["size"] for f in files),
                  "storage": dict(groups), "hashed_this_run": hashed, "hashes_reused": reused,
                  "hashed_bytes_this_run": spent,
                  "pending_hashes": sum(f["kind"] == "evidence" and not f["sha256"]
                                        for f in files),
                  "removed_paths": removed_paths, "sources": sources,
                  "removed_sources": removed_sources, "controls": controls,
                  "controls_changed": controls_changed,
                  "review_queue": [s["source_id"] for s in sources
                                   if s["change"] != "unchanged" or s["layer"] == "pending"
                                   or not s["hashes_complete"] or s["review_required"]]
                                  + removed_sources,
                  "duplicate_groups": duplicate_groups,
                  "links": [f["path"] for f in files if f["kind"] == "link"],
                  "largest_files": sorted(files, key=lambda f: -f["size"])[:20],
                  "evidence_snapshot": snapshot, "registry_counts": registry["counts"],
                  "errors": errors,
                  "claim_boundary": "Index only. No files moved/deleted; no measurements "
                                    "merged; no new scientific admission or predictions."}
        with db:
            db.execute("DELETE FROM files")
            db.executemany("INSERT INTO files VALUES (:path,:size,:mtime_ns,:kind,:source_id,:sha256)",
                           files)
        (output / "report.json").write_text(encode(report), encoding="utf-8")
        (output / "registry-current.json").write_text(encode(registry), encoding="utf-8")
        lines = ["# 数据池管理报告", "",
                 f"文件：{len(files):,}；总占用：{report['total_bytes'] / 1024**3:.2f} GiB。", "",
                 "| 存储类别 | 文件数 | GiB |", "|---|---:|---:|"]
        lines += [f"| {kind} | {g['files']:,} | {g['bytes'] / 1024**3:.2f} |"
                  for kind, g in sorted(groups.items())]
        lines += ["", f"本轮计算 {hashed} 个哈希，复用 {reused} 个；"
                  f"仍待计算 {report['pending_hashes']} 个。", "",
                  f"当前注册表接受记录：{registry['counts']['accepted_records']}。",
                  f"证据快照：{snapshot or '未完成，不能视为冻结版本'}。", "",
                  "层级标记针对来源目录；接受来源内仍可能包含被排除的行。",
                  "缓存可在相应任务停止后清理；环境、权重及原始证据保留。",
                  "相同哈希仅表示文件内容相同，不代表测量重复，也不自动删除。", "",
                  f"单独记录 {len(report['links'])} 个链接；临时缓存与链接不纳入证据快照。", "",
                  "## 来源索引", "", "| 来源 | 层级 | 变化 | MiB |",
                  "|---|---|---|---:|"]
        lines += [f"| {s['source_id']} | {s['layer']} | {s['change']} | "
                  f"{s['bytes'] / 1024**2:.2f} |" for s in sources]
        lines += ["", "## 归档候选", "",
                  "以下目录保留原路径和全部证据，仅列为归档候选；未移动或删除。", ""]
        lines += [f"- `{POOL}/{s['source_id']}`" for s in sources
                  if s["layer"] == "archive_candidate"]
        (output / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return report
    finally:
        db.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--hash-budget-mib", type=int, default=1024)
    parser.add_argument("--verify-hashes", action="store_true",
                        help="Rehash evidence instead of trusting unchanged size/mtime")
    parser.add_argument("--reviewed-source", action="append", default=[],
                        help="Acknowledge a source ONLY after completing its evidence review")
    args = parser.parse_args()
    if args.hash_budget_mib < 0:
        parser.error("hash budget must be nonnegative")
    report = manage(args.root, args.hash_budget_mib * 1024**2, args.verify_hashes,
                    tuple(args.reviewed_source))
    summary = {k: report[k] for k in ("total_files", "total_bytes", "storage",
                                    "hashed_this_run", "hashes_reused", "pending_hashes",
                                    "evidence_snapshot", "registry_counts")}
    summary["error_count"] = len(report["errors"])
    print(encode(summary))
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
