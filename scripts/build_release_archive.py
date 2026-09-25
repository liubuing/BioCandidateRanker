"""Build a deterministic release archive for public deposit.

The archive contains every Git-tracked file plus each frozen evidence artifact
named by the release manifest. Entries are written in sorted order with fixed
timestamps, so rebuilding the same commit reproduces byte-identical output.

A RELEASE_MANIFEST.json is embedded in the archive and written alongside it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_RELEASE_MANIFEST = Path("configs/software_implementation_release_v6.json")
DEFAULT_LICENSE_INVENTORY = Path("configs/data_license_inventory.json")
MANIFEST_MEMBER_NAME = "RELEASE_MANIFEST.json"
ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
EXCLUDED_PREFIXES = (".git/", ".zcode/")
EXCLUDED_SUFFIXES = (".pyc", ".pyo")
EXCLUDED_PARTS = ("__pycache__",)

ROLE_TRACKED = "tracked_source"
ROLE_FROZEN = "frozen_evidence"

# Absolute local paths that must not appear in a double-blind deposit.
LOCAL_PATH_PATTERNS = (
    re.compile(r"[A-Za-z]:\\\\?[^\"'\s,;)]{3,}"),
    re.compile(r"/(?:home|Users)/[A-Za-z0-9._-]+[^\"'\s,;)]*"),
)


def _run_git(root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return completed.stdout


def _excluded(relative: str) -> bool:
    if relative.startswith(EXCLUDED_PREFIXES):
        return True
    if relative.endswith(EXCLUDED_SUFFIXES):
        return True
    return any(part in EXCLUDED_PARTS for part in Path(relative).parts)


def list_tracked_files(root: Path) -> list[str]:
    """Return sorted repo-relative paths of every Git-tracked file."""
    output = _run_git(root, "ls-files", "-z")
    if not output:
        return []
    names = [name for name in output.split("\0") if name]
    return sorted(name for name in names if not _excluded(name.replace("\\", "/")))


def source_date(root: Path) -> str:
    """Derive a stable build timestamp from the commit, not the wall clock."""
    committed = _run_git(root, "log", "-1", "--format=%cI")
    if committed and committed.strip():
        return committed.strip()
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def head_commit(root: Path) -> str | None:
    commit = _run_git(root, "rev-parse", "HEAD")
    return commit.strip() if commit and commit.strip() else None


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(filename=name, date_time=ZIP_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 0
    info.external_attr = 0o644 << 16
    return info


def collect_entries(root: Path, release_manifest: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    """Resolve archive members and flag frozen evidence missing from disk."""
    entries: list[dict[str, Any]] = []
    seen: set[str] = set()
    for relative in list_tracked_files(root):
        path = root / relative
        if not path.is_file():
            continue
        entries.append({"path": relative, "role": ROLE_TRACKED})
        seen.add(relative)

    missing: list[str] = []
    for frozen in release_manifest.get("frozen_files", []):
        relative = frozen["path"].replace("\\", "/")
        if relative in seen:
            continue
        path = root / relative
        if not path.is_file():
            missing.append(relative)
            continue
        entries.append({"path": relative, "role": ROLE_FROZEN})
        seen.add(relative)

    entries.sort(key=lambda entry: entry["path"])
    return entries, sorted(missing)


def scan_local_paths(
    entries: list[dict[str, Any]],
    root: Path,
    max_matches_per_file: int = 5,
) -> list[dict[str, Any]]:
    """Find absolute local filesystem paths in text members.

    A double-blind deposit must not carry the authors' home or project paths, and
    frozen audit files frequently echo the machine that produced them.
    """
    findings: list[dict[str, Any]] = []
    for entry in entries:
        if entry["path"].endswith((".pt", ".pth", ".ckpt", ".lmdb", ".zip")):
            continue
        try:
            text = (root / entry["path"]).read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        matches: list[str] = []
        for pattern in LOCAL_PATH_PATTERNS:
            for match in pattern.finditer(text):
                if len(matches) >= max_matches_per_file:
                    break
                matches.append(match.group(0))
            if len(matches) >= max_matches_per_file:
                break
        if matches:
            findings.append({"path": entry["path"], "matches": matches})
    return findings


REDACTION = "<redacted-local-path>"


def redact_payload(payload: bytes) -> tuple[bytes, int]:
    """Replace absolute local paths in a text member with a neutral marker.

    Returns the transformed payload and the number of replacements. Binary
    members are returned untouched.
    """
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload, 0
    replaced = 0
    for pattern in LOCAL_PATH_PATTERNS:
        text, count = pattern.subn(REDACTION, text)
        replaced += count
    return text.encode("utf-8"), replaced


def build_redact_map(
    entries: list[dict[str, Any]],
    root: Path,
) -> tuple[dict[str, bytes], list[dict[str, Any]]]:
    """Redacted payload for every member that embeds a local path."""
    redact_map: dict[str, bytes] = {}
    redactions: list[dict[str, Any]] = []
    for entry in entries:
        if entry["path"].endswith((".pt", ".pth", ".ckpt", ".lmdb", ".zip")):
            continue
        payload = (root / entry["path"]).read_bytes()
        transformed, replaced = redact_payload(payload)
        if replaced:
            redact_map[entry["path"]] = transformed
            redactions.append({"path": entry["path"], "replaced": replaced})
    return redact_map, redactions


def member_bytes(entry: dict[str, Any], root: Path, redact_map: dict[str, bytes] | None) -> bytes:
    if redact_map and entry["path"] in redact_map:
        return redact_map[entry["path"]]
    return (root / entry["path"]).read_bytes()


def check_license_ownership(
    entries: list[dict[str, Any]],
    license_inventory: dict[str, Any] | None,
) -> tuple[list[str], list[dict[str, Any]]]:
    """Confirm every archived path belongs to a corpus that may ship.

    A frozen artifact carrying measurement data from a non-redistributable corpus
    must not enter a public archive, so ownership is checked before writing.
    """
    if license_inventory is None:
        return [], []

    owned_prefixes: list[tuple[str, dict[str, Any]]] = []
    for corpus in license_inventory.get("corpora", []):
        for prefix in corpus.get("release_paths", []):
            owned_prefixes.append((prefix, corpus))
        for derived in corpus.get("derived_artifacts_in_release", []):
            owned_prefixes.append((derived["path"], corpus))

    violations: list[str] = []
    ownership: list[dict[str, Any]] = []
    for entry in entries:
        path = entry["path"]
        match = None
        for prefix, corpus in owned_prefixes:
            if path == prefix or path.startswith(prefix):
                if match is None or len(prefix) > len(match[0]):
                    match = (prefix, corpus)
        if match is None:
            # Source, docs, and tests are the project's own work, not a corpus.
            # Only artifacts/ can carry third-party data, so only artifacts/
            # must be attributable.
            ownership.append(
                {
                    "path": path,
                    "corpus_id": "project-source" if not path.startswith("artifacts/") else None,
                }
            )
            continue
        corpus = match[1]
        ownership.append({"path": path, "corpus_id": corpus["id"]})
        blocked = corpus.get("redistribution") in {"not_verified", "blocked_pending_permission"}
        if blocked and not corpus.get("redistribute_in_release", False):
            derived = next(
                (
                    item for item in corpus.get("derived_artifacts_in_release", [])
                    if item["path"] == path
                ),
                None,
            )
            if derived and not derived.get("contains_measurements", False):
                continue
            violations.append(
                f"{path} belongs to {corpus['id']}, which may not ship "
                f"({corpus.get('blocking_reason', 'redistribution unresolved')})"
            )
    return violations, ownership


def build_manifest(
    root: Path,
    release_manifest_path: Path,
    release_manifest: dict[str, Any],
    entries: list[dict[str, Any]],
    missing: list[str],
    redact_map: dict[str, bytes] | None = None,
) -> dict[str, Any]:
    # Resolve against root, not the process working directory, so --root works
    # from anywhere.
    resolved_manifest = (
        release_manifest_path if release_manifest_path.is_absolute()
        else root / release_manifest_path
    )
    files: list[dict[str, Any]] = []
    total = 0
    for entry in entries:
        payload = member_bytes(entry, root, redact_map)
        total += len(payload)
        files.append(
            {
                "path": entry["path"],
                "role": entry["role"],
                "sha256": sha256_bytes(payload),
                "size_bytes": len(payload),
            }
        )

    # Every frozen path is checked, including ones that are also Git-tracked;
    # a tracked file that drifted from the freeze must not pass silently.
    frozen_index = {
        item["path"].replace("\\", "/"): item for item in release_manifest.get("frozen_files", [])
    }
    verified: list[dict[str, Any]] = []
    for record in files:
        expected = frozen_index.get(record["path"])
        if expected is None:
            continue
        verified.append(
            {
                "path": record["path"],
                "role": record["role"],
                "sha256_matches_manifest": record["sha256"] == expected["sha256"],
                "size_matches_manifest": record["size_bytes"] == expected["size_bytes"],
            }
        )

    return {
        "schema_version": 1,
        "release_id": release_manifest.get("release_id"),
        "package_version": release_manifest.get("package_version"),
        "built_from_commit": head_commit(root),
        "built_on": source_date(root),
        "builder": "scripts/build_release_archive.py",
        "source_manifest": release_manifest_path.as_posix(),
        "source_manifest_sha256": sha256_bytes(resolved_manifest.read_bytes()),
        "claim_boundary": release_manifest.get("claim_boundary"),
        "file_count": len(files),
        "total_size_bytes": total,
        "frozen_evidence_verified": verified,
        "missing_frozen_evidence": missing,
        "files": files,
    }


def write_archive(
    archive_path: Path,
    root: Path,
    entries: list[dict[str, Any]],
    manifest: dict[str, Any],
    redact_map: dict[str, bytes] | None = None,
) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for entry in entries:
            archive.writestr(_zip_info(entry["path"]), member_bytes(entry, root, redact_map))
        payload = json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8")
        archive.writestr(_zip_info(MANIFEST_MEMBER_NAME), payload)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=DEFAULT_RELEASE_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/release-archive"))
    parser.add_argument("--license-inventory", type=Path, default=DEFAULT_LICENSE_INVENTORY)
    parser.add_argument(
        "--skip-license-check",
        action="store_true",
        help="build without the corpus ownership check (local inspection only)",
    )
    parser.add_argument(
        "--fail-on-local-paths",
        action="store_true",
        help="refuse to deposit while any member embeds an absolute local path",
    )
    parser.add_argument(
        "--redact-paths",
        action="store_true",
        help="replace absolute local paths in the archived bytes; writes a "
             "-redacted suffix archive, leaving the integrity archive untouched",
    )
    parser.add_argument(
        "--require-frozen-evidence",
        action="store_true",
        help="fail instead of recording frozen evidence that is absent from disk",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    manifest_path = (root / args.manifest).resolve() if not args.manifest.is_absolute() else args.manifest
    if not manifest_path.is_file():
        print(json.dumps({"valid": False, "issues": [f"missing manifest: {manifest_path}"]}, indent=2))
        return 1

    release_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    release_id = release_manifest.get("release_id", "unreleased")

    inventory_path = (
        args.license_inventory if args.license_inventory.is_absolute()
        else root / args.license_inventory
    )
    license_inventory = (
        json.loads(inventory_path.read_text(encoding="utf-8"))
        if inventory_path.is_file()
        else None
    )
    if license_inventory is None and not args.skip_license_check:
        print(
            json.dumps(
                {
                    "valid": False,
                    "issues": [f"missing license inventory: {args.license_inventory.as_posix()}"],
                },
                indent=2,
            )
        )
        return 1

    entries, missing = collect_entries(root, release_manifest)
    if missing and args.require_frozen_evidence:
        print(
            json.dumps(
                {
                    "valid": False,
                    "issues": [f"missing frozen evidence: {path}" for path in missing],
                },
                indent=2,
            )
        )
        return 1

    violations, ownership = check_license_ownership(entries, license_inventory)
    if violations:
        print(
            json.dumps(
                {
                    "valid": False,
                    "issues": violations,
                    "remedy": (
                        "resolve the corpus terms, drop the artifact from the release "
                        "manifest, or record it as a derived artifact that carries no "
                        "measurements"
                    ),
                },
                indent=2,
            )
        )
        return 1

    redact_map: dict[str, bytes] | None = None
    redactions: list[dict[str, Any]] = []
    if args.redact_paths:
        redact_map, redactions = build_redact_map(entries, root)

    manifest = build_manifest(
        root, args.manifest, release_manifest, entries, missing, redact_map=redact_map
    )
    unowned = sorted(item["path"] for item in ownership if item["corpus_id"] is None)
    local_paths = scan_local_paths(entries, root)
    manifest["license_ownership"] = {
        "inventory": args.license_inventory.as_posix() if license_inventory else None,
        "unowned_paths": unowned,
        "by_corpus": {},
    }
    for item in ownership:
        if item["corpus_id"]:
            manifest["license_ownership"]["by_corpus"].setdefault(item["corpus_id"], []).append(
                item["path"]
            )
    manifest["local_path_findings"] = local_paths
    if args.redact_paths:
        # Post-redaction guard: the deposit bytes must be clean. Frozen-evidence
        # hash mismatches are expected for redacted members and the verifier
        # tolerates exactly those paths.
        residual = []
        for entry in entries:
            payload = member_bytes(entry, root, redact_map)
            try:
                text = payload.decode("utf-8")
            except UnicodeDecodeError:
                continue
            if any(pattern.search(text) for pattern in LOCAL_PATH_PATTERNS):
                residual.append(entry["path"])
        manifest["redacted"] = True
        manifest["redactions"] = redactions
        manifest["pre_redaction_findings_count"] = len(local_paths)
        manifest["local_path_findings"] = []
        manifest["residual_local_paths"] = residual
        local_paths = []
        if residual:
            print(
                json.dumps(
                    {
                        "valid": False,
                        "issues": [f"redaction left local paths in: {path}" for path in residual],
                    },
                    indent=2,
                )
            )
            return 1

    if local_paths and args.fail_on_local_paths:
        print(
            json.dumps(
                {
                    "valid": False,
                    "issues": [
                        f"{item['path']} embeds local paths: {', '.join(item['matches'])}"
                        for item in local_paths
                    ],
                    "remedy": "redact the field before deposit, or exclude the file",
                },
                indent=2,
            )
        )
        return 1

    if args.dry_run:
        print(
            json.dumps(
                {
                    "dry_run": True,
                    "release_id": release_id,
                    "file_count": manifest["file_count"],
                    "total_size_bytes": manifest["total_size_bytes"],
                    "missing_frozen_evidence": missing,
                },
                indent=2,
            )
        )
        return 0

    output_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    suffix = "-redacted" if args.redact_paths else ""
    archive_path = output_dir / f"{release_id}{suffix}.zip"
    sidecar_path = output_dir / f"{release_id}{suffix}-manifest.json"

    write_archive(archive_path, root, entries, manifest, redact_map)
    sidecar_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "valid": True,
                "release_id": release_id,
                "archive": archive_path.relative_to(root).as_posix(),
                "archive_sha256": sha256_bytes(archive_path.read_bytes()),
                "archive_size_bytes": archive_path.stat().st_size,
                "manifest": sidecar_path.relative_to(root).as_posix(),
                "file_count": manifest["file_count"],
                "missing_frozen_evidence": missing,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
