from __future__ import annotations

import hashlib
import csv
import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Union


PathLike = Union[str, Path]


class ManifestVerificationError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class FileManifest:
    sha256: str
    size_bytes: int
    row_count: int

    def __post_init__(self) -> None:
        digest = self.sha256.lower()
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("sha256 must be a 64-character hexadecimal digest")
        if self.size_bytes < 0 or self.row_count < 0:
            raise ValueError("size_bytes and row_count must be non-negative")
        object.__setattr__(self, "sha256", digest)


def _json_row_count(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ManifestVerificationError(f"cannot count rows in {path}: {exc}") from exc

    if isinstance(payload, list):
        return len(payload)
    if isinstance(payload, dict):
        for key in ("data", "records", "rows", "candidates", "candidate_ids"):
            if isinstance(payload.get(key), list):
                return len(payload[key])
        partitions = payload.get("partition_candidate_ids")
        if isinstance(partitions, dict) and all(isinstance(rows, list) for rows in partitions.values()):
            return sum(len(rows) for rows in partitions.values())
    raise ManifestVerificationError(
        "JSON row counting requires an array or a data/records/rows array"
    )


def _row_count(path: Path) -> int:
    if path.suffix.lower() in {".tsv", ".txt", ".xls"}:
        try:
            with path.open("r", encoding="utf-8") as handle:
                return sum(1 for line in handle if line.strip())
        except (OSError, UnicodeError) as exc:
            raise ManifestVerificationError(f"cannot count rows in {path}: {exc}") from exc
    if path.suffix.lower() == ".zip":
        try:
            with zipfile.ZipFile(path) as archive:
                rows = 0
                members = [
                    name for name in archive.namelist()
                    if "/experiments/" in name.replace("\\", "/")
                    and name.lower().endswith(".csv")
                    and not name.startswith("__MACOSX/")
                ]
                for name in members:
                    with archive.open(name) as raw:
                        handle = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
                        reader = csv.reader(handle)
                        next(reader, None)
                        rows += sum(1 for row in reader if any(cell.strip() for cell in row))
                return rows
        except (OSError, UnicodeError, csv.Error, zipfile.BadZipFile) as exc:
            raise ManifestVerificationError(f"cannot count rows in {path}: {exc}") from exc
    if path.suffix.lower() == ".csv":
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.reader(handle)
                next(reader, None)
                return sum(1 for row in reader if any(cell.strip() for cell in row))
        except (OSError, UnicodeError, csv.Error) as exc:
            raise ManifestVerificationError(f"cannot count rows in {path}: {exc}") from exc
    if path.suffix.lower() in {".jsonl", ".ndjson"}:
        try:
            with path.open("r", encoding="utf-8") as handle:
                return sum(1 for line in handle if line.strip())
        except (OSError, UnicodeError) as exc:
            raise ManifestVerificationError(f"cannot count rows in {path}: {exc}") from exc
    return _json_row_count(path)


def verify_manifest(path: PathLike, manifest: FileManifest) -> None:
    source = Path(path)
    if not source.is_file():
        raise ManifestVerificationError(f"source is not a regular file: {source}")

    digest = hashlib.sha256()
    size = 0
    try:
        with source.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                digest.update(chunk)
    except OSError as exc:
        raise ManifestVerificationError(f"cannot read source {source}: {exc}") from exc

    actual_digest = digest.hexdigest()
    rows = _row_count(source)
    mismatches = []
    if actual_digest != manifest.sha256:
        mismatches.append(f"sha256 expected {manifest.sha256}, got {actual_digest}")
    if size != manifest.size_bytes:
        mismatches.append(f"size expected {manifest.size_bytes}, got {size}")
    if rows != manifest.row_count:
        mismatches.append(f"rows expected {manifest.row_count}, got {rows}")
    if mismatches:
        raise ManifestVerificationError("manifest verification failed: " + "; ".join(mismatches))


def build_manifest(path: PathLike) -> FileManifest:
    """Compute an identity manifest without modifying the source asset."""
    source = Path(path)
    if not source.is_file():
        raise ManifestVerificationError(f"source is not a regular file: {source}")
    digest = hashlib.sha256()
    size = 0
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return FileManifest(digest.hexdigest(), size, _row_count(source))
