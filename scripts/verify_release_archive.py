"""Verify a release archive against its embedded manifest and the release freeze.

Checks performed:
  1. every archive member matches the size and SHA256 recorded in
     RELEASE_MANIFEST.json;
  2. no archive member is unlisted, and no listed file is absent;
  3. every frozen-evidence member still matches the release manifest's own
     SHA256 and size records;
  4. the archive records the source-manifest hash of the release it claims.

Exits non-zero when any check fails, so this can gate a deposit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

DEFAULT_RELEASE_MANIFEST = Path("configs/software_implementation_release_v6.json")
MANIFEST_MEMBER_NAME = "RELEASE_MANIFEST.json"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def verify_archive(
    archive_path: Path,
    release_manifest_path: Path | None,
) -> list[str]:
    issues: list[str] = []

    if not archive_path.is_file():
        return [f"missing archive: {archive_path}"]

    with zipfile.ZipFile(archive_path) as archive:
        names = archive.namelist()
        if MANIFEST_MEMBER_NAME not in names:
            return [f"archive has no {MANIFEST_MEMBER_NAME}"]

        try:
            embedded = json.loads(archive.read(MANIFEST_MEMBER_NAME).decode("utf-8"))
        except json.JSONDecodeError as error:
            return [f"{MANIFEST_MEMBER_NAME} is not valid JSON: {error}"]

        recorded = {entry["path"]: entry for entry in embedded.get("files", [])}

        actual_members = [name for name in names if name != MANIFEST_MEMBER_NAME]
        unlisted = sorted(set(actual_members) - set(recorded))
        for name in unlisted:
            issues.append(f"member not listed in manifest: {name}")

        absent = sorted(set(recorded) - set(actual_members))
        for name in absent:
            issues.append(f"listed file absent from archive: {name}")

        for name in sorted(set(actual_members) & set(recorded)):
            payload = archive.read(name)
            entry = recorded[name]
            if len(payload) != entry["size_bytes"]:
                issues.append(f"size mismatch: {name}")
            if sha256_bytes(payload) != entry["sha256"]:
                issues.append(f"SHA256 mismatch: {name}")

        if embedded.get("file_count") != len(recorded):
            issues.append(
                f"file_count {embedded.get('file_count')} != recorded entries {len(recorded)}"
            )

        for record in embedded.get("frozen_evidence_verified", []):
            if not record.get("sha256_matches_manifest"):
                issues.append(f"frozen evidence SHA256 differs from release manifest: {record['path']}")
            if not record.get("size_matches_manifest"):
                issues.append(f"frozen evidence size differs from release manifest: {record['path']}")

    if release_manifest_path is not None:
        if not release_manifest_path.is_file():
            issues.append(f"missing release manifest: {release_manifest_path}")
        else:
            release_manifest = json.loads(release_manifest_path.read_text(encoding="utf-8"))
            expected_hash = sha256_bytes(release_manifest_path.read_bytes())
            if embedded.get("source_manifest_sha256") != expected_hash:
                issues.append(
                    "archive was built from a different release-manifest snapshot "
                    f"({embedded.get('source_manifest_sha256')} != {expected_hash})"
                )
            if embedded.get("release_id") != release_manifest.get("release_id"):
                issues.append(
                    f"release_id mismatch: archive {embedded.get('release_id')!r} vs "
                    f"manifest {release_manifest.get('release_id')!r}"
                )
            if embedded.get("missing_frozen_evidence"):
                issues.append(
                    "archive is missing frozen evidence: "
                    + ", ".join(embedded["missing_frozen_evidence"])
                )

    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=DEFAULT_RELEASE_MANIFEST)
    parser.add_argument(
        "--no-release-manifest",
        action="store_true",
        help="verify archive internal consistency only, without the release freeze",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    archive_path = args.archive if args.archive.is_absolute() else root / args.archive
    release_manifest_path = None
    if not args.no_release_manifest:
        release_manifest_path = (
            args.manifest if args.manifest.is_absolute() else root / args.manifest
        )

    issues = verify_archive(archive_path, release_manifest_path)
    print(json.dumps({"valid": not issues, "issues": issues}, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
