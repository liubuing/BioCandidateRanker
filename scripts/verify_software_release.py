"""Verify every file identity in the software implementation release manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path


def verify_release(root: Path, manifest_path: Path) -> list[str]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    issues: list[str] = []
    pyproject = root / "pyproject.toml"
    if pyproject.is_file():
        project = tomllib.loads(pyproject.read_text(encoding="utf-8")).get("project", {})
        package_version = project.get("version")
        declared = manifest.get("package_version")
        if package_version and declared is not None and package_version != declared:
            issues.append(
                f"package_version mismatch: manifest {declared!r} vs pyproject {package_version!r}")
    for entry in manifest["frozen_files"]:
        path = root / entry["path"]
        if not path.is_file():
            issues.append(f"missing: {entry['path']}")
            continue
        payload = path.read_bytes()
        if len(payload) != entry["size_bytes"]:
            issues.append(f"size mismatch: {entry['path']}")
        digest = hashlib.sha256(payload).hexdigest()
        if digest != entry["sha256"]:
            issues.append(f"SHA256 mismatch: {entry['path']}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("configs/software_implementation_release_v1.json"),
    )
    args = parser.parse_args()
    issues = verify_release(args.root, args.manifest)
    print(json.dumps({"valid": not issues, "issues": issues}, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
