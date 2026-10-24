"""Verify every file identity recorded in the software release manifest.

Reads the release manifest JSON (default
configs/software_implementation_release_v11.json), cross-checks its
package_version against pyproject.toml when both are present, and re-hashes
every frozen_files entry against the given repository root. Prints the JSON
verdict {"valid": bool, "issues": [...]} and exits non-zero on any issue so it
can gate release and deposit steps; those two output keys are a test contract.

中文：逐文件校验软件发布清单记录的身份：读发布清单（默认 v10），交叉核对
package_version 与 pyproject.toml，逐文件重算大小与 SHA256；输出 JSON 的
"valid"/"issues" 键是测试契约，不得更改。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path


def verify_release(root: Path, manifest_path: Path) -> list[str]:
    """Return the release's integrity issues; an empty list means the freeze holds.
中文：返回该发布的完整性问题；空列表即冻结仍然成立。
"""
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
        default=Path("configs/software_implementation_release_v11.json"),
    )
    args = parser.parse_args()
    issues = verify_release(args.root, args.manifest)
    print(json.dumps({"valid": not issues, "issues": issues}, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
