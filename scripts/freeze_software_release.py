"""Re-freeze the software implementation release at the current working tree.

Creates the next release manifest from its predecessor by recomputing every
frozen file's size and hash, bumping the release id, and recording which files
actually changed. The frozen file list is inherited unchanged, so a re-freeze
re-pins content without quietly widening or narrowing the scope.

Reports the drift against the predecessor before writing, so the reason for a
new release is recorded rather than implied.

Usage:
    python scripts/freeze_software_release.py \
        --predecessor configs/software_implementation_release_v6.json \
        --release-id software-implementation-v7 \
        --note "Fixes the MMseqs2 subprocess decoding ..." --write
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def binary_identity(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def verify_against(manifest: dict[str, Any], root: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Return the drift entries and the paths missing from disk."""
    drift: list[dict[str, Any]] = []
    missing: list[str] = []
    for entry in manifest["frozen_files"]:
        path = root / entry["path"]
        if not path.is_file():
            missing.append(entry["path"])
            continue
        current = binary_identity(path)
        if current["sha256"] != entry["sha256"] or current["size_bytes"] != entry["size_bytes"]:
            drift.append(
                {
                    "path": entry["path"],
                    "previous_sha256": entry["sha256"],
                    "current_sha256": current["sha256"],
                    "previous_size_bytes": entry["size_bytes"],
                    "current_size_bytes": current["size_bytes"],
                }
            )
    return drift, missing


def rebuild(
    predecessor: dict[str, Any],
    release_id: str,
    note: str,
    drift: list[dict[str, Any]],
    root: Path,
) -> dict[str, Any]:
    changed = {entry["path"] for entry in drift}
    frozen_files = []
    for entry in predecessor["frozen_files"]:
        path = root / entry["path"]
        identity = binary_identity(path)
        record: dict[str, Any] = {"path": entry["path"], **identity}
        if entry["path"] in changed:
            record["changed_in_this_release"] = True
        frozen_files.append(record)

    next_manifest = {
        key: value for key, value in predecessor.items()
        if key not in {"frozen_files", "release_id", "predecessor", "scope", "claim_boundary"}
    }
    next_manifest.update(
        {
            "release_id": release_id,
            "predecessor": predecessor["release_id"],
            "frozen_on": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "scope": (
                f"Re-freeze of the working tree inherited from {predecessor['release_id']}. "
                f"Frozen file list unchanged ({len(frozen_files)} files); content re-pinned. "
                f"{len(changed)} file(s) changed: "
                + (", ".join(sorted(changed)) if changed else "none")
                + ". "
                + note.strip()
            ),
            "claim_boundary": predecessor.get(
                "claim_boundary", "Frozen identity record; no new scientific claim."
            ),
            "frozen_files": frozen_files,
        }
    )
    return next_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--predecessor", type=Path, required=True)
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--note", default="")
    parser.add_argument("--write", action="store_true", help="write the manifest; default is a dry run")
    args = parser.parse_args()

    root = args.root.resolve()
    predecessor_path = (
        args.predecessor if args.predecessor.is_absolute() else root / args.predecessor
    )
    if not predecessor_path.is_file():
        print(json.dumps({"valid": False, "issues": [f"missing predecessor: {args.predecessor}"]}, indent=2))
        return 1

    predecessor = json.loads(predecessor_path.read_text(encoding="utf-8"))
    drift, missing = verify_against(predecessor, root)

    report: dict[str, Any] = {
        "release_id": args.release_id,
        "predecessor": predecessor["release_id"],
        "frozen_file_count": len(predecessor["frozen_files"]),
        "drift": drift,
        "missing": missing,
    }

    if missing:
        report["valid"] = False
        report["remedy"] = "a frozen file is absent from disk; restore it or revise the scope deliberately"
        print(json.dumps(report, indent=2))
        return 1

    next_manifest = rebuild(predecessor, args.release_id, args.note, drift, root)
    output_path = (
        args.output if args.output is not None
        else Path("configs") / f"{args.release_id.replace('-', '_')}.json".replace(
            "software_implementation_v", "software_implementation_release_v"
        )
    )
    if not output_path.is_absolute():
        output_path = root / output_path

    report["output"] = output_path.relative_to(root).as_posix()
    report["written"] = bool(args.write)

    if args.write:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(next_manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
