"""Exercise the release archive builder and verifier on a synthetic tree."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(module_name: str):
    spec = importlib.util.spec_from_file_location(
        module_name, ROOT / "scripts" / f"{module_name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = _load("build_release_archive")
verifier = _load("verify_release_archive")

TRACKED = "src/example.py"
TRACKED_PAYLOAD = b"print('example')\n"
EVIDENCE = "artifacts/evidence.json"
EVIDENCE_PAYLOAD = b'{"rows": 3}\n'
RELEASE_MANIFEST = "configs/tiny_release.json"


def _digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@pytest.fixture
def synthetic_root(tmp_path: Path):
    """A miniature repo whose release manifest freezes one tracked and one ignored file."""
    (tmp_path / "src").mkdir()
    (tmp_path / TRACKED).write_bytes(TRACKED_PAYLOAD)
    (tmp_path / "artifacts").mkdir()
    (tmp_path / EVIDENCE).write_bytes(EVIDENCE_PAYLOAD)
    (tmp_path / "configs").mkdir()
    release = {
        "release_id": "tiny-release",
        "package_version": "0.0.1",
        "claim_boundary": "synthetic fixture",
        "frozen_files": [
            {"path": EVIDENCE, "sha256": _digest(EVIDENCE_PAYLOAD), "size_bytes": len(EVIDENCE_PAYLOAD)},
            {"path": TRACKED, "sha256": _digest(TRACKED_PAYLOAD), "size_bytes": len(TRACKED_PAYLOAD)},
        ],
    }
    (tmp_path / RELEASE_MANIFEST).write_text(json.dumps(release), encoding="utf-8")
    return tmp_path, release


@pytest.fixture
def patched(monkeypatch):
    """Pin the tracked-file list and build timestamp so archives are reproducible."""
    monkeypatch.setattr(builder, "list_tracked_files", lambda root: [TRACKED, RELEASE_MANIFEST])
    monkeypatch.setattr(builder, "source_date", lambda root: "2026-01-01T00:00:00+00:00")
    monkeypatch.setattr(builder, "head_commit", lambda root: "0" * 40)


def _build(root: Path, release: dict, output_dir: Path) -> Path:
    entries, missing = builder.collect_entries(root, release)
    manifest = builder.build_manifest(root, Path(RELEASE_MANIFEST), release, entries, missing)
    archive_path = output_dir / "tiny-release.zip"
    builder.write_archive(archive_path, root, entries, manifest)
    return archive_path


def test_round_trip_verifies(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    archive = _build(root, release, tmp_path / "out")
    release_path = root / RELEASE_MANIFEST
    assert verifier.verify_archive(archive, release_path) == []


def test_frozen_tracked_file_is_verified_not_skipped(synthetic_root, patched, tmp_path):
    """A Git-tracked file that drifted from the freeze must be reported, not ignored."""
    root, release = synthetic_root
    release["frozen_files"][1]["sha256"] = "0" * 64
    archive = _build(root, release, tmp_path / "out")

    with zipfile.ZipFile(archive) as handle:
        embedded = json.loads(handle.read(builder.MANIFEST_MEMBER_NAME).decode("utf-8"))
    checked = {record["path"] for record in embedded["frozen_evidence_verified"]}
    assert checked == {TRACKED, EVIDENCE}

    with zipfile.ZipFile(archive) as handle:
        records = {entry["path"]: entry for entry in json.loads(
            handle.read(builder.MANIFEST_MEMBER_NAME).decode("utf-8"))["frozen_evidence_verified"]}
    assert records[TRACKED]["role"] == builder.ROLE_TRACKED
    assert records[TRACKED]["sha256_matches_manifest"] is False
    assert records[EVIDENCE]["sha256_matches_manifest"] is True


def test_rebuild_is_byte_identical(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    first = _build(root, release, tmp_path / "one")
    second = _build(root, release, tmp_path / "two")
    assert first.read_bytes() == second.read_bytes()


def test_tampered_member_is_detected(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    archive = _build(root, release, tmp_path / "out")

    rewritten = tmp_path / "tampered.zip"
    with zipfile.ZipFile(archive) as source, zipfile.ZipFile(rewritten, "w") as target:
        for name in source.namelist():
            payload = source.read(name)
            if name == TRACKED:
                payload = b"print('tampered')\n"
            target.writestr(name, payload)

    issues = verifier.verify_archive(rewritten, root / RELEASE_MANIFEST)
    assert any("SHA256 mismatch" in issue for issue in issues)


def test_missing_frozen_evidence_is_recorded(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    (root / EVIDENCE).unlink()
    archive = _build(root, release, tmp_path / "out")

    with zipfile.ZipFile(archive) as handle:
        embedded = json.loads(handle.read(builder.MANIFEST_MEMBER_NAME).decode("utf-8"))
    assert embedded["missing_frozen_evidence"] == [EVIDENCE]

    issues = verifier.verify_archive(archive, root / RELEASE_MANIFEST)
    assert any("missing frozen evidence" in issue for issue in issues)


def test_stale_source_manifest_is_detected(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    archive = _build(root, release, tmp_path / "out")

    release["claim_boundary"] = "edited after the archive was built"
    (root / RELEASE_MANIFEST).write_text(json.dumps(release), encoding="utf-8")

    issues = verifier.verify_archive(archive, root / RELEASE_MANIFEST)
    assert any("different release-manifest snapshot" in issue for issue in issues)


def test_internal_consistency_check_without_release_manifest(synthetic_root, patched, tmp_path):
    root, release = synthetic_root
    archive = _build(root, release, tmp_path / "out")
    assert verifier.verify_archive(archive, None) == []


@pytest.mark.parametrize(
    "relative",
    [".git/HEAD", ".zcode/plans/x.md", "src/__pycache__/a.cpython-311.pyc", "src/a.pyc"],
)
def test_excluded_paths_are_not_archived(relative):
    assert builder._excluded(relative) is True


@pytest.mark.parametrize("relative", ["src/a.py", "artifacts/x.json", "docs/notes.md"])
def test_source_paths_are_archived(relative):
    assert builder._excluded(relative) is False
