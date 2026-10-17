"""Storage governance tests use synthetic files and never modify scientific evidence."""

import os
from pathlib import Path

import pytest

from scripts.manage_data_pool import OUTPUT, POOL, manage, storage_kind


@pytest.fixture
def pool(tmp_path, monkeypatch):
    from biocandidate import temporal_registry

    monkeypatch.setattr(temporal_registry, "build_registry", lambda *_: {
        "counts": {"accepted_records": 1},
        "sources": [{"source_id": "sample", "accepted_evidence_final": True,
                     "csv_accepted_records": 1, "disposition": "accepted_pool_source"}],
    })
    for name in ("configs/temporal_absolute_kinetics_protocol.json",
                 "configs/software_implementation_release_v7.json",
                 "artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json"):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")
    source = tmp_path / POOL / "sample"
    source.mkdir(parents=True)
    (source / "first.csv").write_bytes(b"abc")
    (source / "copy.csv").write_bytes(b"abc")
    return tmp_path


def test_incremental_snapshot_and_duplicates(pool):
    first = manage(pool, 10000)
    frozen = pool / first["evidence_snapshot"]
    content = frozen.read_bytes()
    assert first["review_queue"] == ["sample"]
    assert first["duplicate_groups"][0]["redundant_bytes"] == 3
    second = manage(pool, 10000)
    assert second["hashed_this_run"] == 0
    assert second["hashes_reused"] == first["hashed_this_run"]
    assert second["review_queue"] == []
    assert second["evidence_snapshot"] == first["evidence_snapshot"]
    path = pool / POOL / "sample/first.csv"
    path.write_bytes(b"abcd")
    third = manage(pool, 10000)
    assert third["hashed_this_run"] == 1
    assert third["review_queue"] == ["sample"]
    assert third["evidence_snapshot"] != first["evidence_snapshot"]
    assert frozen.read_bytes() == content
    assert path.read_bytes() == b"abcd"
    assert manage(pool, 10000)["review_queue"] == ["sample"]
    acknowledged = manage(pool, 10000, reviewed=("sample",))
    assert acknowledged["sources"][0]["review_required"] is False
    assert manage(pool, 10000)["review_queue"] == []


def test_budget_and_full_verification(pool):
    report = manage(pool, 0)
    assert report["pending_hashes"] > 0
    assert report["evidence_snapshot"] is None
    first = manage(pool, 10000)
    path = pool / POOL / "sample/first.csv"
    stat = path.stat()
    path.write_bytes(b"xyz")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    verified = manage(pool, 10000, verify=True)
    assert verified["evidence_snapshot"] != first["evidence_snapshot"]
    assert verified["review_queue"] == ["sample"]


def test_control_change_and_source_removal(pool):
    manage(pool, 10000)
    (pool / "configs/temporal_absolute_kinetics_protocol.json").write_text('{"v": 2}')
    report = manage(pool, 10000)
    assert report["controls_changed"]
    assert report["review_queue"] == ["sample"]
    for path in (pool / POOL / "sample").iterdir():
        path.unlink()
    report = manage(pool, 10000)
    assert report["removed_sources"] == ["sample"]
    assert report["review_queue"] == ["sample"]
    assert len(report["removed_paths"]) == 2
    assert (pool / OUTPUT / "index.sqlite").exists()


def test_storage_categories():
    assert storage_kind("artifacts/model/pip-cache/abc.body") == "cache"
    assert storage_kind("artifacts/model/venv/lib/test.so") == "environment"
    assert storage_kind(f"{POOL}/sample/cluster_tmp/file") == "cache"
    assert storage_kind(f"{POOL}/sample/raw.csv") == "evidence"
    assert storage_kind("artifacts/training/best.pt") == "research_output"
    assert Path(OUTPUT).parts[0] == "artifacts"


def test_layers_do_not_admit_unregistered_sources(pool, monkeypatch):
    from biocandidate import temporal_registry

    monkeypatch.setattr(temporal_registry, "build_registry", lambda *_: {
        "counts": {"accepted_records": 0},
        "sources": [{"source_id": "sample", "csv_accepted_records": 0,
                     "disposition": "finalized_zero_candidate_source"}],
    })
    unknown = pool / POOL / "unregistered"
    unknown.mkdir()
    (unknown / "raw.txt").write_text("raw")
    report = manage(pool, 10000)
    layers = {s["source_id"]: s["layer"] for s in report["sources"]}
    assert layers == {"sample": "archive_candidate", "unregistered": "pending"}
    assert (pool / POOL / "sample/first.csv").read_bytes() == b"abc"
