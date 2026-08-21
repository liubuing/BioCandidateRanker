import json
from pathlib import Path

from scripts.verify_software_release import verify_release


ROOT = Path(__file__).resolve().parents[1]


def test_software_release_manifest_is_internally_consistent():
    manifest_path = ROOT / "configs/software_implementation_release_v4.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["release_id"] == "software-implementation-v4"
    assert manifest["predecessor"] == "software-implementation-v3"
    assert manifest["temporal_pool"]["ready"] is False
    assert manifest["temporal_pool"]["predictions_permitted"] is False
    assert verify_release(ROOT, manifest_path) == []


def test_v1_v2_v3_release_manifests_exist_as_historical_records():
    for path, release_id, predecessor in (
        ("configs/software_implementation_release_v1.json", "software-implementation-v1", None),
        ("configs/software_implementation_release_v2.json", "software-implementation-v2", "software-implementation-v1"),
        ("configs/software_implementation_release_v3.json", "software-implementation-v3", "software-implementation-v2"),
    ):
        manifest = json.loads((ROOT / path).read_text(encoding="utf-8"))
        assert manifest["release_id"] == release_id
        assert "frozen_files" in manifest
        if predecessor is not None:
            assert manifest["predecessor"] == predecessor


def test_external_data_blockers_are_explicit():
    blockers = json.loads(
        (ROOT / "configs/external_data_blockers.json").read_text(encoding="utf-8")
    )
    by_id = {item["id"]: item for item in blockers["workstreams"]}
    required = {
        "absolute_kinetics_expansion",
        "governed_km_training",
        "governed_absolute_activity_training",
        "experimental_flux_training",
        "prospective_candidate_ranking",
    }
    assert required == set(by_id)
    assert all(by_id[item]["status"] == "blocked_pending_external_data" for item in required)
