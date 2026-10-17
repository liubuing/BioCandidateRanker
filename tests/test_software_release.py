import json
from pathlib import Path

from scripts.verify_software_release import verify_release


ROOT = Path(__file__).resolve().parents[1]


# The closed baseline. Bump both when a new release is frozen; the tests below
# then verify the new one, and the superseded release joins the historical list.
CURRENT_RELEASE_ID = "software-implementation-v7"
CURRENT_RELEASE_MANIFEST = "configs/software_implementation_release_v7.json"


def test_current_release_manifest_is_internally_consistent():
    manifest_path = ROOT / CURRENT_RELEASE_MANIFEST
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["release_id"] == CURRENT_RELEASE_ID
    assert manifest["temporal_pool"]["ready"] is False
    assert manifest["temporal_pool"]["predictions_permitted"] is False
    assert verify_release(ROOT, manifest_path) == []


def test_current_release_does_not_widen_the_frozen_scope():
    """A re-freeze re-pins content; it must not change which files are frozen."""
    manifest = json.loads((ROOT / CURRENT_RELEASE_MANIFEST).read_text(encoding="utf-8"))
    # release id "software-implementation-v6" -> "software_implementation_release_v6.json"
    predecessor_slug = manifest["predecessor"].replace("-", "_")
    predecessor_file = predecessor_slug.replace("implementation_", "implementation_release_")
    predecessor = json.loads(
        (ROOT / "configs" / f"{predecessor_file}.json").read_text(encoding="utf-8")
    )
    assert [entry["path"] for entry in manifest["frozen_files"]] == [
        entry["path"] for entry in predecessor["frozen_files"]
    ]


def test_v1_through_v6_release_manifests_exist_as_historical_records():
    for path, release_id, predecessor in (
        ("configs/software_implementation_release_v1.json", "software-implementation-v1", None),
        ("configs/software_implementation_release_v2.json", "software-implementation-v2", "software-implementation-v1"),
        ("configs/software_implementation_release_v3.json", "software-implementation-v3", "software-implementation-v2"),
        ("configs/software_implementation_release_v4.json", "software-implementation-v4", "software-implementation-v3"),
        ("configs/software_implementation_release_v5.json", "software-implementation-v5", "software-implementation-v4"),
        ("configs/software_implementation_release_v6.json", "software-implementation-v6", "software-implementation-v5"),
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
