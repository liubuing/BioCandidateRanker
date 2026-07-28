import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "artifacts/project-state-archive-2026-07-27.sha256.json"


def test_project_state_archive_files_match_receipt():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for entry in receipt["files"]:
        payload = (ROOT / entry["path"]).read_bytes()
        assert len(payload) == entry["size_bytes"]
        assert hashlib.sha256(payload).hexdigest() == entry["sha256"]


def test_project_state_archive_preserves_claim_boundary():
    archive = json.loads(
        (ROOT / "artifacts/project-state-archive-2026-07-27.json").read_text(
            encoding="utf-8"
        )
    )
    assert archive["scientific_state"] == "blocked_pending_external_data"
    assert archive["temporal_absolute_kinetics_pool"]["model_predictions_permitted"] is False
    assert archive["external_request_state"]["messages_sent"] == 0
    assert archive["external_request_state"]["prospective_campaign_collaborator_supplied"] is False
