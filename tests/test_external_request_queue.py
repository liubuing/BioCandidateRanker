import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
QUEUE = ROOT / "artifacts/external/external-request-queue.json"
TRACKER = ROOT / "artifacts/external/external-response-tracker.csv"


def test_external_request_queue_does_not_claim_messages_were_sent():
    queue = json.loads(QUEUE.read_text(encoding="utf-8"))
    assert queue["messages_sent"] == 0
    assert len(queue["requests"]) == 4
    assert all(item["sent_at"] is None for item in queue["requests"])
    assert all(item["message_reference"] is None for item in queue["requests"])

    for request in queue["requests"]:
        assert (ROOT / request["draft"]).is_file()
        assert all((ROOT / attachment).exists() for attachment in request["attachments"])


def test_response_tracker_matches_request_queue():
    queue = json.loads(QUEUE.read_text(encoding="utf-8"))
    with TRACKER.open(newline="", encoding="utf-8") as handle:
        tracker = list(csv.DictReader(handle))
    assert {row["request_id"] for row in tracker} == {
        item["request_id"] for item in queue["requests"]
    }
    assert all(row["status"] in {"not_sent", "recipient_required"} for row in tracker)
