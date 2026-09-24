"""Reconcile the external-request queue, the response tracker, and the drafts.

The project has three records of its outreach and they drifted apart:

  * artifacts/external/external-request-queue.json  - structured state
  * artifacts/external/external-response-tracker.csv - flat state
  * docs/external-requests/*.md                      - the drafts themselves

This script reports the union and flags any draft that no queue entry tracks,
any queue entry whose draft is missing, and any disagreement about send status.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

DEFAULT_QUEUE = Path("artifacts/external/external-request-queue.json")
DEFAULT_TRACKER = Path("artifacts/external/external-response-tracker.csv")
DEFAULT_DRAFTS = Path("docs/external-requests")

SENDER_FIELDS = ("sender_name", "sender_organization", "sender_contact")
UNFILLED_MARKERS = ("[", "TO BE IDENTIFIED")

# The queue and the CSV use different vocabularies for the same state. Only a
# genuine contradiction should be reported, so equivalent labels are grouped.
COMPATIBLE_STATES = (
    {"ready_to_send_after_sender_fields", "not_sent", "ready_to_send"},
    {"blocked_user_must_supply_recipient", "recipient_required", "blocked"},
    {"sent", "sent"},
    {"responded", "replied", "response_received", "responded"},
    {"declined", "declined", "not_provided"},
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tracker(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["request_id"]: row for row in csv.DictReader(handle)}


def draft_recipient(path: Path) -> str | None:
    """First recipient-looking line, which is how every draft in this repo opens."""
    for line in path.read_text(encoding="utf-8").splitlines()[:6]:
        stripped = line.strip()
        for prefix in ("To:", "Recipient:"):
            if stripped.startswith(prefix):
                return stripped[len(prefix):].strip()
    return None


def placeholder_recipient(value: str | None) -> bool:
    if not value:
        return True
    return any(marker in value for marker in UNFILLED_MARKERS)


def states_agree(left: str | None, right: str | None) -> bool:
    """True when two labels sit in the same compatible-state group."""
    if left == right:
        return True
    return any(left in group and right in group for group in COMPATIBLE_STATES)


def reconcile(root: Path) -> dict[str, Any]:
    issues: list[str] = []
    queue_path = root / DEFAULT_QUEUE
    tracker_path = root / DEFAULT_TRACKER
    drafts_dir = root / DEFAULT_DRAFTS

    queue = load_json(queue_path) if queue_path.is_file() else {"requests": [], "messages_sent": 0}
    if not queue_path.is_file():
        issues.append(f"missing queue: {DEFAULT_QUEUE.as_posix()}")

    tracker = load_tracker(tracker_path) if tracker_path.is_file() else {}
    if not tracker_path.is_file():
        issues.append(f"missing tracker: {DEFAULT_TRACKER.as_posix()}")

    draft_paths = sorted(drafts_dir.glob("*.md")) if drafts_dir.is_dir() else []
    if not draft_paths:
        issues.append(f"no drafts under {DEFAULT_DRAFTS.as_posix()}")

    queued_drafts = {
        Path(request["draft"]).resolve(): request
        for request in queue.get("requests", [])
        if request.get("draft")
    }
    queued_ids = {request["request_id"] for request in queue.get("requests", [])}

    entries: list[dict[str, Any]] = []
    untracked: list[str] = []

    for draft in draft_paths:
        resolved = draft.resolve()
        request = queued_drafts.get(resolved)
        relative = draft.relative_to(root).as_posix()
        recipient = draft_recipient(draft)
        entry = {
            "draft": relative,
            "recipient": recipient,
            "recipient_is_placeholder": placeholder_recipient(recipient),
            "request_id": request["request_id"] if request else None,
            "queue_status": request.get("status") if request else None,
            "tracked_in_queue": request is not None,
            "tracked_in_csv": bool(request and request["request_id"] in tracker),
        }
        if request:
            csv_row = tracker.get(request["request_id"], {})
            entry["csv_status"] = csv_row.get("status")
            entry["sent_at"] = request.get("sent_at") or csv_row.get("sent_at") or None
            entry["response_status"] = request.get("response_status") or csv_row.get("acceptance_status")
            if csv_row and csv_row.get("status") and request.get("status"):
                if not states_agree(csv_row["status"], request["status"]):
                    issues.append(
                        f"{request['request_id']}: queue status {request['status']!r} disagrees "
                        f"with tracker status {csv_row['status']!r}"
                    )
        else:
            untracked.append(relative)
            issues.append(
                f"{relative} has no queue entry, so its send state is unrecorded"
            )
        entries.append(entry)

    orphan_queue_entries = []
    for request in queue.get("requests", []):
        draft = Path(request["draft"]).resolve() if request.get("draft") else None
        if draft and not draft.is_file():
            orphan_queue_entries.append(request["request_id"])
            issues.append(
                f"{request['request_id']}: queued draft is missing on disk ({request['draft']})"
            )

    csv_only = sorted(set(tracker) - queued_ids)
    for request_id in csv_only:
        issues.append(f"{request_id}: present in the tracker but absent from the queue")

    sent = [entry for entry in entries if entry.get("sent_at")]
    ready = [entry for entry in entries if entry.get("queue_status") == "ready_to_send_after_sender_fields"]
    blocked = [entry for entry in entries if entry.get("queue_status") == "blocked_user_must_supply_recipient"]
    placeholders = [entry["draft"] for entry in entries if entry["recipient_is_placeholder"]]

    return {
        "schema_version": 1,
        "queue_path": DEFAULT_QUEUE.as_posix(),
        "queue_generated_on": queue.get("generated_on"),
        "messages_sent_recorded": queue.get("messages_sent", 0),
        "draft_count": len(entries),
        "queued_count": len(queued_ids),
        "tracked_count": sum(1 for entry in entries if entry["tracked_in_queue"]),
        "sent_count": len(sent),
        "ready_to_send_count": len(ready),
        "blocked_count": len(blocked),
        "placeholder_recipient_count": len(placeholders),
        "entries": entries,
        "untracked_drafts": untracked,
        "orphan_queue_entries": orphan_queue_entries,
        "issues": issues,
        "headline": (
            "No request has been sent. "
            f"{len(placeholders)} of {len(entries)} drafts still carry a placeholder recipient, "
            "so the entire external-data workstream is idle."
        ),
        "claim_boundary": (
            "Outreach bookkeeping only. Draft text is not evidence of contact, and no reply "
            "is treated as data."
        ),
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = ["# External request status", ""]
    lines.append(report["headline"])
    lines.append("")
    lines.append("| Draft | Recipient | Queue status | Tracked | Sent |")
    lines.append("|---|---|---|---|---|")
    for entry in report["entries"]:
        recipient = entry["recipient"] or "(missing)"
        if entry["recipient_is_placeholder"]:
            recipient += " *(placeholder)*"
        lines.append(
            f"| `{entry['draft']}` | {recipient} | {entry.get('queue_status') or 'unrecorded'} | "
            f"{'yes' if entry['tracked_in_queue'] else 'no'} | {entry.get('sent_at') or 'no'} |"
        )
    lines.append("")
    lines.append(f"Ready to send: {report['ready_to_send_count']}")
    lines.append(f"Blocked on a recipient: {report['blocked_count']}")
    lines.append(f"Placeholder recipients: {report['placeholder_recipient_count']}")
    lines.append("")
    if report["issues"]:
        lines.append("## Issues")
        lines.append("")
        for issue in report["issues"]:
            lines.append(f"- {issue}")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--markdown", type=Path, default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    report = reconcile(root)

    if args.markdown:
        markdown_path = args.markdown if args.markdown.is_absolute() else root / args.markdown
        markdown_path.parent.mkdir(parents=True, exist_ok=True)
        markdown_path.write_text(render_markdown(report), encoding="utf-8")
        report["markdown"] = args.markdown.as_posix()

    if args.json_out:
        json_path = args.json_out if args.json_out.is_absolute() else root / args.json_out
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
