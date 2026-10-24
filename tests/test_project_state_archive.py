"""Pins the 2026-07-27 project-state archive: every receipted file still matches
its recorded SHA256 (live pointers excepted), RESUME_HERE.md names a current
release manifest, and the archived claim boundary stays honestly blocked.

中文：钉住 2026-07-27 项目状态归档：回执文件仍与记录的 SHA256 一致、RESUME_HERE 指向当前 release 清单、主张边界保持诚实阻塞。
"""

import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "artifacts/project-state-archive-2026-07-27.sha256.json"

# RESUME_HERE.md is a live pointer: later releases intentionally update it to name the
# current manifest. Its 2026-07-27 identity remains in the receipt and in git history,
# so the live file is checked for a current release pointer rather than a frozen hash.
LIVE_POINTERS = {"RESUME_HERE.md"}

pytestmark = pytest.mark.skipif(
    not RECEIPT.is_file(),
    reason="local archive receipt under gitignored artifacts/ is required",
)


def test_project_state_archive_files_match_receipt():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for entry in receipt["files"]:
        if entry["path"] in LIVE_POINTERS:
            continue
        payload = (ROOT / entry["path"]).read_bytes()
        assert len(payload) == entry["size_bytes"]
        assert hashlib.sha256(payload).hexdigest() == entry["sha256"]


def test_live_resume_pointer_names_a_release_manifest():
    resume = (ROOT / "RESUME_HERE.md").read_text(encoding="utf-8")
    assert "software_implementation_release_v" in resume


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
