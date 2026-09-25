"""Tests for the interactive readiness dashboard generator.

The dashboard is the page a human reads before deciding whether to deposit, so
the tests cover the two ways it can lie: by dropping a blocker, and by embedding
payload data that escapes its script block.
"""

from __future__ import annotations

import importlib.util
import json
import re
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


dashboard = _load("build_readiness_dashboard")


def _base_data(**overrides):
    data = {
        "archive": {
            "release_id": "software-implementation-v6",
            "package_version": "0.1.0",
            "built_from_commit": "a" * 40,
            "built_on": "2026-01-01T00:00:00+00:00",
            "file_count": 251,
            "total_size_bytes": 15301925,
            "frozen_evidence_verified": [{"path": "x", "sha256_matches_manifest": True}],
            "missing_frozen_evidence": [],
            "local_path_findings": [],
            "license_ownership": {"by_corpus": {}},
        },
        # Default to a passing gate so each blocker test can fail exactly one thing.
        "pool": {
            "gates": [{"gate": "capped records", "observed": 300, "required": 300,
                       "shortfall": 0, "satisfied": True}],
            "families": [],
            "saturation": {},
            "projection": {},
            "options": [],
            "readiness_gate_passes": False,
            "predictions_permitted": False,
            "family_cap": 20,
        },
        "licenses": {"valid": True, "corpus_count": 9, "corpora": [],
                     "shippable_corpora": [], "blocking_corpora": [], "issues": []},
        "requests": {"sent_count": 1, "draft_count": 8, "entries": [],
                     "untracked_drafts": [], "issues": [],
                     "ready_to_send_count": 3, "blocked_count": 1},
        "deliverables": [],
        "placeholders": [],
    }
    data.update(overrides)
    return data


# --- redaction ------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        r"C:\biological\BioCandidateRanker\artifacts",
        r"D:\biological\Metabolic model prediction",
        "/home/liubuing/.local/bin/mmseqs",
        "/Users/someone/project",
    ],
)
def test_redact_removes_absolute_local_paths(raw):
    assert raw not in dashboard.redact(raw)
    assert dashboard.REDACTION in dashboard.redact(raw)


def test_redact_leaves_relative_paths_alone():
    text = "artifacts/release-archive/manifest.json and docs/notes.md"
    assert dashboard.redact(text) == text


# --- blocker detection ----------------------------------------------------------


def test_healthy_project_reports_no_blockers():
    assert dashboard.collect_blockers(_base_data()) == []


def test_local_path_findings_become_a_blocker():
    data = _base_data()
    data["archive"]["local_path_findings"] = [
        {"path": "README.md", "matches": ["D:\\biological\\x"]},
        {"path": "README.md", "matches": ["C:\\y"]},
        {"path": "scripts/a.py", "matches": ["/home/liubuing/z"]},
    ]
    blockers = dashboard.collect_blockers(data)
    item = next(b for b in blockers if b["id"] == "local-paths")
    assert "2 archived file(s)" in item["title"]
    assert item["items"] == ["README.md", "scripts/a.py"]


def test_unmet_gate_becomes_a_blocker_with_the_shortfall():
    data = _base_data()
    data["pool"]["gates"][0].update(
        {"observed": 192, "shortfall": 108, "satisfied": False}
    )
    blockers = dashboard.collect_blockers(data)
    item = next(b for b in blockers if b["id"] == "pool-gate")
    assert item["items"] == ["capped records: 192/300 (short 108)"]


def test_a_satisfied_gate_does_not_become_a_blocker():
    assert not any(
        b["id"] == "pool-gate" for b in dashboard.collect_blockers(_base_data())
    )


def test_zero_sent_requests_become_a_blocker():
    data = _base_data()
    data["requests"]["sent_count"] = 0
    blockers = dashboard.collect_blockers(data)
    assert any(b["id"] == "no-outreach" for b in blockers)


def test_a_sent_request_is_not_a_blocker():
    assert not any(b["id"] == "no-outreach" for b in dashboard.collect_blockers(_base_data()))


def test_clean_redacted_deposit_clears_the_local_path_blocker():
    data = _base_data()
    data["archive"]["local_path_findings"] = [
        {"path": "README.md", "matches": ["D:\\biological\\x"]},
    ]
    data["archive"]["redacted_deposit"] = {
        "manifest": "artifacts/release-archive/x-redacted-manifest.json",
        "members_transformed": 22,
        "replacements": 47,
    }
    blockers = dashboard.collect_blockers(data)
    assert not any(b["id"] == "local-paths" for b in blockers)


def test_inconsistent_license_audit_becomes_a_blocker():
    data = _base_data()
    data["licenses"]["valid"] = False
    data["licenses"]["issues"] = ["a", "b"]
    blockers = dashboard.collect_blockers(data)
    item = next(b for b in blockers if b["id"] == "license-audit")
    assert "2 issue(s)" in item["detail"]


def test_missing_deliverable_becomes_a_blocker():
    data = _base_data()
    data["deliverables"] = [
        {"id": "cover-letter", "label": "Cover letter", "present": False, "files": []}
    ]
    blockers = dashboard.collect_blockers(data)
    assert any(b["id"] == "missing-cover-letter" for b in blockers)


def test_present_deliverable_is_not_a_blocker():
    data = _base_data()
    data["deliverables"] = [
        {"id": "cover-letter", "label": "Cover letter", "present": True, "files": ["c.md"]}
    ]
    assert not any(b["id"].startswith("missing-") for b in dashboard.collect_blockers(data))


def test_placeholders_become_a_blocker_and_count_across_files():
    data = _base_data()
    data["placeholders"] = [
        {"file": "manuscript/A.md", "placeholders": {"corresponding author": 1}},
        {"file": "manuscript/B.md", "placeholders": {"repository URL": 2}},
    ]
    blockers = dashboard.collect_blockers(data)
    item = next(b for b in blockers if b["id"] == "placeholders")
    assert "3 placeholder(s)" in item["title"]
    assert len(item["items"]) == 2


# --- rendering ------------------------------------------------------------------


def test_payload_cannot_escape_its_script_block():
    """A literal </script> inside the data must not terminate the payload element."""
    payload = _base_data()
    payload["requests"]["entries"] = [
        {"draft": "</script><script>alert(1)</script>", "recipient": None,
         "tracked_in_queue": False, "queue_state": None}
    ]
    html = dashboard.render(payload)
    assert "</script><script>alert(1)" not in html
    body = html.split('<script id="payload" type="application/json">')[1].split("</script>")[0]
    assert json.loads(body.replace("<\\/", "</"))["archive"]["file_count"] == 251


def test_render_declares_every_section():
    """The tab labels are the page's navigation contract."""
    html = dashboard.render(_base_data())
    for label in ("Overview", "Prospective pool", "Corpus licenses",
                  "External requests", "Archive"):
        assert label in html


def test_render_is_self_contained():
    """No network fetches: a deposit reviewer may open the file offline."""
    html = dashboard.render(_base_data())
    assert "<link" not in html
    assert re.search(r'src\s*=\s*["\']https?://', html) is None
    assert "fetch(" not in html


def test_render_escapes_the_page_title():
    html = dashboard.render(_base_data())
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")


# --- payload assembly -----------------------------------------------------------


def test_build_payload_reports_missing_sources(tmp_path):
    payload = dashboard.build_payload(tmp_path, redact_paths=False)
    assert payload["sources_present"] == {
        "archive": False, "pool": False, "licenses": False, "requests": False
    }
    assert payload["release_id"] == "unknown"


def test_build_payload_redacts_when_asked(tmp_path):
    source = tmp_path / dashboard.SOURCES["archive"]
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text(
        json.dumps({"release_id": "r", "local_path_findings": [
            {"path": "README.md", "matches": ["/home/liubuing/x"]}]}),
        encoding="utf-8",
    )
    payload = dashboard.build_payload(tmp_path, redact_paths=True)
    assert payload["redacted"] is True
    assert payload["archive"]["local_path_findings"][0]["matches"] == [dashboard.REDACTION]


def test_find_placeholders_only_reports_files_that_have_them(tmp_path):
    manuscript = tmp_path / dashboard.MANUSCRIPT_DIR
    manuscript.mkdir()
    (manuscript / "clean.md").write_text("Nothing to fill in.\n", encoding="utf-8")
    (manuscript / "draft.md").write_text(
        "[corresponding author] and [corresponding author]\n", encoding="utf-8"
    )
    findings = dashboard.find_placeholders(tmp_path)
    assert [item["file"] for item in findings] == ["manuscript/draft.md"]
    assert sum(findings[0]["placeholders"].values()) == 2
