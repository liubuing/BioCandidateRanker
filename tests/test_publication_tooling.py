"""Tests for the publication governance tooling.

These tools decide what may be deposited and what may be claimed, so the tests
focus on the guards: a blocked corpus must not ship, a gate must not silently
pass, and an untracked draft must not look recorded.
"""

from __future__ import annotations

import importlib.util
import json
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


dashboard = _load("prospective_pool_dashboard")
license_audit = _load("audit_development_corpus_license")
requests_status = _load("external_request_status")
figure1 = _load("figure1_data_task_map")
archive = _load("build_release_archive")

AUDIT_FIXTURE = {
    "generated_on": "2026-07-27",
    "readiness_gate_passes": False,
    "parameters": {"maximum_records_per_family": 20},
    "counts": {
        "accepted_sources": 18,
        "records_after_global_family_cap": 192,
        "global_families": 25,
        "substrates_after_global_family_cap": 54,
    },
    "family_record_counts_before_cap": {
        "family-0001": 20,
        "family-0002": 8,
        "family-0003": 3,
    },
}

THRESHOLDS = {
    "records_after_family_cap": 192,
    "required_records": 300,
    "global_mmseqs_families": 25,
    "required_global_mmseqs_families": 30,
    "unique_substrates": 54,
    "required_unique_substrates": 50,
    "predictions_permitted": False,
}


# --- prospective pool dashboard -------------------------------------------------


def test_gates_read_live_audit_counts_under_audit_key_names():
    """The audit names its counts differently from the release manifest."""
    gates = dashboard.evaluate_gates(AUDIT_FIXTURE, THRESHOLDS)
    observed = {gate["gate"]: gate for gate in gates}
    assert observed["capped records"]["observed"] == 192
    assert observed["capped records"]["observed_source"] == "audit"
    assert observed["global MMseqs families"]["observed"] == 25
    assert observed["unique substrates"]["observed"] == 54


def test_all_three_gates_are_evaluated():
    """A missing gate would silently understate the shortfall."""
    gates = dashboard.evaluate_gates(AUDIT_FIXTURE, THRESHOLDS)
    assert len(gates) == 3


def test_shortfall_and_satisfaction():
    gates = {gate["gate"]: gate for gate in dashboard.evaluate_gates(AUDIT_FIXTURE, THRESHOLDS)}
    assert gates["capped records"]["shortfall"] == 108
    assert gates["capped records"]["satisfied"] is False
    assert gates["global MMseqs families"]["shortfall"] == 5
    assert gates["unique substrates"]["satisfied"] is True
    assert gates["unique substrates"]["shortfall"] == 0


def test_gates_fall_back_to_the_manifest_when_the_audit_lacks_a_count():
    thin = {"counts": {"global_families": 25}, "parameters": {"maximum_records_per_family": 20}}
    gates = {gate["gate"]: gate for gate in dashboard.evaluate_gates(thin, THRESHOLDS)}
    assert gates["capped records"]["observed_source"] == "release_manifest"
    assert gates["capped records"]["observed"] == 192


def test_family_headroom_counts_capped_families_as_full():
    rows = dashboard.family_table(AUDIT_FIXTURE, 20)
    at_cap = [row for row in rows if row["at_cap"]]
    assert [row["family"] for row in at_cap] == ["family-0001"]
    assert sum(row["headroom"] for row in rows) == 0 + 12 + 17


def test_projection_only_counts_unmet_gates():
    gates = dashboard.evaluate_gates(AUDIT_FIXTURE, THRESHOLDS)
    projection = dashboard.source_projection(AUDIT_FIXTURE, gates)
    assert projection["records_needed"] == 108
    assert projection["families_needed"] == 5
    assert projection["records_per_accepted_source"] == pytest.approx(10.67, abs=0.01)


# --- license audit --------------------------------------------------------------


def _inventory(**overrides):
    entry = {
        "id": "corpus-a",
        "role": "training",
        "license_name": "unknown",
        "license_basis": "unresolved upstream terms",
        "redistribution": "not_verified",
        "redistribute_in_release": False,
        "blocking_reason": "terms unresolved",
        "evidence": ["evidence.json"],
    }
    entry.update(overrides)
    return {
        "policy": {
            "allowed_redistribution_states": [
                "permitted_with_attribution",
                "permitted_unconditionally",
                "not_verified",
                "blocked_pending_permission",
                "not_applicable_reference_only",
            ]
        },
        "corpora": [entry],
    }


@pytest.fixture
def evidence_root(tmp_path):
    (tmp_path / "evidence.json").write_text("{}", encoding="utf-8")
    return tmp_path


def test_license_audit_accepts_a_consistent_entry(evidence_root):
    report = license_audit.audit(evidence_root, _inventory(), None, cross_checks=())
    assert report["valid"] is True
    assert report["blocking_corpora"] == ["corpus-a"]


def test_license_audit_flags_missing_evidence(evidence_root):
    report = license_audit.audit(
        evidence_root, _inventory(evidence=["absent.json"]), None, cross_checks=()
    )
    assert report["valid"] is False
    assert any("missing evidence path" in issue for issue in report["issues"])


def test_license_audit_flags_shipping_a_blocked_corpus(evidence_root):
    report = license_audit.audit(
        evidence_root, _inventory(redistribute_in_release=True), None, cross_checks=()
    )
    assert report["valid"] is False
    assert any("shippable in a release" in issue for issue in report["issues"])


def test_license_audit_flags_a_blocked_corpus_without_a_reason(evidence_root):
    inventory = _inventory()
    del inventory["corpora"][0]["blocking_reason"]
    report = license_audit.audit(evidence_root, inventory, None, cross_checks=())
    assert any("no blocking_reason" in issue for issue in report["issues"])


def test_license_audit_rejects_an_unknown_redistribution_state(evidence_root):
    report = license_audit.audit(
        evidence_root, _inventory(redistribution="probably_fine"), None, cross_checks=()
    )
    assert report["valid"] is False
    assert any("not in the allowed set" in issue for issue in report["issues"])


def test_license_audit_flags_duplicate_corpus_ids(evidence_root):
    inventory = _inventory()
    inventory["corpora"].append(dict(inventory["corpora"][0]))
    report = license_audit.audit(evidence_root, inventory, None, cross_checks=())
    assert any("duplicate corpus id" in issue for issue in report["issues"])


# --- archive corpus ownership ---------------------------------------------------


def test_archive_refuses_a_measurement_bearing_blocked_corpus(tmp_path):
    entries = [{"path": "artifacts/blocked/table.csv", "role": archive.ROLE_FROZEN}]
    inventory = {
        "corpora": [
            {
                "id": "blocked",
                "redistribution": "not_verified",
                "redistribute_in_release": False,
                "blocking_reason": "unresolved",
                "release_paths": ["artifacts/blocked/"],
            }
        ]
    }
    violations, _ = archive.check_license_ownership(entries, inventory)
    assert len(violations) == 1
    assert "may not ship" in violations[0]


def test_archive_allows_a_declared_measurement_free_derivative():
    entries = [{"path": "artifacts/blocked/split.json", "role": archive.ROLE_FROZEN}]
    inventory = {
        "corpora": [
            {
                "id": "blocked",
                "redistribution": "not_verified",
                "redistribute_in_release": False,
                "blocking_reason": "unresolved",
                "release_paths": [],
                "derived_artifacts_in_release": [
                    {
                        "path": "artifacts/blocked/split.json",
                        "kind": "identifier mapping",
                        "contains_measurements": False,
                    }
                ],
            }
        ]
    }
    violations, ownership = archive.check_license_ownership(entries, inventory)
    assert violations == []
    assert ownership[0]["corpus_id"] == "blocked"


def test_archive_reports_an_unattributable_artifact():
    entries = [{"path": "artifacts/mystery/data.csv", "role": archive.ROLE_FROZEN}]
    _, ownership = archive.check_license_ownership(entries, {"corpora": []})
    assert ownership[0]["corpus_id"] is None


def test_archive_treats_project_source_as_self_owned():
    entries = [{"path": "src/biocandidate/model.py", "role": archive.ROLE_TRACKED}]
    violations, ownership = archive.check_license_ownership(entries, {"corpora": []})
    assert violations == []
    assert ownership[0]["corpus_id"] == "project-source"


def test_local_path_scan_finds_a_windows_and_a_posix_path(tmp_path):
    (tmp_path / "audit.json").write_text(
        json.dumps({"source": "C:\\biological\\Secret", "tool": "/home/liubuing/bin/x"}),
        encoding="utf-8",
    )
    entries = [{"path": "audit.json", "role": archive.ROLE_FROZEN}]
    findings = archive.scan_local_paths(entries, tmp_path)
    assert len(findings) == 1
    assert len(findings[0]["matches"]) == 2


def test_local_path_scan_skips_binary_members(tmp_path):
    (tmp_path / "weights.pt").write_bytes(b"\x00C:\\biological\x00")
    entries = [{"path": "weights.pt", "role": archive.ROLE_FROZEN}]
    assert archive.scan_local_paths(entries, tmp_path) == []


# --- external request reconciliation --------------------------------------------


def test_compatible_state_labels_are_not_reported_as_a_conflict():
    assert requests_status.states_agree("ready_to_send_after_sender_fields", "not_sent")
    assert requests_status.states_agree(
        "blocked_user_must_supply_recipient", "recipient_required"
    )
    assert not requests_status.states_agree("ready_to_send_after_sender_fields", "sent")


@pytest.mark.parametrize(
    "value,expected",
    [
        ("Mark C. Hall <mchall@purdue.edu>", False),
        ("`[COLLABORATOR NAME]`", True),
        ("[USER MUST SUPPLY A REAL EXPERIMENTAL COLLABORATOR]", True),
        (None, True),
    ],
)
def test_placeholder_recipient_detection(value, expected):
    assert requests_status.placeholder_recipient(value) is expected


def test_reconcile_flags_a_draft_with_no_queue_entry(tmp_path):
    drafts = tmp_path / requests_status.DEFAULT_DRAFTS
    drafts.mkdir(parents=True)
    (drafts / "ORPHAN_REQUEST.md").write_text("# Orphan\n\nTo: Real Person\n", encoding="utf-8")

    queue = tmp_path / requests_status.DEFAULT_QUEUE
    queue.parent.mkdir(parents=True, exist_ok=True)
    queue.write_text(json.dumps({"requests": [], "messages_sent": 0}), encoding="utf-8")

    tracker = tmp_path / requests_status.DEFAULT_TRACKER
    tracker.write_text(
        "request_id,recipient,status,sent_at,acceptance_status\n", encoding="utf-8"
    )

    report = requests_status.reconcile(tmp_path)
    assert report["untracked_drafts"] == ["docs/external-requests/ORPHAN_REQUEST.md"]
    assert report["tracked_count"] == 0


# --- figure 1 -------------------------------------------------------------------


def test_figure_svg_is_well_formed_and_covers_every_corpus():
    svg = figure1.build_svg(figure1.CORPORA)
    assert svg.startswith("<svg")
    assert svg.rstrip().endswith("</svg>")
    for corpus in figure1.CORPORA:
        assert corpus["label"] in svg


def test_figure_has_one_lane_per_corpus():
    """One y position per corpus is what keeps labels from colliding."""
    lanes = [figure1._lane_y(index) for index in range(len(figure1.CORPORA))]
    assert len(set(lanes)) == len(figure1.CORPORA)


def test_figure_marks_the_unscored_pool_as_never_scored():
    dashed = [c for c in figure1.CORPORA if c["use"] == "never scored"]
    assert [c["id"] for c in dashed] == ["temporal-pool"]
    assert "stroke-dasharray" in figure1.build_svg(figure1.CORPORA)


def test_figure_legend_lists_each_corpus_with_provenance():
    legend = figure1.build_legend(figure1.CORPORA)
    for corpus in figure1.CORPORA:
        assert corpus["license"] in legend
