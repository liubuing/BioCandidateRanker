import hashlib
import json

from biocandidate import experimental_flux_audit as audit_module


FIXTURE = (
    b"GlcV\tGlcV\t6.3\t6.58\t100\tJouhten2008\tr_1714\tD-glucose exchange\teq\n"
    b"G6P\tF6P\t88\t88\t88\tJouhten2008\tr_0467\tisomerase\teq\n"
    b"other_cs\tout\t1\t0.037\t3.7\tNissen04\n"
)


def test_inspection_never_guesses_missing_mappings_or_column_semantics():
    result = audit_module.inspect_flux_table(FIXTURE, {"r_1714", "r_0467"})
    assert result["has_header"] is False
    assert result["column_width_counts"] == {"6": 1, "9": 2}
    assert result["mapped_row_count"] == 2
    assert result["unmapped_row_count"] == 1
    assert result["all_declared_mapping_ids_in_model"] is True
    assert result["candidate_numeric_columns_without_declared_semantics"] == [2, 3, 4]


def test_audit_emits_only_blocker_when_governance_is_incomplete(tmp_path, monkeypatch):
    model = tmp_path / "Yeast-MetaTwin.yml"
    model.write_text('- reactions:\n    - id: "r_1714"\n    - id: "r_0467"\n', encoding="utf-8")
    license_payload = b"Creative Commons Attribution 4.0 International\n"

    def request(url):
        payload = license_payload if url.endswith("LICENSE.md?ref=v9.1.0") else FIXTURE
        return payload, {"etag": '"fixture"'}

    monkeypatch.setattr(audit_module, "_request", request)
    output = tmp_path / "audit" / "blocker.json"
    result = audit_module.audit(model, output)
    saved = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == saved["status"] == "blocked"
    assert saved["decision"] == "skip_training_and_checkpoint"
    assert "No FBA solution" in saved["claim_boundary"]
    assert saved["source"]["sha256"] == hashlib.sha256(FIXTURE).hexdigest()
    assert "original_flux_unit" in saved["blockers"][0]["missing_fields"]
    assert "uncertainty" in saved["blockers"][0]["missing_fields"]
    assert sorted(path.name for path in output.parent.iterdir()) == ["blocker.json"]
