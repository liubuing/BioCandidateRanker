import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "acquire_sabio_rk_km.py"
SPEC = importlib.util.spec_from_file_location("acquire_sabio_rk_km", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_frozen_sabio_blocker_fails_closed() -> None:
    source = ROOT / "artifacts" / "external" / "governed-km-sabio-rk"
    blocker = json.loads((source / "blocker.json").read_text(encoding="ascii"))
    assert blocker["status"] == "blocked_skip_training"
    assert blocker["decision"] == {
        "checkpoint_emitted": False,
        "mmseqs_run": False,
        "normalized_rows_emitted": 0,
        "reason_codes": [
            "PERMISSION_SCOPE_NOT_GOVERNED",
            "EXACT_ASSAYED_SEQUENCE_NOT_EXPOSED",
        ],
        "split_emitted": False,
        "training_performed": False,
    }
    assert blocker["positive_api_evidence"]["stable_sample_entry_id"] == 39181
    assert blocker["positive_api_evidence"]["sample_original_km"] == {
        "normalized_value_molar": 1.6e-06,
        "species_key": "1 | 1D-myo-Inositol 1,3,4,5-tetrakisphosphate | Substrate",
        "unit": "\u00b5M",
        "value": 1.6,
    }
    assert all(blocker["permission_evidence"]["exact_terms_phrases"].values())
    for name in MODULE.PROHIBITED_OUTPUTS:
        assert not (source / name).exists()


def test_legacy_404_redirects_were_not_retried() -> None:
    blocker = json.loads((
        ROOT / "artifacts" / "external" / "governed-km-sabio-rk" / "blocker.json"
    ).read_text(encoding="ascii"))
    attempts = blocker["attempted_endpoints"]
    assert len(attempts) == 3
    assert all(item["location"] == "/ui/404" for item in attempts)
    assert all(item["retry_count"] == 0 for item in attempts)
    assert all(item["disposition"] == "recorded_and_not_retried" for item in attempts)
