from __future__ import annotations

import argparse
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "external" / "governed-km-sabio-rk"
BASE = "https://sabiork.h-its.org"
QUERY = "ParameterType:Km AND UniProtKB_AC:*"
LEGACY_ENDPOINTS = (
    "/sabioRestWebServices/kineticlawsExportTsv?q=EntryID:1",
    "/sabioRestWebServices/searchKineticLaws/kinlaws?q=EntryID:1",
    "/sabioRestWebServices/kineticlaws/1",
)
PROHIBITED_OUTPUTS = (
    "normalized-kinetics.jsonl",
    "source-manifest.json",
    "mmseqs-family-split.json",
    "log10-km-smoke.pt",
)


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ANN001
        return None


def request(url: str, *, allow_redirects: bool = True) -> tuple[int, dict[str, str], bytes]:
    """Make one request. In particular, HTTP 404 is returned and never retried."""
    opener = urllib.request.build_opener() if allow_redirects else urllib.request.build_opener(
        NoRedirect
    )
    try:
        with opener.open(url, timeout=30) as response:
            return response.status, dict(response.headers.items()), response.read()
    except urllib.error.HTTPError as exc:
        return exc.code, dict(exc.headers.items()), exc.read()


def json_request(path: str, params: dict[str, str] | None = None) -> tuple[int, Any]:
    url = BASE + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    status, _, body = request(url)
    return status, json.loads(body) if body else None


def acquire_evidence() -> dict[str, Any]:
    attempts = []
    for path in LEGACY_ENDPOINTS:
        status, headers, _ = request(BASE + path, allow_redirects=False)
        attempts.append({
            "method": "GET",
            "url": BASE + path,
            "status": status,
            "location": headers.get("Location"),
            "retry_count": 0,
            "disposition": "recorded_and_not_retried" if headers.get("Location") == "/ui/404"
            else "recorded",
        })

    count_status, count = json_request("/api/ft/proxy-select", {
        "q": QUERY, "rows": "0", "wt": "json",
    })
    native_status, native_count = json_request("/api/ft/proxy-select", {
        "q": QUERY + " AND IsRecombinant:false", "rows": "0", "wt": "json",
    })
    sample_status, sample = json_request("/api/sabio/kinlaw-entry/39181/entity-website")
    sequence_status, sequence = json_request("/api/sabio/uniprot/Q9UNW1/entity-website")
    asset_status, _, asset = request(BASE + "/assets/index-CK1ttORE.js")
    asset_text = asset.decode("utf-8")
    terms_phrases = (
        "non-exclusive and non-transferable license",
        "Non-Commercial Purpose only",
        "solely for internal non-commercial research and academic purposes",
        "Users will cite SABIO-RK in publications or presentations",
    )
    return {
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "attempted_endpoints": attempts,
        "current_api": {
            "query_url": BASE + "/api/ft/proxy-select",
            "query": QUERY,
            "status": count_status,
            "matching_entries": count.get("response", {}).get("numFound"),
            "nonrecombinant_query_status": native_status,
            "nonrecombinant_matching_entries": native_count.get("response", {}).get("numFound"),
            "sample_entry_url": BASE + "/api/sabio/kinlaw-entry/39181/entity-website",
            "sample_entry_status": sample_status,
            "sample_entry": sample.get("data") if sample else None,
            "uniprot_entity_status": sequence_status,
            "uniprot_entity": sequence.get("data") if sequence else None,
        },
        "terms_bundle": {
            "url": BASE + "/assets/index-CK1ttORE.js",
            "status": asset_status,
            "size_bytes": len(asset),
            "sha256": hashlib.sha256(asset).hexdigest(),
            "required_phrases_present": {
                phrase: phrase in asset_text for phrase in terms_phrases
            },
        },
    }


def build_blocker(evidence: dict[str, Any]) -> dict[str, Any]:
    api = evidence["current_api"]
    sample = api.get("sample_entry") or {}
    km = [
        parameter for parameter in sample.get("kineticlaw", {}).get("parameter", [])
        if parameter.get("parameter_type", {}).get("name") == "Km"
    ]
    return {
        "schema_version": 1,
        "source_id": "sabio-rk-public-api-2026-07-27",
        "task": "governed row-level log10_km training",
        "status": "blocked_skip_training",
        "decision": {
            "normalized_rows_emitted": 0,
            "mmseqs_run": False,
            "split_emitted": False,
            "training_performed": False,
            "checkpoint_emitted": False,
            "reason_codes": [
                "PERMISSION_SCOPE_NOT_GOVERNED",
                "EXACT_ASSAYED_SEQUENCE_NOT_EXPOSED",
            ],
        },
        "scope": "Targeted SABIO-RK public API and SABIO-RK terms only; no broad search.",
        "retrieved_at": evidence["retrieved_at"],
        "attempted_endpoints": evidence["attempted_endpoints"],
        "positive_api_evidence": {
            "km_uniprot_entries_reported": api.get("matching_entries"),
            "km_uniprot_nonrecombinant_entries_reported": api.get(
                "nonrecombinant_matching_entries"
            ),
            "stable_sample_entry_id": sample.get("id"),
            "stable_sample_reaction_id": sample.get("reaction", {}).get("id"),
            "sample_original_km": ({
                "value": km[0].get("start_value"),
                "unit": km[0].get("unit", {}).get("name"),
                "normalized_value_molar": km[0].get("n_start_value"),
                "species_key": km[0].get("species", {}).get("species_key"),
            } if km else None),
            "sample_pubmed_id": sample.get("publication", {}).get("pubmed_id"),
            "sample_assay_conditions": sample.get("experimental_conditions"),
            "sample_uniprot_accessions": [
                protein.get("uniprot_id")
                for protein in sample.get("enzyme_description", {}).get("proteins", [])
            ],
            "terms_bundle": evidence["terms_bundle"],
        },
        "permission_evidence": {
            "terms_url": BASE + "/ui/terms",
            "observed_terms": "Non-Commercial Purpose License",
            "allowed_scope": "internal non-commercial research and academic purposes only",
            "transferability": "non-transferable",
            "commercial_use": "requires a separately negotiated HITS license",
            "attribution": "SABIO-RK citation required when extracted data are used",
            "exact_terms_phrases": evidence["terms_bundle"]["required_phrases_present"],
            "conclusion": "No project/user purpose declaration or permission was available to bind this restricted license to training and redistribution of row-level derived artifacts.",
        },
        "blocking_gates": [
            {
                "gate": "permission_for_intended_training_and_artifact_use",
                "reason_code": "PERMISSION_SCOPE_NOT_GOVERNED",
                "detail": "The public grant is purpose-restricted and non-transferable, not a generally reusable open-data license. The workspace does not declare an eligible licensee, exclusively internal non-commercial use, or authorization to redistribute derived row-level training artifacts.",
                "resolution": "Obtain written HITS permission covering the intended project, model training, derived artifacts/checkpoints, and planned sharing, or bind an authorized non-commercial-only governance declaration approved by the data controller.",
            },
            {
                "gate": "exact_assayed_sequence_or_unambiguous_mapping",
                "reason_code": "EXACT_ASSAYED_SEQUENCE_NOT_EXPOSED",
                "detail": "The API exposes UniProt accessions and construct descriptors, but its UniProt entity endpoint exposes metadata rather than sequence. Canonical UniProt sequence retrieval would not resolve tags, truncations, mutations, processing, or other assayed-construct differences.",
                "resolution": "Acquire exact final assayed chains per entry or a source-backed accession/isoform mapping that explicitly proves the canonical chain was assayed unchanged.",
            },
        ],
        "claim_boundary": "No UniKP Km pickle, prediction, inferred unit, normalized label, family split, MMseqs result, or checkpoint was used or produced.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in PROHIBITED_OUTPUTS:
        (args.output_dir / name).unlink(missing_ok=True)
    blocker = build_blocker(acquire_evidence())
    (args.output_dir / "blocker.json").write_text(
        json.dumps(blocker, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
