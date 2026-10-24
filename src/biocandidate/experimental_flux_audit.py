"""Fail-closed audit of the external anaerobic yeast flux table (yeast-GEM v9.1.0).

Fetches the two pinned upstream release files, records their byte identities,
and audits the headerless flux table against the local Yeast-MetaTwin model.
The source carries no header, no units, and no condition/strain metadata, so
the audit always concludes ``status: blocked`` and enumerates the missing
provenance fields; no experimental-flux dataset, split, manifest, or checkpoint
is ever produced. The written report is the fail-closed artifact that must
exist before any experimental-flux work can be attempted.

中文：对外部厌氧酵母通量表（yeast-GEM v9.1.0）做 fail-closed 审计。抓取两个固定
版本的上游发布文件，记录其字节身份，并对照本地 Yeast-MetaTwin 模型审计该无表头
通量表。数据源没有表头、单位，也没有条件/菌株元数据，因此审计结论恒为
``status: blocked``，并列出缺失的来源凭证字段；绝不产出任何实验通量数据集、
切分、清单或检查点。写出的报告是任何实验通量工作启动前必须存在的 fail-closed
工件。
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path


SOURCE_TAG = "v9.1.0"
REPOSITORY_URL = "https://github.com/SysBioChalmers/yeast-GEM"
API_ROOT = "https://api.github.com/repos/SysBioChalmers/yeast-GEM"
FLUX_PATH = "data/physiology/flux_data_anaerobic.tsv"
LICENSE_PATH = "LICENSE.md"
# Upstream git blob SHA-1 values of the pinned release files, recorded verbatim
# in the report next to the fetched bytes' own SHA256; a change upstream
# invalidates the recorded provenance even though this audit never re-derives
# the blob hashes itself.
EXPECTED_BLOBS = {
    FLUX_PATH: "ba9166c4a2e0319cf2001e82218d282a1d337b0e",
    LICENSE_PATH: "d7244c0d623e30a760d098e4ed829f4fc37ac262",
}
ATTEMPTED_URLS = (
    f"{API_ROOT}/contents/{FLUX_PATH}?ref={SOURCE_TAG}",
    f"{API_ROOT}/contents/{LICENSE_PATH}?ref={SOURCE_TAG}",
    f"{API_ROOT}/releases/tags/{SOURCE_TAG}",
    "https://www.ebi.ac.uk/metabolights/ws/studies?query=%2213C%22%20AND%20yeast",
    "https://www.ebi.ac.uk/biostudies/api/v1/search?query=yeast%2013C%20flux&size=20",
    "https://zenodo.org/api/records?q=%22Saccharomyces%20cerevisiae%22%20%2213C-MFA%22&size=20",
    "https://fairdomhub.org/search.json?q=yeast%2013C%20flux",
    "https://www.omicsdi.org/ws/dataset/search?query=organism:%22Saccharomyces%20cerevisiae%22%20AND%20%2213C-MFA%22&size=20",
)
REQUIRED_FIELDS = (
    "stable_row_id", "experimental_evidence_type", "strain_and_genotype",
    "biological_replicate", "medium_and_feed", "oxygenation", "temperature",
    "growth_phase", "sampling_time", "original_flux_unit", "uncertainty",
    "measurement_method", "reaction_direction", "sign_convention",
    "biomass_normalization", "source_row_citation",
)


def _request(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github.raw+json",
            "User-Agent": "BioCandidateRanker-targeted-flux-audit/1",
        },
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read(), {key.lower(): value for key, value in response.headers.items()}


def _model_reaction_ids(model_path: Path) -> set[str]:
    text = model_path.read_text(encoding="utf-8")
    return set(re.findall(r'^\s*- id: ["\']?(r_\d+)["\']?\s*$', text, re.MULTILINE))


def inspect_flux_table(payload: bytes, model_reaction_ids: set[str]) -> dict:
    """Structural audit of the headerless TSV.

    Measures row widths, dataset label counts, coverage of explicit ``r_*``
    reaction mappings against the model, and numeric columns whose semantics
    (units, normalization) the source never declares.

    中文：对无表头 TSV 做结构审计。测量行宽、数据集标签计数、显式 ``r_*`` 反应
    映射相对模型的覆盖率，以及源文件从未声明语义（单位、归一化方式）的数值列。
    """
    text = payload.decode("utf-8-sig")
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise ValueError("source flux table is empty")
    widths = Counter(len(row) for row in rows)
    dataset_counts = Counter(row[5] if len(row) > 5 and row[5] else "<missing>" for row in rows)
    mapped = [row[6] for row in rows if len(row) > 6 and re.fullmatch(r"r_\d+", row[6])]
    unique_mapped = set(mapped)
    positive_columns = [
        index for index in (2, 3, 4)
        if any(len(row) > index and _is_positive_number(row[index]) for row in rows)
    ]
    return {
        "row_count": len(rows),
        "column_width_counts": {str(key): value for key, value in sorted(widths.items())},
        "has_header": False,
        "dataset_counts": dict(sorted(dataset_counts.items())),
        "mapped_row_count": len(mapped),
        "unmapped_row_count": len(rows) - len(mapped),
        "unique_mapped_reaction_count": len(unique_mapped),
        "mapped_ids_absent_from_yeast_metatwin": sorted(unique_mapped - model_reaction_ids),
        "all_declared_mapping_ids_in_model": unique_mapped <= model_reaction_ids,
        "candidate_numeric_columns_without_declared_semantics": positive_columns,
    }


def _is_positive_number(value: str) -> bool:
    try:
        return float(value) > 0
    except ValueError:
        return False


def audit(model_path: str | Path, output_path: str | Path) -> dict:
    """Fetch the pinned files, audit the table, and write the fail-closed report.

    ``status: blocked`` is unconditional: the source cannot supply the required
    provenance fields, so the written report is the audit's only artifact and no
    dataset, split, or checkpoint is ever produced.

    中文：抓取固定版本的文件，审计通量表，并写出 fail-closed 报告。
    ``status: blocked`` 是无条件的：数据源无法提供所需的来源凭证字段，因此写出的
    报告是本审计唯一的工件，绝不产出任何数据集、切分或检查点。
    """
    model = Path(model_path).resolve()
    if not model.is_file():
        raise RuntimeError(f"Yeast-MetaTwin model asset is not a file: {model}")
    source_url = ATTEMPTED_URLS[0]
    license_url = ATTEMPTED_URLS[1]
    payload, source_headers = _request(source_url)
    license_payload, license_headers = _request(license_url)
    table = inspect_flux_table(payload, _model_reaction_ids(model))
    blockers = [
        {
            "scope": "all_rows",
            "missing_fields": list(REQUIRED_FIELDS),
            "detail": (
                "The source has no header or data dictionary. Numeric columns cannot be assigned "
                "absolute versus glucose-normalized semantics or units without guessing."
            ),
        },
        {
            "scope": "unmapped_rows",
            "missing_fields": ["yeast_metatwin_reaction_id", "mapping_confidence"],
            "affected_rows": table["unmapped_row_count"],
            "detail": "Rows without an explicit r_* identifier are excluded; name matching is prohibited.",
        },
        {
            "scope": "dataset",
            "missing_fields": ["leakage_safe_condition_split"],
            "detail": (
                "Condition, strain, and replicate identities are absent, so groups cannot be "
                "partitioned without leakage."
            ),
        },
    ]
    result = {
        "schema_version": 1,
        "status": "blocked",
        "task": "experimental_log10_flux",
        "retrieved_on": date.today().isoformat(),
        "decision": "skip_training_and_checkpoint",
        "claim_boundary": "No FBA solution was used as a label; no experimental model was trained.",
        "source": {
            "repository": REPOSITORY_URL,
            "release": SOURCE_TAG,
            "dataset_path": FLUX_PATH,
            "stable_file_url": source_url,
            "git_blob_sha1": EXPECTED_BLOBS[FLUX_PATH],
            "sha256": hashlib.sha256(payload).hexdigest(),
            "size_bytes": len(payload),
            "response_etag": source_headers.get("etag"),
            "license": "CC-BY-4.0",
            "license_path": LICENSE_PATH,
            "license_url": license_url,
            "license_git_blob_sha1": EXPECTED_BLOBS[LICENSE_PATH],
            "license_sha256": hashlib.sha256(license_payload).hexdigest(),
            "license_response_etag": license_headers.get("etag"),
            "citations_as_encoded_by_source": sorted(table["dataset_counts"]),
        },
        "model": {
            "path": str(model),
            "sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
        },
        "table_audit": table,
        "blockers": blockers,
        "attempted_source_urls": list(ATTEMPTED_URLS),
        "skipped_404_urls": [],
        "prohibited_inputs": ["FBA solutions", "model predictions", "name-guessed mappings"],
        "outputs_not_created": [
            "experimental_flux.jsonl", "split.json", "manifest.json",
            "experimental_smoke.pt", "experimental_smoke_metrics.json",
        ],
    }
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
