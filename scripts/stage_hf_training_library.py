"""Stage the Hugging Face training-library dataset package.

Collects only artifacts whose redistribution is permitted by
configs/data_license_inventory.json, redacts absolute local paths with the
same rules as build_release_archive.py, converts the frozen homology-cold
split to Parquet, and verifies no staged text file still embeds a local path.

The UniKP/DLKcat measurement rows themselves are deliberately not staged:
their corpus is marked `not_verified` and blocks public release.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_release_archive import LOCAL_PATH_PATTERNS, redact_payload  # noqa: E402

BUILD = ROOT / "artifacts" / "hf-upload-build"
STAGE = BUILD / "dataset"

# repo-relative source -> staged destination (relative to STAGE)
TEXT_FILES = {
    "artifacts/homology-final/homology_split.json": "data/homology_split.json",
    "artifacts/unikp_audit.json": "data/unikp_audit.json",
    "artifacts/unikp-km-training-blocker-audit.json": "data/unikp-km-training-blocker-audit.json",
    "artifacts/absolute-activity-training-blocker-audit.json": "data/absolute-activity-training-blocker-audit.json",
    "artifacts/full-homology-baseline/test_metrics.json": "data/full-homology-baseline/test_metrics.json",
    "artifacts/external/enzengdb-v1/audit.json": "data/enzengdb-v1/audit.json",
    "artifacts/external/enzengdb-v1/audit-homology-cold.json": "data/enzengdb-v1/audit-homology-cold.json",
    "artifacts/external/enzengdb-v1/selection.json": "data/enzengdb-v1/selection.json",
    "artifacts/external/enzengdb-v1/selection-homology-cold.json": "data/enzengdb-v1/selection-homology-cold.json",
    "artifacts/external/enzengdb-v1/rank-split/rank_split.json": "data/enzengdb-v1/rank-split/rank_split.json",
    "artifacts/external/enzengdb-v1/rank-split/campaign_representatives.fasta": "data/enzengdb-v1/rank-split/campaign_representatives.fasta",
    "artifacts/external/enzengdb-v1/ranker-three-seed-summary.json": "data/enzengdb-v1/ranker-three-seed-summary.json",
    "artifacts/external/enzengdb-v1/zero-shot-kcat-proxy.json": "data/enzengdb-v1/zero-shot-kcat-proxy.json",
    "artifacts/external/absolute-kinetics-screen/dryad-4964723/selection.json": "data/lunzer-imdh/selection.json",
    "artifacts/external/absolute-kinetics-screen/dryad-4964723/homology-audit.json": "data/lunzer-imdh/homology-audit.json",
    "artifacts/external/absolute-kinetics-screen/dryad-4964723/biocandidate-full-evaluation.json": "data/lunzer-imdh/biocandidate-full-evaluation.json",
    "artifacts/external/absolute-kinetics-screen/dryad-4964723/dlkcat-evaluation.json": "data/lunzer-imdh/dlkcat-evaluation.json",
    "artifacts/external/temporal-absolute-kinetics/registry-readiness-audit.json": "data/temporal/registry-readiness-audit.json",
    "artifacts/external/temporal-global-family-audit/readiness-audit.json": "data/temporal/readiness-audit.json",
    "artifacts/external/temporal-global-family-audit/global-accepted-components.fasta": "data/temporal/global-accepted-components.fasta",
    "artifacts/external/temporal-global-family-audit/global-accepted-mapping.csv": "data/temporal/global-accepted-mapping.csv",
    "artifacts/uncertainty-calibration-v2/grouped-cv.json": "data/uncertainty-calibration-v2/grouped-cv.json",
    "configs/unikp_source_manifest.json": "configs/unikp_source_manifest.json",
    "configs/data_license_inventory.json": "configs/data_license_inventory.json",
    "configs/temporal_absolute_kinetics_protocol.json": "configs/temporal_absolute_kinetics_protocol.json",
}

RELEASE_FILES = {
    "artifacts/hf-upload-build/software-implementation-v10-redacted.zip": "release/software-implementation-v10-redacted.zip",
    "artifacts/hf-upload-build/software-implementation-v10-redacted-manifest.json": "release/software-implementation-v10-redacted-manifest.json",
}


def main() -> int:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)

    copied = []
    for src_rel, dst_rel in {**TEXT_FILES, **RELEASE_FILES}.items():
        src = ROOT / src_rel
        dst = STAGE / dst_rel
        if not src.is_file():
            print(f"MISSING: {src_rel}")
            return 1
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        copied.append(dst_rel)

    # Redact local paths in every staged text member.
    redacted = {}
    for rel in copied:
        path = STAGE / rel
        if path.suffix in {".zip", ".pt"}:
            continue
        payload = path.read_bytes()
        transformed, count = redact_payload(payload)
        if count:
            path.write_bytes(transformed)
            redacted[rel] = count

    # Parquet view of the frozen homology-cold split.
    split = json.loads((STAGE / "data/homology_split.json").read_text(encoding="utf-8"))
    import pandas as pd

    frame = pd.DataFrame(split["rows"])[["source_row", "sequence_id", "cluster", "split"]]
    frame = frame.sort_values("source_row").reset_index(drop=True)
    frame.to_parquet(STAGE / "data/homology_split.parquet", index=False)

    # Residual local-path guard over every staged text file.
    residual = []
    for path in STAGE.rglob("*"):
        if not path.is_file() or path.suffix in {".zip", ".pt"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if any(pattern.search(text) for pattern in LOCAL_PATH_PATTERNS):
            residual.append(str(path.relative_to(STAGE)))

    summary = {
        "staged_root": str(STAGE),
        "file_count": len(copied) + 1,
        "split_rows": len(frame),
        "split_counts": frame["split"].value_counts().to_dict(),
        "redacted_files": redacted,
        "residual_local_paths": residual,
    }
    print(json.dumps(summary, indent=2))
    return 1 if residual else 0


if __name__ == "__main__":
    raise SystemExit(main())
