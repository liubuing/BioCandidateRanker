"""Export the frozen homology-cold split as per-partition CSVs for SOTA comparison.

Reads the external UniKP corpus through the same adapter and median-aggregation
pipeline used by every reported BioCandidateRanker run, applies the frozen
homology split manifest with source-identity validation, and writes one CSV per
partition plus a SHA256 manifest under artifacts/external/sota-homology-cold/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from biocandidate.data.homology import apply_split_manifest
from biocandidate.data.unikp import aggregate_pair_measurements, read_unikp_json


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "artifacts" / "external" / "sota-homology-cold"
EXPECTED = {
    "sha256": "13643b0b36374f8d3f64d8b014882cf1b3b58946eeaae2b9dcd59e8b2c2d6719",
    "size_bytes": 12132719,
    "row_count": 17010,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument(
        "--split-manifest", type=Path,
        default=ROOT / "artifacts" / "homology-final" / "homology_split.json",
    )
    args = parser.parse_args()

    corpus = args.corpus
    payload = corpus.read_bytes()
    identity = {
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
        "row_count": len(json.loads(payload)),
    }
    if identity["sha256"] != EXPECTED["sha256"] or identity["size_bytes"] != EXPECTED["size_bytes"]:
        raise ValueError("corpus identity does not match the frozen protocol")

    read = read_unikp_json(corpus, source_rows=None)
    if len(read.records) != 16838:
        raise ValueError(f"unexpected accepted row count {len(read.records)}")
    source_identity = {
        "row_count": identity["row_count"],
        "sha256": identity["sha256"],
        "size_bytes": identity["size_bytes"],
    }
    split_records = apply_split_manifest(
        read.records, args.split_manifest, source_identity=source_identity
    )

    aggregated = aggregate_pair_measurements(split_records)
    counts = {"train": 0, "validation": 0, "test": 0}
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    data_dir = OUTPUT_ROOT / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, dict[str, object]] = {}

    for split in ("train", "validation", "test"):
        rows = [record for record in aggregated if record.split == split]
        counts[split] = len(rows)
        target = data_dir / f"{split}.csv"
        with target.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["pair_id", "sequence", "smiles", "log10_kcat", "n_measurements"]
            )
            for index, record in enumerate(rows):
                writer.writerow([
                    f"{split}-{index:05d}",
                    record.sequence,
                    record.substrate_smiles,
                    repr(record.log10_kcat),
                    record.replicate_count,
                ])
        written[split] = {
            "path": str(target.relative_to(ROOT)),
            "sha256": sha256(target),
            "size_bytes": target.stat().st_size,
            "rows": len(rows),
        }

    manifest = {
        "schema_version": 1,
        "protocol_id": "sota-homology-cold-comparison-v1",
        "generated_on": "2026-09-16",
        "corpus_identity": identity,
        "split_manifest": {
            "path": str(args.split_manifest.relative_to(ROOT)),
            "sha256": sha256(args.split_manifest),
        },
        "aggregation": "median log10_kcat per (sequence, canonical SMILES) pair",
        "expected_counts": {"train": 13157, "validation": 1640, "test": 1646},
        "actual_counts": counts,
        "partitions": written,
    }
    if counts != manifest["expected_counts"]:
        raise ValueError(f"aggregated counts {counts} differ from frozen expectations")
    (OUTPUT_ROOT / "data-manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest["actual_counts"], indent=2))


if __name__ == "__main__":
    main()
