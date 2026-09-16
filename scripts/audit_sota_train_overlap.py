"""Audit overlap between published predictors' training corpora and the frozen test split.

Mode A (as-published checkpoints) is leakage-favorable because the UniKP/DLKcat training
corpus subsumes this project's development corpus. This script quantifies that overlap
exactly: how many frozen test pairs/sequences appear verbatim in each predictor's
published training corpus. Governed by
configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True,
                        help="UniKP/DLKcat published training corpus (raw JSON)")
    parser.add_argument("--test-csv", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/data/test.csv")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/overlap-audit.json")
    args = parser.parse_args()

    rows = json.loads(args.corpus.read_text(encoding="utf-8"))
    corpus_sequences = set()
    corpus_pairs = set()
    for row in rows:
        sequence = str(row.get("Sequence", "")).strip()
        smiles = str(row.get("Smiles", "")).strip()
        if sequence:
            corpus_sequences.add(sequence)
            corpus_pairs.add((sequence, smiles))

    with args.test_csv.open(newline="", encoding="utf-8") as handle:
        test_rows = list(csv.DictReader(handle))
    test_sequences = {row["sequence"] for row in test_rows}
    test_pairs = {(row["sequence"], row["smiles"]) for row in test_rows}

    exact_sequence_hits = sum(1 for row in test_rows if row["sequence"] in corpus_sequences)
    exact_pair_hits = sum(1 for row in test_rows if (row["sequence"], row["smiles"]) in corpus_pairs)

    receipt = {
        "schema_version": 1,
        "protocol_id": "sota-homology-cold-comparison-v1",
        "generated_on": "2026-09-16",
        "claim": "As-published UniKP/DLKcat checkpoints have seen the frozen test rows in training; Mode A numbers are a leakage-favorable upper reference only.",
        "corpus": {
            "path": str(args.corpus),
            "sha256": sha256(args.corpus),
            "rows": len(rows),
            "unique_sequences": len(corpus_sequences),
            "unique_pairs": len(corpus_pairs),
        },
        "test_partition": {
            "path": str(args.test_csv),
            "sha256": sha256(args.test_csv),
            "rows": len(test_rows),
            "unique_sequences": len(test_sequences),
            "unique_pairs": len(test_pairs),
        },
        "exact_overlap": {
            "test_rows_with_sequence_in_corpus": exact_sequence_hits,
            "test_rows_with_pair_in_corpus": exact_pair_hits,
            "sequence_overlap_fraction": round(exact_sequence_hits / len(test_rows), 4),
            "pair_overlap_fraction": round(exact_pair_hits / len(test_rows), 4),
        },
        "interpretation": "Any Mode A evaluation on this test partition is in-corpus recall, not generalization. Only Mode B (retrained on the frozen train partition) supports comparison.",
    }
    args.out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt["exact_overlap"], indent=2))


if __name__ == "__main__":
    main()
