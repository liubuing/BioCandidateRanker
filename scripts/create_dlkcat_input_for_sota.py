"""Create DLKcat-format inputs from the frozen homology-cold split CSVs.

Tokenization (radius-2 fingerprints, 3-gram protein words, log2 labels) is taken
verbatim from the pinned upstream preprocessing module so the input format matches
the published DLKcat pipeline exactly. Dictionaries are built over the union of all
partitions (inputs only, no labels), and one npy set is written per partition.
Governed by configs/sota_homology_cold_comparison_protocol.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import pickle
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "artifacts" / "external" / "sota-homology-cold" / "upstream-dlkcat"
RADIUS = 2
NGRAM = 3
LOG10_OF_2 = math.log10(2)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_upstream_preprocess() -> object:
    spec = importlib.util.spec_from_file_location(
        "dlkcat_upstream_preprocess", UPSTREAM / "preprocess_all.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/data")
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "artifacts/external/sota-homology-cold/dlkcat-input")
    args = parser.parse_args()

    upstream = load_upstream_preprocess()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    partitions = {}
    for split in ("train", "validation", "test"):
        path = args.data_dir / f"{split}.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            partitions[split] = list(csv.DictReader(handle))

    compounds, adjacencies, proteins, regression = [], [], [], []
    index = {"train": 0, "validation": 0, "test": 0}
    for split in ("train", "validation", "test"):
        for row in partitions[split]:
            if "." in row["smiles"]:
                raise ValueError(f"multi-component SMILES in frozen split: {row['smiles']}")
            mol = upstream.Chem.MolFromSmiles(row["smiles"])
            if mol is None:
                raise ValueError(f"RDKit failed to parse frozen split SMILES: {row['smiles']}")
            mol = upstream.Chem.AddHs(mol)
            atoms = upstream.create_atoms(mol)
            i_jbond_dict = upstream.create_ijbonddict(mol)
            fingerprints = upstream.extract_fingerprints(atoms, i_jbond_dict, RADIUS)
            compounds.append(fingerprints)
            adjacencies.append(upstream.create_adjacency(mol))
            proteins.append(upstream.split_sequence(row["sequence"], NGRAM))
            log10_kcat = float(row["log10_kcat"])
            regression.append(np.array([log10_kcat / LOG10_OF_2]))
            index[split] += 1
        print(f"{split}: {index[split]} rows encoded")

    for split in ("train", "validation", "test"):
        split_dir = args.out_dir / split
        split_dir.mkdir(parents=True, exist_ok=True)
        start = 0 if split == "train" else index["train"]
        if split == "validation":
            start = index["train"]
        elif split == "test":
            start = index["train"] + index["validation"]
        def save_ragged(name, values):
            np.save(split_dir / name, np.asarray(values, dtype=object))

        save_ragged("compounds", compounds[start:start + index[split]])
        save_ragged("adjacencies", adjacencies[start:start + index[split]])
        save_ragged("proteins", proteins[start:start + index[split]])
        np.save(split_dir / "regression",
                np.asarray([float(item[0]) for item in regression[start:start + index[split]]],
                           dtype=np.float64).reshape(-1, 1))

    for name, dictionary in (
        ("fingerprint_dict", upstream.fingerprint_dict),
        ("atom_dict", upstream.atom_dict),
        ("bond_dict", upstream.bond_dict),
        ("edge_dict", upstream.edge_dict),
        ("sequence_dict", upstream.word_dict),
    ):
        with (args.out_dir / f"{name}.pickle").open("wb") as handle:
            pickle.dump(dict(dictionary), handle)

    manifest = {
        "schema_version": 1,
        "protocol_id": "sota-homology-cold-comparison-v1",
        "radius": RADIUS,
        "ngram": NGRAM,
        "label_encoding": "log2(kcat), upstream convention",
        "rows": index,
        "vocabulary_sizes": {
            "n_fingerprint": len(upstream.fingerprint_dict),
            "n_word": len(upstream.word_dict),
            "n_atom": len(upstream.atom_dict),
            "n_bond": len(upstream.bond_dict),
            "n_edge": len(upstream.edge_dict),
        },
        "sources": {
            "data-manifest.json": sha256(ROOT / "artifacts/external/sota-homology-cold/data-manifest.json"),
            "upstream_preprocess_all.py": sha256(UPSTREAM / "preprocess_all.py"),
        },
        "outputs": {
            f"{split}/{name}": sha256(args.out_dir / split / f"{name}.npy")
            for split in index for name in ("compounds", "adjacencies", "proteins", "regression")
        },
    }
    (args.out_dir / "input-manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest["rows"], indent=2))
    print(json.dumps(manifest["vocabulary_sizes"], indent=2))


if __name__ == "__main__":
    main()
