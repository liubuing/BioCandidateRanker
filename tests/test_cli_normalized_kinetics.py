import json
import tempfile
import unittest
from pathlib import Path

from biocandidate.cli import build_parser, evaluate_command, train_command
from biocandidate.data import build_manifest


def _label(value, unit):
    return {
        "value": value,
        "unit": unit,
        "evidence": "direct",
        "provenance": {"doi": "10.1000/synthetic"},
    }


class NormalizedKineticsCliTest(unittest.TestCase):
    def test_parser_defaults_to_unikp_and_accepts_normalized_tasks(self):
        parser = build_parser()
        unikp = parser.parse_args(["train", "--data", "data.json", "--manifest", "m.json"])
        self.assertEqual(unikp.data_format, "unikp")
        self.assertIsNone(unikp.tasks)

        normalized = parser.parse_args([
            "evaluate", "--checkpoint", "model.pt", "--data", "data.jsonl",
            "--data-format", "normalized-kinetics", "--tasks", "log10_km",
            "log10_kcat_per_km", "--manifest", "m.json", "--split-manifest", "s.json",
            "--output", "result.json",
        ])
        self.assertEqual(normalized.data_format, "normalized-kinetics")
        self.assertEqual(normalized.tasks, ["log10_km", "log10_kcat_per_km"])

    def test_parser_rejects_activity_and_unknown_data_format(self):
        parser = build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args([
                "train", "--data", "data.jsonl", "--manifest", "m.json",
                "--tasks", "log10_activity",
            ])
        with self.assertRaises(SystemExit):
            parser.parse_args([
                "train", "--data", "data.jsonl", "--manifest", "m.json",
                "--data-format", "other",
            ])

    def test_train_smoke_and_sparse_partition_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "kinetics.jsonl"
            rows = []
            for index, split in enumerate(("train", "train", "validation", "validation")):
                rows.append({
                    "schema_version": 1,
                    "sequence": "ACDE" + "F" * index,
                    "substrate_smiles": "CCO",
                    "source_dataset": "synthetic",
                    "source_row": index,
                    "split": split,
                    "labels": {"log10_kcat": _label(index + 1, "s^-1")},
                })
            data.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
            manifest = root / "manifest.json"
            identity = build_manifest(data)
            manifest.write_text(json.dumps({
                "identity": {
                    "sha256": identity.sha256,
                    "size_bytes": identity.size_bytes,
                    "row_count": identity.row_count,
                }
            }), encoding="utf-8")
            split_manifest = root / "split.json"
            split_manifest.write_text(json.dumps({
                "format_version": 1,
                "rows": [
                    {"source_row": index, "split": split}
                    for index, split in enumerate(("train", "train", "validation", "validation"))
                ],
            }), encoding="utf-8")
            parser = build_parser()
            common = [
                "train", "--data", str(data), "--data-format", "normalized-kinetics",
                "--tasks", "log10_kcat", "--manifest", str(manifest), "--epochs", "1",
                "--split-manifest", str(split_manifest), "--conflict-policy", "median",
                "--batch-size", "2", "--d-model", "8", "--num-heads", "2",
                "--protein-layers", "1", "--molecule-layers", "1", "--fusion-layers", "1",
                "--chunk-size", "8", "--context-buckets", "8", "--output-dir", str(root),
            ]
            args = parser.parse_args(common)
            train_command(args)
            self.assertTrue((root / "latest.pt").is_file())

            output = root / "evaluation.json"
            evaluate_args = parser.parse_args([
                "evaluate", "--checkpoint", str(root / "latest.pt"), "--data", str(data),
                "--data-format", "normalized-kinetics", "--tasks", "log10_kcat",
                "--manifest", str(manifest), "--split-manifest", str(split_manifest),
                "--partition", "validation", "--batch-size", "2", "--output", str(output),
            ])
            evaluate_command(evaluate_args)
            evaluation = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(evaluation["model_metrics"]["log10_kcat"]["count"], 2)

            args = parser.parse_args(common)
            args.tasks = ["log10_km"]
            with self.assertRaisesRegex(ValueError, "zero train observations"):
                train_command(args)


if __name__ == "__main__":
    unittest.main()
