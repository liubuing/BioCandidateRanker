import tempfile
import unittest
from pathlib import Path

from biocandidate.data.lunzer import MUTATION_POSITIONS, REFERENCE_SEQUENCE, read_lunzer_tsv
from biocandidate.cli import build_parser
from biocandidate.data import build_manifest


class LunzerAdapterTest(unittest.TestCase):
    def test_requires_complete_frozen_landscape(self):
        header = "236\t289\t290\t296\t337\t341\tlnKmNAD\tlnKmNADP\tlnNAD\tlnNADP\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "landscape.xls"
            path.write_text(header + "D\tD\tI\tA\tG\tR\t1\t2\t3\t4\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "512 genotypes"):
                read_lunzer_tsv(path)
            self.assertEqual(build_manifest(path).row_count, 2)

    def test_reference_positions_match_declared_wildtype(self):
        self.assertEqual(
            "".join(REFERENCE_SEQUENCE[position - 1] for position in MUTATION_POSITIONS),
            "DDIAGR",
        )

    def test_cli_registers_audit_and_evaluation_commands(self):
        parser = build_parser()
        audit = parser.parse_args(["audit-lunzer", "--data", "source.xls"])
        self.assertEqual(audit.func.__name__, "audit_lunzer_command")
        evaluate = parser.parse_args([
            "evaluate-lunzer", "--checkpoint", "model.pt", "--data", "source.xls",
            "--selection-manifest", "selection.json", "--output", "metrics.json",
        ])
        self.assertEqual(evaluate.func.__name__, "evaluate_lunzer_command")
        dlkcat = parser.parse_args([
            "evaluate-lunzer-dlkcat-output", "--data", "source.xls",
            "--selection-manifest", "selection.json", "--predictions", "output.tsv",
            "--dlkcat-checkpoint", "saved_model", "--output", "metrics.json",
        ])
        self.assertEqual(dlkcat.func.__name__, "evaluate_lunzer_dlkcat_output_command")


if __name__ == "__main__":
    unittest.main()
