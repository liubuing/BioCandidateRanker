import csv
import tempfile
import unittest
from pathlib import Path

from biocandidate.data import (
    EnzEngDBObservation,
    EnzymeSubstrateRecord,
    campaign_rank_records,
    campaign_representative_records,
    read_enzengdb_directory,
    select_primary_reactant,
)
from biocandidate.evaluation import campaign_ranking_metrics


class EnzEngDBAdapterTest(unittest.TestCase):
    def test_selects_largest_concrete_carbon_reactant(self):
        selected = select_primary_reactant("O.CCO.CCCCC>>CCCCCO")
        self.assertEqual(selected, "CCCCC")

    def test_rejects_wildcard_only_reaction(self):
        with self.assertRaisesRegex(ValueError, "no concrete"):
            select_primary_reactant("[*:1]C.[*:2]N>>[*:1]CN[*:2]")

    def test_reads_valid_rows_and_audits_invalid_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "campaign.csv"
            fields = [
                "id", "amino_acid_substitutions", "aa_sequence",
                "fitness_value", "reaction_smiles",
            ]
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow({
                    "id": "parent", "amino_acid_substitutions": "#PARENT#",
                    "aa_sequence": "ACDE", "fitness_value": "1.5",
                    "reaction_smiles": "CCO.O>>CC=O",
                })
                writer.writerow({
                    "id": "bad", "amino_acid_substitutions": "A1*",
                    "aa_sequence": "AC*E", "fitness_value": "2.0",
                    "reaction_smiles": "CCO>>CC=O",
                })
            result = read_enzengdb_directory(directory)
        self.assertEqual(result.audit.total_rows, 2)
        self.assertEqual(result.audit.accepted_rows, 1)
        self.assertEqual(result.audit.campaigns, 1)
        observation = result.observations[0]
        self.assertEqual(observation.record.enzyme_type, "parent")
        self.assertEqual(observation.record.substrate_smiles, "CCO")
        self.assertIsNone(observation.record.log10_activity)
        self.assertIn("non-canonical", result.audit.rejected[0].reason)


def observation(campaign, fitness, sequence, candidate_id):
    return EnzEngDBObservation(
        record=EnzymeSubstrateRecord(
            sequence=sequence,
            substrate_smiles="CCO",
            organism="",
            ec="",
            enzyme_type="variant",
            candidate_id=candidate_id,
        ),
        campaign_id=campaign,
        fitness_value=fitness,
        endpoint_name="fitness",
        source_file=f"{campaign}.csv",
    )


class CampaignRecordTest(unittest.TestCase):
    def test_average_tie_percentile_ranks(self):
        observations = [
            observation("a", 1.0, "AAA", "low"),
            observation("a", 2.0, "AAC", "tie-1"),
            observation("a", 2.0, "AAG", "tie-2"),
            observation("a", 4.0, "AAT", "high"),
            observation("single", 8.0, "ACC", "only"),
        ]
        ranked = campaign_rank_records(observations)
        self.assertEqual([item.activity_rank for item in ranked], [0.0, 0.5, 0.5, 1.0, 0.5])
        self.assertEqual([item.campaign_group for item in ranked], ["a"] * 4 + ["single"])

    def test_representative_is_actual_consensus_nearest_and_deterministic(self):
        observations = [
            observation("b", 0.0, "AAAA", "wrong-length"),
            observation("a", 0.0, "AAA", "first"),
            observation("a", 0.0, "AAC", "consensus"),
            observation("a", 0.0, "ACC", "last"),
            observation("b", 0.0, "CCC", "modal-1"),
            observation("b", 0.0, "CCA", "modal-2"),
        ]
        forward = campaign_representative_records(observations)
        reverse = campaign_representative_records(reversed(observations))
        self.assertEqual(forward, reverse)
        self.assertEqual([item.campaign_group for item in forward], ["a", "b"])
        self.assertEqual(forward[0].candidate_id, "consensus")
        self.assertEqual(len(forward[1].sequence), 3)


class CampaignRankingMetricsTest(unittest.TestCase):
    def test_perfect_order_has_perfect_metrics(self):
        metrics = campaign_ranking_metrics([1, 2, 3, 4], [10, 20, 30, 40])
        self.assertAlmostEqual(metrics["spearman"], 1.0)
        self.assertAlmostEqual(metrics["pairwise_accuracy"], 1.0)
        self.assertAlmostEqual(metrics["top_10pct_recall"], 1.0)

    def test_reverse_order_is_negative(self):
        metrics = campaign_ranking_metrics([4, 3, 2, 1], [10, 20, 30, 40])
        self.assertAlmostEqual(metrics["spearman"], -1.0)
        self.assertAlmostEqual(metrics["pairwise_accuracy"], 0.0)

    def test_prediction_ties_receive_random_expectation(self):
        metrics = campaign_ranking_metrics([1, 1, 1, 1], [10, 20, 30, 40])
        self.assertIsNone(metrics["spearman"])
        self.assertAlmostEqual(metrics["pairwise_accuracy"], 0.5)
        self.assertAlmostEqual(metrics["top_10pct_enrichment"], 1.0)


if __name__ == "__main__":
    unittest.main()
