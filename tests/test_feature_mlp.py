import unittest

import numpy as np
import torch

from biocandidate.feature_mlp import FeatureMLP, standardize_features, train_feature_mlp


class FeatureMLPTest(unittest.TestCase):
    def test_standardization_uses_train_statistics_and_handles_constant_columns(self):
        train = np.asarray([[1.0, 2.0], [3.0, 2.0]], dtype=np.float32)
        validation = np.asarray([[5.0, 2.0]], dtype=np.float32)
        scaled, mean, stddev = standardize_features(train, validation, validation)
        self.assertTrue(np.allclose(mean, [2.0, 2.0]))
        self.assertTrue(np.allclose(stddev, [1.0, 1.0]))
        self.assertTrue(np.allclose(scaled["validation"], [[3.0, 0.0]]))

    def test_training_selects_finite_checkpoint(self):
        x = np.arange(24, dtype=np.float32).reshape(12, 2)
        y = x[:, 0] * 0.1
        features = {"train": x[:8], "validation": x[8:10], "test": x[10:]}
        labels = {"train": y[:8], "validation": y[8:10], "test": y[10:]}
        model, metrics, _, _ = train_feature_mlp(
            features, labels, seed=7, epochs=2, batch_size=4, learning_rate=0.01,
            weight_decay=0.0, hidden_dimensions=(4,), dropout=0.0,
            device=torch.device("cpu"),
        )
        self.assertIsInstance(model, FeatureMLP)
        self.assertIn(metrics["selected_epoch"], (1, 2))
        self.assertTrue(np.isfinite(metrics["test_metrics"]["rmse"]))


if __name__ == "__main__":
    unittest.main()
