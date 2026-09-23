import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location('paired', Path(__file__).resolve().parents[1]/'scripts/paired_family_analysis.py')
paired = importlib.util.module_from_spec(spec)
spec.loader.exec_module(paired)

class PairedFamilyTests(unittest.TestCase):
    def test_whole_cluster_weights_match_explicit_rows(self):
        y = np.array([0., 1., 3., 5.]); p = np.array([1., 1., 2., 6.])
        groups = np.array([0, 0, 1, 1])
        actual = paired.metrics_by_draw([p], y, groups, np.array([[2., 1.]]))[0]
        idx = [0, 1, 0, 1, 2, 3]
        expected = [np.sqrt(np.mean((p[idx]-y[idx])**2)), np.mean(abs(p[idx]-y[idx])), np.corrcoef(p[idx], y[idx])[0,1]]
        np.testing.assert_allclose(actual, expected)

    def test_alignment_and_label_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'p.csv'
            path.write_text('pair_id,p,y\nb,4,2\na,3,1\n')
            np.testing.assert_equal(paired.aligned(path,['a','b'],np.array([1.,2.]),'p','y'),[3,4])
            with self.assertRaises(ValueError):
                paired.aligned(path,['a','b'],np.array([1.,9.]),'p','y')

    def test_seed_average_is_not_prediction_ensemble(self):
        y = np.array([0.,1.,2.]); predictions = [y-1, y+1]
        actual = paired.metrics_by_draw(predictions,y,np.arange(3),np.ones((1,3)))
        self.assertAlmostEqual(actual[0,0], 1.)

if __name__ == '__main__':
    unittest.main()
