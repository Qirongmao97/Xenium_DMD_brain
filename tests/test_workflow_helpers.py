"""Smoke tests for data handling without Scanpy, Seurat, or study data."""

import gzip
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd
from scipy.io import mmread

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from allen_reference import class_sample
from export_seurat_matrix import export_matrix
from isoform_coexpression import LONG, SHORT, summarize
from isoform_expression import sample_tests, summarize as sample_mean
from skipping_differential_expression import matched_draw, negative_pools


class WorkflowHelpers(unittest.TestCase):
    def test_codetection_denominator(self):
        columns = LONG + SHORT
        values = pd.DataFrame(
            [[1, 0, 0, 1, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
            columns=columns,
        )
        meta = pd.DataFrame({"Region": ["CA1"] * 4, "Class": ["Glut"] * 4})
        out = summarize(values, meta).loc[("CA1", "Glut")]
        self.assertEqual(out["Coexpression_Either"], 3)
        self.assertAlmostEqual(out["percent_coexpressed"], 100 / 3)
        self.assertAlmostEqual(out["percent_A"] + out["percent_B"] + out["percent_coexpressed"], 100)
        self.assertEqual(out["short_only_display"], 0)

    def test_sample_means_do_not_pool_cells_across_replicates(self):
        df = pd.DataFrame(
            {"Sample": ["w1", "w1", "w2", "m1", "m2"],
             "Group": ["WT", "WT", "WT", "mdx52", "mdx52"],
             "Dp427c": [1, 3, 4, 7, 9]}
        )
        means = sample_mean(df, [])
        result = sample_tests(means, ["Dp427c"]).iloc[0]
        self.assertEqual(result["n_WT"], 2)
        self.assertEqual(result["n_mdx52"], 2)
        self.assertAlmostEqual(result["mdx52_minus_WT"], 5)

    def test_control_matches_disjoint_and_prefer_exact_tier(self):
        obs = pd.DataFrame({"Class": ["A", "B", "A", "A", "C"], "Region": ["r1", "r2", "r1", "r2", "r2"]})
        negatives = np.array([2, 3, 4])
        usage = np.array([0, 0, 0, 1, 0])
        indices, diagnostic = matched_draw(
            {("A", "r1"): 1, ("B", "r2"): 1},
            negative_pools(obs, negatives, ["Class", "Region"]),
            negative_pools(obs, negatives, "Class"), negative_pools(obs, negatives, "Region"),
            negatives, usage, 1, np.random.default_rng(42),
        )
        np.testing.assert_array_equal(indices, np.array([2, 4]))
        self.assertEqual(diagnostic, dict(exact=1, class_only=0, region_only=1, global_pool=0))

    def test_capped_controls_can_fail_without_changing_usage(self):
        usage = np.array([0, 1])
        indices, _ = matched_draw(
            {("A", "r"): 1}, {("A", "r"): np.array([1])}, {}, {},
            np.array([1]), usage, 1, np.random.default_rng(42),
        )
        self.assertIsNone(indices)
        np.testing.assert_array_equal(usage, np.array([0, 1]))

    def test_matrix_market_orientation_and_layer(self):
        data = SimpleNamespace(
            X=np.array([[1, 2, 3], [4, 5, 6]]),
            layers={"counts": np.array([[1, 0, 2], [0, 0, 7]])},
            obs_names=pd.Index(["cell_a", "cell_b"]), var_names=pd.Index(["f1", "f2", "f3"]),
        )
        with TemporaryDirectory() as tmp:
            out = Path(tmp)
            export_matrix(data, out, layer="counts")
            with gzip.open(out / "matrix.mtx.gz", "rb") as handle:
                np.testing.assert_array_equal(mmread(handle).toarray(), data.layers["counts"].T)
            with gzip.open(out / "barcodes.tsv.gz", "rt") as handle:
                self.assertEqual(handle.read().splitlines(), list(data.obs_names))

    def test_class_sampling_uses_only_10xv2(self):
        metadata = pd.DataFrame(
            {"library_method": ["10Xv2"] * 20 + ["10Xv3"] * 10, "class": ["A"] * 10 + ["B"] * 20}
        )
        test = class_sample(metadata, 0.1, 42)
        self.assertEqual(test["class"].value_counts().to_dict(), {"A": 1, "B": 1})
        pd.testing.assert_frame_equal(test, class_sample(metadata, 0.1, 42))


if __name__ == "__main__":
    unittest.main()
