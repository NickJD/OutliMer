import importlib.util
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
DATA = os.path.join(ROOT, "tests", "data")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from OutliMer import classification


def _require_pandas():
    if importlib.util.find_spec("pandas") is None:
        raise unittest.SkipTest("pandas is not installed")
    import pandas as pd
    return pd


class ClassificationTests(unittest.TestCase):
    def test_combine_scores_ties_are_not_all_anomalies(self):
        pd = _require_pandas()
        scores = pd.Series([0.0, 0.0], index=["s1", "s2"])

        combined = classification.combine_scores(scores, scores, scores)

        self.assertEqual(combined["anomaly_score"].tolist(), [0.0, 0.0])

        flagged = classification.flag_anomalies(combined, contamination=0.5)
        self.assertFalse(flagged["is_anomaly"].any())

    def test_anomaly_flags_follow_contamination(self):
        pd = _require_pandas()
        ranked = pd.DataFrame(
            {"anomaly_score": [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]},
            index=[f"s{i}" for i in range(10)],
        )
        flagged = classification.flag_anomalies(ranked, contamination=0.2)

        self.assertEqual(int(flagged["is_anomaly"].sum()), 2)
        self.assertEqual(flagged.index[flagged["is_anomaly"]].tolist(), ["s0", "s1"])

    def test_report_mode_does_not_require_union_csv(self):
        pd = _require_pandas()
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "kmer_report.csv")
            pd.DataFrame({
                "sample": ["s1", "s2"],
                "n_hashes": [100, 100],
                "n_in_db": [95, 20],
                "pct_in_db": [0.95, 0.20],
                "n_new_hashes": [5, 80],
            }).to_csv(report_path, index=False)

            rc = classification.main([
                "--mode", "report",
                "--report-csv", report_path,
                "--out-dir", tmpdir,
            ])

            self.assertEqual(rc, 0)
            self.assertTrue(
                os.path.exists(
                    os.path.join(tmpdir, "outliers_report_from_reportcsv.csv")
                )
            )

    def test_enrichment_uses_explicit_foreground_and_background(self):
        _require_pandas()
        old_loader = classification._load_fisher_exact

        def fake_loader():
            def fake_fisher(table):
                foreground_present = table[0][0]
                background_present = table[1][0]
                pvalue = 0.01 if foreground_present > background_present else 1.0
                return 1.0, pvalue
            return fake_fisher

        classification._load_fisher_exact = fake_loader
        try:
            result = classification.compute_enrichment(
                {
                    "s1": {1: 2, 2: 1},
                    "s2": {1: 1},
                    "s3": {2: 1},
                    "s4": {3: 1},
                },
                foreground_samples=["s1", "s2"],
                background_samples=["s3", "s4"],
            )
        finally:
            classification._load_fisher_exact = old_loader

        self.assertEqual(result.index[0], 1)
        self.assertIn("p_adj", result.columns)
        self.assertEqual(result.loc[1, "count_in_foreground"], 2)
        self.assertEqual(result.loc[1, "count_in_background"], 0)

    def test_metadata_query_selects_samples(self):
        pd = _require_pandas()
        metadata = pd.DataFrame({
            "status": ["case", "control", "case"],
        }, index=["s1", "s2", "s3"])

        selected = classification.select_samples_by_query(
            metadata, "status == 'case'"
        )

        self.assertEqual(selected, ["s1", "s3"])

    def test_metadata_loader_rejects_duplicate_samples(self):
        _require_pandas()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "metadata.tsv")
            with open(path, "w") as fh:
                fh.write("sample\tstatus\ns1\tcase\ns1\tcontrol\n")
            with self.assertRaisesRegex(ValueError, "duplicate metadata"):
                classification.load_sample_metadata(path)

    def test_union_loader_rejects_invalid_matrices(self):
        _require_pandas()
        with tempfile.TemporaryDirectory() as tmpdir:
            negative = os.path.join(tmpdir, "negative.csv")
            with open(negative, "w") as fh:
                fh.write("hash,s1\n1,-1\n")
            with self.assertRaisesRegex(ValueError, "non-negative"):
                classification.load_union_csv(negative)

            duplicate = os.path.join(tmpdir, "duplicate.csv")
            with open(duplicate, "w") as fh:
                fh.write("hash,s1,s1\n1,1,2\n")
            with self.assertRaisesRegex(ValueError, "duplicate union CSV"):
                classification.load_union_csv(duplicate)

            malformed = os.path.join(tmpdir, "malformed.csv")
            with open(malformed, "w") as fh:
                fh.write("hash,s1\n1,abc\n")
            with self.assertRaisesRegex(ValueError, "numeric"):
                classification.load_union_csv(malformed)

    def test_compare_union_matrices(self):
        pd = _require_pandas()
        with tempfile.TemporaryDirectory() as tmpdir:
            baseline = os.path.join(tmpdir, "baseline.csv")
            query = os.path.join(tmpdir, "query.csv")
            pd.DataFrame({"b1": [1, 1, 0]}, index=[1, 2, 3]).to_csv(
                baseline, index_label="hash"
            )
            pd.DataFrame({"q1": [1, 0, 1]}, index=[1, 2, 3]).to_csv(
                query, index_label="hash"
            )

            compared = classification.compare_union_matrices(baseline, query)

        self.assertEqual(int(compared.loc["q1", "n_hashes"]), 2)
        self.assertEqual(int(compared.loc["q1", "n_in_baseline"]), 1)
        self.assertAlmostEqual(float(compared.loc["q1", "pct_in_baseline"]), 0.5)

    def test_multiqc_and_rocrate_writers(self):
        pd = _require_pandas()
        ranked = pd.DataFrame({"anomaly_score": [0.8]}, index=["s1"])
        with tempfile.TemporaryDirectory() as tmpdir:
            mqc = classification.write_multiqc_custom_content(
                tmpdir, ranked, "anomaly_score"
            )
            crate = classification.write_ro_crate(
                tmpdir,
                "test run",
                {"input": "input.csv"},
                {"multiqc": mqc},
                {"mode": "test"},
            )

            self.assertTrue(os.path.exists(mqc))
            self.assertTrue(os.path.exists(crate))
            with open(mqc) as fh:
                payload = __import__("json").load(fh)

        self.assertEqual(payload["id"], "outlimer_anomaly_scores")
        self.assertIn("s1", payload["data"])

    def test_sparse_npz_writer_when_scipy_available(self):
        if importlib.util.find_spec("scipy") is None:
            raise unittest.SkipTest("scipy is not installed")
        pd = _require_pandas()
        df = pd.DataFrame({"s1": [1, 0], "s2": [0, 2]}, index=[11, 22])
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "matrix.npz")
            classification.write_sparse_npz(df, out)

            self.assertTrue(os.path.exists(out))
            self.assertTrue(os.path.exists(out + ".labels.json"))

    def test_sample_qc_metrics_counts_unique_hashes(self):
        pd = _require_pandas()
        df = pd.DataFrame({
            "s1": [5, 0, 1],
            "s2": [0, 2, 1],
            "s3": [0, 0, 1],
        }, index=[11, 22, 33])
        ranked = pd.DataFrame({"anomaly_score": [0.8, 0.1, 0.2]},
                              index=["s1", "s2", "s3"])

        qc = classification.compute_sample_qc_metrics(df, ranked)

        self.assertEqual(int(qc.loc["s1", "n_hashes"]), 2)
        self.assertEqual(int(qc.loc["s1", "unique_hashes"]), 1)
        self.assertAlmostEqual(float(qc.loc["s1", "pct_unique"]), 0.5)
        self.assertAlmostEqual(float(qc.loc["s1", "anomaly_score"]), 0.8)

    def test_feature_matrix_normalises_depth_and_selects_variability(self):
        pd = _require_pandas()
        df = pd.DataFrame({
            "s1": [1000, 10, 0],
            "s2": [2000, 20, 0],
            "s3": [1000, 10, 100],
        }, index=[11, 22, 33])

        _, _, all_log = classification.prepare_feature_matrix(df, top_M=0)
        selected, _, _ = classification.prepare_feature_matrix(df, top_M=1)

        self.assertAlmostEqual(float(all_log.loc["s1", 11]),
                               float(all_log.loc["s2", 11]))
        self.assertEqual(selected.columns.tolist(), [33])

    def test_feature_matrix_rejects_empty_samples(self):
        pd = _require_pandas()
        df = pd.DataFrame({"s1": [1, 2], "empty": [0, 0]}, index=[11, 22])

        with self.assertRaisesRegex(ValueError, "empty"):
            classification.prepare_feature_matrix(df, top_M=0)

    def test_mean_jaccard_excludes_self_distance(self):
        pd = _require_pandas()
        binary = pd.DataFrame({1: [1, 0], 2: [0, 1]}, index=["s1", "s2"])

        mean = classification.compute_mean_jaccard_distance(binary)

        self.assertEqual(mean.to_dict(), {"s1": 1.0, "s2": 1.0})

    def test_truth_fixture_separates_composition_from_depth(self):
        _require_pandas()
        df = classification.load_union_csv(
            os.path.join(DATA, "classifier_truth.csv")
        )
        _, binary, normalised = classification.prepare_feature_matrix(df, top_M=0)
        combined = classification.flag_anomalies(
            classification.combine_scores(
                classification.compute_mean_jaccard_distance(binary),
                classification.compute_isolation_forest(normalised, 0.1),
                classification.compute_lof(normalised, 0.1),
            ),
            contamination=0.1,
        )

        self.assertEqual(combined.index[0], "composition_outlier")
        self.assertTrue(bool(combined.loc["composition_outlier", "is_anomaly"]))
        self.assertFalse(bool(combined.loc["depth_control", "is_anomaly"]))

    def test_diagnostic_plot_writers_when_dependencies_available(self):
        required = ["matplotlib", "seaborn", "sklearn", "scipy"]
        missing = [name for name in required if importlib.util.find_spec(name) is None]
        if missing:
            raise unittest.SkipTest(
                "plot dependencies are not installed: " + ", ".join(missing)
            )
        pd = _require_pandas()
        df = pd.DataFrame({
            "s1": [5, 0, 1, 0, 4],
            "s2": [0, 2, 1, 0, 0],
            "s3": [0, 0, 1, 7, 0],
            "s4": [1, 1, 0, 0, 0],
        }, index=[11, 22, 33, 44, 55])
        X, X_binary, X_log = classification.prepare_feature_matrix(df, top_M=0)
        scores = pd.Series([0.7, 0.2, 0.9, 0.1], index=X.index)
        combined = classification.combine_scores(scores, scores, scores)
        qc = classification.compute_sample_qc_metrics(df, combined)
        enrichment = pd.DataFrame({
            "pvalue": [0.001, 0.02, 0.5],
            "p_adj": [0.003, 0.03, 0.5],
        }, index=[11, 22, 33])

        with tempfile.TemporaryDirectory() as tmpdir:
            paths = {
                "pca": os.path.join(tmpdir, "pca.png"),
                "mds": os.path.join(tmpdir, "mds.png"),
                "qc": os.path.join(tmpdir, "qc.png"),
                "cluster": os.path.join(tmpdir, "cluster.png"),
                "driver": os.path.join(tmpdir, "driver.png"),
                "enrichment": os.path.join(tmpdir, "enrichment.png"),
            }
            classification.compute_pca_plot(X_log, paths["pca"], combined)
            classification.compute_jaccard_mds_plot(X_binary, paths["mds"], combined)
            classification.compute_anomaly_qc_plot(qc, paths["qc"])
            classification.compute_clustered_jaccard_heatmap(X_binary, paths["cluster"])
            classification.compute_driver_hash_heatmap(df, combined, paths["driver"])
            classification.compute_enrichment_plot(enrichment, paths["enrichment"])

            for path in paths.values():
                self.assertTrue(os.path.exists(path), path)
                self.assertGreater(os.path.getsize(path), 0)


if __name__ == "__main__":
    unittest.main()
