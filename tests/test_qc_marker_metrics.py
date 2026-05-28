from __future__ import annotations

import unittest

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from nico_wrapper.qc import (
    marker_annotation_qc,
    marker_de_recovery_metrics,
    marker_logfc_metrics,
    marker_score_metrics,
    marker_set_scores,
    marker_set_summary,
)


class MarkerMetricTests(unittest.TestCase):
    def _adata(self, *, sparse_x: bool = False) -> AnnData:
        x = np.array(
            [
                [5, 4, 0, 0, 1, 1],
                [6, 5, 0, 0, 1, 1],
                [5, 4, 0, 0, 1, 1],
                [0, 0, 5, 4, 1, 1],
                [0, 0, 6, 5, 1, 1],
                [0, 0, 5, 4, 1, 1],
            ],
            dtype=float,
        )
        if sparse_x:
            x = sparse.csr_matrix(x)
        return AnnData(
            x,
            obs=pd.DataFrame({"cell_type": ["A", "A", "A", "B", "B", "B"]}),
            var=pd.DataFrame(index=["g1", "g2", "g3", "g4", "g5", "g6"]),
        )

    def test_summary_and_logfc_edge_cases(self) -> None:
        adata = self._adata()
        markers = {"A": ["g1", "g1", "missing"], "B": ["g3", "g4"], "C": ["g2"]}

        summary = marker_set_summary(adata, markers, labels=["B", "A", "C"])
        self.assertEqual(summary["cell_type"].tolist(), ["B", "A", "C"])
        self.assertEqual(summary.loc[1, "n_markers_provided"], 2)
        self.assertEqual(summary.loc[1, "n_markers_used"], 1)
        self.assertEqual(summary.loc[1, "missing_marker_names"], ["missing"])
        self.assertEqual(summary.loc[2, "n_cells"], 0)

        logfc, per_marker = marker_logfc_metrics(
            adata,
            markers,
            labels=["A", "B", "C"],
            return_per_marker=True,
        )
        self.assertGreater(logfc.loc[0, "mean_marker_logFC"], 0)
        self.assertGreater(logfc.loc[1, "mean_marker_logFC"], 0)
        self.assertTrue(np.isnan(logfc.loc[2, "mean_marker_logFC"]))
        self.assertEqual(per_marker["cell_type"].tolist(), ["A", "B", "B", "C"])
        self.assertTrue(
            per_marker.loc[per_marker["cell_type"] == "C", "marker_logFC"].isna().all()
        )

    def test_sparse_logfc_matches_dense(self) -> None:
        markers = {"A": ["g1", "g2"], "B": ["g3", "g4"]}
        dense = marker_logfc_metrics(self._adata(), markers)
        sparse_result = marker_logfc_metrics(self._adata(sparse_x=True), markers)
        np.testing.assert_allclose(
            dense["mean_marker_logFC"].to_numpy(),
            sparse_result["mean_marker_logFC"].to_numpy(),
        )

    def test_raw_count_like_x_is_rejected_unless_raw_is_selected(self) -> None:
        raw = np.array(
            [
                [500, 400, 0, 0],
                [600, 500, 0, 0],
                [0, 0, 700, 600],
                [0, 0, 800, 700],
            ],
            dtype=float,
        )
        adata = AnnData(
            raw,
            obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
            var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
        )
        markers = {"A": ["g1", "g2"], "B": ["g3", "g4"]}

        with self.assertRaisesRegex(
            ValueError, "adata.X appears to contain raw counts"
        ):
            marker_logfc_metrics(adata, markers)

        normalized = adata.copy()
        normalized.X = np.log1p(raw)
        result = marker_logfc_metrics(normalized, markers)
        self.assertEqual(result["n_markers_used"].tolist(), [2, 2])

        raw_selected = normalized.copy()
        raw_selected.raw = AnnData(
            raw,
            obs=raw_selected.obs.copy(),
            var=raw_selected.var.copy(),
        )
        raw_result = marker_logfc_metrics(raw_selected, markers, use_raw=True)
        self.assertEqual(raw_result["n_markers_used"].tolist(), [2, 2])

    def test_raw_count_like_layer_is_rejected_for_scores(self) -> None:
        import scanpy as sc

        original_score_genes = sc.tl.score_genes

        def fake_score_genes(adata, gene_list, *, score_name, **kwargs):
            adata.obs[score_name] = np.zeros(adata.n_obs)

        sc.tl.score_genes = fake_score_genes
        try:
            adata = self._adata()
            adata.layers["counts"] = np.array(
                [
                    [500, 400, 0, 0, 100, 100],
                    [600, 500, 0, 0, 100, 100],
                    [500, 400, 0, 0, 100, 100],
                    [0, 0, 700, 600, 100, 100],
                    [0, 0, 800, 700, 100, 100],
                    [0, 0, 700, 600, 100, 100],
                ],
                dtype=float,
            )
            with self.assertRaisesRegex(
                ValueError, "adata.layers\\['counts'\\] appears to contain raw counts"
            ):
                marker_set_scores(adata, {"A": ["g1"], "B": ["g3"]}, layer="counts")
        finally:
            sc.tl.score_genes = original_score_genes

    def test_plain_string_marker_set_is_rejected(self) -> None:
        adata = self._adata()
        with self.assertRaisesRegex(ValueError, "not a single string"):
            marker_set_summary(adata, {"A": "g1"})

    def test_gene_symbol_resolution_ignores_ambiguous_symbols(self) -> None:
        adata = AnnData(
            np.ones((2, 4)),
            obs=pd.DataFrame({"cell_type": ["A", "B"]}),
            var=pd.DataFrame(
                {"symbol": ["S1", "DUP", "DUP", "S4"]},
                index=["ENS1", "ENS2", "ENS3", "ENS4"],
            ),
        )
        summary = marker_set_summary(
            adata,
            {"A": ["S1", "DUP", "ENS3"]},
            gene_symbols_key="symbol",
        )
        self.assertEqual(summary.loc[0, "n_markers_provided"], 3)
        self.assertEqual(summary.loc[0, "n_markers_used"], 2)
        self.assertEqual(summary.loc[0, "missing_marker_names"], ["DUP"])

    def test_scores_do_not_mutate_unless_requested(self) -> None:
        import scanpy as sc

        original_score_genes = sc.tl.score_genes

        def fake_score_genes(adata, gene_list, *, score_name, **kwargs):
            positions = [adata.var_names.get_loc(gene) for gene in gene_list]
            adata.obs[score_name] = np.asarray(adata.X[:, positions]).mean(axis=1)

        sc.tl.score_genes = fake_score_genes
        try:
            adata = self._adata()
            markers = {"A": ["g1"], "B": ["g3"]}
            scores = marker_set_scores(adata, markers)
            self.assertEqual(scores.columns.tolist(), ["A", "B"])
            self.assertFalse(
                any(column.startswith("marker_score__") for column in adata.obs)
            )

            marker_set_scores(adata, markers, copy_scores_to_obs=True)
            self.assertIn("marker_score__A", adata.obs)
            self.assertIn("marker_score__B", adata.obs)

            metrics = marker_score_metrics(adata, markers)
            self.assertEqual(metrics["correct_top_score_fraction"].tolist(), [1.0, 1.0])
        finally:
            sc.tl.score_genes = original_score_genes

    def test_scores_handle_duplicate_var_names_without_mutating_them(self) -> None:
        import scanpy as sc

        original_score_genes = sc.tl.score_genes

        def fake_score_genes(adata, gene_list, *, score_name, **kwargs):
            self.assertTrue(adata.var_names.is_unique)
            positions = [adata.var_names.get_loc(gene) for gene in gene_list]
            adata.obs[score_name] = np.asarray(adata.X[:, positions]).mean(axis=1)

        sc.tl.score_genes = fake_score_genes
        try:
            adata = AnnData(
                np.array([[1, 10, 0], [0, 0, 2]], dtype=float),
                obs=pd.DataFrame({"cell_type": ["A", "B"]}),
                var=pd.DataFrame(index=["g1", "g1", "g2"]),
            )
            self.assertFalse(adata.var_names.is_unique)
            scores = marker_set_scores(adata, {"A": ["g1"], "B": ["g2"]})
            self.assertEqual(scores["A"].tolist(), [1.0, 0.0])
            self.assertFalse(adata.var_names.is_unique)
        finally:
            sc.tl.score_genes = original_score_genes

    def test_de_recovery_and_final_column_order(self) -> None:
        too_few = AnnData(
            np.array([[1, 0], [0, 1], [0, 2]], dtype=float),
            obs=pd.DataFrame({"cell_type": ["A", "B", "B"]}),
            var=pd.DataFrame(index=["g1", "g2"]),
        )
        too_few_summary = marker_de_recovery_metrics(
            too_few,
            {"A": ["g1"], "B": ["g2"]},
            method="t-test",
        )
        self.assertTrue(too_few_summary["fraction_markers_significant_DE"].isna().all())

        duplicate_var_adata = AnnData(
            np.array([[1, 10], [1, 9], [1, 0], [1, 0]], dtype=float),
            obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
            var=pd.DataFrame(index=["g1", "g1"]),
        )
        duplicate_de = marker_de_recovery_metrics(
            duplicate_var_adata,
            {"A": ["g1"], "B": []},
            method="t-test",
            top_n=1,
        )
        self.assertEqual(duplicate_de.loc[0, "marker_recall_top1_DE"], 0.0)

        adata = self._adata()
        markers = {"A": ["g1", "g2"], "B": ["g3", "g4"]}
        de_summary, de_table = marker_de_recovery_metrics(
            adata,
            markers,
            method="t-test",
            top_n=2,
            return_de_table=True,
        )
        self.assertIn("marker_recall_top2_DE", de_summary)
        self.assertIn("is_marker_for_group", de_table)
        self.assertEqual(de_summary["n_markers_used"].tolist(), [2, 2])

        import scanpy as sc

        original_score_genes = sc.tl.score_genes

        def fake_score_genes(adata, gene_list, *, score_name, **kwargs):
            positions = [adata.var_names.get_loc(gene) for gene in gene_list]
            adata.obs[score_name] = np.asarray(adata.X[:, positions]).mean(axis=1)

        sc.tl.score_genes = fake_score_genes
        try:
            final = marker_annotation_qc(
                adata,
                markers,
                score_ctrl_size=1,
                score_n_bins=2,
                de_method="t-test",
                de_top_n=2,
            )
        finally:
            sc.tl.score_genes = original_score_genes

        self.assertEqual(
            final.columns.tolist(),
            [
                "cell_type",
                "n_cells",
                "n_markers_provided",
                "n_markers_used",
                "mean_marker_logFC",
                "fraction_markers_logFC_gt_0.25",
                "correct_top_score_fraction",
                "median_score_margin",
                "fraction_markers_significant_DE",
                "marker_recall_top2_DE",
            ],
        )


if __name__ == "__main__":
    unittest.main()
