"""Isolated GSEApy binding and real NiCo pathway filtering regressions."""

from dataclasses import asdict
from importlib.metadata import version
import os
from pathlib import Path
import pickle
import tempfile
from types import FunctionType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch, sentinel

import numpy as np
import pandas as pd

from nico_wrapper.covariation.config import CovariationReportConfig
from nico_wrapper.covariation.reports import _adapt_nico_pathway_analysis, run_pathway_enrichment
from nico_wrapper.covariation.results import CovariationResult


def _make_probe(gseapy_namespace):
    """Give the probe its own GSEApy global, defaults, and a nonempty closure."""

    background = object()
    token_default = object()

    def probe(token=token_default, *, organism="Mouse", gene_list, gene_sets="ProbeDB", cutoff=0.03):
        # This global is supplied in the private globals dictionary below.
        result = gseapy.enrichr(
            token, organism=organism, gene_list=gene_list, gene_sets=gene_sets,
            cutoff=cutoff, background=background,
        )
        gseapy.dotplot(result, marker=background)
        return result

    original = FunctionType(
        probe.__code__, {"gseapy": gseapy_namespace, "other_global": object()},
        probe.__name__, probe.__defaults__, probe.__closure__,
    )
    original.__kwdefaults__ = probe.__kwdefaults__
    return original, background


class GseaAdapterTests(unittest.TestCase):
    def setUp(self):
        self.gseapy = SimpleNamespace(enrichr=Mock(), dotplot=Mock(), untouched=sentinel.attribute)
        self.original, self.background = _make_probe(self.gseapy)
        self.original_globals = self.original.__globals__
        self.globals_before = dict(self.original_globals)
        self.gseapy_before = dict(vars(self.gseapy))
        self.adapted = _adapt_nico_pathway_analysis(self.original)

    def assert_shared_bindings_unchanged(self):
        self.assertIs(self.original.__globals__, self.original_globals)
        self.assertEqual(self.original.__globals__, self.globals_before)
        self.assertEqual(vars(self.gseapy), self.gseapy_before)

    def test_code_defaults_closure_and_other_bindings_are_preserved(self):
        self.assertIs(self.adapted.__code__, self.original.__code__)
        self.assertEqual(self.adapted.__name__, self.original.__name__)
        self.assertIs(self.adapted.__defaults__, self.original.__defaults__)
        self.assertIs(self.adapted.__kwdefaults__, self.original.__kwdefaults__)
        self.assertIsNotNone(self.original.__closure__)
        self.assertIs(self.adapted.__closure__, self.original.__closure__)
        self.assertIsNot(self.adapted.__globals__, self.original_globals)
        private_gseapy = self.adapted.__globals__["gseapy"]
        self.assertIsNot(private_gseapy, self.gseapy)
        self.assertIsNot(private_gseapy.enrichr, self.gseapy.enrichr)
        self.assertIs(private_gseapy.dotplot, self.gseapy.dotplot)
        self.assertIs(private_gseapy.untouched, self.gseapy.untouched)
        self.assertIs(self.adapted.__globals__["other_global"], self.original_globals["other_global"])

        def enrichr(*args, **kwargs):
            self.assert_shared_bindings_unchanged()
            return sentinel.enriched

        self.gseapy.enrichr.side_effect = enrichr
        genes = ["Rps3", "MT-CO1", "MixedCaseGene"]
        result = self.adapted(gene_list=genes)
        self.assertIs(result, sentinel.enriched)
        self.gseapy.enrichr.assert_called_once_with(
            self.original.__defaults__[0], organism="mouse", gene_list=genes,
            gene_sets="ProbeDB", cutoff=0.03, background=self.background,
        )
        self.gseapy.dotplot.assert_called_once_with(sentinel.enriched, marker=self.background)
        self.assert_shared_bindings_unchanged()

    def test_organism_is_the_only_translated_argument(self):
        for upstream, expected in (("Mouse", "mouse"), ("Human", "human")):
            with self.subTest(organism=upstream):
                self.gseapy.enrichr.reset_mock()
                self.gseapy.dotplot.reset_mock()
                genes = ["Rps3", "RPL5", "mt-Co1", "MT-CO1", "MixedCaseGene"]
                databases = ["GO_Biological_Process_2021", "Reactome_2016"]
                token = object()

                def enrichr(*args, **kwargs):
                    self.assert_shared_bindings_unchanged()
                    self.assertIs(args[0], token)
                    self.assertIs(kwargs["gene_list"], genes)
                    self.assertIs(kwargs["gene_sets"], databases)
                    return sentinel.enriched

                self.gseapy.enrichr.side_effect = enrichr
                result = self.adapted(
                    token, organism=upstream, gene_list=genes, gene_sets=databases, cutoff=0.17
                )
                self.assertIs(result, sentinel.enriched)
                self.gseapy.enrichr.assert_called_once_with(
                    token, organism=expected, gene_list=genes, gene_sets=databases,
                    cutoff=0.17, background=self.background,
                )
                self.gseapy.dotplot.assert_called_once_with(sentinel.enriched, marker=self.background)
                self.assert_shared_bindings_unchanged()
                # A call outside the adapter still supplies title case unchanged.
                self.original(token, organism=upstream, gene_list=genes, gene_sets=databases, cutoff=0.17)
                self.assertEqual(self.gseapy.enrichr.call_args.kwargs["organism"], upstream)
                self.assert_shared_bindings_unchanged()

    def test_enrichment_exception_propagates_without_changing_shared_bindings(self):
        error = RuntimeError("enrichment request failed")

        def fail(*args, **kwargs):
            self.assert_shared_bindings_unchanged()
            self.assertEqual(kwargs["organism"], "human")
            raise error

        self.gseapy.enrichr.side_effect = fail
        with self.assertRaises(RuntimeError) as raised:
            self.adapted(organism="Human", gene_list=["ACTB"])
        self.assertIs(raised.exception, error)
        self.gseapy.enrichr.assert_called_once()
        self.gseapy.dotplot.assert_not_called()
        self.assert_shared_bindings_unchanged()


class RealNicoPathwayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import matplotlib

        matplotlib.use("Agg")
        from nico import Covariations as scov

        cls.scov = scov

    def test_real_filtering_order_and_barplot_for_both_species_and_inclusion_settings(self):
        self.assertEqual(version("nico-sc-sp"), "1.8.0")
        self.assertEqual(version("gseapy"), "1.3.1")
        for organism, excluded, ordinary in (
            ("mouse", ["Rps3", "Rpl5", "mt-Co1"], ["Actb", "Gapdh", "Sox2", "Pax6", "Slc1a3"]),
            ("human", ["RPS3", "RPL5", "MT-CO1"], ["ACTB", "GAPDH", "SOX2", "PAX6", "SLC1A3"]),
        ):
            for include in (False, True):
                with self.subTest(organism=organism, include=include), tempfile.TemporaryDirectory() as temporary:
                    directory = Path(temporary) / "covariations_R0_F1"
                    directory.mkdir()
                    # Scrambled rows ensure the observed order comes from NiCo's ranking.
                    genes = np.array([
                        ordinary[2], excluded[1], ordinary[0], excluded[2],
                        ordinary[4], excluded[0], ordinary[1], ordinary[3],
                    ])
                    correlations = np.array([0.6, 0.9, 0.8, 0.85, 0.5, 0.95, 0.7, 0.55]).reshape(8, 1)
                    record = (
                        correlations, np.ones((3, 1)), genes, np.ones(8),
                        np.ones(8), correlations.copy(), 1.0,
                    )
                    factors_pickle = directory / "factors_info.p"
                    with factors_pickle.open("wb") as handle:
                        pickle.dump(({0: record}, {}, {}), handle)
                    artifact_before = factors_pickle.read_bytes()
                    state = SimpleNamespace(
                        covariation_dir=str(directory) + os.sep,
                        spatialcell_unique_clustername=["Cell A"],
                        spatialcell_unique_clusterid=[0],
                        no_of_pc=1,
                    )
                    state_before = dict(vars(state))
                    result = CovariationResult(
                        output_dir=Path(temporary), covariation_dir=directory,
                        radius=0, radius_tag="0", n_factors=1, modality="double",
                        factors_pickle=factors_pickle,
                        feature_matrix_npz=directory / "Principal_component_feature_matrix.npz",
                        regression_dir=directory / "Regression_outputs", nico_result=state,
                    )
                    config = CovariationReportConfig(
                        organism=organism, include_rps_rpl_mt_genes=include,
                        pathway_top_genes=4, pathway_databases=("ProbeDB",),
                        choose_celltypes=("Cell A",), choose_factors_id=(1,),
                        pathway_plot_as="barplot", saveas="png", dpi=72, show=False,
                    )
                    config_before = asdict(config)
                    expected_genes = [*excluded, ordinary[0]] if include else ordinary[:4]
                    original_pathway = self.scov.pathway_analysis
                    original_globals = dict(original_pathway.__globals__)
                    original_gseapy = self.scov.gseapy
                    original_enrichr = original_gseapy.enrichr

                    def enrichr(*, gene_list, organism, gene_sets, cutoff):
                        self.assertEqual(organism, config.organism)
                        self.assertEqual(gene_list, expected_genes)
                        self.assertEqual(gene_sets, "ProbeDB")
                        self.assertEqual(cutoff, 0.05)
                        self.assertIs(self.scov.pathway_analysis, original_pathway)
                        self.assertIs(self.scov.gseapy, original_gseapy)
                        self.assertIs(original_pathway.__globals__["gseapy"], original_gseapy)
                        return SimpleNamespace(res2d=pd.DataFrame({
                            "Adjusted P-value": [0.001], "Odds Ratio": [4.0], "Term": ["Test pathway"],
                        }))

                    # Only the enrichment request is mocked; NiCo filtering and plotting run unchanged.
                    with patch.object(original_gseapy, "enrichr", side_effect=enrichr) as request:
                        figures = run_pathway_enrichment(result, config=config)
                        request.assert_called_once()
                        self.assertEqual(request.call_args.kwargs["organism"], organism)
                        self.assertEqual(request.call_args.kwargs["gene_list"], expected_genes)
                    self.assertEqual(len(figures), 1)
                    self.assertEqual(figures[0].parent, directory / "Pathway_figures")
                    self.assertEqual(figures[0].suffix, ".png")
                    self.assertTrue(figures[0].is_file())
                    self.assertGreater(figures[0].stat().st_size, 0)
                    self.assertEqual(asdict(config), config_before)
                    self.assertEqual(vars(state), state_before)
                    self.assertEqual(factors_pickle.read_bytes(), artifact_before)
                    self.assertIs(self.scov.pathway_analysis, original_pathway)
                    self.assertEqual(original_pathway.__globals__, original_globals)
                    self.assertIs(original_gseapy.enrichr, original_enrichr)


if __name__ == "__main__":
    unittest.main()
