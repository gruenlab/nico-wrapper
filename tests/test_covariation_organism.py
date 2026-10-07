"""Public organism contract and NiCo report-call assembly regressions."""

from dataclasses import FrozenInstanceError, asdict
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import DEFAULT, Mock, patch, sentinel

import pandas as pd
from typer.testing import CliRunner

from nico_wrapper.covariation import cli, reports
from nico_wrapper.covariation.config import CovariationReportConfig
from nico_wrapper.covariation.results import CovariationReportOutputs, CovariationResult


class OrganismConfigTests(unittest.TestCase):
    def test_default_and_supported_organisms_remain_lowercase_and_frozen(self):
        self.assertEqual(CovariationReportConfig().organism, "mouse")
        for organism in ("mouse", "human"):
            with self.subTest(organism=organism):
                config = CovariationReportConfig(organism=organism)
                self.assertEqual(config.organism, organism)
                self.assertEqual(asdict(config)["organism"], organism)
                with self.assertRaises(FrozenInstanceError):
                    config.organism = "human"

    def test_invalid_organisms_are_rejected_without_normalization(self):
        for organism in ("Mouse", "Human", "zebrafish"):
            with self.subTest(organism=organism):
                with self.assertRaises(ValueError) as raised:
                    CovariationReportConfig(organism=organism)
                self.assertEqual(
                    str(raised.exception),
                    f"organism must be 'mouse' or 'human' (lowercase); got {organism!r}.",
                )


class PathwayCliTests(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def test_default_and_explicit_organisms_are_forwarded_lowercase(self):
        for options, organism in (([], "mouse"), (["--organism", "mouse"], "mouse"), (["--organism", "human"], "human")):
            with self.subTest(options=options), patch.object(
                cli, "load_covariation_result", return_value=sentinel.result
            ) as load, patch.object(cli, "run_pathway_enrichment", return_value=(Path("figure.png"),)) as run:
                result = self.runner.invoke(cli.app, ["pathway", "--output-dir", "existing-output", *options])
                self.assertEqual(result.exit_code, 0, result.output)
                load.assert_called_once_with(Path("existing-output"), radius="0", n_factors=3, load_state=True)
                run.assert_called_once()
                self.assertIs(run.call_args.args[0], sentinel.result)
                config = run.call_args.kwargs["config"]
                self.assertIsInstance(config, CovariationReportConfig)
                self.assertEqual(config.organism, organism)
                self.assertIn("pathway_figure: figure.png", result.output)

    def test_invalid_input_fails_before_loading_or_enrichment(self):
        for organism in ("Mouse", "Human", "zebrafish"):
            with self.subTest(organism=organism), patch.object(
                cli, "load_covariation_result"
            ) as load, patch.object(cli, "run_pathway_enrichment") as run:
                result = self.runner.invoke(
                    cli.app, ["pathway", "--output-dir", "missing-output", "--organism", organism]
                )
                self.assertEqual(result.exit_code, 1, result.output)
                self.assertIn(
                    f"pathway error: organism must be 'mouse' or 'human' (lowercase); got {organism!r}.",
                    result.output,
                )
                load.assert_not_called()
                run.assert_not_called()

    def test_standalone_and_umbrella_help_show_lowercase_values_and_default(self):
        from nico_wrapper.cli import app as umbrella_app

        for app, command in ((cli.app, ["pathway"]), (umbrella_app, ["covariation", "pathway"])):
            with self.subTest(command=command):
                result = self.runner.invoke(
                    app, [*command, "--help"], color=False, env={"COLUMNS": "160", "TERM": "dumb", "NO_COLOR": "1"}
                )
                self.assertEqual(result.exit_code, 0, result.output)
                organism_line = next(line for line in result.output.splitlines() if "--organism" in line)
                self.assertIn("mouse or human (lowercase)", organism_line)
                self.assertIn("[default: mouse]", organism_line)


class ReportBoundaryTests(unittest.TestCase):
    def setUp(self):
        from nico import Covariations as scov

        self.scov = scov
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        self.state = SimpleNamespace(spatialcell_unique_clustername=["Cell A"], no_of_pc=3)
        self.result = CovariationResult(
            output_dir=directory,
            covariation_dir=directory,
            radius=0,
            radius_tag="0",
            n_factors=3,
            modality="double",
            factors_pickle=directory / "factors_info.p",
            feature_matrix_npz=directory / "Principal_component_feature_matrix.npz",
            regression_dir=directory / "Regression_outputs",
            nico_result=self.state,
        )

    def test_all_four_boundaries_preserve_selection_filtering_and_config(self):
        for organism in ("mouse", "human"):
            for include in (False, True):
                with self.subTest(organism=organism, include=include):
                    config = CovariationReportConfig(
                        organism=organism,
                        include_rps_rpl_mt_genes=include,
                        choose_celltypes=("Cell A",),
                        choose_factors_id=(2,),
                        top_genes_per_factor=7,
                        pathway_top_genes=4,
                        pathway_databases=("ProbeDB",),
                        pathway_plot_as="dotplot",
                        correlation_with_spearman=False,
                        positively_correlated=False,
                        saveas="png",
                        dpi=72,
                        transparent=True,
                        show=True,
                    )
                    before = asdict(config)
                    pathway = Mock()
                    with patch.multiple(
                        self.scov,
                        extract_and_plot_top_genes_from_chosen_factor_in_celltype=DEFAULT,
                        plot_top_genes_for_a_given_celltype_from_all_factors=DEFAULT,
                        plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors=DEFAULT,
                    ) as calls, patch.object(reports, "_adapt_nico_pathway_analysis", return_value=pathway) as adapt:
                        extract = calls["extract_and_plot_top_genes_from_chosen_factor_in_celltype"]
                        extract.return_value = pd.DataFrame({"gene": ["UnchangedGene"]})
                        table = reports.extract_top_genes(self.result, cell_type="Cell A", factor_id=2, config=config)
                        reports.plot_top_genes_all_factors(self.result, config=config)
                        reports.plot_top_genes_pair(
                            self.result, celltype_pair=("Cell A", "Cell B"), factor_ids=(2, 1), config=config
                        )
                        reports.run_pathway_enrichment(self.result, config=config)
                        common = dict(
                            organism=organism.capitalize(),
                            rps_rpl_mt_genes_included=include,
                            correlation_with_spearman=False,
                            saveas="png",
                            dpi=72,
                            transparent_mode=True,
                            showit=True,
                        )
                        extract.assert_called_once_with(
                            self.state, choose_celltype="Cell A", choose_factor_id=2,
                            top_NOG=7, positively_correlated=False, **common,
                        )
                        calls["plot_top_genes_for_a_given_celltype_from_all_factors"].assert_called_once_with(
                            self.state, choose_celltypes=["Cell A"], top_NOG=7, **common,
                        )
                        calls["plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors"].assert_called_once_with(
                            self.state, choose_interacting_celltype_pair=["Cell A", "Cell B"],
                            visualize_factors_id=[2, 1], top_NOG=7, **common,
                        )
                        adapt.assert_called_once_with(self.scov.pathway_analysis)
                        pathway.assert_called_once_with(
                            self.state, NOG_pathway=4, choose_celltypes=["Cell A"], choose_factors_id=[2],
                            database=["ProbeDB"], display_plot_as="dotplot", savefigure=True,
                            positively_correlated=False, **common,
                        )
                    self.assertEqual(asdict(config), before)
                    self.assertIs(self.result.nico_result, self.state)
                    self.assertEqual(table["gene"].tolist(), ["UnchangedGene"])
                    manifest = reports.write_report_manifest(
                        self.result, outputs=CovariationReportOutputs(), config=config, kinds=("pathway",)
                    )
                    self.assertEqual(json.loads(manifest.read_text())["config"]["organism"], organism)
                    self.assertEqual(asdict(config), before)

    def test_all_factor_extraction_and_dispatcher_use_corrected_boundaries(self):
        for organism in ("mouse", "human"):
            with self.subTest(organism=organism):
                config = CovariationReportConfig(
                    organism=organism, kinds=("top-genes-all-factors", "pathway"),
                    choose_celltypes=("Cell A",), choose_factors_id=(2,), include_rps_rpl_mt_genes=False,
                )
                pathway = Mock()
                with patch.object(
                    self.scov, "extract_and_plot_top_genes_from_chosen_factor_in_celltype",
                    return_value=pd.DataFrame({"gene": ["UnchangedGene"]}),
                ) as extract, patch.object(
                    self.scov, "plot_top_genes_for_a_given_celltype_from_all_factors"
                ) as plot, patch.object(reports, "_adapt_nico_pathway_analysis", return_value=pathway):
                    outputs = reports.generate_covariation_reports(self.result, config=config)
                    for call in (extract, plot, pathway):
                        call.assert_called_once()
                        self.assertEqual(call.call_args.kwargs["organism"], organism.capitalize())
                        self.assertIs(call.call_args.kwargs["rps_rpl_mt_genes_included"], False)
                    self.assertEqual(extract.call_args.kwargs["choose_factor_id"], 2)
                    self.assertEqual(len(outputs.top_gene_outputs), 1)
                    self.assertTrue(outputs.top_gene_outputs[0].is_file())
                    manifest = json.loads(outputs.report_manifest_json.read_text())
                    self.assertEqual(manifest["config"]["organism"], organism)
                    self.assertEqual(config.organism, organism)


if __name__ == "__main__":
    unittest.main()
