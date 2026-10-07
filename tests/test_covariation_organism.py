"""Report CLI settings, public organism contract, and NiCo call assembly."""

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


class ReportCliTests(unittest.TestCase):
    cases = (
        ("pathway", (), "run_pathway_enrichment"),
        ("reports", (), "generate_covariation_reports"),
        ("top-genes", ("--cell-type", "Cell A", "--factor-id", "2"), "extract_top_genes"),
        ("top-genes", ("--cell-type", "Cell A", "--all-factors"), "plot_top_genes_all_factors"),
        (
            "top-genes",
            ("--cell-type", "Cell A", "--factor-id", "2", "--pair-cell-type", "Cell B", "--pair-factor-id", "1"),
            "plot_top_genes_pair",
        ),
    )

    def setUp(self):
        self.runner = CliRunner()
        self.load = self.enterContext(patch.object(cli, "load_covariation_result", return_value=sentinel.result))
        self.calls = self.enterContext(patch.multiple(
            cli,
            extract_top_genes=DEFAULT,
            plot_top_genes_all_factors=DEFAULT,
            plot_top_genes_pair=DEFAULT,
            run_pathway_enrichment=DEFAULT,
            generate_covariation_reports=DEFAULT,
        ))
        table = pd.DataFrame({"gene": ["UnchangedGene"]})
        table.attrs["output_path"] = "genes.tsv"
        self.calls["extract_top_genes"].return_value = table
        for target in ("plot_top_genes_all_factors", "plot_top_genes_pair", "run_pathway_enrichment"):
            self.calls[target].return_value = (Path("figure.png"),)
        self.calls["generate_covariation_reports"].return_value = CovariationReportOutputs()

    def invoke_report(self, command, options, target):
        self.load.reset_mock()
        for call in self.calls.values():
            call.reset_mock()
        result = self.runner.invoke(cli.app, [command, "--output-dir", "existing-output", *options])
        self.assertEqual(result.exit_code, 0, result.output)
        self.load.assert_called_once_with(Path("existing-output"), radius="0", n_factors=3, load_state=True)
        for name, call in self.calls.items():
            if name != target:
                call.assert_not_called()
        call = self.calls[target]
        call.assert_called_once()
        self.assertIs(call.call_args.args[0], sentinel.result)
        config = call.call_args.kwargs["config"]
        self.assertIsInstance(config, CovariationReportConfig)
        if command == "pathway":
            self.assertIn("pathway_figure: figure.png", result.output)
        return config

    def test_inclusion_and_mouse_remain_defaults_for_every_report_path(self):
        for command, options, target in self.cases:
            with self.subTest(command=command, target=target):
                config = self.invoke_report(command, options, target)
                self.assertIs(config.include_rps_rpl_mt_genes, True)
                self.assertEqual(config.organism, "mouse")
                if command == "reports":
                    self.assertEqual(config.kinds, CovariationReportConfig().kinds)
                    self.assertNotIn("top-genes-all-factors", config.kinds)
                    self.assertNotIn("pathway", config.kinds)

    def test_explicit_inclusion_exclusion_and_species_reach_every_report_path(self):
        for command, options, target in self.cases:
            for organism in ("mouse", "human"):
                for flag, include in (("--include-rps-rpl-mt-genes", True), ("--exclude-rps-rpl-mt-genes", False)):
                    with self.subTest(command=command, target=target, organism=organism, include=include):
                        selected = (*options, "--organism", organism, flag)
                        if command == "reports":
                            selected = (*selected, "--kind", "top-genes-all-factors", "--kind", "pathway")
                        config = self.invoke_report(command, selected, target)
                        self.assertIs(config.include_rps_rpl_mt_genes, include)
                        self.assertEqual(config.organism, organism)
                        if command == "reports":
                            self.assertEqual(config.kinds, ("top-genes-all-factors", "pathway"))

    def test_all_top_gene_paths_preserve_other_settings_and_selections(self):
        for command, options, target in self.cases:
            if command != "top-genes":
                continue
            with self.subTest(target=target):
                config = self.invoke_report(command, (
                    *options, "--organism", "human", "--exclude-rps-rpl-mt-genes",
                    "--top-n", "7", "--plot-format", "png", "--negative", "--show",
                ), target)
                self.assertEqual(config, CovariationReportConfig(
                    organism="human", include_rps_rpl_mt_genes=False,
                    top_genes_per_factor=7, saveas="png", positively_correlated=False, show=True,
                    choose_celltypes=("Cell A",) if target == "plot_top_genes_all_factors" else (),
                ))
                kwargs = self.calls[target].call_args.kwargs
                if target == "extract_top_genes":
                    self.assertEqual(kwargs["cell_type"], "Cell A")
                    self.assertEqual(kwargs["factor_id"], 2)
                elif target == "plot_top_genes_pair":
                    self.assertEqual(kwargs["celltype_pair"], ("Cell A", "Cell B"))
                    self.assertEqual(kwargs["factor_ids"], (2, 1))

    def test_invalid_organisms_fail_before_loading_or_report_execution(self):
        for command, options, target in self.cases[:3]:
            for organism in ("Mouse", "Human", "zebrafish"):
                with self.subTest(command=command, organism=organism):
                    result = self.runner.invoke(cli.app, [
                        command, "--output-dir", "missing-output", *options, "--organism", organism,
                    ])
                    self.assertEqual(result.exit_code, 1, result.output)
                    error_prefix = "report" if command == "reports" else command
                    self.assertIn(
                        f"{error_prefix} error: organism must be 'mouse' or 'human' (lowercase); got {organism!r}.",
                        result.output,
                    )
                    self.load.assert_not_called()
                    for call in self.calls.values():
                        call.assert_not_called()

    def test_standalone_and_umbrella_help_expose_options_defaults_and_scope(self):
        from nico_wrapper.cli import app as umbrella_app

        for command in ("pathway", "top-genes", "reports"):
            for app, prefix in ((cli.app, ()), (umbrella_app, ("covariation",))):
                with self.subTest(command=command, prefix=prefix):
                    result = self.runner.invoke(
                        app, [*prefix, command, "--help"], color=False,
                        env={"COLUMNS": "240", "TERM": "dumb", "NO_COLOR": "1"},
                    )
                    self.assertEqual(result.exit_code, 0, result.output)
                    organism_line = next(line for line in result.output.splitlines() if "--organism" in line)
                    self.assertIn("mouse or human (lowercase)", organism_line)
                    self.assertIn("[default: mouse]", organism_line)
                    filter_line = next(line for line in result.output.splitlines() if "--include-rps-rpl-mt-genes" in line)
                    self.assertIn("--exclude-rps-rpl-mt-genes", filter_line)
                    self.assertIn("[default: include-rps-rpl-mt-genes]", filter_line)
                    self.assertIn("not model fitting", filter_line)
                    if command == "reports":
                        self.assertIn("only top-genes-all-factors and pathway reports", filter_line)


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
