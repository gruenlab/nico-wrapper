"""Optional report wrappers for NiCo covariation outputs."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass, replace
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any, Literal

import pandas as pd

from .config import CovariationReportConfig
from .io import read_covariation_state
from .results import CovariationReportOutputs, CovariationResult


def require_nico_covariation_state(result: CovariationResult):
    """Return raw NiCo covariation state or fail with a clear message."""

    if result.nico_result is not None:
        return result.nico_result
    if result.state_pickle is not None and result.state_pickle.exists():
        return read_covariation_state(result.state_pickle)
    raise ValueError(
        "This covariation report requires NiCo covariation state. "
        "Load with load_state=True or re-run with persist_state=True."
    )


def generate_covariation_reports(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> CovariationReportOutputs:
    """Generate selected optional covariation reports/plots.

    Kinds needing specific cell-type/factor choices are exposed as focused
    functions/CLI commands rather than through this broad dispatcher.
    """

    kinds = set(config.kinds)
    if "all" in kinds:
        kinds = {
            "regression-circleplots",
            "regression-heatmaps",
            "pvalue-sizebar",
            "gene-correlation-excel",
            "lr-summary",
            "lr-plots",
            "factor-gene-heatmaps",
            "top-genes-all-factors",
            "feature-matrix",
        }

    supported = {
        "regression-circleplots",
        "regression-heatmaps",
        "pvalue-sizebar",
        "gene-correlation-excel",
        "lr-summary",
        "lr-plots",
        "factor-gene-heatmaps",
        "top-genes-all-factors",
        "pathway",
        "feature-matrix",
        "sc-umap",
        "spatial-umap",
    }
    unsupported = kinds.difference(supported)
    if unsupported:
        raise ValueError("Unsupported report kind(s): " + ", ".join(sorted(unsupported)))

    circleplots: tuple[Path, ...] = ()
    heatmaps: tuple[Path, ...] = ()
    other: tuple[Path, ...] = ()
    gene_excel: Path | None = None
    lr_xlsx: Path | None = None
    lr_txt: Path | None = None
    factor_heatmaps: tuple[Path, ...] = ()
    lr_plots: tuple[Path, ...] = ()
    pathway_figures: tuple[Path, ...] = ()
    top_gene_outputs: tuple[Path, ...] = ()
    top_gene_plots: tuple[Path, ...] = ()
    feature_matrix_plot: Path | None = None
    umap_plots: tuple[Path, ...] = ()

    if "gene-correlation-excel" in kinds:
        gene_excel = export_gene_correlations(result)
    if "regression-circleplots" in kinds:
        circleplots = plot_regression_covariations(result, kind="circle", config=config)
    if "regression-heatmaps" in kinds:
        heatmaps = plot_regression_covariations(result, kind="heatmap", config=config)
    if "pvalue-sizebar" in kinds:
        other = (*other, *print_pvalue_sizebar(result, config=config))
    if "lr-summary" in kinds:
        lr_xlsx, lr_txt = export_ligand_receptor_summary(result, config=config)
    if "lr-plots" in kinds:
        lr_plots = plot_ligand_receptor_interactions(result, config=config)
    if "factor-gene-heatmaps" in kinds:
        factor_heatmaps = plot_factor_gene_heatmaps(result, config=config)
    if "top-genes-all-factors" in kinds:
        top_gene_outputs = extract_top_genes_all_factors(result, config=config)
        top_gene_plots = plot_top_genes_all_factors(result, config=config)
    if "pathway" in kinds:
        pathway_figures = run_pathway_enrichment(result, config=config)
    if "feature-matrix" in kinds:
        feature_matrix_plot = plot_feature_matrix(result, config=config)
    for kind, modality in (("sc-umap", "sc"), ("spatial-umap", "spatial")):
        if kind in kinds:
            umap_plots = (
                *umap_plots,
                plot_factor_umap(
                    result,
                    modality=modality,  # type: ignore[arg-type]
                    celltype_pair=config.choose_celltypes,
                    factor_ids=config.choose_factors_id,
                    config=config,
                ),
            )

    outputs = CovariationReportOutputs(
        regression_circleplots=circleplots,
        regression_heatmaps=heatmaps,
        gene_correlation_excel=gene_excel,
        factor_gene_heatmaps=factor_heatmaps,
        ligand_receptor_summary_xlsx=lr_xlsx,
        regression_summary_txt=lr_txt,
        ligand_receptor_plots=lr_plots,
        pathway_figures=pathway_figures,
        top_gene_outputs=top_gene_outputs,
        top_gene_plots=top_gene_plots,
        umap_plots=umap_plots,
        feature_matrix_plot=feature_matrix_plot,
        other_outputs=other,
    )
    manifest = write_report_manifest(result, outputs=outputs, config=config, kinds=tuple(sorted(kinds)))
    return replace(outputs, report_manifest_json=manifest)


def export_gene_correlations(result: CovariationResult) -> Path:
    """Write NiCo's gene-correlation Excel workbook."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    scov.make_excel_sheet_for_gene_correlation(state)
    return result.covariation_dir / "gene_correlation.xlsx"


def plot_regression_covariations(
    result: CovariationResult,
    *,
    kind: Literal["circle", "heatmap"] = "circle",
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Generate regression covariation plots and return created files."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    before = _snapshot(result.regression_dir, config.saveas)
    common = {
        "choose_celltypes": list(config.choose_celltypes),
        "saveas": config.saveas,
        "transparent_mode": config.transparent,
        "showit": config.show,
        "dpi": config.dpi,
    }
    if kind == "circle":
        scov.plot_significant_regression_covariations_as_circleplot(
            state,
            pvalue_cutoff=config.pvalue_cutoff,
            **common,
        )
    else:
        scov.plot_significant_regression_covariations_as_heatmap(state, **common)
    return _created_or_existing(result.regression_dir, config.saveas, before)


def print_pvalue_sizebar(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Generate NiCo's p-value sizebar legend."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    before = _snapshot(result.regression_dir, config.saveas)
    scov.print_pvalue_sizebar(
        state,
        saveas=config.saveas,
        showit=config.show,
        transparent_mode=config.transparent,
        dpi=config.dpi,
    )
    path = result.regression_dir / f"pvalue_cirlce_sizebar.{config.saveas}"
    return _created_or_known(path, result.regression_dir, config.saveas, before)


def export_ligand_receptor_summary(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, Path]:
    """Write NiCo's ligand-receptor workbook and regression summary text."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    scov.save_LR_interactions_in_excelsheet_and_regression_summary_in_textfile_for_interacting_cell_types(
        state,
        pvalueCutoff=config.pvalue_cutoff,
        correlation_with_spearman=config.correlation_with_spearman,
        Ligand_Factor_thres=config.ligand_factor_threshold,
        Receptor_Factor_thres=config.receptor_factor_threshold,
        Ligand_proportion_of_cells_expressed_thres=config.ligand_population_threshold,
        Receptor_proportion_of_cells_expressed_thres=config.receptor_population_threshold,
        number_of_top_genes_to_print=config.number_of_top_genes_to_print,
    )
    return (
        result.covariation_dir / "Lig_and_Rec_enrichment_in_interacting_celltypes.xlsx",
        result.covariation_dir / "Regression_summary.txt",
    )


def plot_ligand_receptor_interactions(
    result: CovariationResult,
    *,
    choose_interacting_celltype_pair: tuple[str, ...] = (),
    choose_factors_id: tuple[int, ...] = (),
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Generate ligand-receptor interaction plots."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    directories = (
        result.covariation_dir / "Plot_ligand_receptor_in_niche",
        result.covariation_dir / "Plot_ligand_receptor_in_niche_cc_vs_nc",
        result.covariation_dir / "Plot_ligand_receptor_in_niche_nc_vs_cc",
    )
    before = set().union(*(_snapshot(directory, config.saveas) for directory in directories))
    scov.find_LR_interactions_in_interacting_cell_types(
        state,
        choose_interacting_celltype_pair=list(choose_interacting_celltype_pair),
        choose_factors_id=list(choose_factors_id),
        pvalueCutoff=config.pvalue_cutoff,
        dpi=config.dpi,
        correlation_with_spearman=config.correlation_with_spearman,
        Ligand_Factor_thres=config.ligand_factor_threshold,
        Receptor_Factor_thres=config.receptor_factor_threshold,
        Ligand_proportion_of_cells_expressed_thres=config.ligand_population_threshold,
        Receptor_proportion_of_cells_expressed_thres=config.receptor_population_threshold,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
    )
    after = set().union(*(_snapshot(directory, config.saveas) for directory in directories))
    created = after.difference(before)
    return tuple(sorted(created or after))


def plot_factor_gene_heatmaps(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Generate cosine/Spearman factor-gene heatmaps."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    output_dir = result.covariation_dir / "NMF_output"
    before = _snapshot(output_dir, config.saveas)
    scov.plot_cosine_and_spearman_correlation_to_factors(
        state,
        choose_celltypes=list(config.choose_celltypes),
        NOG_Fa=config.top_genes_per_factor,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
        dpi=config.dpi,
    )
    return _created_or_existing(output_dir, config.saveas, before)


def extract_top_genes(
    result: CovariationResult,
    *,
    cell_type: str,
    factor_id: int,
    output_path: str | Path | None = None,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> pd.DataFrame:
    """Extract and plot top genes for one cell type/factor.

    The returned DataFrame has ``attrs['output_path']`` set to the wrapper TSV
    that was written.
    """

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    table = scov.extract_and_plot_top_genes_from_chosen_factor_in_celltype(
        state,
        choose_celltype=cell_type,
        choose_factor_id=factor_id,
        top_NOG=config.top_genes_per_factor,
        rps_rpl_mt_genes_included=config.include_rps_rpl_mt_genes,
        organism=config.organism,
        correlation_with_spearman=config.correlation_with_spearman,
        positively_correlated=config.positively_correlated,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
        dpi=config.dpi,
    )
    if not isinstance(table, pd.DataFrame):
        raise ValueError(f"No top-gene table was returned for cell type {cell_type!r}, factor {factor_id}.")
    path = Path(output_path) if output_path is not None else _top_genes_path(result, cell_type, factor_id, config)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)
    table.attrs["output_path"] = str(path)
    return table


def extract_top_genes_all_factors(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Extract top-gene TSVs for all selected cell types/factors."""

    state = require_nico_covariation_state(result)
    cell_types = tuple(config.choose_celltypes) or tuple(map(str, state.spatialcell_unique_clustername))
    factors = tuple(config.choose_factors_id) or tuple(range(1, int(state.no_of_pc) + 1))
    paths: list[Path] = []
    for cell_type in cell_types:
        for factor_id in factors:
            table = extract_top_genes(result, cell_type=cell_type, factor_id=factor_id, config=config)
            output_path = table.attrs.get("output_path")
            if output_path is not None:
                paths.append(Path(output_path))
    return tuple(paths)


def plot_top_genes_all_factors(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Plot top genes from all factors for selected cell types."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    output_dir = result.covariation_dir / "dotplots"
    before = _snapshot(output_dir, config.saveas)
    scov.plot_top_genes_for_a_given_celltype_from_all_factors(
        state,
        choose_celltypes=list(config.choose_celltypes),
        top_NOG=config.top_genes_per_factor,
        rps_rpl_mt_genes_included=config.include_rps_rpl_mt_genes,
        organism=config.organism,
        correlation_with_spearman=config.correlation_with_spearman,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
        dpi=config.dpi,
    )
    return _created_or_existing(output_dir, config.saveas, before)


def plot_top_genes_pair(
    result: CovariationResult,
    *,
    celltype_pair: tuple[str, str],
    factor_ids: tuple[int, int],
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Plot top genes for a pair of cell types and selected factors."""

    from nico import Covariations as scov

    output_dir = result.covariation_dir / "dotplots"
    before = _snapshot(output_dir, config.saveas)
    state = require_nico_covariation_state(result)
    scov.plot_top_genes_for_pair_of_celltypes_from_two_chosen_factors(
        state,
        choose_interacting_celltype_pair=list(celltype_pair),
        visualize_factors_id=list(factor_ids),
        top_NOG=config.top_genes_per_factor,
        dpi=config.dpi,
        organism=config.organism,
        rps_rpl_mt_genes_included=config.include_rps_rpl_mt_genes,
        correlation_with_spearman=config.correlation_with_spearman,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
    )
    return _created_or_existing(output_dir, config.saveas, before)


def plot_feature_matrix(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> Path:
    """Plot NiCo's weighted factor-neighborhood feature matrix."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    scov.plot_feature_matrices(
        state,
        showit=config.show,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        dpi=config.dpi,
    )
    return result.covariation_dir / f"Feature_matrix_PC.{config.saveas}"


def run_pathway_enrichment(
    result: CovariationResult,
    *,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Run NiCo pathway enrichment and return generated pathway figures."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    output_dir = result.covariation_dir / "Pathway_figures"
    before = _snapshot(output_dir, config.saveas)
    scov.pathway_analysis(
        state,
        NOG_pathway=config.pathway_top_genes,
        choose_factors_id=list(config.choose_factors_id),
        correlation_with_spearman=config.correlation_with_spearman,
        saveas=config.saveas,
        savefigure=True,
        positively_correlated=config.positively_correlated,
        rps_rpl_mt_genes_included=config.include_rps_rpl_mt_genes,
        choose_celltypes=list(config.choose_celltypes),
        organism=config.organism,
        database=list(config.pathway_databases),
        display_plot_as=config.pathway_plot_as,
        showit=config.show,
        transparent_mode=config.transparent,
        dpi=config.dpi,
    )
    return _created_or_existing(output_dir, config.saveas, before)


def plot_factor_umap(
    result: CovariationResult,
    *,
    modality: Literal["sc", "spatial"],
    celltype_pair: tuple[str, ...],
    factor_ids: tuple[int, ...],
    umap_key: str = "X_umap",
    config: CovariationReportConfig = CovariationReportConfig(),
) -> Path:
    """Plot factor loadings on scRNA-seq or spatial UMAP coordinates."""

    from nico import Covariations as scov

    if not celltype_pair or not factor_ids:
        raise ValueError("celltype_pair and factor_ids must be non-empty.")
    if len(celltype_pair) not in {1, 2} or len(factor_ids) not in {1, 2}:
        raise ValueError("NiCo UMAP plotting supports one or two cell types/factors.")
    state = require_nico_covariation_state(result)
    if modality == "sc":
        scov.visualize_factors_in_scRNAseq_umap(
            state,
            choose_interacting_celltype_pair=list(celltype_pair),
            visualize_factors_id=list(factor_ids),
            umap_tag=umap_key,
            saveas=config.saveas,
            transparent_mode=config.transparent,
            showit=config.show,
            dpi=config.dpi,
        )
        return result.covariation_dir / f"scRNAseq_factors_in_umap.{config.saveas}"
    scov.visualize_factors_in_spatial_umap(
        state,
        choose_interacting_celltype_pair=list(celltype_pair),
        visualize_factors_id=list(factor_ids),
        umap_tag=umap_key,
        saveas=config.saveas,
        transparent_mode=config.transparent,
        showit=config.show,
        dpi=config.dpi,
    )
    return result.covariation_dir / f"spatial_factors_in_umap.{config.saveas}"


def plot_colocalized_factors(
    result: CovariationResult,
    *,
    central_cell_type: str,
    neighbor_cell_type: str,
    central_factor_id: int,
    neighbor_factor_id: int,
    include_bar: bool = True,
    include_violin: bool = False,
    config: CovariationReportConfig = CovariationReportConfig(),
) -> tuple[Path, ...]:
    """Plot colocalized central/neighbor factor loadings."""

    from nico import Covariations as scov

    state = require_nico_covariation_state(result)
    output_dir = result.covariation_dir / "colocalization"
    before = _snapshot(output_dir, config.saveas)
    output_cc, output_nc, cc_not_colocalized, nc_not_colocalized = scov.visualization_of_colocalized_celltype_factors_as_scatterplot(
        state,
        CC_name=central_cell_type,
        NC_name=neighbor_cell_type,
        CC_factor_id=central_factor_id,
        NC_factor_id=neighbor_factor_id,
        saveas=config.saveas,
        showit=config.show,
        transparent_mode=config.transparent,
        dpi=config.dpi,
    )
    if include_bar:
        scov.visualization_of_colocalized_celltype_factors_as_bar_violin_plot(
            state,
            CC_name=central_cell_type,
            NC_name=neighbor_cell_type,
            CC_factor_id=central_factor_id,
            NC_factor_id=neighbor_factor_id,
            CC_unique_colocalized_loadings=output_cc[1],
            NC_unique_colocalized_loadings=output_nc[1],
            CC_not_colocalized_loadings=cc_not_colocalized,
            NC_not_colocalized_loadings=nc_not_colocalized,
            visualize_as="BarPlot",
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
    if include_violin:
        scov.visualization_of_colocalized_celltype_factors_as_bar_violin_plot(
            state,
            CC_name=central_cell_type,
            NC_name=neighbor_cell_type,
            CC_factor_id=central_factor_id,
            NC_factor_id=neighbor_factor_id,
            CC_unique_colocalized_loadings=output_cc[1],
            NC_unique_colocalized_loadings=output_nc[1],
            CC_not_colocalized_loadings=cc_not_colocalized,
            NC_not_colocalized_loadings=nc_not_colocalized,
            visualize_as="ViolinPlot",
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
    return _created_or_existing(output_dir, config.saveas, before)


def write_report_manifest(
    result: CovariationResult,
    *,
    outputs: CovariationReportOutputs,
    config: CovariationReportConfig,
    kinds: tuple[str, ...],
    output_path: str | Path | None = None,
) -> Path:
    """Write a manifest for optional report outputs."""

    path = Path(output_path) if output_path is not None else result.covariation_dir / "report_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "covariation_dir": str(result.covariation_dir),
        "radius": result.radius,
        "n_factors": result.n_factors,
        "kinds": list(kinds),
        "outputs": _jsonify(asdict(outputs)),
        "config": _jsonify(asdict(config) if is_dataclass(config) else config),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _top_genes_path(
    result: CovariationResult,
    cell_type: str,
    factor_id: int,
    config: CovariationReportConfig,
) -> Path:
    direction = "pos" if config.positively_correlated else "neg"
    return result.covariation_dir / f"top_genes_{_safe_stem(cell_type)}_Fa{factor_id}_{direction}.tsv"


def _safe_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "cell_type"


def _jsonify(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_jsonify(item) for item in value]
    if isinstance(value, list):
        return [_jsonify(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonify(item) for key, item in value.items()}
    return value


def _snapshot(directory: Path, suffix: str) -> set[Path]:
    if not directory.exists():
        return set()
    return set(directory.rglob(f"*.{suffix}"))


def _created_or_existing(directory: Path, suffix: str, before: set[Path]) -> tuple[Path, ...]:
    if not directory.exists():
        return ()
    after = set(directory.rglob(f"*.{suffix}"))
    created = after.difference(before)
    return tuple(sorted(created or after))


def _created_or_known(path: Path, directory: Path, suffix: str, before: set[Path]) -> tuple[Path, ...]:
    created = _created_or_existing(directory, suffix, before)
    if created:
        return created
    return (path,) if path.exists() else ()
