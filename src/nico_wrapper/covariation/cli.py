"""Typer app for NiCo niche covariation analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

from anndata import read_h5ad
import typer

from .config import CovariationConfig, CovariationReportConfig, DEFAULT_RIDGE_ALPHAS
from .io import planned_covariation_paths
from .pipeline import export_regression_table, load_covariation_result, run_covariation
from .reports import (
    extract_top_genes,
    generate_covariation_reports,
    plot_colocalized_factors,
    plot_factor_umap,
    plot_feature_matrix,
    plot_ligand_receptor_interactions,
    plot_top_genes_all_factors,
    plot_top_genes_pair,
    run_pathway_enrichment,
)
from .results import CovariationReportOutputs, CovariationResult
from .validation import ValidationError, resolve_niche_context, validate_covariation_inputs
from nico_wrapper.niche.io import normalize_radius

app = typer.Typer(
    name="nico-covariation",
    help="Run NiCo niche covariation analysis after niche interactions.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Run NiCo niche covariation analysis after niche interactions."""


@app.command("run")
def run(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir containing niche artifacts.")],
    ref_dir: Annotated[Optional[Path], typer.Option("--ref-dir", help="Reference input dir for double modality.")] = None,
    spatial_dir: Annotated[Optional[Path], typer.Option("--spatial-dir", help="Spatial/query input dir for double modality.")] = None,
    radius: Annotated[str, typer.Option("--radius", help="Radius matching nico-niche artifacts.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors per cell type.")] = 3,
    modality: Annotated[str, typer.Option("--modality", help="double or single.")] = "double",
    factorization: Annotated[str, typer.Option("--factorization", help="inmf or nmf-transfer.")] = "inmf",
    ref_label_key: Annotated[str, typer.Option("--ref-label-key", help="Reference .obs label column.")] = "cluster",
    annotated_h5ad_name: Annotated[str, typer.Option("--annotated-h5ad-name", help="Annotated h5ad name for single modality.")] = "nico_celltype_annotation.h5ad",
    ref_original_counts_name: Annotated[str, typer.Option("--ref-original-counts-name", help="Original reference h5ad filename.")] = "Original_counts.h5ad",
    ref_sct_name: Annotated[str, typer.Option("--ref-sct-name", help="Normalized reference h5ad filename.")] = "sct_singleCell.h5ad",
    spatial_sct_name: Annotated[str, typer.Option("--spatial-sct-name", help="Normalized spatial h5ad filename.")] = "sct_spatial.h5ad",
    ligand_receptor_db: Annotated[Optional[Path], typer.Option("--ligand-receptor-db", help="Ligand-receptor DB path; auto-detected when possible.")] = None,
    ridge_alpha: Annotated[Optional[list[float]], typer.Option("--ridge-alpha", help="RidgeCV alpha. May be repeated.")] = None,
    logistic_coef_cutoff: Annotated[float, typer.Option("--logistic-coef-cutoff", help="Niche logistic coefficient cutoff.")] = 0.0,
    ridge_coef_cutoff: Annotated[float, typer.Option("--ridge-coef-cutoff", help="Ridge coefficient cutoff used by exports/reports.")] = 0.0,
    expression_population_cutoff: Annotated[float, typer.Option("--expression-population-cutoff", help="Expression threshold for population fractions.")] = 0.0,
    seed: Annotated[int, typer.Option("--seed", help="Random seed passed to NiCo.")] = 541,
    persist_state: Annotated[bool, typer.Option("--persist-state/--no-persist-state", help="Write raw NiCo state pickle.")] = True,
    export_regression: Annotated[bool, typer.Option("--export-regression/--no-export-regression", help="Write regression coefficient TSV.")] = True,
    write_manifest: Annotated[bool, typer.Option("--write-manifest/--no-manifest", help="Write manifest JSON.")] = True,
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Replace existing outputs.")] = False,
) -> None:
    """Run core covariation analysis."""

    config = _make_config(
        radius=radius,
        n_factors=n_factors,
        modality=modality,
        factorization=factorization,
        ref_label_key=ref_label_key,
        annotated_h5ad_name=annotated_h5ad_name,
        ref_original_counts_name=ref_original_counts_name,
        ref_sct_name=ref_sct_name,
        spatial_sct_name=spatial_sct_name,
        ligand_receptor_db=ligand_receptor_db,
        ridge_alpha=ridge_alpha,
        logistic_coef_cutoff=logistic_coef_cutoff,
        ridge_coef_cutoff=ridge_coef_cutoff,
        expression_population_cutoff=expression_population_cutoff,
        seed=seed,
        persist_state=persist_state,
        export_regression=export_regression,
        write_manifest=write_manifest,
        overwrite=overwrite,
    )
    try:
        result = run_covariation(output_dir=output_dir, ref_dir=ref_dir, spatial_dir=spatial_dir, config=config)
    except ValidationError as exc:
        typer.secho(f"validation error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    _echo_result(result)


@app.command("validate")
def validate(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir containing niche artifacts.")],
    ref_dir: Annotated[Optional[Path], typer.Option("--ref-dir", help="Reference input dir for double modality.")] = None,
    spatial_dir: Annotated[Optional[Path], typer.Option("--spatial-dir", help="Spatial/query input dir for double modality.")] = None,
    radius: Annotated[str, typer.Option("--radius", help="Radius matching nico-niche artifacts.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors per cell type.")] = 3,
    modality: Annotated[str, typer.Option("--modality", help="double or single.")] = "double",
    factorization: Annotated[str, typer.Option("--factorization", help="inmf or nmf-transfer.")] = "inmf",
    ref_label_key: Annotated[str, typer.Option("--ref-label-key", help="Reference .obs label column.")] = "cluster",
    annotated_h5ad_name: Annotated[str, typer.Option("--annotated-h5ad-name", help="Annotated h5ad name for single modality.")] = "nico_celltype_annotation.h5ad",
    ref_original_counts_name: Annotated[str, typer.Option("--ref-original-counts-name", help="Original reference h5ad filename.")] = "Original_counts.h5ad",
    ref_sct_name: Annotated[str, typer.Option("--ref-sct-name", help="Normalized reference h5ad filename.")] = "sct_singleCell.h5ad",
    spatial_sct_name: Annotated[str, typer.Option("--spatial-sct-name", help="Normalized spatial h5ad filename.")] = "sct_spatial.h5ad",
    ligand_receptor_db: Annotated[Optional[Path], typer.Option("--ligand-receptor-db", help="Ligand-receptor DB path.")] = None,
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Allow existing planned outputs.")] = False,
) -> None:
    """Validate covariation inputs without running NiCo."""

    config = _make_config(
        radius=radius,
        n_factors=n_factors,
        modality=modality,
        factorization=factorization,
        ref_label_key=ref_label_key,
        annotated_h5ad_name=annotated_h5ad_name,
        ref_original_counts_name=ref_original_counts_name,
        ref_sct_name=ref_sct_name,
        spatial_sct_name=spatial_sct_name,
        ligand_receptor_db=ligand_receptor_db,
        overwrite=overwrite,
    )
    try:
        validate_covariation_inputs(output_dir=output_dir, ref_dir=ref_dir, spatial_dir=spatial_dir, config=config)
        niche = resolve_niche_context(niche_result=None, output_dir=output_dir, radius=config.radius)
    except ValidationError as exc:
        typer.secho(f"validation error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    paths = planned_covariation_paths(output_dir, radius=config.radius, n_factors=config.n_factors)
    typer.echo(f"validated: {output_dir}")
    typer.echo(f"radius: {paths.radius_tag}")
    typer.echo(f"n_factors: {config.n_factors}")
    typer.echo(f"modality: {config.modality}")
    typer.echo(f"cell_types: {len(niche.cell_type_names or {})}")
    if config.modality == "double" and ref_dir is not None and spatial_dir is not None:
        typer.echo(f"shared_genes: {_shared_raw_genes(ref_dir / config.ref_sct_name, spatial_dir / config.spatial_sct_name)}")
    typer.echo(f"planned_covariation_dir: {paths.covariation_dir}")


@app.command("artifacts")
def artifacts(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value to inspect.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    check: Annotated[bool, typer.Option("--check/--no-check", help="Exit nonzero when required artifacts are missing.")] = False,
) -> None:
    """Print expected covariation artifacts and existence status."""

    paths = planned_covariation_paths(output_dir, radius=radius, n_factors=n_factors)
    all_paths = {
        "factors_pickle": paths.factors_pickle,
        "feature_matrix_npz": paths.feature_matrix_npz,
        "regression_dir": paths.regression_dir,
        "state_pickle": paths.state_pickle,
        "regression_tsv": paths.regression_tsv,
        "manifest_json": paths.manifest_json,
    }
    missing_required = False
    for name, path in all_paths.items():
        ok = path.exists() and (path.is_dir() if name == "regression_dir" else path.is_file())
        if name in {"factors_pickle", "feature_matrix_npz", "regression_dir"} and not ok:
            missing_required = True
        typer.echo(f"{name}: {path} [{'ok' if ok else 'missing'}]")
    if check and missing_required:
        raise typer.Exit(1)


@app.command("export")
def export(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    kind: Annotated[str, typer.Option("--kind", help="Export kind; currently only regression.")] = "regression",
    output: Annotated[Optional[Path], typer.Option("--output", help="Output TSV path.")] = None,
    ridge_coef_cutoff: Annotated[Optional[float], typer.Option("--ridge-coef-cutoff", help="Override coefficient cutoff.")] = None,
) -> None:
    """Export stable tables from an existing result."""

    if kind != "regression":
        typer.secho("export error: only --kind regression is supported", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        path = export_regression_table(result, output_path=output, ridge_coef_cutoff=ridge_coef_cutoff)
    except (ValidationError, ValueError) as exc:
        typer.secho(f"export error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"regression_tsv: {path}")


@app.command("reports")
def reports(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    kind: Annotated[Optional[list[str]], typer.Option("--kind", help="Report kind. May be repeated.")] = None,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    dpi: Annotated[int, typer.Option("--dpi", help="Plot DPI.")] = 300,
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
    transparent: Annotated[bool, typer.Option("--transparent/--opaque", help="Transparent plot background.")] = False,
    pvalue_cutoff: Annotated[float, typer.Option("--pvalue-cutoff", help="P-value cutoff.")] = 0.05,
    cell_type: Annotated[Optional[list[str]], typer.Option("--cell-type", help="Restrict to selected cell types. May be repeated.")] = None,
    factor_id: Annotated[Optional[list[int]], typer.Option("--factor-id", help="Restrict to selected factor IDs. May be repeated.")] = None,
) -> None:
    """Generate optional reports from an existing covariation result."""

    config = CovariationReportConfig(
        kinds=tuple(kind) if kind else CovariationReportConfig().kinds,
        saveas=plot_format,
        dpi=dpi,
        show=show,
        transparent=transparent,
        pvalue_cutoff=pvalue_cutoff,
        choose_celltypes=tuple(cell_type or ()),
        choose_factors_id=tuple(factor_id or ()),
    )
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        outputs = generate_covariation_reports(result, config=config)
    except (ValidationError, ValueError) as exc:
        typer.secho(f"report error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    _echo_report_outputs(outputs)


@app.command("top-genes")
def top_genes(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    cell_type: Annotated[str, typer.Option("--cell-type", help="Cell type to inspect.")],
    factor_id: Annotated[Optional[int], typer.Option("--factor-id", help="1-based factor ID.")] = None,
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    top_n: Annotated[int, typer.Option("--top-n", help="Number of genes to export/plot.")] = 30,
    output: Annotated[Optional[Path], typer.Option("--output", help="Output TSV path.")] = None,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    positive: Annotated[bool, typer.Option("--positive/--negative", help="Use positive or negative factor correlations.")] = True,
    all_factors: Annotated[bool, typer.Option("--all-factors/--single-factor", help="Plot top genes across all factors.")] = False,
    pair_cell_type: Annotated[Optional[str], typer.Option("--pair-cell-type", help="Second cell type for paired top-gene plot.")] = None,
    pair_factor_id: Annotated[Optional[int], typer.Option("--pair-factor-id", help="Second factor ID for paired top-gene plot.")] = None,
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Extract top genes for one cell type/factor."""

    config = CovariationReportConfig(
        saveas=plot_format,
        show=show,
        top_genes_per_factor=top_n,
        positively_correlated=positive,
    )
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        if pair_cell_type is not None:
            if factor_id is None or pair_factor_id is None:
                raise ValueError("--pair-cell-type requires --factor-id and --pair-factor-id.")
            for path in plot_top_genes_pair(
                result,
                celltype_pair=(cell_type, pair_cell_type),
                factor_ids=(factor_id, pair_factor_id),
                config=config,
            ):
                typer.echo(f"top_genes_plot: {path}")
        elif all_factors:
            config = CovariationReportConfig(
                saveas=plot_format,
                show=show,
                top_genes_per_factor=top_n,
                positively_correlated=positive,
                choose_celltypes=(cell_type,),
            )
            for path in plot_top_genes_all_factors(result, config=config):
                typer.echo(f"top_genes_plot: {path}")
        else:
            if factor_id is None:
                raise ValueError("--factor-id is required unless --all-factors is used.")
            table = extract_top_genes(result, cell_type=cell_type, factor_id=factor_id, output_path=output, config=config)
            typer.echo(f"top_genes_tsv: {table.attrs.get('output_path')}")
    except (ValidationError, ValueError) as exc:
        typer.secho(f"top-genes error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc


@app.command("lr")
def lr(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    central_cell_type: Annotated[Optional[str], typer.Option("--central-cell-type", help="Central cell type for focused LR plots.")] = None,
    neighbor_cell_type: Annotated[Optional[str], typer.Option("--neighbor-cell-type", help="Neighbor cell type for focused LR plots.")] = None,
    central_factor_id: Annotated[Optional[int], typer.Option("--central-factor-id", help="Central factor for focused LR plots.")] = None,
    neighbor_factor_id: Annotated[Optional[int], typer.Option("--neighbor-factor-id", help="Neighbor factor for focused LR plots.")] = None,
    summary: Annotated[bool, typer.Option("--summary/--no-summary", help="Write LR summary workbook/text.")] = True,
    plots: Annotated[bool, typer.Option("--plots/--no-plots", help="Generate LR plots.")] = True,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    pvalue_cutoff: Annotated[float, typer.Option("--pvalue-cutoff", help="P-value cutoff.")] = 0.05,
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Generate ligand-receptor summaries and/or plots."""

    config = CovariationReportConfig(saveas=plot_format, pvalue_cutoff=pvalue_cutoff, show=show)
    pair = tuple(x for x in (central_cell_type, neighbor_cell_type) if x)
    factor_ids = tuple(x for x in (central_factor_id, neighbor_factor_id) if x is not None)
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        if summary:
            outputs = generate_covariation_reports(result, config=CovariationReportConfig(kinds=("lr-summary",), pvalue_cutoff=pvalue_cutoff))
            _echo_report_outputs(outputs)
        if plots:
            for path in plot_ligand_receptor_interactions(result, choose_interacting_celltype_pair=pair, choose_factors_id=factor_ids, config=config):
                typer.echo(f"ligand_receptor_plot: {path}")
    except (ValidationError, ValueError) as exc:
        typer.secho(f"lr error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc


@app.command("pathway")
def pathway(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    cell_type: Annotated[Optional[list[str]], typer.Option("--cell-type", help="Restrict to selected cell types. May be repeated.")] = None,
    factor_id: Annotated[Optional[list[int]], typer.Option("--factor-id", help="Restrict to selected factors. May be repeated.")] = None,
    top_genes: Annotated[int, typer.Option("--top-genes", help="Top genes per factor for enrichment.")] = 50,
    database: Annotated[Optional[list[str]], typer.Option("--database", help="Enrichr database. May be repeated.")] = None,
    organism: Annotated[str, typer.Option("--organism", help="Mouse or Human.")] = "Mouse",
    plot_as: Annotated[str, typer.Option("--plot-as", help="barplot or dotplot.")] = "barplot",
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Run optional pathway enrichment. May require network access to Enrichr."""

    config = CovariationReportConfig(
        saveas=plot_format,
        show=show,
        choose_celltypes=tuple(cell_type or ()),
        choose_factors_id=tuple(factor_id or ()),
        pathway_top_genes=top_genes,
        pathway_databases=tuple(database) if database else CovariationReportConfig().pathway_databases,
        organism=organism,  # type: ignore[arg-type]
        pathway_plot_as=plot_as,  # type: ignore[arg-type]
    )
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        for path in run_pathway_enrichment(result, config=config):
            typer.echo(f"pathway_figure: {path}")
    except (ValidationError, ValueError) as exc:
        typer.secho(f"pathway error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc


@app.command("umap")
def umap(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    cell_type: Annotated[list[str], typer.Option("--cell-type", help="Cell type. Repeat once or twice.")],
    factor_id: Annotated[list[int], typer.Option("--factor-id", help="1-based factor ID. Repeat once or twice.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    modality: Annotated[str, typer.Option("--modality", help="sc or spatial.")] = "spatial",
    umap_key: Annotated[str, typer.Option("--umap-key", help="AnnData obsm UMAP key.")] = "X_umap",
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Plot factor loadings on UMAP coordinates."""

    config = CovariationReportConfig(saveas=plot_format, show=show)
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        path = plot_factor_umap(
            result,
            modality=modality,  # type: ignore[arg-type]
            celltype_pair=tuple(cell_type),
            factor_ids=tuple(factor_id),
            umap_key=umap_key,
            config=config,
        )
    except (ValidationError, ValueError, KeyError) as exc:
        typer.secho(f"umap error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"umap_plot: {path}")


@app.command("feature-matrix")
def feature_matrix(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Plot the weighted factor-neighborhood feature matrix."""

    config = CovariationReportConfig(saveas=plot_format, show=show)
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        path = plot_feature_matrix(result, config=config)
    except (ValidationError, ValueError) as exc:
        typer.secho(f"feature-matrix error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"feature_matrix_plot: {path}")


@app.command("colocalize")
def colocalize(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Base NiCo output dir.")],
    central_cell_type: Annotated[str, typer.Option("--central-cell-type", help="Central cell type.")],
    neighbor_cell_type: Annotated[str, typer.Option("--neighbor-cell-type", help="Neighbor cell type.")],
    central_factor_id: Annotated[int, typer.Option("--central-factor-id", help="Central factor ID.")],
    neighbor_factor_id: Annotated[int, typer.Option("--neighbor-factor-id", help="Neighbor factor ID.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value.")] = "0",
    n_factors: Annotated[int, typer.Option("--n-factors", help="Number of latent factors.")] = 3,
    bar: Annotated[bool, typer.Option("--bar/--no-bar", help="Also create bar plot.")] = True,
    violin: Annotated[bool, typer.Option("--violin/--no-violin", help="Also create violin plot.")] = False,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    show: Annotated[bool, typer.Option("--show/--no-show", help="Keep figures open.")] = False,
) -> None:
    """Plot colocalized central/neighbor factor loadings."""

    config = CovariationReportConfig(saveas=plot_format, show=show)
    try:
        result = load_covariation_result(output_dir, radius=radius, n_factors=n_factors, load_state=True)
        paths = plot_colocalized_factors(
            result,
            central_cell_type=central_cell_type,
            neighbor_cell_type=neighbor_cell_type,
            central_factor_id=central_factor_id,
            neighbor_factor_id=neighbor_factor_id,
            include_bar=bar,
            include_violin=violin,
            config=config,
        )
    except (ValidationError, ValueError) as exc:
        typer.secho(f"colocalize error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    for path in paths:
        typer.echo(f"colocalization_plot: {path}")


def _make_config(
    *,
    radius: str,
    n_factors: int,
    modality: str,
    factorization: str,
    ref_label_key: str,
    annotated_h5ad_name: str = "nico_celltype_annotation.h5ad",
    ref_original_counts_name: str = "Original_counts.h5ad",
    ref_sct_name: str = "sct_singleCell.h5ad",
    spatial_sct_name: str = "sct_spatial.h5ad",
    ligand_receptor_db: Path | None = None,
    ridge_alpha: list[float] | None = None,
    logistic_coef_cutoff: float = 0.0,
    ridge_coef_cutoff: float = 0.0,
    expression_population_cutoff: float = 0.0,
    seed: int = 541,
    persist_state: bool = True,
    export_regression: bool = True,
    write_manifest: bool = True,
    overwrite: bool = False,
) -> CovariationConfig:
    radius_value, _ = normalize_radius(radius)
    return CovariationConfig(
        radius=radius_value,
        n_factors=n_factors,
        modality=modality,  # type: ignore[arg-type]
        factorization=factorization,  # type: ignore[arg-type]
        ref_label_key=ref_label_key,
        annotated_h5ad_name=annotated_h5ad_name,
        ref_original_counts_name=ref_original_counts_name,
        ref_sct_name=ref_sct_name,
        spatial_sct_name=spatial_sct_name,
        ligand_receptor_db=ligand_receptor_db,
        ridge_alphas=tuple(ridge_alpha) if ridge_alpha else DEFAULT_RIDGE_ALPHAS,
        logistic_coef_cutoff=logistic_coef_cutoff,
        ridge_coef_cutoff=ridge_coef_cutoff,
        expression_population_cutoff=expression_population_cutoff,
        seed=seed,
        persist_state=persist_state,
        export_regression_table=export_regression,
        write_manifest=write_manifest,
        overwrite=overwrite,
    )


def _shared_raw_genes(ref_sct: Path, spatial_sct: Path) -> int:
    ref = read_h5ad(ref_sct, backed="r")
    spatial = read_h5ad(spatial_sct, backed="r")
    try:
        if ref.raw is None or spatial.raw is None:
            return 0
        return len(set(ref.raw.var_names).intersection(map(str, spatial.raw.var_names)))
    finally:
        ref.file.close()
        spatial.file.close()


def _echo_result(result: CovariationResult) -> None:
    typer.echo(f"covariation_dir: {result.covariation_dir}")
    typer.echo(f"factors_pickle: {result.factors_pickle}")
    typer.echo(f"feature_matrix_npz: {result.feature_matrix_npz}")
    typer.echo(f"regression_dir: {result.regression_dir}")
    if result.state_pickle is not None:
        typer.echo(f"state_pickle: {result.state_pickle}")
    if result.regression_tsv is not None:
        typer.echo(f"regression_tsv: {result.regression_tsv}")
    if result.manifest_json is not None:
        typer.echo(f"manifest_json: {result.manifest_json}")


def _echo_report_outputs(outputs: CovariationReportOutputs) -> None:
    for path in outputs.regression_circleplots:
        typer.echo(f"regression_circleplot: {path}")
    for path in outputs.regression_heatmaps:
        typer.echo(f"regression_heatmap: {path}")
    if outputs.gene_correlation_excel:
        typer.echo(f"gene_correlation_excel: {outputs.gene_correlation_excel}")
    for path in outputs.factor_gene_heatmaps:
        typer.echo(f"factor_gene_heatmap: {path}")
    if outputs.ligand_receptor_summary_xlsx:
        typer.echo(f"ligand_receptor_summary_xlsx: {outputs.ligand_receptor_summary_xlsx}")
    if outputs.regression_summary_txt:
        typer.echo(f"regression_summary_txt: {outputs.regression_summary_txt}")
    for path in outputs.ligand_receptor_plots:
        typer.echo(f"ligand_receptor_plot: {path}")
    for path in outputs.pathway_figures:
        typer.echo(f"pathway_figure: {path}")
    for path in outputs.top_gene_outputs:
        typer.echo(f"top_genes_tsv: {path}")
    for path in outputs.top_gene_plots:
        typer.echo(f"top_genes_plot: {path}")
    for path in outputs.umap_plots:
        typer.echo(f"umap_plot: {path}")
    if outputs.feature_matrix_plot:
        typer.echo(f"feature_matrix_plot: {outputs.feature_matrix_plot}")
    if outputs.report_manifest_json:
        typer.echo(f"report_manifest_json: {outputs.report_manifest_json}")
    for path in outputs.other_outputs:
        typer.echo(f"report: {path}")


if __name__ == "__main__":
    app()
