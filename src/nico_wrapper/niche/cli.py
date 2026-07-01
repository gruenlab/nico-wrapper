"""Typer command line app for NiCo niche interaction analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

from anndata import read_h5ad
import typer

from .config import (
    InteractionModelConfig,
    NeighborhoodConfig,
    NicheInteractionConfig,
    NichePlotConfig,
    ProximityConfig,
)
from .io import normalize_radius, planned_artifact_paths
from .pipeline import export_interaction_table, load_niche_result, run_niche_interactions
from .plotting import plot_niche_result
from .proximity import run_proximity_analysis
from .results import NicheInteractionResult
from .validation import ValidationError, validate_covariation_artifacts, validate_niche_inputs

app = typer.Typer(
    name="nico-niche",
    help="Run NiCo niche interaction analysis after label transfer.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Run NiCo niche interaction analysis after label transfer."""


@app.command("run")
def run(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing transferred labels and receiving niche outputs.")],
    anndata_filename: Annotated[
        str,
        typer.Option("--anndata-filename", help="Annotated spatial AnnData file name under output-dir."),
    ] = "nico_celltype_annotation.h5ad",
    label_key: Annotated[str, typer.Option("--label-key", help=".obs column containing transferred labels.")] = "nico_ct",
    spatial_key: Annotated[str, typer.Option("--spatial-key", help=".obsm key containing spatial coordinates.")] = "spatial",
    radius: Annotated[str, typer.Option("--radius", help="0 for Delaunay; positive value for fixed-radius neighbors.")] = "0",
    epsilon_threshold: Annotated[
        float,
        typer.Option("--epsilon-threshold", help="Maximum Delaunay edge length when radius is 0."),
    ] = 100.0,
    exclude_cell_type: Annotated[
        Optional[list[str]],
        typer.Option("--exclude-cell-type", help="Additional cell type to exclude. May be repeated."),
    ] = None,
    k_fold: Annotated[int, typer.Option("--k-fold", help="Stratified cross-validation folds.")] = 5,
    n_repeats: Annotated[int, typer.Option("--n-repeats", help="Repeated CV rounds after C selection.")] = 1,
    seed: Annotated[int, typer.Option("--seed", help="Random seed passed to NiCo.")] = 36851234,
    n_jobs: Annotated[int, typer.Option("--n-jobs", help="sklearn parallel jobs; -1 uses all CPUs.")] = -1,
    c_value: Annotated[
        Optional[list[float]],
        typer.Option("--c-value", help="Candidate sklearn C value. May be repeated."),
    ] = None,
    make_plots: Annotated[bool, typer.Option("--make-plots/--no-plots", help="Generate selected plots after analysis.")] = False,
    plot_kind: Annotated[
        Optional[list[str]],
        typer.Option("--plot-kind", help="Plot kind to generate. May be repeated."),
    ] = None,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format, e.g. pdf or png.")] = "pdf",
    interaction_cutoff: Annotated[
        float,
        typer.Option("--interaction-cutoff", help="Positive normalized coefficient cutoff for plots/tables."),
    ] = 0.1,
    proximity: Annotated[bool, typer.Option("--proximity/--no-proximity", help="Also run proximity analysis.")] = False,
    n_permutations: Annotated[int, typer.Option("--n-permutations", help="Proximity random permutations.")] = 1000,
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Replace existing niche outputs.")] = False,
) -> None:
    """Run niche interaction analysis in one go."""

    radius_value, _ = _parse_radius(radius)
    config = NicheInteractionConfig(
        anndata_filename=anndata_filename,
        label_key=label_key,
        spatial_key=spatial_key,
        neighborhood=NeighborhoodConfig(
            radius=radius_value,
            epsilon_threshold=epsilon_threshold,
            additional_excluded_cell_types=tuple(exclude_cell_type or ()),
        ),
        model=InteractionModelConfig(
            k_fold=k_fold,
            n_repeats=n_repeats,
            seed=seed,
            n_jobs=n_jobs,
            c_values=tuple(c_value) if c_value else None,
        ),
        plots=NichePlotConfig(
            enabled=make_plots,
            kinds=tuple(plot_kind) if plot_kind else ("confusion", "coefficients", "scores", "graph"),
            saveas=plot_format,
            interaction_cutoff=interaction_cutoff,
        ),
        proximity=ProximityConfig(enabled=proximity, n_permutations=n_permutations, saveas=plot_format),
        overwrite=overwrite,
    )
    try:
        result = run_niche_interactions(output_dir=output_dir, config=config)
    except ValidationError as exc:
        typer.secho(f"validation error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc

    _echo_result(result)


@app.command("validate")
def validate(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing transferred labels.")],
    anndata_filename: Annotated[str, typer.Option("--anndata-filename", help="Annotated spatial AnnData file name.")] = "nico_celltype_annotation.h5ad",
    label_key: Annotated[str, typer.Option("--label-key", help=".obs label column.")] = "nico_ct",
    spatial_key: Annotated[str, typer.Option("--spatial-key", help=".obsm coordinate key.")] = "spatial",
    radius: Annotated[str, typer.Option("--radius", help="Radius value to validate planned artifacts for.")] = "0",
    epsilon_threshold: Annotated[float, typer.Option("--epsilon-threshold", help="Delaunay edge threshold.")] = 100.0,
    exclude_cell_type: Annotated[
        Optional[list[str]],
        typer.Option("--exclude-cell-type", help="Additional cell type to exclude. May be repeated."),
    ] = None,
    k_fold: Annotated[int, typer.Option("--k-fold", help="Cross-validation folds.")] = 5,
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Allow existing planned outputs.")] = False,
) -> None:
    """Validate niche interaction inputs without running NiCo."""

    radius_value, radius_tag = _parse_radius(radius)
    config = NicheInteractionConfig(
        anndata_filename=anndata_filename,
        label_key=label_key,
        spatial_key=spatial_key,
        neighborhood=NeighborhoodConfig(
            radius=radius_value,
            epsilon_threshold=epsilon_threshold,
            additional_excluded_cell_types=tuple(exclude_cell_type or ()),
        ),
        model=InteractionModelConfig(k_fold=k_fold),
        overwrite=overwrite,
    )
    try:
        validate_niche_inputs(output_dir, config=config)
    except ValidationError as exc:
        typer.secho(f"validation error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc

    adata = read_h5ad(output_dir / anndata_filename)
    cell_types = adata.obs[label_key].astype(str).nunique()
    typer.echo(f"validated: {output_dir / anndata_filename}")
    typer.echo(f"cells: {adata.n_obs}")
    typer.echo(f"cell_types: {cell_types}")
    typer.echo(f"radius: {radius_tag}")
    typer.echo(f"planned_prediction_dir: {output_dir / 'niche_prediction_linear'}")


@app.command("artifacts")
def artifacts(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing niche outputs.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value to inspect.")] = "0",
    anndata_filename: Annotated[str, typer.Option("--anndata-filename", help="Annotated spatial AnnData file name.")] = "nico_celltype_annotation.h5ad",
    check: Annotated[bool, typer.Option("--check/--no-check", help="Exit nonzero when covariation artifacts are missing.")] = False,
) -> None:
    """Print expected niche artifacts and their existence status."""

    radius_value, _ = _parse_radius(radius)
    paths = planned_artifact_paths(output_dir, radius=radius_value, anndata_filename=anndata_filename)
    all_paths = {
        "used_cell_types_tsv": paths.used_cell_types_tsv,
        "used_clusters_csv": paths.used_clusters_csv,
        "neighbors_pickle": paths.neighbors_pickle,
        "distances_pickle": paths.distances_pickle,
        "normalized_neighborhood_npz": paths.normalized_neighborhood_npz,
        "classifier_matrices_npz": paths.classifier_matrices_npz,
        "manifest_json": paths.manifest_json,
        "metrics_tsv": paths.metrics_tsv,
        "interactions_tsv": paths.interactions_tsv,
    }
    for name, path in all_paths.items():
        status = "ok" if path.exists() else "missing"
        typer.echo(f"{name}: {path} [{status}]")

    try:
        ready = validate_covariation_artifacts(output_dir, radius=radius_value, anndata_filename=anndata_filename)
    except ValidationError as exc:
        ready = False
        if check:
            typer.secho(f"covariation_ready: false ({exc})", fg=typer.colors.RED, err=True)
            raise typer.Exit(1) from exc
    typer.echo(f"covariation_ready: {str(ready).lower()}")


@app.command("export")
def export(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing niche outputs.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value to export.")] = "0",
    output: Annotated[Optional[Path], typer.Option("--output", help="Output TSV path.")] = None,
    cutoff: Annotated[float, typer.Option("--cutoff", help="Positive normalized coefficient cutoff.")] = 0.0,
    include_self_edges: Annotated[
        bool,
        typer.Option("--include-self-edges/--exclude-self-edges", help="Include source=target interactions."),
    ] = True,
) -> None:
    """Export a directed interaction coefficient table from existing artifacts."""

    radius_value, _ = _parse_radius(radius)
    try:
        result = load_niche_result(output_dir, radius=radius_value)
        path = export_interaction_table(result, cutoff=cutoff, include_self_edges=include_self_edges, output_path=output)
    except (ValidationError, ValueError) as exc:
        typer.secho(f"export error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"interactions_tsv: {path}")


@app.command("plot")
def plot(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing niche outputs.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value to plot.")] = "0",
    kind: Annotated[Optional[list[str]], typer.Option("--kind", help="Plot kind. May be repeated.")] = None,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
    interaction_cutoff: Annotated[float, typer.Option("--interaction-cutoff", help="Graph cutoff.")] = 0.1,
    choose_cell_type: Annotated[
        Optional[list[str]],
        typer.Option("--choose-cell-type", help="Cell type to include in top-coefficients plot. May be repeated. Omit for all."),
    ] = None,
) -> None:
    """Generate plots from existing niche artifacts."""

    radius_value, _ = _parse_radius(radius)
    try:
        result = load_niche_result(output_dir, radius=radius_value)
        paths = plot_niche_result(
            result,
            config=NichePlotConfig(
                enabled=True,
                kinds=tuple(kind) if kind else ("confusion", "coefficients", "scores", "graph"),
                saveas=plot_format,
                interaction_cutoff=interaction_cutoff,
                choose_celltypes=tuple(choose_cell_type) if choose_cell_type else (),
            ),
        )
    except (ValidationError, ValueError) as exc:
        typer.secho(f"plot error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    for path in paths:
        typer.echo(f"plot: {path}")


@app.command("proximity")
def proximity_command(
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory containing niche outputs.")],
    radius: Annotated[str, typer.Option("--radius", help="Radius tag/value for proximity.")] = "0",
    n_permutations: Annotated[int, typer.Option("--n-permutations", help="Number of randomized label permutations.")] = 1000,
    seed: Annotated[Optional[int], typer.Option("--seed", help="Optional random seed for permutations.")] = None,
    plot_format: Annotated[str, typer.Option("--plot-format", help="Plot file format.")] = "pdf",
) -> None:
    """Run observed-vs-randomized cell-type proximity analysis."""

    radius_value, _ = _parse_radius(radius)
    try:
        result = load_niche_result(output_dir, radius=radius_value)
        prox = run_proximity_analysis(
            result,
            config=ProximityConfig(n_permutations=n_permutations, seed=seed, saveas=plot_format),
        )
    except (ValidationError, ValueError) as exc:
        typer.secho(f"proximity error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    if prox.observed_tsv:
        typer.echo(f"observed_tsv: {prox.observed_tsv}")
    if prox.ratio_tsv:
        typer.echo(f"ratio_tsv: {prox.ratio_tsv}")
    if prox.plot_path:
        typer.echo(f"plot: {prox.plot_path}")


def _parse_radius(value: str) -> tuple[int | float, str]:
    try:
        return normalize_radius(value)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _echo_result(result: NicheInteractionResult) -> None:
    typer.echo(f"annotated_h5ad: {result.annotated_h5ad}")
    typer.echo(f"used_cell_types_tsv: {result.used_cell_types_tsv}")
    typer.echo(f"used_clusters_csv: {result.used_clusters_csv}")
    typer.echo(f"neighbors_pickle: {result.neighbors_pickle}")
    typer.echo(f"distances_pickle: {result.distances_pickle}")
    typer.echo(f"normalized_neighborhood_npz: {result.normalized_neighborhood_npz}")
    typer.echo(f"classifier_matrices_npz: {result.classifier_matrices_npz}")
    if result.metrics_tsv is not None:
        typer.echo(f"metrics_tsv: {result.metrics_tsv}")
    if result.interactions_tsv is not None:
        typer.echo(f"interactions_tsv: {result.interactions_tsv}")
    if result.manifest_json is not None:
        typer.echo(f"manifest_json: {result.manifest_json}")
    if result.selected_c is not None:
        typer.echo(f"selected_c: {result.selected_c}")
    typer.echo(f"covariation_ready: {str(result.covariation_ready).lower()}")


if __name__ == "__main__":
    app()
