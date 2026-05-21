"""Typer command line app for NiCo preprocessing."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

from .config import NiCoSCTransformConfig, PearsonResidualsConfig
from .io import load_reference_from_sparse, load_spatial_query_from_csv
from .pipeline import preprocess_nico_inputs

app = typer.Typer(
    name="nico-preprocess",
    help="Convert raw reference/spatial inputs and build NiCo-ready preprocessing files.",
    no_args_is_help=True,
)


NormalizationName = Annotated[
    str,
    typer.Option(
        "--normalization",
        help="Normalization method: 'nico-sctransform' or 'pearson-residuals'.",
    ),
]


def _normalization_config(
    normalization: str,
    *,
    sct_min_cells: int,
    sct_gmean_eps: float,
    sct_n_genes: int | None,
    sct_n_cells: int | None,
    sct_bin_size: int,
    sct_bw_adjust: float,
    pearson_theta: float,
    pearson_clip: float | None,
    pearson_check_values: bool,
    pearson_layer: str | None,
) -> NiCoSCTransformConfig | PearsonResidualsConfig:
    """Create a normalization config from CLI options.

    Method-specific options are intentionally accepted by the build command even
    when inactive. The selected ``--normalization`` value determines which group
    is used.
    """

    if normalization == "nico-sctransform":
        return NiCoSCTransformConfig(
            min_cells=sct_min_cells,
            gmean_eps=sct_gmean_eps,
            n_genes=sct_n_genes,
            n_cells=sct_n_cells,
            bin_size=sct_bin_size,
            bw_adjust=sct_bw_adjust,
        )
    if normalization == "pearson-residuals":
        return PearsonResidualsConfig(
            theta=pearson_theta,
            clip=pearson_clip,
            check_values=pearson_check_values,
            layer=pearson_layer,
        )
    raise typer.BadParameter("Expected 'nico-sctransform' or 'pearson-residuals'.")


@app.command("convert-reference-sparse")
def convert_reference_sparse(
    counts: Annotated[Path, typer.Option("--counts", help="Sparse triplet count table.")],
    genes: Annotated[Path, typer.Option("--genes", help="Gene metadata CSV.")],
    barcodes: Annotated[Path, typer.Option("--barcodes", help="Cell barcode CSV.")],
    output: Annotated[Path, typer.Option("--output", help="Output reference .h5ad path.")],
    sep: Annotated[str, typer.Option("--sep", help="Separator used by the sparse count table.")] = " ",
    header: Annotated[Optional[int], typer.Option("--header", help="Header row for the sparse count table.")] = 1,
    one_based_indices: Annotated[
        bool,
        typer.Option("--one-based-indices/--zero-based-indices", help="Whether sparse indices are one-based."),
    ] = True,
) -> None:
    """Convert sparse triplet reference files to a raw reference h5ad."""

    adata = load_reference_from_sparse(
        counts_path=counts,
        genes_path=genes,
        barcodes_path=barcodes,
        sep=sep,
        header=header,
        one_based_indices=one_based_indices,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)


@app.command("convert-spatial-csv")
def convert_spatial_csv(
    counts: Annotated[Path, typer.Option("--counts", help="Spatial count CSV.")],
    coordinates: Annotated[Path, typer.Option("--coordinates", help="Spatial coordinate CSV.")],
    output: Annotated[Path, typer.Option("--output", help="Output spatial .h5ad path.")],
    counts_orientation: Annotated[
        str,
        typer.Option("--counts-orientation", help="'genes_by_cells' or 'cells_by_genes'."),
    ] = "genes_by_cells",
    barcode_col: Annotated[str, typer.Option("--barcode-col", help="Coordinate barcode column name or index.")] = "0",
    coordinate_cols: Annotated[
        Optional[list[str]],
        typer.Option("--coordinate-col", help="Coordinate column name/index. May be repeated."),
    ] = None,
    spatial_key: Annotated[str, typer.Option("--spatial-key", help="Key to store coordinates in .obsm.")] = "spatial",
    reorder_coordinates: Annotated[
        bool,
        typer.Option("--reorder-coordinates/--no-reorder-coordinates", help="Reorder coordinates by barcode."),
    ] = True,
) -> None:
    """Convert spatial/Xenium count and coordinate CSVs to a raw spatial h5ad."""

    adata = load_spatial_query_from_csv(
        counts_path=counts,
        coordinates_path=coordinates,
        counts_orientation=counts_orientation,  # type: ignore[arg-type]
        barcode_col=_parse_column_selector(barcode_col),
        coordinate_cols=[_parse_column_selector(col) for col in coordinate_cols] if coordinate_cols else None,
        spatial_key=spatial_key,
        reorder_coordinates=reorder_coordinates,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)


@app.command("build")
def build(
    reference: Annotated[Path, typer.Option("--reference", help="Raw reference .h5ad path.")],
    spatial: Annotated[Path, typer.Option("--spatial", help="Raw spatial/query .h5ad path.")],
    ref_out_dir: Annotated[Path, typer.Option("--ref-out-dir", help="Directory for reference NiCo outputs.")],
    spatial_out_dir: Annotated[Path, typer.Option("--spatial-out-dir", help="Directory for spatial NiCo outputs.")],
    normalization: NormalizationName = "nico-sctransform",
    spatial_key: Annotated[str, typer.Option("--spatial-key", help="Key containing spatial coordinates in .obsm.")] = "spatial",
    ref_label_key: Annotated[str, typer.Option("--ref-label-key", help="Reference cell type/cluster column.")] = "cluster",
    min_cell_counts: Annotated[int, typer.Option("--min-cell-counts", help="Minimum total counts per retained cell.")] = 5,
    min_gene_cells: Annotated[int, typer.Option("--min-gene-cells", help="Minimum cells per retained gene.")] = 1,
    gene_space: Annotated[str, typer.Option("--gene-space", help="'shared' or 'reference_all'.")] = "shared",
    spatial_n_pcs: Annotated[int, typer.Option("--spatial-n-pcs", help="PCs for spatial neighbors.")] = 30,
    leiden_resolutions: Annotated[
        Optional[list[float]],
        typer.Option("--leiden-resolution", help="Leiden resolution. May be repeated."),
    ] = None,
    make_reference_umap: Annotated[
        bool,
        typer.Option("--make-reference-umap/--no-make-reference-umap", help="Compute reference UMAP in Original_counts."),
    ] = True,
    random_state: Annotated[int, typer.Option("--random-state", help="Random seed for Scanpy steps.")] = 0,
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Replace existing output files.")] = False,
    # NiCo SCTransform-specific options.
    sct_min_cells: Annotated[int, typer.Option("--sct-min-cells", help="SCTransform-only: min cells per gene.")] = 1,
    sct_gmean_eps: Annotated[float, typer.Option("--sct-gmean-eps", help="SCTransform-only: geometric mean epsilon.")] = 1.0,
    sct_n_genes: Annotated[Optional[int], typer.Option("--sct-n-genes", help="SCTransform-only: genes sampled for fitting.")] = 500,
    sct_n_cells: Annotated[Optional[int], typer.Option("--sct-n-cells", help="SCTransform-only: cells sampled for fitting.")] = None,
    sct_bin_size: Annotated[int, typer.Option("--sct-bin-size", help="SCTransform-only: genes per fitting bin.")] = 500,
    sct_bw_adjust: Annotated[float, typer.Option("--sct-bw-adjust", help="SCTransform-only: bandwidth adjustment.")] = 3.0,
    # Pearson residual-specific options.
    pearson_theta: Annotated[float, typer.Option("--pearson-theta", help="Pearson-only: overdispersion theta.")] = 100.0,
    pearson_clip: Annotated[Optional[float], typer.Option("--pearson-clip", help="Pearson-only: clipping threshold.")] = None,
    pearson_check_values: Annotated[
        bool,
        typer.Option("--pearson-check-values/--no-pearson-check-values", help="Pearson-only: validate count values."),
    ] = True,
    pearson_layer: Annotated[Optional[str], typer.Option("--pearson-layer", help="Pearson-only: layer to normalize.")] = None,
) -> None:
    """Build NiCo-ready h5ad files from raw reference and spatial h5ad inputs."""

    norm_config = _normalization_config(
        normalization,
        sct_min_cells=sct_min_cells,
        sct_gmean_eps=sct_gmean_eps,
        sct_n_genes=sct_n_genes,
        sct_n_cells=sct_n_cells,
        sct_bin_size=sct_bin_size,
        sct_bw_adjust=sct_bw_adjust,
        pearson_theta=pearson_theta,
        pearson_clip=pearson_clip,
        pearson_check_values=pearson_check_values,
        pearson_layer=pearson_layer,
    )
    outputs = preprocess_nico_inputs(
        reference_h5ad=reference,
        spatial_h5ad=spatial,
        ref_out_dir=ref_out_dir,
        spatial_out_dir=spatial_out_dir,
        spatial_key=spatial_key,
        ref_label_key=ref_label_key,
        normalization=norm_config,
        min_cell_counts=min_cell_counts,
        min_gene_cells=min_gene_cells,
        gene_space=gene_space,  # type: ignore[arg-type]
        spatial_n_pcs=spatial_n_pcs,
        leiden_resolutions=leiden_resolutions or [0.4, 0.5],
        make_reference_umap=make_reference_umap,
        random_state=random_state,
        overwrite=overwrite,
    )
    for name, path in outputs.items():
        typer.echo(f"{name}: {path}")


def _parse_column_selector(value: str) -> str | int:
    """Parse a CLI column selector as an integer index when possible."""

    try:
        return int(value)
    except ValueError:
        return value


if __name__ == "__main__":
    app()
