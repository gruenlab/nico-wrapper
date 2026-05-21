"""Typer command line app for NiCo label transfer."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

from .config import AnchorConfig, AnnotationConfig, LabelTransferConfig, TieStrategy
from .pipeline import run_label_transfer

app = typer.Typer(
    name="nico-transfer",
    help="Transfer reference labels onto NiCo-preprocessed spatial/query data.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Transfer reference labels onto NiCo-preprocessed spatial/query data."""


@app.command("run")
def run(
    ref_dir: Annotated[
        Path,
        typer.Option(
            "--ref-dir",
            help="Directory from nico-preprocess containing Original_counts.h5ad and sct_singleCell.h5ad.",
        ),
    ],
    spatial_dir: Annotated[
        Path,
        typer.Option("--spatial-dir", help="Directory from nico-preprocess containing sct_spatial.h5ad."),
    ],
    output_dir: Annotated[Path, typer.Option("--output-dir", help="Directory for NiCo transfer outputs.")],
    annotation_dir: Annotated[
        Optional[Path],
        typer.Option("--annotation-dir", help="Directory for annotation intermediates; defaults to output-dir/annotations."),
    ] = None,
    ref_label_key: Annotated[
        str,
        typer.Option("--ref-label-key", help="Reference .obs column containing cell labels."),
    ] = "cluster",
    spatial_cluster_key: Annotated[
        str,
        typer.Option("--spatial-cluster-key", help="Spatial .obs guide cluster column, e.g. leiden0.5."),
    ] = "leiden0.5",
    neighbors: Annotated[int, typer.Option("--neighbors", help="K for MNN anchors and spatial KNN graph.")] = 50,
    n_pcs: Annotated[int, typer.Option("--n-pcs", help="Number of PCs for the transfer space.")] = 50,
    minkowski_order: Annotated[
        int,
        typer.Option("--minkowski-order", help="Minkowski distance order; 2 is Euclidean."),
    ] = 2,
    dispersion_cutoff: Annotated[
        float,
        typer.Option("--dispersion-cutoff", help="Anchor pruning cutoff across spatial guide clusters."),
    ] = 0.15,
    iterations: Annotated[int, typer.Option("--iterations", help="Label propagation iterations.")] = 3,
    tie_strategy: Annotated[
        TieStrategy,
        typer.Option("--tie-strategy", help="How to resolve ambiguous labels: majority or weighted."),
    ] = TieStrategy.MAJORITY,
    spatial_sct_filename: Annotated[
        str,
        typer.Option("--spatial-sct-filename", help="Normalized spatial/query AnnData file name."),
    ] = "sct_spatial.h5ad",
    sc_sct_filename: Annotated[
        str,
        typer.Option("--sc-sct-filename", help="Normalized single-cell/reference AnnData file name."),
    ] = "sct_singleCell.h5ad",
    sc_full_filename: Annotated[
        str,
        typer.Option("--sc-full-filename", help="Full/original reference AnnData file name."),
    ] = "Original_counts.h5ad",
    output_h5ad_name: Annotated[
        str,
        typer.Option("--output-h5ad-name", help="Annotated spatial AnnData output file name."),
    ] = "nico_celltype_annotation.h5ad",
    output_label_key: Annotated[
        str,
        typer.Option("--output-label-key", help=".obs column for transferred labels in the output AnnData."),
    ] = "nico_ct",
    overwrite: Annotated[bool, typer.Option("--overwrite/--no-overwrite", help="Replace existing outputs.")] = False,
    cleanup_intermediate: Annotated[
        bool,
        typer.Option("--cleanup-intermediate/--keep-intermediate", help="Remove temporary anchor files after success."),
    ] = False,
) -> None:
    """Run label transfer on outputs produced by nico-preprocess."""

    config = LabelTransferConfig(
        anchors=AnchorConfig(
            neighbors=neighbors,
            n_pcs=n_pcs,
            minkowski_order=minkowski_order,
            spatial_sct_filename=spatial_sct_filename,
            sc_sct_filename=sc_sct_filename,
            sc_full_filename=sc_full_filename,
        ),
        annotation=AnnotationConfig(
            ref_label_key=ref_label_key,
            spatial_cluster_key=spatial_cluster_key,
            dispersion_cutoff=dispersion_cutoff,
            iterations=iterations,
            tie_strategy=tie_strategy,
            output_label_key=output_label_key,
            output_h5ad_name=output_h5ad_name,
        ),
        overwrite=overwrite,
        cleanup_intermediate=cleanup_intermediate,
    )
    outputs = run_label_transfer(
        ref_dir=ref_dir,
        spatial_dir=spatial_dir,
        output_dir=output_dir,
        annotation_dir=annotation_dir,
        config=config,
    )

    typer.echo(f"annotated_h5ad: {outputs.annotated_h5ad}")
    if outputs.anchors_npz is not None:
        typer.echo(f"anchors_npz: {outputs.anchors_npz}")
    for path in outputs.iteration_cluster_csvs:
        typer.echo(f"iteration_cluster_csv: {path}")
    for path in outputs.iteration_celltype_csvs:
        typer.echo(f"iteration_celltype_csv: {path}")


if __name__ == "__main__":
    app()
