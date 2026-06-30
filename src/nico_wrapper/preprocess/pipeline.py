"""Joint NiCo preprocessing pipeline."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import scanpy as sc
from anndata import AnnData

from .config import NiCoSCTransformConfig, NormalizationConfig, PearsonResidualsConfig
from .validation import (
    ValidationError,
    require_file,
    validate_expected_files_exist,
    validate_joint_adata,
    validate_leiden_keys,
    validate_no_nan_or_inf,
    validate_normalization_config,
    validate_normalization_layers,
    validate_output_files,
    validate_reference_adata,
    validate_spatial_adata,
)

OutputPaths = dict[str, Path]


def preprocess_nico_inputs(
    reference_h5ad: str | Path,
    spatial_h5ad: str | Path,
    *,
    ref_out_dir: str | Path,
    spatial_out_dir: str | Path,
    spatial_key: str = "spatial",
    ref_label_key: str = "cluster",
    normalization: NormalizationConfig = NiCoSCTransformConfig(),
    ref_min_cell_counts: int = 5,
    spatial_min_cell_counts: int = 5,
    ref_min_gene_cells: int = 1,
    spatial_min_gene_cells: int = 1,
    gene_space: Literal["shared", "reference_all"] = "shared",
    spatial_n_pcs: int = 30,
    leiden_resolutions: Sequence[float] = (0.4, 0.5),
    make_reference_umap: bool = True,
    random_state: int = 0,
    overwrite: bool = False,
) -> OutputPaths:
    """Build NiCo-ready preprocessing outputs from reference and spatial h5ad files.

    This is the required joint step after any optional input conversion. It
    validates both datasets, filters cells/genes, aligns modalities, performs the
    selected normalization strategy, computes spatial PCA/neighbors/UMAP/Leiden,
    and writes the file names expected by downstream NiCo annotation functions.

    Parameters
    ----------
    reference_h5ad
        Path to a reference scRNA-seq ``.h5ad`` file with raw counts in ``.X``
        or in the layer selected by the normalization config.
    spatial_h5ad
        Path to a spatial/Xenium query ``.h5ad`` file with raw counts and spatial
        coordinates under ``.obsm[spatial_key]``.
    ref_out_dir
        Directory where reference outputs are written. The expected outputs are
        ``Original_counts.h5ad`` and ``sct_singleCell.h5ad``.
    spatial_out_dir
        Directory where spatial outputs are written. The expected output is
        ``sct_spatial.h5ad``.
    spatial_key
        Key under ``spatial.obsm`` containing cell coordinates.
    ref_label_key
        Reference cell type / cluster column in ``reference.obs``. This is not
        transformed here, but validating it early helps catch downstream NiCo
        annotation failures.
    normalization
        Normalization configuration. Use ``NiCoSCTransformConfig`` for the
        notebook/paper-style SCTransform path or ``PearsonResidualsConfig`` for
        Scanpy analytic Pearson residuals.
    ref_min_cell_counts
        Minimum total counts required for a cell to be retained before modality
        alignment in reference data.
    spatial_min_cell_counts
        Minimum total counts required for a cell to be retained before modality
        alignment in spatial data.
    ref_min_gene_cells
        Minimum number of cells in which a gene must be observed before being
        retained before modality alignment in reference data.
    spatial_min_gene_cells
        Minimum number of cells in which a gene must be observed before being
        retained before modality alignment in spatial data.
    gene_space
        ``"shared"`` subsets both modalities to common genes before
        normalization. ``"reference_all"`` preserves all reference genes while
        subsetting only spatial data to shared genes; this mirrors one variant in
        the original notebook, but ``"shared"`` is the safer default.
    spatial_n_pcs
        Number of principal components used for spatial neighbor graph
        construction.
    leiden_resolutions
        Leiden clustering resolutions to compute on the normalized spatial data.
        Results should be stored as ``leiden{resolution}``, e.g. ``leiden0.5``.
    make_reference_umap
        Whether to compute normalized/log reference PCA/neighbors/UMAP in
        ``Original_counts.h5ad`` for downstream visualization/reference checks.
    random_state
        Random seed passed to stochastic Scanpy steps where applicable.
    overwrite
        Whether existing output files may be replaced.

    Returns
    -------
    dict[str, Path]
        Mapping with at least ``original_counts``, ``sct_single_cell``, and
        ``sct_spatial`` output paths.
    """

    _validate_preprocessing_arguments(
        reference_h5ad=reference_h5ad,
        spatial_h5ad=spatial_h5ad,
        ref_out_dir=ref_out_dir,
        spatial_out_dir=spatial_out_dir,
        normalization=normalization,
        ref_min_cell_counts=ref_min_cell_counts,
        spatial_min_cell_counts=spatial_min_cell_counts,
        ref_min_gene_cells=ref_min_gene_cells,
        spatial_min_gene_cells=spatial_min_gene_cells,
        gene_space=gene_space,
        spatial_n_pcs=spatial_n_pcs,
        leiden_resolutions=leiden_resolutions,
        overwrite=overwrite,
    )

    reference = sc.read_h5ad(reference_h5ad)
    spatial = sc.read_h5ad(spatial_h5ad)

    source_layer = (
        normalization.layer
        if isinstance(normalization, PearsonResidualsConfig)
        else None
    )
    source_layer = (
        normalization.layer
        if isinstance(normalization, PearsonResidualsConfig)
        else None
    )
    validate_normalization_layers(reference, spatial, normalization)
    validate_reference_adata(reference, ref_label_key=ref_label_key, layer=source_layer)
    validate_spatial_adata(spatial, spatial_key=spatial_key, layer=source_layer)

    reference = reference.copy()
    spatial = spatial.copy()

    _filter_counts(
        reference,
        min_cell_counts=ref_min_cell_counts,
        min_gene_cells=ref_min_gene_cells,
    )
    _filter_counts(
        spatial,
        min_cell_counts=spatial_min_cell_counts,
        min_gene_cells=spatial_min_gene_cells,
    )

    validate_reference_adata(reference, ref_label_key=ref_label_key, layer=source_layer)
    validate_spatial_adata(spatial, spatial_key=spatial_key, layer=source_layer)
    validate_joint_adata(reference, spatial, gene_space=gene_space)

    ref_out_dir = Path(ref_out_dir)
    spatial_out_dir = Path(spatial_out_dir)
    ref_out_dir.mkdir(parents=True, exist_ok=True)
    spatial_out_dir.mkdir(parents=True, exist_ok=True)

    output_paths: OutputPaths = {
        "original_counts": ref_out_dir / "Original_counts.h5ad",
        "sct_single_cell": ref_out_dir / "sct_singleCell.h5ad",
        "sct_spatial": spatial_out_dir / "sct_spatial.h5ad",
    }

    original_counts = _make_original_counts(
        reference, make_reference_umap=make_reference_umap, random_state=random_state
    )
    original_counts.write_h5ad(output_paths["original_counts"])

    reference_common, spatial_common = _align_gene_space(
        reference, spatial, gene_space=gene_space
    )
    # Re-filter after gene-space alignment. Some cells can pass whole-transcriptome
    # filtering but have too few/zero counts in the shared gene space, which can
    # make SCTransform/Pearson residual normalization unstable.
    _filter_counts(
        reference_common,
        min_cell_counts=ref_min_cell_counts,
        min_gene_cells=ref_min_gene_cells,
    )
    _filter_counts(
        spatial_common,
        min_cell_counts=spatial_min_cell_counts,
        min_gene_cells=spatial_min_gene_cells,
    )
    validate_joint_adata(reference_common, spatial_common, gene_space=gene_space)

    reference_normalized, spatial_normalized = _normalize_modalities(
        reference_common,
        spatial_common,
        normalization=normalization,
        spatial_key=spatial_key,
    )

    validate_no_nan_or_inf(reference_normalized, label="Normalized reference")
    validate_no_nan_or_inf(spatial_normalized, label="Normalized spatial")

    _compute_spatial_embedding_and_clusters(
        spatial_normalized,
        n_pcs=spatial_n_pcs,
        leiden_resolutions=leiden_resolutions,
        random_state=random_state,
    )
    validate_leiden_keys(spatial_normalized, leiden_resolutions)

    reference_normalized.write_h5ad(output_paths["sct_single_cell"])
    spatial_normalized.write_h5ad(output_paths["sct_spatial"])
    validate_expected_files_exist(output_paths.values())
    return output_paths


def _validate_preprocessing_arguments(
    *,
    reference_h5ad: str | Path,
    spatial_h5ad: str | Path,
    ref_out_dir: str | Path,
    spatial_out_dir: str | Path,
    normalization: NormalizationConfig,
    ref_min_cell_counts: int,
    spatial_min_cell_counts: int,
    ref_min_gene_cells: int,
    spatial_min_gene_cells: int,
    gene_space: str,
    spatial_n_pcs: int,
    leiden_resolutions: Sequence[float],
    overwrite: bool,
) -> None:
    require_file(reference_h5ad, label="Reference h5ad file")
    require_file(spatial_h5ad, label="Spatial h5ad file")
    validate_normalization_config(normalization)
    if (ref_min_cell_counts < 0) or (spatial_min_cell_counts < 0):
        raise ValidationError(
            "ref_min_cell_counts and spatial_min_cell_counts must be >= 0."
        )
    if (ref_min_gene_cells < 0) or (spatial_min_gene_cells < 0):
        raise ValidationError(
            "ref_min_gene_cells and spatial_min_gene_cells must be >= 0."
        )
    if spatial_n_pcs <= 0:
        raise ValidationError("spatial_n_pcs must be > 0.")
    if not leiden_resolutions:
        raise ValidationError("At least one Leiden resolution is required.")
    if any(resolution <= 0 for resolution in leiden_resolutions):
        raise ValidationError("Leiden resolutions must be > 0.")
    if gene_space not in {"shared", "reference_all"}:
        raise ValidationError("gene_space must be 'shared' or 'reference_all'.")

    ref_out_dir = Path(ref_out_dir)
    spatial_out_dir = Path(spatial_out_dir)
    validate_output_files(
        [
            ref_out_dir / "Original_counts.h5ad",
            ref_out_dir / "sct_singleCell.h5ad",
            spatial_out_dir / "sct_spatial.h5ad",
        ],
        overwrite=overwrite,
    )


def _filter_counts(
    adata: AnnData, *, min_cell_counts: int, min_gene_cells: int
) -> None:
    if min_cell_counts > 0:
        sc.pp.filter_cells(adata, min_counts=min_cell_counts)
    if min_gene_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_gene_cells)
    if adata.n_obs == 0:
        raise ValidationError("Filtering removed all cells.")
    if adata.n_vars == 0:
        raise ValidationError("Filtering removed all genes.")


def _make_original_counts(
    reference: AnnData, *, make_reference_umap: bool, random_state: int
) -> AnnData:
    original_counts = reference.copy()
    original_counts.raw = original_counts.copy()

    if make_reference_umap:
        analysis = original_counts.copy()
        sc.pp.normalize_total(analysis)
        sc.pp.log1p(analysis)
        _run_pca_neighbors_umap(analysis, requested_n_pcs=50, random_state=random_state)
        for key in ("X_pca", "X_umap"):
            if key in analysis.obsm:
                original_counts.obsm[key] = analysis.obsm[key]
        for key in ("pca", "neighbors", "umap"):
            if key in analysis.uns:
                original_counts.uns[key] = analysis.uns[key]
        for key in ("PCs",):
            if key in analysis.varm:
                original_counts.varm[key] = analysis.varm[key]
        for key in ("connectivities", "distances"):
            if key in analysis.obsp:
                original_counts.obsp[key] = analysis.obsp[key]

    return original_counts


def _align_gene_space(
    reference: AnnData,
    spatial: AnnData,
    *,
    gene_space: Literal["shared", "reference_all"],
) -> tuple[AnnData, AnnData]:
    reference_gene_to_index = {
        gene: idx for idx, gene in enumerate(reference.var_names)
    }
    spatial_indices: list[int] = []
    reference_indices: list[int] = []
    for spatial_idx, gene in enumerate(spatial.var_names):
        reference_idx = reference_gene_to_index.get(gene)
        if reference_idx is not None:
            spatial_indices.append(spatial_idx)
            reference_indices.append(reference_idx)

    if not spatial_indices:
        raise ValidationError(
            "No shared genes found between reference and spatial data."
        )

    spatial_common = spatial[:, spatial_indices].copy()
    if gene_space == "shared":
        reference_common = reference[:, reference_indices].copy()
    elif gene_space == "reference_all":
        reference_common = reference.copy()
    else:
        raise ValidationError("gene_space must be 'shared' or 'reference_all'.")

    return reference_common, spatial_common


def _normalize_modalities(
    reference: AnnData,
    spatial: AnnData,
    *,
    normalization: NormalizationConfig,
    spatial_key: str,
) -> tuple[AnnData, AnnData]:
    if isinstance(normalization, NiCoSCTransformConfig):
        return _normalize_with_nico_sctransform(
            reference, spatial, normalization=normalization, spatial_key=spatial_key
        )
    if isinstance(normalization, PearsonResidualsConfig):
        return _normalize_with_pearson_residuals(
            reference, spatial, normalization=normalization
        )
    raise ValidationError(
        f"Unsupported normalization config type: {type(normalization).__name__}."
    )


def _normalize_with_nico_sctransform(
    reference: AnnData,
    spatial: AnnData,
    *,
    normalization: NiCoSCTransformConfig,
    spatial_key: str,
) -> tuple[AnnData, AnnData]:
    try:
        from nico import Annotations as sann
    except ImportError as exc:  # pragma: no cover - environment-specific fallback
        raise ImportError(
            "NiCo SCTransform normalization requires the 'nico' package."
        ) from exc

    reference_raw = reference.copy()
    spatial_raw = spatial.copy()

    reference_normalized = sann.SCTransform(
        reference,
        min_cells=normalization.min_cells,
        gmean_eps=normalization.gmean_eps,
        n_genes=normalization.n_genes,
        n_cells=normalization.n_cells,
        bin_size=normalization.bin_size,
        bw_adjust=normalization.bw_adjust,
        inplace=False,
    )
    spatial_normalized = sann.SCTransform(
        spatial,
        min_cells=normalization.min_cells,
        gmean_eps=normalization.gmean_eps,
        n_genes=normalization.n_genes,
        n_cells=normalization.n_cells,
        bin_size=normalization.bin_size,
        bw_adjust=normalization.bw_adjust,
        inplace=False,
    )

    reference_normalized.raw = reference_raw
    spatial_normalized.raw = spatial_raw
    if spatial_key in spatial_raw.obsm:
        spatial_normalized.obsm[spatial_key] = spatial_raw.obsm[spatial_key].copy()
    return reference_normalized, spatial_normalized


def _normalize_with_pearson_residuals(
    reference: AnnData,
    spatial: AnnData,
    *,
    normalization: PearsonResidualsConfig,
) -> tuple[AnnData, AnnData]:
    reference_normalized = reference.copy()
    spatial_normalized = spatial.copy()
    reference_normalized.raw = reference.copy()
    spatial_normalized.raw = spatial.copy()

    if normalization.layer is not None:
        reference_normalized.X = reference_normalized.layers[normalization.layer].copy()
        spatial_normalized.X = spatial_normalized.layers[normalization.layer].copy()

    for adata, label in (
        (reference_normalized, "reference"),
        (spatial_normalized, "spatial"),
    ):
        with warnings.catch_warnings():
            warnings.simplefilter("default")
            sc.experimental.pp.normalize_pearson_residuals(
                adata,
                theta=normalization.theta,
                clip=normalization.clip,
                check_values=normalization.check_values,
                layer=None,
                inplace=True,
            )
        if normalization.layer is not None:
            adata.uns.setdefault("pearson_residuals_normalization", {})[
                "source_layer"
            ] = normalization.layer
        adata.uns.setdefault("nico_wrapper_preprocess", {})["normalization"] = (
            f"pearson_residuals:{label}"
        )

    return reference_normalized, spatial_normalized


def _compute_spatial_embedding_and_clusters(
    spatial: AnnData,
    *,
    n_pcs: int,
    leiden_resolutions: Sequence[float],
    random_state: int,
) -> None:
    _run_pca_neighbors_umap(spatial, requested_n_pcs=n_pcs, random_state=random_state)
    for resolution in leiden_resolutions:
        sc.tl.leiden(
            spatial,
            resolution=resolution,
            key_added=f"leiden{resolution}",
            random_state=random_state,
        )


def _run_pca_neighbors_umap(
    adata: AnnData, *, requested_n_pcs: int, random_state: int
) -> int:
    if n_comps < 1:
        raise ValidationError(
            "PCA requires at least two cells and two genes after filtering/alignment. "
            f"Observed {adata.n_obs} cells and {adata.n_vars} genes."
        )
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)
    sc.pp.neighbors(adata, n_pcs=n_comps)
    sc.tl.umap(adata, random_state=random_state)
    return n_comps
