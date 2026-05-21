"""Input conversion and loading helpers for NiCo preprocessing."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

from .validation import require_file, validate_sparse_reference_inputs, validate_spatial_csv_inputs, validate_spatial_adata


IndexColumn = str | int


def load_reference_from_sparse(
    counts_path: str | Path,
    genes_path: str | Path,
    barcodes_path: str | Path,
    *,
    annotation_files: dict[str, str | Path] | None = None,
    embedding_files: dict[str, str | Path] | None = None,
    sep: str = " ",
    header: int | None = 1,
    one_based_indices: bool = True,
) -> AnnData:
    """Load a scRNA-seq reference from a sparse triplet count representation.

    The expected count table has three columns: gene index, cell index, and count.
    This function infers the matrix shape from ``genes_path`` and
    ``barcodes_path`` and constructs a sparse ``AnnData`` object with cells as
    observations and genes as variables.

    Parameters
    ----------
    counts_path
        Path to the sparse count table containing ``gene, cell, count`` rows.
    genes_path
        CSV file containing gene identifiers/names. The implementation uses the
        second column as ``adata.var_names`` when present, otherwise the first.
    barcodes_path
        CSV file containing cell barcodes. The implementation uses the second
        column as ``adata.obs_names`` when present, otherwise the first.
    annotation_files
        Optional mapping from desired ``adata.obs`` column names to CSV files.
        Each file is expected to contain one annotation value per cell. The third
        column is used when present, otherwise the last column.
    embedding_files
        Optional mapping from desired ``adata.obsm`` keys, e.g. ``"X_umap"`` or
        ``"X_tsne"``, to CSV files containing embeddings. Columns 2 and 3 are
        used when present, otherwise the last two columns.
    sep
        Separator used by ``counts_path``.
    header
        Header row passed to ``pandas.read_csv`` for ``counts_path``. The sparse
        tutorial uses ``header=1``.
    one_based_indices
        Whether sparse matrix indices are one-based, as in MatrixMarket-style
        exports. If true, indices are shifted down by one.

    Returns
    -------
    AnnData
        Reference data with raw counts in ``.X``.
    """

    validate_sparse_reference_inputs(
        counts_path=counts_path,
        genes_path=genes_path,
        barcodes_path=barcodes_path,
        sep=sep,
        header=header,
        one_based_indices=one_based_indices,
    )

    genes = pd.read_csv(genes_path)
    barcodes = pd.read_csv(barcodes_path)
    gene_names = _infer_names(genes)
    cell_names = _infer_names(barcodes)

    counts = pd.read_csv(counts_path, sep=sep, header=header)
    gene_index = pd.to_numeric(counts.iloc[:, 0], errors="raise").to_numpy(dtype=np.int64)
    cell_index = pd.to_numeric(counts.iloc[:, 1], errors="raise").to_numpy(dtype=np.int64)
    values = pd.to_numeric(counts.iloc[:, 2], errors="raise").to_numpy(dtype=float)

    if one_based_indices:
        gene_index = gene_index - 1
        cell_index = cell_index - 1

    matrix = sparse.coo_matrix(
        (values, (cell_index, gene_index)),
        shape=(len(cell_names), len(gene_names)),
    ).tocsr()
    matrix.sum_duplicates()

    adata = AnnData(matrix)
    adata.obs_names = pd.Index(cell_names).astype(str)
    adata.var_names = pd.Index(gene_names).astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    for key, path in (annotation_files or {}).items():
        adata.obs[key] = _read_annotation_column(path, expected_rows=adata.n_obs)

    for key, path in (embedding_files or {}).items():
        adata.obsm[key] = _read_embedding_matrix(path, expected_rows=adata.n_obs)

    return adata


def load_spatial_query_from_csv(
    counts_path: str | Path,
    coordinates_path: str | Path,
    *,
    counts_orientation: Literal["genes_by_cells", "cells_by_genes"] = "genes_by_cells",
    barcode_col: IndexColumn = 0,
    coordinate_cols: Sequence[IndexColumn] | None = None,
    spatial_key: str = "spatial",
    reorder_coordinates: bool = True,
) -> AnnData:
    """Load a spatial/Xenium query dataset from count and coordinate CSV files.

    Parameters
    ----------
    counts_path
        CSV count matrix. By default rows are genes and columns are cells, as in
        the tutorial's ``gene_by_cell.csv`` file.
    coordinates_path
        CSV file containing cell barcodes and spatial centroid coordinates.
    counts_orientation
        Orientation of the count matrix. ``"genes_by_cells"`` means the matrix
        will be transposed after loading so cells become observations.
        ``"cells_by_genes"`` means cells are already rows.
    barcode_col
        Column in ``coordinates_path`` containing cell barcodes. Can be a column
        name or positional index.
    coordinate_cols
        Columns in ``coordinates_path`` containing spatial coordinates. If
        ``None``, all columns except ``barcode_col`` are used.
    spatial_key
        Key under ``adata.obsm`` where coordinates are stored.
    reorder_coordinates
        If true, reorder coordinate rows to match ``adata.obs_names`` using the
        barcode column. If false, require coordinates to already be in the same
        order as the count matrix.

    Returns
    -------
    AnnData
        Spatial query data with raw counts in ``.X`` and coordinates in
        ``.obsm[spatial_key]``.
    """

    validate_spatial_csv_inputs(
        counts_path=counts_path,
        coordinates_path=coordinates_path,
        counts_orientation=counts_orientation,
        barcode_col=barcode_col,
        coordinate_cols=coordinate_cols,
        reorder_coordinates=reorder_coordinates,
    )

    counts = pd.read_csv(counts_path, index_col=0)
    numeric_counts = counts.apply(pd.to_numeric, errors="raise")

    if counts_orientation == "genes_by_cells":
        matrix = sparse.csr_matrix(numeric_counts.to_numpy(dtype=float).T)
        obs_names = counts.columns.astype(str)
        var_names = counts.index.astype(str)
    elif counts_orientation == "cells_by_genes":
        matrix = sparse.csr_matrix(numeric_counts.to_numpy(dtype=float))
        obs_names = counts.index.astype(str)
        var_names = counts.columns.astype(str)
    else:
        raise ValueError("counts_orientation must be 'genes_by_cells' or 'cells_by_genes'.")

    adata = AnnData(matrix)
    adata.obs_names = pd.Index(obs_names).astype(str)
    adata.var_names = pd.Index(var_names).astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    coordinates = pd.read_csv(coordinates_path)
    barcode_values = _select_column(coordinates, barcode_col).astype(str)
    resolved_coordinate_cols = _resolve_coordinate_cols(coordinates, barcode_col, coordinate_cols)
    coordinate_frame = coordinates.loc[:, resolved_coordinate_cols].apply(pd.to_numeric, errors="raise")
    coordinate_frame.index = barcode_values

    if reorder_coordinates:
        coordinate_frame = coordinate_frame.loc[adata.obs_names]
    elif not np.array_equal(coordinate_frame.index.to_numpy(), adata.obs_names.to_numpy()):
        raise ValueError("Coordinate rows are not in the same order as count matrix cells.")

    adata.obsm[spatial_key] = coordinate_frame.to_numpy(dtype=float)
    validate_spatial_adata(adata, spatial_key=spatial_key)
    return adata


def load_spatial_query_from_h5ad(
    path: str | Path,
    *,
    spatial_key: str = "spatial",
    require_spatial: bool = True,
) -> AnnData:
    """Load a spatial/Xenium query dataset from an existing ``.h5ad`` file.

    This is primarily a thin wrapper around ``scanpy.read_h5ad`` that optionally
    validates the presence and shape of spatial coordinates.

    Parameters
    ----------
    path
        Path to the spatial ``.h5ad`` file.
    spatial_key
        Key under ``adata.obsm`` containing spatial coordinates.
    require_spatial
        If true, raise an error when ``spatial_key`` is missing or has the wrong
        number of rows.

    Returns
    -------
    AnnData
        Loaded spatial query data.
    """

    require_file(path, label="Spatial h5ad file")
    adata = sc.read_h5ad(path)
    if require_spatial:
        validate_spatial_adata(adata, spatial_key=spatial_key, check_counts=False)
    return adata


def _infer_names(frame: pd.DataFrame) -> pd.Index:
    column_position = 1 if frame.shape[1] > 1 else 0
    return pd.Index(frame.iloc[:, column_position].astype(str))


def _read_annotation_column(path: str | Path, *, expected_rows: int) -> pd.Series:
    require_file(path, label="Annotation file")
    frame = pd.read_csv(path)
    if frame.shape[0] != expected_rows:
        raise ValueError(f"Annotation file {path} has {frame.shape[0]} rows; expected {expected_rows}.")
    column_position = 2 if frame.shape[1] > 2 else frame.shape[1] - 1
    return frame.iloc[:, column_position].to_numpy()


def _read_embedding_matrix(path: str | Path, *, expected_rows: int) -> np.ndarray:
    require_file(path, label="Embedding file")
    frame = pd.read_csv(path)
    if frame.shape[0] != expected_rows:
        raise ValueError(f"Embedding file {path} has {frame.shape[0]} rows; expected {expected_rows}.")
    if frame.shape[1] >= 3:
        values = frame.iloc[:, [1, 2]]
    elif frame.shape[1] >= 2:
        values = frame.iloc[:, -2:]
    else:
        raise ValueError(f"Embedding file {path} must contain at least two columns.")
    return values.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)


def _select_column(frame: pd.DataFrame, selector: IndexColumn) -> pd.Series:
    if isinstance(selector, int):
        return frame.iloc[:, selector]
    return frame.loc[:, selector]


def _resolve_coordinate_cols(
    frame: pd.DataFrame,
    barcode_col: IndexColumn,
    coordinate_cols: Sequence[IndexColumn] | None,
) -> list[str]:
    if coordinate_cols is None:
        barcode_name = frame.columns[barcode_col] if isinstance(barcode_col, int) else barcode_col
        return [str(col) for col in frame.columns if col != barcode_name]

    resolved: list[str] = []
    for col in coordinate_cols:
        if isinstance(col, int):
            resolved.append(str(frame.columns[col]))
        else:
            resolved.append(col)
    return resolved
