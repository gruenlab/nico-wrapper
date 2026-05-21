"""Validation helpers for NiCo preprocessing inputs and outputs.

This module intentionally contains validation logic only. It does not convert
files, mutate AnnData objects, normalize data, or write preprocessing outputs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence
import math
import warnings

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse

from .config import NiCoSCTransformConfig, NormalizationConfig, PearsonResidualsConfig


ColumnSelector = str | int


class ValidationError(ValueError):
    """Raised when preprocessing inputs violate required assumptions."""


def require_file(path: str | Path, *, label: str = "file") -> Path:
    """Validate that ``path`` exists and is a regular file."""

    resolved = Path(path)
    if not resolved.exists():
        raise ValidationError(f"{label} does not exist: {resolved}")
    if not resolved.is_file():
        raise ValidationError(f"{label} is not a file: {resolved}")
    return resolved


def validate_output_files(paths: Sequence[str | Path], *, overwrite: bool = False) -> None:
    """Validate output file paths before writing.

    Parameters
    ----------
    paths
        Output files that will be written.
    overwrite
        If false, existing output files are rejected.
    """

    for raw_path in paths:
        path = Path(raw_path)
        if path.exists() and not overwrite:
            raise ValidationError(f"Output file already exists: {path}. Use overwrite=True to replace it.")
        parent = path.parent if path.parent != Path("") else Path(".")
        if parent.exists() and not parent.is_dir():
            raise ValidationError(f"Output parent exists but is not a directory: {parent}")


def validate_sparse_reference_inputs(
    counts_path: str | Path,
    genes_path: str | Path,
    barcodes_path: str | Path,
    *,
    sep: str = " ",
    header: int | None = 1,
    one_based_indices: bool = True,
) -> None:
    """Validate sparse-triplet reference input files.

    The count table is expected to contain at least three columns representing
    gene index, cell index, and count. Matrix dimensions are inferred from the
    number of rows in ``genes_path`` and ``barcodes_path``.
    """

    counts_path = require_file(counts_path, label="Sparse reference counts file")
    genes_path = require_file(genes_path, label="Reference genes file")
    barcodes_path = require_file(barcodes_path, label="Reference barcodes file")

    genes = pd.read_csv(genes_path)
    barcodes = pd.read_csv(barcodes_path)
    if genes.shape[0] == 0:
        raise ValidationError("Reference genes file contains no genes.")
    if barcodes.shape[0] == 0:
        raise ValidationError("Reference barcodes file contains no barcodes.")

    gene_names = _infer_name_column_values(genes, label="gene")
    barcode_names = _infer_name_column_values(barcodes, label="barcode")
    _warn_if_duplicate_names(gene_names, label="gene names")
    _warn_if_duplicate_names(barcode_names, label="cell barcodes")

    counts = pd.read_csv(counts_path, sep=sep, header=header)
    if counts.shape[1] < 3:
        raise ValidationError(
            f"Sparse reference counts file must contain at least 3 columns; found {counts.shape[1]}."
        )
    if counts.shape[0] == 0:
        raise ValidationError("Sparse reference counts file contains no nonzero entries.")

    gene_index = pd.to_numeric(counts.iloc[:, 0], errors="coerce")
    cell_index = pd.to_numeric(counts.iloc[:, 1], errors="coerce")
    values = pd.to_numeric(counts.iloc[:, 2], errors="coerce")

    if gene_index.isna().any() or cell_index.isna().any():
        raise ValidationError("Sparse reference gene/cell indices must be numeric.")
    if values.isna().any():
        raise ValidationError("Sparse reference counts must be numeric.")
    if not _is_integer_like(gene_index.to_numpy()) or not _is_integer_like(cell_index.to_numpy()):
        raise ValidationError("Sparse reference gene/cell indices must be integer-like.")
    if (values < 0).any():
        raise ValidationError("Sparse reference counts contain negative values.")

    gene_index_arr = gene_index.to_numpy(dtype=np.int64)
    cell_index_arr = cell_index.to_numpy(dtype=np.int64)
    if one_based_indices:
        min_gene, min_cell = 1, 1
        max_gene, max_cell = genes.shape[0], barcodes.shape[0]
    else:
        min_gene, min_cell = 0, 0
        max_gene, max_cell = genes.shape[0] - 1, barcodes.shape[0] - 1

    if gene_index_arr.min() < min_gene or gene_index_arr.max() > max_gene:
        raise ValidationError(
            "Sparse reference gene index out of bounds: "
            f"observed [{gene_index_arr.min()}, {gene_index_arr.max()}], "
            f"allowed [{min_gene}, {max_gene}]."
        )
    if cell_index_arr.min() < min_cell or cell_index_arr.max() > max_cell:
        raise ValidationError(
            "Sparse reference cell index out of bounds: "
            f"observed [{cell_index_arr.min()}, {cell_index_arr.max()}], "
            f"allowed [{min_cell}, {max_cell}]."
        )


def validate_spatial_csv_inputs(
    counts_path: str | Path,
    coordinates_path: str | Path,
    *,
    counts_orientation: Literal["genes_by_cells", "cells_by_genes"] = "genes_by_cells",
    barcode_col: ColumnSelector = 0,
    coordinate_cols: Sequence[ColumnSelector] | None = None,
    reorder_coordinates: bool = True,
) -> None:
    """Validate spatial count and coordinate CSV inputs.

    The count CSV is assumed to have an identifier column in the first position:
    gene IDs for ``genes_by_cells`` input or cell IDs for ``cells_by_genes``
    input.
    """

    counts_path = require_file(counts_path, label="Spatial counts file")
    coordinates_path = require_file(coordinates_path, label="Spatial coordinates file")

    if counts_orientation not in {"genes_by_cells", "cells_by_genes"}:
        raise ValidationError("counts_orientation must be 'genes_by_cells' or 'cells_by_genes'.")

    counts = pd.read_csv(counts_path, index_col=0)
    if counts.shape[0] == 0 or counts.shape[1] == 0:
        raise ValidationError("Spatial count matrix must contain at least one row and one column.")

    numeric_counts = counts.apply(pd.to_numeric, errors="coerce")
    if numeric_counts.isna().any().any():
        raise ValidationError("Spatial count matrix contains non-numeric values.")
    if (numeric_counts < 0).any().any():
        raise ValidationError("Spatial count matrix contains negative values.")

    if counts_orientation == "genes_by_cells":
        cell_ids = pd.Index(counts.columns.astype(str))
    else:
        cell_ids = pd.Index(counts.index.astype(str))
    if cell_ids.empty:
        raise ValidationError("Spatial count matrix contains no cells.")
    _warn_if_duplicate_names(cell_ids, label="spatial count cell IDs")

    coordinates = pd.read_csv(coordinates_path)
    if coordinates.shape[0] == 0:
        raise ValidationError("Spatial coordinate file contains no rows.")

    raw_barcode_series = _select_column(coordinates, barcode_col, label="barcode_col")
    if raw_barcode_series.isna().any() or (raw_barcode_series.astype(str) == "").any():
        raise ValidationError("Spatial coordinate barcode column contains missing or empty values.")
    barcode_series = raw_barcode_series.astype(str)
    _warn_if_duplicate_names(pd.Index(barcode_series), label="spatial coordinate barcodes")

    resolved_coordinate_cols = _resolve_coordinate_cols(coordinates, barcode_col, coordinate_cols)
    if len(resolved_coordinate_cols) < 2:
        raise ValidationError("Spatial coordinates must contain at least two coordinate columns.")

    coord_values = coordinates.loc[:, resolved_coordinate_cols].apply(pd.to_numeric, errors="coerce")
    if coord_values.isna().any().any():
        raise ValidationError("Spatial coordinate columns must be numeric.")
    if not np.isfinite(coord_values.to_numpy(dtype=float)).all():
        raise ValidationError("Spatial coordinate columns contain non-finite values.")

    coordinate_barcodes = pd.Index(barcode_series.astype(str))
    missing = cell_ids.difference(coordinate_barcodes)
    extra = coordinate_barcodes.difference(cell_ids)
    if len(missing) or len(extra):
        message_parts = []
        if len(missing):
            message_parts.append(f"missing coordinates for {len(missing)} count cells")
        if len(extra):
            message_parts.append(f"{len(extra)} coordinate barcodes are not present in counts")
        raise ValidationError("Spatial barcode mismatch: " + "; ".join(message_parts) + ".")

    if not reorder_coordinates and not np.array_equal(cell_ids.to_numpy(), coordinate_barcodes.to_numpy()):
        raise ValidationError(
            "Spatial coordinates are not in the same order as count matrix cells. "
            "Use reorder_coordinates=True to reorder by barcode."
        )


def validate_reference_adata(
    adata: AnnData,
    *,
    ref_label_key: str = "cluster",
    layer: str | None = None,
    require_label: bool = True,
    check_counts: bool = True,
    warn_non_integer: bool = True,
) -> None:
    """Validate a raw reference AnnData object for the build step."""

    _validate_basic_adata(adata, label="Reference")
    if require_label and ref_label_key not in adata.obs:
        raise ValidationError(f"Reference .obs is missing required label column: {ref_label_key!r}.")
    if check_counts:
        matrix = _get_matrix(adata, layer=layer, label="Reference")
        _validate_count_matrix(matrix, label="Reference", warn_non_integer=warn_non_integer)


def validate_spatial_adata(
    adata: AnnData,
    *,
    spatial_key: str = "spatial",
    layer: str | None = None,
    check_counts: bool = True,
    warn_non_integer: bool = True,
) -> None:
    """Validate a raw spatial/query AnnData object for the build step."""

    _validate_basic_adata(adata, label="Spatial")
    if spatial_key not in adata.obsm:
        raise ValidationError(f"Spatial AnnData is missing coordinates in .obsm[{spatial_key!r}].")
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2:
        raise ValidationError(f"Spatial coordinates .obsm[{spatial_key!r}] must be a 2D matrix.")
    if coords.shape[0] != adata.n_obs:
        raise ValidationError(
            f"Spatial coordinates row count ({coords.shape[0]}) does not match number of cells ({adata.n_obs})."
        )
    if coords.shape[1] < 2:
        raise ValidationError("Spatial coordinates must have at least two columns.")
    if not np.issubdtype(coords.dtype, np.number):
        raise ValidationError("Spatial coordinates must be numeric.")
    if not np.isfinite(coords).all():
        raise ValidationError("Spatial coordinates contain non-finite values.")
    if check_counts:
        matrix = _get_matrix(adata, layer=layer, label="Spatial")
        _validate_count_matrix(matrix, label="Spatial", warn_non_integer=warn_non_integer)


def validate_joint_adata(
    reference: AnnData,
    spatial: AnnData,
    *,
    min_shared_genes: int = 1,
    gene_space: Literal["shared", "reference_all"] = "shared",
) -> None:
    """Validate assumptions involving both modalities."""

    if gene_space not in {"shared", "reference_all"}:
        raise ValidationError("gene_space must be 'shared' or 'reference_all'.")
    shared = reference.var_names.intersection(spatial.var_names)
    if len(shared) < min_shared_genes:
        raise ValidationError(
            f"Only {len(shared)} shared genes found between reference and spatial data; "
            f"minimum required is {min_shared_genes}."
        )


def validate_normalization_config(config: NormalizationConfig) -> None:
    """Validate a normalization configuration object."""

    if isinstance(config, NiCoSCTransformConfig):
        if config.min_cells < 1:
            raise ValidationError("SCTransform min_cells must be >= 1.")
        if config.gmean_eps <= 0:
            raise ValidationError("SCTransform gmean_eps must be > 0.")
        if config.n_genes is not None and config.n_genes <= 0:
            raise ValidationError("SCTransform n_genes must be positive or None.")
        if config.n_cells is not None and config.n_cells <= 0:
            raise ValidationError("SCTransform n_cells must be positive or None.")
        if config.bin_size <= 0:
            raise ValidationError("SCTransform bin_size must be > 0.")
        if config.bw_adjust <= 0:
            raise ValidationError("SCTransform bw_adjust must be > 0.")
        return

    if isinstance(config, PearsonResidualsConfig):
        if not math.isfinite(config.theta) and config.theta != float("inf"):
            raise ValidationError("Pearson theta must be finite or infinity.")
        if config.theta <= 0:
            raise ValidationError("Pearson theta must be > 0.")
        if config.clip is not None and config.clip < 0:
            raise ValidationError("Pearson clip must be >= 0 or None.")
        return

    raise ValidationError(f"Unsupported normalization config type: {type(config).__name__}.")


def validate_normalization_layers(
    reference: AnnData,
    spatial: AnnData,
    config: NormalizationConfig,
) -> None:
    """Validate normalization-specific data-layer requirements."""

    if isinstance(config, PearsonResidualsConfig) and config.layer is not None:
        if config.layer not in reference.layers:
            raise ValidationError(f"Reference is missing Pearson layer {config.layer!r}.")
        if config.layer not in spatial.layers:
            raise ValidationError(f"Spatial data is missing Pearson layer {config.layer!r}.")


def validate_no_nan_or_inf(adata: AnnData, *, label: str, layer: str | None = None) -> None:
    """Validate that a processed matrix contains no NaN or infinite values."""

    matrix = _get_matrix(adata, layer=layer, label=label)
    data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if np.isnan(data).any():
        raise ValidationError(f"{label} matrix contains NaN values.")
    if not np.isfinite(data).all():
        raise ValidationError(f"{label} matrix contains infinite values.")


def validate_leiden_keys(adata: AnnData, resolutions: Sequence[float]) -> None:
    """Validate that requested Leiden clustering keys exist in ``adata.obs``."""

    missing = [f"leiden{resolution}" for resolution in resolutions if f"leiden{resolution}" not in adata.obs]
    if missing:
        raise ValidationError(f"Spatial AnnData is missing Leiden columns: {', '.join(missing)}.")


def validate_expected_files_exist(paths: Sequence[str | Path]) -> None:
    """Validate that expected output files were created."""

    for path in paths:
        require_file(path, label="Expected output file")


def _validate_basic_adata(adata: AnnData, *, label: str) -> None:
    if adata.n_obs == 0:
        raise ValidationError(f"{label} AnnData contains no cells.")
    if adata.n_vars == 0:
        raise ValidationError(f"{label} AnnData contains no genes.")
    if not adata.obs_names.is_unique:
        warnings.warn(f"{label} cell names are not unique; implementation should make them unique.", UserWarning)
    if not adata.var_names.is_unique:
        warnings.warn(f"{label} gene names are not unique; implementation should make them unique.", UserWarning)


def _get_matrix(adata: AnnData, *, layer: str | None, label: str):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise ValidationError(f"{label} AnnData is missing layer {layer!r}.")
    return adata.layers[layer]


def _validate_count_matrix(matrix, *, label: str, warn_non_integer: bool) -> None:
    if matrix is None:
        raise ValidationError(f"{label} count matrix is missing.")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValidationError(f"{label} count matrix is empty.")

    data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if data.size == 0:
        return
    if np.isnan(data).any():
        raise ValidationError(f"{label} count matrix contains NaN values.")
    if not np.isfinite(data).all():
        raise ValidationError(f"{label} count matrix contains infinite values.")
    if (data < 0).any():
        raise ValidationError(f"{label} count matrix contains negative values.")
    if warn_non_integer and not _is_integer_like(data):
        warnings.warn(
            f"{label} count matrix contains non-integer values; normalization expects raw counts.",
            UserWarning,
        )


def _infer_name_column_values(frame: pd.DataFrame, *, label: str) -> pd.Index:
    if frame.shape[1] == 0:
        raise ValidationError(f"Cannot infer {label} names from an empty table.")
    column_position = 1 if frame.shape[1] > 1 else 0
    raw_values = frame.iloc[:, column_position]
    if raw_values.isna().any() or (raw_values.astype(str) == "").any():
        raise ValidationError(f"Inferred {label} names contain missing or empty values.")
    return pd.Index(raw_values.astype(str))


def _warn_if_duplicate_names(values: pd.Index | Sequence[object], *, label: str) -> None:
    index = pd.Index(values)
    if not index.is_unique:
        duplicated_count = int(index.duplicated().sum())
        warnings.warn(f"Duplicate {label} found ({duplicated_count} duplicate rows).", UserWarning)


def _select_column(frame: pd.DataFrame, selector: ColumnSelector, *, label: str) -> pd.Series:
    if isinstance(selector, int):
        if selector < 0 or selector >= frame.shape[1]:
            raise ValidationError(f"{label} index {selector} is out of bounds for {frame.shape[1]} columns.")
        return frame.iloc[:, selector]
    if selector not in frame.columns:
        raise ValidationError(f"{label} {selector!r} not found in columns: {list(frame.columns)!r}.")
    return frame.loc[:, selector]


def _resolve_coordinate_cols(
    frame: pd.DataFrame,
    barcode_col: ColumnSelector,
    coordinate_cols: Sequence[ColumnSelector] | None,
) -> list[str]:
    if coordinate_cols is None:
        if isinstance(barcode_col, int):
            if barcode_col < 0 or barcode_col >= frame.shape[1]:
                raise ValidationError(f"barcode_col index {barcode_col} is out of bounds for {frame.shape[1]} columns.")
            barcode_name = frame.columns[barcode_col]
        else:
            barcode_name = barcode_col
        return [col for col in frame.columns if col != barcode_name]

    resolved: list[str] = []
    for col in coordinate_cols:
        if isinstance(col, int):
            if col < 0 or col >= frame.shape[1]:
                raise ValidationError(f"coordinate column index {col} is out of bounds for {frame.shape[1]} columns.")
            resolved.append(str(frame.columns[col]))
        else:
            if col not in frame.columns:
                raise ValidationError(f"coordinate column {col!r} not found in columns: {list(frame.columns)!r}.")
            resolved.append(col)
    return resolved


def _is_integer_like(values: np.ndarray) -> bool:
    values = np.asarray(values)
    if values.size == 0:
        return True
    if not np.issubdtype(values.dtype, np.number):
        return False
    finite_values = values[np.isfinite(values)]
    return bool(np.all(np.equal(finite_values, np.floor(finite_values))))
