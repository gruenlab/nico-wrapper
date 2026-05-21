"""Validation helpers for NiCo label-transfer inputs and outputs.

This module validates the assumptions needed by NiCo's annotation functions. It
intentionally does not run anchor discovery, perform label transfer, mutate
AnnData objects, create directories, or write output files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import warnings

from anndata import AnnData, read_h5ad
import numpy as np
from scipy import sparse

from .config import LabelTransferConfig, TieStrategy


class ValidationError(ValueError):
    """Raised when NiCo label-transfer inputs violate required assumptions."""


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
        Output files that will be written by the transfer workflow.
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


def validate_label_transfer_config(config: LabelTransferConfig) -> None:
    """Validate scalar label-transfer configuration values.

    This checks values that are independent of the input AnnData sizes. Shape-
    dependent checks, such as whether ``n_pcs`` is feasible, are performed by
    :func:`validate_label_transfer_inputs` after loading the AnnData files.
    """

    anchors = config.anchors
    annotation = config.annotation

    if anchors.neighbors < 2:
        raise ValidationError("neighbors must be >= 2 because NiCo expects iterable KNN results.")
    if anchors.n_pcs <= 0:
        raise ValidationError("n_pcs must be > 0.")
    if anchors.minkowski_order < 1:
        raise ValidationError("minkowski_order must be >= 1.")

    _validate_plain_filename(anchors.spatial_sct_filename, label="spatial_sct_filename")
    _validate_plain_filename(anchors.sc_sct_filename, label="sc_sct_filename")
    _validate_plain_filename(anchors.sc_full_filename, label="sc_full_filename")
    _validate_plain_filename(annotation.output_h5ad_name, label="output_h5ad_name")

    _validate_nonempty_key(annotation.ref_label_key, label="ref_label_key")
    _validate_nonempty_key(annotation.spatial_cluster_key, label="spatial_cluster_key")
    _validate_nonempty_key(annotation.output_label_key, label="output_label_key")

    if not 0 <= annotation.dispersion_cutoff <= 1:
        raise ValidationError("dispersion_cutoff must be between 0 and 1.")
    if annotation.iterations < 1:
        raise ValidationError("iterations must be >= 1.")
    tie_strategy = annotation.tie_strategy.value if isinstance(annotation.tie_strategy, TieStrategy) else annotation.tie_strategy
    if tie_strategy not in {TieStrategy.MAJORITY.value, TieStrategy.WEIGHTED.value}:
        raise ValidationError("tie_strategy must be 'majority' or 'weighted'.")


def validate_label_transfer_inputs(
    ref_dir: str | Path,
    spatial_dir: str | Path,
    *,
    output_dir: str | Path,
    annotation_dir: str | Path | None = None,
    config: LabelTransferConfig = LabelTransferConfig(),
    warn_dense_memory_gb: float | None = 8.0,
) -> None:
    """Validate files and AnnData schemas for the complete transfer workflow.

    Parameters
    ----------
    ref_dir
        Directory containing the preprocessed reference files. The file names
        are taken from ``config.anchors``.
    spatial_dir
        Directory containing the preprocessed spatial/query file. The file name
        is taken from ``config.anchors``.
    output_dir
        Planned output directory for the annotated AnnData and NiCo artifacts.
    annotation_dir
        Planned directory for annotation intermediates. If omitted, validation
        assumes ``output_dir / "annotations"``.
    config
        Label-transfer configuration to validate against the data.
    warn_dense_memory_gb
        If not ``None``, warn when NiCo's known dense matrix conversions for the
        shared gene space are estimated to require at least this many GiB.

    Raises
    ------
    ValidationError
        If any required file, key, shape, output path, or cross-file invariant is
        invalid.

    Notes
    -----
    ``n_pcs`` is intentionally a hard error when too large. The transfer space
    should not be silently clamped because that changes annotation behavior.
    """

    validate_label_transfer_config(config)

    ref_dir = Path(ref_dir)
    spatial_dir = Path(spatial_dir)
    output_dir = Path(output_dir)
    resolved_annotation_dir = Path(annotation_dir) if annotation_dir is not None else output_dir / "annotations"

    _validate_directory_path(ref_dir, label="ref_dir", must_exist=True)
    _validate_directory_path(spatial_dir, label="spatial_dir", must_exist=True)
    _validate_directory_path(output_dir, label="output_dir", must_exist=False)
    _validate_directory_path(resolved_annotation_dir, label="annotation_dir", must_exist=False)

    anchors = config.anchors
    annotation = config.annotation

    sc_full_path = require_file(ref_dir / anchors.sc_full_filename, label="Full/original reference AnnData")
    sc_sct_path = require_file(ref_dir / anchors.sc_sct_filename, label="Normalized reference AnnData")
    spatial_sct_path = require_file(spatial_dir / anchors.spatial_sct_filename, label="Normalized spatial AnnData")

    planned_outputs = _planned_output_files(
        output_dir=output_dir,
        annotation_dir=resolved_annotation_dir,
        neighbors=anchors.neighbors,
        iterations=annotation.iterations,
        output_h5ad_name=annotation.output_h5ad_name,
    )
    validate_output_files(planned_outputs, overwrite=config.overwrite)

    sc_full = _read_h5ad(sc_full_path, label="Full/original reference AnnData")
    sc_sct = _read_h5ad(sc_sct_path, label="Normalized reference AnnData")
    spatial_sct = _read_h5ad(spatial_sct_path, label="Normalized spatial AnnData")

    validate_full_reference_adata(sc_full, ref_label_key=annotation.ref_label_key)
    validate_normalized_reference_adata(sc_sct)
    validate_normalized_spatial_adata(
        spatial_sct,
        spatial_cluster_key=annotation.spatial_cluster_key,
        output_label_key=annotation.output_label_key,
        overwrite=config.overwrite,
    )
    validate_joint_transfer_adata(
        sc_full=sc_full,
        sc_sct=sc_sct,
        spatial_sct=spatial_sct,
        neighbors=anchors.neighbors,
        n_pcs=anchors.n_pcs,
        warn_dense_memory_gb=warn_dense_memory_gb,
    )


def validate_full_reference_adata(adata: AnnData, *, ref_label_key: str = "cluster") -> None:
    """Validate ``Original_counts.h5ad`` for label transfer.

    The reference labels must be present and safe for NiCo's manual CSV-writing
    and internal sentinel values.
    """

    _validate_basic_adata(adata, label="Full/original reference")
    _validate_matrix(adata.X, label="Full/original reference .X")

    if ref_label_key not in adata.obs:
        raise ValidationError(f"Full/original reference .obs is missing label column: {ref_label_key!r}.")

    labels = adata.obs[ref_label_key]
    if labels.isna().any():
        raise ValidationError(f"Reference label column {ref_label_key!r} contains missing values.")

    label_values = labels.astype(str)
    if (label_values.str.len() == 0).any():
        raise ValidationError(f"Reference label column {ref_label_key!r} contains empty labels.")

    reserved = {"NM", "xxxx"}
    observed_reserved = sorted(reserved.intersection(set(label_values)))
    if observed_reserved:
        raise ValidationError(
            f"Reference label column {ref_label_key!r} contains NiCo-reserved labels: {observed_reserved}."
        )

    if label_values.str.contains("_a#d_", regex=False).any():
        raise ValidationError(f"Reference labels may not contain NiCo's internal delimiter '_a#d_'.")
    if label_values.str.contains(r"[,\r\n]", regex=True).any():
        raise ValidationError(
            f"Reference labels may not contain commas or newlines because NiCo writes annotation CSVs manually."
        )


def validate_normalized_reference_adata(adata: AnnData) -> None:
    """Validate the normalized reference AnnData used for anchor discovery."""

    _validate_basic_adata(adata, label="Normalized reference")
    _validate_matrix(adata.X, label="Normalized reference .X")
    if not _matrix_has_nonzero(adata.X):
        raise ValidationError("Normalized reference .X contains no nonzero values.")


def validate_normalized_spatial_adata(
    adata: AnnData,
    *,
    spatial_cluster_key: str = "leiden0.5",
    output_label_key: str = "nico_ct",
    overwrite: bool = False,
) -> None:
    """Validate the normalized spatial/query AnnData used for transfer."""

    _validate_basic_adata(adata, label="Normalized spatial")
    _validate_matrix(adata.X, label="Normalized spatial .X")
    if not _matrix_has_nonzero(adata.X):
        raise ValidationError("Normalized spatial .X contains no nonzero values.")

    if spatial_cluster_key not in adata.obs:
        raise ValidationError(f"Normalized spatial .obs is missing guide cluster column: {spatial_cluster_key!r}.")

    clusters = adata.obs[spatial_cluster_key]
    if clusters.isna().any():
        raise ValidationError(f"Spatial guide cluster column {spatial_cluster_key!r} contains missing values.")
    try:
        sorted(list(np.unique(clusters.to_numpy())))
    except TypeError as exc:
        raise ValidationError(
            f"Spatial guide cluster column {spatial_cluster_key!r} contains values that NiCo cannot sort."
        ) from exc

    if output_label_key in adata.obs and not overwrite:
        raise ValidationError(
            f"Output label key {output_label_key!r} already exists in normalized spatial .obs. "
            "Use overwrite=True to replace it."
        )

    if "spatial" in adata.obsm:
        _validate_spatial_coordinates(adata, spatial_key="spatial")


def validate_joint_transfer_adata(
    *,
    sc_full: AnnData,
    sc_sct: AnnData,
    spatial_sct: AnnData,
    neighbors: int,
    n_pcs: int,
    warn_dense_memory_gb: float | None = 8.0,
) -> None:
    """Validate cross-file invariants needed by NiCo label transfer."""

    if neighbors > spatial_sct.n_obs:
        raise ValidationError(
            f"neighbors ({neighbors}) must be <= number of normalized spatial cells ({spatial_sct.n_obs})."
        )
    if neighbors > sc_sct.n_obs:
        raise ValidationError(
            f"neighbors ({neighbors}) must be <= number of normalized reference cells ({sc_sct.n_obs})."
        )

    missing_reference_cells = sc_sct.obs_names.difference(sc_full.obs_names)
    if len(missing_reference_cells):
        raise ValidationError(
            "Normalized reference contains cells absent from full/original reference AnnData: "
            f"{len(missing_reference_cells)} missing."
        )

    shared_genes = sc_sct.var_names.intersection(spatial_sct.var_names)
    if len(shared_genes) == 0:
        raise ValidationError("Normalized reference and spatial AnnData files share no genes.")

    max_pcs = min(sc_sct.n_obs, len(shared_genes)) - 1
    if n_pcs > max_pcs:
        raise ValidationError(
            "n_pcs is too large for PCA on the normalized reference shared-gene matrix: "
            f"requested {n_pcs}, maximum allowed is {max_pcs} "
            f"(reference cells={sc_sct.n_obs}, shared genes={len(shared_genes)})."
        )

    if not _shared_gene_matrix_has_nonzero(sc_sct, shared_genes):
        raise ValidationError("Normalized reference has no nonzero values in the shared gene space.")
    if not _shared_gene_matrix_has_nonzero(spatial_sct, shared_genes):
        raise ValidationError("Normalized spatial data has no nonzero values in the shared gene space.")

    if warn_dense_memory_gb is not None:
        estimated_bytes = (sc_sct.n_obs + spatial_sct.n_obs) * len(shared_genes) * 8
        estimated_gib = estimated_bytes / 1024**3
        if estimated_gib >= warn_dense_memory_gb:
            warnings.warn(
                "NiCo densifies normalized matrices in the shared gene space during anchor discovery. "
                f"Estimated dense array memory is {estimated_gib:.2f} GiB "
                f"({sc_sct.n_obs} reference cells, {spatial_sct.n_obs} spatial cells, "
                f"{len(shared_genes)} shared genes).",
                UserWarning,
                stacklevel=2,
            )


def _planned_output_files(
    *,
    output_dir: Path,
    annotation_dir: Path,
    neighbors: int,
    iterations: int,
    output_h5ad_name: str,
) -> list[Path]:
    files = [
        output_dir / output_h5ad_name,
        annotation_dir / f"anchors_data_{neighbors}.npz",
        annotation_dir / "final_sct_sc.h5ad",
        annotation_dir / "final_sct_sp.h5ad",
    ]
    for iteration in range(1, iterations + 1):
        files.append(annotation_dir / f"{iteration}_nico_annotation_cluster.csv")
        files.append(annotation_dir / f"{iteration}_nico_annotation_ct_name.csv")
    return files


def _validate_directory_path(path: Path, *, label: str, must_exist: bool) -> None:
    if path.exists() and not path.is_dir():
        raise ValidationError(f"{label} exists but is not a directory: {path}")
    if must_exist and not path.exists():
        raise ValidationError(f"{label} does not exist: {path}")
    parent = path.parent if path.parent != Path("") else Path(".")
    if parent.exists() and not parent.is_dir():
        raise ValidationError(f"{label} parent exists but is not a directory: {parent}")


def _validate_plain_filename(value: str, *, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{label} must be a non-empty file name.")
    path = Path(value)
    if path.is_absolute() or path.name != value:
        raise ValidationError(f"{label} must be a plain file name, not a path: {value!r}.")


def _validate_nonempty_key(value: str, *, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{label} must be a non-empty string.")


def _read_h5ad(path: Path, *, label: str) -> AnnData:
    try:
        return read_h5ad(path)
    except Exception as exc:  # pragma: no cover - depends on h5ad backend failure modes
        raise ValidationError(f"Could not read {label} at {path}: {exc}") from exc


def _validate_basic_adata(adata: AnnData, *, label: str) -> None:
    if adata.n_obs == 0:
        raise ValidationError(f"{label} AnnData contains no cells.")
    if adata.n_vars == 0:
        raise ValidationError(f"{label} AnnData contains no genes.")
    if not adata.obs_names.is_unique:
        raise ValidationError(f"{label} AnnData .obs_names must be unique.")
    if not adata.var_names.is_unique:
        raise ValidationError(f"{label} AnnData .var_names must be unique.")
    if adata.obs_names.isna().any():
        raise ValidationError(f"{label} AnnData .obs_names contains missing values.")
    if adata.var_names.isna().any():
        raise ValidationError(f"{label} AnnData .var_names contains missing values.")


def _validate_matrix(matrix: object, *, label: str) -> None:
    if sparse.issparse(matrix):
        if not np.issubdtype(matrix.dtype, np.number):
            raise ValidationError(f"{label} must have a numeric dtype; found {matrix.dtype}.")
        if not np.isfinite(matrix.data).all():
            raise ValidationError(f"{label} contains NaN or infinite values.")
        return

    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValidationError(f"{label} must be a 2D matrix.")
    if not np.issubdtype(array.dtype, np.number):
        raise ValidationError(f"{label} must have a numeric dtype; found {array.dtype}.")
    if not np.isfinite(array).all():
        raise ValidationError(f"{label} contains NaN or infinite values.")


def _matrix_has_nonzero(matrix: object) -> bool:
    if sparse.issparse(matrix):
        return bool(matrix.nnz and np.any(matrix.data != 0))
    return bool(np.any(np.asarray(matrix) != 0))


def _shared_gene_matrix_has_nonzero(adata: AnnData, shared_genes: Sequence[str]) -> bool:
    indices = adata.var_names.get_indexer(shared_genes)
    if np.any(indices < 0):
        return False
    matrix = adata.X[:, indices]
    return _matrix_has_nonzero(matrix)


def _validate_spatial_coordinates(adata: AnnData, *, spatial_key: str) -> None:
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
        raise ValidationError("Spatial coordinates contain NaN or infinite values.")
