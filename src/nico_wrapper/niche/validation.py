"""Validation helpers for NiCo niche interaction inputs and outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import pickle

from anndata import AnnData, read_h5ad
import numpy as np
import pandas as pd

from .config import NicheInteractionConfig
from .io import (
    core_output_paths,
    load_classifier_matrices,
    load_used_cell_types,
    normalize_radius,
    planned_artifact_paths,
    required_covariation_paths,
    wrapper_output_paths,
)


class ValidationError(ValueError):
    """Raised when niche interaction inputs violate wrapper assumptions."""


def validate_niche_config(config: NicheInteractionConfig) -> None:
    """Validate scalar niche interaction configuration values.

    This function does not read input data. Shape- and schema-dependent checks
    are performed by :func:`validate_niche_inputs`.
    """

    _validate_plain_filename(config.anndata_filename, label="anndata_filename")
    _validate_nonempty_string(config.label_key, label="label_key")
    _validate_nonempty_string(config.spatial_key, label="spatial_key")

    neighborhood = config.neighborhood
    try:
        normalize_radius(neighborhood.radius)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    if not np.isfinite(neighborhood.epsilon_threshold) or neighborhood.epsilon_threshold <= 0:
        raise ValidationError("epsilon_threshold must be a positive finite number.")
    for cell_type in neighborhood.additional_excluded_cell_types:
        _validate_nonempty_string(cell_type, label="additional_excluded_cell_types entry")

    model = config.model
    if model.k_fold < 2:
        raise ValidationError("k_fold must be >= 2.")
    if model.n_repeats < 1:
        raise ValidationError("n_repeats must be >= 1.")
    if model.c_values is not None:
        if not model.c_values:
            raise ValidationError("c_values must not be empty when provided.")
        for value in model.c_values:
            if not np.isfinite(value) or value <= 0:
                raise ValidationError("all c_values must be positive finite numbers.")

    plots = config.plots
    if not plots.saveas.strip():
        raise ValidationError("plots.saveas must be non-empty.")
    if plots.dpi <= 0:
        raise ValidationError("plots.dpi must be > 0.")
    if not np.isfinite(plots.interaction_cutoff):
        raise ValidationError("plots.interaction_cutoff must be finite.")

    proximity = config.proximity
    if proximity.n_permutations < 1:
        raise ValidationError("proximity.n_permutations must be >= 1.")
    if not np.isfinite(proximity.observed_threshold):
        raise ValidationError("proximity.observed_threshold must be finite.")
    if not proximity.saveas.strip():
        raise ValidationError("proximity.saveas must be non-empty.")


def validate_niche_inputs(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> None:
    """Validate input AnnData, configuration, and planned output files.

    Parameters
    ----------
    output_dir
        Directory containing the annotated spatial AnnData and receiving NiCo
        niche interaction artifacts.
    config
        Niche analysis configuration.

    Raises
    ------
    ValidationError
        If a required file, AnnData key, scalar config value, or output path is
        invalid.
    """

    validate_niche_config(config)
    paths = planned_artifact_paths(
        output_dir,
        radius=config.neighborhood.radius,
        anndata_filename=config.anndata_filename,
    )

    _validate_directory_path(paths.output_dir, label="output_dir", must_exist=True)
    require_file(paths.annotated_h5ad, label="Annotated spatial AnnData")
    validate_output_files(
        [*core_output_paths(paths), *wrapper_output_paths(paths)],
        overwrite=config.overwrite,
    )

    adata = _read_h5ad(paths.annotated_h5ad, label="Annotated spatial AnnData")
    validate_annotated_spatial_adata(
        adata,
        label_key=config.label_key,
        spatial_key=config.spatial_key,
        k_fold=config.model.k_fold,
        additional_excluded_cell_types=config.neighborhood.additional_excluded_cell_types,
    )


def validate_annotated_spatial_adata(
    adata: AnnData,
    *,
    label_key: str = "nico_ct",
    spatial_key: str = "spatial",
    k_fold: int = 5,
    additional_excluded_cell_types: Sequence[str] = (),
) -> None:
    """Validate an annotated spatial AnnData object for niche analysis."""

    if adata.n_obs == 0:
        raise ValidationError("Annotated spatial AnnData contains no cells.")
    if not adata.obs_names.is_unique:
        raise ValidationError("Annotated spatial AnnData .obs_names must be unique.")
    if adata.obs_names.isna().any():
        raise ValidationError("Annotated spatial AnnData .obs_names contains missing values.")

    if label_key not in adata.obs:
        raise ValidationError(f"Annotated spatial .obs is missing label column: {label_key!r}.")
    labels = adata.obs[label_key]
    if labels.isna().any():
        raise ValidationError(f"Label column {label_key!r} contains missing values.")
    label_values = labels.astype(str)
    if (label_values.str.len() == 0).any():
        raise ValidationError(f"Label column {label_key!r} contains empty labels.")
    try:
        sorted(list(np.unique(label_values.to_numpy())))
    except TypeError as exc:
        raise ValidationError(f"Label column {label_key!r} contains values NiCo cannot sort.") from exc

    if spatial_key not in adata.obsm:
        raise ValidationError(f"Annotated spatial .obsm is missing coordinates: {spatial_key!r}.")
    _validate_spatial_coordinates(adata, spatial_key=spatial_key)

    excluded = {"NM", *additional_excluded_cell_types}
    counts = label_values.value_counts()
    retained_counts = counts[~counts.index.isin(excluded)]
    retained_counts = retained_counts[retained_counts >= 5]
    if len(retained_counts) < 3:
        raise ValidationError(
            "At least three cell types with >=5 cells must remain after excluding "
            f"{sorted(excluded)!r}; found {len(retained_counts)}. "
            "NiCo's interaction module is not robust for binary-class ROC output."
        )
    too_small_for_cv = retained_counts[retained_counts < k_fold]
    if len(too_small_for_cv):
        details = ", ".join(f"{name}={count}" for name, count in too_small_for_cv.items())
        raise ValidationError(
            f"Every retained cell type needs at least k_fold={k_fold} cells before modeling; too small: {details}."
        )


def validate_covariation_artifacts(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
) -> bool:
    """Validate that downstream NiCo covariation artifacts exist and are readable."""

    paths = planned_artifact_paths(output_dir, radius=radius, anndata_filename=anndata_filename)
    for path in required_covariation_paths(paths):
        require_file(path, label="Covariation prerequisite")
        if path.stat().st_size == 0:
            raise ValidationError(f"Covariation prerequisite is empty: {path}")
    try:
        matrices = load_classifier_matrices(paths.classifier_matrices_npz)
        cell_types = load_used_cell_types(paths.used_cell_types_tsv)
        clusters = pd.read_csv(paths.used_clusters_csv)
        with paths.neighbors_pickle.open("rb") as handle:
            neighbors = pickle.load(handle)
        with paths.distances_pickle.open("rb") as handle:
            distances = pickle.load(handle)
    except Exception as exc:
        raise ValidationError(f"Could not read covariation prerequisite artifacts: {exc}") from exc

    if len(clusters) != len(neighbors):
        raise ValidationError(
            "used_Clusters row count does not match neighbor-list length: "
            f"{len(clusters)} rows vs {len(neighbors)} neighbor entries."
        )
    if len(neighbors) != len(distances):
        raise ValidationError(
            "Neighbor and distance pickle lengths differ: "
            f"{len(neighbors)} neighbors vs {len(distances)} distance entries."
        )
    for index, (neighbor_row, distance_row) in enumerate(zip(neighbors, distances)):
        if len(neighbor_row) != len(distance_row):
            raise ValidationError(f"Neighbor and distance entry lengths differ at row {index}.")
    coef = np.asarray(matrices["coef"])
    if coef.ndim != 2:
        raise ValidationError(f"Classifier coefficient matrix must be 2D; found shape {coef.shape}.")
    if coef.shape[0] > len(cell_types) or coef.shape[1] > len(cell_types):
        raise ValidationError(
            "Classifier coefficient shape is inconsistent with used_CT.txt: "
            f"coef shape={coef.shape}, used cell types={len(cell_types)}."
        )
    return True


def validate_output_files(paths: Sequence[str | Path], *, overwrite: bool = False) -> None:
    """Validate planned output files before writing."""

    for raw_path in paths:
        path = Path(raw_path)
        if path.exists() and not overwrite:
            raise ValidationError(f"Output file already exists: {path}. Use overwrite=True to replace it.")
        parent = path.parent if path.parent != Path("") else Path(".")
        if parent.exists() and not parent.is_dir():
            raise ValidationError(f"Output parent exists but is not a directory: {parent}")


def require_file(path: str | Path, *, label: str = "file") -> Path:
    """Validate that ``path`` exists and is a regular file."""

    resolved = Path(path)
    if not resolved.exists():
        raise ValidationError(f"{label} does not exist: {resolved}")
    if not resolved.is_file():
        raise ValidationError(f"{label} is not a file: {resolved}")
    return resolved


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


def _validate_nonempty_string(value: str, *, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{label} must be a non-empty string.")


def _read_h5ad(path: Path, *, label: str) -> AnnData:
    try:
        return read_h5ad(path)
    except Exception as exc:  # pragma: no cover - h5ad backend-specific
        raise ValidationError(f"Could not read {label} at {path}: {exc}") from exc


def _validate_spatial_coordinates(adata: AnnData, *, spatial_key: str) -> None:
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2:
        raise ValidationError(f"Spatial coordinates .obsm[{spatial_key!r}] must be a 2D matrix.")
    if coords.shape[0] != adata.n_obs:
        raise ValidationError(
            f"Spatial coordinates row count ({coords.shape[0]}) does not match number of cells ({adata.n_obs})."
        )
    if coords.shape[1] not in {2, 3}:
        raise ValidationError(
            f"Spatial coordinates must have exactly 2 or 3 columns for NiCo interactions; found {coords.shape[1]}."
        )
    if not np.issubdtype(coords.dtype, np.number):
        raise ValidationError("Spatial coordinates must be numeric.")
    if not np.isfinite(coords).all():
        raise ValidationError("Spatial coordinates contain NaN or infinite values.")
