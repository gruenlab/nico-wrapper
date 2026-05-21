"""Validation helpers for NiCo covariation inputs and outputs."""

from __future__ import annotations

from pathlib import Path

from anndata import AnnData, read_h5ad
import numpy as np
import pandas as pd

from nico_wrapper.niche.io import load_used_cell_types, normalize_radius
from nico_wrapper.niche.pipeline import load_niche_result
from nico_wrapper.niche.results import NicheInteractionResult
from nico_wrapper.niche.validation import validate_covariation_artifacts as validate_niche_covariation_artifacts

from .config import CovariationConfig
from .io import (
    load_feature_matrix,
    planned_covariation_paths,
    required_covariation_output_paths,
    resolve_ligand_receptor_db,
)


class ValidationError(ValueError):
    """Raised when covariation inputs violate wrapper assumptions."""


def validate_covariation_config(config: CovariationConfig) -> None:
    """Validate scalar covariation configuration values."""

    try:
        normalize_radius(config.radius)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    if config.n_factors < 1:
        raise ValidationError("n_factors must be >= 1.")
    if config.modality not in {"double", "single"}:
        raise ValidationError("modality must be 'double' or 'single'.")
    if config.factorization not in {"inmf", "nmf-transfer"}:
        raise ValidationError("factorization must be 'inmf' or 'nmf-transfer'.")
    _validate_nonempty_string(config.ref_label_key, label="ref_label_key")
    for label, value in (
        ("annotated_h5ad_name", config.annotated_h5ad_name),
        ("ref_original_counts_name", config.ref_original_counts_name),
        ("ref_sct_name", config.ref_sct_name),
        ("spatial_sct_name", config.spatial_sct_name),
    ):
        _validate_plain_filename(value, label=label)
    if config.ligand_receptor_db is not None:
        require_file(config.ligand_receptor_db, label="Ligand-receptor DB")
    if not config.ridge_alphas:
        raise ValidationError("ridge_alphas must not be empty.")
    for value in config.ridge_alphas:
        if not np.isfinite(value) or value <= 0:
            raise ValidationError("all ridge_alphas must be positive finite numbers.")
    for label, value in (
        ("logistic_coef_cutoff", config.logistic_coef_cutoff),
        ("ridge_coef_cutoff", config.ridge_coef_cutoff),
        ("expression_population_cutoff", config.expression_population_cutoff),
        ("shap_cluster_cutoff", config.shap_cluster_cutoff),
    ):
        if not np.isfinite(value):
            raise ValidationError(f"{label} must be finite.")
    if not isinstance(config.seed, int):
        raise ValidationError("seed must be an integer.")
    if config.shap_analysis:
        raise ValidationError("shap_analysis=True is not supported by the wrapper yet because upstream NiCo's shap import is disabled.")


def resolve_niche_context(
    *,
    niche_result: NicheInteractionResult | None,
    output_dir: str | Path | None,
    radius: int | float | str,
) -> NicheInteractionResult:
    """Resolve and validate the niche result used as covariation input."""

    radius_value, _ = normalize_radius(radius)
    if niche_result is None and output_dir is None:
        raise ValidationError("Provide either niche_result or output_dir.")
    if niche_result is not None:
        if output_dir is not None and Path(output_dir).resolve() != niche_result.output_dir.resolve():
            raise ValidationError(f"output_dir {Path(output_dir)} does not match niche_result.output_dir {niche_result.output_dir}.")
        niche_radius, _ = normalize_radius(niche_result.radius)
        if niche_radius != radius_value:
            raise ValidationError(f"config radius {radius_value!r} does not match niche_result radius {niche_radius!r}.")
        resolved = niche_result
    else:
        if output_dir is None:
            raise ValidationError("Provide either niche_result or output_dir.")
        try:
            resolved = load_niche_result(output_dir, radius=radius_value)
        except Exception as exc:
            raise ValidationError(f"Could not load niche result for covariation: {exc}") from exc

    try:
        validate_niche_covariation_artifacts(resolved.output_dir, radius=resolved.radius)
    except Exception as exc:
        raise ValidationError(f"Invalid niche artifacts for covariation: {exc}") from exc
    return resolved


def validate_covariation_inputs(
    *,
    niche_result: NicheInteractionResult | None = None,
    output_dir: str | Path | None = None,
    ref_dir: str | Path | None = None,
    spatial_dir: str | Path | None = None,
    config: CovariationConfig = CovariationConfig(),
) -> None:
    """Validate all prerequisites for a covariation run."""

    validate_covariation_config(config)
    resolved_niche = resolve_niche_context(niche_result=niche_result, output_dir=output_dir, radius=config.radius)
    paths = planned_covariation_paths(resolved_niche.output_dir, radius=config.radius, n_factors=config.n_factors)
    validate_output_paths(paths, config=config)
    try:
        resolve_ligand_receptor_db(config.ligand_receptor_db)
    except FileNotFoundError as exc:
        raise ValidationError(str(exc)) from exc

    if config.modality == "double":
        if ref_dir is None:
            raise ValidationError("ref_dir is required when modality='double'.")
        if spatial_dir is None:
            raise ValidationError("spatial_dir is required when modality='double'.")
        validate_double_modality_inputs(ref_dir=ref_dir, spatial_dir=spatial_dir, niche_result=resolved_niche, config=config)
    else:
        validate_single_modality_inputs(output_dir=resolved_niche.output_dir, config=config)


def validate_double_modality_inputs(
    *,
    ref_dir: str | Path,
    spatial_dir: str | Path,
    niche_result: NicheInteractionResult,
    config: CovariationConfig,
) -> None:
    """Validate reference and spatial expression files for double modality."""

    ref_path = Path(ref_dir)
    spatial_path = Path(spatial_dir)
    ref_original = require_file(ref_path / config.ref_original_counts_name, label="Original reference AnnData")
    ref_sct = require_file(ref_path / config.ref_sct_name, label="Normalized reference AnnData")
    spatial_sct = require_file(spatial_path / config.spatial_sct_name, label="Normalized spatial AnnData")

    ref_original_adata = _read_h5ad(ref_original, label="Original reference AnnData")
    ref_sct_adata = _read_h5ad(ref_sct, label="Normalized reference AnnData")
    spatial_sct_adata = _read_h5ad(spatial_sct, label="Normalized spatial AnnData")
    try:
        _validate_nonempty_adata(ref_original_adata, label="Original reference AnnData")
        _validate_nonempty_adata(ref_sct_adata, label="Normalized reference AnnData")
        _validate_nonempty_adata(spatial_sct_adata, label="Normalized spatial AnnData")
        if config.ref_label_key not in ref_original_adata.obs:
            raise ValidationError(f"Original reference .obs is missing label column: {config.ref_label_key!r}.")
        if ref_original_adata.raw is None:
            raise ValidationError("Original reference AnnData must have .raw set for NiCo covariation.")
        if ref_sct_adata.raw is None:
            raise ValidationError("Normalized reference AnnData must have .raw set for NiCo covariation.")
        if spatial_sct_adata.raw is None:
            raise ValidationError("Normalized spatial AnnData must have .raw set for NiCo covariation.")
        shared_genes = set(ref_sct_adata.raw.var_names).intersection(map(str, spatial_sct_adata.raw.var_names))
        if len(shared_genes) < config.n_factors:
            raise ValidationError(f"Shared raw gene count ({len(shared_genes)}) must be at least n_factors={config.n_factors}.")

        missing_ref_cells = [name for name in ref_sct_adata.obs_names if name not in ref_original_adata.obs_names]
        if missing_ref_cells:
            preview = ", ".join(map(str, missing_ref_cells[:5]))
            raise ValidationError(f"Original reference counts are missing normalized reference cells: {preview}")
        ref_labels = ref_original_adata.obs[config.ref_label_key].reindex(ref_sct_adata.obs_names).astype(str)
        required_names = set((niche_result.cell_type_names or load_used_cell_types(niche_result.used_cell_types_tsv)).values())
        missing = sorted(required_names.difference(set(ref_labels)))
        if missing:
            raise ValidationError("Spatial cell types are missing from normalized reference labels: " + ", ".join(missing))
        small_ref = ref_labels.value_counts().reindex(sorted(required_names), fill_value=0)
        small_ref = small_ref[small_ref < config.n_factors]
        if len(small_ref):
            details = ", ".join(f"{name}={int(count)}" for name, count in small_ref.items())
            raise ValidationError(f"Reference cell types need at least n_factors={config.n_factors} cells: {details}.")

        spatial_counts = _spatial_cluster_counts(niche_result.used_clusters_csv)
        small_spatial = {ct_id: count for ct_id, count in spatial_counts.items() if count < config.n_factors}
        if small_spatial:
            names = niche_result.cell_type_names or load_used_cell_types(niche_result.used_cell_types_tsv)
            details = ", ".join(f"{names.get(ct_id, str(ct_id))}={count}" for ct_id, count in sorted(small_spatial.items()))
            raise ValidationError(f"Spatial cell types need at least n_factors={config.n_factors} cells: {details}.")
    finally:
        _close_if_backed(ref_original_adata)
        _close_if_backed(ref_sct_adata)
        _close_if_backed(spatial_sct_adata)


def validate_single_modality_inputs(*, output_dir: str | Path, config: CovariationConfig) -> None:
    """Validate annotated spatial AnnData for single modality."""

    path = require_file(Path(output_dir) / config.annotated_h5ad_name, label="Annotated spatial AnnData")
    adata = _read_h5ad(path, label="Annotated spatial AnnData")
    try:
        _validate_nonempty_adata(adata, label="Annotated spatial AnnData")
        if config.ref_label_key not in adata.obs:
            raise ValidationError(f"Annotated spatial .obs is missing label column: {config.ref_label_key!r}.")
        if adata.raw is None:
            raise ValidationError("Annotated spatial AnnData must have .raw set for single-modality covariation.")
    finally:
        _close_if_backed(adata)


def validate_output_paths(paths, *, config: CovariationConfig) -> None:
    """Validate planned covariation outputs against overwrite policy."""

    planned = list(required_covariation_output_paths(paths))
    if config.persist_state:
        planned.append(paths.state_pickle)
    if config.write_manifest:
        planned.append(paths.manifest_json)
    if config.export_regression_table:
        planned.append(paths.regression_tsv)
    for path in planned:
        if path.exists() and not config.overwrite:
            raise ValidationError(f"Output file already exists: {path}. Use overwrite=True to replace it.")
        parent = path.parent if path.parent != Path("") else Path(".")
        if parent.exists() and not parent.is_dir():
            raise ValidationError(f"Output parent exists but is not a directory: {parent}")


def validate_covariation_artifacts(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
) -> bool:
    """Validate that completed covariation artifacts exist and are readable."""

    paths = planned_covariation_paths(output_dir, radius=radius, n_factors=n_factors)
    return validate_covariation_outputs(paths)


def validate_covariation_outputs(paths) -> bool:
    """Validate completed covariation core outputs."""

    missing: list[Path] = []
    empty: list[Path] = []
    for path in required_covariation_output_paths(paths):
        if not path.exists() or not path.is_file():
            missing.append(path)
        elif path.stat().st_size == 0:
            empty.append(path)
    if missing:
        raise ValidationError("Missing covariation artifacts: " + ", ".join(map(str, missing)))
    if empty:
        raise ValidationError("Empty covariation artifacts: " + ", ".join(map(str, empty)))
    if not paths.regression_dir.exists() or not paths.regression_dir.is_dir():
        raise ValidationError(f"Regression output directory is missing: {paths.regression_dir}")
    try:
        load_feature_matrix(paths.feature_matrix_npz)
    except Exception as exc:
        raise ValidationError(f"Could not read covariation feature matrix: {exc}") from exc
    return True


def require_file(path: str | Path, *, label: str = "file") -> Path:
    """Validate that ``path`` exists and is a regular file."""

    resolved = Path(path)
    if not resolved.exists():
        raise ValidationError(f"{label} does not exist: {resolved}")
    if not resolved.is_file():
        raise ValidationError(f"{label} is not a file: {resolved}")
    return resolved


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
        return read_h5ad(path, backed="r")
    except Exception as exc:  # pragma: no cover - h5ad backend-specific
        raise ValidationError(f"Could not read {label} at {path}: {exc}") from exc


def _validate_nonempty_adata(adata: AnnData, *, label: str) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValidationError(f"{label} must contain nonzero cells and genes.")
    if not adata.obs_names.is_unique:
        raise ValidationError(f"{label} .obs_names must be unique.")
    if not adata.var_names.is_unique:
        raise ValidationError(f"{label} .var_names must be unique.")


def _spatial_cluster_counts(path: Path) -> dict[int, int]:
    table = pd.read_csv(path)
    if table.shape[1] < 2:
        raise ValidationError(f"used_Clusters file must contain at least two columns: {path}")
    values = pd.to_numeric(table.iloc[:, 1], errors="raise").astype(int)
    return {int(key): int(value) for key, value in values.value_counts().items()}


def _close_if_backed(adata: AnnData) -> None:
    file_obj = getattr(adata, "file", None)
    close = getattr(file_obj, "close", None)
    if close is not None:
        close()
