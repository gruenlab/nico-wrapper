"""Public API for NiCo covariation analysis."""

from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
import shutil
from typing import Any

from anndata import read_h5ad
import numpy as np
import pandas as pd

from nico_wrapper.niche.io import load_used_cell_types, normalize_radius
from nico_wrapper.niche.pipeline import load_niche_result
from nico_wrapper.niche.results import NicheInteractionResult

from .config import CovariationConfig
from .io import (
    planned_covariation_paths,
    read_covariation_manifest,
    read_covariation_state,
    resolve_ligand_receptor_db,
    write_covariation_manifest,
    write_covariation_state,
)
from .results import CovariationResult
from .validation import (
    ValidationError,
    resolve_niche_context,
    validate_covariation_inputs,
    validate_covariation_outputs,
)


def run_covariation(
    *,
    niche_result: NicheInteractionResult | None = None,
    output_dir: str | Path | None = None,
    ref_dir: str | Path | None = None,
    spatial_dir: str | Path | None = None,
    config: CovariationConfig = CovariationConfig(),
) -> CovariationResult:
    """Run NiCo niche covariation analysis."""

    validate_covariation_inputs(
        niche_result=niche_result,
        output_dir=output_dir,
        ref_dir=ref_dir,
        spatial_dir=spatial_dir,
        config=config,
    )
    resolved_niche = resolve_niche_context(niche_result=niche_result, output_dir=output_dir, radius=config.radius)
    paths = planned_covariation_paths(resolved_niche.output_dir, radius=config.radius, n_factors=config.n_factors)

    upstream_ref_dir = ref_dir
    upstream_spatial_dir = spatial_dir
    upstream_ref_original_counts = config.ref_original_counts_name
    if config.modality == "double":
        upstream_ref_dir, upstream_spatial_dir, upstream_ref_original_counts = _prepare_upstream_input_dirs(
            niche_result=resolved_niche,
            ref_dir=ref_dir,
            spatial_dir=spatial_dir,
            config=config,
        )

    nico_result = _run_upstream_gene_covariation_analysis(
        niche_result=resolved_niche,
        ref_dir=upstream_ref_dir,
        spatial_dir=upstream_spatial_dir,
        ref_original_counts_name=upstream_ref_original_counts,
        config=config,
    )
    _validate_nico_result(nico_result)

    result = build_covariation_result(
        resolved_niche.output_dir,
        radius=config.radius,
        n_factors=config.n_factors,
        modality=config.modality,
        niche_result=resolved_niche,
        nico_result=nico_result,
    )
    result = replace(
        result,
        state_pickle=result.state_pickle if config.persist_state else None,
        regression_tsv=result.regression_tsv if config.export_regression_table else None,
        manifest_json=None,
    )

    if config.persist_state:
        state_path = write_covariation_state(nico_result, paths.state_pickle)
        result = replace(result, state_pickle=state_path)
    if config.export_regression_table:
        table_path = export_regression_table(result, ridge_coef_cutoff=config.ridge_coef_cutoff)
        result = replace(result, regression_tsv=table_path)
    if config.write_manifest:
        manifest_path = write_covariation_manifest(result, config, paths.manifest_json)
        result = replace(result, manifest_json=manifest_path)
    return result


def build_covariation_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
    modality: str = "double",
    niche_result: NicheInteractionResult | None = None,
    nico_result: Any | None = None,
) -> CovariationResult:
    """Build a stable covariation result from artifacts and metadata."""

    paths = planned_covariation_paths(output_dir, radius=radius, n_factors=n_factors)
    validate_covariation_outputs(paths)
    cell_type_names: dict[int, str] | None = None
    if niche_result is not None:
        cell_type_names = niche_result.cell_type_names or load_used_cell_types(niche_result.used_cell_types_tsv)
    else:
        used_ct = paths.output_dir / "used_CT.txt"
        if used_ct.exists():
            cell_type_names = load_used_cell_types(used_ct)

    return CovariationResult(
        output_dir=paths.output_dir,
        covariation_dir=paths.covariation_dir,
        radius=paths.radius,
        radius_tag=paths.radius_tag,
        n_factors=paths.n_factors,
        modality=modality,
        factors_pickle=paths.factors_pickle,
        feature_matrix_npz=paths.feature_matrix_npz,
        regression_dir=paths.regression_dir,
        state_pickle=paths.state_pickle if paths.state_pickle.exists() else None,
        manifest_json=paths.manifest_json if paths.manifest_json.exists() else None,
        regression_tsv=paths.regression_tsv if paths.regression_tsv.exists() else None,
        niche_result=niche_result,
        cell_type_names=cell_type_names,
        nico_result=nico_result,
    )


def load_covariation_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
    load_state: bool = True,
) -> CovariationResult:
    """Load an existing covariation result from disk."""

    paths = planned_covariation_paths(output_dir, radius=radius, n_factors=n_factors)
    modality = "double"
    cell_type_names: dict[int, str] | None = None
    if paths.manifest_json.exists():
        manifest = read_covariation_manifest(paths.manifest_json)
        modality = str(manifest.get("modality", modality))
        raw_names = manifest.get("cell_type_names")
        if isinstance(raw_names, dict):
            cell_type_names = {int(key): str(value) for key, value in raw_names.items()}

    niche_result: NicheInteractionResult | None
    try:
        niche_result = load_niche_result(paths.output_dir, radius=paths.radius)
    except Exception:
        niche_result = None

    nico_result = None
    if load_state and paths.state_pickle.exists():
        nico_result = read_covariation_state(paths.state_pickle)

    result = build_covariation_result(
        paths.output_dir,
        radius=paths.radius,
        n_factors=paths.n_factors,
        modality=modality,
        niche_result=niche_result,
        nico_result=nico_result,
    )
    if cell_type_names is not None:
        result = replace(result, cell_type_names=cell_type_names)
    return result


def export_regression_table(
    result: CovariationResult,
    *,
    output_path: str | Path | None = None,
    ridge_coef_cutoff: float | None = None,
) -> Path:
    """Export NiCo ridge-regression coefficients to a long TSV table."""

    nico_result = _require_regression_state(result)
    save_reg_coef = getattr(nico_result, "save_reg_coef", None)
    if save_reg_coef is None:
        raise ValueError("Regression export requires covariation state containing NiCo save_reg_coef.")
    cutoff = _resolve_ridge_cutoff(nico_result, ridge_coef_cutoff)
    cell_type_names = result.cell_type_names or _cell_type_names_from_state_or_disk(result, nico_result)
    name_to_id = {name: cell_id for cell_id, name in cell_type_names.items()}

    rows: list[dict[str, object]] = []
    for central_id_raw, payload in dict(save_reg_coef).items():
        central_id = int(central_id_raw)
        if len(payload) < 10:
            raise ValueError(f"Unexpected save_reg_coef payload for cell type {central_id}: expected at least 10 entries.")
        coef = np.asarray(payload[0], dtype=float)
        alpha = np.asarray(payload[2], dtype=float)
        xlabel = [str(value) for value in np.asarray(payload[3]).tolist()]
        score = np.asarray(payload[4], dtype=float)
        pvalue = np.asarray(payload[7], dtype=float)
        pve = np.asarray(payload[8], dtype=float)
        if coef.ndim != 2:
            raise ValueError(f"Coefficient matrix for cell type {central_id} must be 2D; found {coef.shape}.")
        if pvalue.shape != coef.shape:
            raise ValueError(f"P-value matrix for cell type {central_id} must match coefficient shape {coef.shape}; found {pvalue.shape}.")
        n_factors = coef.shape[0]
        if coef.shape[1] % n_factors != 0:
            raise ValueError(f"Coefficient columns for cell type {central_id} are not divisible by n_factors={n_factors}.")
        n_neighbors = coef.shape[1] // n_factors
        largest_abs = float(np.max(np.abs(coef))) if coef.size else 0.0
        normalized_coef = coef / largest_abs if largest_abs else np.full_like(coef, np.nan, dtype=float)
        for central_factor_index in range(n_factors):
            for neighbor_index in range(n_neighbors):
                neighbor_name = xlabel[neighbor_index] if neighbor_index < len(xlabel) else str(neighbor_index)
                neighbor_id = name_to_id.get(neighbor_name)
                for neighbor_factor_index in range(n_factors):
                    col = neighbor_index * n_factors + neighbor_factor_index
                    coefficient = float(coef[central_factor_index, col])
                    normalized_coefficient = float(normalized_coef[central_factor_index, col])
                    rows.append(
                        {
                            "central_cell_type_id": central_id,
                            "central_cell_type": cell_type_names.get(central_id, str(central_id)),
                            "central_factor": central_factor_index + 1,
                            "neighbor_cell_type_id": neighbor_id,
                            "neighbor_cell_type": neighbor_name,
                            "neighbor_factor": neighbor_factor_index + 1,
                            "neighbor_block_index": neighbor_index,
                            "logistic_score": _safe_index(score, neighbor_index),
                            "ridge_coefficient": coefficient,
                            "normalized_ridge_coefficient": normalized_coefficient,
                            "abs_normalized_ridge_coefficient": abs(normalized_coefficient)
                            if np.isfinite(normalized_coefficient)
                            else float("nan"),
                            "p_value": float(pvalue[central_factor_index, col]),
                            "selected_alpha": _safe_index(alpha, central_factor_index),
                            "explained_variance": _safe_index(pve, central_factor_index),
                            "passes_ridge_cutoff": bool(np.isfinite(normalized_coefficient) and abs(normalized_coefficient) >= cutoff),
                            "radius": result.radius,
                            "n_factors": n_factors,
                        }
                    )

    destination = Path(output_path) if output_path is not None else result.covariation_dir / "regression_coefficients.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(destination, sep="\t", index=False)
    return destination


def _run_upstream_gene_covariation_analysis(
    *,
    niche_result: NicheInteractionResult,
    ref_dir: str | Path | None,
    spatial_dir: str | Path | None,
    ref_original_counts_name: str,
    config: CovariationConfig,
) -> Any:
    try:
        from nico import Covariations as scov
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo covariation analysis requires the 'nico' package.") from exc

    radius_value, _ = normalize_radius(config.radius)
    lr_db = resolve_ligand_receptor_db(config.ligand_receptor_db)
    return scov.gene_covariation_analysis(
        Radius=radius_value,
        output_niche_prediction_dir=_as_nico_dir(niche_result.output_dir),
        refpath=_as_nico_dir(ref_dir) if ref_dir is not None else "./inputRef/",
        quepath=_as_nico_dir(spatial_dir) if spatial_dir is not None else "./inputQuery/",
        ref_cluster_tag=config.ref_label_key,
        ref_original_counts=ref_original_counts_name,
        LRdbFilename=str(lr_db),
        iNMFmode=(config.factorization == "inmf"),
        no_of_factors=config.n_factors,
        shap_analysis=config.shap_analysis,
        shap_cluster_cutoff=config.shap_cluster_cutoff,
        cutoff_to_count_exp_cell_population=config.expression_population_cutoff,
        seed=config.seed,
        spatial_integration_modality=config.modality,
        anndata_object_name=config.annotated_h5ad_name,
        lambda_c=list(config.ridge_alphas),
        coeff_cutoff_for_rid_reg=config.ridge_coef_cutoff,
        logistic_coef_cutoff=config.logistic_coef_cutoff,
    )


def _prepare_upstream_input_dirs(
    *,
    niche_result: NicheInteractionResult,
    ref_dir: str | Path | None,
    spatial_dir: str | Path | None,
    config: CovariationConfig,
) -> tuple[Path, Path, str]:
    """Stage NiCo-compatible input filenames/cell sets when upstream requires it.

    NiCo hardcodes ``sct_singleCell.h5ad`` and ``sct_spatial.h5ad`` and assumes
    ``Original_counts.h5ad`` has the same reference cells as ``sct_singleCell``.
    The wrapper accepts configurable filenames and subsets the original-counts
    file into an internal staging directory when preprocessing retained fewer
    normalized reference cells.
    """

    if ref_dir is None or spatial_dir is None:
        raise ValidationError("ref_dir and spatial_dir are required for double-modality covariation.")

    ref_path = Path(ref_dir)
    spatial_path = Path(spatial_dir)
    ref_sct = ref_path / config.ref_sct_name
    ref_original = ref_path / config.ref_original_counts_name
    spatial_sct = spatial_path / config.spatial_sct_name

    needs_ref_stage = config.ref_sct_name != "sct_singleCell.h5ad" or config.ref_original_counts_name != "Original_counts.h5ad"
    needs_spatial_stage = config.spatial_sct_name != "sct_spatial.h5ad"
    ref_obs = _obs_names(ref_sct)
    original_obs = _obs_names(ref_original)
    if ref_obs != original_obs:
        needs_ref_stage = True

    if not needs_ref_stage and not needs_spatial_stage:
        return ref_path, spatial_path, config.ref_original_counts_name

    stage_base = niche_result.output_dir / ".covariation_inputs"
    staged_ref = stage_base / "ref"
    staged_spatial = stage_base / "spatial"
    staged_ref.mkdir(parents=True, exist_ok=True)
    staged_spatial.mkdir(parents=True, exist_ok=True)

    if needs_ref_stage:
        _link_or_copy(ref_sct, staged_ref / "sct_singleCell.h5ad")
        _write_original_subset(ref_original, ref_obs, staged_ref / "Original_counts.h5ad", overwrite=config.overwrite)
    else:
        staged_ref = ref_path

    if needs_spatial_stage:
        _link_or_copy(spatial_sct, staged_spatial / "sct_spatial.h5ad")
    else:
        staged_spatial = spatial_path

    return staged_ref, staged_spatial, "Original_counts.h5ad"


def _as_nico_dir(path: str | Path) -> str:
    """Return a directory string with a trailing separator for NiCo."""

    text = os.fspath(Path(path))
    return text if text.endswith(os.sep) else text + os.sep


def _validate_nico_result(nico_result: Any) -> None:
    missing = [name for name in ("save_reg_coef", "covariation_dir", "pc_of_sp_clusterid") if not hasattr(nico_result, name)]
    if missing:
        raise RuntimeError("NiCo covariation result is missing expected fields: " + ", ".join(missing))


def _require_regression_state(result: CovariationResult) -> Any:
    if result.nico_result is not None:
        return result.nico_result
    if result.state_pickle is not None and result.state_pickle.exists():
        return read_covariation_state(result.state_pickle)
    raise ValueError(
        "Regression export requires covariation state containing NiCo save_reg_coef. "
        "Re-run with persist_state=True or load with load_state=True."
    )


def _resolve_ridge_cutoff(nico_result: Any, explicit: float | None) -> float:
    if explicit is not None:
        return float(explicit)
    if hasattr(nico_result, "coeff_cutoff_for_rid_reg"):
        return float(nico_result.coeff_cutoff_for_rid_reg)
    return 0.0


def _cell_type_names_from_state_or_disk(result: CovariationResult, nico_result: Any) -> dict[int, str]:
    if hasattr(nico_result, "spatialcell_unique_clusterid") and hasattr(nico_result, "spatialcell_unique_clustername"):
        ids = np.asarray(nico_result.spatialcell_unique_clusterid).tolist()
        names = np.asarray(nico_result.spatialcell_unique_clustername).tolist()
        return {int(cell_id): str(name) for cell_id, name in zip(ids, names)}
    used_ct = result.output_dir / "used_CT.txt"
    if used_ct.exists():
        return load_used_cell_types(used_ct)
    return {}


def _obs_names(path: Path) -> tuple[str, ...]:
    adata = read_h5ad(path, backed="r")
    try:
        return tuple(map(str, adata.obs_names))
    finally:
        adata.file.close()


def _write_original_subset(source: Path, obs_names: tuple[str, ...], destination: Path, *, overwrite: bool) -> None:
    if destination.exists() and not overwrite and _obs_names(destination) == obs_names:
        return
    original_obs = set(_obs_names(source))
    missing = [name for name in obs_names if name not in original_obs]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValidationError(f"Original reference counts are missing normalized reference cells: {preview}")
    adata = read_h5ad(source)
    subset = adata[list(obs_names), :].copy()
    destination.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(destination)


def _link_or_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    try:
        destination.symlink_to(source.resolve())
    except OSError:
        shutil.copy2(source, destination)


def _safe_index(values: np.ndarray, index: int) -> float:
    flat = np.ravel(values)
    if index >= flat.size:
        return float("nan")
    value = flat[index]
    return float(value) if np.isfinite(value) else float("nan")
