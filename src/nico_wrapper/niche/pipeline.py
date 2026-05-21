"""Public API for NiCo niche interaction analysis."""

from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
from typing import Any

import numpy as np

from .config import NicheInteractionConfig
from .io import (
    load_used_cell_types,
    metrics_from_score,
    normalize_radius,
    planned_artifact_paths,
    read_manifest,
    read_metrics_table,
    required_covariation_paths,
    write_interaction_table,
    write_manifest,
    write_metrics_table,
)
from .results import NicheArtifactPaths, NicheInteractionResult
from .validation import validate_covariation_artifacts, validate_niche_inputs


def run_niche_interactions(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> NicheInteractionResult:
    """Run complete NiCo niche interaction analysis after label transfer.

    Parameters
    ----------
    output_dir
        Directory containing the annotated spatial AnnData produced by
        ``nico-transfer``. The same directory receives NiCo's niche interaction
        artifacts.
    config
        Niche interaction configuration controlling input keys, neighborhood
        construction, model cross-validation, sidecar exports, and overwrite
        behavior.

    Returns
    -------
    NicheInteractionResult
        Stable wrapper result with paths to NiCo-compatible artifacts and
        wrapper sidecars.
    """

    validate_niche_inputs(output_dir, config=config)
    resolved_output_dir = Path(output_dir)

    nico_result = _run_upstream_spatial_neighborhood_analysis(
        output_dir=resolved_output_dir,
        config=config,
    )

    result = build_niche_result(
        resolved_output_dir,
        radius=config.neighborhood.radius,
        anndata_filename=config.anndata_filename,
        nico_result=nico_result,
    )

    metrics_path: Path | None = None
    if result.metrics:
        metrics_path = write_metrics_table(result.metrics, result.prediction_dir / f"metrics_{result.radius_tag}.tsv", radius=result.radius)
        result = replace(result, metrics_tsv=metrics_path)

    interactions_path: Path | None = None
    if config.export_interactions:
        interactions_path = export_interaction_table(result)
        result = replace(result, interactions_tsv=interactions_path)

    if config.write_manifest:
        manifest_path = write_manifest(result, config, result.prediction_dir / f"niche_manifest_{result.radius_tag}.json")
        result = replace(result, manifest_json=manifest_path)

    if config.plots.enabled:
        from .plotting import plot_niche_result

        plot_niche_result(result, config=config.plots)

    if config.proximity.enabled:
        from .proximity import run_proximity_analysis

        run_proximity_analysis(result, config=config.proximity)

    return result


def build_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
    nico_result: Any | None = None,
) -> NicheInteractionResult:
    """Build a :class:`NicheInteractionResult` from artifacts and metadata.

    Parameters
    ----------
    output_dir
        Base NiCo output directory containing niche artifacts.
    radius
        Radius value used for the niche analysis.
    anndata_filename
        Annotated AnnData filename under ``output_dir``.
    nico_result
        Optional raw NiCo result returned by
        ``spatial_neighborhood_analysis``. When present, selected ``C``, metrics,
        classes, and cell-type names are captured in the stable wrapper result.
    """

    paths = planned_artifact_paths(output_dir, radius=radius, anndata_filename=anndata_filename)
    _require_created_artifacts(paths)

    selected_c: float | None = None
    metrics: dict[str, tuple[float, float]] | None = None
    classes: tuple[int, ...] | None = None
    cell_type_names: dict[int, str] | None = None

    if nico_result is not None:
        if hasattr(nico_result, "lambda_c"):
            selected_c = float(nico_result.lambda_c)
        if hasattr(nico_result, "score"):
            metrics = metrics_from_score(nico_result.score)
        if hasattr(nico_result, "classes"):
            classes = tuple(int(value) for value in np.asarray(nico_result.classes).tolist())
        if hasattr(nico_result, "nameOfCellType"):
            cell_type_names = {int(key): str(value) for key, value in dict(nico_result.nameOfCellType).items()}

    if cell_type_names is None:
        cell_type_names = load_used_cell_types(paths.used_cell_types_tsv)
    if classes is None:
        classes = tuple(sorted(cell_type_names))

    covariation_ready = validate_covariation_artifacts(
        paths.output_dir,
        radius=paths.radius,
        anndata_filename=anndata_filename,
    )

    manifest = paths.manifest_json if paths.manifest_json.exists() else None
    metrics_tsv = paths.metrics_tsv if paths.metrics_tsv.exists() else None
    interactions_tsv = paths.interactions_tsv if paths.interactions_tsv.exists() else None

    return NicheInteractionResult(
        output_dir=paths.output_dir,
        annotated_h5ad=paths.annotated_h5ad,
        prediction_dir=paths.prediction_dir,
        radius=paths.radius,
        radius_tag=paths.radius_tag,
        used_cell_types_tsv=paths.used_cell_types_tsv,
        used_clusters_csv=paths.used_clusters_csv,
        neighbors_pickle=paths.neighbors_pickle,
        distances_pickle=paths.distances_pickle,
        normalized_neighborhood_npz=paths.normalized_neighborhood_npz,
        classifier_matrices_npz=paths.classifier_matrices_npz,
        manifest_json=manifest,
        metrics_tsv=metrics_tsv,
        interactions_tsv=interactions_tsv,
        selected_c=selected_c,
        metrics=metrics,
        cell_type_names=cell_type_names,
        classes=classes,
        covariation_ready=covariation_ready,
        nico_result=nico_result,
    )


def load_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
) -> NicheInteractionResult:
    """Load an existing niche interaction result from disk.

    The loader prefers the wrapper manifest when present and falls back to
    discovering standard NiCo artifacts. It does not rerun analysis and therefore
    cannot reconstruct in-memory-only NiCo fields such as ROC curves.
    """

    paths = planned_artifact_paths(output_dir, radius=radius, anndata_filename=anndata_filename)
    selected_c: float | None = None
    metrics: dict[str, tuple[float, float]] | None = None
    classes: tuple[int, ...] | None = None

    if paths.manifest_json.exists():
        manifest = read_manifest(paths.manifest_json)
        selected = manifest.get("model", {}).get("selected_c")
        selected_c = float(selected) if selected is not None else None
        raw_metrics = manifest.get("metrics")
        if isinstance(raw_metrics, dict):
            metrics = {str(key): (float(value[0]), float(value[1])) for key, value in raw_metrics.items()}
        raw_classes = manifest.get("classes")
        if raw_classes is not None:
            classes = tuple(int(value) for value in raw_classes)

    result = build_niche_result(
        output_dir,
        radius=paths.radius,
        anndata_filename=anndata_filename,
        nico_result=None,
    )
    if metrics is None and paths.metrics_tsv.exists():
        metrics = read_metrics_table(paths.metrics_tsv)
    return replace(result, selected_c=selected_c, metrics=metrics, classes=classes or result.classes)


def export_interaction_table(
    result: NicheInteractionResult,
    *,
    cutoff: float = 0.0,
    normalized_cutoff: float | None = None,
    include_self_edges: bool = True,
    output_path: str | Path | None = None,
) -> Path:
    """Export logistic-regression coefficients as a directed interaction table.

    Parameters
    ----------
    result
        Niche interaction result containing classifier matrix and cell-type
        artifacts.
    cutoff
        Positive normalized coefficient cutoff used when
        ``normalized_cutoff`` is not supplied.
    normalized_cutoff
        Explicit positive normalized coefficient cutoff.
    include_self_edges
        Whether to include source/target self-interactions in the table.
    output_path
        Optional output path. Defaults to ``interactions_{radius}.tsv`` under
        the niche prediction directory.
    """

    destination = Path(output_path) if output_path is not None else result.prediction_dir / f"interactions_{result.radius_tag}.tsv"
    return write_interaction_table(
        result,
        destination,
        cutoff=cutoff,
        normalized_cutoff=normalized_cutoff,
        include_self_edges=include_self_edges,
    )


def summarize_niche_result(result: NicheInteractionResult) -> dict[str, object]:
    """Return a compact dictionary summary of a niche interaction result."""

    return {
        "output_dir": result.output_dir,
        "annotated_h5ad": result.annotated_h5ad,
        "prediction_dir": result.prediction_dir,
        "radius": result.radius,
        "n_cell_types": len(result.cell_type_names or {}),
        "selected_c": result.selected_c,
        "covariation_ready": result.covariation_ready,
        "classifier_matrices_npz": result.classifier_matrices_npz,
        "interactions_tsv": result.interactions_tsv,
        "metrics_tsv": result.metrics_tsv,
        "manifest_json": result.manifest_json,
    }


def _run_upstream_spatial_neighborhood_analysis(
    *,
    output_dir: Path,
    config: NicheInteractionConfig,
) -> Any:
    try:
        from nico import Interactions as sint
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo niche interaction analysis requires the 'nico' package.") from exc

    model = config.model
    neighborhood = config.neighborhood
    radius_value, _ = normalize_radius(neighborhood.radius)
    c_values = list(model.c_values) if model.c_values is not None else list(np.power(2.0, np.arange(-12, 12)))

    return sint.spatial_neighborhood_analysis(
        output_nico_dir=_as_nico_dir(output_dir),
        anndata_object_name=config.anndata_filename,
        spatial_cluster_tag=config.label_key,
        spatial_coordinate_tag=config.spatial_key,
        Radius=radius_value,
        n_repeats=model.n_repeats,
        K_fold=model.k_fold,
        seed=model.seed,
        n_jobs=model.n_jobs,
        lambda_c_ranges=c_values,
        epsilonThreshold=neighborhood.epsilon_threshold,
        removed_CTs_before_finding_CT_CT_interactions=list(neighborhood.additional_excluded_cell_types),
    )


def _as_nico_dir(path: str | Path) -> str:
    """Return a directory string with a trailing separator for NiCo."""

    text = os.fspath(Path(path))
    return text if text.endswith(os.sep) else text + os.sep


def _require_created_artifacts(paths: NicheArtifactPaths) -> None:
    missing: list[Path] = []
    empty: list[Path] = []
    for path in (*required_covariation_paths(paths), paths.normalized_neighborhood_npz):
        if not path.exists() or not path.is_file():
            missing.append(path)
        elif path.stat().st_size == 0:
            empty.append(path)
    if missing:
        raise RuntimeError("NiCo niche analysis did not create expected artifacts: " + ", ".join(map(str, missing)))
    if empty:
        raise RuntimeError("NiCo niche analysis created empty artifacts: " + ", ".join(map(str, empty)))
