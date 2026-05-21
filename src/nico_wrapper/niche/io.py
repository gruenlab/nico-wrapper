"""I/O helpers for NiCo niche interaction artifacts."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .config import NicheInteractionConfig
from .results import NicheArtifactPaths, NicheInteractionResult

METRIC_NAMES: tuple[str, ...] = (
    "accuracy",
    "macro_f1",
    "macro_precision",
    "macro_recall",
    "micro_f1",
    "micro_precision",
    "micro_recall",
    "weighted_f1",
    "weighted_precision",
    "weighted_recall",
    "cohen_kappa",
    "log_loss",
    "matthews_corrcoef",
    "hamming_loss",
    "zero_one_loss",
)

REQUIRED_CLASSIFIER_KEYS: tuple[str, ...] = ("cmn", "coef", "cmn_std", "coef_std", "CTFeatures")


def normalize_radius(value: int | float | str) -> tuple[int | float, str]:
    """Normalize a user radius value and return ``(value, file_tag)``.

    Integer-like values are normalized to integers so that ``0``, ``0.0``, and
    ``"0"`` all produce NiCo-compatible artifact names such as
    ``neighbors_0.p``.
    """

    if isinstance(value, bool):
        raise ValueError("radius must be numeric, not boolean.")
    if isinstance(value, int):
        if value < 0:
            raise ValueError("radius must be non-negative.")
        return value, str(value)
    if isinstance(value, float):
        if not np.isfinite(value) or value < 0:
            raise ValueError("radius must be a finite non-negative number.")
        if value.is_integer():
            integer = int(value)
            return integer, str(integer)
        return value, _format_float_tag(value)
    text = str(value).strip()
    if not text:
        raise ValueError("radius must be non-empty.")
    try:
        parsed = float(text)
    except ValueError as exc:
        raise ValueError(f"radius must be numeric; received {value!r}.") from exc
    if not np.isfinite(parsed) or parsed < 0:
        raise ValueError("radius must be a finite non-negative number.")
    if parsed.is_integer():
        integer = int(parsed)
        return integer, str(integer)
    return parsed, _format_float_tag(parsed)


def planned_artifact_paths(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
) -> NicheArtifactPaths:
    """Return all standard artifact paths for a niche interaction run."""

    resolved_output_dir = Path(output_dir)
    radius_value, radius_tag = normalize_radius(radius)
    prediction_dir = resolved_output_dir / "niche_prediction_linear"
    return NicheArtifactPaths(
        output_dir=resolved_output_dir,
        annotated_h5ad=resolved_output_dir / anndata_filename,
        prediction_dir=prediction_dir,
        radius=radius_value,
        radius_tag=radius_tag,
        used_cell_types_tsv=resolved_output_dir / "used_CT.txt",
        used_clusters_csv=resolved_output_dir / f"used_Clusters{radius_tag}.csv",
        neighbors_pickle=resolved_output_dir / f"neighbors_{radius_tag}.p",
        distances_pickle=resolved_output_dir / f"distances_{radius_tag}.p",
        normalized_neighborhood_npz=prediction_dir / f"normalized_spatial_neighborhood_{radius_tag}.npz",
        classifier_matrices_npz=prediction_dir / f"classifier_matrices_{radius_tag}.npz",
        manifest_json=prediction_dir / f"niche_manifest_{radius_tag}.json",
        metrics_tsv=prediction_dir / f"metrics_{radius_tag}.tsv",
        interactions_tsv=prediction_dir / f"interactions_{radius_tag}.tsv",
    )


def required_covariation_paths(paths: NicheArtifactPaths) -> tuple[Path, ...]:
    """Return files required by downstream NiCo covariation."""

    return (
        paths.used_cell_types_tsv,
        paths.used_clusters_csv,
        paths.neighbors_pickle,
        paths.distances_pickle,
        paths.classifier_matrices_npz,
    )


def core_output_paths(paths: NicheArtifactPaths) -> tuple[Path, ...]:
    """Return files written by upstream NiCo's main niche analysis."""

    return (*required_covariation_paths(paths), paths.normalized_neighborhood_npz)


def wrapper_output_paths(paths: NicheArtifactPaths) -> tuple[Path, ...]:
    """Return wrapper sidecar files that may be written after analysis."""

    return (paths.manifest_json, paths.metrics_tsv, paths.interactions_tsv)


def load_used_cell_types(path: str | Path) -> dict[int, str]:
    """Load NiCo's ``used_CT.txt`` as ``{cell_type_id: cell_type_name}``."""

    table = pd.read_csv(path, sep="\t", header=None, comment="#")
    if table.shape[1] < 2:
        raise ValueError(f"used_CT file must contain at least two tab-separated columns: {path}")
    result: dict[int, str] = {}
    for _, row in table.iterrows():
        result[int(row.iloc[0])] = str(row.iloc[1])
    return result


def load_classifier_matrices(path: str | Path) -> Mapping[str, np.ndarray]:
    """Load and validate NiCo's classifier matrix ``.npz`` artifact."""

    data = np.load(path, allow_pickle=True)
    missing = [key for key in REQUIRED_CLASSIFIER_KEYS if key not in data]
    if missing:
        raise ValueError(f"Classifier matrix file is missing keys {missing}: {path}")
    return {key: data[key] for key in data.files}


def metrics_from_score(score: Any) -> dict[str, tuple[float, float]]:
    """Convert NiCo's ``score`` array to a metric mapping."""

    array = np.asarray(score, dtype=float)
    if array.ndim != 2 or array.shape[1] < 2:
        raise ValueError(f"Expected score array with shape (n_metrics, 2); found {array.shape}.")
    metrics: dict[str, tuple[float, float]] = {}
    for index, name in enumerate(METRIC_NAMES):
        if index >= array.shape[0]:
            break
        metrics[name] = (float(array[index, 0]), float(array[index, 1]))
    return metrics


def write_metrics_table(
    metrics: Mapping[str, tuple[float, float]],
    output_path: str | Path,
    *,
    radius: int | float,
) -> Path:
    """Write metric means and standard deviations to a TSV file."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"metric": name, "mean": values[0], "std": values[1], "radius": radius}
        for name, values in metrics.items()
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def read_metrics_table(path: str | Path) -> dict[str, tuple[float, float]]:
    """Read a wrapper ``metrics_{radius}.tsv`` file."""

    table = pd.read_csv(path, sep="\t")
    required = {"metric", "mean", "std"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Metrics table is missing columns {sorted(missing)}: {path}")
    return {str(row.metric): (float(row.mean), float(row.std)) for row in table.itertuples(index=False)}


def write_interaction_table(
    result: NicheInteractionResult,
    output_path: str | Path,
    *,
    cutoff: float = 0.0,
    normalized_cutoff: float | None = None,
    include_self_edges: bool = True,
) -> Path:
    """Write the logistic-regression coefficient matrix as a directed edge table.

    In NiCo's public linear path, coefficient rows are central/predicted cell
    types and columns are neighboring/source cell types. The exported direction
    is therefore ``source_cell_type -> target_cell_type``.
    """

    matrices = load_classifier_matrices(result.classifier_matrices_npz)
    coef = np.asarray(matrices["coef"], dtype=float)
    coef_std = np.asarray(matrices.get("coef_std", np.zeros_like(coef)), dtype=float)
    features = [str(value) for value in np.asarray(matrices["CTFeatures"]).tolist()]

    if coef.ndim != 2:
        raise ValueError(f"Expected 2D coefficient matrix; found shape {coef.shape}.")
    if any(" " in feature for feature in features):
        raise ValueError("Cross-term CTFeatures are not supported by the interaction table exporter.")

    cell_type_names = result.cell_type_names or load_used_cell_types(result.used_cell_types_tsv)
    target_ids = result.classes if result.classes is not None and len(result.classes) == coef.shape[0] else tuple(range(coef.shape[0]))
    source_ids = [_feature_to_cell_type_id(feature, fallback=index) for index, feature in enumerate(features[: coef.shape[1]])]

    max_abs = float(np.max(np.abs(coef))) if coef.size else 0.0
    norm_threshold = cutoff if normalized_cutoff is None else normalized_cutoff
    rows: list[dict[str, object]] = []
    for row_index, target_id in enumerate(target_ids):
        for col_index, source_id in enumerate(source_ids):
            if row_index >= coef.shape[0] or col_index >= coef.shape[1]:
                continue
            is_self = int(source_id) == int(target_id)
            if is_self and not include_self_edges:
                continue
            value = float(coef[row_index, col_index])
            std_value = float(coef_std[row_index, col_index]) if coef_std.shape == coef.shape else float("nan")
            normalized = value / max_abs if max_abs else float("nan")
            rows.append(
                {
                    "source_cell_type": cell_type_names.get(int(source_id), str(source_id)),
                    "target_cell_type": cell_type_names.get(int(target_id), str(target_id)),
                    "source_cell_type_id": int(source_id),
                    "target_cell_type_id": int(target_id),
                    "coefficient": value,
                    "coefficient_std": std_value,
                    "normalized_coefficient": normalized,
                    "abs_normalized_coefficient": abs(normalized) if np.isfinite(normalized) else float("nan"),
                    "is_self_interaction": is_self,
                    "passes_positive_cutoff": bool(np.isfinite(normalized) and normalized > norm_threshold),
                    "radius": result.radius,
                }
            )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def write_manifest(result: NicheInteractionResult, config: NicheInteractionConfig, output_path: str | Path) -> Path:
    """Write a JSON manifest describing a niche interaction run."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "wrapper_version": _package_version("nico-wrapper"),
        "nico_package": _package_version("nico-sc-sp"),
        "output_dir": str(result.output_dir),
        "annotated_h5ad": str(result.annotated_h5ad),
        "prediction_dir": str(result.prediction_dir),
        "label_key": config.label_key,
        "spatial_key": config.spatial_key,
        "radius": result.radius,
        "radius_tag": result.radius_tag,
        "epsilon_threshold": config.neighborhood.epsilon_threshold,
        "excluded_cell_types": ["NM", *config.neighborhood.additional_excluded_cell_types],
        "model": {
            "k_fold": config.model.k_fold,
            "n_repeats": config.model.n_repeats,
            "seed": config.model.seed,
            "n_jobs": config.model.n_jobs,
            "c_values": list(config.model.c_values) if config.model.c_values is not None else None,
            "selected_c": result.selected_c,
        },
        "artifacts": {
            "used_cell_types_tsv": str(result.used_cell_types_tsv),
            "used_clusters_csv": str(result.used_clusters_csv),
            "neighbors_pickle": str(result.neighbors_pickle),
            "distances_pickle": str(result.distances_pickle),
            "normalized_neighborhood_npz": str(result.normalized_neighborhood_npz),
            "classifier_matrices_npz": str(result.classifier_matrices_npz),
            "metrics_tsv": str(result.metrics_tsv) if result.metrics_tsv else None,
            "interactions_tsv": str(result.interactions_tsv) if result.interactions_tsv else None,
        },
        "covariation_ready": result.covariation_ready,
        "metrics": result.metrics,
        "cell_type_names": result.cell_type_names,
        "classes": result.classes,
        "config": asdict(config),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def read_manifest(path: str | Path) -> dict[str, Any]:
    """Read a niche manifest JSON file."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def _feature_to_cell_type_id(feature: str, *, fallback: int) -> int:
    if feature.startswith("x"):
        suffix = feature[1:]
        if suffix.isdigit():
            return int(suffix)
    return fallback


def _format_float_tag(value: float) -> str:
    return str(value)


def _package_version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None
