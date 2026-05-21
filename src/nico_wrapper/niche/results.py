"""Result dataclasses for NiCo niche interaction analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class NicheArtifactPaths:
    """Expected file locations for one niche interaction run."""

    output_dir: Path
    annotated_h5ad: Path
    prediction_dir: Path
    radius: int | float
    radius_tag: str
    used_cell_types_tsv: Path
    used_clusters_csv: Path
    neighbors_pickle: Path
    distances_pickle: Path
    normalized_neighborhood_npz: Path
    classifier_matrices_npz: Path
    manifest_json: Path
    metrics_tsv: Path
    interactions_tsv: Path


@dataclass(frozen=True)
class NicheInteractionResult:
    """Stable result object returned by the niche interaction wrapper.

    The object exposes the NiCo-compatible artifact paths required by downstream
    covariation analysis plus wrapper sidecars such as metrics and interaction
    tables. ``nico_result`` is intentionally optional because it is only
    available immediately after a fresh upstream NiCo run.
    """

    output_dir: Path
    annotated_h5ad: Path
    prediction_dir: Path
    radius: int | float
    radius_tag: str
    used_cell_types_tsv: Path
    used_clusters_csv: Path
    neighbors_pickle: Path
    distances_pickle: Path
    normalized_neighborhood_npz: Path
    classifier_matrices_npz: Path
    manifest_json: Path | None = None
    metrics_tsv: Path | None = None
    interactions_tsv: Path | None = None
    selected_c: float | None = None
    metrics: dict[str, tuple[float, float]] | None = None
    cell_type_names: dict[int, str] | None = None
    classes: tuple[int, ...] | None = None
    covariation_ready: bool = False
    nico_result: Any | None = None


@dataclass(frozen=True)
class ProximityResult:
    """Result from observed-vs-randomized cell-type proximity analysis."""

    observed_tsv: Path | None
    ratio_tsv: Path | None
    plot_path: Path | None
    observed: dict[str, float]
    ratio: dict[str, float]
