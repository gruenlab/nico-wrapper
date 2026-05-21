"""Configuration objects for NiCo label transfer."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any


class TieStrategy(str, Enum):
    """Strategies for resolving ambiguous spatial-cell label assignments."""

    MAJORITY = "majority"
    WEIGHTED = "weighted"


@dataclass(frozen=True)
class AnchorConfig:
    """Configuration for finding reference/query anchors.

    Parameters
    ----------
    neighbors
        Number of nearest neighbors used both for mutual-nearest-neighbor anchor
        discovery and for the spatial KNN graph used during label propagation.
    n_pcs
        Number of principal components used for the shared transfer space.
        NiCo computes PCs from the normalized reference data and projects both
        reference and spatial data into that space.
    minkowski_order
        Minkowski distance order used by ``scipy.spatial.cKDTree``. ``2`` is
        Euclidean distance and matches NiCo's default.
    spatial_sct_filename
        File name of the normalized spatial/query AnnData under ``spatial_dir``.
    sc_sct_filename
        File name of the normalized single-cell/reference AnnData under
        ``ref_dir``.
    sc_full_filename
        File name of the full/original reference AnnData under ``ref_dir``.
        This file must contain the reference labels in ``.obs``.
    """

    neighbors: int = 50
    n_pcs: int = 50
    minkowski_order: int = 2
    spatial_sct_filename: str = "sct_spatial.h5ad"
    sc_sct_filename: str = "sct_singleCell.h5ad"
    sc_full_filename: str = "Original_counts.h5ad"


@dataclass(frozen=True)
class AnnotationConfig:
    """Configuration for NiCo-based label transfer and propagation.

    Parameters
    ----------
    ref_label_key
        Column in ``Original_counts.h5ad.obs`` containing reference cell labels.
    spatial_cluster_key
        Column in ``sct_spatial.h5ad.obs`` containing spatial guide clusters.
        NiCo uses this to prune noisy anchors and constrain propagation. The
        default follows NiCo's own default.
    dispersion_cutoff
        Cutoff used by NiCo when pruning anchors distributed across spatial
        clusters.
    iterations
        Number of iterative rounds used to propagate labels from anchored cells
        to unmapped spatial cells.
    tie_strategy
        How ambiguous assignments should be resolved. ``majority`` maps to
        NiCo's ``resolved_tie_issue_with_weighted_nearest_neighbor='No'``;
        ``weighted`` maps to ``'Yes'``.
    output_label_key
        Column name to write into the annotated spatial AnnData ``.obs``.
    output_h5ad_name
        File name for the annotated spatial AnnData written under ``output_dir``.
    """

    ref_label_key: str = "cluster"
    spatial_cluster_key: str = "leiden0.5"
    dispersion_cutoff: float = 0.15
    iterations: int = 3
    tie_strategy: TieStrategy = TieStrategy.MAJORITY
    output_label_key: str = "nico_ct"
    output_h5ad_name: str = "nico_celltype_annotation.h5ad"


@dataclass(frozen=True)
class LabelTransferConfig:
    """Top-level configuration for the complete label-transfer workflow."""

    anchors: AnchorConfig = field(default_factory=AnchorConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    overwrite: bool = False
    cleanup_intermediate: bool = False


@dataclass(frozen=True)
class AnchorResult:
    """Result metadata from anchor discovery.

    ``nico_result`` is intentionally opaque for now: NiCo returns a
    ``types.SimpleNamespace`` with in-memory AnnData objects and path metadata.
    The wrapper gives callers stable path fields while preserving the original
    NiCo object for the later annotation step.
    """

    nico_result: Any
    output_dir: Path
    annotation_dir: Path
    anchors_npz: Path


@dataclass(frozen=True)
class AnnotationResult:
    """Result metadata from label propagation before saving the final AnnData."""

    nico_result: Any
    output_dir: Path
    annotation_dir: Path
    iteration_cluster_csvs: tuple[Path, ...]
    iteration_celltype_csvs: tuple[Path, ...]


@dataclass(frozen=True)
class TransferOutputs:
    """Paths produced by the complete label-transfer workflow."""

    output_dir: Path
    annotation_dir: Path
    annotated_h5ad: Path
    anchors_npz: Path | None
    iteration_cluster_csvs: tuple[Path, ...]
    iteration_celltype_csvs: tuple[Path, ...]
