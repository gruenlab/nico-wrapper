"""Utilities for NiCo niche / spatial cell-type interaction analysis."""

from .config import (
    InteractionModelConfig,
    NeighborhoodConfig,
    NicheInteractionConfig,
    NichePlotConfig,
    ProximityConfig,
)
from .io import normalize_radius, planned_artifact_paths
from .pipeline import (
    build_niche_result,
    export_interaction_table,
    load_niche_result,
    run_niche_interactions,
    summarize_niche_result,
)
from .results import NicheArtifactPaths, NicheInteractionResult, ProximityResult
from .validation import (
    ValidationError,
    validate_annotated_spatial_adata,
    validate_covariation_artifacts,
    validate_niche_config,
    validate_niche_inputs,
)

__all__ = [
    "InteractionModelConfig",
    "NeighborhoodConfig",
    "NicheInteractionConfig",
    "NichePlotConfig",
    "ProximityConfig",
    "NicheArtifactPaths",
    "NicheInteractionResult",
    "ProximityResult",
    "ValidationError",
    "normalize_radius",
    "planned_artifact_paths",
    "validate_annotated_spatial_adata",
    "validate_covariation_artifacts",
    "validate_niche_config",
    "validate_niche_inputs",
    "build_niche_result",
    "export_interaction_table",
    "load_niche_result",
    "run_niche_interactions",
    "summarize_niche_result",
]
