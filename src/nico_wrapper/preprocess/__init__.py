"""Utilities for preparing NiCo reference and spatial query inputs."""

from .config import NiCoSCTransformConfig, PearsonResidualsConfig
from .io import load_reference_from_sparse, load_spatial_query_from_csv, load_spatial_query_from_h5ad
from .pipeline import preprocess_nico_inputs
from .validation import ValidationError

__all__ = [
    "ValidationError",
    "NiCoSCTransformConfig",
    "PearsonResidualsConfig",
    "load_reference_from_sparse",
    "load_spatial_query_from_csv",
    "load_spatial_query_from_h5ad",
    "preprocess_nico_inputs",
]
