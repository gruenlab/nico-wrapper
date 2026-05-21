"""Utilities for NiCo label transfer from reference to spatial/query data."""

from .config import (
    AnchorConfig,
    AnchorResult,
    AnnotationConfig,
    AnnotationResult,
    LabelTransferConfig,
    TieStrategy,
    TransferOutputs,
)
from .pipeline import find_anchors, run_label_transfer, save_transfer_result, transfer_labels
from .validation import ValidationError, validate_label_transfer_config, validate_label_transfer_inputs

__all__ = [
    "AnchorConfig",
    "AnchorResult",
    "AnnotationConfig",
    "AnnotationResult",
    "LabelTransferConfig",
    "TieStrategy",
    "TransferOutputs",
    "ValidationError",
    "validate_label_transfer_config",
    "validate_label_transfer_inputs",
    "find_anchors",
    "run_label_transfer",
    "save_transfer_result",
    "transfer_labels",
]
