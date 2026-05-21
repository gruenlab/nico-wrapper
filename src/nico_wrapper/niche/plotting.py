"""Plotting wrappers for NiCo niche interaction results."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from .config import NichePlotConfig
from .results import NicheInteractionResult

SUPPORTED_PLOT_KINDS: tuple[str, ...] = (
    "confusion",
    "coefficients",
    "scores",
    "graph",
    "roc",
    "predicted-probabilities",
    "top-coefficients",
    "all",
)


def plot_niche_result(result: NicheInteractionResult, *, config: NichePlotConfig = NichePlotConfig()) -> list[Path]:
    """Generate selected plots for a niche interaction result.

    Parameters
    ----------
    result
        Result returned by :func:`nico_wrapper.niche.run_niche_interactions` or
        :func:`nico_wrapper.niche.load_niche_result`.
    config
        Plotting configuration.

    Returns
    -------
    list[pathlib.Path]
        Expected output paths for plots requested by ``config``. Some upstream
        NiCo plotting functions compute filenames internally; the returned paths
        are best-effort expected locations.
    """

    try:
        from nico import Interactions as sint
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo plotting requires the 'nico' package.") from exc

    namespace = _to_nico_plot_namespace(result)
    kinds = _expand_plot_kinds(config.kinds)
    paths: list[Path] = []

    if "confusion" in kinds:
        sint.plot_confusion_matrix(
            namespace,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"Confusing_matrix_R{result.radius_tag}.{config.saveas}")

    if "coefficients" in kinds:
        sint.plot_coefficient_matrix(
            namespace,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"weight_matrix_R{result.radius_tag}.{config.saveas}")

    if "scores" in kinds:
        _require_fresh_or_metric_data(result, kind="scores")
        sint.plot_evaluation_scores(
            namespace,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"scores_{result.radius_tag}.{config.saveas}")

    if "graph" in kinds:
        graph_fn = (
            sint.plot_niche_interactions_with_edge_weight
            if config.graph_edge_labels
            else sint.plot_niche_interactions_without_edge_weight
        )
        graph_fn(
            namespace,
            niche_cutoff=config.interaction_cutoff,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        name = "with" if config.graph_edge_labels else "without"
        paths.append(result.prediction_dir / f"Niche_interactions_{name}_edge_weights_R{result.radius_tag}.{config.saveas}")

    if "top-coefficients" in kinds:
        sint.find_interacting_cell_types(
            namespace,
            celltype_niche_interaction_cutoff=config.interaction_cutoff,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"TopCoeff_R{result.radius_tag}")

    if "roc" in kinds:
        if result.nico_result is None:
            raise ValueError("ROC plots require a fresh result with in-memory ROC arrays from NiCo.")
        sint.plot_roc_results(
            namespace,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"ROC_R{result.radius_tag}.{config.saveas}")

    if "predicted-probabilities" in kinds:
        if result.nico_result is None:
            raise ValueError("Predicted-probability plots require a fresh result with in-memory arrays from NiCo.")
        sint.plot_predicted_probabilities(
            namespace,
            saveas=config.saveas,
            showit=config.show,
            transparent_mode=config.transparent,
            dpi=config.dpi,
        )
        paths.append(result.prediction_dir / f"predicted_probability_R{result.radius_tag}.{config.saveas}")

    return paths


def _to_nico_plot_namespace(result: NicheInteractionResult) -> SimpleNamespace:
    if result.nico_result is not None:
        return result.nico_result
    if result.metrics is None:
        score = None
    else:
        from .io import METRIC_NAMES
        import numpy as np

        score = np.array([result.metrics.get(name, (float("nan"), float("nan"))) for name in METRIC_NAMES])
    return SimpleNamespace(
        outputdir=_as_trailing_dir(result.output_dir),
        fout=str(result.classifier_matrices_npz),
        niche_pred_outdir=_as_trailing_dir(result.prediction_dir),
        nameOfCellType=result.cell_type_names or {},
        Radius=result.radius,
        BothLinearAndCrossTerms=1,
        classes=list(result.classes or sorted((result.cell_type_names or {}).keys())),
        lambda_c=result.selected_c,
        score=score,
    )


def _expand_plot_kinds(kinds: tuple[str, ...]) -> set[str]:
    requested = {kind.strip().lower() for kind in kinds}
    unknown = requested.difference(SUPPORTED_PLOT_KINDS)
    if unknown:
        raise ValueError(f"Unsupported plot kind(s): {sorted(unknown)}")
    if "all" in requested:
        return {"confusion", "coefficients", "scores", "graph", "top-coefficients"}
    return requested


def _require_fresh_or_metric_data(result: NicheInteractionResult, *, kind: str) -> None:
    if result.nico_result is None and result.metrics is None:
        raise ValueError(f"{kind} plots require metric data from a fresh run or manifest.")


def _as_trailing_dir(path: Path) -> str:
    text = str(path)
    return text if text.endswith("/") else text + "/"
