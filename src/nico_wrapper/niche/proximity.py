"""Cell-type proximity analysis wrappers for NiCo niche results."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import ProximityConfig
from .plotting import _to_nico_plot_namespace
from .results import NicheInteractionResult, ProximityResult


def run_proximity_analysis(
    result: NicheInteractionResult,
    *,
    config: ProximityConfig = ProximityConfig(),
) -> ProximityResult:
    """Run NiCo's observed-vs-randomized cell-type proximity analysis.

    Parameters
    ----------
    result
        Existing niche interaction result with neighbor and cell-type artifacts.
    config
        Proximity analysis configuration.

    Returns
    -------
    ProximityResult
        Observed pair values, observed/randomized ratios, and output paths.
    """

    try:
        from nico import Interactions as sint
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo proximity analysis requires the 'nico' package.") from exc

    namespace = _to_nico_plot_namespace(result)
    rng_state: Any | None = None
    if config.seed is not None:
        rng_state = np.random.get_state()
        np.random.seed(config.seed)
    try:
        observed, ratio = sint.visualization_of_top_celltype_proximity_pairs(
            namespace,
            saveas=config.saveas,
            showit=False,
            n_rand_permute=config.n_permutations,
            Observed_Threshold=config.observed_threshold,
            remove_self_pairs=config.remove_self_pairs,
            celltype_proximity_as_counts=config.as_counts,
        )
    finally:
        if rng_state is not None:
            np.random.set_state(rng_state)

    observed_tsv = _write_pair_table(
        observed,
        result.prediction_dir / f"proximity_observed_{result.radius_tag}.tsv",
        value_name="observed",
    )
    ratio_tsv = _write_pair_table(
        ratio,
        result.prediction_dir / f"proximity_ratio_{result.radius_tag}.tsv",
        value_name="observed_over_randomized",
    )
    plot_path = result.prediction_dir / f"celltype_proximity_pairs{result.radius_tag}.{config.saveas}"
    return ProximityResult(
        observed_tsv=observed_tsv,
        ratio_tsv=ratio_tsv,
        plot_path=plot_path if plot_path.exists() else None,
        observed={str(key): float(value) for key, value in observed.items()},
        ratio={str(key): float(value) for key, value in ratio.items()},
    )


def _write_pair_table(values: dict[str, float], path: Path, *, value_name: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for pair, value in values.items():
        central, neighbor = _split_pair(str(pair))
        rows.append(
            {
                "pair": pair,
                "central_cell_type": central,
                "neighbor_cell_type": neighbor,
                value_name: float(value),
            }
        )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _split_pair(pair: str) -> tuple[str, str]:
    if ":" not in pair:
        return pair, ""
    central, neighbor = pair.split(":", 1)
    return central, neighbor
