"""Configuration objects for NiCo preprocessing."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NiCoSCTransformConfig:
    """Configuration for NiCo's SCTransform-like normalization.

    This corresponds to ``nico.Annotations.SCTransform`` / ``Annotations.SCTransform``
    as used in the original data preparation notebooks.

    Parameters
    ----------
    min_cells
        Minimum number of cells in which a gene must have counts before being
        retained by SCTransform.
    gmean_eps
        Small offset used when computing geometric means.
    n_genes
        Number of genes sampled for model fitting. ``None`` uses all eligible
        genes.
    n_cells
        Number of cells sampled for model fitting. ``None`` uses all cells.
    bin_size
        Number of genes per fitting bin.
    bw_adjust
        Bandwidth adjustment used by the SCTransform fitting procedure.
    """

    min_cells: int = 1
    gmean_eps: float = 1.0
    n_genes: int | None = 500
    n_cells: int | None = None
    bin_size: int = 500
    bw_adjust: float = 3.0


@dataclass(frozen=True)
class PearsonResidualsConfig:
    """Configuration for Scanpy analytic Pearson residual normalization.

    This corresponds to ``scanpy.experimental.pp.normalize_pearson_residuals``.

    Parameters
    ----------
    theta
        Negative-binomial overdispersion parameter. Larger values imply less
        overdispersion; ``float("inf")`` corresponds approximately to a Poisson
        model.
    clip
        Residual clipping threshold. ``None`` uses Scanpy's default
        ``sqrt(n_obs)``. ``float("inf")`` disables clipping.
    check_values
        Whether Scanpy should check that the input contains raw integer counts.
    layer
        Optional layer to normalize instead of ``.X``.
    """

    theta: float = 100.0
    clip: float | None = None
    check_values: bool = True
    layer: str | None = None


NormalizationConfig = NiCoSCTransformConfig | PearsonResidualsConfig
