"""Configuration objects for NiCo niche interaction analysis."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class NeighborhoodConfig:
    """Configuration for spatial neighborhood construction.

    Parameters
    ----------
    radius
        Neighborhood mode. ``0`` selects NiCo's Delaunay-neighborhood mode;
        positive values select fixed-radius neighbors.
    epsilon_threshold
        Maximum Delaunay edge length used when ``radius == 0``.
    additional_excluded_cell_types
        Cell-type labels to exclude in addition to NiCo's implicit ``"NM"``
        exclusion.
    """

    radius: int | float = 0
    epsilon_threshold: float = 100.0
    additional_excluded_cell_types: tuple[str, ...] = ()


@dataclass(frozen=True)
class InteractionModelConfig:
    """Configuration for NiCo's logistic-regression interaction model.

    Parameters
    ----------
    k_fold
        Number of stratified cross-validation folds.
    n_repeats
        Number of repeated cross-validation rounds after hyperparameter
        selection.
    seed
        Random seed passed to NiCo.
    n_jobs
        Number of parallel sklearn jobs. ``-1`` uses all available CPUs.
    c_values
        Candidate inverse regularization strengths for sklearn logistic
        regression. If ``None``, NiCo's default grid is used.
    """

    k_fold: int = 5
    n_repeats: int = 1
    seed: int = 36851234
    n_jobs: int = -1
    c_values: tuple[float, ...] | None = None


@dataclass(frozen=True)
class NichePlotConfig:
    """Configuration for optional niche interaction plots.

    Parameters
    ----------
    enabled
        Whether plots should be generated after the main analysis.
    kinds
        Plot kinds to generate. Supported values are implemented in
        :mod:`nico_wrapper.niche.plotting`.
    saveas
        File extension used by upstream NiCo plotting functions, for example
        ``"pdf"`` or ``"png"``.
    dpi
        Plot resolution.
    transparent
        Whether plots should use transparent backgrounds.
    show
        Whether upstream plotting functions should keep figures open.
    interaction_cutoff
        Positive normalized coefficient cutoff used for graph-style plots.
    graph_edge_labels
        Whether graph plots should include edge-weight labels.
    choose_celltypes
        Cell types to plot for the ``"top-coefficients"`` kind. An empty tuple
        (default) plots all cell types, matching the original NiCo behaviour.
    """

    enabled: bool = False
    kinds: tuple[str, ...] = ("confusion", "coefficients", "scores", "graph")
    saveas: str = "pdf"
    dpi: int = 300
    transparent: bool = False
    show: bool = False
    interaction_cutoff: float = 0.1
    graph_edge_labels: bool = False
    choose_celltypes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ProximityConfig:
    """Configuration for observed-vs-randomized cell-type proximity analysis.

    Parameters
    ----------
    enabled
        Whether proximity analysis should run as part of the main pipeline.
    n_permutations
        Number of random label permutations used to build the null distribution.
    observed_threshold
        Minimum observed value for a pair to be included in the output.
    remove_self_pairs
        Whether to exclude same cell-type pairs from the analysis.
    as_counts
        If ``True``, co-localisation is measured as raw cell counts; if
        ``False``, proportions are used.
    seed
        Optional random seed for reproducible permutations.
    saveas
        Plot file extension, e.g. ``"pdf"`` or ``"png"``.
    show
        Whether to display the proximity plot inline (calls ``plt.show()``
        via NiCo internally).  Set ``True`` when running in a Jupyter
        notebook; keep ``False`` for non-interactive / CLI use.
    """

    enabled: bool = False
    n_permutations: int = 1000
    observed_threshold: float = 0.0
    remove_self_pairs: bool = True
    as_counts: bool = True
    seed: int | None = None
    saveas: str = "pdf"
    show: bool = False


@dataclass(frozen=True)
class NicheInteractionConfig:
    """Top-level configuration for complete niche interaction analysis.

    Parameters
    ----------
    anndata_filename
        Annotated spatial AnnData filename under ``output_dir``.
    label_key
        ``.obs`` column containing spatial cell-type labels.
    spatial_key
        ``.obsm`` key containing spatial coordinates.
    neighborhood
        Spatial neighborhood configuration.
    model
        Logistic-regression model configuration.
    plots
        Optional plotting configuration.
    proximity
        Optional proximity-analysis configuration.
    overwrite
        Whether existing output artifacts may be replaced.
    write_manifest
        Whether to write a JSON manifest sidecar.
    export_interactions
        Whether to export a TSV interaction coefficient table.
    """

    anndata_filename: str = "nico_celltype_annotation.h5ad"
    label_key: str = "nico_ct"
    spatial_key: str = "spatial"
    neighborhood: NeighborhoodConfig = field(default_factory=NeighborhoodConfig)
    model: InteractionModelConfig = field(default_factory=InteractionModelConfig)
    plots: NichePlotConfig = field(default_factory=NichePlotConfig)
    proximity: ProximityConfig = field(default_factory=ProximityConfig)
    overwrite: bool = False
    write_manifest: bool = True
    export_interactions: bool = True
