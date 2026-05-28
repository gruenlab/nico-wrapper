"""Marker-based annotation quality-control metric APIs.

This module intentionally contains public method stubs only. Implementations will
be added after the API is finalized.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

from anndata import AnnData
import pandas as pd

MarkerSets = Mapping[str, Sequence[str]]



def marker_set_summary(
    adata: AnnData,
    markers: MarkerSets,
    *,
    label_key: str = "cell_type",
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Summarize available marker genes and cell counts per annotated cell type.

    Returns one row per cell type with the number of cells, number of provided
    markers, number of markers found in the AnnData object, and optionally
    missing marker names.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    label_key
        Column in ``adata.obs`` containing cell type annotations.
    layer
        Optional layer name to use when checking available genes/data source.
    use_raw
        Whether to use ``adata.raw`` instead of ``adata.X``.
    gene_symbols_key
        Optional column in ``adata.var`` or ``adata.raw.var`` containing gene
        symbols if markers are not keyed by ``var_names``.
    labels
        Optional subset/order of cell types to evaluate.

    Returns
    -------
    pandas.DataFrame
        Per-cell-type marker availability summary.
    """
    raise NotImplementedError


def marker_logfc_metrics(
    adata: AnnData,
    markers: MarkerSets,
    *,
    label_key: str = "cell_type",
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
    pseudocount: float = 1e-9,
    logfc_threshold: float = 0.25,
    min_cells: int = 1,
    return_per_marker: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Compute one-vs-rest marker expression enrichment metrics per cell type.

    For each marker gene of each cell type, computes log2 fold change between
    cells annotated as that type and all remaining cells. Summarizes marker
    enrichment as mean marker logFC and fraction of markers above a configurable
    logFC threshold.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    label_key
        Column in ``adata.obs`` containing cell type annotations.
    layer
        Optional layer to use for expression values.
    use_raw
        Whether to use ``adata.raw``.
    gene_symbols_key
        Optional gene-symbol column in ``adata.var``.
    labels
        Optional subset/order of cell types to evaluate.
    pseudocount
        Small value added to group means before computing log fold change.
    logfc_threshold
        Threshold used for the enriched-marker fraction.
    min_cells
        Minimum number of cells required for a cell type to be evaluated.
    return_per_marker
        If True, also return a per-marker table with mean expression and logFC.

    Returns
    -------
    pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
        Summary metrics, and optionally per-marker logFC values.
    """
    raise NotImplementedError


def marker_set_scores(
    adata: AnnData,
    markers: MarkerSets,
    *,
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
    ctrl_size: int = 50,
    gene_pool: Sequence[str] | None = None,
    n_bins: int = 25,
    random_state: int | None = 0,
    score_prefix: str = "marker_score__",
    copy_scores_to_obs: bool = False,
) -> pd.DataFrame:
    """Compute marker-set/module scores for each marker set.

    Intended to use Scanpy's ``score_genes`` and return a cell-by-cell-type
    score matrix where rows are cells and columns are marker-set labels.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    layer
        Optional layer to use for expression values.
    use_raw
        Whether to use ``adata.raw``.
    gene_symbols_key
        Optional gene-symbol column in ``adata.var``.
    labels
        Optional subset/order of marker sets to score.
    ctrl_size
        Number of control genes sampled per expression bin by Scanpy.
    gene_pool
        Optional background gene pool for control gene selection.
    n_bins
        Number of expression bins used by Scanpy.
    random_state
        Random seed for reproducible control-gene sampling.
    score_prefix
        Prefix used if scores are written to ``adata.obs``.
    copy_scores_to_obs
        If True, store computed scores in ``adata.obs``.

    Returns
    -------
    pandas.DataFrame
        Score matrix with cells as rows and marker-set labels as columns.
    """
    raise NotImplementedError


def marker_score_metrics(
    adata: AnnData,
    markers: MarkerSets,
    *,
    label_key: str = "cell_type",
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
    ctrl_size: int = 50,
    gene_pool: Sequence[str] | None = None,
    n_bins: int = 25,
    random_state: int | None = 0,
    include_median_assigned_score: bool = True,
) -> pd.DataFrame:
    """Compute marker-set agreement metrics from per-cell marker scores.

    For each cell, identifies the highest-scoring marker set. For each annotated
    cell type, reports the fraction of cells whose assigned annotation is the
    top-scoring marker set and the median score margin against the best
    alternative marker set.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    label_key
        Column in ``adata.obs`` containing assigned annotations.
    layer
        Optional expression layer.
    use_raw
        Whether to use ``adata.raw``.
    gene_symbols_key
        Optional gene-symbol column in ``adata.var``.
    labels
        Optional subset/order of cell types to evaluate.
    ctrl_size
        Control gene count for Scanpy ``score_genes``.
    gene_pool
        Optional background gene pool.
    n_bins
        Number of expression bins for control gene matching.
    random_state
        Random seed for module score calculation.
    include_median_assigned_score
        Whether to include the median assigned marker score.

    Returns
    -------
    pandas.DataFrame
        Per-cell-type marker-score QC metrics.
    """
    raise NotImplementedError


def marker_de_recovery_metrics(
    adata: AnnData,
    markers: MarkerSets,
    *,
    label_key: str = "cell_type",
    layer: str | None = None,
    use_raw: bool | None = None,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
    method: Literal["wilcoxon", "t-test", "t-test_overestim_var"] = "wilcoxon",
    alpha: float = 0.05,
    logfc_threshold: float = 0.25,
    top_n: int = 50,
    rank_by: Literal["scores", "logfoldchanges", "pvals_adj"] = "scores",
    corr_method: Literal[
        "benjamini-hochberg",
        "bonferroni",
    ] = "benjamini-hochberg",
    tie_correct: bool = False,
    key_added: str = "rank_genes_groups_marker_qc",
    return_de_table: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Compute one-vs-rest differential-expression marker recovery metrics.

    Intended to run Scanpy ``rank_genes_groups`` using the annotation column as
    groups and ``reference="rest"``. For each cell type, calculates the fraction
    of provided markers recovered as significant upregulated DE genes and the
    recall of markers among the top N ranked DE genes.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    label_key
        Column in ``adata.obs`` containing cell type annotations.
    layer
        Optional expression layer for DE.
    use_raw
        Whether Scanpy should use ``adata.raw``.
    gene_symbols_key
        Optional gene-symbol column in ``adata.var``.
    labels
        Optional subset/order of cell types to evaluate.
    method
        Differential-expression test passed to Scanpy.
    alpha
        Adjusted p-value threshold for significant marker recovery.
    logfc_threshold
        Minimum log fold change for significant upregulated DE genes.
    top_n
        Number of top DE genes used for marker recall.
    rank_by
        Field used to rank DE genes for top-N marker recall.
    corr_method
        Multiple-testing correction method passed to Scanpy.
    tie_correct
        Whether to use tie correction for Wilcoxon.
    key_added
        Key used in ``adata.uns`` for Scanpy DE results.
    return_de_table
        If True, also return a long-form DE results table.

    Returns
    -------
    pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
        Per-cell-type DE recovery metrics, and optionally long-form DE results.
    """
    raise NotImplementedError


def marker_annotation_qc(
    adata: AnnData,
    markers: MarkerSets,
    *,
    label_key: str = "cell_type",
    layer: str | None = None,
    use_raw: bool = False,
    gene_symbols_key: str | None = None,
    labels: Sequence[str] | None = None,
    logfc_pseudocount: float = 1e-9,
    logfc_threshold: float = 0.25,
    score_ctrl_size: int = 50,
    score_gene_pool: Sequence[str] | None = None,
    score_n_bins: int = 25,
    score_random_state: int | None = 0,
    de_method: Literal["wilcoxon", "t-test", "t-test_overestim_var"] = "wilcoxon",
    de_alpha: float = 0.05,
    de_logfc_threshold: float = 0.25,
    de_top_n: int = 50,
    include_optional_score_metrics: bool = False,
) -> pd.DataFrame:
    """Compute the full marker-based annotation QC table.

    Combines marker logFC enrichment, marker-set score agreement, and
    one-vs-rest DE marker recovery into one compact per-cell-type table.

    Parameters
    ----------
    adata
        Annotated data matrix.
    markers
        Mapping from cell type name to marker gene names.
    label_key
        Column in ``adata.obs`` containing cell type annotations.
    layer
        Optional expression layer used across metrics.
    use_raw
        Whether to use ``adata.raw``.
    gene_symbols_key
        Optional gene-symbol column in ``adata.var``.
    labels
        Optional subset/order of cell types to evaluate.
    logfc_pseudocount
        Pseudocount for marker logFC calculations.
    logfc_threshold
        Threshold for enriched-marker fraction.
    score_ctrl_size
        Control gene count for marker-set scoring.
    score_gene_pool
        Optional background gene pool for marker-set scoring.
    score_n_bins
        Number of expression bins for marker-set scoring.
    score_random_state
        Random seed for marker-set scoring.
    de_method
        Differential-expression method.
    de_alpha
        Adjusted p-value threshold for significant marker recovery.
    de_logfc_threshold
        LogFC threshold for significant marker recovery.
    de_top_n
        Top-N DE genes used for marker recall.
    include_optional_score_metrics
        Whether to include optional score summaries such as median assigned
        score.

    Returns
    -------
    pandas.DataFrame
        Final marker-based annotation QC table.
    """
    raise NotImplementedError
