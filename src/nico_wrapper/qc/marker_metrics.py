"""Marker-based annotation quality-control metric APIs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse

MarkerSets = Mapping[str, Sequence[str]]


@dataclass(frozen=True)
class ResolvedMarkerSet:
    """Marker-set resolution against one AnnData variable namespace."""

    provided: list[str]
    used_input_names: list[str]
    used_var_names: list[str]
    used_indices: list[int]
    missing: list[str]


# ---------------------------------------------------------------------------
# Validation and label helpers


def _validate_label_key(adata: AnnData, label_key: str) -> pd.Series:
    if label_key not in adata.obs:
        raise KeyError(
            f"Annotation label key {label_key!r} is not present in adata.obs. "
            f"Available columns include: {list(adata.obs.columns[:10])!r}"
        )
    return adata.obs[label_key]


def _selected_labels(markers: MarkerSets, labels: Sequence[str] | None) -> list[str]:
    selected = list(markers.keys()) if labels is None else list(labels)
    duplicates = [label for label in dict.fromkeys(selected) if selected.count(label) > 1]
    if duplicates:
        raise ValueError(f"labels contains duplicate entries: {duplicates!r}")

    missing = [label for label in selected if label not in markers]
    if missing:
        raise KeyError(f"No marker set provided for requested labels: {missing!r}")
    return selected


def _count_cells_by_label(obs_labels: pd.Series, labels: Sequence[str]) -> dict[str, int]:
    return {label: int((obs_labels == label).sum()) for label in labels}


# ---------------------------------------------------------------------------
# Expression source helpers


def _validate_expression_source(
    adata: AnnData,
    layer: str | None,
    use_raw: bool | None,
    *,
    reject_layer_with_raw: bool = True,
) -> None:
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True was requested, but adata.raw is None.")
        if layer is not None and reject_layer_with_raw:
            raise ValueError(
                "layer and use_raw=True cannot be used together because "
                "adata.raw has no layers."
            )
    if layer is not None and layer not in adata.layers:
        raise KeyError(
            f"Layer {layer!r} is not present in adata.layers. "
            f"Available layers: {list(adata.layers.keys())!r}"
        )


def _get_matrix_var_names_and_var(
    adata: AnnData,
    *,
    layer: str | None,
    use_raw: bool,
):
    _validate_expression_source(adata, layer, use_raw)
    if use_raw:
        # _validate_expression_source guarantees raw is present.
        return adata.raw.X, pd.Index(adata.raw.var_names), adata.raw.var
    if layer is not None:
        return adata.layers[layer], pd.Index(adata.var_names), adata.var
    return adata.X, pd.Index(adata.var_names), adata.var


def _effective_scanpy_use_raw(
    adata: AnnData,
    *,
    layer: str | None,
    use_raw: bool | None,
) -> bool:
    if layer is not None:
        if use_raw is True:
            raise ValueError("Scanpy cannot use both a layer and use_raw=True.")
        _validate_expression_source(adata, layer, False)
        return False

    if use_raw is None:
        return adata.raw is not None

    _validate_expression_source(adata, None, use_raw)
    return bool(use_raw)


# ---------------------------------------------------------------------------
# Marker and gene resolution helpers


def _deduplicate_preserve_order(values: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    deduplicated: list[str] = []
    for value in values:
        item = str(value)
        if item in seen:
            continue
        seen.add(item)
        deduplicated.append(item)
    return deduplicated


def _build_gene_lookup(
    var_names: pd.Index,
    var: pd.DataFrame,
    gene_symbols_key: str | None,
) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for index, var_name in enumerate(var_names):
        lookup.setdefault(str(var_name), index)

    if gene_symbols_key is None:
        return lookup

    if gene_symbols_key not in var:
        raise KeyError(
            f"gene_symbols_key {gene_symbols_key!r} is not present in var. "
            f"Available columns include: {list(var.columns[:10])!r}"
        )

    symbols = var[gene_symbols_key]
    symbol_counts = symbols[symbols.notna()].astype(str).value_counts()
    for index, raw_symbol in enumerate(symbols):
        if pd.isna(raw_symbol):
            continue
        symbol = str(raw_symbol)
        if symbol_counts.get(symbol, 0) != 1:
            continue
        # Keep exact var_name matches authoritative.
        lookup.setdefault(symbol, index)

    return lookup


def _resolve_marker_set(
    marker_genes: Sequence[str],
    lookup: Mapping[str, int],
    var_names: pd.Index,
) -> ResolvedMarkerSet:
    if isinstance(marker_genes, str):
        raise ValueError(
            "Marker sets must be sequences of gene names, not a single string. "
            f"Got {marker_genes!r}."
        )
    provided = _deduplicate_preserve_order(marker_genes)
    used_input_names: list[str] = []
    used_var_names: list[str] = []
    used_indices: list[int] = []
    missing: list[str] = []
    seen_indices: set[int] = set()

    for marker in provided:
        index = lookup.get(marker)
        if index is None:
            missing.append(marker)
            continue
        if index in seen_indices:
            # Two input names can resolve to the same variable (e.g. var_name and
            # symbol). Count/use that variable once for computation.
            continue
        seen_indices.add(index)
        used_input_names.append(marker)
        used_var_names.append(str(var_names[index]))
        used_indices.append(index)

    return ResolvedMarkerSet(
        provided=provided,
        used_input_names=used_input_names,
        used_var_names=used_var_names,
        used_indices=used_indices,
        missing=missing,
    )


def _resolve_marker_sets(
    markers: MarkerSets,
    labels: Sequence[str],
    var_names: pd.Index,
    var: pd.DataFrame,
    gene_symbols_key: str | None,
) -> dict[str, ResolvedMarkerSet]:
    lookup = _build_gene_lookup(var_names, var, gene_symbols_key)
    return {
        label: _resolve_marker_set(markers[label], lookup, var_names)
        for label in labels
    }


def _resolve_gene_pool(
    gene_pool: Sequence[str] | None,
    lookup: Mapping[str, int],
    var_names: pd.Index,
) -> list[str] | None:
    if gene_pool is None:
        return None

    resolved = _resolve_marker_set(gene_pool, lookup, var_names).used_var_names
    if not resolved:
        raise ValueError("gene_pool was supplied but none of its genes resolved.")
    return resolved


def _make_unique_index_strings(index: pd.Index, *, join: str = "-") -> pd.Index:
    """Return a string index with duplicate names made unique like AnnData."""
    used: set[str] = set()
    counts: dict[str, int] = {}
    values: list[str] = []
    for raw_value in index:
        value = str(raw_value)
        if value not in used:
            used.add(value)
            counts[value] = 0
            values.append(value)
            continue

        counts[value] = counts.get(value, 0) + 1
        candidate = f"{value}{join}{counts[value]}"
        while candidate in used:
            counts[value] += 1
            candidate = f"{value}{join}{counts[value]}"
        used.add(candidate)
        values.append(candidate)
    return pd.Index(values)


# ---------------------------------------------------------------------------
# Numeric helpers


def _column_means(matrix, row_mask: np.ndarray, col_indices: Sequence[int]) -> np.ndarray:
    if len(col_indices) == 0:
        return np.asarray([], dtype=float)

    mask = np.asarray(row_mask, dtype=bool)
    subset = matrix[mask, :][:, list(col_indices)]
    means = subset.mean(axis=0)
    if sparse.issparse(means):
        means = means.A
    return np.asarray(means, dtype=float).ravel()


def _safe_fraction(numerator: int, denominator: int) -> float:
    return np.nan if denominator == 0 else float(numerator) / float(denominator)


def _format_float_for_column(value: float) -> str:
    return f"{value:g}"


def _logfc_fraction_column(logfc_threshold: float) -> str:
    return f"fraction_markers_logFC_gt_{_format_float_for_column(logfc_threshold)}"


def _topn_recall_column(top_n: int) -> str:
    return f"marker_recall_top{top_n}_DE"


def _empty_de_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "cell_type",
            "rank",
            "names",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "is_marker_for_group",
            "is_significant_upregulated",
        ]
    )


# ---------------------------------------------------------------------------
# Private computation helpers used by public APIs


def _marker_set_summary_from_resolved(
    labels: Sequence[str],
    obs_labels: pd.Series,
    resolved: Mapping[str, ResolvedMarkerSet],
) -> pd.DataFrame:
    counts = _count_cells_by_label(obs_labels, labels)
    rows = []
    for label in labels:
        marker_set = resolved[label]
        rows.append(
            {
                "cell_type": label,
                "n_cells": counts[label],
                "n_markers_provided": len(marker_set.provided),
                "n_markers_used": len(marker_set.used_var_names),
                "missing_marker_names": list(marker_set.missing),
            }
        )
    return pd.DataFrame(rows)


def _marker_logfc_metrics_from_resolved(
    matrix,
    labels: Sequence[str],
    obs_labels: pd.Series,
    resolved: Mapping[str, ResolvedMarkerSet],
    *,
    pseudocount: float,
    logfc_threshold: float,
    min_cells: int,
    return_per_marker: bool,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    if pseudocount <= 0:
        raise ValueError("pseudocount must be greater than 0.")
    if min_cells < 1:
        raise ValueError("min_cells must be at least 1.")

    threshold_column = _logfc_fraction_column(logfc_threshold)
    summary_rows: list[dict[str, object]] = []
    per_marker_rows: list[dict[str, object]] = []

    for label in labels:
        target_mask = np.asarray(obs_labels == label, dtype=bool)
        rest_mask = ~target_mask
        n_cells = int(target_mask.sum())
        n_rest_cells = int(rest_mask.sum())
        marker_set = resolved[label]
        n_markers_used = len(marker_set.used_indices)

        means = np.full(n_markers_used, np.nan, dtype=float)
        rest_means = np.full(n_markers_used, np.nan, dtype=float)
        marker_logfc = np.full(n_markers_used, np.nan, dtype=float)

        evaluable = (
            n_cells >= min_cells
            and n_rest_cells > 0
            and n_markers_used > 0
        )
        if evaluable:
            means = _column_means(matrix, target_mask, marker_set.used_indices)
            rest_means = _column_means(matrix, rest_mask, marker_set.used_indices)
            marker_logfc = np.log2((means + pseudocount) / (rest_means + pseudocount))
            mean_logfc = float(np.nanmean(marker_logfc))
            fraction = float(np.sum(marker_logfc > logfc_threshold) / n_markers_used)
        else:
            mean_logfc = np.nan
            fraction = np.nan

        summary_rows.append(
            {
                "cell_type": label,
                "n_cells": n_cells,
                "n_rest_cells": n_rest_cells,
                "n_markers_provided": len(marker_set.provided),
                "n_markers_used": n_markers_used,
                "mean_marker_logFC": mean_logfc,
                threshold_column: fraction,
            }
        )

        if return_per_marker:
            for marker, var_name, mean, rest_mean, logfc in zip(
                marker_set.used_input_names,
                marker_set.used_var_names,
                means,
                rest_means,
                marker_logfc,
                strict=True,
            ):
                per_marker_rows.append(
                    {
                        "cell_type": label,
                        "marker": marker,
                        "var_name": var_name,
                        "mean_expression": float(mean) if not np.isnan(mean) else np.nan,
                        "rest_mean_expression": (
                            float(rest_mean) if not np.isnan(rest_mean) else np.nan
                        ),
                        "marker_logFC": float(logfc) if not np.isnan(logfc) else np.nan,
                        "is_logFC_gt_threshold": (
                            bool(logfc > logfc_threshold)
                            if not np.isnan(logfc)
                            else np.nan
                        ),
                    }
                )

    summary = pd.DataFrame(summary_rows)
    if not return_per_marker:
        return summary
    per_marker = pd.DataFrame(
        per_marker_rows,
        columns=[
            "cell_type",
            "marker",
            "var_name",
            "mean_expression",
            "rest_mean_expression",
            "marker_logFC",
            "is_logFC_gt_threshold",
        ],
    )
    return summary, per_marker


def _compute_marker_scores(
    adata: AnnData,
    markers: MarkerSets,
    labels: Sequence[str],
    *,
    layer: str | None,
    use_raw: bool,
    gene_symbols_key: str | None,
    ctrl_size: int,
    gene_pool: Sequence[str] | None,
    n_bins: int,
    random_state: int | None,
    score_prefix: str,
    copy_scores_to_obs: bool,
) -> pd.DataFrame:
    import scanpy as sc

    _validate_expression_source(adata, layer, use_raw)
    _, var_names, var = _get_matrix_var_names_and_var(adata, layer=layer, use_raw=use_raw)
    resolved = _resolve_marker_sets(markers, labels, var_names, var, gene_symbols_key)
    lookup = _build_gene_lookup(var_names, var, gene_symbols_key)
    resolved_gene_pool_indices: list[int] | None = None
    if gene_pool is not None:
        resolved_gene_pool_set = _resolve_marker_set(gene_pool, lookup, var_names)
        resolved_gene_pool_indices = resolved_gene_pool_set.used_indices
        if not resolved_gene_pool_indices:
            raise ValueError("gene_pool was supplied but none of its genes resolved.")

    score_use_raw = use_raw
    score_layer = layer
    if use_raw:
        # Scanpy's score_genes requires a unique gene index. raw.var_names cannot
        # be safely changed in-place, so score against a temporary raw AnnData.
        work_adata = adata.raw.to_adata()
        work_adata.obs = adata.obs.copy()
        score_use_raw = False
        score_layer = None
    elif copy_scores_to_obs and var_names.is_unique:
        work_adata = adata
    else:
        work_adata = adata.copy()

    score_var_names = pd.Index(work_adata.var_names)
    if not score_var_names.is_unique:
        score_var_names = _make_unique_index_strings(score_var_names)
        work_adata.var_names = score_var_names

    resolved_gene_pool = (
        None
        if resolved_gene_pool_indices is None
        else [str(score_var_names[index]) for index in resolved_gene_pool_indices]
    )
    scores = pd.DataFrame(index=adata.obs_names)

    for label in labels:
        marker_set = resolved[label]
        score_name = f"{score_prefix}{label}"
        if not marker_set.used_indices:
            score_values = np.full(adata.n_obs, np.nan, dtype=float)
            if copy_scores_to_obs:
                adata.obs[score_name] = score_values
        else:
            gene_list = [str(score_var_names[index]) for index in marker_set.used_indices]
            try:
                sc.tl.score_genes(
                    work_adata,
                    gene_list=gene_list,
                    score_name=score_name,
                    ctrl_size=ctrl_size,
                    gene_pool=resolved_gene_pool,
                    n_bins=n_bins,
                    random_state=random_state,
                    use_raw=score_use_raw,
                    layer=score_layer,
                    copy=False,
                )
            except RuntimeError as error:
                if "No control genes found" not in str(error):
                    raise
                score_values = np.full(adata.n_obs, np.nan, dtype=float)
                if copy_scores_to_obs:
                    adata.obs[score_name] = score_values
            else:
                score_values = work_adata.obs[score_name].to_numpy(dtype=float)
                if copy_scores_to_obs and work_adata is not adata:
                    adata.obs[score_name] = score_values
        scores[label] = score_values

    return scores


def _marker_score_metrics_from_scores(
    scores: pd.DataFrame,
    labels: Sequence[str],
    obs_labels: pd.Series,
    *,
    include_median_assigned_score: bool,
) -> pd.DataFrame:
    numeric_scores = scores.loc[:, list(labels)].apply(pd.to_numeric, errors="coerce")

    non_nan_rows = numeric_scores.notna().any(axis=1)
    top_labels = pd.Series(pd.NA, index=numeric_scores.index, dtype="object")
    if non_nan_rows.any():
        top_labels.loc[non_nan_rows] = numeric_scores.loc[non_nan_rows].idxmax(axis=1)

    rows: list[dict[str, object]] = []
    for label in labels:
        cell_mask = np.asarray(obs_labels == label, dtype=bool)
        n_cells = int(cell_mask.sum())
        row: dict[str, object] = {"cell_type": label, "n_cells": n_cells}

        assigned_unavailable = (
            n_cells == 0
            or label not in numeric_scores
            or numeric_scores.loc[cell_mask, label].isna().all()
        )
        if assigned_unavailable:
            row["correct_top_score_fraction"] = np.nan
            row["median_score_margin"] = np.nan
            if include_median_assigned_score:
                row["median_assigned_marker_score"] = np.nan
            rows.append(row)
            continue

        assigned_scores = numeric_scores.loc[cell_mask, label]
        row["correct_top_score_fraction"] = float(
            (top_labels.loc[cell_mask] == label).sum() / n_cells
        )

        if len(labels) <= 1:
            row["median_score_margin"] = np.nan
        else:
            alternative_labels = [other for other in labels if other != label]
            alternative_scores = numeric_scores.loc[cell_mask, alternative_labels]
            best_alternative = alternative_scores.max(axis=1, skipna=True)
            margins = assigned_scores - best_alternative
            row["median_score_margin"] = (
                float(np.nanmedian(margins.to_numpy(dtype=float)))
                if margins.notna().any()
                else np.nan
            )

        if include_median_assigned_score:
            row["median_assigned_marker_score"] = (
                float(np.nanmedian(assigned_scores.to_numpy(dtype=float)))
                if assigned_scores.notna().any()
                else np.nan
            )
        rows.append(row)

    return pd.DataFrame(rows)


def _de_var_namespace(
    adata: AnnData,
    *,
    layer: str | None,
    effective_use_raw: bool,
) -> tuple[pd.Index, pd.DataFrame]:
    if effective_use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True was requested, but adata.raw is None.")
        return pd.Index(adata.raw.var_names), adata.raw.var
    _validate_expression_source(adata, layer, False)
    return pd.Index(adata.var_names), adata.var


def _extract_rank_genes_groups_df(
    work_adata: AnnData,
    groups: Sequence[str],
    key_added: str,
) -> pd.DataFrame:
    import scanpy as sc

    tables: list[pd.DataFrame] = []
    for group in groups:
        group_df = sc.get.rank_genes_groups_df(work_adata, group=group, key=key_added)
        group_df = group_df.copy()
        if "group" in group_df:
            group_df = group_df.rename(columns={"group": "cell_type"})
        if "cell_type" not in group_df:
            group_df.insert(0, "cell_type", group)
        tables.append(group_df)
    if not tables:
        return _empty_de_table()
    return pd.concat(tables, ignore_index=True)


def _marker_de_recovery_metrics_from_resolved(
    adata: AnnData,
    labels: Sequence[str],
    obs_labels: pd.Series,
    resolved: Mapping[str, ResolvedMarkerSet],
    *,
    layer: str | None,
    effective_use_raw: bool,
    gene_symbols_key: str | None,
    source_var: pd.DataFrame,
    method: Literal["wilcoxon", "t-test", "t-test_overestim_var"],
    alpha: float,
    logfc_threshold: float,
    top_n: int,
    rank_by: Literal["scores", "logfoldchanges", "pvals_adj"],
    corr_method: Literal["benjamini-hochberg", "bonferroni"],
    tie_correct: bool,
    key_added: str,
    return_de_table: bool,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    import scanpy as sc

    if not 0 <= alpha <= 1:
        raise ValueError("alpha must be between 0 and 1, inclusive.")
    if top_n < 1:
        raise ValueError("top_n must be at least 1.")
    if rank_by not in {"scores", "logfoldchanges", "pvals_adj"}:
        raise ValueError("rank_by must be one of 'scores', 'logfoldchanges', or 'pvals_adj'.")

    counts = _count_cells_by_label(obs_labels, labels)
    recall_column = _topn_recall_column(top_n)
    observed_groups = pd.Index(obs_labels.dropna().astype(str).unique())
    total_cells = int(adata.n_obs)
    present_selected_labels = [
        label
        for label in labels
        if counts[label] >= 2 and total_cells - counts[label] >= 2
    ]

    de_df = _empty_de_table()
    de_var_names = pd.Index(source_var.index.astype(str))
    if len(observed_groups) >= 2 and present_selected_labels:
        if effective_use_raw:
            work_adata = adata.raw.to_adata()
            work_adata.obs = adata.obs.copy()
            de_use_raw = False
            de_layer = None
        else:
            work_adata = adata.copy()
            de_use_raw = False
            de_layer = layer

        de_var_names = pd.Index(work_adata.var_names.astype(str))
        if not de_var_names.is_unique:
            de_var_names = _make_unique_index_strings(de_var_names)
            work_adata.var_names = de_var_names

        sc.tl.rank_genes_groups(
            work_adata,
            groupby=obs_labels.name,
            groups=present_selected_labels,
            reference="rest",
            method=method,
            corr_method=corr_method,
            tie_correct=tie_correct,
            key_added=key_added,
            use_raw=de_use_raw,
            layer=de_layer,
            copy=False,
        )
        de_df = _extract_rank_genes_groups_df(work_adata, present_selected_labels, key_added)
        de_df = de_df.reset_index(drop=True)
        de_df["rank"] = de_df.groupby("cell_type", sort=False).cumcount() + 1

        if gene_symbols_key is not None and gene_symbols_key in source_var:
            symbol_lookup = {
                str(var_name): source_var[gene_symbols_key].iloc[index]
                for index, var_name in enumerate(de_var_names)
                if pd.notna(source_var[gene_symbols_key].iloc[index])
            }
            de_df[gene_symbols_key] = de_df["names"].astype(str).map(symbol_lookup)

        marker_lookup = {
            label: {str(de_var_names[index]) for index in resolved[label].used_indices}
            for label in labels
        }
        de_df["is_marker_for_group"] = [
            str(name) in marker_lookup.get(str(cell_type), set())
            for cell_type, name in zip(de_df["cell_type"], de_df["names"], strict=False)
        ]
        de_df["is_significant_upregulated"] = (
            (pd.to_numeric(de_df.get("pvals_adj"), errors="coerce") < alpha)
            & (pd.to_numeric(de_df.get("logfoldchanges"), errors="coerce") > logfc_threshold)
        )

    summary_rows: list[dict[str, object]] = []
    for label in labels:
        marker_set = resolved[label]
        n_markers_used = len(marker_set.used_var_names)
        group_de = de_df.loc[de_df["cell_type"].astype(str) == str(label)].copy()

        if n_markers_used == 0 or group_de.empty:
            fraction_significant = np.nan
            topn_recall = np.nan
        else:
            markers_used = set(marker_set.used_var_names)
            significant_genes = set(
                group_de.loc[group_de["is_significant_upregulated"], "names"].astype(str)
            )
            fraction_significant = _safe_fraction(
                len(markers_used.intersection(significant_genes)),
                n_markers_used,
            )

            ascending = rank_by == "pvals_adj"
            sorted_group = group_de.sort_values(
                by=rank_by,
                ascending=ascending,
                na_position="last",
                kind="mergesort",
            )
            top_genes = set(sorted_group.head(top_n)["names"].astype(str))
            topn_recall = _safe_fraction(len(markers_used.intersection(top_genes)), n_markers_used)

        summary_rows.append(
            {
                "cell_type": label,
                "n_cells": counts[label],
                "n_markers_provided": len(marker_set.provided),
                "n_markers_used": n_markers_used,
                "fraction_markers_significant_DE": fraction_significant,
                recall_column: topn_recall,
            }
        )

    summary = pd.DataFrame(summary_rows)
    if return_de_table:
        return summary, de_df.reset_index(drop=True)
    return summary


# ---------------------------------------------------------------------------
# Public API


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
    """Summarize marker availability and annotation cell counts.

    Returns one row per selected marker-set label with ``cell_type``,
    ``n_cells``, ``n_markers_provided``, ``n_markers_used``, and
    ``missing_marker_names``. Marker genes are resolved against the selected
    expression source and output rows preserve selected label order.
    """
    obs_labels = _validate_label_key(adata, label_key)
    selected = _selected_labels(markers, labels)
    _, var_names, var = _get_matrix_var_names_and_var(adata, layer=layer, use_raw=use_raw)
    resolved = _resolve_marker_sets(markers, selected, var_names, var, gene_symbols_key)
    return _marker_set_summary_from_resolved(selected, obs_labels, resolved)


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
    """Compute one-vs-rest marker expression enrichment per cell type.

    For each selected label, compares target-cell and rest-cell marker means and
    reports mean marker log2 fold-change plus the fraction above
    ``logfc_threshold``. If requested, also returns one row per resolved marker.
    Unevaluable comparisons are represented with ``NaN`` metrics.
    """
    obs_labels = _validate_label_key(adata, label_key)
    selected = _selected_labels(markers, labels)
    matrix, var_names, var = _get_matrix_var_names_and_var(adata, layer=layer, use_raw=use_raw)
    resolved = _resolve_marker_sets(markers, selected, var_names, var, gene_symbols_key)
    return _marker_logfc_metrics_from_resolved(
        matrix,
        selected,
        obs_labels,
        resolved,
        pseudocount=pseudocount,
        logfc_threshold=logfc_threshold,
        min_cells=min_cells,
        return_per_marker=return_per_marker,
    )


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
    """Compute Scanpy marker-set/module scores for each marker set.

    Returns a cell-by-marker-set score matrix indexed by ``adata.obs_names``.
    The original AnnData is not mutated unless ``copy_scores_to_obs=True``, in
    which case score columns named ``f"{score_prefix}{label}"`` are written to
    ``adata.obs``.
    """
    selected = _selected_labels(markers, labels)
    return _compute_marker_scores(
        adata,
        markers,
        selected,
        layer=layer,
        use_raw=use_raw,
        gene_symbols_key=gene_symbols_key,
        ctrl_size=ctrl_size,
        gene_pool=gene_pool,
        n_bins=n_bins,
        random_state=random_state,
        score_prefix=score_prefix,
        copy_scores_to_obs=copy_scores_to_obs,
    )


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
    """Summarize marker-score agreement with assigned annotations.

    Computes marker-set scores internally, identifies each cell's top-scoring
    selected marker set, and reports per-label top-score agreement and median
    assigned-vs-best-alternative score margin. The original AnnData is not
    modified.
    """
    obs_labels = _validate_label_key(adata, label_key)
    selected = _selected_labels(markers, labels)
    scores = _compute_marker_scores(
        adata,
        markers,
        selected,
        layer=layer,
        use_raw=use_raw,
        gene_symbols_key=gene_symbols_key,
        ctrl_size=ctrl_size,
        gene_pool=gene_pool,
        n_bins=n_bins,
        random_state=random_state,
        score_prefix="__marker_qc_score__",
        copy_scores_to_obs=False,
    )
    return _marker_score_metrics_from_scores(
        scores,
        selected,
        obs_labels,
        include_median_assigned_score=include_median_assigned_score,
    )


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

    Runs Scanpy ``rank_genes_groups`` on a temporary AnnData object and reports
    per-label fractions of resolved markers recovered as significant
    upregulated DE genes and among the top ``top_n`` ranked DE genes. Optionally
    returns the annotated long-form DE table.
    """
    obs_labels = _validate_label_key(adata, label_key)
    selected = _selected_labels(markers, labels)
    effective_use_raw = _effective_scanpy_use_raw(adata, layer=layer, use_raw=use_raw)
    var_names, var = _de_var_namespace(adata, layer=layer, effective_use_raw=effective_use_raw)
    resolved = _resolve_marker_sets(markers, selected, var_names, var, gene_symbols_key)
    return _marker_de_recovery_metrics_from_resolved(
        adata,
        selected,
        obs_labels,
        resolved,
        layer=layer,
        effective_use_raw=effective_use_raw,
        gene_symbols_key=gene_symbols_key,
        source_var=var,
        method=method,
        alpha=alpha,
        logfc_threshold=logfc_threshold,
        top_n=top_n,
        rank_by=rank_by,
        corr_method=corr_method,
        tie_correct=tie_correct,
        key_added=key_added,
        return_de_table=return_de_table,
    )


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
    """Compute the compact marker-based annotation QC table.

    Composes the private marker availability, logFC, marker-score, and DE
    recovery computations into the final per-cell-type table recommended by
    ``METRIC_QC.md``. Output order follows ``labels`` or marker mapping order.
    """
    obs_labels = _validate_label_key(adata, label_key)
    selected = _selected_labels(markers, labels)
    matrix, var_names, var = _get_matrix_var_names_and_var(adata, layer=layer, use_raw=use_raw)
    resolved = _resolve_marker_sets(markers, selected, var_names, var, gene_symbols_key)

    summary = _marker_set_summary_from_resolved(selected, obs_labels, resolved).drop(
        columns=["missing_marker_names"]
    )
    logfc = _marker_logfc_metrics_from_resolved(
        matrix,
        selected,
        obs_labels,
        resolved,
        pseudocount=logfc_pseudocount,
        logfc_threshold=logfc_threshold,
        min_cells=1,
        return_per_marker=False,
    ).drop(columns=["n_cells", "n_markers_provided", "n_markers_used"])
    scores = _compute_marker_scores(
        adata,
        markers,
        selected,
        layer=layer,
        use_raw=use_raw,
        gene_symbols_key=gene_symbols_key,
        ctrl_size=score_ctrl_size,
        gene_pool=score_gene_pool,
        n_bins=score_n_bins,
        random_state=score_random_state,
        score_prefix="__marker_qc_score__",
        copy_scores_to_obs=False,
    )
    score_metrics = _marker_score_metrics_from_scores(
        scores,
        selected,
        obs_labels,
        include_median_assigned_score=include_optional_score_metrics,
    ).drop(columns=["n_cells"])

    effective_use_raw = _effective_scanpy_use_raw(adata, layer=layer, use_raw=use_raw)
    de_var_names, de_var = _de_var_namespace(
        adata,
        layer=layer,
        effective_use_raw=effective_use_raw,
    )
    de_resolved = (
        resolved
        if de_var_names.equals(var_names)
        else _resolve_marker_sets(markers, selected, de_var_names, de_var, gene_symbols_key)
    )
    de_metrics = _marker_de_recovery_metrics_from_resolved(
        adata,
        selected,
        obs_labels,
        de_resolved,
        layer=layer,
        effective_use_raw=effective_use_raw,
        gene_symbols_key=gene_symbols_key,
        source_var=de_var,
        method=de_method,
        alpha=de_alpha,
        logfc_threshold=de_logfc_threshold,
        top_n=de_top_n,
        rank_by="scores",
        corr_method="benjamini-hochberg",
        tie_correct=False,
        key_added="rank_genes_groups_marker_qc",
        return_de_table=False,
    ).drop(columns=["n_cells", "n_markers_provided", "n_markers_used"])

    result = summary.merge(logfc, on="cell_type", how="left")
    result = result.merge(score_metrics, on="cell_type", how="left")
    result = result.merge(de_metrics, on="cell_type", how="left")

    columns = [
        "cell_type",
        "n_cells",
        "n_markers_provided",
        "n_markers_used",
        "mean_marker_logFC",
        _logfc_fraction_column(logfc_threshold),
        "correct_top_score_fraction",
        "median_score_margin",
    ]
    if include_optional_score_metrics:
        columns.append("median_assigned_marker_score")
    columns.extend(
        [
            "fraction_markers_significant_DE",
            _topn_recall_column(de_top_n),
        ]
    )
    return result.loc[:, columns]
