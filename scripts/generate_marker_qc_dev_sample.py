#!/usr/bin/env python3
"""Generate marker-based QC tables for the bundled follicular dev sample.

Example:
    uv run python scripts/generate_marker_qc_dev_sample.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from nico_wrapper.qc import (
    load_marker_sets_json,
    marker_annotation_qc,
    marker_de_recovery_metrics,
    marker_logfc_metrics,
    marker_score_metrics,
    marker_set_scores,
    marker_set_summary,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata",
        type=Path,
        default=Path("dev-samples/sce_follicular.h5ad"),
        help="Input AnnData .h5ad file.",
    )
    parser.add_argument(
        "--markers",
        type=Path,
        default=Path("dev-samples/marker_genes.json"),
        help="Marker-set JSON file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("dev-samples/qc-output"),
        help="Directory for generated QC CSV files.",
    )
    parser.add_argument(
        "--label-key",
        default="celltype",
        help="adata.obs column containing annotations.",
    )
    parser.add_argument(
        "--layer",
        default="logcounts",
        help="Expression layer to use. Pass an empty string to use adata.X.",
    )
    parser.add_argument(
        "--de-method",
        default="wilcoxon",
        choices=["wilcoxon", "t-test", "t-test_overestim_var"],
        help="Scanpy rank_genes_groups method.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Top-N DE genes for marker recall.",
    )
    parser.add_argument(
        "--include-obs-score-columns",
        action="store_true",
        help="Also copy marker score columns into the loaded AnnData object before writing scores.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    layer = args.layer or None
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.adata)
    # The dev sample contains duplicate categorical observation and variable
    # names. Convert to plain string indexes before making them unique; unique
    # names avoid pandas/Scanpy alignment errors while preserving first exact
    # marker matches such as CD19, MS4A1, and TRAC.
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    markers = load_marker_sets_json(args.markers, allow_empty_marker_sets=True)
    labels = list(markers)

    summary = marker_set_summary(
        adata,
        markers,
        label_key=args.label_key,
        layer=layer,
        labels=labels,
    )
    summary.to_csv(args.output_dir / "marker_set_summary.csv", index=False)

    logfc_summary, logfc_per_marker = marker_logfc_metrics(
        adata,
        markers,
        label_key=args.label_key,
        layer=layer,
        labels=labels,
        return_per_marker=True,
    )
    logfc_summary.to_csv(args.output_dir / "marker_logfc_summary.csv", index=False)
    logfc_per_marker.to_csv(args.output_dir / "marker_logfc_per_marker.csv", index=False)

    scores = marker_set_scores(
        adata,
        markers,
        layer=layer,
        labels=labels,
        copy_scores_to_obs=args.include_obs_score_columns,
    )
    scores.to_csv(args.output_dir / "marker_scores.csv")

    score_metrics = marker_score_metrics(
        adata,
        markers,
        label_key=args.label_key,
        layer=layer,
        labels=labels,
    )
    score_metrics.to_csv(args.output_dir / "marker_score_metrics.csv", index=False)

    de_summary, de_table = marker_de_recovery_metrics(
        adata,
        markers,
        label_key=args.label_key,
        layer=layer,
        labels=labels,
        method=args.de_method,
        top_n=args.top_n,
        return_de_table=True,
    )
    de_summary.to_csv(args.output_dir / "marker_de_recovery_summary.csv", index=False)
    de_table.to_csv(args.output_dir / "marker_de_table.csv", index=False)

    final_qc = marker_annotation_qc(
        adata,
        markers,
        label_key=args.label_key,
        layer=layer,
        labels=labels,
        de_method=args.de_method,
        de_top_n=args.top_n,
    )
    final_qc.to_csv(args.output_dir / "marker_annotation_qc.csv", index=False)

    print(f"Wrote marker QC outputs to {args.output_dir}")
    print(final_qc.to_string(index=False))


if __name__ == "__main__":
    main()
