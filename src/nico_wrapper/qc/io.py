"""Input/output helpers for marker-based QC."""

from __future__ import annotations

import json
from os import PathLike
from pathlib import Path


def load_marker_sets_json(
    path: str | PathLike[str],
    *,
    encoding: str = "utf-8",
    strip_whitespace: bool = True,
    drop_empty_genes: bool = True,
    deduplicate: bool = True,
    allow_empty_marker_sets: bool = False,
) -> dict[str, list[str]]:
    """Load marker sets from a JSON file.

    Expects a top-level JSON object mapping cell type names to lists of marker
    gene names, for example ``{"B cell": ["MS4A1", "CD79A"]}``.

    Parameters
    ----------
    path
        Path to the marker-set JSON file.
    encoding
        File encoding used when reading the JSON file.
    strip_whitespace
        Whether to strip leading/trailing whitespace from cell type and gene
        names.
    drop_empty_genes
        Whether to drop empty gene names after optional whitespace stripping.
    deduplicate
        Whether to remove duplicate genes within each marker set while
        preserving their first-seen order.
    allow_empty_marker_sets
        Whether marker sets with no genes after parsing are allowed.

    Returns
    -------
    dict[str, list[str]]
        Parsed marker sets keyed by cell type name.

    Raises
    ------
    ValueError
        If the JSON structure does not match ``dict[str, list[str]]`` or if an
        empty marker set is encountered while ``allow_empty_marker_sets`` is
        False.
    """
    marker_path = Path(path)
    with marker_path.open("r", encoding=encoding) as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            "Marker-set JSON must contain a top-level object mapping cell type "
            "names to lists of marker gene names."
        )

    marker_sets: dict[str, list[str]] = {}
    for raw_cell_type, raw_genes in data.items():
        if not isinstance(raw_cell_type, str):
            raise ValueError(
                "Marker-set JSON keys must be strings; found key "
                f"{raw_cell_type!r} of type {type(raw_cell_type).__name__}."
            )

        cell_type = raw_cell_type.strip() if strip_whitespace else raw_cell_type
        if not cell_type:
            raise ValueError("Marker-set JSON contains an empty cell type name.")
        if cell_type in marker_sets:
            raise ValueError(
                "Marker-set JSON contains duplicate cell type names after "
                f"normalization: {cell_type!r}."
            )

        if not isinstance(raw_genes, list):
            raise ValueError(
                f"Marker set for {cell_type!r} must be a list of gene names; "
                f"found {type(raw_genes).__name__}."
            )

        genes: list[str] = []
        seen: set[str] = set()
        for index, raw_gene in enumerate(raw_genes):
            if not isinstance(raw_gene, str):
                raise ValueError(
                    f"Marker gene at index {index} for {cell_type!r} must be "
                    f"a string; found {type(raw_gene).__name__}."
                )

            gene = raw_gene.strip() if strip_whitespace else raw_gene
            if not gene and drop_empty_genes:
                continue
            if deduplicate:
                if gene in seen:
                    continue
                seen.add(gene)
            genes.append(gene)

        if not genes and not allow_empty_marker_sets:
            raise ValueError(f"Marker set for {cell_type!r} is empty.")

        marker_sets[cell_type] = genes

    return marker_sets
