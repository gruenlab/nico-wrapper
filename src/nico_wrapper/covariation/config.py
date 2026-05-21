"""Configuration objects for NiCo covariation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

DEFAULT_RIDGE_ALPHAS: tuple[float, ...] = tuple(float(x) for x in np.power(2.0, np.arange(-10, 10)))


@dataclass(frozen=True)
class CovariationConfig:
    """Configuration for core NiCo niche covariation analysis."""

    radius: int | float = 0
    n_factors: int = 3

    modality: Literal["double", "single"] = "double"
    factorization: Literal["inmf", "nmf-transfer"] = "inmf"

    ref_label_key: str = "cluster"
    annotated_h5ad_name: str = "nico_celltype_annotation.h5ad"
    ref_original_counts_name: str = "Original_counts.h5ad"
    ref_sct_name: str = "sct_singleCell.h5ad"
    spatial_sct_name: str = "sct_spatial.h5ad"

    ligand_receptor_db: Path | None = None

    ridge_alphas: tuple[float, ...] = DEFAULT_RIDGE_ALPHAS
    logistic_coef_cutoff: float = 0.0
    ridge_coef_cutoff: float = 0.0
    expression_population_cutoff: float = 0.0

    seed: int = 541
    shap_analysis: bool = False
    shap_cluster_cutoff: float = 0.5

    overwrite: bool = False
    write_manifest: bool = True
    persist_state: bool = True
    export_regression_table: bool = True


@dataclass(frozen=True)
class CovariationReportConfig:
    """Configuration for optional covariation report generation."""

    enabled: bool = False
    kinds: tuple[str, ...] = (
        "regression-circleplots",
        "regression-heatmaps",
        "gene-correlation-excel",
        "lr-summary",
    )

    saveas: str = "pdf"
    dpi: int = 300
    transparent: bool = False
    show: bool = False

    choose_celltypes: tuple[str, ...] = ()
    choose_factors_id: tuple[int, ...] = ()

    pvalue_cutoff: float = 0.05
    regression_plot_cutoff: float | None = None

    top_genes_per_factor: int = 30
    correlation_with_spearman: bool = True
    positively_correlated: bool = True
    include_rps_rpl_mt_genes: bool = True
    organism: Literal["Mouse", "Human"] = "Mouse"

    ligand_factor_threshold: float = 0.2
    receptor_factor_threshold: float = 0.2
    ligand_population_threshold: float = 0.2
    receptor_population_threshold: float = 0.2
    number_of_top_genes_to_print: int = 20

    pathway_databases: tuple[str, ...] = (
        "GO_Biological_Process_2021",
        "BioPlanet_2019",
        "Reactome_2016",
    )
    pathway_top_genes: int = 50
    pathway_plot_as: Literal["barplot", "dotplot"] = "barplot"
