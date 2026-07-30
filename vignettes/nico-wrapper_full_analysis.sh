#!/bin/bash
#SBATCH --job-name=nico-wrapper_full
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=jobID_%j.log


# ── Input datasets ─────────────────────────────────────────────────────────────
REFERENCE_H5AD="<path_to_ref_data/scRNAseq_data.h5ad>"
SPATIAL_H5AD="<path_to_spatial_data/spatial_data.h5ad>"

# ── Key column names ───────────────────────────────────────────────────────────
# SPATIAL_KEY="spatial"   # .obsm key holding XY coordinates in the spatial file
# REF_LABEL_KEY="cluster" # .obs column holding cell-type labels in the reference data

# ── Cell types and factors of interest (used in covariation reports) ───────────
# Central (CC) = Stem/TA; Neighbor (NC) = Paneth; both use factor 1.
CC_NAME="Stem/TA"
NC_NAME="Paneth"
CC_FACTOR_ID=1
NC_FACTOR_ID=1


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 1 — Core pipeline
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 1a. Preprocess and normalize reference and spatial datasets
uv run nico-wrapper preprocess build \
  --reference $REFERENCE_H5AD \
  --spatial $SPATIAL_H5AD \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery

# 1b. Transfer cell-type labels from reference to spatial cells
uv run nico-wrapper transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis

# 1c. Identify spatial niche interactions
uv run nico-wrapper niche run \
  --output-dir nico_analysis

# 1d. Run niche covariation analysis
uv run nico-wrapper covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 2 — Export stable tables
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Export directed niche interaction coefficient table (TSV)
uv run nico-wrapper niche export \
  --output-dir nico_analysis \
  --cutoff 0.0 \
  --include-self-edges

# Export regression coefficient table (TSV)
uv run nico-wrapper covariation export \
  --output-dir nico_analysis \
  --kind regression


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 3 — Niche interaction plots
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Confusion matrix, coefficient matrix, evaluation scores, interaction graph,
# and top interaction coefficients per cell type (all at once)
uv run nico-wrapper niche plot \
  --output-dir nico_analysis \
  --kind confusion \
  --kind coefficients \
  --kind scores \
  --kind graph \
  --kind top-coefficients \
  --plot-format png \
  --interaction-cutoff 0.1


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 4 — Proximity analysis (observed vs. randomized co-localization)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

uv run nico-wrapper niche proximity \
  --output-dir nico_analysis \
  --n-permutations 1000 \
  --seed 42 \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 5 — Covariation reports (broad)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 5a. Regression circle plots + heatmaps for Stem/TA, and p-value sizebar legend
uv run nico-wrapper covariation reports \
  --output-dir nico_analysis \
  --kind regression-circleplots \
  --kind regression-heatmaps \
  --kind pvalue-sizebar \
  --cell-type "Stem/TA" \
  --pvalue-cutoff 0.05 \
  --plot-format png \
  --dpi 300

# 5b. Gene-correlation Excel workbook and LR summary workbook (all cell types)
uv run nico-wrapper covariation reports \
  --output-dir nico_analysis \
  --kind gene-correlation-excel \
  --kind lr-summary

# 5c. Factor-gene heatmaps (cosine + Spearman correlation) for Paneth
uv run nico-wrapper covariation reports \
  --output-dir nico_analysis \
  --kind factor-gene-heatmaps \
  --cell-type "Paneth" \
  --plot-format png \
  --dpi 300

# 5d. Top genes across all factors as dot plots for Paneth and Stem/TA
uv run nico-wrapper covariation reports \
  --output-dir nico_analysis \
  --kind top-genes-all-factors \
  --cell-type "Paneth" \
  --cell-type "Stem/TA" \
  --plot-format png \
  --dpi 150

# 5e. Weighted factor-neighborhood feature matrix
uv run nico-wrapper covariation feature-matrix \
  --output-dir nico_analysis \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 6 — Top genes for selected cell type / factor
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 6a. Extract top genes for Stem/TA factor 1 (negatively correlated; saves TSV + dotplot)
uv run nico-wrapper covariation top-genes \
  --output-dir nico_analysis \
  --cell-type "Stem/TA" \
  --factor-id $NC_FACTOR_ID \
  --top-n 20 \
  --negative \
  --plot-format png

# 6b. Paired top-gene dot plot: Stem/TA (factor 1) vs Paneth (factor 1)
uv run nico-wrapper covariation top-genes \
  --output-dir nico_analysis \
  --cell-type "Stem/TA" \
  --factor-id $NC_FACTOR_ID \
  --pair-cell-type "Paneth" \
  --pair-factor-id $CC_FACTOR_ID \
  --top-n 20 \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 7 — Factor UMAP plots
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 7a. Spatial UMAP: Stem/TA (factor 1) and Paneth (factor 1)
uv run nico-wrapper covariation umap \
  --output-dir nico_analysis \
  --modality spatial \
  --cell-type "Stem/TA" --factor-id $NC_FACTOR_ID \
  --cell-type "Paneth"  --factor-id $CC_FACTOR_ID \
  --plot-format png

# 7b. scRNA-seq UMAP: Stem/TA (factor 1) and Paneth (factor 1)
uv run nico-wrapper covariation umap \
  --output-dir nico_analysis \
  --modality sc \
  --cell-type "Stem/TA" --factor-id $NC_FACTOR_ID \
  --cell-type "Paneth"  --factor-id $CC_FACTOR_ID \
  --plot-format png

# 7c. Spatial UMAP: single cell type — Stem/TA, factor 1
uv run nico-wrapper covariation umap \
  --output-dir nico_analysis \
  --modality spatial \
  --cell-type "Stem/TA" --factor-id $NC_FACTOR_ID \
  --plot-format png

# 7d. scRNA-seq UMAP: single cell type — Stem/TA, factor 1
uv run nico-wrapper covariation umap \
  --output-dir nico_analysis \
  --modality sc \
  --cell-type "Stem/TA" --factor-id $NC_FACTOR_ID \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 8 — Colocalization factor plots (colocalized vs. non-colocalized cells)
# Central = Stem/TA (factor 1); Neighbor = Paneth (factor 1)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 8a. Default (spatial scatter only — no bar, no violin)
uv run nico-wrapper covariation colocalize \
  --output-dir nico_analysis \
  --central-cell-type "$CC_NAME" \
  --neighbor-cell-type "$NC_NAME" \
  --central-factor-id $CC_FACTOR_ID \
  --neighbor-factor-id $NC_FACTOR_ID \
  --no-bar \
  --no-violin \
  --plot-format png

# 8b. Bar plot: factor loadings of colocalized vs non-colocalized instances
uv run nico-wrapper covariation colocalize \
  --output-dir nico_analysis \
  --central-cell-type "$CC_NAME" \
  --neighbor-cell-type "$NC_NAME" \
  --central-factor-id $CC_FACTOR_ID \
  --neighbor-factor-id $NC_FACTOR_ID \
  --bar \
  --no-violin \
  --plot-format png

# 8c. Violin plot: factor loadings of colocalized vs non-colocalized instances
uv run nico-wrapper covariation colocalize \
  --output-dir nico_analysis \
  --central-cell-type "$CC_NAME" \
  --neighbor-cell-type "$NC_NAME" \
  --central-factor-id $CC_FACTOR_ID \
  --neighbor-factor-id $NC_FACTOR_ID \
  --no-bar \
  --violin \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 9 — Ligand-receptor interaction analysis
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 9a. Focused pair: Stem/TA → Paneth (central factor 1, neighbor factor 1)
uv run nico-wrapper covariation lr \
  --output-dir nico_analysis \
  --central-cell-type "Stem/TA" \
  --neighbor-cell-type "Paneth" \
  --central-factor-id $CC_FACTOR_ID \
  --neighbor-factor-id $NC_FACTOR_ID \
  --pvalue-cutoff 0.05 \
  --plot-format png

# 9b. All significant interaction partners of Paneth (no factor filter)
uv run nico-wrapper covariation lr \
  --output-dir nico_analysis \
  --central-cell-type "Paneth" \
  --pvalue-cutoff 0.05 \
  --plot-format png


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 10 — Pathway enrichment analysis (requires network access to Enrichr)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 10a. Stem/TA factor 1 — dotplot (positively correlated genes, top 50)
uv run nico-wrapper covariation pathway \
  --output-dir nico_analysis \
  --cell-type "Stem/TA" \
  --factor-id $CC_FACTOR_ID \
  --top-genes 50 \
  --database "GO_Biological_Process_2021" \
  --organism Mouse \
  --plot-as dotplot \
  --plot-format png

# 10b. Stem/TA factor 1 — barplot
uv run nico-wrapper covariation pathway \
  --output-dir nico_analysis \
  --cell-type "Stem/TA" \
  --factor-id $CC_FACTOR_ID \
  --top-genes 50 \
  --database "GO_Biological_Process_2021" \
  --organism Mouse \
  --plot-as barplot \
  --plot-format png
