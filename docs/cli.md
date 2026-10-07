# CLI usage

From PyPI, run the umbrella CLI directly with `uvx`:

```bash
uvx --python 3.12 nico-wrapper --help
uvx --python 3.12 nico-wrapper preprocess --help
uvx --python 3.12 nico-wrapper transfer --help
uvx --python 3.12 nico-wrapper niche --help
uvx --python 3.12 nico-wrapper covariation --help
```

If your default interpreter is already Python 3.12, you can omit `--python 3.12`.

From the repository root, use `uv run` unless your `.venv` is activated:

```bash
uv run nico-wrapper --help
uv run nico-preprocess --help
uv run nico-transfer --help
uv run nico-niche --help
uv run nico-covariation --help
```

The umbrella and standalone interfaces execute the same Typer applications. The standalone command names remain supported.

## Overall workflow

The diagram shows the intended high-level flow. Initial user inputs appear at the start; intermediate NiCo artifacts are summarized only by major workflow stage.

```text
Initial inputs
--------------
- Reference scRNA-seq counts / AnnData
- Reference cell labels, default obs["cluster"]
- Spatial or Xenium counts / AnnData
- Spatial coordinates, default obsm["spatial"]

        |
        v
+-------------------------------+
| nico-preprocess               |
| optional conversion + build   |
+-------------------------------+
        |
        v
+-------------------------------+
| nico-transfer run             |
| reference labels -> spatial   |
+-------------------------------+
        |
        v
+-------------------------------+
| nico-niche run                |
| spatial cell-type interaction |
+-------------------------------+
        |
        v
+-------------------------------+
| nico-covariation run          |
| latent-factor covariation     |
+-------------------------------+
        |
        v
Final outputs
-------------
- Annotated spatial AnnData
- NiCo niche interaction artifacts
- Covariation factors, regression outputs, and optional reports
- Wrapper manifests and TSV sidecars where available
```

## General patterns

Use either the umbrella interface:

```bash
uvx --python 3.12 nico-wrapper <group> <command> [OPTIONS]
```

or an installed/source-checkout standalone command:

```bash
uv run <standalone-command> <command> [OPTIONS]
```

The command groups map directly onto the existing standalone commands:

| Umbrella command | Standalone command | Purpose |
|---|---|---|
| `nico-wrapper preprocess` | `nico-preprocess` | Prepare NiCo-ready input files |
| `nico-wrapper transfer` | `nico-transfer` | Transfer reference labels onto spatial cells |
| `nico-wrapper niche` | `nico-niche` | Run and inspect niche interaction analysis |
| `nico-wrapper covariation` | `nico-covariation` | Run and inspect covariation analysis and reports |

Default paths used throughout the examples:

```text
inputRef/       preprocessed reference files
inputQuery/     preprocessed spatial/query files
nico_analysis/  transfer, niche, and covariation outputs
```

## Recommended minimal workflow

If your reference and spatial inputs are already raw `.h5ad` files:

```bash
uv run nico-preprocess build \
  --reference data/reference_raw.h5ad \
  --spatial data/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery

uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis

uv run nico-niche run \
  --output-dir nico_analysis

uv run nico-covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

If your inputs are CSV/sparse files, run the relevant `nico-preprocess convert-*` command first.

# `nico-preprocess`

Preprocessing has two optional conversion commands and one required build command.

```text
convert-reference-sparse   optional sparse reference files -> raw reference h5ad
convert-spatial-csv        optional spatial count/coordinate CSV -> raw spatial h5ad
build                      raw h5ad inputs -> NiCo-ready outputs
```

## `nico-preprocess convert-reference-sparse`

Converts reference scRNA-seq sparse triplet files into a raw `.h5ad` file.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--counts` | yes | none | Sparse triplet count table. |
| `--genes` | yes | none | Gene metadata CSV. |
| `--barcodes` | yes | none | Cell barcode CSV. |
| `--output` | yes | none | Output raw reference `.h5ad`. |
| `--sep` | no | space | Separator used by the sparse count table. |
| `--header` | no | `1` | Header row passed to pandas for the count table. |
| `--one-based-indices / --zero-based-indices` | no | one-based | Whether sparse row/column indices are one-based. |

### Minimal example

```bash
uv run nico-preprocess convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad
```

### Outputs

```text
work/reference_raw.h5ad
```

The resulting AnnData stores cells in `.obs_names`, genes in `.var_names`, and counts in `.X`.

## `nico-preprocess convert-spatial-csv`

Converts a spatial count matrix plus coordinate table into a raw spatial `.h5ad` file.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--counts` | yes | none | Spatial count CSV. |
| `--coordinates` | yes | none | Spatial coordinate CSV. |
| `--output` | yes | none | Output raw spatial `.h5ad`. |
| `--counts-orientation` | no | `genes_by_cells` | Use `genes_by_cells` or `cells_by_genes`. |
| `--barcode-col` | no | `0` | Coordinate barcode column name or index. |
| `--coordinate-col` | no | all non-barcode columns | Coordinate column name/index. May be repeated. |
| `--spatial-key` | no | `spatial` | `.obsm` key used for coordinates. |
| `--reorder-coordinates / --no-reorder-coordinates` | no | reorder | Reorder coordinates to match count-matrix cells. |

### Minimal example

```bash
uv run nico-preprocess convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad
```

### Outputs

```text
work/spatial_raw.h5ad
```

Coordinates are written to `.obsm["spatial"]` by default.

## `nico-preprocess build`

Builds NiCo-ready reference and spatial files from raw `.h5ad` inputs.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--reference` | yes | none | Raw reference `.h5ad`. |
| `--spatial` | yes | none | Raw spatial/query `.h5ad`. |
| `--ref-out-dir` | yes | none | Output directory for reference NiCo files. |
| `--spatial-out-dir` | yes | none | Output directory for spatial NiCo files. |
| `--normalization` | no | `nico-sctransform` | `nico-sctransform` or `pearson-residuals`. |
| `--spatial-key` | no | `spatial` | Spatial coordinate key in `.obsm`. |
| `--ref-label-key` | no | `cluster` | Reference label column in `.obs`. |
| `--ref_min-cell-counts` | no | `5` | Minimum total counts per retained cell in reference data. |
| `--min-cell-counts` | no | `5` | Minimum total counts per retained cell in spatial data. |
| `--ref_min-gene-cells` | no | `1` | Minimum cells per retained gene in reference data. |
| `--spatial_min-gene-cells` | no | `1` | Minimum cells per retained gene in spatial data. |
| `--gene-space` | no | `shared` | `shared` or `reference_all`. |
| `--spatial-n-pcs` | no | `30` | PCs for spatial neighbor graph construction. |
| `--leiden-resolution` | no | `0.4`, `0.5` | Leiden resolution. May be repeated. |
| `--make-reference-umap / --no-make-reference-umap` | no | make | Compute reference UMAP in `Original_counts.h5ad`. |
| `--random-state` | no | `0` | Random seed for Scanpy steps. |
| `--overwrite / --no-overwrite` | no | no overwrite | Replace existing outputs. |
| `--sct-min-cells` | no | `1` | SCTransform-only minimum cells per gene. |
| `--sct-gmean-eps` | no | `1.0` | SCTransform geometric mean epsilon. |
| `--sct-n-genes` | no | `500` | SCTransform genes sampled for fitting. |
| `--sct-n-cells` | no | all cells | SCTransform cells sampled for fitting. |
| `--sct-bin-size` | no | `500` | SCTransform genes per fitting bin. |
| `--sct-bw-adjust` | no | `3.0` | SCTransform bandwidth adjustment. |
| `--pearson-theta` | no | `100.0` | Pearson residual overdispersion theta. |
| `--pearson-clip` | no | Scanpy default | Pearson residual clipping threshold. |
| `--pearson-check-values / --no-pearson-check-values` | no | check | Validate Pearson residual count values. |
| `--pearson-layer` | no | `.X` | Layer to normalize instead of `.X`. |

### Minimal example

```bash
uv run nico-preprocess build \
  --reference work/reference_raw.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery
```

### Outputs

```text
inputRef/Original_counts.h5ad
inputRef/sct_singleCell.h5ad
inputQuery/sct_spatial.h5ad
```

`nico-transfer` expects these files by default.

# `nico-transfer`

`nico-transfer` wraps NiCo anchor discovery and label propagation.

## `nico-transfer run`

Transfers labels from preprocessed reference cells onto preprocessed spatial/query cells.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--ref-dir` | yes | none | Directory containing `Original_counts.h5ad` and `sct_singleCell.h5ad`. |
| `--spatial-dir` | yes | none | Directory containing `sct_spatial.h5ad`. |
| `--output-dir` | yes | none | Output directory for transfer outputs. |
| `--annotation-dir` | no | `<output-dir>/annotations` | Directory for annotation intermediates. |
| `--ref-label-key` | no | `cluster` | Reference `.obs` label column. |
| `--spatial-cluster-key` | no | `leiden0.5` | Spatial guide cluster column. |
| `--neighbors` | no | `50` | K for MNN anchors and spatial KNN graph. |
| `--n-pcs` | no | `50` | PCs for the transfer space. |
| `--minkowski-order` | no | `2` | Minkowski distance order; `2` is Euclidean. |
| `--dispersion-cutoff` | no | `0.15` | Anchor pruning cutoff across spatial guide clusters. |
| `--iterations` | no | `3` | Label propagation iterations. |
| `--tie-strategy` | no | `majority` | `majority` or `weighted`. |
| `--spatial-sct-filename` | no | `sct_spatial.h5ad` | Normalized spatial file name under `spatial_dir`. |
| `--sc-sct-filename` | no | `sct_singleCell.h5ad` | Normalized reference file name under `ref_dir`. |
| `--sc-full-filename` | no | `Original_counts.h5ad` | Full/original reference file name under `ref_dir`. |
| `--output-h5ad-name` | no | `nico_celltype_annotation.h5ad` | Annotated spatial AnnData file name. |
| `--output-label-key` | no | `nico_ct` | `.obs` column for transferred labels. |
| `--overwrite / --no-overwrite` | no | no overwrite | Replace existing outputs. |
| `--cleanup-intermediate / --keep-intermediate` | no | keep | Remove temporary anchor files after success. |

### Minimal example

```bash
uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis
```

### Outputs

```text
nico_analysis/nico_celltype_annotation.h5ad
nico_analysis/annotations/anchors_data_50.npz
nico_analysis/annotations/*_nico_annotation_cluster.csv
nico_analysis/annotations/*_nico_annotation_ct_name.csv
```

The final label column is `.obs["nico_ct"]` by default.

# `nico-niche`

`nico-niche` runs and inspects spatial cell-type interaction analysis after label transfer. The main command is `run`; the other commands validate inputs or export/plot existing artifacts.

## `nico-niche run`

Runs NiCo niche interaction analysis.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing transferred labels and receiving niche outputs. |
| `--anndata-filename` | no | `nico_celltype_annotation.h5ad` | Annotated spatial file under `output_dir`. |
| `--label-key` | no | `nico_ct` | `.obs` column containing transferred labels. |
| `--spatial-key` | no | `spatial` | `.obsm` key containing coordinates. |
| `--radius` | no | `0` | `0` for Delaunay; positive value for fixed-radius neighbors. |
| `--epsilon-threshold` | no | `100.0` | Maximum Delaunay edge length when radius is `0`. |
| `--exclude-cell-type` | no | none | Additional cell type to exclude. May be repeated. |
| `--k-fold` | no | `5` | Stratified cross-validation folds. |
| `--n-repeats` | no | `1` | Repeated CV rounds after C selection. |
| `--seed` | no | `36851234` | Random seed passed to NiCo. |
| `--n-jobs` | no | `-1` | sklearn parallel jobs. |
| `--c-value` | no | NiCo grid | Candidate sklearn C value. May be repeated. |
| `--make-plots / --no-plots` | no | no plots | Generate selected plots after analysis. |
| `--plot-kind` | no | built-in defaults | Plot kind. May be repeated. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--interaction-cutoff` | no | `0.1` | Positive normalized coefficient cutoff for graph-style plots. |
| `--proximity / --no-proximity` | no | no proximity | Also run proximity analysis. |
| `--n-permutations` | no | `1000` | Proximity random permutations. |
| `--overwrite / --no-overwrite` | no | no overwrite | Replace existing niche outputs. |

### Minimal example

```bash
uv run nico-niche run \
  --output-dir nico_analysis
```

### Outputs

Important outputs include:

```text
nico_analysis/used_CT.txt
nico_analysis/used_Clusters0.csv
nico_analysis/neighbors_0.p
nico_analysis/distances_0.p
nico_analysis/niche_prediction_linear/classifier_matrices_0.npz
nico_analysis/niche_prediction_linear/niche_manifest_0.json
nico_analysis/niche_prediction_linear/metrics_0.tsv
nico_analysis/niche_prediction_linear/interactions_0.tsv
```

## `nico-niche validate`

Validates niche inputs without running NiCo.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing transferred labels. |
| `--anndata-filename` | no | `nico_celltype_annotation.h5ad` | Annotated spatial AnnData file. |
| `--label-key` | no | `nico_ct` | `.obs` label column. |
| `--spatial-key` | no | `spatial` | `.obsm` coordinate key. |
| `--radius` | no | `0` | Radius value used to plan artifacts. |
| `--epsilon-threshold` | no | `100.0` | Delaunay edge threshold. |
| `--exclude-cell-type` | no | none | Additional cell type to exclude. May be repeated. |
| `--k-fold` | no | `5` | Cross-validation folds. |
| `--overwrite / --no-overwrite` | no | no overwrite | Allow existing planned outputs. |

### Minimal example

```bash
uv run nico-niche validate --output-dir nico_analysis
```

## `nico-niche artifacts`

Prints expected niche artifacts and whether they exist.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing niche outputs. |
| `--radius` | no | `0` | Radius tag/value to inspect. |
| `--anndata-filename` | no | `nico_celltype_annotation.h5ad` | Annotated spatial AnnData file name. |
| `--check / --no-check` | no | no check | Exit nonzero when covariation artifacts are missing. |

### Minimal example

```bash
uv run nico-niche artifacts --output-dir nico_analysis
```

## `nico-niche export`

Exports a directed cell-type interaction coefficient table from existing niche artifacts.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing niche outputs. |
| `--radius` | no | `0` | Radius tag/value to export. |
| `--output` | no | default interactions path | Output TSV path. |
| `--cutoff` | no | `0.0` | Positive normalized coefficient cutoff. |
| `--include-self-edges / --exclude-self-edges` | no | include | Include source=target interactions. |

### Minimal example

```bash
uv run nico-niche export --output-dir nico_analysis
```

## `nico-niche plot`

Generates plots from existing niche artifacts.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing niche outputs. |
| `--radius` | no | `0` | Radius tag/value to plot. |
| `--kind` | no | built-in defaults | Plot kind. May be repeated. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--interaction-cutoff` | no | `0.1` | Graph cutoff. |
| `--choose-cell-type` | no | all cell types | Cell type to include in `top-coefficients` plots. May be repeated. |

### Minimal example

```bash
uv run nico-niche plot --output-dir nico_analysis
```

## `nico-niche proximity`

Runs observed-vs-randomized cell-type proximity analysis from an existing niche result.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Directory containing niche outputs. |
| `--radius` | no | `0` | Radius tag/value for proximity. |
| `--n-permutations` | no | `1000` | Number of randomized label permutations. |
| `--seed` | no | none | Optional random seed for permutations. |
| `--plot-format` | no | `pdf` | Plot file format. |

### Minimal example

```bash
uv run nico-niche proximity --output-dir nico_analysis
```

# `nico-covariation`

`nico-covariation` runs NiCo latent-factor covariation after niche analysis and provides focused export/report commands. The core analysis is `run`; report commands operate on existing covariation artifacts.

## Report-only gene filtering

The `top-genes`, `pathway`, and `reports` commands include ribosomal/mitochondrial gene symbols by default (`--include-rps-rpl-mt-genes`). Use `--exclude-rps-rpl-mt-genes` to remove candidate symbols with these exact, case-sensitive prefixes:

| `--organism` | Excluded prefixes |
|---|---|
| `mouse` (default) | `Rps`, `Rpl`, `mt-` |
| `human` | `RPS`, `RPL`, `MT-` |

Set `--organism` to match your gene symbols. NiCo filters previously computed gene–factor associations before selecting the top N genes, replacing excluded genes with the next-ranked eligible genes. This applies to single-factor, all-factor, and paired top-gene outputs, and to gene lists submitted to Enrichr for pathway enrichment.

Filtering does not change normalization, NMF/iNMF factors, gene–factor correlations, or ridge regression. Reuse existing fitted results and rerun only the affected reports. This is symbol-prefix filtering: Ensembl IDs do not match these prefixes, and it does not exclude all genes involved in mitochondrial biology.

## `nico-covariation run`

Runs core covariation analysis.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory containing niche artifacts. |
| `--ref-dir` | for double modality | none | Reference input directory. |
| `--spatial-dir` | for double modality | none | Spatial/query input directory. |
| `--radius` | no | `0` | Radius matching `nico-niche` artifacts. |
| `--n-factors` | no | `3` | Number of latent factors per cell type. |
| `--modality` | no | `double` | `double` or `single`. |
| `--factorization` | no | `inmf` | `inmf` or `nmf-transfer`. |
| `--ref-label-key` | no | `cluster` | Reference `.obs` label column. |
| `--annotated-h5ad-name` | no | `nico_celltype_annotation.h5ad` | Annotated h5ad name for single modality. |
| `--ref-original-counts-name` | no | `Original_counts.h5ad` | Original reference h5ad filename. |
| `--ref-sct-name` | no | `sct_singleCell.h5ad` | Normalized reference h5ad filename. |
| `--spatial-sct-name` | no | `sct_spatial.h5ad` | Normalized spatial h5ad filename. |
| `--ligand-receptor-db` | no | auto-detected | Ligand-receptor database path. |
| `--ridge-alpha` | no | built-in grid | RidgeCV alpha. May be repeated. |
| `--logistic-coef-cutoff` | no | `0.0` | Niche logistic coefficient cutoff. |
| `--ridge-coef-cutoff` | no | `0.0` | Ridge coefficient cutoff for exports/reports. |
| `--expression-population-cutoff` | no | `0.0` | Expression threshold for population fractions. |
| `--seed` | no | `541` | Random seed passed to NiCo. |
| `--persist-state / --no-persist-state` | no | persist | Write raw NiCo state pickle. |
| `--export-regression / --no-export-regression` | no | export | Write regression coefficient TSV. |
| `--write-manifest / --no-manifest` | no | write | Write manifest JSON. |
| `--overwrite / --no-overwrite` | no | no overwrite | Replace existing outputs. |

### Minimal example

```bash
uv run nico-covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

### Outputs

```text
nico_analysis/covariations_R0_F3/factors_info.p
nico_analysis/covariations_R0_F3/Principal_component_feature_matrix.npz
nico_analysis/covariations_R0_F3/Regression_outputs/
nico_analysis/covariations_R0_F3/covariation_state.pkl
nico_analysis/covariations_R0_F3/regression_coefficients.tsv
nico_analysis/covariations_R0_F3/covariation_manifest.json
```

## `nico-covariation validate`

Validates covariation inputs without running NiCo.

### Parameters

Same core input/planning parameters as `run`, except analysis-only export flags are omitted:

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory containing niche artifacts. |
| `--ref-dir` | for double modality | none | Reference input directory. |
| `--spatial-dir` | for double modality | none | Spatial/query input directory. |
| `--radius` | no | `0` | Radius matching niche artifacts. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--modality` | no | `double` | `double` or `single`. |
| `--factorization` | no | `inmf` | `inmf` or `nmf-transfer`. |
| `--ref-label-key` | no | `cluster` | Reference label column. |
| `--annotated-h5ad-name` | no | `nico_celltype_annotation.h5ad` | Annotated h5ad name. |
| `--ref-original-counts-name` | no | `Original_counts.h5ad` | Original reference filename. |
| `--ref-sct-name` | no | `sct_singleCell.h5ad` | Normalized reference filename. |
| `--spatial-sct-name` | no | `sct_spatial.h5ad` | Normalized spatial filename. |
| `--ligand-receptor-db` | no | auto-detected | Ligand-receptor database path. |
| `--overwrite / --no-overwrite` | no | no overwrite | Allow existing planned outputs. |

### Minimal example

```bash
uv run nico-covariation validate \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

## `nico-covariation artifacts`

Prints expected covariation artifacts and whether they exist.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value to inspect. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--check / --no-check` | no | no check | Exit nonzero when required artifacts are missing. |

### Minimal example

```bash
uv run nico-covariation artifacts --output-dir nico_analysis
```

## `nico-covariation export`

Exports stable tables from an existing result. Currently only regression export is supported.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--kind` | no | `regression` | Export kind; currently only `regression`. |
| `--output` | no | default regression path | Output TSV path. |
| `--ridge-coef-cutoff` | no | result/default cutoff | Override coefficient cutoff. |

### Minimal example

```bash
uv run nico-covariation export --output-dir nico_analysis
```

## `nico-covariation reports`

Generates optional report bundles from an existing covariation result.

Gene exclusion affects only the `top-genes-all-factors` and `pathway` report kinds, not every output in a bundle. The default bundle contains neither affected kind.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--kind` | no | report defaults | Report kind. May be repeated. |
| `--include-rps-rpl-mt-genes / --exclude-rps-rpl-mt-genes` | no | include | Retain or exclude ribosomal/mitochondrial symbols in `top-genes-all-factors` and `pathway` only. |
| `--organism` | no | `mouse` | `mouse` or `human` (lowercase); selects species-specific symbol prefixes. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--dpi` | no | `300` | Plot DPI. |
| `--show / --no-show` | no | no show | Keep figures open. |
| `--transparent / --opaque` | no | opaque | Transparent plot background. |
| `--pvalue-cutoff` | no | `0.05` | P-value cutoff. |
| `--cell-type` | no | all | Restrict to selected cell types. May be repeated. |
| `--factor-id` | no | all | Restrict to selected factor IDs. May be repeated. |

### Minimal example

```bash
uv run nico-covariation reports --output-dir nico_analysis
```

To exclude matching mouse symbols from all-factor top-gene reports:

```bash
uv run nico-covariation reports \
  --output-dir nico_analysis \
  --kind top-genes-all-factors \
  --organism mouse \
  --exclude-rps-rpl-mt-genes
```

## `nico-covariation top-genes`

Extracts or plots top genes for a cell type/factor.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--cell-type` | yes | none | Cell type to inspect. |
| `--factor-id` | for single-factor export | none | 1-based factor ID. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--top-n` | no | `30` | Number of genes to export/plot. |
| `--include-rps-rpl-mt-genes / --exclude-rps-rpl-mt-genes` | no | include | Retain or exclude ribosomal/mitochondrial symbols before top-gene selection. |
| `--organism` | no | `mouse` | `mouse` or `human` (lowercase); selects species-specific symbol prefixes. |
| `--output` | no | default path | Output TSV path. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--positive / --negative` | no | positive | Use positive or negative factor correlations. |
| `--all-factors / --single-factor` | no | single factor | Plot top genes across all factors. |
| `--pair-cell-type` | no | none | Second cell type for paired top-gene plot. |
| `--pair-factor-id` | no | none | Second factor ID for paired plot. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation top-genes \
  --output-dir nico_analysis \
  --cell-type APCs \
  --factor-id 1
```

For human symbols, using the umbrella entry point:

```bash
uv run nico-wrapper covariation top-genes \
  --output-dir nico_analysis \
  --cell-type APCs \
  --factor-id 1 \
  --organism human \
  --exclude-rps-rpl-mt-genes
```

The same species and inclusion options apply with `--all-factors` or `--pair-cell-type`/`--pair-factor-id`.

## `nico-covariation lr`

Generates ligand-receptor summaries and/or focused plots.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--central-cell-type` | no | all/summary | Central cell type for focused plots. |
| `--neighbor-cell-type` | no | all/summary | Neighbor cell type for focused plots. |
| `--central-factor-id` | no | all/summary | Central factor for focused plots. |
| `--neighbor-factor-id` | no | all/summary | Neighbor factor for focused plots. |
| `--summary / --no-summary` | no | summary | Write LR summary workbook/text. |
| `--plots / --no-plots` | no | plots | Generate LR plots. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--pvalue-cutoff` | no | `0.05` | P-value cutoff. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation lr --output-dir nico_analysis
```

## `nico-covariation pathway`

Runs optional pathway enrichment from covariation factors. This may require Enrichr/network access.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--cell-type` | no | all | Restrict to selected cell types. May be repeated. |
| `--factor-id` | no | all | Restrict to selected factors. May be repeated. |
| `--top-genes` | no | `50` | Top genes per factor for enrichment. |
| `--database` | no | GO/BioPlanet/Reactome defaults | Enrichr database. May be repeated. |
| `--organism` | no | `mouse` | `mouse` or `human` (lowercase). |
| `--include-rps-rpl-mt-genes / --exclude-rps-rpl-mt-genes` | no | include | Retain or exclude ribosomal/mitochondrial symbols before selecting genes for enrichment. |
| `--plot-as` | no | `barplot` | `barplot` or `dotplot`. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation pathway \
  --output-dir nico_analysis \
  --cell-type APCs \
  --factor-id 1
```

To exclude matching human symbols from enrichment gene lists:

```bash
uv run nico-covariation pathway \
  --output-dir nico_analysis \
  --organism human \
  --exclude-rps-rpl-mt-genes
```

## `nico-covariation umap`

Plots factor loadings on UMAP coordinates.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--cell-type` | yes | none | Cell type. Repeat once or twice. |
| `--factor-id` | yes | none | 1-based factor ID. Repeat once or twice. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--modality` | no | `spatial` | `sc` or `spatial`. |
| `--umap-key` | no | `X_umap` | AnnData `.obsm` UMAP key. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation umap \
  --output-dir nico_analysis \
  --cell-type APCs \
  --factor-id 1
```

## `nico-covariation feature-matrix`

Plots the weighted factor-neighborhood feature matrix.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation feature-matrix --output-dir nico_analysis
```

## `nico-covariation colocalize`

Plots colocalized central/neighbor factor loadings. With NiCo 1.8.0, the scatterplot uses neighbor loadings on the x-axis and central loadings on the y-axis, shows colocalized and non-colocalized central cells, and fits the regression across all central cells. Bar and violin plots compare colocalized and non-colocalized loadings for both cell types.

### Parameters

| Parameter | Required | Default | Description |
|---|---:|---|---|
| `--output-dir` | yes | none | Base NiCo output directory. |
| `--central-cell-type` | yes | none | Central cell type. |
| `--neighbor-cell-type` | yes | none | Neighbor cell type. |
| `--central-factor-id` | yes | none | Central factor ID. |
| `--neighbor-factor-id` | yes | none | Neighbor factor ID. |
| `--radius` | no | `0` | Radius tag/value. |
| `--n-factors` | no | `3` | Number of latent factors. |
| `--axis-log-scale / --no-axis-log-scale` | no | no axis log scale | Use log-scaled axes in the scatter plot. |
| `--bar / --no-bar` | no | bar | Also create bar plot. |
| `--violin / --no-violin` | no | no violin | Also create violin plot. |
| `--plot-format` | no | `pdf` | Plot file format. |
| `--show / --no-show` | no | no show | Keep figures open. |

### Minimal example

```bash
uv run nico-covariation colocalize \
  --output-dir nico_analysis \
  --central-cell-type APCs \
  --neighbor-cell-type DCs \
  --central-factor-id 1 \
  --neighbor-factor-id 1
```

## Notes on intermediate outputs

NiCo creates several intermediate files for anchors, niche interactions, neighborhoods, and covariation. The wrapper preserves the important ones for downstream steps and adds some more user-friendly sidecars. For basic usage, you normally only need to track the command inputs and the top-level output directories shown above.
