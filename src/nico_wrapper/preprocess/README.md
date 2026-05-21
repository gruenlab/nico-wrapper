# NiCo Preprocessing CLI

This directory contains a small, standalone preprocessing tool for preparing reference scRNA-seq and spatial/Xenium query data for the downstream NiCo workflow.

The tool has two optional conversion commands and one required joint build command:

```text
convert-reference-sparse   # optional: sparse reference files -> raw reference h5ad
convert-spatial-csv        # optional: spatial count/coordinate CSV files -> raw spatial h5ad
build                      # required: reference h5ad + spatial h5ad -> NiCo-ready outputs
```

Run commands from the repository root with:

```bash
python -m nico_wrapper.preprocess.cli --help
```

## Overall workflow

If both reference and spatial query data are already available as `.h5ad` files, skip the conversion commands and run only `build`.

If one or both inputs are CSV/sparse exports, first convert them to raw `.h5ad`, then run `build`. The conversion commands create parent directories for their output paths when needed.

```bash
# Optional: only needed for sparse reference input
python -m nico_wrapper.preprocess.cli convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad

# Optional: only needed for spatial CSV input
python -m nico_wrapper.preprocess.cli convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad

# Required joint preprocessing step
python -m nico_wrapper.preprocess.cli build \
  --reference work/reference_raw.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery
```

## Command: `convert-reference-sparse`

Converts a reference scRNA-seq dataset stored as sparse triplet files into a raw `.h5ad` file.

Example:

```bash
python -m nico_wrapper.preprocess.cli convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad
```

### Assumed input format

`--counts` is expected to be a sparse triplet table with three columns:

```text
gene_index  cell_index  count
```

By default, indices are assumed to be one-based, matching MatrixMarket-style exports. Use `--zero-based-indices` if the file already uses Python-style zero-based indexing.

`--genes` should contain gene identifiers/names.

`--barcodes` should contain cell barcodes.

The output `.h5ad` has:

- cells in `.obs_names`
- genes in `.var_names`
- raw counts in `.X` as a sparse float matrix

The float matrix dtype is intentional. NiCo's SCTransform implementation modifies sparse matrix data in place during log/geometric-mean calculations; integer sparse matrices can silently truncate intermediate values and lead to downstream `NaN`/`inf` failures.

### Main options

```text
--counts PATH              sparse triplet count table
--genes PATH               gene metadata CSV
--barcodes PATH            cell barcode CSV
--output PATH              output raw reference h5ad
--sep TEXT                 separator for the sparse count table; default: space
--header INTEGER           header row passed to pandas; default: 1
--one-based-indices / --zero-based-indices
```

## Command: `convert-spatial-csv`

Converts spatial/Xenium count and coordinate CSV files into a raw spatial `.h5ad` file.

Example for the tutorial-style `gene_by_cell.csv`:

```bash
python -m nico_wrapper.preprocess.cli convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad \
  --counts-orientation genes_by_cells
```

### Assumed input format

`--counts` is a dense count matrix CSV.

By default, it is assumed to be oriented as:

```text
rows    = genes
columns = cells
```

This corresponds to `--counts-orientation genes_by_cells`. The loader transposes the matrix so that the resulting `AnnData` has cells as observations and genes as variables.

If the CSV already has cells as rows and genes as columns, use:

```bash
--counts-orientation cells_by_genes
```

`--coordinates` is expected to contain one barcode column and one or more coordinate columns. Coordinates are stored in:

```python
adata.obsm["spatial"]
```

by default.

The command can reorder coordinate rows to match the count matrix cell order using barcodes. This is enabled by default. Barcode matching is strict: every count-matrix cell must have coordinates, and extra coordinate barcodes are rejected.

The output `.h5ad` stores counts as a sparse float matrix. This is intentional for compatibility with the default NiCo SCTransform normalization path.

### Main options

```text
--counts PATH                         spatial count CSV
--coordinates PATH                    coordinate CSV
--output PATH                         output raw spatial h5ad
--counts-orientation TEXT             genes_by_cells or cells_by_genes
--barcode-col TEXT                    coordinate barcode column name or index; default: 0
--coordinate-col TEXT                 coordinate column name/index; may be repeated
--spatial-key TEXT                    obsm key for coordinates; default: spatial
--reorder-coordinates / --no-reorder-coordinates
```

## Command: `build`

Runs the joint NiCo preprocessing step on a raw reference `.h5ad` and a raw spatial `.h5ad`.

Example:

```bash
python -m nico_wrapper.preprocess.cli build \
  --reference inputRef/input_ref.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery \
  --normalization nico-sctransform \
  --leiden-resolution 0.4 \
  --leiden-resolution 0.5
```

### What `build` does

The build step is where the two modalities are processed together. It is responsible for:

1. reading the reference and spatial `.h5ad` files
2. validating required fields
3. filtering low-count cells and low-support genes on the raw full gene space
4. saving full filtered reference counts as `Original_counts.h5ad`
5. finding shared genes between reference and spatial data
6. aligning modalities to the selected gene space
7. filtering again after gene-space alignment
8. normalizing both modalities with the selected strategy
9. computing PCA/neighbors/UMAP on the normalized spatial data
10. computing Leiden clusters for the requested resolutions
11. writing NiCo-ready output files

### Required inputs

Reference `.h5ad`:

- cells in `.obs_names`
- genes in `.var_names`
- raw counts in `.X`, unless a normalization strategy explicitly uses a layer
- reference cell type/cluster column in `.obs`, default: `cluster`

Spatial `.h5ad`:

- cells in `.obs_names`
- genes in `.var_names`
- raw counts in `.X`, unless a normalization strategy explicitly uses a layer
- spatial coordinates in `.obsm["spatial"]` by default

The spatial coordinate matrix must have one row per spatial cell. Coordinates may have two or more dimensions; the tutorial input includes `Xcoord`, `Ycoord`, and `Zcoord`, so `.obsm["spatial"]` is `n_cells x 3`.

### Outputs

The build command writes the files expected by downstream NiCo functions:

```text
<ref-out-dir>/Original_counts.h5ad
<ref-out-dir>/sct_singleCell.h5ad
<spatial-out-dir>/sct_spatial.h5ad
```

`Original_counts.h5ad` contains the reference count data and a `.raw` copy. If `--make-reference-umap` is enabled, it also contains a standard Scanpy normalized/log-transformed PCA/neighbors/UMAP representation for reference inspection.

`sct_singleCell.h5ad` contains the normalized reference data.

`sct_spatial.h5ad` contains the normalized spatial data, spatial coordinates, UMAP, neighbors, and Leiden cluster columns such as:

```text
leiden0.4
leiden0.5
```

These Leiden columns can be passed later to NiCo as `guiding_spatial_cluster_resolution_tag`, for example `leiden0.4`.

## Normalization strategies

The build command supports two normalization strategies:

```text
--normalization nico-sctransform
--normalization pearson-residuals
```

The selected normalization is applied to both the reference and spatial data after filtering and gene-space alignment. The build step filters once before alignment and again after alignment. The second pass is important because cells can pass whole-transcriptome filtering but have too few or zero counts in the shared gene space, which can make SCTransform or Pearson residual normalization unstable.

### `nico-sctransform`

This is the default and mirrors the original NiCo notebook/paper path using NiCo's SCTransform-like implementation.

Note: SCTransform is sensitive to sparse integer matrices because its implementation performs in-place transformations on sparse `.data`. The converters therefore store count matrices as sparse float matrices. If you provide your own `.h5ad`, make sure raw counts are non-negative and preferably stored with a float-compatible dtype.

```bash
--normalization nico-sctransform
```

Method-specific options:

```text
--sct-min-cells INTEGER      SCTransform-only: min cells per gene; default: 1
--sct-gmean-eps FLOAT        SCTransform-only: geometric mean epsilon; default: 1.0
--sct-n-genes INTEGER        SCTransform-only: genes sampled for fitting; default: 500
--sct-n-cells INTEGER        SCTransform-only: cells sampled for fitting; default: all cells
--sct-bin-size INTEGER       SCTransform-only: genes per fitting bin; default: 500
--sct-bw-adjust FLOAT        SCTransform-only: bandwidth adjustment; default: 3.0
```

These options are ignored when `--normalization pearson-residuals` is selected.

### `pearson-residuals`

This uses Scanpy's analytic Pearson residual normalization:

```python
scanpy.experimental.pp.normalize_pearson_residuals
```

Example:

```bash
python -m nico_wrapper.preprocess.cli build \
  --reference inputRef/input_ref.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery \
  --normalization pearson-residuals \
  --pearson-theta 100
```

Method-specific options:

```text
--pearson-theta FLOAT                         Pearson-only: overdispersion theta; default: 100
--pearson-clip FLOAT                          Pearson-only: clipping threshold; default: sqrt(n_obs)
--pearson-check-values / --no-pearson-check-values
--pearson-layer TEXT                          Pearson-only: layer to normalize instead of .X
```

These options are ignored when `--normalization nico-sctransform` is selected.

## Shared gene handling

By default, build uses:

```bash
--gene-space shared
```

This subsets both reference and spatial data to genes shared by the two modalities before normalization. This is the safest default and follows the intent of the notebooks. After subsetting, the build step repeats count/gene filtering to remove cells or genes that become uninformative in the shared gene space.

An alternative mode is available:

```bash
--gene-space reference_all
```

This preserves all reference genes while still aligning spatial data to shared genes. This mirrors one variant in the original notebook, but `shared` is recommended unless there is a specific reason to preserve the full reference gene set.

## Typical tutorial-style usage

For the original tutorial where the reference is already `input_ref.h5ad` and the spatial query is `gene_by_cell.csv` plus `tissue_positions_list.csv`:

```bash
mkdir -p work

python -m nico_wrapper.preprocess.cli convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad \
  --counts-orientation genes_by_cells

python -m nico_wrapper.preprocess.cli build \
  --reference inputRef/input_ref.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery \
  --normalization nico-sctransform \
  --leiden-resolution 0.4 \
  --leiden-resolution 0.5
```

For sparse-format reference input:

```bash
mkdir -p work

python -m nico_wrapper.preprocess.cli convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad

python -m nico_wrapper.preprocess.cli convert-spatial-csv \
  --counts inputQuery/counts.csv \
  --coordinates inputQuery/tissue_pos.csv \
  --output work/spatial_raw.h5ad \
  --counts-orientation cells_by_genes

python -m nico_wrapper.preprocess.cli build \
  --reference work/reference_raw.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery
```

## Current implementation status

The conversion commands and build pipeline are implemented. The tutorial-style input archives `inputRef.zip` and `inputQuery.zip` have been smoke-tested with both `pearson-residuals` and the default `nico-sctransform` normalization path.

Two implementation details from testing are worth keeping in mind:

1. Converted count matrices are written as sparse float matrices, not sparse integer matrices, to avoid SCTransform truncation/`NaN` issues.
2. The build step filters both before and after shared-gene alignment, because shared-gene subsetting can create zero/low-count cells even when the original full matrix passed filtering.
