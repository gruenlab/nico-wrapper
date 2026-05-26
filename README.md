This is a test commit
# nico-wrapper

> [!WARNING]
> This is a vibe-coded proof of concept for a NiCo wrapper. It can change at any moment, may contain bugs, and should not be treated as a stable production interface.

`nico-wrapper` provides small Python APIs and Typer CLIs around the NiCo (`nico-sc-sp`) preprocessing and label-transfer workflow for reference scRNA-seq and spatial/Xenium query data.

The package exposes four tools:

- `nico-preprocess`: convert raw inputs and build NiCo-ready AnnData files.
- `nico-transfer`: transfer reference labels onto preprocessed spatial/query data.
- `nico-niche`: run NiCo niche / spatial cell-type interaction analysis after transfer.
- `nico-covariation`: run NiCo latent-factor covariation analysis after niche interactions.

Most plotting is optional and disabled by default.

## Installation

This project is managed with [`uv`](https://docs.astral.sh/uv/) and requires Python 3.11.

```bash
uv venv --python 3.11
uv sync
```

For platform-specific system dependency notes, especially for `pygraphviz`, see [`docs/setup.md`](docs/setup.md).

Run the CLIs from the repository root with:

```bash
uv run nico-preprocess --help
uv run nico-transfer --help
uv run nico-niche --help
uv run nico-covariation --help
```

Equivalent module entry points are also available:

```bash
uv run python -m nico_wrapper.preprocess.cli --help
uv run python -m nico_wrapper.transfer.cli --help
uv run python -m nico_wrapper.niche.cli --help
uv run python -m nico_wrapper.covariation.cli --help
```

## End-to-end workflow

If both reference and spatial query inputs are already `.h5ad` files, skip the conversion commands and run only `nico-preprocess build` followed by `nico-transfer run`.

```bash
# Optional: sparse reference files -> raw reference h5ad
uv run nico-preprocess convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad

# Optional: spatial count/coordinate CSV files -> raw spatial h5ad
uv run nico-preprocess convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad \
  --counts-orientation genes_by_cells

# Required: raw h5ad inputs -> NiCo-ready preprocessing outputs
uv run nico-preprocess build \
  --reference work/reference_raw.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery

# Label transfer from reference to spatial/query cells
uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis

# Niche interaction analysis on transferred labels
uv run nico-niche run \
  --output-dir nico_analysis

# Factor covariation analysis after niche interactions
uv run nico-covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

## `nico-preprocess`

`nico-preprocess` prepares reference and spatial/query data for downstream NiCo annotation.

It provides three commands:

```text
convert-reference-sparse   optional: sparse reference files -> raw reference h5ad
convert-spatial-csv        optional: spatial count/coordinate CSV files -> raw spatial h5ad
build                      required: reference h5ad + spatial h5ad -> NiCo-ready outputs
```

### `convert-reference-sparse`

Converts a reference scRNA-seq dataset stored as sparse triplet files into a raw `.h5ad` file.

```bash
uv run nico-preprocess convert-reference-sparse \
  --counts inputRef/counts.csv \
  --genes inputRef/genes.csv \
  --barcodes inputRef/barcodes.csv \
  --output work/reference_raw.h5ad
```

Expected count format:

```text
gene_index  cell_index  count
```

Main options:

```text
--counts PATH              sparse triplet count table
--genes PATH               gene metadata CSV
--barcodes PATH            cell barcode CSV
--output PATH              output raw reference h5ad
--sep TEXT                 separator for the sparse count table; default: space
--header INTEGER           header row passed to pandas; default: 1
--one-based-indices / --zero-based-indices
                           default: one-based indices
```

The output `.h5ad` stores cells in `.obs_names`, genes in `.var_names`, and raw counts in `.X` as a sparse float matrix. The float dtype avoids in-place truncation issues in NiCo's SCTransform implementation.

### `convert-spatial-csv`

Converts spatial/Xenium count and coordinate CSV files into a raw spatial `.h5ad` file.

```bash
uv run nico-preprocess convert-spatial-csv \
  --counts inputQuery/gene_by_cell.csv \
  --coordinates inputQuery/tissue_positions_list.csv \
  --output work/spatial_raw.h5ad \
  --counts-orientation genes_by_cells
```

By default, the count matrix is interpreted as genes by cells and transposed so cells become observations. Use `--counts-orientation cells_by_genes` if cells are already rows.

Coordinates are stored in `.obsm["spatial"]` by default. Barcode matching is strict: every count-matrix cell must have coordinates, and extra coordinate barcodes are rejected.

Main options:

```text
--counts PATH                         spatial count CSV
--coordinates PATH                    coordinate CSV
--output PATH                         output raw spatial h5ad
--counts-orientation TEXT             genes_by_cells or cells_by_genes; default: genes_by_cells
--barcode-col TEXT                    coordinate barcode column name or index; default: 0
--coordinate-col TEXT                 coordinate column name/index; may be repeated; default: all non-barcode columns
--spatial-key TEXT                    obsm key for coordinates; default: spatial
--reorder-coordinates / --no-reorder-coordinates
                                      default: reorder coordinates
```

### `build`

Runs the joint preprocessing step on raw reference and spatial `.h5ad` files.

```bash
uv run nico-preprocess build \
  --reference inputRef/input_ref.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery \
  --normalization nico-sctransform \
  --leiden-resolution 0.4 \
  --leiden-resolution 0.5
```

The build step:

1. reads and validates reference/spatial `.h5ad` files;
2. filters low-count cells and low-support genes;
3. writes filtered reference counts as `Original_counts.h5ad`;
4. aligns the modalities to the selected gene space;
5. normalizes both modalities;
6. computes spatial PCA/neighbors/UMAP and Leiden clusters;
7. writes NiCo-ready output files.

Required reference input:

- cells in `.obs_names`
- genes in `.var_names`
- raw counts in `.X`, unless a normalization strategy uses a layer
- reference label column in `.obs`, default: `cluster`

Required spatial input:

- cells in `.obs_names`
- genes in `.var_names`
- raw counts in `.X`, unless a normalization strategy uses a layer
- coordinates in `.obsm["spatial"]` by default

Main build options and defaults:

```text
--spatial-key TEXT                    default: spatial
--ref-label-key TEXT                  default: cluster
--min-cell-counts INTEGER             default: 5
--min-gene-cells INTEGER              default: 1
--gene-space TEXT                     shared or reference_all; default: shared
--spatial-n-pcs INTEGER               default: 30
--leiden-resolution FLOAT             may be repeated; default: 0.4 and 0.5
--make-reference-umap / --no-make-reference-umap
                                      default: make reference UMAP
--random-state INTEGER                default: 0
--overwrite / --no-overwrite          default: no overwrite
```

Outputs:

```text
<ref-out-dir>/Original_counts.h5ad
<ref-out-dir>/sct_singleCell.h5ad
<spatial-out-dir>/sct_spatial.h5ad
```

`sct_spatial.h5ad` contains normalized spatial data, coordinates, UMAP, neighbors, and Leiden columns such as `leiden0.4` and `leiden0.5`. These can be used later as the spatial guide cluster key for transfer.

### Preprocessing normalization

Supported strategies:

```text
--normalization nico-sctransform     default; uses NiCo's SCTransform-like implementation
--normalization pearson-residuals    uses scanpy.experimental.pp.normalize_pearson_residuals
```

NiCo SCTransform-specific options:

```text
--sct-min-cells INTEGER      default: 1
--sct-gmean-eps FLOAT        default: 1.0
--sct-n-genes INTEGER        default: 500
--sct-n-cells INTEGER        default: all cells
--sct-bin-size INTEGER       default: 500
--sct-bw-adjust FLOAT        default: 3.0
```

Pearson residual-specific options:

```text
--pearson-theta FLOAT       overdispersion theta; default: 100
--pearson-clip FLOAT        clipping threshold; default: Scanpy default sqrt(n_obs)
--pearson-check-values / --no-pearson-check-values
                            validate count values; default: check values
--pearson-layer TEXT        layer to normalize instead of .X; default: .X
```

By default, build uses `--gene-space shared`, which subsets both modalities to shared genes before normalization. `--gene-space reference_all` preserves all reference genes while still aligning spatial data to shared genes.

## `nico-transfer`

`nico-transfer` wraps NiCo's label-transfer/annotation workflow after preprocessing.

It expects these files:

```text
<ref-dir>/Original_counts.h5ad
<ref-dir>/sct_singleCell.h5ad
<spatial-dir>/sct_spatial.h5ad
```

Run transfer with:

```bash
uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis
```

A more explicit example:

```bash
uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis \
  --ref-label-key cluster \
  --spatial-cluster-key leiden0.5 \
  --neighbors 50 \
  --n-pcs 50 \
  --minkowski-order 2 \
  --dispersion-cutoff 0.15 \
  --iterations 3 \
  --tie-strategy majority \
  --output-label-key nico_ct \
  --output-h5ad-name nico_celltype_annotation.h5ad \
  --keep-intermediate
```

### Transfer mechanism

NiCo transfer has two stages:

1. Anchor discovery with `nico.Annotations.find_anchor_cells_between_ref_and_query`.
2. Label propagation with `nico.Annotations.nico_based_annotation`.

The wrapper validates inputs and outputs, preserves access to lower-level APIs, and saves transferred labels into an annotated spatial AnnData object.

### Transfer options

Core inputs:

```text
--ref-dir PATH             Directory with Original_counts.h5ad and sct_singleCell.h5ad
--spatial-dir PATH         Directory with sct_spatial.h5ad
--output-dir PATH          Directory for NiCo transfer outputs
--annotation-dir PATH      Optional annotation/intermediate directory; default: output-dir/annotations
```

Label-transfer behavior:

```text
--ref-label-key TEXT       Reference .obs label column; default: cluster
--spatial-cluster-key TEXT Spatial guide cluster column; default: leiden0.5
--neighbors INTEGER        K for anchors and spatial KNN graph; default: 50
--n-pcs INTEGER            PCs for transfer space; default: 50
--minkowski-order INTEGER  Distance order; default: 2
--dispersion-cutoff FLOAT  Anchor pruning cutoff; default: 0.15
--iterations INTEGER       Propagation iterations; default: 3
--tie-strategy TEXT        majority or weighted; default: majority
```

File-name overrides:

```text
--spatial-sct-filename TEXT  default: sct_spatial.h5ad
--sc-sct-filename TEXT       default: sct_singleCell.h5ad
--sc-full-filename TEXT      default: Original_counts.h5ad
--output-h5ad-name TEXT      default: nico_celltype_annotation.h5ad
--output-label-key TEXT      default: nico_ct
```

Operational flags:

```text
--overwrite / --no-overwrite                  default: no overwrite
--cleanup-intermediate / --keep-intermediate  default: keep intermediates
```

Defaults follow NiCo where practical. In particular, the default spatial guide cluster key is `leiden0.5`, and intermediates are kept by default.

### Transfer outputs

The complete transfer workflow returns and/or prints paths for:

```text
<output-dir>/<output-h5ad-name>
<annotation-dir>/anchors_data_<neighbors>.npz
<annotation-dir>/1_nico_annotation_cluster.csv
<annotation-dir>/1_nico_annotation_ct_name.csv
...
```

The annotated AnnData contains transferred labels in:

```python
adata.obs[output_label_key]
```

with `output_label_key="nico_ct"` by default.

### Transfer Python API

The stable API is exported from `nico_wrapper.transfer`:

```python
from nico_wrapper.transfer import (
    AnchorConfig,
    AnnotationConfig,
    LabelTransferConfig,
    TieStrategy,
    run_label_transfer,
    validate_label_transfer_inputs,
)

config = LabelTransferConfig(
    annotation=AnnotationConfig(
        ref_label_key="cluster",
        spatial_cluster_key="leiden0.5",
        tie_strategy=TieStrategy.MAJORITY,
    ),
)

validate_label_transfer_inputs(
    ref_dir="inputRef",
    spatial_dir="inputQuery",
    output_dir="nico_analysis",
    config=config,
)

outputs = run_label_transfer(
    ref_dir="inputRef",
    spatial_dir="inputQuery",
    output_dir="nico_analysis",
    config=config,
)
```

Lower-level APIs are also exported:

```python
from nico_wrapper.transfer import find_anchors, transfer_labels, save_transfer_result
```

## `nico-niche`

`nico-niche` wraps NiCo's spatial neighborhood / niche interaction analysis after label transfer.

It expects an annotated spatial AnnData under the transfer output directory:

```text
<output-dir>/nico_celltype_annotation.h5ad
```

with transferred labels in `.obs["nico_ct"]` and spatial coordinates in `.obsm["spatial"]` by default.

Run the default Delaunay-neighborhood analysis with:

```bash
uv run nico-niche run \
  --output-dir nico_analysis
```

A more explicit example:

```bash
uv run nico-niche run \
  --output-dir nico_analysis \
  --anndata-filename nico_celltype_annotation.h5ad \
  --label-key nico_ct \
  --spatial-key spatial \
  --radius 0 \
  --epsilon-threshold 100 \
  --k-fold 5 \
  --n-repeats 1 \
  --seed 36851234 \
  --n-jobs -1
```

`--radius 0` uses NiCo's Delaunay-neighborhood mode. Positive radius values use fixed-radius neighborhoods.

Useful supporting commands:

```bash
uv run nico-niche validate --output-dir nico_analysis
uv run nico-niche artifacts --output-dir nico_analysis --radius 0
uv run nico-niche export --output-dir nico_analysis --radius 0
```

Core outputs for downstream covariation are:

```text
<output-dir>/used_CT.txt
<output-dir>/used_Clusters0.csv
<output-dir>/neighbors_0.p
<output-dir>/distances_0.p
<output-dir>/niche_prediction_linear/classifier_matrices_0.npz
```

The wrapper also writes easier-to-consume sidecars when available:

```text
<output-dir>/niche_prediction_linear/niche_manifest_0.json
<output-dir>/niche_prediction_linear/metrics_0.tsv
<output-dir>/niche_prediction_linear/interactions_0.tsv
```

### Niche Python API

```python
from nico_wrapper.niche import run_niche_interactions

result = run_niche_interactions(output_dir="nico_analysis")
print(result.classifier_matrices_npz)
print(result.covariation_ready)
```

For custom settings:

```python
from nico_wrapper.niche import (
    InteractionModelConfig,
    NeighborhoodConfig,
    NicheInteractionConfig,
    run_niche_interactions,
)

config = NicheInteractionConfig(
    neighborhood=NeighborhoodConfig(radius=0, epsilon_threshold=100),
    model=InteractionModelConfig(k_fold=5, n_repeats=1),
)

result = run_niche_interactions("nico_analysis", config=config)
```

## `nico-covariation`

`nico-covariation` wraps NiCo's latent-factor covariation analysis after `nico-niche` has produced neighborhood artifacts.

Default double-modality run:

```bash
uv run nico-covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

For custom reference label keys, pass the same reference label key used during preprocessing/transfer, for example `--ref-label-key CellType`.

Useful commands:

```bash
uv run nico-covariation validate --output-dir nico_analysis --ref-dir inputRef --spatial-dir inputQuery
uv run nico-covariation artifacts --output-dir nico_analysis --radius 0 --n-factors 3
uv run nico-covariation export --output-dir nico_analysis --kind regression
uv run nico-covariation reports --output-dir nico_analysis --kind regression-circleplots --kind regression-heatmaps
uv run nico-covariation top-genes --output-dir nico_analysis --cell-type APCs --factor-id 1
uv run nico-covariation lr --output-dir nico_analysis --central-cell-type APCs --neighbor-cell-type DCs --central-factor-id 1 --neighbor-factor-id 1
uv run nico-covariation umap --output-dir nico_analysis --modality spatial --cell-type APCs --factor-id 1
uv run nico-covariation pathway --output-dir nico_analysis --cell-type APCs --factor-id 1  # may require Enrichr/network access
```

Core outputs:

```text
<output-dir>/covariations_R0_F3/factors_info.p
<output-dir>/covariations_R0_F3/Principal_component_feature_matrix.npz
<output-dir>/covariations_R0_F3/Regression_outputs/
<output-dir>/covariations_R0_F3/covariation_state.pkl
<output-dir>/covariations_R0_F3/regression_coefficients.tsv
<output-dir>/covariations_R0_F3/covariation_manifest.json
```

Additional focused report commands are available for ligand-receptor plots, top genes, pathway enrichment, UMAP factor plots, feature-matrix plots, and colocalization. Report generation writes `report_manifest.json` when using the `reports` dispatcher. The state pickle can be large because NiCo keeps report inputs in memory. Disable it with `--no-persist-state` if you only need core artifacts and the regression TSV. NiCo's ridge p-values are an upstream OLS-style approximation and are not corrected by the wrapper.

## Implementation notes

- Converted count matrices are stored as sparse float matrices to avoid SCTransform truncation and `NaN`/`inf` issues.
- The build step filters both before and after shared-gene alignment because shared-gene subsetting can create zero/low-count cells.
- The default output file names match NiCo's downstream annotation expectations.
