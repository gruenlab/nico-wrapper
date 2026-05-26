# Python API

This page documents the public API functions defined in the package's `pipeline.py` files. Functions are listed by their position in the overall workflow rather than alphabetically.

## Import conventions

You can import from the package-level module when the symbol is exported:

```python
from nico_wrapper.preprocess import preprocess_nico_inputs
from nico_wrapper.transfer import run_label_transfer
from nico_wrapper.niche import run_niche_interactions
from nico_wrapper.covariation import run_covariation
```

or from the concrete `pipeline.py` module:

```python
from nico_wrapper.preprocess.pipeline import preprocess_nico_inputs
```

Configuration and result dataclasses live beside the pipeline modules, for example `nico_wrapper.transfer.AnchorConfig` or `nico_wrapper.covariation.CovariationResult`.

## API stability

This repository is currently a proof of concept. The functions below are intended to be the structured public pipeline API, but signatures and behavior may still change.

# 1. Preprocessing

## `preprocess_nico_inputs`

Fully qualified import:

```python
from nico_wrapper.preprocess.pipeline import preprocess_nico_inputs
```

Builds NiCo-ready reference and spatial files from raw `.h5ad` inputs. This is the required preprocessing step before label transfer.

### Signature

```python
preprocess_nico_inputs(
    reference_h5ad: str | Path,
    spatial_h5ad: str | Path,
    *,
    ref_out_dir: str | Path,
    spatial_out_dir: str | Path,
    spatial_key: str = "spatial",
    ref_label_key: str = "cluster",
    normalization: NormalizationConfig = NiCoSCTransformConfig(),
    min_cell_counts: int = 5,
    min_gene_cells: int = 1,
    gene_space: Literal["shared", "reference_all"] = "shared",
    spatial_n_pcs: int = 30,
    leiden_resolutions: Sequence[float] = (0.4, 0.5),
    make_reference_umap: bool = True,
    random_state: int = 0,
    overwrite: bool = False,
) -> dict[str, Path]
```

### Parameters

| Parameter | Description |
|---|---|
| `reference_h5ad` | Raw reference scRNA-seq `.h5ad`. |
| `spatial_h5ad` | Raw spatial/query `.h5ad`. |
| `ref_out_dir` | Directory for reference outputs. |
| `spatial_out_dir` | Directory for spatial outputs. |
| `spatial_key` | `.obsm` key containing spatial coordinates. |
| `ref_label_key` | Reference label column in `.obs`. |
| `normalization` | `NiCoSCTransformConfig()` or `PearsonResidualsConfig()`. |
| `min_cell_counts` | Minimum total counts per retained cell. |
| `min_gene_cells` | Minimum cells per retained gene. |
| `gene_space` | `"shared"` or `"reference_all"`. |
| `spatial_n_pcs` | PCs for spatial neighbors/UMAP/Leiden. |
| `leiden_resolutions` | Spatial Leiden resolutions to compute. |
| `make_reference_umap` | Whether to compute reference UMAP in `Original_counts.h5ad`. |
| `random_state` | Random seed for Scanpy steps. |
| `overwrite` | Whether existing output files may be replaced. |

### Returns

A mapping with at least:

```text
original_counts -> <ref_out_dir>/Original_counts.h5ad
sct_single_cell -> <ref_out_dir>/sct_singleCell.h5ad
sct_spatial     -> <spatial_out_dir>/sct_spatial.h5ad
```

### Minimal example

```python
from nico_wrapper.preprocess import preprocess_nico_inputs

outputs = preprocess_nico_inputs(
    "work/reference_raw.h5ad",
    "work/spatial_raw.h5ad",
    ref_out_dir="inputRef",
    spatial_out_dir="inputQuery",
)
```

### Related CLI command

```bash
uv run nico-preprocess build --reference work/reference_raw.h5ad --spatial work/spatial_raw.h5ad --ref-out-dir inputRef --spatial-out-dir inputQuery
```

# 2. Label transfer

## `run_label_transfer`

Fully qualified import:

```python
from nico_wrapper.transfer.pipeline import run_label_transfer
```

Runs the complete NiCo label-transfer workflow: anchor discovery, label propagation, and saving the annotated spatial AnnData.

### Signature

```python
run_label_transfer(
    ref_dir: str | Path,
    spatial_dir: str | Path,
    *,
    output_dir: str | Path,
    annotation_dir: str | Path | None = None,
    config: LabelTransferConfig = LabelTransferConfig(),
) -> TransferOutputs
```

### Parameters

| Parameter | Description |
|---|---|
| `ref_dir` | Directory with `Original_counts.h5ad` and `sct_singleCell.h5ad`. |
| `spatial_dir` | Directory with `sct_spatial.h5ad`. |
| `output_dir` | Directory for transfer outputs. |
| `annotation_dir` | Optional directory for anchor/annotation intermediates. Defaults to `<output_dir>/annotations`. |
| `config` | `LabelTransferConfig` controlling anchor and annotation behavior. |

### Returns

`TransferOutputs` with paths to the annotated h5ad, anchor file if kept, and per-iteration CSVs.

### Minimal example

```python
from nico_wrapper.transfer import run_label_transfer

outputs = run_label_transfer(
    "inputRef",
    "inputQuery",
    output_dir="nico_analysis",
)
print(outputs.annotated_h5ad)
```

### Related CLI command

```bash
uv run nico-transfer run --ref-dir inputRef --spatial-dir inputQuery --output-dir nico_analysis
```

## `find_anchors`

Fully qualified import:

```python
from nico_wrapper.transfer.pipeline import find_anchors
```

Runs only NiCo's mutual-nearest-neighbor anchor discovery between the preprocessed reference and spatial data.

### Signature

```python
find_anchors(
    ref_dir: str | Path,
    spatial_dir: str | Path,
    *,
    output_dir: str | Path,
    annotation_dir: str | Path | None = None,
    config: AnchorConfig = AnchorConfig(),
    overwrite: bool = False,
) -> AnchorResult
```

### Parameters

| Parameter | Description |
|---|---|
| `ref_dir` | Directory containing reference input files named by `config`. |
| `spatial_dir` | Directory containing spatial input files named by `config`. |
| `output_dir` | Base output directory passed to NiCo. |
| `annotation_dir` | Optional intermediate output directory. |
| `config` | `AnchorConfig` with neighbors, PCs, distance order, and filenames. |
| `overwrite` | Whether existing anchor/intermediate files may be replaced. |

### Returns

`AnchorResult`, including the opaque upstream NiCo result and the anchor `.npz` path.

### Minimal example

```python
from nico_wrapper.transfer import find_anchors

anchors = find_anchors(
    "inputRef",
    "inputQuery",
    output_dir="nico_analysis",
)
```

## `transfer_labels`

Fully qualified import:

```python
from nico_wrapper.transfer.pipeline import transfer_labels
```

Runs label propagation after anchor discovery.

### Signature

```python
transfer_labels(
    anchors: AnchorResult,
    *,
    config: AnnotationConfig = AnnotationConfig(),
) -> AnnotationResult
```

### Parameters

| Parameter | Description |
|---|---|
| `anchors` | `AnchorResult` returned by `find_anchors`. |
| `config` | `AnnotationConfig` with label keys, propagation iterations, tie strategy, and output naming. |

### Returns

`AnnotationResult`, including the upstream NiCo annotation object and generated per-iteration CSV paths.

### Minimal example

```python
from nico_wrapper.transfer import find_anchors, transfer_labels

anchors = find_anchors("inputRef", "inputQuery", output_dir="nico_analysis")
annotation = transfer_labels(anchors)
```

## `save_transfer_result`

Fully qualified import:

```python
from nico_wrapper.transfer.pipeline import save_transfer_result
```

Saves the annotated spatial AnnData from a label-propagation result.

### Signature

```python
save_transfer_result(
    annotation: AnnotationResult,
    *,
    output_dir: str | Path | None = None,
    output_h5ad_name: str = "nico_celltype_annotation.h5ad",
    output_label_key: str = "nico_ct",
    overwrite: bool = False,
) -> Path
```

### Parameters

| Parameter | Description |
|---|---|
| `annotation` | `AnnotationResult` returned by `transfer_labels`. |
| `output_dir` | Optional output directory; defaults to `annotation.output_dir`. |
| `output_h5ad_name` | Output file name for the annotated AnnData. |
| `output_label_key` | `.obs` column receiving transferred labels. |
| `overwrite` | Whether to replace an existing file or label column. |

### Returns

Path to the saved annotated spatial `.h5ad`.

### Minimal example

```python
from nico_wrapper.transfer import find_anchors, save_transfer_result, transfer_labels

anchors = find_anchors("inputRef", "inputQuery", output_dir="nico_analysis")
annotation = transfer_labels(anchors)
path = save_transfer_result(annotation)
```

# 3. Niche interaction analysis

## `run_niche_interactions`

Fully qualified import:

```python
from nico_wrapper.niche.pipeline import run_niche_interactions
```

Runs complete NiCo spatial niche / cell-type interaction analysis after label transfer.

### Signature

```python
run_niche_interactions(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> NicheInteractionResult
```

### Parameters

| Parameter | Description |
|---|---|
| `output_dir` | Directory containing the annotated spatial AnnData and receiving niche artifacts. |
| `config` | `NicheInteractionConfig` with input keys, neighborhood mode, model settings, plots, proximity, and exports. |

### Returns

`NicheInteractionResult` with artifact paths, metrics, cell-type names, and an optional raw NiCo result.

### Minimal example

```python
from nico_wrapper.niche import run_niche_interactions

result = run_niche_interactions("nico_analysis")
print(result.classifier_matrices_npz)
```

### Related CLI command

```bash
uv run nico-niche run --output-dir nico_analysis
```

## `build_niche_result`

Fully qualified import:

```python
from nico_wrapper.niche.pipeline import build_niche_result
```

Builds a stable `NicheInteractionResult` from existing artifacts and optional in-memory NiCo metadata.

### Signature

```python
build_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
    nico_result: Any | None = None,
) -> NicheInteractionResult
```

### Parameters

| Parameter | Description |
|---|---|
| `output_dir` | Base NiCo output directory containing niche artifacts. |
| `radius` | Radius value/tag used for the niche run. |
| `anndata_filename` | Annotated spatial AnnData filename under `output_dir`. |
| `nico_result` | Optional raw result returned by NiCo's spatial neighborhood analysis. |

### Returns

`NicheInteractionResult` assembled from known artifact paths.

### Minimal example

```python
from nico_wrapper.niche import build_niche_result

result = build_niche_result("nico_analysis", radius=0)
```

## `load_niche_result`

Fully qualified import:

```python
from nico_wrapper.niche.pipeline import load_niche_result
```

Loads an existing niche interaction result from disk without rerunning NiCo.

### Signature

```python
load_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
) -> NicheInteractionResult
```

### Parameters

| Parameter | Description |
|---|---|
| `output_dir` | Base output directory. |
| `radius` | Radius value/tag to load. |
| `anndata_filename` | Annotated spatial AnnData filename. |

### Returns

`NicheInteractionResult` populated from the manifest when present, otherwise from standard artifact discovery.

### Minimal example

```python
from nico_wrapper.niche import load_niche_result

result = load_niche_result("nico_analysis")
```

## `export_interaction_table`

Fully qualified import:

```python
from nico_wrapper.niche.pipeline import export_interaction_table
```

Exports logistic-regression coefficients from a niche result as a directed interaction TSV.

### Signature

```python
export_interaction_table(
    result: NicheInteractionResult,
    *,
    cutoff: float = 0.0,
    normalized_cutoff: float | None = None,
    include_self_edges: bool = True,
    output_path: str | Path | None = None,
) -> Path
```

### Parameters

| Parameter | Description |
|---|---|
| `result` | Existing `NicheInteractionResult`. |
| `cutoff` | Positive normalized coefficient cutoff used when `normalized_cutoff` is not supplied. |
| `normalized_cutoff` | Explicit positive normalized coefficient cutoff. |
| `include_self_edges` | Whether to include source/target self-interactions. |
| `output_path` | Optional destination path. Defaults to the niche prediction directory. |

### Returns

Path to the exported interaction TSV.

### Minimal example

```python
from nico_wrapper.niche import export_interaction_table, load_niche_result

result = load_niche_result("nico_analysis")
path = export_interaction_table(result)
```

### Related CLI command

```bash
uv run nico-niche export --output-dir nico_analysis
```

## `summarize_niche_result`

Fully qualified import:

```python
from nico_wrapper.niche.pipeline import summarize_niche_result
```

Returns a compact dictionary summary of a niche result.

### Signature

```python
summarize_niche_result(result: NicheInteractionResult) -> dict[str, object]
```

### Parameters

| Parameter | Description |
|---|---|
| `result` | `NicheInteractionResult` to summarize. |

### Returns

Dictionary with key paths, radius, number of cell types, selected C, covariation readiness, and sidecar paths.

### Minimal example

```python
from nico_wrapper.niche import load_niche_result, summarize_niche_result

summary = summarize_niche_result(load_niche_result("nico_analysis"))
```

# 4. Covariation analysis

## `run_covariation`

Fully qualified import:

```python
from nico_wrapper.covariation.pipeline import run_covariation
```

Runs NiCo niche covariation analysis after niche interaction artifacts exist.

### Signature

```python
run_covariation(
    *,
    niche_result: NicheInteractionResult | None = None,
    output_dir: str | Path | None = None,
    ref_dir: str | Path | None = None,
    spatial_dir: str | Path | None = None,
    config: CovariationConfig = CovariationConfig(),
) -> CovariationResult
```

### Parameters

| Parameter | Description |
|---|---|
| `niche_result` | Optional existing niche result. |
| `output_dir` | Output directory containing niche artifacts; used when `niche_result` is not supplied. |
| `ref_dir` | Reference directory required for double-modality covariation. |
| `spatial_dir` | Spatial/query directory required for double-modality covariation. |
| `config` | `CovariationConfig` with radius, factors, modality, filenames, cutoffs, and export settings. |

### Returns

`CovariationResult` with covariation artifact paths and optional raw NiCo state.

### Minimal example

```python
from nico_wrapper.covariation import run_covariation

result = run_covariation(
    output_dir="nico_analysis",
    ref_dir="inputRef",
    spatial_dir="inputQuery",
)
print(result.covariation_dir)
```

### Related CLI command

```bash
uv run nico-covariation run --output-dir nico_analysis --ref-dir inputRef --spatial-dir inputQuery
```

## `build_covariation_result`

Fully qualified import:

```python
from nico_wrapper.covariation.pipeline import build_covariation_result
```

Builds a stable covariation result from existing artifacts and optional metadata.

### Signature

```python
build_covariation_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
    modality: str = "double",
    niche_result: NicheInteractionResult | None = None,
    nico_result: Any | None = None,
) -> CovariationResult
```

### Parameters

| Parameter | Description |
|---|---|
| `output_dir` | Base output directory. |
| `radius` | Radius value/tag matching niche artifacts. |
| `n_factors` | Number of factors used for the run. |
| `modality` | Covariation modality, usually `"double"` or `"single"`. |
| `niche_result` | Optional related `NicheInteractionResult`. |
| `nico_result` | Optional raw upstream NiCo result. |

### Returns

`CovariationResult` assembled from known artifact paths.

### Minimal example

```python
from nico_wrapper.covariation import build_covariation_result

result = build_covariation_result("nico_analysis", radius=0, n_factors=3)
```

## `load_covariation_result`

Fully qualified import:

```python
from nico_wrapper.covariation.pipeline import load_covariation_result
```

Loads an existing covariation result from disk.

### Signature

```python
load_covariation_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
    load_state: bool = True,
) -> CovariationResult
```

### Parameters

| Parameter | Description |
|---|---|
| `output_dir` | Base output directory. |
| `radius` | Radius value/tag. |
| `n_factors` | Number of latent factors. |
| `load_state` | Whether to load the persisted raw NiCo state pickle when available. |

### Returns

`CovariationResult` populated from artifacts, manifest, and optionally state pickle.

### Minimal example

```python
from nico_wrapper.covariation import load_covariation_result

result = load_covariation_result("nico_analysis")
```

## `export_regression_table`

Fully qualified import:

```python
from nico_wrapper.covariation.pipeline import export_regression_table
```

Exports NiCo ridge-regression coefficients to a long TSV table.

### Signature

```python
export_regression_table(
    result: CovariationResult,
    *,
    output_path: str | Path | None = None,
    ridge_coef_cutoff: float | None = None,
) -> Path
```

### Parameters

| Parameter | Description |
|---|---|
| `result` | `CovariationResult` with in-memory or persisted NiCo state. |
| `output_path` | Optional output TSV path. |
| `ridge_coef_cutoff` | Optional override for the coefficient cutoff. |

### Returns

Path to the exported regression coefficient TSV.

### Minimal example

```python
from nico_wrapper.covariation import export_regression_table, load_covariation_result

result = load_covariation_result("nico_analysis", load_state=True)
path = export_regression_table(result)
```

### Related CLI command

```bash
uv run nico-covariation export --output-dir nico_analysis
```
