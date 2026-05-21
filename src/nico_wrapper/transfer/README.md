# NiCo Transfer CLI/API

`nico_wrapper.transfer` is intended to wrap NiCo's label-transfer/annotation workflow in a small public Python API and Typer CLI.

This is the step after `nico-preprocess build`. It assumes the preprocessing output directories already contain:

```text
<ref-dir>/Original_counts.h5ad
<ref-dir>/sct_singleCell.h5ad
<spatial-dir>/sct_spatial.h5ad
```

Plotting is intentionally out of scope for this package. It can be added later as a separate visualization command/module.

## Current status

The public API, configuration objects, validation helpers, and Typer CLI are implemented. Plotting remains intentionally out of scope.

## Mechanism to be wrapped

NiCo's annotation flow has two main stages:

1. **Anchor discovery** via `nico.Annotations.find_anchor_cells_between_ref_and_query`
   - reads normalized reference and spatial AnnData objects;
   - finds shared genes;
   - computes a shared PCA transfer space;
   - finds mutual-nearest-neighbor anchors;
   - builds a spatial KNN graph;
   - writes intermediate files under the annotation directory.

2. **Label transfer / propagation** via `nico.Annotations.nico_based_annotation`
   - reads reference labels from `Original_counts.h5ad.obs[ref_label_key]`;
   - reads spatial guide clusters from `sct_spatial.h5ad.obs[spatial_cluster_key]`;
   - prunes noisy anchors using the dispersion cutoff;
   - propagates labels from anchors to neighboring spatial cells;
   - writes one pair of annotation CSVs per iteration;
   - produces labels that will be saved into an annotated spatial AnnData.

## Public API

The intended stable API is exported from `nico_wrapper.transfer`:

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

Lower-level APIs are also exported for advanced use:

```python
from nico_wrapper.transfer import find_anchors, transfer_labels, save_transfer_result
```

## CLI

Run from the repository root with:

```bash
python -m nico_wrapper.transfer.cli --help
```

Main command:

```bash
python -m nico_wrapper.transfer.cli run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis
```

A more explicit example:

```bash
python -m nico_wrapper.transfer.cli run \
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

## CLI options

Core inputs:

```text
--ref-dir PATH             Directory with Original_counts.h5ad and sct_singleCell.h5ad
--spatial-dir PATH         Directory with sct_spatial.h5ad
--output-dir PATH          Directory for NiCo transfer outputs
--annotation-dir PATH      Optional annotation/intermediate directory
```

Label-transfer behavior:

```text
--ref-label-key TEXT       Reference .obs label column; default: cluster
--spatial-cluster-key TEXT Spatial guide cluster column; default: leiden0.5
--neighbors INTEGER        K for anchors and spatial KNN graph; default: 50; must be >= 2
--n-pcs INTEGER            PCs for transfer space; default: 50
--minkowski-order INTEGER  Distance order; default: 2
--dispersion-cutoff FLOAT  Anchor pruning cutoff; default: 0.15
--iterations INTEGER       Propagation iterations; default: 3
--tie-strategy TEXT        majority or weighted; default: majority
```

File-name escape hatches:

```text
--spatial-sct-filename TEXT  default: sct_spatial.h5ad
--sc-sct-filename TEXT       default: sct_singleCell.h5ad
--sc-full-filename TEXT      default: Original_counts.h5ad
--output-h5ad-name TEXT      default: nico_celltype_annotation.h5ad
--output-label-key TEXT      default: nico_ct
```

Operational flags:

```text
--overwrite / --no-overwrite
--cleanup-intermediate / --keep-intermediate
```

Defaults follow NiCo where practical. In particular, the default spatial guide cluster key is `leiden0.5`, and intermediates are kept by default.

## Outputs

The complete workflow returns and/or prints paths for:

```text
<output-dir>/<output-h5ad-name>
<annotation-dir>/anchors_data_<neighbors>.npz
<annotation-dir>/1_nico_annotation_cluster.csv
<annotation-dir>/1_nico_annotation_ct_name.csv
...
```

The annotated AnnData should contain transferred labels in:

```python
adata.obs[output_label_key]
```

with `output_label_key="nico_ct"` by default.
