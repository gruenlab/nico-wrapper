# Starting from pre-annotated spatial data

If you already have a spatial `.h5ad` file with cell-type labels (from Seurat, cellTypist, manual annotation, etc.) you can skip `nico-preprocess build` and `nico-transfer run` and start directly from `nico-niche run`.

However, `nico-covariation run` still requires the normalized expression files produced by `nico-preprocess build`. How much of the preprocessing you need to run therefore depends on which steps of the analysis you want to perform.

---

## Which steps can be skipped?

| Goal | Minimum steps required |
|---|---|
| Niche interaction analysis only | Prepare annotated h5ad → `nico-niche run` |
| Niche + covariation analysis | Prepare annotated h5ad + run preprocess on raw counts → `nico-niche run` → `nico-covariation run` |

---

## Option A — Niche interaction only (no covariation)

### 1. Prepare the annotated spatial AnnData

Your input `.h5ad` must satisfy the following requirements for `nico-niche run`:

| Requirement | Details |
|---|---|
| Cell-type labels | A column in `.obs` containing string cell-type labels. Default column name expected: `"nico_ct"`. Pass `--label-key <your_column>` to use a different name. |
| Spatial coordinates | A key in `.obsm` containing a 2-column XY array. Default key: `"spatial"`. Pass `--spatial-key <your_key>` if different. |
| Unique cell names | `.obs_names` must be unique and non-null. |
| Minimum cell counts | At least 3 cell types with ≥ 5 cells each (after excluding the `"NM"` label, which NiCo reserves for low-quality cells). Each retained cell type also needs at least `k_fold` (default 5) cells. |

Save the file at the location `nico-niche run` expects (by default `nico_analysis/nico_celltype_annotation.h5ad`):

```python
import scanpy as sc

adata = sc.read_h5ad("your_annotated_spatial.h5ad")

# Rename your label column to the default if needed
adata.obs["nico_ct"] = adata.obs["your_celltype_column"]

# Ensure spatial coordinates are under the right key
# adata.obsm["spatial"] = adata.obsm["your_spatial_key"]

adata.write_h5ad("nico_analysis/nico_celltype_annotation.h5ad")
```

### 2. Run niche interaction analysis

```bash
uv run nico-niche run \
  --output-dir nico_analysis \
  --label-key nico_ct \
  --spatial-key spatial
```

---

## Option B — Niche + covariation analysis

Covariation requires normalized expression matrices for both modalities (produced by `nico-preprocess build`) in addition to the annotated spatial file. These files are needed by NiCo's joint iNMF factorization.

### 1. Run preprocessing on the raw counts

Even if you have your own normalized data, you still need to run `nico-preprocess build` on the **raw counts** to produce the files `nico-covariation run` expects:

```
inputRef/Original_counts.h5ad   ← raw reference counts + .raw slot
inputRef/sct_singleCell.h5ad    ← normalized reference (with .raw)
inputQuery/sct_spatial.h5ad     ← normalized spatial (with .raw)
```

> **Important:** `nico-covariation run` reads raw gene counts from the `.raw` slot of these files. The normalization needs to be done jointly on the shared gene space — this is why pre-normalized values cannot simply be reused. See [the preprocessing documentation](api.md#preprocess_nico_inputs) for details.

```bash
uv run nico-preprocess build \
  --reference inputRef_raw/reference_raw.h5ad \
  --spatial inputQuery_raw/spatial_raw.h5ad \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery \
  --no-make-reference-umap   # optional: skip reference UMAP if not needed
```

### 2. Place the annotated spatial AnnData

Put your pre-annotated file where `nico-niche run` expects it, ensuring the label and coordinate keys match (see Option A step 1 above):

```bash
cp your_annotated_spatial.h5ad nico_analysis/nico_celltype_annotation.h5ad
```

Or prepare it in Python as shown in Option A.

### 3. Run niche interaction analysis

```bash
uv run nico-niche run \
  --output-dir nico_analysis \
  --label-key nico_ct \
  --spatial-key spatial
```

### 4. Run covariation analysis

```bash
uv run nico-covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery
```

---

## Key label column names

| Step | Default `.obs` column | Override flag |
|---|---|---|
| `nico-niche run` | `nico_ct` | `--label-key` |
| `nico-covariation run` (reference) | `cluster` | `--ref-label-key` |

The reference label key (`--ref-label-key`) must match the column in `Original_counts.h5ad` that holds the original cell-type annotations — this is used by covariation to align reference cells with the spatial cell types identified by niche analysis.
