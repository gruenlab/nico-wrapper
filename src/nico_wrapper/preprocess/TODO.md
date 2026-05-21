# TODO

## Add raw Xenium-to-AnnData assembly command

Add a `nico-preprocess` CLI command for assembling raw 10x Xenium output directories into a raw spatial/query `.h5ad` that can be passed to `nico-preprocess build`.

Proposed command shape:

```bash
python -m nico_wrapper.preprocess.cli convert-xenium \
  --xenium-dir /path/to/xenium_query \
  --output work/spatial_raw.h5ad \
  --spatial-key spatial
```

Implementation instructions:

1. Read the Xenium feature-count matrix from:

   ```text
   <xenium-dir>/cell_feature_matrix.h5
   ```

   using `scanpy.read_10x_h5(...)`.

2. Ensure gene names are unique with `adata.var_names_make_unique()`.

3. Convert `.X` to a float-compatible matrix, preferably sparse `float32`, because downstream NiCo SCTransform modifies sparse matrix data in place and integer sparse matrices can cause truncation/`NaN` issues.

4. Read cell metadata from:

   ```text
   <xenium-dir>/cells.csv.gz
   ```

   Required columns:

   ```text
   cell_id
   x_centroid
   y_centroid
   ```

5. Match `cells.csv.gz[cell_id]` exactly to `adata.obs_names` from the feature matrix.

   - Error if any matrix cell is missing from the cells table.
   - Extra cells in the cells table may be ignored with a warning.
   - Reorder the cells table to match `adata.obs_names` before joining metadata.

6. Join all useful cell metadata columns into `adata.obs`, excluding the duplicate `cell_id` column after it has been used as the index.

7. Store spatial coordinates in:

   ```python
   adata.obsm[spatial_key] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
   ```

8. Store lightweight provenance, for example:

   ```python
   adata.uns["source_xenium_dir"] = str(xenium_dir)
   adata.uns["spatial_coordinate_columns"] = ["x_centroid", "y_centroid"]
   ```

9. Validate the assembled object with `validate_spatial_adata(adata, spatial_key=spatial_key)` before writing.

10. Write the final raw spatial AnnData to `--output` and create parent directories as needed.

The resulting file should be usable as:

```bash
python -m nico_wrapper.preprocess.cli build \
  --reference /path/to/reference.h5ad \
  --spatial work/spatial_raw.h5ad \
  --ref-out-dir output/preprocess_ref \
  --spatial-out-dir output/preprocess_query
```
