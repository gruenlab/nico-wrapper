# Plan: Public API and CLI for NiCo niche interaction analysis

This document proposes a step-by-step implementation plan for a public wrapper API and CLI around NiCo's niche / spatial cell-type interaction module.

The goal is to make niche interaction analysis easy to run after label transfer and to generate all artifacts required by downstream functionality, especially a future wrapper around NiCo 1 covariation analysis.

Primary reference: [`NICHE.md`](NICHE.md).

---

## 1. Context and workflow position

The intended high-level workflow is:

```text
1. nico-preprocess
   Raw reference/spatial data -> NiCo-ready h5ad files

2. nico-transfer
   Preprocessed reference/spatial data -> annotated spatial AnnData

3. nico-niche
   Annotated spatial AnnData -> niche interaction artifacts

4. nico-covariation, future work
   Niche interaction artifacts + expression data -> covariation analysis
```

The niche wrapper should assume that label transfer has already produced an annotated spatial AnnData file, by default:

```text
<nico-output-dir>/nico_celltype_annotation.h5ad
```

with transferred labels in:

```python
adata.obs["nico_ct"]
```

and spatial coordinates in:

```python
adata.obsm["spatial"]
```

The wrapper's core job is to call NiCo's upstream interaction analysis safely, validate inputs and outputs, normalize paths, prevent accidental artifact corruption, and expose a clean Python API and CLI.

---

## 2. Design goals

### 2.1 User-facing goals

1. Provide a one-command CLI for the common use case:

   ```bash
   nico-niche run --output-dir nico_analysis
   ```

2. Provide a simple Python API:

   ```python
   from nico_wrapper.niche import run_niche_interactions

   result = run_niche_interactions(output_dir="nico_analysis")
   ```

3. Preserve compatibility with NiCo covariation by writing NiCo's expected artifact names and locations.

4. Give clear validation errors before running expensive upstream NiCo code.

5. Make outputs discoverable through a typed result object and a manifest file.

6. Keep the default path aligned with current NiCo covariation expectations: linear niche prediction only.

### 2.2 Engineering goals

1. Wrap upstream NiCo behavior without exposing NiCo's raw `SimpleNamespace` as the stable public API.
2. Avoid surprising path behavior from upstream string concatenation.
3. Avoid accidental overwrite or inconsistent artifacts.
4. Keep configuration explicit and typed.
5. Make future covariation wrapper integration straightforward.
6. Separate analysis, plotting, artifact inspection, and proximity functionality.

---

## 3. Proposed package structure

Add a new package:

```text
src/nico_wrapper/niche/
├── __init__.py
├── cli.py
├── config.py
├── pipeline.py
├── plotting.py
├── proximity.py
├── results.py
├── validation.py
└── io.py
```

Suggested responsibilities:

| File | Responsibility |
|---|---|
| `config.py` | Dataclasses/enums for user configuration. |
| `results.py` | Stable result dataclasses and metric/interaction table schemas. |
| `validation.py` | Pure validation helpers for inputs, config, and output paths. |
| `pipeline.py` | Main API functions that call upstream NiCo and write wrapper artifacts. |
| `plotting.py` | Optional public plotting wrappers around NiCo plotting functions. |
| `proximity.py` | Optional proximity analysis helpers/wrappers. |
| `io.py` | Manifest loading/writing, artifact discovery, table export. |
| `cli.py` | Typer CLI app. |
| `__init__.py` | Public API exports. |

Add a console entry point to `pyproject.toml`:

```toml
[project.scripts]
nico-preprocess = "nico_wrapper.preprocess.cli:app"
nico-transfer = "nico_wrapper.transfer.cli:app"
nico-niche = "nico_wrapper.niche.cli:app"
```

---

## 4. Core artifact contract

The wrapper must preserve these NiCo-compatible files because future covariation analysis expects them:

```text
<output-dir>/used_CT.txt
<output-dir>/used_Clusters{Radius}.csv
<output-dir>/neighbors_{Radius}.p
<output-dir>/distances_{Radius}.p
<output-dir>/niche_prediction_linear/classifier_matrices_{Radius}.npz
```

The upstream niche analysis also writes:

```text
<output-dir>/niche_prediction_linear/normalized_spatial_neighborhood_{Radius}.npz
```

The wrapper should add sidecar artifacts that are easier for users and future wrappers to consume:

```text
<output-dir>/niche_prediction_linear/niche_manifest_{Radius}.json
<output-dir>/niche_prediction_linear/metrics_{Radius}.tsv
<output-dir>/niche_prediction_linear/interactions_{Radius}.tsv
```

Optional later sidecars:

```text
<output-dir>/niche_prediction_linear/proximity_pairs_{Radius}.tsv
<output-dir>/niche_prediction_linear/proximity_ratios_{Radius}.tsv
```

---

## 5. Important NiCo behaviors to wrap carefully

The upstream function is:

```python
nico.Interactions.spatial_neighborhood_analysis(...)
```

Important behaviors from `NICHE.md`:

1. `output_nico_dir` is string-concatenated with filenames. The wrapper must pass a string with a trailing path separator.
2. `Radius == 0` means Delaunay neighborhoods, not zero-radius search.
3. Nonzero `Radius` means fixed-radius neighbors.
4. Delaunay edges are filtered by `epsilonThreshold` using strict `<`.
5. The output subdirectory is always `niche_prediction_linear/` on the public path.
6. `NM` is always excluded internally.
7. Additional exclusions are exact cell-type string matches.
8. Cell types with fewer than 5 cells are silently dropped upstream.
9. Zero-neighbor cells produce NaN features and are removed before modeling.
10. The selected logistic-regression `C` is returned in memory but not saved by upstream.
11. Metrics are returned in memory but not saved by upstream.
12. The fitted model is not saved.
13. ROC and predicted probabilities correspond only to the final CV split.
14. `used_CT.txt` is not radius-specific and can become inconsistent if multiple runs with different filters share one output directory.
15. Covariation currently expects `niche_prediction_linear/`, so cross-term support should not be exposed as a default stable path.

---

## 6. Public Python API design

### 6.1 Main imports

The intended public API should be:

```python
from nico_wrapper.niche import (
    NeighborhoodConfig,
    InteractionModelConfig,
    NichePlotConfig,
    ProximityConfig,
    NicheInteractionConfig,
    NicheInteractionResult,
    ValidationError,
    validate_niche_inputs,
    run_niche_interactions,
    load_niche_result,
    validate_covariation_artifacts,
    export_interaction_table,
)
```

### 6.2 Main one-shot function

```python
def run_niche_interactions(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> NicheInteractionResult:
    """Run complete NiCo niche interaction analysis on transferred labels."""
```

Typical usage:

```python
from nico_wrapper.niche import run_niche_interactions

result = run_niche_interactions(output_dir="nico_analysis")
print(result.classifier_matrices_npz)
print(result.covariation_ready)
```

Explicit usage:

```python
from nico_wrapper.niche import (
    InteractionModelConfig,
    NeighborhoodConfig,
    NicheInteractionConfig,
    NichePlotConfig,
    run_niche_interactions,
)

config = NicheInteractionConfig(
    anndata_filename="nico_celltype_annotation.h5ad",
    label_key="nico_ct",
    spatial_key="spatial",
    neighborhood=NeighborhoodConfig(
        radius=0,
        epsilon_threshold=100,
        additional_excluded_cell_types=("LowQuality",),
    ),
    model=InteractionModelConfig(
        k_fold=5,
        n_repeats=1,
        seed=36851234,
        n_jobs=-1,
    ),
    plots=NichePlotConfig(enabled=True),
    overwrite=False,
)

result = run_niche_interactions("nico_analysis", config=config)
```

### 6.3 Lower-level API functions

```python
def validate_niche_inputs(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> None:
    """Validate annotated AnnData, config values, and planned outputs."""
```

```python
def load_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
) -> NicheInteractionResult:
    """Load an existing niche result from manifest or discovered artifacts."""
```

```python
def validate_covariation_artifacts(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
) -> bool:
    """Check that all artifacts required by NiCo covariation exist."""
```

```python
def export_interaction_table(
    result: NicheInteractionResult,
    *,
    cutoff: float = 0.0,
    normalized_cutoff: float | None = None,
    include_self_edges: bool = True,
    output_path: str | Path | None = None,
) -> Path:
    """Export logistic-regression coefficients as a directed interaction table."""
```

```python
def summarize_niche_result(result: NicheInteractionResult) -> dict[str, object]:
    """Return a lightweight summary useful for CLI printing and notebooks."""
```

---

## 7. Configuration dataclasses

### 7.1 `NeighborhoodConfig`

```python
@dataclass(frozen=True)
class NeighborhoodConfig:
    """Neighborhood construction configuration."""

    radius: int | float = 0
    epsilon_threshold: float = 100.0
    additional_excluded_cell_types: tuple[str, ...] = ()
```

Field details:

| Field | Meaning |
|---|---|
| `radius` | `0` means Delaunay; positive values mean fixed-radius search. |
| `epsilon_threshold` | Maximum Delaunay edge length. Used only when `radius == 0`. |
| `additional_excluded_cell_types` | Extra cell-type labels to exclude in addition to upstream NiCo's implicit `NM`. |

Implementation note: keep the public field name lowercase `radius`, but pass it to upstream as `Radius`.

### 7.2 `InteractionModelConfig`

```python
@dataclass(frozen=True)
class InteractionModelConfig:
    """Logistic-regression model and CV configuration."""

    k_fold: int = 5
    n_repeats: int = 1
    seed: int = 36851234
    n_jobs: int = -1
    c_values: tuple[float, ...] | None = None
```

Field details:

| Field | Meaning |
|---|---|
| `k_fold` | Stratified CV folds. Each retained class needs at least this many modeled cells. |
| `n_repeats` | Number of repeated CV rounds after hyperparameter selection. |
| `seed` | Seed passed to upstream NiCo. |
| `n_jobs` | sklearn parallelism. |
| `c_values` | Candidate sklearn inverse regularization strengths. If `None`, use NiCo default `2**np.arange(-12, 12)`. |

### 7.3 `NichePlotConfig`

```python
@dataclass(frozen=True)
class NichePlotConfig:
    """Optional plotting configuration."""

    enabled: bool = False
    kinds: tuple[str, ...] = ("confusion", "coefficients", "scores", "graph")
    saveas: str = "pdf"
    dpi: int = 300
    transparent: bool = False
    show: bool = False
    interaction_cutoff: float = 0.1
    graph_edge_labels: bool = False
```

Supported `kinds`:

```text
confusion
coefficients
scores
roc
predicted-probabilities
top-coefficients
graph
all
```

Initial implementation can support only a safe subset and leave the rest for a later milestone.

### 7.4 `ProximityConfig`

```python
@dataclass(frozen=True)
class ProximityConfig:
    """Observed-vs-randomized cell-type proximity configuration."""

    enabled: bool = False
    n_permutations: int = 1000
    observed_threshold: float = 0.0
    remove_self_pairs: bool = True
    as_counts: bool = True
    seed: int | None = None
    saveas: str = "pdf"
```

Implementation note: upstream proximity randomization is unseeded. If `seed` is provided, the wrapper should set NumPy's random seed around the upstream call and restore prior state if practical.

### 7.5 `NicheInteractionConfig`

```python
@dataclass(frozen=True)
class NicheInteractionConfig:
    """Top-level niche interaction analysis config."""

    anndata_filename: str = "nico_celltype_annotation.h5ad"
    label_key: str = "nico_ct"
    spatial_key: str = "spatial"
    neighborhood: NeighborhoodConfig = field(default_factory=NeighborhoodConfig)
    model: InteractionModelConfig = field(default_factory=InteractionModelConfig)
    plots: NichePlotConfig = field(default_factory=NichePlotConfig)
    proximity: ProximityConfig = field(default_factory=ProximityConfig)
    overwrite: bool = False
    write_manifest: bool = True
    export_interactions: bool = True
```

---

## 8. Result dataclasses

### 8.1 `NicheInteractionResult`

```python
@dataclass(frozen=True)
class NicheInteractionResult:
    """Stable wrapper result for one niche interaction run."""

    output_dir: Path
    annotated_h5ad: Path
    prediction_dir: Path

    radius: int | float
    radius_tag: str

    used_cell_types_tsv: Path
    used_clusters_csv: Path
    neighbors_pickle: Path
    distances_pickle: Path
    normalized_neighborhood_npz: Path
    classifier_matrices_npz: Path

    manifest_json: Path | None = None
    metrics_tsv: Path | None = None
    interactions_tsv: Path | None = None

    selected_c: float | None = None
    metrics: dict[str, tuple[float, float]] | None = None
    cell_type_names: dict[int, str] | None = None
    classes: tuple[int, ...] | None = None
    covariation_ready: bool = False

    nico_result: Any | None = None
```

Notes:

1. `nico_result` should preserve access to the upstream `SimpleNamespace` when available from a fresh run.
2. Loaded results from disk may not have `nico_result`, `selected_c`, or in-memory ROC/probability arrays unless wrapper sidecars were written.
3. `covariation_ready` should be true only if all required downstream covariation artifacts exist.

### 8.2 Interaction table schema

Write to:

```text
<niche-prediction-dir>/interactions_{Radius}.tsv
```

Columns:

```text
source_cell_type
target_cell_type
source_cell_type_id
target_cell_type_id
coefficient
coefficient_std
normalized_coefficient
abs_normalized_coefficient
is_self_interaction
passes_positive_cutoff
radius
```

Interpretation:

- `source_cell_type`: neighboring cell type / feature column.
- `target_cell_type`: central predicted cell type / coefficient row.
- Directed interaction is `source_cell_type -> target_cell_type`.
- `normalized_coefficient` should be coefficient divided by the global maximum absolute coefficient, matching NiCo graph interpretation.

### 8.3 Metrics table schema

Write to:

```text
<niche-prediction-dir>/metrics_{Radius}.tsv
```

Columns:

```text
metric
mean
std
radius
```

Metric order from upstream score array:

| Index | Metric |
|---:|---|
| 0 | accuracy |
| 1 | macro_f1 |
| 2 | macro_precision |
| 3 | macro_recall |
| 4 | micro_f1 |
| 5 | micro_precision |
| 6 | micro_recall |
| 7 | weighted_f1 |
| 8 | weighted_precision |
| 9 | weighted_recall |
| 10 | cohen_kappa |
| 11 | log_loss |
| 12 | matthews_corrcoef |
| 13 | hamming_loss |
| 14 | zero_one_loss |

### 8.4 Manifest schema

Write to:

```text
<niche-prediction-dir>/niche_manifest_{Radius}.json
```

Suggested keys:

```json
{
  "schema_version": 1,
  "created_at": "2026-05-21T00:00:00Z",
  "wrapper_version": "0.1.0",
  "nico_package": "nico-sc-sp==1.6.0",
  "output_dir": "nico_analysis",
  "annotated_h5ad": "nico_analysis/nico_celltype_annotation.h5ad",
  "label_key": "nico_ct",
  "spatial_key": "spatial",
  "radius": 0,
  "radius_tag": "0",
  "epsilon_threshold": 100,
  "excluded_cell_types": ["NM"],
  "model": {
    "k_fold": 5,
    "n_repeats": 1,
    "seed": 36851234,
    "n_jobs": -1,
    "c_values": null,
    "selected_c": 1.0
  },
  "artifacts": {
    "used_cell_types_tsv": "nico_analysis/used_CT.txt",
    "used_clusters_csv": "nico_analysis/used_Clusters0.csv",
    "neighbors_pickle": "nico_analysis/neighbors_0.p",
    "distances_pickle": "nico_analysis/distances_0.p",
    "normalized_neighborhood_npz": "nico_analysis/niche_prediction_linear/normalized_spatial_neighborhood_0.npz",
    "classifier_matrices_npz": "nico_analysis/niche_prediction_linear/classifier_matrices_0.npz",
    "metrics_tsv": "nico_analysis/niche_prediction_linear/metrics_0.tsv",
    "interactions_tsv": "nico_analysis/niche_prediction_linear/interactions_0.tsv"
  },
  "covariation_ready": true
}
```

---

## 9. CLI design

Add a new Typer app:

```text
nico-niche
```

### 9.1 CLI command overview

```text
nico-niche run          Run niche interaction analysis end to end.
nico-niche validate     Validate inputs/config without running analysis.
nico-niche artifacts    Print/check produced artifacts for a given radius.
nico-niche plot         Generate plots from an existing result.
nico-niche export       Export interaction coefficient table from existing result.
nico-niche proximity    Run observed-vs-randomized proximity analysis.
```

Initial MVP can implement only:

```text
run
validate
artifacts
```

Then add `plot`, `export`, and `proximity` in later phases.

### 9.2 `nico-niche run`

Minimal usage:

```bash
uv run nico-niche run \
  --output-dir nico_analysis
```

Explicit usage:

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
  --n-jobs -1 \
  --make-plots \
  --plot-format pdf \
  --interaction-cutoff 0.1 \
  --overwrite
```

Suggested options:

| Option | Default | Meaning |
|---|---:|---|
| `--output-dir PATH` | required | Directory containing transferred annotated h5ad and receiving outputs. |
| `--anndata-filename TEXT` | `nico_celltype_annotation.h5ad` | Annotated spatial AnnData filename. |
| `--label-key TEXT` | `nico_ct` | `.obs` label key. |
| `--spatial-key TEXT` | `spatial` | `.obsm` coordinate key. |
| `--radius TEXT` | `0` | `0` means Delaunay; positive value means fixed-radius search. |
| `--epsilon-threshold FLOAT` | `100` | Delaunay edge length threshold. |
| `--exclude-cell-type TEXT` | repeatable | Additional cell types to exclude. |
| `--k-fold INTEGER` | `5` | CV folds. |
| `--n-repeats INTEGER` | `1` | Repeated CV rounds. |
| `--seed INTEGER` | `36851234` | Random seed. |
| `--n-jobs INTEGER` | `-1` | sklearn parallelism. |
| `--c-value FLOAT` | repeatable | Candidate sklearn C values. If omitted, NiCo defaults. |
| `--make-plots / --no-plots` | `--no-plots` | Generate selected plots after analysis. |
| `--plot-kind TEXT` | repeatable | Plot kinds to generate. |
| `--plot-format TEXT` | `pdf` | Plot file format. |
| `--interaction-cutoff FLOAT` | `0.1` | Positive normalized coefficient cutoff for graphs/tables. |
| `--proximity / --no-proximity` | `--no-proximity` | Also run proximity analysis. |
| `--n-permutations INTEGER` | `1000` | Proximity random permutations. |
| `--overwrite / --no-overwrite` | `--no-overwrite` | Replace existing outputs. |

Important CLI parsing note:

- Parse `--radius` as text, then convert it to `int` when it represents an integer.
- This prevents accidental filenames like `neighbors_0.0.p` when downstream tools expect `neighbors_0.p`.

Example parser behavior:

| CLI input | Internal value | Radius tag |
|---|---:|---|
| `--radius 0` | `0` | `0` |
| `--radius 50` | `50` | `50` |
| `--radius 50.5` | `50.5` | `50.5` |

### 9.3 `nico-niche validate`

Usage:

```bash
uv run nico-niche validate \
  --output-dir nico_analysis \
  --label-key nico_ct \
  --spatial-key spatial \
  --radius 0
```

Responsibilities:

1. Load and validate annotated AnnData.
2. Validate label and coordinate keys.
3. Validate model/neighborhood scalar config.
4. Validate output paths and overwrite policy.
5. Print a clear success message if valid.

Output example:

```text
validated: nico_analysis/nico_celltype_annotation.h5ad
cells: 12345
cell_types: 17
radius: 0
planned_prediction_dir: nico_analysis/niche_prediction_linear
```

### 9.4 `nico-niche artifacts`

Usage:

```bash
uv run nico-niche artifacts \
  --output-dir nico_analysis \
  --radius 0
```

Responsibilities:

1. Discover expected files for a radius.
2. Print paths and existence status.
3. Exit nonzero if required artifacts are missing when `--check` is used.

Output example:

```text
used_cell_types_tsv: nico_analysis/used_CT.txt [ok]
used_clusters_csv: nico_analysis/used_Clusters0.csv [ok]
neighbors_pickle: nico_analysis/neighbors_0.p [ok]
distances_pickle: nico_analysis/distances_0.p [ok]
classifier_matrices_npz: nico_analysis/niche_prediction_linear/classifier_matrices_0.npz [ok]
normalized_neighborhood_npz: nico_analysis/niche_prediction_linear/normalized_spatial_neighborhood_0.npz [ok]
covariation_ready: true
```

### 9.5 `nico-niche export`

Usage:

```bash
uv run nico-niche export \
  --output-dir nico_analysis \
  --radius 0 \
  --cutoff 0.1 \
  --output nico_analysis/niche_prediction_linear/interactions_0.tsv
```

Responsibilities:

1. Load `classifier_matrices_{Radius}.npz`.
2. Load cell-type names from `used_CT.txt`.
3. Convert coefficient matrix into a directed edge table.
4. Write TSV.

### 9.6 `nico-niche plot`

Usage:

```bash
uv run nico-niche plot \
  --output-dir nico_analysis \
  --radius 0 \
  --kind all \
  --plot-format pdf
```

Supported kinds:

```text
confusion
coefficients
scores
roc
predicted-probabilities
top-coefficients
graph
all
```

Implementation note:

- Some upstream plot functions need the in-memory `SimpleNamespace` returned by the original analysis, especially ROC and predicted-probability plots.
- For loaded results, the wrapper may only be able to support plots based on persisted artifacts unless the manifest includes enough fields.
- MVP should prioritize persisted-artifact plots: confusion matrix, coefficient matrix, interaction graph, top coefficients if possible.

### 9.7 `nico-niche proximity`

Usage:

```bash
uv run nico-niche proximity \
  --output-dir nico_analysis \
  --radius 0 \
  --n-permutations 1000 \
  --seed 42 \
  --plot-format pdf
```

Responsibilities:

1. Load existing neighbor/profile artifacts.
2. Run observed-vs-randomized proximity analysis.
3. Save plot and optional TSV sidecars.

Implementation note:

- This should be separate from main logistic-regression niche interaction modeling because it is a different analysis, even though it uses the same neighbor graph.

---

## 10. Validation plan

Create `src/nico_wrapper/niche/validation.py`.

### 10.1 Config validation

Implement:

```python
def validate_niche_config(config: NicheInteractionConfig) -> None:
    ...
```

Checks:

1. `anndata_filename` is a plain filename, not an absolute path or nested path.
2. `label_key` is a non-empty string.
3. `spatial_key` is a non-empty string.
4. `radius >= 0`.
5. `epsilon_threshold > 0`.
6. `k_fold >= 2`.
7. `n_repeats >= 1`.
8. `seed` is an integer.
9. `c_values`, if provided, are positive finite numbers.
10. `additional_excluded_cell_types` contains non-empty strings.
11. Plot format is among supported values, e.g. `pdf`, `png`, `svg`.
12. Proximity permutation count is non-negative or positive depending on command semantics.

### 10.2 Input AnnData validation

Implement:

```python
def validate_annotated_spatial_adata(
    adata: AnnData,
    *,
    label_key: str = "nico_ct",
    spatial_key: str = "spatial",
    k_fold: int = 5,
    additional_excluded_cell_types: Sequence[str] = (),
) -> None:
    ...
```

Checks:

1. `adata.n_obs > 0`.
2. `adata.obs_names` are unique and non-missing.
3. `label_key in adata.obs`.
4. `adata.obs[label_key]` has no missing values.
5. Labels are non-empty after casting to string.
6. Unique labels can be sorted by NiCo.
7. `spatial_key in adata.obsm`.
8. Coordinates are a 2D numeric matrix.
9. Coordinate rows equal `adata.n_obs`.
10. Coordinate columns are either 2 or 3 for upstream NiCo compatibility.
11. Coordinates are finite.
12. After excluding `NM` and additional excluded labels, at least two cell types remain.
13. After applying the known upstream minimum-count filter of 5 cells per cell type, at least two cell types remain.
14. Every retained class should have at least `k_fold` cells before modeling.

Note:

- Upstream later drops zero-neighbor cells from modeling. The wrapper cannot fully know this without building neighborhoods, but it can catch obvious class-size problems early.

### 10.3 Output path validation

Implement:

```python
def planned_niche_artifacts(
    output_dir: str | Path,
    *,
    radius: int | float | str,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
) -> NicheArtifactPaths:
    ...
```

Checks:

1. `output_dir` exists or can be created.
2. `output_dir` is a directory if it exists.
3. `niche_prediction_linear` path is either absent or a directory.
4. Planned output files do not already exist unless `overwrite=True`.
5. Be especially careful with `used_CT.txt`, because it is not radius-specific.

Planned core files:

```text
<output-dir>/used_CT.txt
<output-dir>/used_Clusters{radius_tag}.csv
<output-dir>/neighbors_{radius_tag}.p
<output-dir>/distances_{radius_tag}.p
<output-dir>/niche_prediction_linear/normalized_spatial_neighborhood_{radius_tag}.npz
<output-dir>/niche_prediction_linear/classifier_matrices_{radius_tag}.npz
```

Planned wrapper sidecars:

```text
<output-dir>/niche_prediction_linear/niche_manifest_{radius_tag}.json
<output-dir>/niche_prediction_linear/metrics_{radius_tag}.tsv
<output-dir>/niche_prediction_linear/interactions_{radius_tag}.tsv
```

### 10.4 Covariation artifact validation

Implement:

```python
def validate_covariation_artifacts(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
) -> bool:
    ...
```

Required files:

```text
used_CT.txt
used_Clusters{Radius}.csv
neighbors_{Radius}.p
distances_{Radius}.p
niche_prediction_linear/classifier_matrices_{Radius}.npz
```

Additional checks:

1. Files exist and are non-empty.
2. `classifier_matrices_{Radius}.npz` contains keys:

   ```text
   cmn
   coef
   cmn_std
   coef_std
   CTFeatures
   ```

3. `used_CT.txt` can be parsed into integer ID and name columns.
4. Coefficient matrix shape agrees with number of used cell types where possible.

---

## 11. Pipeline implementation plan

Create `src/nico_wrapper/niche/pipeline.py`.

### 11.1 Main function skeleton

```python
def run_niche_interactions(
    output_dir: str | Path,
    *,
    config: NicheInteractionConfig = NicheInteractionConfig(),
) -> NicheInteractionResult:
    validate_niche_inputs(output_dir, config=config)

    resolved_output_dir = Path(output_dir)
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    nico_result = _run_upstream_spatial_neighborhood_analysis(
        output_dir=resolved_output_dir,
        config=config,
    )

    result = _build_result_from_upstream(
        output_dir=resolved_output_dir,
        config=config,
        nico_result=nico_result,
    )

    if config.export_interactions:
        interactions_tsv = export_interaction_table(result)
        result = replace(result, interactions_tsv=interactions_tsv)

    metrics_tsv = export_metrics_table(result, nico_result=nico_result)
    result = replace(result, metrics_tsv=metrics_tsv)

    if config.write_manifest:
        manifest_json = write_niche_manifest(result, config=config)
        result = replace(result, manifest_json=manifest_json)

    if config.plots.enabled:
        plot_niche_result(result, config=config.plots)

    if config.proximity.enabled:
        run_proximity_analysis(result, config=config.proximity)

    return result
```

### 11.2 Upstream call wrapper

Implement:

```python
def _run_upstream_spatial_neighborhood_analysis(
    *,
    output_dir: Path,
    config: NicheInteractionConfig,
) -> Any:
    try:
        from nico import Interactions as sint
    except ImportError as exc:
        raise ImportError("NiCo niche analysis requires the 'nico' package.") from exc

    model = config.model
    neighborhood = config.neighborhood

    c_values = (
        list(model.c_values)
        if model.c_values is not None
        else list(np.power(2.0, np.arange(-12, 12)))
    )

    return sint.spatial_neighborhood_analysis(
        output_nico_dir=_as_nico_dir(output_dir),
        anndata_object_name=config.anndata_filename,
        spatial_cluster_tag=config.label_key,
        spatial_coordinate_tag=config.spatial_key,
        Radius=neighborhood.radius,
        n_repeats=model.n_repeats,
        K_fold=model.k_fold,
        seed=model.seed,
        n_jobs=model.n_jobs,
        lambda_c_ranges=c_values,
        epsilonThreshold=neighborhood.epsilon_threshold,
        removed_CTs_before_finding_CT_CT_interactions=list(
            neighborhood.additional_excluded_cell_types
        ),
    )
```

Use helper:

```python
def _as_nico_dir(path: str | Path) -> str:
    text = os.fspath(Path(path))
    return text if text.endswith(os.sep) else text + os.sep
```

### 11.3 Result construction

Implement:

```python
def build_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    anndata_filename: str = "nico_celltype_annotation.h5ad",
    nico_result: Any | None = None,
) -> NicheInteractionResult:
    ...
```

Responsibilities:

1. Resolve `radius_tag`.
2. Construct expected paths.
3. Read selected C, classes, metrics, and cell-type names from `nico_result` if available.
4. Load cell-type names from `used_CT.txt` if `nico_result` is absent.
5. Validate covariation artifacts.
6. Return a `NicheInteractionResult`.

### 11.4 Metrics export

Implement:

```python
def export_metrics_table(
    result: NicheInteractionResult,
    *,
    nico_result: Any | None = None,
    output_path: str | Path | None = None,
) -> Path | None:
    ...
```

Behavior:

- If `nico_result` has `score`, export all 15 metrics.
- If no scores are available, return `None` or raise depending on caller mode.

### 11.5 Interaction export

Implement:

```python
def export_interaction_table(
    result: NicheInteractionResult,
    *,
    cutoff: float = 0.0,
    normalized_cutoff: float | None = None,
    include_self_edges: bool = True,
    output_path: str | Path | None = None,
) -> Path:
    ...
```

Steps:

1. Load `classifier_matrices_npz`.
2. Extract `coef` and `coef_std`.
3. Load cell-type names from result or `used_CT.txt`.
4. Compute global maximum absolute coefficient.
5. For each target class row and source feature column, emit one row.
6. Normalize coefficients by global max absolute coefficient.
7. Mark whether the edge passes positive cutoff.
8. Optionally skip self-edges.
9. Write TSV.

Important:

- In the public linear path, feature columns map directly to cell-type IDs.
- If future cross-term features appear, the exporter should detect `CTFeatures` containing spaces and either raise a clear unsupported error or encode them separately.

### 11.6 Manifest write/load

Implement:

```python
def write_niche_manifest(
    result: NicheInteractionResult,
    *,
    config: NicheInteractionConfig,
    output_path: str | Path | None = None,
) -> Path:
    ...
```

```python
def load_niche_result(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
) -> NicheInteractionResult:
    ...
```

Loading behavior:

1. Prefer manifest if present.
2. Fall back to artifact discovery.
3. Parse cell types from `used_CT.txt`.
4. Parse metrics/interactions sidecars if present.
5. Do not require `nico_result`.

---

## 12. Plotting implementation plan

Create `src/nico_wrapper/niche/plotting.py`.

### 12.1 Initial plotting support

Initial safe wrappers:

```python
def plot_confusion_matrix(result: NicheInteractionResult, *, config: NichePlotConfig) -> Path:
    ...
```

```python
def plot_coefficient_matrix(result: NicheInteractionResult, *, config: NichePlotConfig) -> Path:
    ...
```

```python
def plot_interaction_graph(result: NicheInteractionResult, *, config: NichePlotConfig) -> Path:
    ...
```

```python
def plot_niche_result(result: NicheInteractionResult, *, config: NichePlotConfig) -> list[Path]:
    ...
```

### 12.2 Handling upstream plot functions

Some upstream plot functions expect the original `SimpleNamespace` returned by `spatial_neighborhood_analysis()`. The wrapper can support two modes:

1. If `result.nico_result is not None`, call upstream plot functions directly.
2. If loading from disk, reconstruct a minimal namespace where possible.

Minimal namespace fields needed by common plots:

```text
fout
niche_pred_outdir
Radius
classes
nameOfCellType
lambda_c
score
```

Possible helper:

```python
def to_nico_plot_namespace(result: NicheInteractionResult) -> SimpleNamespace:
    ...
```

Limitations:

- `plot_predicted_probabilities` needs final split arrays not saved by upstream.
- `plot_roc_results` needs ROC arrays not saved by upstream.
- These should be supported only immediately after a fresh run unless the wrapper later persists those arrays.

### 12.3 Plot support priority

Phase 1:

1. Confusion matrix.
2. Coefficient matrix.
3. Evaluation scores if metrics are available.
4. Interaction graph from persisted coefficient matrix.

Phase 2:

1. Top coefficients by central cell type.
2. Proximity plots.
3. ROC and predicted probabilities only for fresh runs or if persisted.

---

## 13. Proximity implementation plan

Create `src/nico_wrapper/niche/proximity.py`.

### 13.1 Public API

```python
def compute_celltype_proximity_pairs(
    result: NicheInteractionResult,
    *,
    config: ProximityConfig = ProximityConfig(),
) -> ProximityResult:
    ...
```

```python
def plot_celltype_proximity_pairs(
    result: NicheInteractionResult,
    *,
    config: ProximityConfig = ProximityConfig(),
) -> list[Path]:
    ...
```

```python
def run_proximity_analysis(
    result: NicheInteractionResult,
    *,
    config: ProximityConfig = ProximityConfig(),
) -> ProximityResult:
    ...
```

### 13.2 Result object

```python
@dataclass(frozen=True)
class ProximityResult:
    observed_tsv: Path | None
    ratio_tsv: Path | None
    plot_path: Path | None
    observed: dict[str, float]
    ratio: dict[str, float]
```

### 13.3 Seed handling

Upstream randomization uses `np.random.permutation` without a seed. If `config.seed` is provided:

1. Save NumPy global RNG state.
2. Set `np.random.seed(config.seed)`.
3. Call upstream function.
4. Restore previous RNG state in `finally`.

This makes CLI proximity runs reproducible without permanently mutating global RNG state.

---

## 14. Step-by-step implementation phases

### Phase 0: Confirm scope and decisions

Before implementation, decide:

1. Should `nico-niche run` make plots by default?
   - Recommendation: default `--no-plots` for speed and reproducibility; users can enable `--make-plots`.
2. Should proximity run by default?
   - Recommendation: no; keep it separate because it is distinct from logistic-regression interactions and can be slow.
3. Should the wrapper allow multiple runs with different filters in the same `output_dir`?
   - Recommendation: discourage this; protect with overwrite checks and manifest warnings because `used_CT.txt` is not radius-specific.
4. Should cross-term niche prediction be exposed?
   - Recommendation: no for now; downstream covariation expects linear artifacts.

### Phase 1: Add package skeleton and config/results

Files:

```text
src/nico_wrapper/niche/__init__.py
src/nico_wrapper/niche/config.py
src/nico_wrapper/niche/results.py
```

Tasks:

1. Create dataclasses:
   - `NeighborhoodConfig`
   - `InteractionModelConfig`
   - `NichePlotConfig`
   - `ProximityConfig`
   - `NicheInteractionConfig`
   - `NicheInteractionResult`
   - optional `NicheArtifactPaths`
2. Add `__all__` exports in `__init__.py`.
3. Keep defaults aligned with upstream NiCo.

Acceptance criteria:

- `from nico_wrapper.niche import NicheInteractionConfig` works.
- Dataclasses instantiate with defaults.

### Phase 2: Implement validation

File:

```text
src/nico_wrapper/niche/validation.py
```

Tasks:

1. Add `ValidationError` or reuse/export the existing transfer validation error style.
2. Implement `require_file` or reuse an existing common helper if one is later extracted.
3. Implement config validation.
4. Implement annotated AnnData validation.
5. Implement output artifact planning and overwrite validation.
6. Implement covariation artifact validation.

Acceptance criteria:

- Missing annotated h5ad raises a clear error.
- Missing label key raises a clear error.
- Missing spatial key raises a clear error.
- Invalid coordinates raise a clear error.
- Existing outputs are rejected unless `overwrite=True`.

### Phase 3: Implement core pipeline

File:

```text
src/nico_wrapper/niche/pipeline.py
```

Tasks:

1. Implement `_as_nico_dir`.
2. Implement radius parsing/tag helper.
3. Implement upstream call wrapper.
4. Implement `run_niche_interactions`.
5. Implement result construction from paths and upstream namespace.
6. Implement checks that upstream created all expected artifacts.

Acceptance criteria:

- `run_niche_interactions(output_dir="nico_analysis")` calls upstream NiCo with normalized paths.
- Result contains all expected paths.
- Missing upstream-created files raise clear runtime errors.
- `covariation_ready` is true only when required files exist.

### Phase 4: Implement metrics and interaction exports

Files:

```text
src/nico_wrapper/niche/io.py
src/nico_wrapper/niche/pipeline.py
```

Tasks:

1. Implement `load_used_cell_types`.
2. Implement `export_metrics_table`.
3. Implement `export_interaction_table`.
4. Implement robust loading of `classifier_matrices_{Radius}.npz`.
5. Implement manifest writing.
6. Implement `load_niche_result`.

Acceptance criteria:

- Fresh run writes `metrics_{Radius}.tsv` when `nico_result.score` exists.
- Fresh run writes `interactions_{Radius}.tsv`.
- Fresh run writes `niche_manifest_{Radius}.json`.
- `load_niche_result` can reconstruct a result without re-running NiCo.

### Phase 5: Implement CLI MVP

File:

```text
src/nico_wrapper/niche/cli.py
```

Update:

```text
pyproject.toml
README.md
```

Tasks:

1. Add Typer app.
2. Implement `run` command.
3. Implement `validate` command.
4. Implement `artifacts` command.
5. Add console entry point.
6. Add README section for `nico-niche`.

Acceptance criteria:

- `uv run nico-niche --help` works.
- `uv run nico-niche run --help` works.
- `uv run nico-niche validate --output-dir ...` validates only.
- `uv run nico-niche artifacts --output-dir ... --radius 0` reports artifact status.

### Phase 6: Implement plotting wrappers

File:

```text
src/nico_wrapper/niche/plotting.py
```

Tasks:

1. Build minimal NiCo plot namespace from `NicheInteractionResult`.
2. Implement plot dispatcher.
3. Implement `nico-niche plot` command.
4. Support `--kind` repeated option and `--kind all`.
5. Document limitations for ROC and predicted-probability plots.

Acceptance criteria:

- Confusion and coefficient plots can be generated from a loaded result.
- Fresh run with `--make-plots` generates selected plots.
- Unsupported plot kinds fail with clear errors rather than obscure attribute errors.

### Phase 7: Implement proximity wrapper

File:

```text
src/nico_wrapper/niche/proximity.py
```

CLI:

```text
nico-niche proximity
```

Tasks:

1. Wrap upstream proximity function.
2. Add seed handling.
3. Write observed and ratio TSV sidecars if practical.
4. Add CLI command.
5. Document that proximity is separate from logistic-regression niche interactions.

Acceptance criteria:

- Existing niche artifacts can be used to run proximity analysis.
- `--seed` makes the randomized baseline reproducible.
- Outputs are printed clearly.

### Phase 8: Testing

Add tests under:

```text
tests/niche/
```

Suggested tests:

1. Config defaults instantiate.
2. Invalid radius fails.
3. Invalid `k_fold` fails.
4. Missing h5ad fails.
5. Missing label key fails.
6. Missing spatial key fails.
7. Invalid coordinate shape fails.
8. Existing outputs fail without overwrite.
9. Radius tag normalization:
   - `"0" -> "0"`
   - `0 -> "0"`
   - `0.0 -> "0"` or reject ambiguous floats, depending on chosen policy.
10. `used_CT.txt` parser works.
11. Interaction table export works from synthetic classifier matrix.
12. Manifest write/load round trip works.
13. CLI help commands work through Typer runner.

Optional integration test:

- Create a tiny synthetic AnnData with enough cells/classes and run upstream if feasible. This may be slow or fragile due to Delaunay/logistic regression; keep as optional.

### Phase 9: Documentation

Update:

```text
README.md
docs/setup.md
```

Add:

1. Workflow section including `nico-niche` between transfer and future covariation.
2. CLI examples.
3. Python API examples.
4. Explanation of radius modes.
5. List of generated files.
6. Covariation readiness contract.
7. Known limitations inherited from NiCo.

Suggested README example:

```bash
uv run nico-transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis

uv run nico-niche run \
  --output-dir nico_analysis \
  --radius 0 \
  --epsilon-threshold 100

uv run nico-niche artifacts \
  --output-dir nico_analysis \
  --radius 0
```

---

## 15. Detailed CLI implementation sketch

```python
"""Typer command line app for NiCo niche interaction analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

from .config import (
    InteractionModelConfig,
    NeighborhoodConfig,
    NicheInteractionConfig,
    NichePlotConfig,
    ProximityConfig,
)
from .pipeline import run_niche_interactions, load_niche_result
from .validation import validate_niche_inputs, validate_covariation_artifacts

app = typer.Typer(
    name="nico-niche",
    help="Run NiCo niche interaction analysis after label transfer.",
    no_args_is_help=True,
)


@app.command("run")
def run(...):
    config = NicheInteractionConfig(...)
    result = run_niche_interactions(output_dir=output_dir, config=config)
    typer.echo(f"annotated_h5ad: {result.annotated_h5ad}")
    typer.echo(f"used_cell_types_tsv: {result.used_cell_types_tsv}")
    typer.echo(f"used_clusters_csv: {result.used_clusters_csv}")
    typer.echo(f"neighbors_pickle: {result.neighbors_pickle}")
    typer.echo(f"distances_pickle: {result.distances_pickle}")
    typer.echo(f"classifier_matrices_npz: {result.classifier_matrices_npz}")
    if result.metrics_tsv:
        typer.echo(f"metrics_tsv: {result.metrics_tsv}")
    if result.interactions_tsv:
        typer.echo(f"interactions_tsv: {result.interactions_tsv}")
    if result.manifest_json:
        typer.echo(f"manifest_json: {result.manifest_json}")
    typer.echo(f"covariation_ready: {str(result.covariation_ready).lower()}")
```

---

## 16. Risk register and mitigation

| Risk | Cause | Mitigation |
|---|---|---|
| Path bugs | Upstream concatenates strings. | Always pass trailing-separator directory strings to NiCo. |
| Artifact corruption | `used_CT.txt` is not radius-specific. | Strict overwrite validation and manifest warnings. |
| Multiple config runs in same output dir | Radius-specific files plus shared `used_CT.txt`. | Recommend separate output dirs for different exclusions/labels. |
| Upstream warnings suppression | Importing Interactions changes global warning behavior. | Delay import until analysis call; document side effect. |
| Delaunay failure | Degenerate geometry or too few points. | Pre-validate point count; let scipy errors surface with context. |
| No neighbors | Threshold/radius too strict. | Optional preflight or catch runtime errors and suggest larger radius/epsilon. |
| CV class-size failure | Too few cells per class. | Validate retained cell-type counts before run. |
| Hyperparameter search may run long | Upstream repeats until best C repeats. | Document; later maybe add timeout or custom implementation. |
| Metrics not persisted upstream | Upstream only returns metrics in memory. | Export wrapper `metrics_{Radius}.tsv` immediately after fresh run. |
| Selected C not persisted upstream | Upstream only returns it in memory. | Save in manifest. |
| ROC/probability plots unavailable after reload | Upstream doesn't save arrays. | Support only on fresh result or add optional persistence later. |
| Graphviz missing | Network plots use Graphviz/pydot. | Make plotting optional and fail with clear install message. |

---

## 17. Open design questions

1. Should `nico-niche run` generate plots by default?
   - Recommendation: no.

2. Should proximity analysis be part of `run` by default?
   - Recommendation: no.

3. Should the wrapper save ROC and predicted-probability arrays?
   - Recommendation: not in MVP; consider later if users need reproducible diagnostic plots from loaded results.

4. Should we expose cross-term models?
   - Recommendation: no, because downstream covariation expects linear artifacts.

5. Should we make a separate output directory per niche run?
   - Recommendation: keep upstream-compatible default layout, but consider a future `run_id` or `analysis_dir` mode if users need many configurations.

6. Should `radius=0.0` be accepted?
   - Recommendation: normalize to integer `0` to preserve downstream artifact names.

7. Should `NM` be exposed as configurable?
   - Upstream always excludes it. The wrapper should document this rather than pretending it is optional.

---

## 18. MVP definition

The minimum useful implementation is:

1. `NicheInteractionConfig`, `NeighborhoodConfig`, `InteractionModelConfig`, `NicheInteractionResult`.
2. `validate_niche_inputs`.
3. `run_niche_interactions` wrapping `nico.Interactions.spatial_neighborhood_analysis`.
4. Result path construction and artifact existence checks.
5. Metrics export from fresh runs.
6. Interaction table export from classifier matrix.
7. Manifest write/load.
8. CLI commands:
   - `nico-niche run`
   - `nico-niche validate`
   - `nico-niche artifacts`
9. README documentation.

Everything else can follow in later increments.

---

## 19. Future covariation integration notes

The future covariation wrapper should be able to consume `NicheInteractionResult` or the manifest directly.

Potential future API:

```python
from nico_wrapper.niche import load_niche_result
from nico_wrapper.covariation import run_covariation_analysis

niche = load_niche_result("nico_analysis", radius=0)

cov = run_covariation_analysis(
    niche_result=niche,
    ref_dir="inputRef",
    spatial_dir="inputQuery",
)
```

Because NiCo covariation expects the base output directory, the niche result should expose:

```python
niche.output_dir
niche.radius
```

and guarantee these files exist:

```text
niche.output_dir / "used_CT.txt"
niche.output_dir / f"used_Clusters{radius}.csv"
niche.output_dir / f"neighbors_{radius}.p"
niche.output_dir / f"distances_{radius}.p"
niche.output_dir / "niche_prediction_linear" / f"classifier_matrices_{radius}.npz"
```

This means the niche wrapper should prioritize artifact compatibility over a cleaner but incompatible directory layout.

---

## 20. Summary

The niche wrapper should provide a clear bridge between label transfer and covariation:

```bash
nico-transfer run --output-dir nico_analysis
nico-niche run --output-dir nico_analysis
nico-niche artifacts --output-dir nico_analysis --radius 0
```

It should improve NiCo's raw interaction module by adding:

1. typed configs;
2. preflight validation;
3. stable result objects;
4. overwrite protection;
5. manifest and table sidecars;
6. covariation artifact checks;
7. optional plotting/proximity commands.

The most important implementation constraint is to preserve NiCo's downstream artifact contract exactly, especially:

```text
niche_prediction_linear/classifier_matrices_{Radius}.npz
```

and the base-directory files:

```text
used_CT.txt
used_Clusters{Radius}.csv
neighbors_{Radius}.p
distances_{Radius}.p
```
