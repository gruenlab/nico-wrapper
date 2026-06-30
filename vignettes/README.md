# Vignettes

This directory contains worked examples that illustrate how to run the full `nico-wrapper` analysis using different interfaces. They complement the reference documentation in [`../docs/api.md`](../docs/api.md) and [`../docs/cli.md`](../docs/cli.md).

## Files

| File | Interface | Description |
|------|-----------|-------------|
| `nico-wrapper_python_full_analysis.ipynb` | Python API | Interactive notebook covering the complete analysis pipeline using `nico-wrapper` Python functions. |
| `nico-wrapper_full_analysis.sh` | CLI | Shell script running the full analysis end-to-end via CLI commands, including covariation reports and exports. SLURM-ready for HPC batch submission. |
| `nico-wrapper_core_pipeline.sh` | CLI | Minimal shell script with just the four core pipeline steps (`nico-preprocess`, `nico-transfer`, `nico-niche`, `nico-covariation`). A lean starting point for scripted runs. |
| `nico-wrapper_python_data_exploration.ipynb` | CLI + Python API | Hybrid workflow: core pipeline steps are run via CLI commands, followed by interactive result exploration and plotting using Python API functions. |

## Documentation

Full package documentation lives in the [`../docs/`](../docs/) directory:

- [`install.md`](../docs/install.md) — installation instructions
- [`cli.md`](../docs/cli.md) — reference for all CLI commands and options
- [`api.md`](../docs/api.md) — Python API reference

