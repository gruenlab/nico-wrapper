# nico-wrapper

`nico-wrapper` provides Python APIs and Typer command-line tools for running the main NiCo ([`nico-sc-sp`](https://github.com/ankitbioinfo/nico_tutorial)) workflow on reference single-cell RNA-seq data and spatial/Xenium query data.

[`nico-wrapper`](https://github.com/gruenlab/nico-wrapper) is being developed at the [Grün Lab](https://www.med.uni-wuerzburg.de/en/systemimmunologie/research/quantitative-single-cell-biology-of-the-immune-system-gruen-lab/) at the Institute for Systems Immunology under the leadership of Dominic Grün.

It covers the practical pipeline around NiCo:

1. preprocessing raw reference and spatial data into NiCo-ready AnnData files;
2. transferring reference labels onto spatial/query cells;
3. running NiCo niche / spatial cell-type interaction analysis;
4. running NiCo latent-factor covariation analysis and selected exports/reports.

## Relationship to NiCo

This project is a wrapper around the upstream [`nico-sc-sp`](https://github.com/ankitbioinfo/nico_tutorial) package, not a replacement for NiCo and not a fork of the NiCo method.

NiCo still performs the core scientific computations, including anchor discovery, label propagation, niche interaction modeling, and covariation analysis. `nico-wrapper` adds:

- repository-local CLI entry points;
- typed configuration objects and pipeline functions;
- input/output validation before expensive NiCo calls;
- predictable output locations;
- small sidecar exports such as manifests, metrics tables, and coefficient tables where available.

When you need methodological details, cite or consult NiCo. When you need a structured way to execute the workflow from this repository, use `nico-wrapper`.

## Documentation

- [Installation](https://github.com/gruenlab/nico-wrapper/blob/main/docs/install.md)
- [CLI usage](https://github.com/gruenlab/nico-wrapper/blob/main/docs/cli.md)
- [Python API](https://github.com/gruenlab/nico-wrapper/blob/main/docs/api.md)

## Installation

Version 1.0.0 supports Python 3.11 on Linux.

### 1. Install uv

Install `uv` with Homebrew, your system package manager, or the official standalone installer. See the [`uv` installation documentation](https://docs.astral.sh/uv/getting-started/installation/) for the available methods.

### 2. Run nico-wrapper with uvx

Run `nico-wrapper` directly in an isolated environment without installing it permanently:

```bash
uvx --python 3.11 nico-wrapper --help
uvx --python 3.11 nico-wrapper preprocess --help
```

`uvx` installs and caches the package. The explicit Python request makes `uv` select or download Python 3.11; if your default interpreter is already Python 3.11, you can omit `--python 3.11`.

### 3. Install nico-wrapper

Install `nico-wrapper` as a persistent command-line tool:

```bash
uv tool install --python 3.11 nico-wrapper
nico-wrapper --help
```

Alternatively, install it in a project-local virtual environment:

```bash
uv venv --python 3.11
source .venv/bin/activate
uv pip install nico-wrapper
```

## Quick start

With `uvx`, use the umbrella command groups:

```bash
uvx --python 3.11 nico-wrapper preprocess --help
uvx --python 3.11 nico-wrapper transfer --help
uvx --python 3.11 nico-wrapper niche --help
uvx --python 3.11 nico-wrapper covariation --help
```

A minimal end-to-end CLI skeleton is documented in [CLI usage](https://github.com/gruenlab/nico-wrapper/blob/main/docs/cli.md#recommended-minimal-workflow).

## Main command groups

| Command | Purpose |
|---|---|
| `nico-wrapper preprocess` | Convert/build NiCo-ready reference and spatial inputs |
| `nico-wrapper transfer` | Transfer reference labels onto spatial/query cells |
| `nico-wrapper niche` | Run spatial niche interaction analysis |
| `nico-wrapper covariation` | Run latent-factor covariation analysis and reports |

## Vignettes: start here

The [`vignettes/`](https://github.com/gruenlab/nico-wrapper/tree/main/vignettes) directory contains worked examples for different starting points:

| Goal | Start here |
|---|---|
| Run the minimal CLI workflow | [Core CLI pipeline](vignettes/nico-wrapper_core_pipeline.sh) |
| Run the complete CLI analysis with exports and plots | [Full CLI analysis](vignettes/nico-wrapper_full_analysis.sh) |
| Work interactively through the full Python API | [Full Python workflow](vignettes/nico-wrapper_python_full_analysis.ipynb) |
| Explore results from an existing run | [Interactive result exploration](vignettes/nico-wrapper_python_data_exploration.ipynb) |
| Use existing spatial cell-type annotations | [Starting from pre-annotated spatial data](vignettes/starting_from_annotated_spatial.md) |

New CLI users should begin with the **core pipeline**, which covers preprocessing, label transfer, niche analysis, and covariation analysis. The full CLI example adds exports, diagnostic plots, proximity analysis, covariation reports, ligand–receptor analysis, and pathway enrichment.

See the [CLI documentation](docs/cli.md) for command options and the [Python API documentation](docs/api.md) for programmatic usage.

## License and citation

`nico-wrapper` is licensed under the [MIT License](https://github.com/gruenlab/nico-wrapper/blob/main/LICENSE).

NiCo is used as an external library dependency. When using this wrapper, cite the upstream NiCo method:

> Agrawal A, Thomann S, Basu S, Grün D. NiCo identifies extrinsic drivers of cell state modulation by niche covariation analysis. *Nature Communications*. 2024;15:10628. [doi:10.1038/s41467-024-54973-w](https://doi.org/10.1038/s41467-024-54973-w).

## Contributors

See [CONTRIBUTORS.md](https://github.com/gruenlab/nico-wrapper/blob/main/CONTRIBUTORS.md).
