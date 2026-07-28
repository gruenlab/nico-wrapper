# nico-wrapper

> [!WARNING]
> This is a proof-of-concept wrapper around NiCo. It can change at any moment, may contain bugs, and should not be treated as a stable production interface.

`nico-wrapper` provides Python APIs and Typer command-line tools for running the main NiCo (`nico-sc-sp`) workflow on reference single-cell RNA-seq data and spatial/Xenium query data.

It covers the practical pipeline around NiCo:

1. preprocessing raw reference and spatial data into NiCo-ready AnnData files;
2. transferring reference labels onto spatial/query cells;
3. running NiCo niche / spatial cell-type interaction analysis;
4. running NiCo latent-factor covariation analysis and selected exports/reports.

## Relationship to NiCo

This project is a wrapper around the upstream [`nico-sc-sp`](https://pypi.org/project/nico-sc-sp/) package, not a replacement for NiCo and not a fork of the NiCo method.

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

```bash
uv venv --python 3.11
uv pip install nico-wrapper
source .venv/bin/activate
```

Alternatively, install into an existing Python 3.11 environment:

```bash
python -m pip install nico-wrapper
```

## Quick start

```bash
nico-preprocess --help
nico-transfer --help
nico-niche --help
nico-covariation --help
```

A minimal end-to-end CLI skeleton is documented in [CLI usage](https://github.com/gruenlab/nico-wrapper/blob/main/docs/cli.md#recommended-minimal-workflow).

## Main command groups

```text
nico-preprocess     Convert/build NiCo-ready reference and spatial inputs
nico-transfer       Transfer reference labels onto spatial/query cells
nico-niche          Run spatial niche interaction analysis
nico-covariation    Run latent-factor covariation analysis and reports
```

## License and citation

`nico-wrapper` is licensed under the [MIT License](https://github.com/gruenlab/nico-wrapper/blob/main/LICENSE).

NiCo is used as an external library dependency. When using this wrapper, cite the upstream NiCo method:

> Agrawal A, Thomann S, Basu S, Grün D. NiCo identifies extrinsic drivers of cell state modulation by niche covariation analysis. *Nature Communications*. 2024;15:10628. [doi:10.1038/s41467-024-54973-w](https://doi.org/10.1038/s41467-024-54973-w).

## Contributors

See [CONTRIBUTORS.md](https://github.com/gruenlab/nico-wrapper/blob/main/CONTRIBUTORS.md).
