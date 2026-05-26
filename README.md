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

- [Installation](docs/install.md)
- [CLI usage](docs/cli.md)
- [Python API](docs/api.md)

## Quick start

```bash
uv venv --python 3.11
uv sync

uv run nico-preprocess --help
uv run nico-transfer --help
uv run nico-niche --help
uv run nico-covariation --help
```

A minimal end-to-end CLI skeleton is documented in [CLI usage](docs/cli.md#recommended-minimal-workflow).

## Main command groups

```text
nico-preprocess     Convert/build NiCo-ready reference and spatial inputs
nico-transfer       Transfer reference labels onto spatial/query cells
nico-niche          Run spatial niche interaction analysis
nico-covariation    Run latent-factor covariation analysis and reports
```

## License and citation

See this repository's license for wrapper code. For the underlying method and implementation, follow the citation guidance of upstream NiCo / `nico-sc-sp`.
