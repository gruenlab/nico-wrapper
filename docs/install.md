# Installation

`nico-wrapper` 1.0.0 supports **Python 3.11 on Linux**. The project metadata declares `>=3.11,<3.12`.

## Install from PyPI

With `uv`:

```bash
uv venv --python 3.11
uv pip install nico-wrapper
source .venv/bin/activate
```

Or install into an existing Python 3.11 environment:

```bash
python -m pip install nico-wrapper
```

Check the installed commands:

```bash
nico-preprocess --help
nico-transfer --help
nico-niche --help
nico-covariation --help
```

## Install from source

Use this path for local development or unreleased changes:

```bash
# Install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"

# Clone the public repository
git clone https://github.com/gruenlab/nico-wrapper.git
cd nico-wrapper

# Create the supported environment and install the locked dependencies
uv venv --python 3.11
uv sync
```

You can either activate the environment:

```bash
source .venv/bin/activate
```

or run commands from the repository through `uv run`.

## Requirements

- Linux
- Python 3.11
- `uv` for the source-development workflow
- a working C/C++ build toolchain if a scientific dependency has no compatible wheel
- the upstream `nico-sc-sp==1.6.0` package, installed automatically

Most dependencies should install from wheels on supported Python 3.11 Linux systems. System package names vary by distribution if a local build is required.

## NiCo dependency

This project depends on the upstream `nico-sc-sp` package through `pyproject.toml`:

```text
nico-sc-sp==1.6.0
```

You do not need to clone NiCo separately. The wrapper calls upstream NiCo modules internally, including:

- `nico.Annotations` for SCTransform, anchor discovery, and label transfer;
- `nico.Interactions` for spatial niche interaction analysis;
- `nico.Covariations` for latent-factor covariation analysis.

## Verify the installation

Check the package and core imports:

```bash
python - <<'PY'
from importlib.metadata import version

import nico
import nico_wrapper
import scanpy

print("nico-wrapper:", version("nico-wrapper"))
print("nico-sc-sp:", version("nico-sc-sp"))
print("scanpy:", scanpy.__version__)
print("imports: ok")
PY
```

When working from a source checkout without activating the environment, prefix the command with `uv run`.

Equivalent module entry points are also available:

```bash
python -m nico_wrapper.preprocess.cli --help
python -m nico_wrapper.transfer.cli --help
python -m nico_wrapper.niche.cli --help
python -m nico_wrapper.covariation.cli --help
```

## Starting input expectations

For the main workflow you need:

- a reference scRNA-seq dataset with cells as observations and genes as variables;
- a reference label column, default `.obs["cluster"]`;
- a spatial/query dataset with cells/spots as observations and genes as variables;
- spatial coordinates, default `.obsm["spatial"]`;
- raw count values in `.X`, unless you intentionally configure a normalization layer.

The preprocessing CLI can also convert two common raw formats:

- sparse reference triplet files: counts, genes, and barcodes;
- spatial count CSV plus coordinate CSV.

See [CLI usage](cli.md) for command-level details.

## Troubleshooting

### A scientific dependency attempts a source build

First confirm that the environment uses the supported Python version:

```bash
python --version
```

With `uv`, recreate the environment if necessary:

```bash
rm -rf .venv
uv venv --python 3.11
uv sync
```

If no wheel exists for your Linux distribution, install its standard compiler toolchain and the development libraries named by the failing package, then retry.

### CLI command not found

If you installed into `.venv`, activate it:

```bash
source .venv/bin/activate
which nico-preprocess
```

From a source checkout you can instead run:

```bash
uv run nico-preprocess --help
```
