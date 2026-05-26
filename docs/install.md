# Installation

`nico-wrapper` is managed with [`uv`](https://docs.astral.sh/uv/) and requires **Python 3.11**. The project metadata declares `>=3.11,<3.12`.

## Requirements

- Python 3.11
- `uv`
- A working C/C++ build toolchain for packages that may compile locally
- Graphviz development headers for `pygraphviz`
- The Python dependencies pinned in `pyproject.toml` / `uv.lock`, including `nico-sc-sp==1.6.0`

## Quick install from source

```bash
# Install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"

# Clone this repository
git clone <repo-url> nico-wrapper
cd nico-wrapper

# Create the supported Python environment and install dependencies
uv venv --python 3.11
uv sync
```

You can either activate the environment:

```bash
source .venv/bin/activate
```

or run commands through `uv run` without activation.

## System dependencies

Most dependencies should install from wheels. `pygraphviz` needs Graphviz C libraries and headers. Install these before `uv sync` if your platform does not already provide them.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  pkg-config \
  python3.11 \
  python3.11-dev \
  graphviz \
  graphviz-dev
```

If Python 3.11 is not available from your OS repositories, let `uv` install it:

```bash
uv python install 3.11
uv venv --python 3.11
uv sync
```

### RHEL / Rocky Linux / AlmaLinux / CentOS Stream / Fedora

```bash
sudo dnf install -y \
  gcc \
  gcc-c++ \
  make \
  redhat-rpm-config \
  pkgconf-pkg-config \
  python3.11 \
  python3.11-devel \
  graphviz \
  graphviz-devel
```

## NiCo dependency

This repository depends on the upstream `nico-sc-sp` package through `pyproject.toml`:

```text
nico-sc-sp==1.6.0
```

`uv sync` installs NiCo into the local environment. You do not need to clone NiCo separately for normal wrapper usage.

The wrapper calls upstream NiCo Python modules internally, for example:

- `nico.Annotations` for SCTransform, anchor discovery, and label transfer;
- `nico.Interactions` for spatial niche interaction analysis;
- `nico.Covariations` for latent-factor covariation analysis.

## Verifying the installation

Check imports and key compiled dependencies:

```bash
uv run python - <<'PY'
import sys
import scanpy
import pygraphviz
import llvmlite
import numba
import nico

print(sys.version)
print('scanpy:', scanpy.__version__)
print('pygraphviz:', pygraphviz.__version__)
print('llvmlite:', llvmlite.__version__)
print('numba:', numba.__version__)
print('nico import ok')
PY
```

Check CLI entry points:

```bash
uv run nico-preprocess --help
uv run nico-transfer --help
uv run nico-niche --help
uv run nico-covariation --help
```

Equivalent module entry points are also available:

```bash
uv run python -m nico_wrapper.preprocess.cli --help
uv run python -m nico_wrapper.transfer.cli --help
uv run python -m nico_wrapper.niche.cli --help
uv run python -m nico_wrapper.covariation.cli --help
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

### `pygraphviz` cannot find Graphviz headers

Typical errors include `graphviz/cgraph.h: No such file or directory` or `cannot find -lgvc`.

Install Graphviz development headers, then retry:

```bash
# Ubuntu / Debian
sudo apt install graphviz graphviz-dev

# RHEL / Fedora variants
sudo dnf install graphviz graphviz-devel

uv sync --reinstall-package pygraphviz
```

If Graphviz is installed in a non-standard location, expose compiler flags:

```bash
export CFLAGS="$(pkg-config --cflags libgvc)"
export LDFLAGS="$(pkg-config --libs libgvc)"
uv sync --reinstall-package pygraphviz
```

### `llvmlite` / `numba` attempts a source build

The pinned `numba` / `llvmlite` versions should normally install from wheels on standard Python 3.11 Linux setups.

First confirm Python 3.11:

```bash
uv run python --version
```

Then retry:

```bash
uv sync --reinstall-package llvmlite
```

Only if no wheel is available for your platform, install a compatible LLVM package and set `LLVM_CONFIG`.

### CLI command not found

Run from the repository root with `uv run`:

```bash
uv run nico-preprocess --help
```

If using an activated shell, make sure the environment is active:

```bash
source .venv/bin/activate
which nico-preprocess
```
