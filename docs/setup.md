# nico-wrapper

This project wraps/uses `nico-sc-sp` and related single-cell analysis tooling. It is managed with [`uv`](https://docs.astral.sh/uv/) and requires **Python 3.11** (`pyproject.toml` declares `>=3.11,<3.12`).

## Quick start

```bash
# 1. Install uv if it is not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Restart your shell or add uv to PATH for the current shell
export PATH="$HOME/.local/bin:$PATH"

# 3. Clone the repository and enter it
git clone <repo-url> nico-wrapper
cd nico-wrapper

# 4. Create the virtual environment explicitly
uv venv --python 3.11

# 5. Install the project dependencies from pyproject.toml/uv.lock
uv sync

# 6. Optional: activate the environment for interactive shell use
source .venv/bin/activate
```

You can also run commands without activating the environment:

```bash
uv run python -c "import scanpy, pygraphviz, llvmlite; print('install ok')"
uv run jupyter lab
```

## System dependencies

Most Python packages install from prebuilt wheels. A few packages, especially `pygraphviz` and occasionally `llvmlite`, may need system development headers and compilers on Linux.

Install the OS packages below **before** running `uv sync`.

### Ubuntu / Debian

For Ubuntu 22.04/24.04 or recent Debian:

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

If Python 3.11 is not available from your distribution repositories, let `uv` download and manage Python automatically:

```bash
uv python install 3.11
uv venv --python 3.11
uv sync
```

Optional packages useful when `llvmlite` has to be built from source instead of using a wheel:

```bash
sudo apt install -y llvm-14 llvm-14-dev
export LLVM_CONFIG=/usr/bin/llvm-config-14
```

Normally this is **not required** because `llvmlite` should install from a wheel for Python 3.11 on common Linux architectures.

### RHEL, Rocky Linux, AlmaLinux, CentOS Stream, Fedora

Install compilers, Python headers, pkg-config, and Graphviz headers:

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

On RHEL-compatible distributions, `graphviz-devel` may require the CodeReady Builder/CRB repository.

RHEL 9:

```bash
sudo subscription-manager repos --enable codeready-builder-for-rhel-9-$(arch)-rpms
```

RHEL 8:

```bash
sudo subscription-manager repos --enable codeready-builder-for-rhel-8-$(arch)-rpms
```

Rocky Linux / AlmaLinux / CentOS Stream 9:

```bash
sudo dnf install -y dnf-plugins-core
sudo dnf config-manager --set-enabled crb
```

Rocky Linux / AlmaLinux 8:

```bash
sudo dnf install -y dnf-plugins-core
sudo dnf config-manager --set-enabled powertools || sudo dnf config-manager --set-enabled PowerTools
```

Fedora usually does not need extra repositories:

```bash
sudo dnf install -y graphviz graphviz-devel python3-devel gcc gcc-c++ make pkgconf-pkg-config
```

If Python 3.11 is not packaged for your RHEL variant, use `uv`'s managed Python instead:

```bash
uv python install 3.11
uv venv --python 3.11
uv sync
```

Optional packages useful only if `llvmlite` must be built from source:

```bash
sudo dnf install -y llvm14 llvm14-devel || sudo dnf install -y llvm-devel
export LLVM_CONFIG=$(command -v llvm-config-14 || command -v llvm-config)
```

Prefer the prebuilt `llvmlite` wheel when available. Building `llvmlite` from source requires an LLVM version compatible with the pinned `numba`/`llvmlite` versions.

## Installing uv

Recommended Linux installer:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
uv --version
```

Alternative installation methods:

```bash
# With pipx
pipx install uv

# With an existing Python/pip installation
python -m pip install --user uv
```

## Creating and using the virtual environment

From the repository root:

```bash
# Create .venv with Python 3.11
uv venv --python 3.11

# Install locked dependencies and the local project
uv sync

# Activate if desired
source .venv/bin/activate
```

If you want `uv sync` to recreate the environment from scratch:

```bash
rm -rf .venv
uv venv --python 3.11
uv sync
```

## Notes for `pygraphviz`

`pygraphviz` compiles against the system Graphviz C libraries. If installation fails with errors such as `graphviz/cgraph.h: No such file or directory` or `cannot find -lgvc`, install the Graphviz development package for your OS:

- Ubuntu/Debian: `sudo apt install graphviz graphviz-dev`
- RHEL/Fedora variants: `sudo dnf install graphviz graphviz-devel`

Then retry:

```bash
uv sync --reinstall-package pygraphviz
```

If Graphviz is installed in a non-standard location, expose its compiler/linker flags before retrying:

```bash
export CFLAGS="$(pkg-config --cflags libgvc)"
export LDFLAGS="$(pkg-config --libs libgvc)"
uv sync --reinstall-package pygraphviz
```

## Notes for `llvmlite` / `numba`

This project pins compatible versions of `numba` and `llvmlite` through `pyproject.toml`/`uv.lock`. On standard x86_64 Linux with Python 3.11, `llvmlite` should install from a prebuilt wheel.

If `llvmlite` attempts a source build and fails:

1. Confirm you are using Python 3.11:

   ```bash
   uv run python --version
   ```

2. Upgrade/retry with wheels enabled:

   ```bash
   uv sync --reinstall-package llvmlite
   ```

3. Only if no wheel is available for your platform, install a compatible LLVM development package and set `LLVM_CONFIG` as shown in the Ubuntu/RHEL sections above.

## Verify the installation

```bash
uv run python - <<'PY'
import sys
import scanpy
import pygraphviz
import llvmlite
import numba

print(sys.version)
print('scanpy:', scanpy.__version__)
print('pygraphviz:', pygraphviz.__version__)
print('llvmlite:', llvmlite.__version__)
print('numba:', numba.__version__)
PY
```

If this command completes without import errors, the environment is ready.
