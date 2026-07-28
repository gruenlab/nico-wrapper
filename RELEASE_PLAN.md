# PyPI Release Plan

## Current snapshot

- **Audit date:** 2026-07-28
- **Target distribution:** `nico-wrapper`
- **Branch:** `main` at `9ee200a` (`feat: add license statement`)
- **Remote state:** local `main` is one commit ahead of `origin/main`
- **Committed version:** `0.1.0`
- **Prepared release version:** `1.0.0` in the working tree
- **Verdict:** **Packaging and metadata preparation are complete.** Remaining owner actions are making the GitHub repository public, later committing/pushing/tagging the candidate, setting up manual index credentials, rehearsing on TestPyPI, and uploading to PyPI.

## Decisions

- [x] Institutional publication authorization is taken care of.
- [x] Use MIT for the wrapper; NiCo remains an external dependency.
- [x] Use `1.0.0` as the first public version.
- [x] Use an explicit Hatch allowlist for artifact contents.
- [x] Provide full PyPI-facing metadata.
- [x] Advertise Python 3.11 on Linux only.
- [x] List Alexander Dallmann as author and list Alexander Dallmann and Christian Eger as maintainers.
- [x] Use the public `gruenlab/nico-wrapper` GitHub repository for homepage, source, documentation, and issues.
- [x] Use a minimal direct runtime dependency declaration and rely on pinned `nico-sc-sp==1.6.0` for its scientific stack.
- [x] Do not add product tests or continuous-integration gates for this project.
- [x] Publish manually with API tokens rather than GitHub Actions/Trusted Publishing.
- [ ] Commit, push, and tag the final release later.

## Accomplished

### Package and documentation

- [x] Added wrapper modules for preprocessing, label transfer, niche analysis, and covariation analysis under `src/nico_wrapper/`.
- [x] Added four console entry points: `nico-preprocess`, `nico-transfer`, `nico-niche`, and `nico-covariation`.
- [x] Added installation, CLI, and Python API documentation under `docs/`.
- [x] Added a proof-of-concept/alpha warning and documented that NiCo remains the external scientific implementation.
- [x] Documented the Python 3.11/Linux support contract.
- [x] Added PyPI installation instructions using `uv` and `pip`.
- [x] Replaced relative README documentation links with public GitHub URLs.
- [x] Added the upstream NiCo citation and DOI.
- [x] Added `CONTRIBUTORS.md` listing Alexander Dallmann and Christian Eger, with a public README link.
- [x] Removed obsolete Graphviz/`pygraphviz` installation requirements after removing the unused direct dependency.

### License and release identity

- [x] Added the canonical MIT `LICENSE`, naming `Julius-Maximilians-Universität Würzburg` as copyright holder.
- [x] Added PEP 639 metadata: `license = "MIT"`, `license-files = ["LICENSE"]`, and `hatchling>=1.27`.
- [x] Updated the README license statement and NiCo relationship clarification.
- [x] Prepared version `1.0.0` in `pyproject.toml` and `uv.lock`.
- [x] Removed the tracked `src/nico_wrapper/preprocess/TODO.md` from the package. Its exact content is preserved as an intentionally unversioned root `TODO.md` through `.git/info/exclude`.

### Artifact selection

`pyproject.toml` now defines:

- [x] `packages = ["src/nico_wrapper"]` for the wheel.
- [x] An sdist `only-include` list containing `src`, `docs`, `README.md`, `LICENSE`, `CONTRIBUTORS.md`, and `pyproject.toml`.
- [x] Defense-in-depth exclusions for caches, virtual environments, build outputs, `dev-samples`, `scripts`, `.h5ad` files, and subagent artifacts.
- [x] Skipping traversal of excluded directories.

Hatch also includes standard sdist control/metadata files such as `.gitignore` and generated `PKG-INFO`. No unversioned root TODO, analysis, dataset, cache, transcript, or script enters either artifact.

### Public metadata

- [x] Added project URLs for homepage, repository, documentation, and issues.
- [x] Added alpha, science/bioinformatics, Python 3.11, and Linux classifiers.
- [x] Added bioinformatics, NiCo, single-cell, spatial-transcriptomics, and Xenium keywords.
- [x] Added maintainer metadata for Alexander Dallmann and Christian Eger.
- [x] Wheel metadata reports Metadata 2.4, `License-Expression: MIT`, and `License-File: LICENSE`.

### Dependency cleanup

The wrapper now declares only its direct runtime contract:

- `nico-sc-sp==1.6.0`
- `anndata>=0.8`
- `numpy==1.26.1`
- `pandas==2.1.1`
- `scanpy==1.11.2`
- `scipy==1.11.3`
- `typer>=0.15.0`

Completed cleanup:

- [x] Added direct `anndata` metadata.
- [x] Removed mandatory `jupyterlab` and `ipykernel`.
- [x] Removed unused direct `pygraphviz` and `llvmlite` declarations.
- [x] Removed redundant declarations already owned by pinned NiCo.
- [x] Regenerated `uv.lock`; the lock resolves 58 packages.

### Validation evidence

Validation was run on the current dirty checkout and on a clean exported candidate containing the same intended release changes:

- [x] Both builds produced byte-identical artifacts.
- [x] The sdist contains 38 files and is 66,977 bytes.
- [x] The wheel contains 34 files and is 76,048 bytes.
- [x] The sdist builds the wheel successfully.
- [x] Both artifacts pass `twine check`.
- [x] Artifact paths contain only the approved package/docs/metadata files.
- [x] Wheel metadata contains the intended version, dependencies, URLs, maintainers, keywords, classifiers, Python requirement, and license.
- [x] A clean Python 3.11 wheel installation resolved 56 packages and passed `uv pip check`.
- [x] The installed package imports and reports version `1.0.0` through `importlib.metadata`.
- [x] All four installed console scripts pass `--help`.
- [x] A no-dependency installation from the sdist builds and imports successfully.
- [x] `uv lock --check` passes.

## Remaining owner actions

### 1. Make the GitHub repository public

The selected metadata points to `https://github.com/gruenlab/nico-wrapper`, but the repository is not currently discoverable without authentication.

- [ ] Make the repository public before TestPyPI/PyPI publication.
- [ ] Verify the homepage, `docs/`, and issues URLs without signing in.
- [ ] Recheck all README links after the final branch is pushed.

### 2. Commit, push, and tag later

Current release work is not yet committed, and unrelated local analysis/development files remain in the worktree.

When ready:

- [ ] Separate or leave untracked unrelated local files.
- [ ] Review and commit `pyproject.toml`, `uv.lock`, `README.md`, `CONTRIBUTORS.md`, `docs/install.md`, `RELEASE_PLAN.md`, and the deletion of `src/nico_wrapper/preprocess/TODO.md`.
- [ ] Push `main`.
- [ ] Create annotated tag `v1.0.0` at the exact release commit.
- [ ] Confirm the tag version equals `pyproject.toml`.

The root `TODO.md` is intentionally local-only and must not be added.

### 3. Prepare manual TestPyPI and PyPI credentials

PyPI and TestPyPI are separate services with separate accounts and tokens.

- [ ] Verify account email and 2FA on both services.
- [ ] Confirm that `nico-wrapper` is available on both indexes.
- [ ] Create an account-scoped API token for the first upload if project-scoped tokens are unavailable before project creation.
- [ ] Store tokens in a password manager/keyring or inject them as environment variables; never commit them or place them directly in shell history.
- [ ] Use `__token__` as the Twine username.
- [ ] After first upload, replace broad bootstrap tokens with project-scoped tokens and revoke the broad tokens.

### 4. Rehearse on TestPyPI and publish

From a clean checkout of the tagged commit:

1. Run the local validation below.
2. Upload the two `dist/` files to TestPyPI:

   ```bash
   uvx twine upload --repository testpypi dist/*
   ```

3. Inspect the rendered page, metadata, links, license, and files.
4. Create a clean environment and install the production dependencies from PyPI. Then install only the wrapper artifact from TestPyPI, which is not a complete dependency mirror:

   ```bash
   uv venv --python 3.11 .testpypi-venv
   uv pip install --python .testpypi-venv/bin/python \
     "nico-sc-sp==1.6.0" "anndata>=0.8" "numpy==1.26.1" \
     "pandas==2.1.1" "scanpy==1.11.2" "scipy==1.11.3" \
     "typer>=0.15.0"
   uv pip install --python .testpypi-venv/bin/python \
     --index-url https://test.pypi.org/simple/ \
     --no-deps nico-wrapper==1.0.0
   ```

5. Run the version/import and four CLI `--help` checks in `.testpypi-venv`:

   ```bash
   .testpypi-venv/bin/python -c \
     'import nico_wrapper; from importlib.metadata import version; assert version("nico-wrapper") == "1.0.0"'
   .testpypi-venv/bin/nico-preprocess --help
   .testpypi-venv/bin/nico-transfer --help
   .testpypi-venv/bin/nico-niche --help
   .testpypi-venv/bin/nico-covariation --help
   ```
6. If correct, upload the exact same artifacts to production PyPI:

   ```bash
   uvx twine upload dist/*
   ```

If an uploaded TestPyPI artifact is wrong, increment the version before retrying; an index does not allow an uploaded filename to be reused.

## Local release validation

No product test suite or continuous-integration gate is required. Run these packaging and CLI checks manually from a clean checkout of the intended commit:

```bash
set -euo pipefail

# Remove outputs from a previous validation run before checking source state.
rm -rf \
  build dist \
  .release-wheel-venv .release-sdist-venv \
  release-sha256.txt

test -z "$(git status --porcelain)"
test "$(git describe --tags --exact-match)" = "v1.0.0"
test "$(grep '^version =' pyproject.toml)" = 'version = "1.0.0"'
uv lock --check

uv run --with build python -m build
uv run --with twine twine check dist/*
tar -tzf dist/nico_wrapper-*.tar.gz
unzip -l dist/nico_wrapper-*.whl
sha256sum dist/* | tee release-sha256.txt

rm -rf .release-wheel-venv
uv venv --python 3.11 .release-wheel-venv
uv pip install --python .release-wheel-venv/bin/python dist/nico_wrapper-*.whl
uv pip check --python .release-wheel-venv/bin/python
.release-wheel-venv/bin/python -c \
  'import nico_wrapper; from importlib.metadata import version; assert version("nico-wrapper") == "1.0.0"'
.release-wheel-venv/bin/nico-preprocess --help
.release-wheel-venv/bin/nico-transfer --help
.release-wheel-venv/bin/nico-niche --help
.release-wheel-venv/bin/nico-covariation --help

rm -rf .release-sdist-venv
uv venv --python 3.11 .release-sdist-venv
uv pip install --python .release-sdist-venv/bin/python \
  --no-deps dist/nico_wrapper-*.tar.gz
.release-sdist-venv/bin/python -c 'import nico_wrapper'

# Keep dist/ and the hash manifest; remove temporary installation environments.
rm -rf .release-wheel-venv .release-sdist-venv
```

Review every `tar`/`unzip` path before upload. Upload only the two intended files in `dist/`; keep `release-sha256.txt` outside `dist/`. A repeated validation run starts by removing these prior outputs.

## Final checklist

Before tagging:

- [x] Institutional MIT publication authorization is handled.
- [x] Artifact selection, metadata, documentation, dependencies, and lockfile are prepared and locally validated.
- [ ] GitHub repository and all package links are publicly reachable.
- [ ] The exact release commit is clean, reviewed, and pushed.

Before production upload:

- [ ] `v1.0.0` points to the intended commit and matches `pyproject.toml`.
- [ ] Local build, metadata, archive, installation, import, and CLI checks pass.
- [ ] TestPyPI rendering, installation, and smoke checks pass.
- [ ] Production upload uses the exact locally validated artifacts.

After publication:

- [ ] Production installation succeeds.
- [ ] The production page, files, metadata, links, and license are correct.
- [ ] Broad bootstrap tokens are replaced with project-scoped tokens where applicable.
