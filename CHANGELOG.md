# Changelog

## 1.1.1 (unreleased)

### Added

- Expose `--include-rps-rpl-mt-genes / --exclude-rps-rpl-mt-genes` in the covariation `pathway`, `top-genes`, and `reports` commands. Inclusion remains the default. Exclusion removes species-specific ribosomal/mitochondrial symbol prefixes before top-gene selection and pathway enrichment; model fitting is unchanged. In report bundles, only `top-genes-all-factors` and `pathway` are affected, and neither is in the default bundle.
- Add `--organism` to covariation `top-genes` and `reports`, accepting lowercase `mouse` (default) or `human`. Invalid organisms fail before loading artifacts or generating reports.

### Fixed

- Accept lowercase `mouse` (default) and `human` organism inputs in covariation reports and the pathway CLI. Add a wrapper-local GSEApy adapter that fixes pathway enrichment casing while preserving NiCo's species-specific gene filtering.
- Preserve the selected organism, gene inclusion flag, and other report settings in `top-genes --all-factors` instead of recreating the configuration.

## 1.1.0

### Changed

- Require Python 3.12 instead of Python 3.11. Linux and macOS remain supported.
- Upgrade the upstream NiCo dependency from `nico-sc-sp==1.6.0` to `nico-sc-sp==1.8.0`.
- Update dependency requirements to `anndata>=0.13.1`, `numpy==2.5.3`, `pandas==3.0.6`, `scanpy==1.12.4`, and `scipy==1.18.1`.
- Adapt colocalization reports to NiCo 1.8.0. Scatterplots place neighbor factor loadings on the x-axis and central factor loadings on the y-axis, show colocalized and non-colocalized central cells, and fit the regression across all central cells. Bar and violin plots compare colocalized and non-colocalized loadings for both cell types.
- Use linear axes by default for colocalization scatterplots. Log-scaled axes can be selected with `axis_log_scale=True` in `plot_colocalized_factors` or `--axis-log-scale` in the `covariation colocalize` CLI command.
- Update the full-analysis and data-exploration notebook vignettes with revised covariation plotting examples.
- Update installation and CLI examples to use Python 3.12.

## 1.0.0

### Added

- Initial release with Python APIs and command-line tools for preprocessing, label transfer, niche analysis, and covariation analysis.
- Umbrella `nico-wrapper` command and standalone `nico-preprocess`, `nico-transfer`, `nico-niche`, and `nico-covariation` commands.
- Bundled ligand–receptor database, installation and usage documentation, and workflow vignettes.
- Python 3.11 support on Linux and macOS.
