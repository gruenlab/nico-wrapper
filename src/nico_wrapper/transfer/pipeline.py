"""Public API for NiCo label transfer."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import os

from anndata import read_h5ad

from .config import (
    AnchorConfig,
    AnchorResult,
    AnnotationConfig,
    AnnotationResult,
    LabelTransferConfig,
    TieStrategy,
    TransferOutputs,
)
from .validation import (
    ValidationError,
    require_file,
    validate_label_transfer_config,
    validate_label_transfer_inputs,
    validate_output_files,
)


def run_label_transfer(
    ref_dir: str | Path,
    spatial_dir: str | Path,
    *,
    output_dir: str | Path,
    annotation_dir: str | Path | None = None,
    config: LabelTransferConfig = LabelTransferConfig(),
) -> TransferOutputs:
    """Run the complete NiCo label-transfer workflow.

    This is the main public API for transferring reference labels onto a
    preprocessed spatial/query dataset. Inputs are expected to be the output
    directories produced by ``nico-preprocess build``.

    Parameters
    ----------
    ref_dir
        Directory containing the preprocessed reference files, by default
        ``Original_counts.h5ad`` and ``sct_singleCell.h5ad``.
    spatial_dir
        Directory containing the preprocessed spatial/query file, by default
        ``sct_spatial.h5ad``.
    output_dir
        Directory where NiCo label-transfer outputs should be written.
    annotation_dir
        Optional directory for intermediate annotation artifacts. If omitted,
        ``output_dir / "annotations"`` is used.
    config
        Label-transfer configuration controlling anchor discovery, label
        propagation, output naming, overwrite behavior, and cleanup behavior.

    Returns
    -------
    TransferOutputs
        Paths to the annotated AnnData, anchors file, and per-iteration CSV
        outputs. ``anchors_npz`` is ``None`` if ``cleanup_intermediate`` removed
        the anchor file after a successful run.
    """

    validate_label_transfer_inputs(
        ref_dir=ref_dir,
        spatial_dir=spatial_dir,
        output_dir=output_dir,
        annotation_dir=annotation_dir,
        config=config,
    )

    resolved_output_dir = Path(output_dir)
    resolved_annotation_dir = Path(annotation_dir) if annotation_dir is not None else resolved_output_dir / "annotations"
    _ensure_directory(resolved_output_dir, label="output_dir")
    _ensure_directory(resolved_annotation_dir, label="annotation_dir")

    anchors = find_anchors(
        ref_dir=ref_dir,
        spatial_dir=spatial_dir,
        output_dir=resolved_output_dir,
        annotation_dir=resolved_annotation_dir,
        config=config.anchors,
        overwrite=config.overwrite,
    )
    annotation = transfer_labels(anchors, config=config.annotation)
    annotated_h5ad = save_transfer_result(
        annotation,
        output_dir=resolved_output_dir,
        output_h5ad_name=config.annotation.output_h5ad_name,
        output_label_key=config.annotation.output_label_key,
        overwrite=config.overwrite,
    )

    anchors_npz: Path | None = anchors.anchors_npz
    if config.cleanup_intermediate:
        _cleanup_intermediate_files(annotation)
        anchors_npz = None

    return TransferOutputs(
        output_dir=resolved_output_dir,
        annotation_dir=resolved_annotation_dir,
        annotated_h5ad=annotated_h5ad,
        anchors_npz=anchors_npz,
        iteration_cluster_csvs=annotation.iteration_cluster_csvs,
        iteration_celltype_csvs=annotation.iteration_celltype_csvs,
    )


def find_anchors(
    ref_dir: str | Path,
    spatial_dir: str | Path,
    *,
    output_dir: str | Path,
    annotation_dir: str | Path | None = None,
    config: AnchorConfig = AnchorConfig(),
    overwrite: bool = False,
) -> AnchorResult:
    """Find mutual-nearest-neighbor anchors between reference and spatial data.

    Parameters
    ----------
    ref_dir
        Directory containing the normalized reference AnnData and full reference
        AnnData. File names are taken from ``config``.
    spatial_dir
        Directory containing the normalized spatial/query AnnData. File names
        are taken from ``config``.
    output_dir
        Base output directory passed through to NiCo.
    annotation_dir
        Optional directory for anchor and annotation intermediates. If omitted,
        ``output_dir / "annotations"`` is used.
    config
        Anchor-finding configuration.
    overwrite
        Whether existing anchor/intermediate files may be replaced.

    Returns
    -------
    AnchorResult
        Wrapper containing NiCo's in-memory anchor object plus stable paths to
        the output and anchors artifact.
    """

    _validate_anchor_config(config)

    ref_dir = Path(ref_dir)
    spatial_dir = Path(spatial_dir)
    output_dir = Path(output_dir)
    resolved_annotation_dir = Path(annotation_dir) if annotation_dir is not None else output_dir / "annotations"

    sc_full_path = require_file(ref_dir / config.sc_full_filename, label="Full/original reference AnnData")
    sc_sct_path = require_file(ref_dir / config.sc_sct_filename, label="Normalized reference AnnData")
    spatial_sct_path = require_file(spatial_dir / config.spatial_sct_filename, label="Normalized spatial AnnData")
    _validate_anchor_shapes(
        sc_full_path=sc_full_path,
        sc_sct_path=sc_sct_path,
        spatial_sct_path=spatial_sct_path,
        neighbors=config.neighbors,
        n_pcs=config.n_pcs,
    )

    _ensure_directory(output_dir, label="output_dir")
    _ensure_directory(resolved_annotation_dir, label="annotation_dir")

    anchors_npz = resolved_annotation_dir / f"anchors_data_{config.neighbors}.npz"
    final_sct_sc = resolved_annotation_dir / "final_sct_sc.h5ad"
    final_sct_sp = resolved_annotation_dir / "final_sct_sp.h5ad"
    validate_output_files([anchors_npz, final_sct_sc, final_sct_sp], overwrite=overwrite)

    try:
        from nico import Annotations as sann
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo label transfer requires the 'nico' package.") from exc

    nico_result = sann.find_anchor_cells_between_ref_and_query(
        refpath=_as_nico_dir(ref_dir),
        quepath=_as_nico_dir(spatial_dir),
        output_annotation_dir=_as_nico_dir(resolved_annotation_dir),
        output_nico_dir=_as_nico_dir(output_dir),
        spatial_sct_anndata_filename=config.spatial_sct_filename,
        sc_sct_anndata_filename=config.sc_sct_filename,
        sc_full_anndata_filename=config.sc_full_filename,
        neigh=config.neighbors,
        no_of_pc=config.n_pcs,
        minkowski_order=config.minkowski_order,
    )

    _require_created_file(anchors_npz, label="NiCo anchors file")
    _require_created_file(final_sct_sc, label="NiCo final normalized reference intermediate")
    _require_created_file(final_sct_sp, label="NiCo final normalized spatial intermediate")

    return AnchorResult(
        nico_result=nico_result,
        output_dir=output_dir,
        annotation_dir=resolved_annotation_dir,
        anchors_npz=anchors_npz,
    )


def transfer_labels(
    anchors: AnchorResult,
    *,
    config: AnnotationConfig = AnnotationConfig(),
) -> AnnotationResult:
    """Propagate reference labels to spatial/query cells using NiCo.

    Parameters
    ----------
    anchors
        Result returned by :func:`find_anchors`.
    config
        Annotation and label-propagation configuration.

    Returns
    -------
    AnnotationResult
        Wrapper containing NiCo's in-memory annotation object and paths to the
        per-iteration annotation CSV files.
    """

    validate_label_transfer_config(LabelTransferConfig(anchors=AnchorConfig(), annotation=config))

    try:
        from nico import Annotations as sann
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise ImportError("NiCo label transfer requires the 'nico' package.") from exc

    nico_result = sann.nico_based_annotation(
        anchors.nico_result,
        ref_cluster_tag=config.ref_label_key,
        across_spatial_clusters_dispersion_cutoff=config.dispersion_cutoff,
        guiding_spatial_cluster_resolution_tag=config.spatial_cluster_key,
        number_of_iteration_to_perform_celltype_annotations=config.iterations,
        resolved_tie_issue_with_weighted_nearest_neighbor=_tie_strategy_for_nico(config.tie_strategy),
    )

    iteration_cluster_csvs = tuple(
        anchors.annotation_dir / f"{iteration}_nico_annotation_cluster.csv"
        for iteration in range(1, config.iterations + 1)
    )
    iteration_celltype_csvs = tuple(
        anchors.annotation_dir / f"{iteration}_nico_annotation_ct_name.csv"
        for iteration in range(1, config.iterations + 1)
    )
    for path in (*iteration_cluster_csvs, *iteration_celltype_csvs):
        _require_created_file(path, label="NiCo annotation CSV")

    if not hasattr(nico_result, "nico_cluster"):
        raise RuntimeError("NiCo annotation did not produce a 'nico_cluster' result.")
    if not hasattr(nico_result, "ad_sp_ori"):
        raise RuntimeError("NiCo annotation did not return the annotated spatial AnnData object ('ad_sp_ori').")

    return AnnotationResult(
        nico_result=nico_result,
        output_dir=anchors.output_dir,
        annotation_dir=anchors.annotation_dir,
        iteration_cluster_csvs=iteration_cluster_csvs,
        iteration_celltype_csvs=iteration_celltype_csvs,
    )


def save_transfer_result(
    annotation: AnnotationResult,
    *,
    output_dir: str | Path | None = None,
    output_h5ad_name: str = "nico_celltype_annotation.h5ad",
    output_label_key: str = "nico_ct",
    overwrite: bool = False,
) -> Path:
    """Save the annotated spatial AnnData object.

    Parameters
    ----------
    annotation
        Result returned by :func:`transfer_labels`.
    output_dir
        Optional output directory. If omitted, use ``annotation.output_dir``.
    output_h5ad_name
        File name for the annotated spatial AnnData.
    output_label_key
        ``.obs`` column to receive the transferred labels. Unlike NiCo's helper,
        this supports keys other than ``"nico_ct"``.
    overwrite
        Whether an existing output file or existing AnnData ``.obs`` column may
        be replaced.

    Returns
    -------
    pathlib.Path
        Path to the saved annotated AnnData object.
    """

    _validate_plain_output_filename(output_h5ad_name, label="output_h5ad_name")
    if not isinstance(output_label_key, str) or not output_label_key.strip():
        raise ValidationError("output_label_key must be a non-empty string.")

    resolved_output_dir = Path(output_dir) if output_dir is not None else annotation.output_dir
    _ensure_directory(resolved_output_dir, label="output_dir")
    output_path = resolved_output_dir / output_h5ad_name
    validate_output_files([output_path], overwrite=overwrite)

    nico_result = annotation.nico_result
    if not hasattr(nico_result, "ad_sp_ori") or not hasattr(nico_result, "nico_cluster"):
        raise RuntimeError("AnnotationResult does not contain NiCo spatial annotations to save.")

    adata = nico_result.ad_sp_ori.copy()
    if output_label_key in adata.obs and not overwrite:
        raise ValidationError(
            f"Output label key {output_label_key!r} already exists in annotated spatial .obs. "
            "Use overwrite=True to replace it."
        )
    if len(nico_result.nico_cluster) != adata.n_obs:
        raise RuntimeError(
            "NiCo annotation length does not match spatial AnnData rows: "
            f"{len(nico_result.nico_cluster)} labels for {adata.n_obs} cells."
        )

    adata.obs[output_label_key] = nico_result.nico_cluster
    adata.write_h5ad(output_path)
    return output_path


def _validate_anchor_config(config: AnchorConfig) -> None:
    validate_label_transfer_config(LabelTransferConfig(anchors=config, annotation=AnnotationConfig()))


def _ensure_directory(path: Path, *, label: str) -> None:
    if path.exists() and not path.is_dir():
        raise ValidationError(f"{label} exists but is not a directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _validate_anchor_shapes(
    *,
    sc_full_path: Path,
    sc_sct_path: Path,
    spatial_sct_path: Path,
    neighbors: int,
    n_pcs: int,
) -> None:
    try:
        sc_full = read_h5ad(sc_full_path)
        sc_sct = read_h5ad(sc_sct_path)
        spatial_sct = read_h5ad(spatial_sct_path)
    except Exception as exc:  # pragma: no cover - h5ad backend-specific
        raise ValidationError(f"Could not read transfer AnnData input: {exc}") from exc

    for label, adata in (
        ("Full/original reference", sc_full),
        ("Normalized reference", sc_sct),
        ("Normalized spatial", spatial_sct),
    ):
        if adata.n_obs == 0:
            raise ValidationError(f"{label} AnnData contains no cells.")
        if adata.n_vars == 0:
            raise ValidationError(f"{label} AnnData contains no genes.")
        if not adata.obs_names.is_unique:
            raise ValidationError(f"{label} AnnData .obs_names must be unique.")
        if not adata.var_names.is_unique:
            raise ValidationError(f"{label} AnnData .var_names must be unique.")

    if neighbors > sc_sct.n_obs:
        raise ValidationError(
            f"neighbors ({neighbors}) must be <= number of normalized reference cells ({sc_sct.n_obs})."
        )
    if neighbors > spatial_sct.n_obs:
        raise ValidationError(
            f"neighbors ({neighbors}) must be <= number of normalized spatial cells ({spatial_sct.n_obs})."
        )

    missing_reference_cells = sc_sct.obs_names.difference(sc_full.obs_names)
    if len(missing_reference_cells):
        raise ValidationError(
            "Normalized reference contains cells absent from full/original reference AnnData: "
            f"{len(missing_reference_cells)} missing."
        )

    shared_genes = sc_sct.var_names.intersection(spatial_sct.var_names)
    if len(shared_genes) == 0:
        raise ValidationError("Normalized reference and spatial AnnData files share no genes.")
    max_pcs = min(sc_sct.n_obs, len(shared_genes)) - 1
    if n_pcs > max_pcs:
        raise ValidationError(
            "n_pcs is too large for PCA on the normalized reference shared-gene matrix: "
            f"requested {n_pcs}, maximum allowed is {max_pcs} "
            f"(reference cells={sc_sct.n_obs}, shared genes={len(shared_genes)})."
        )


def _tie_strategy_for_nico(tie_strategy: TieStrategy | str) -> str:
    value = tie_strategy.value if isinstance(tie_strategy, TieStrategy) else tie_strategy
    if value == TieStrategy.MAJORITY.value:
        return "No"
    if value == TieStrategy.WEIGHTED.value:
        return "Yes"
    raise ValidationError("tie_strategy must be 'majority' or 'weighted'.")


def _as_nico_dir(path: str | Path) -> str:
    """Return a directory string with a trailing separator for NiCo.

    NiCo's annotation functions concatenate directory strings and file names, so
    callers must pass a trailing slash/backslash.
    """

    text = os.fspath(Path(path))
    return text if text.endswith(os.sep) else text + os.sep


def _cleanup_intermediate_files(annotation: AnnotationResult) -> None:
    paths = [
        annotation.annotation_dir / "final_sct_sc.h5ad",
        annotation.annotation_dir / "final_sct_sp.h5ad",
    ]
    fmnn = getattr(annotation.nico_result, "fmnn", None)
    if fmnn is not None:
        paths.append(Path(fmnn))
    for path in _unique_paths(paths):
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def _unique_paths(paths: Iterable[Path]) -> list[Path]:
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in paths:
        resolved = Path(path)
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return unique


def _require_created_file(path: Path, *, label: str) -> None:
    if not path.exists() or not path.is_file():
        raise RuntimeError(f"{label} was not created: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")


def _validate_plain_output_filename(value: str, *, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{label} must be a non-empty file name.")
    path = Path(value)
    if path.is_absolute() or path.name != value:
        raise ValidationError(f"{label} must be a plain file name, not a path: {value!r}.")
