"""I/O helpers for NiCo covariation artifacts."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import pickle
from typing import Any

import numpy as np

from nico_wrapper.niche.io import normalize_radius

from .config import CovariationConfig
from .results import CovariationArtifactPaths, CovariationResult

FEATURE_MATRIX_KEY = "weighted_neighborhood_of_factors_in_niche"


def planned_covariation_paths(
    output_dir: str | Path,
    *,
    radius: int | float | str = 0,
    n_factors: int = 3,
) -> CovariationArtifactPaths:
    """Return all standard artifact paths for a covariation run."""

    resolved_output_dir = Path(output_dir)
    radius_value, radius_tag = normalize_radius(radius)
    covariation_dir = resolved_output_dir / f"covariations_R{radius_tag}_F{n_factors}"
    return CovariationArtifactPaths(
        output_dir=resolved_output_dir,
        covariation_dir=covariation_dir,
        radius=radius_value,
        radius_tag=radius_tag,
        n_factors=n_factors,
        factors_pickle=covariation_dir / "factors_info.p",
        feature_matrix_npz=covariation_dir / "Principal_component_feature_matrix.npz",
        regression_dir=covariation_dir / "Regression_outputs",
        state_pickle=covariation_dir / "covariation_state.pkl",
        manifest_json=covariation_dir / "covariation_manifest.json",
        regression_tsv=covariation_dir / "regression_coefficients.tsv",
    )


def required_covariation_output_paths(paths: CovariationArtifactPaths) -> tuple[Path, ...]:
    """Return core files that prove upstream covariation completed."""

    return (paths.factors_pickle, paths.feature_matrix_npz)


def wrapper_output_paths(paths: CovariationArtifactPaths) -> tuple[Path, ...]:
    """Return wrapper sidecar paths."""

    return (paths.state_pickle, paths.manifest_json, paths.regression_tsv)


def resolve_ligand_receptor_db(path: str | Path | None = None) -> Path:
    """Resolve NiCo's ligand-receptor database path."""

    if path is not None:
        resolved = Path(path)
        if not resolved.exists() or not resolved.is_file():
            raise FileNotFoundError(f"Ligand-receptor DB does not exist or is not a file: {resolved}")
        return resolved

    try:
        import nico
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise FileNotFoundError("Could not import nico to locate bundled NiCoLRdb.txt.") from exc

    package_path = Path(nico.__file__).parent
    candidates = [
        package_path / "utils" / "NiCoLRdb.txt",
        package_path.parent / "utils" / "NiCoLRdb.txt",
        Path.cwd() / "NiCoLRdb.txt",
        Path.cwd() / "utils" / "NiCoLRdb.txt",
    ]
    for parent in Path.cwd().resolve().parents:
        candidates.extend(
            [
                parent / "NiCoLRdb.txt",
                parent / "nico_tutorial" / "NiCoLRdb.txt",
                parent / "nico_tutorial" / "NiCo" / "utils" / "NiCoLRdb.txt",
            ]
        )
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate
    raise FileNotFoundError("Could not locate NiCoLRdb.txt; pass ligand_receptor_db explicitly.")


def write_covariation_state(nico_result: Any, output_path: str | Path) -> Path:
    """Persist the raw NiCo covariation namespace needed by reports."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(nico_result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def read_covariation_state(path: str | Path) -> Any:
    """Read a raw NiCo covariation namespace from a state pickle."""

    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def write_covariation_manifest(
    result: CovariationResult,
    config: CovariationConfig,
    output_path: str | Path,
) -> Path:
    """Write a JSON manifest describing a covariation run."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "wrapper_version": _package_version("nico-wrapper"),
        "nico_package": _package_version("nico-sc-sp"),
        "output_dir": str(result.output_dir),
        "covariation_dir": str(result.covariation_dir),
        "radius": result.radius,
        "radius_tag": result.radius_tag,
        "n_factors": result.n_factors,
        "modality": result.modality,
        "factorization": config.factorization,
        "artifacts": {
            "factors_pickle": str(result.factors_pickle),
            "feature_matrix_npz": str(result.feature_matrix_npz),
            "regression_dir": str(result.regression_dir),
            "state_pickle": str(result.state_pickle) if result.state_pickle else None,
            "regression_tsv": str(result.regression_tsv) if result.regression_tsv else None,
            "manifest_json": str(path),
        },
        "niche": _niche_manifest_payload(result),
        "cell_type_names": result.cell_type_names,
        "config": asdict(config),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


def read_covariation_manifest(path: str | Path) -> dict[str, Any]:
    """Read a covariation manifest JSON file."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_factors_info(path: str | Path) -> tuple[Any, Any, Any]:
    """Load NiCo's ``factors_info.p`` tuple."""

    with Path(path).open("rb") as handle:
        value = pickle.load(handle)
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError(f"Expected factors_info.p to contain a 3-tuple; found {type(value).__name__}.")
    return value


def load_feature_matrix(path: str | Path) -> np.ndarray:
    """Load NiCo's weighted factor neighborhood feature matrix."""

    data = np.load(path, allow_pickle=True)
    if FEATURE_MATRIX_KEY not in data:
        raise ValueError(f"Feature matrix is missing key {FEATURE_MATRIX_KEY!r}: {path}")
    return np.asarray(data[FEATURE_MATRIX_KEY])


def _niche_manifest_payload(result: CovariationResult) -> dict[str, Any] | None:
    niche = result.niche_result
    if niche is None:
        return None
    return {
        "output_dir": str(niche.output_dir),
        "radius": niche.radius,
        "radius_tag": niche.radius_tag,
        "classifier_matrices_npz": str(niche.classifier_matrices_npz),
    }


def _package_version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None
