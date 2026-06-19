from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import anndata as ad
import pandas as pd
import scanpy as sc


SUPPORTED_SINGLE_FILE_SUFFIXES = {
    ".h5ad",
    ".h5",
    ".loom",
    ".csv",
    ".tsv",
    ".txt",
}


def load_input(input_path: str):
    """
    Load single-cell data from either:

    1. A single file:
       - .h5ad
       - 10x .h5
       - .loom
       - .csv/.tsv/.txt expression matrix

    2. A 10x MTX directory:
       - matrix.mtx / matrix.mtx.gz
       - barcodes.tsv / barcodes.tsv.gz
       - features.tsv or genes.tsv

    3. A multi-sample directory:
       - sample_1/sample_1.h5ad
       - sample_2/sample_2.h5ad
       - sample_3/filtered_feature_bc_matrix/
       - or direct files inside one folder

    Multi-sample directories are concatenated into one AnnData object with:
       - sample_id
       - batch_id
       - source_file
       - condition, if inferable or already present
    """

    path = Path(input_path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")

    if path.is_dir():
        if _looks_like_10x_mtx_dir(path):
            print(f"Detected format: 10x_mtx_directory")
            adata = sc.read_10x_mtx(str(path), var_names="gene_symbols", cache=False)
            adata = _finalize_adata(adata)
            return adata

        print(f"Detected format: multi_sample_directory")
        return _load_multi_sample_directory(path)

    detected = _detect_file_format(path)
    print(f"Detected format: {detected}")

    adata = _load_single_input(path)
    adata = _finalize_adata(adata)

    return adata


# ---------------------------------------------------------------------
# Single input loading
# ---------------------------------------------------------------------
def _detect_file_format(path: Path) -> str:
    suffix = path.suffix.lower()

    if suffix == ".h5ad":
        return "h5ad"

    if suffix == ".h5":
        return "10x_h5"

    if suffix == ".loom":
        return "loom"

    if suffix == ".csv":
        return "csv_matrix"

    if suffix == ".tsv":
        return "tsv_matrix"

    if suffix == ".txt":
        return "txt_matrix"

    raise ValueError(
        f"Unsupported input file format: {path}. "
        f"Supported suffixes: {sorted(SUPPORTED_SINGLE_FILE_SUFFIXES)}"
    )


def _load_single_input(path: Path):
    suffix = path.suffix.lower()

    if suffix == ".h5ad":
        return sc.read_h5ad(str(path))

    if suffix == ".h5":
        return sc.read_10x_h5(str(path))

    if suffix == ".loom":
        return sc.read_loom(str(path))

    if suffix in {".csv", ".tsv", ".txt"}:
        return _load_expression_matrix(path)

    raise ValueError(f"Unsupported input file format: {path}")


def _load_expression_matrix(path: Path):
    suffix = path.suffix.lower()

    if suffix == ".csv":
        df = pd.read_csv(path, index_col=0)

    else:
        df = pd.read_csv(path, index_col=0, sep="\t")

    # Heuristic:
    # If rows greatly exceed columns, assume genes x cells and transpose.
    # AnnData expects cells x genes.
    if df.shape[0] > df.shape[1]:
        df = df.T

    adata = ad.AnnData(df)

    return adata


# ---------------------------------------------------------------------
# Multi-sample directory loading
# ---------------------------------------------------------------------
def _load_multi_sample_directory(root: Path):
    sample_inputs = _discover_sample_inputs(root)

    if not sample_inputs:
        raise ValueError(
            f"No supported single-cell inputs found inside directory: {root}"
        )

    adatas = []
    sample_ids = []

    for sample_id, sample_path in sample_inputs:
        print(f"Loading sample: {sample_id} -> {sample_path}")

        if sample_path.is_dir() and _looks_like_10x_mtx_dir(sample_path):
            sample_adata = sc.read_10x_mtx(
                str(sample_path),
                var_names="gene_symbols",
                cache=False,
            )
        else:
            sample_adata = _load_single_input(sample_path)

        sample_adata = _finalize_adata(sample_adata)

        # loaded_sample_id records the file/folder we loaded from.
        # Do not overwrite biological sample_id if it already exists.
        sample_adata.obs["loaded_sample_id"] = str(sample_id)
        sample_adata.obs["source_file"] = str(sample_path)

        if (
            "sample_id" not in sample_adata.obs.columns
            or _obs_column_is_empty(sample_adata, "sample_id")
        ):
            sample_adata.obs["sample_id"] = str(sample_id)

        if (
            "batch_id" not in sample_adata.obs.columns
            or _obs_column_is_empty(sample_adata, "batch_id")
        ):
            sample_adata.obs["batch_id"] = str(sample_id)

        if "condition" not in sample_adata.obs.columns:
            inferred_condition = _infer_condition_from_name(sample_id)

            if inferred_condition is not None:
                sample_adata.obs["condition"] = inferred_condition
                sample_adata.obs["condition_raw"] = str(sample_id)

        adatas.append(sample_adata)
        sample_ids.append(sample_id)

    if len(adatas) == 1:
        return _finalize_adata(adatas[0])

    print(f"Concatenating {len(adatas)} samples...")

    combined = ad.concat(
        adatas,
        join="outer",
        label="loaded_sample_id",
        keys=sample_ids,
        index_unique="-",
        fill_value=0,
    )

    # Preserve explicit sample_id column if available.
    if "sample_id" not in combined.obs.columns:
        combined.obs["sample_id"] = combined.obs["loaded_sample_id"].astype(str)

    if "batch_id" not in combined.obs.columns:
        combined.obs["batch_id"] = combined.obs["sample_id"].astype(str)

    combined = _finalize_adata(combined)

    print(
        f"Combined dataset: {combined.n_obs:,} cells x {combined.n_vars:,} genes "
        f"from {len(sample_ids)} samples"
    )

    return combined


def _discover_sample_inputs(root: Path) -> List[Tuple[str, Path]]:
    """
    Discover sample inputs one level below root.

    Supported patterns:
    - root/sample1.h5ad
    - root/sample2.h5ad
    - root/sample1/filtered_feature_bc_matrix/
    - root/sample2/sample2.h5ad
    """

    discovered: List[Tuple[str, Path]] = []

    # Direct files in root.
    for item in sorted(root.iterdir()):
        if item.name.startswith("."):
            continue

        if item.is_file() and item.suffix.lower() in SUPPORTED_SINGLE_FILE_SUFFIXES:
            sample_id = item.stem
            discovered.append((sample_id, item))

    # Direct subdirectories.
    for item in sorted(root.iterdir()):
        if item.name.startswith("."):
            continue

        if not item.is_dir():
            continue

        if item.name in {"outputs", "runs", "reports", "__pycache__"}:
            continue

        if _looks_like_10x_mtx_dir(item):
            discovered.append((item.name, item))
            continue

        inner = _find_primary_input_inside_directory(item)

        if inner is not None:
            discovered.append((item.name, inner))

    # Deduplicate by path.
    seen = set()
    unique = []

    for sample_id, path in discovered:
        key = str(path.resolve())

        if key in seen:
            continue

        seen.add(key)
        unique.append((sample_id, path))

    return unique


def _find_primary_input_inside_directory(directory: Path) -> Optional[Path]:
    """
    Find the most likely single-cell input inside one sample directory.
    """

    if _looks_like_10x_mtx_dir(directory):
        return directory

    candidates = []

    for item in sorted(directory.iterdir()):
        if item.name.startswith("."):
            continue

        if item.is_file() and item.suffix.lower() in SUPPORTED_SINGLE_FILE_SUFFIXES:
            candidates.append(item)

        if item.is_dir() and _looks_like_10x_mtx_dir(item):
            candidates.append(item)

    if not candidates:
        return None

    priority = {
        ".h5ad": 0,
        ".h5": 1,
        ".loom": 2,
        ".csv": 3,
        ".tsv": 4,
        ".txt": 5,
    }

    def sort_key(path: Path):
        if path.is_dir():
            return -1

        return priority.get(path.suffix.lower(), 99)

    candidates = sorted(candidates, key=sort_key)

    return candidates[0]


# ---------------------------------------------------------------------
# 10x detection
# ---------------------------------------------------------------------
def _looks_like_10x_mtx_dir(path: Path) -> bool:
    if not path.is_dir():
        return False

    names = {p.name for p in path.iterdir() if p.is_file()}

    has_matrix = "matrix.mtx" in names or "matrix.mtx.gz" in names
    has_barcodes = "barcodes.tsv" in names or "barcodes.tsv.gz" in names

    has_features = (
        "features.tsv" in names
        or "features.tsv.gz" in names
        or "genes.tsv" in names
        or "genes.tsv.gz" in names
    )

    return has_matrix and has_barcodes and has_features


# ---------------------------------------------------------------------
# Metadata inference
# ---------------------------------------------------------------------
def _infer_condition_from_name(name: str) -> Optional[str]:
    """
    Conservative condition inference from sample/folder names.
    """

    x = str(name).lower()

    healthy_terms = [
        "healthy",
        "control",
        "ctrl",
        "normal",
        "untreated",
        "vehicle",
        "baseline",
    ]

    covid_terms = [
        "covid",
        "sarscov2",
        "sars_cov_2",
        "infected",
    ]

    ifn_terms = [
        "ifn",
        "ifnb",
        "ifn_beta",
        "stim",
        "stimulated",
        "treated",
    ]

    disease_terms = [
        "disease",
        "case",
        "tumor",
        "tumour",
        "cancer",
        "patient",
    ]

    if any(term in x for term in healthy_terms):
        return "Healthy"

    if any(term in x for term in covid_terms):
        return "COVID"

    if any(term in x for term in ifn_terms):
        return "IFN_beta"

    if any(term in x for term in disease_terms):
        return "Disease"

    return None


# ---------------------------------------------------------------------
# Final cleanup
# ---------------------------------------------------------------------
def _obs_column_is_empty(adata, col: str) -> bool:
    """
    Return True if an obs column is missing or contains only empty/null-like values.
    """

    if col not in adata.obs.columns:
        return True

    values = adata.obs[col]

    if values.isna().all():
        return True

    as_str = values.astype(str).str.strip().str.lower()

    empty_tokens = {"", "nan", "none", "null", "na", "n/a"}

    return bool(as_str.isin(empty_tokens).all())

def _finalize_adata(adata):
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if "feature_name" not in adata.var.columns:
        adata.var["feature_name"] = adata.var_names.astype(str)

    return adata
