from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd
import scanpy as sc
from anndata import AnnData


PathLike = Union[str, Path]


def _make_unique_and_preserve_symbols(adata: AnnData) -> AnnData:
    """
    Standardize AnnData after loading.

    Ensures:
    - obs names are unique
    - var names are unique
    - feature_name exists for gene-symbol mapping later
    """

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if "feature_name" not in adata.var.columns:
        adata.var["feature_name"] = adata.var_names.astype(str)

    return adata


def _read_csv_or_tsv(path: Path) -> AnnData:
    """
    Read a simple expression matrix.

    Expected common layout:
    - rows = genes
    - columns = cells/samples

    OmniSeqAI expects:
    - rows = cells
    - columns = genes

    So this loader transposes the matrix by default.
    """

    sep = "," if path.suffix.lower() == ".csv" else "\t"

    df = pd.read_csv(path, sep=sep, index_col=0)

    if df.empty:
        raise ValueError(f"Expression matrix is empty: {path}")

    # Convert all values to numeric where possible.
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Most expression matrices are genes x cells, so transpose.
    adata = AnnData(df.T)

    adata.var["feature_name"] = adata.var_names.astype(str)

    return adata


def _is_10x_mtx_folder(path: Path) -> bool:
    """
    Detect a 10x-style matrix folder.

    Supports compressed and uncompressed names.
    """

    required_any = [
        ["matrix.mtx", "matrix.mtx.gz"],
        ["barcodes.tsv", "barcodes.tsv.gz"],
        ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"],
    ]

    names = {p.name for p in path.iterdir()}

    for options in required_any:
        if not any(option in names for option in options):
            return False

    return True


def load_input(input_path: PathLike) -> AnnData:
    """
    Load supported single-cell input formats into AnnData.

    Supported:
    - .h5ad
    - 10x .h5
    - 10x mtx folder
    - .loom
    - .csv
    - .tsv
    - .txt

    Returns:
    - AnnData object
    """

    path = Path(input_path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")

    # --------------------------------------------------
    # 10x matrix directory
    # --------------------------------------------------
    if path.is_dir():
        if _is_10x_mtx_folder(path):
            print(f"Detected 10x matrix folder: {path}")

            adata = sc.read_10x_mtx(
                str(path),
                var_names="gene_symbols",
                cache=False,
            )

            return _make_unique_and_preserve_symbols(adata)

        raise ValueError(
            f"Unsupported directory format: {path}\n"
            "Expected a 10x folder containing matrix.mtx, barcodes.tsv, and features.tsv/genes.tsv."
        )

    suffix = path.suffix.lower()

    # --------------------------------------------------
    # AnnData
    # --------------------------------------------------
    if suffix == ".h5ad":
        print("Detected format: h5ad")
        adata = sc.read_h5ad(str(path))
        return _make_unique_and_preserve_symbols(adata)

    # --------------------------------------------------
    # 10x h5
    # --------------------------------------------------
    if suffix == ".h5":
        print("Detected format: 10x h5")

        try:
            adata = sc.read_10x_h5(str(path))
            return _make_unique_and_preserve_symbols(adata)
        except Exception as e:
            raise ValueError(
                f"Could not read .h5 file as 10x h5: {path}\n"
                f"Original error: {type(e).__name__}: {e}"
            )

    # --------------------------------------------------
    # Loom
    # --------------------------------------------------
    if suffix == ".loom":
        print("Detected format: loom")
        adata = sc.read_loom(str(path))
        return _make_unique_and_preserve_symbols(adata)

    # --------------------------------------------------
    # CSV / TSV / TXT expression matrix
    # --------------------------------------------------
    if suffix in {".csv", ".tsv", ".txt"}:
        print(f"Detected format: expression matrix ({suffix})")
        adata = _read_csv_or_tsv(path)
        return _make_unique_and_preserve_symbols(adata)

    raise ValueError(
        f"Unsupported input format: {path}\n"
        "Supported formats: .h5ad, .h5, 10x mtx folder, .loom, .csv, .tsv, .txt"
    )
