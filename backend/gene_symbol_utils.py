from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


KNOWN_BIOLOGY_GENES = {
    "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27",
    "IFI44", "IFI44L", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
    "STAT1", "STAT2", "IRF7", "IRF9", "CXCL10", "DDX58", "HERC6",
    "UBE2L6", "BST2", "RSAD2", "EPSTI1", "GBP1", "GBP2", "GBP4",
    "GBP5", "CCL2", "CCL3", "CCL4", "CCL8", "IL1B", "TNF",
    "NFKBIA", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAP2",
    "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7", "NKG7", "GNLY",
    "GZMB", "PRF1",
}


def is_ensembl_like(value: object) -> bool:
    value = str(value).strip().upper()
    return value.startswith("ENSG") or value.startswith("ENSMUSG")


def is_bad_symbol(value: object) -> bool:
    value = str(value).strip()
    return value == "" or value.lower() in {"nan", "none", "null", "na", "n/a"}


def clean_symbol(symbol: object, fallback: object) -> str:
    symbol = str(symbol).strip()
    fallback = str(fallback).strip()

    if is_bad_symbol(symbol):
        return fallback

    return symbol


def choose_best_symbol_column(var: pd.DataFrame) -> Optional[str]:
    """
    Choose the best gene-symbol column from adata.var.

    Scoring favors:
    - non-empty values
    - non-Ensembl-looking values
    - known biology marker overlap
    """

    candidate_cols = [
        "gene_symbol",
        "feature_name",
        "symbol",
        "gene_name",
        "name",
        "external_gene_name",
        "hgnc_symbol",
    ]

    scored = []

    for col in candidate_cols:
        if col not in var.columns:
            continue

        values = var[col].astype(str).fillna("").str.strip()
        valid = values[
            (values != "")
            & (values.str.lower() != "nan")
            & (values.str.lower() != "none")
            & (values.str.lower() != "null")
        ]

        if len(valid) == 0:
            continue

        non_ensembl_fraction = (~valid.map(is_ensembl_like)).mean()
        known_overlap = valid.str.upper().isin(KNOWN_BIOLOGY_GENES).sum()
        unique_fraction = valid.nunique() / max(len(valid), 1)

        score = (
            10.0 * float(non_ensembl_fraction)
            + 2.0 * float(known_overlap > 0)
            + 1.0 * float(unique_fraction)
        )

        scored.append((score, col, non_ensembl_fraction, known_overlap))

    if not scored:
        return None

    scored = sorted(scored, key=lambda x: x[0], reverse=True)
    return scored[0][1]


def build_gene_symbol_map(adata) -> Tuple[Dict[str, str], Optional[str]]:
    """
    Build a robust map from possible gene IDs to readable gene symbols.

    Uses:
    - adata.var_names
    - possible ID columns in adata.var
    - best available symbol column
    """

    if adata is None:
        return {}, None

    if not hasattr(adata, "var"):
        return {}, None

    var = adata.var.copy()

    symbol_col = choose_best_symbol_column(var)

    if symbol_col is None:
        return {}, None

    symbols = var[symbol_col].astype(str).map(
        lambda x: clean_symbol(x, x)
    )

    symbol_map: Dict[str, str] = {}

    # var_names -> symbols
    for key, symbol in zip(adata.var_names.astype(str), symbols):
        symbol_map[str(key)] = str(symbol)

    # extra possible ID columns -> symbols
    for id_col in [
        "gene_id",
        "feature_id",
        "ensembl_id",
        "id",
        "gene_ids",
    ]:
        if id_col in var.columns:
            for key, symbol in zip(var[id_col].astype(str), symbols):
                symbol_map[str(key)] = str(symbol)

    return symbol_map, symbol_col


def ensure_de_gene_symbols(
    de_results: pd.DataFrame,
    adata=None,
    gene_col: str = "names",
) -> pd.DataFrame:
    """
    Standardize a DE table.

    Guarantees:
    - gene_id: original gene identifier
    - gene_symbol: readable gene symbol where available
    - names: readable gene symbol for plots/reports
    """

    if de_results is None:
        return de_results

    if not hasattr(de_results, "empty") or de_results.empty:
        return de_results

    if gene_col not in de_results.columns:
        return de_results

    df = de_results.copy()

    original_names = df[gene_col].astype(str)

    if "gene_id" not in df.columns:
        df["gene_id"] = original_names

    # Already usable: names are not mostly Ensembl and contain some known genes.
    names_upper = original_names.str.upper()
    known_hits = names_upper.isin(KNOWN_BIOLOGY_GENES).sum()
    ensembl_fraction = original_names.map(is_ensembl_like).mean()

    if ensembl_fraction < 0.25 and known_hits > 0:
        df["gene_symbol"] = original_names
        df["names"] = original_names
        return df

    symbol_map, source_col = build_gene_symbol_map(adata)

    if not symbol_map:
        # Keep IDs, but still expose gene_symbol column for downstream code.
        df["gene_symbol"] = original_names
        df["names"] = original_names
        df.attrs["gene_symbol_source"] = None
        return df

    mapped = original_names.map(symbol_map)

    bad = (
        mapped.isna()
        | (mapped.astype(str).str.strip() == "")
        | (mapped.astype(str).str.lower() == "nan")
        | (mapped.astype(str).str.lower() == "none")
    )

    mapped.loc[bad] = original_names.loc[bad]

    df["gene_symbol"] = mapped.astype(str)
    df["names"] = df["gene_symbol"].astype(str)
    df.attrs["gene_symbol_source"] = source_col

    return df
