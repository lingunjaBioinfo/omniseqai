from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd

from backend.biology_signatures import BIOLOGY_SIGNATURES, normalize_gene_symbol


def _pick_col(df: pd.DataFrame, candidates):
    lower_to_original = {str(c).lower(): c for c in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]

    return None


def load_user_signatures(path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """
    Load user-defined biology signatures from CSV/TSV.

    Required columns:
    - signature_name or signature or pathway or program
    - gene or gene_symbol or symbol

    Optional columns:
    - description
    - expected_direction

    Returns a dictionary compatible with BIOLOGY_SIGNATURES.
    """

    if not path:
        return {}

    p = Path(path)

    if not p.exists():
        raise FileNotFoundError(f"Signature file not found: {p}")

    if p.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(p, sep="\t")
    else:
        df = pd.read_csv(p)

    if df.empty:
        raise ValueError(f"Signature file is empty: {p}")

    sig_col = _pick_col(
        df,
        [
            "signature_name",
            "signature",
            "pathway",
            "program",
            "gene_set",
            "geneset",
        ],
    )

    gene_col = _pick_col(
        df,
        [
            "gene",
            "gene_symbol",
            "symbol",
            "gene_name",
        ],
    )

    desc_col = _pick_col(
        df,
        [
            "description",
            "desc",
            "label",
        ],
    )

    direction_col = _pick_col(
        df,
        [
            "expected_direction",
            "direction",
        ],
    )

    if sig_col is None:
        raise ValueError(
            "Signature file must contain one signature column: "
            "signature_name, signature, pathway, program, gene_set, or geneset."
        )

    if gene_col is None:
        raise ValueError(
            "Signature file must contain one gene column: "
            "gene, gene_symbol, symbol, or gene_name."
        )

    signatures: Dict[str, Dict[str, Any]] = {}

    for _, row in df.iterrows():
        signature_name = str(row[sig_col]).strip()

        if not signature_name or signature_name.lower() == "nan":
            continue

        gene = normalize_gene_symbol(row[gene_col])

        if not gene or gene.lower() == "nan":
            continue

        if signature_name not in signatures:
            description = signature_name.replace("_", " ")

            if desc_col is not None and pd.notna(row.get(desc_col)):
                description = str(row[desc_col]).strip() or description

            expected_direction = "case_up"

            if direction_col is not None and pd.notna(row.get(direction_col)):
                expected_direction = str(row[direction_col]).strip() or "case_up"

            signatures[signature_name] = {
                "description": description,
                "expected_direction": expected_direction,
                "genes": [],
                "source": "user",
            }

        if gene not in signatures[signature_name]["genes"]:
            signatures[signature_name]["genes"].append(gene)

    if not signatures:
        raise ValueError(f"No valid signatures found in: {p}")

    return signatures


def merge_builtin_and_user_signatures(
    user_signature_path: Optional[str] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Return built-in signatures plus optional user-defined signatures.

    If a user signature has the same name as a built-in signature, the user
    version is kept under user_<signature_name> to avoid silent overwriting.
    """

    merged: Dict[str, Dict[str, Any]] = {
        name: {
            **payload,
            "genes": list(payload.get("genes", [])),
            "source": payload.get("source", "built_in"),
        }
        for name, payload in BIOLOGY_SIGNATURES.items()
    }

    user_signatures = load_user_signatures(user_signature_path)

    for name, payload in user_signatures.items():
        final_name = name

        if final_name in merged:
            final_name = f"user_{name}"

        merged[final_name] = payload

    return merged
