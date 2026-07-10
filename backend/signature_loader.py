from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from backend.biology_signatures import BIOLOGY_SIGNATURES, normalize_gene_symbol


SIGNATURE_COLUMNS = [
    "signature_name",
    "signature",
    "pathway",
    "program",
    "gene_set",
    "geneset",
]

GENE_COLUMNS = [
    "gene",
    "gene_symbol",
    "symbol",
    "gene_name",
]

GENE_LIST_COLUMNS = [
    "genes",
    "gene_list",
    "gene_symbols",
    "symbols",
]

DESCRIPTION_COLUMNS = [
    "description",
    "desc",
    "label",
]

DIRECTION_COLUMNS = [
    "expected_direction",
    "direction",
]

VALID_DIRECTIONS = {
    "case_up",
    "up",
    "upregulated",
    "positive",
    "higher",
    "case_down",
    "down",
    "downregulated",
    "negative",
    "lower",
    "either",
    "any",
    "changed",
    "both",
}


def _pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lower_to_original = {str(col).lower().strip(): col for col in df.columns}

    for candidate in candidates:
        key = candidate.lower().strip()

        if key in lower_to_original:
            return lower_to_original[key]

    return None


def _read_signature_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")

    if suffix == ".csv":
        return pd.read_csv(path)

    raise ValueError(
        "Unsupported signature file format. "
        "Use .csv, .tsv, or .txt."
    )


def _split_genes(value: Any) -> List[str]:
    if pd.isna(value):
        return []

    text = str(value).strip()

    if not text or text.lower() == "nan":
        return []

    for sep in [";", "|"]:
        text = text.replace(sep, ",")

    genes = []

    for part in text.split(","):
        gene = normalize_gene_symbol(part)

        if gene and gene.lower() != "nan":
            genes.append(gene)

    return genes


def _clean_signature_name(value: Any) -> str:
    if pd.isna(value):
        return ""

    name = str(value).strip()

    if not name or name.lower() == "nan":
        return ""

    return name


def _clean_description(value: Any, fallback: str) -> str:
    if pd.isna(value):
        return fallback.replace("_", " ")

    description = str(value).strip()

    if not description or description.lower() == "nan":
        return fallback.replace("_", " ")

    return description


def _clean_expected_direction(value: Any) -> str:
    if pd.isna(value):
        return "case_up"

    direction = str(value).strip().lower()

    if not direction or direction == "nan":
        return "case_up"

    if direction not in VALID_DIRECTIONS:
        raise ValueError(
            f"Invalid expected_direction value: {direction}. "
            f"Accepted values are: {', '.join(sorted(VALID_DIRECTIONS))}"
        )

    return direction


def load_user_signatures(path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """
    Load user-defined biology signatures from CSV/TSV/TXT.

    Supported formats:

    1. One gene per row:
       signature_name,gene,description,expected_direction
       my_signature,ISG15,Interferon response,case_up
       my_signature,IFIT1,Interferon response,case_up

    2. Multiple genes in one cell:
       signature_name,genes,description,expected_direction
       my_signature,"ISG15,IFIT1,MX1",Interferon response,case_up

    Required:
    - one signature-name column
    - one gene column OR one gene-list column

    Optional:
    - description
    - expected_direction
    """

    if not path:
        return {}

    p = Path(path)

    if not p.exists():
        raise FileNotFoundError(f"Signature file not found: {p}")

    df = _read_signature_table(p)

    if df.empty:
        raise ValueError(f"Signature file is empty: {p}")

    sig_col = _pick_col(df, SIGNATURE_COLUMNS)
    gene_col = _pick_col(df, GENE_COLUMNS)
    gene_list_col = _pick_col(df, GENE_LIST_COLUMNS)
    desc_col = _pick_col(df, DESCRIPTION_COLUMNS)
    direction_col = _pick_col(df, DIRECTION_COLUMNS)

    if sig_col is None:
        raise ValueError(
            "Signature file must contain one signature column: "
            + ", ".join(SIGNATURE_COLUMNS)
            + "."
        )

    if gene_col is None and gene_list_col is None:
        raise ValueError(
            "Signature file must contain either one single-gene column "
            f"({', '.join(GENE_COLUMNS)}) or one gene-list column "
            f"({', '.join(GENE_LIST_COLUMNS)})."
        )

    signatures: Dict[str, Dict[str, Any]] = {}

    for row_index, row in df.iterrows():
        signature_name = _clean_signature_name(row.get(sig_col))

        if not signature_name:
            continue

        row_genes: List[str] = []

        if gene_col is not None:
            row_genes.extend(_split_genes(row.get(gene_col)))

        if gene_list_col is not None:
            row_genes.extend(_split_genes(row.get(gene_list_col)))

        row_genes = list(dict.fromkeys(row_genes))

        if not row_genes:
            continue

        if signature_name not in signatures:
            description = signature_name.replace("_", " ")

            if desc_col is not None:
                description = _clean_description(
                    row.get(desc_col),
                    fallback=signature_name,
                )

            expected_direction = "case_up"

            if direction_col is not None:
                expected_direction = _clean_expected_direction(
                    row.get(direction_col)
                )

            signatures[signature_name] = {
                "description": description,
                "expected_direction": expected_direction,
                "genes": [],
                "source": "user",
            }

        else:
            if direction_col is not None:
                row_direction = _clean_expected_direction(row.get(direction_col))
                stored_direction = signatures[signature_name]["expected_direction"]

                if row_direction != stored_direction:
                    raise ValueError(
                        "Inconsistent expected_direction values for signature "
                        f"'{signature_name}' near row {row_index + 2}: "
                        f"found '{row_direction}' but earlier rows used "
                        f"'{stored_direction}'."
                    )

        for gene in row_genes:
            if gene not in signatures[signature_name]["genes"]:
                signatures[signature_name]["genes"].append(gene)

    if not signatures:
        raise ValueError(f"No valid signatures found in: {p}")

    empty_signatures = [
        name
        for name, payload in signatures.items()
        if not payload.get("genes")
    ]

    if empty_signatures:
        raise ValueError(
            "The following signatures have no valid genes: "
            + ", ".join(empty_signatures)
        )

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
