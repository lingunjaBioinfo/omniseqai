from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd


def _is_dataframe(obj: Any) -> bool:
    return hasattr(obj, "empty") and hasattr(obj, "columns")


def _pick_col(df: pd.DataFrame, candidates: List[str]):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _safe_int(value) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def _safe_float(value) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def _top_signature(summary: pd.DataFrame):
    if summary is None or summary.empty:
        return None

    if "detected" in summary.columns:
        detected = summary[summary["detected"] == True].copy()
    else:
        detected = summary.copy()

    if detected.empty:
        return None

    sort_cols = []
    ascending = []

    if "n_case_up_hits" in detected.columns:
        sort_cols.append("n_case_up_hits")
        ascending.append(False)

    if "case_up_hit_fraction" in detected.columns:
        sort_cols.append("case_up_hit_fraction")
        ascending.append(False)

    if sort_cols:
        detected = detected.sort_values(sort_cols, ascending=ascending)

    return detected.iloc[0].to_dict()


def _driver_celltypes(summary: pd.DataFrame, signature: str, top_n: int = 8) -> List[str]:
    if summary is None or summary.empty:
        return []

    required = {"signature", "scope", "cell_type", "n_case_up_hits"}

    if not required.issubset(set(summary.columns)):
        return []

    df = summary.copy()
    df = df[
        (df["signature"].astype(str) == str(signature))
        & (df["scope"].astype(str) == "celltype_specific")
    ].copy()

    if df.empty:
        return []

    df["n_case_up_hits"] = pd.to_numeric(
        df["n_case_up_hits"],
        errors="coerce",
    ).fillna(0)

    if "case_up_hit_fraction" in df.columns:
        df["case_up_hit_fraction"] = pd.to_numeric(
            df["case_up_hit_fraction"],
            errors="coerce",
        ).fillna(0)
    else:
        df["case_up_hit_fraction"] = 0.0

    df = df.sort_values(
        ["n_case_up_hits", "case_up_hit_fraction"],
        ascending=[False, False],
    )

    labels = []

    for _, row in df.head(top_n).iterrows():
        labels.append(
            f"{row['cell_type']} ({_safe_int(row['n_case_up_hits'])} genes)"
        )

    return labels


def _representative_genes(
    hits: pd.DataFrame,
    signature: str,
    top_n: int = 12,
) -> List[str]:
    if hits is None or not _is_dataframe(hits) or hits.empty:
        return []

    if "signature" not in hits.columns:
        return []

    df = hits[hits["signature"].astype(str) == str(signature)].copy()

    if df.empty:
        return []

    # Prefer whole-dataset evidence for the representative gene list.
    if "scope" in df.columns:
        whole = df[df["scope"].astype(str) == "whole_dataset"].copy()

        if not whole.empty:
            df = whole

    gene_col = _pick_col(
        df,
        [
            "gene_symbol",
            "gene",
            "names",
            "symbol",
            "gene_name",
            "gene_id",
        ],
    )

    if gene_col is None:
        return []

    logfc_col = _pick_col(
        df,
        [
            "logfoldchanges",
            "logFC",
            "log2FoldChange",
            "avg_log2FC",
        ],
    )

    padj_col = _pick_col(
        df,
        [
            "pvals_adj",
            "padj",
            "p_val_adj",
            "FDR",
            "fdr",
            "qval",
        ],
    )

    df[gene_col] = df[gene_col].astype(str)

    if logfc_col is not None:
        df[logfc_col] = pd.to_numeric(df[logfc_col], errors="coerce").fillna(0)
    else:
        df["_logfc_missing"] = 0
        logfc_col = "_logfc_missing"

    if padj_col is not None:
        df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce").fillna(1)
    else:
        df["_padj_missing"] = 1
        padj_col = "_padj_missing"

    df["_abs_lfc"] = df[logfc_col].abs()

    df = df.sort_values(
        [padj_col, "_abs_lfc"],
        ascending=[True, False],
    )

    genes = []

    for gene in df[gene_col].tolist():
        if gene not in genes:
            genes.append(gene)

        if len(genes) >= top_n:
            break

    return genes


def _unique_additional_signatures(
    summary: pd.DataFrame,
    top_signature: str,
    top_n: int = 5,
) -> List[str]:
    if summary is None or summary.empty:
        return []

    if "signature" not in summary.columns:
        return []

    if "detected" in summary.columns:
        df = summary[summary["detected"] == True].copy()
    else:
        df = summary.copy()

    if df.empty:
        return []

    if "n_case_up_hits" in df.columns:
        df["n_case_up_hits"] = pd.to_numeric(
            df["n_case_up_hits"],
            errors="coerce",
        ).fillna(0)
    else:
        df["n_case_up_hits"] = 0

    if "case_up_hit_fraction" in df.columns:
        df["case_up_hit_fraction"] = pd.to_numeric(
            df["case_up_hit_fraction"],
            errors="coerce",
        ).fillna(0)
    else:
        df["case_up_hit_fraction"] = 0

    df = df.sort_values(
        ["n_case_up_hits", "case_up_hit_fraction"],
        ascending=[False, False],
    )

    signatures = []

    for sig in df["signature"].astype(str).tolist():
        if sig == str(top_signature):
            continue

        if sig not in signatures:
            signatures.append(sig)

        if len(signatures) >= top_n:
            break

    return signatures


def format_biology_section(results: Dict[str, Any]) -> str:
    validation = results.get("biology_validation") or {}

    if not validation:
        return ""

    summary = validation.get("summary")
    hits = validation.get("hits")
    interpretation = validation.get("interpretation", [])

    lines = []
    lines.append("Biology Validation")
    lines.append("------------------")

    if summary is None or not _is_dataframe(summary) or summary.empty:
        lines.append("No biology validation results were available.")
        lines.append("")
        return "\n".join(lines)

    detected = summary.copy()

    if "detected" in detected.columns:
        detected = detected[detected["detected"] == True].copy()

    if detected.empty:
        lines.append("No predefined biological signature reached the detection threshold.")
        lines.append("")
        return "\n".join(lines)

    top = _top_signature(summary)

    if top:
        top_signature = str(top.get("signature", "unknown_signature"))
        top_description = str(top.get("description", ""))
        top_comparison = str(top.get("comparison", "unknown comparison"))
        top_hits = _safe_int(top.get("n_case_up_hits", 0))
        top_expected = _safe_int(top.get("n_expected_genes", 0))
        top_fraction = _safe_float(top.get("case_up_hit_fraction", 0.0))

        lines.append("Biological conclusion:")
        lines.append("")
        lines.append(
            f"- Dominant program: {top_signature}."
        )

        if top_description:
            lines.append(
                f"- Interpretation: {top_description}"
            )

        if top_expected > 0:
            lines.append(
                f"- Evidence strength: {top_hits}/{top_expected} expected genes detected "
                f"({top_fraction:.1%}) in {top_comparison}."
            )
        else:
            lines.append(
                f"- Evidence strength: {top_hits} expected genes detected in {top_comparison}."
            )

        driver_celltypes = _driver_celltypes(summary, top_signature)

        if driver_celltypes:
            lines.append(
                "- Main driver cell types: "
                + ", ".join(driver_celltypes)
                + "."
            )

        representative_genes = _representative_genes(hits, top_signature)

        if representative_genes:
            lines.append(
                "- Representative supporting genes: "
                + ", ".join(representative_genes)
                + "."
            )

        additional = _unique_additional_signatures(summary, top_signature)

        if additional:
            lines.append(
                "- Additional detected programs: "
                + ", ".join(additional)
                + "."
            )

        lines.append("")

    detected = detected.sort_values(
        ["n_case_up_hits", "case_up_hit_fraction"],
        ascending=[False, False],
    )

    lines.append("Detected biological programs:")
    lines.append("")

    for _, row in detected.head(10).iterrows():
        signature = row.get("signature", "unknown_signature")
        description = row.get("description", "")
        n_hits = _safe_int(row.get("n_case_up_hits", 0))
        n_expected = _safe_int(row.get("n_expected_genes", 0))
        scope = row.get("scope", "unknown_scope")
        cell_type = row.get("cell_type", "unknown_cell_type")
        comparison = row.get("comparison", "unknown_comparison")

        lines.append(
            f"- {signature}: {description} "
            f"({n_hits}/{n_expected} expected genes detected; "
            f"scope={scope}; cell_type={cell_type}; comparison={comparison})"
        )

    lines.append("")

    # Older validator-generated interpretation is intentionally omitted here.
    # The Biological conclusion section above is more specific and avoids redundancy.

    return "\n".join(lines)


def append_biology_section(report_text: str, results: Dict[str, Any]) -> str:
    section = format_biology_section(results)

    if not section:
        return report_text

    if "Biology Validation" in report_text:
        return report_text

    insertion_markers = [
        "\nIntegration and Batch QC\n",
        "\nTables Generated\n",
        "\nMethods\n",
    ]

    for marker in insertion_markers:
        if marker in report_text:
            return report_text.replace(marker, "\n" + section + marker, 1)

    return report_text.rstrip() + "\n\n" + section + "\n"
