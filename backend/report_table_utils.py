from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List


def _display_table_path(path: str) -> str:
    """
    Convert internal output paths such as:
        outputs/tables/condition_de/file.csv

    into run-folder-relative report paths:
        tables/condition_de/file.csv
    """

    if not path:
        return ""

    p = Path(path)
    parts = list(p.parts)

    if "tables" in parts:
        idx = parts.index("tables")
        return Path(*parts[idx:]).as_posix()

    return p.as_posix()


def _describe_table(key: str, path: str) -> str:
    """
    Human-readable description for exported table paths.
    """

    key_lower = str(key).lower()
    path_lower = str(path).lower()

    if "run_summary" in key_lower or path_lower.endswith("run_summary.json"):
        return "Run metadata, selected mode, dataset dimensions, and detected metadata columns."

    if "celltype_counts_by_condition" in key_lower or "celltype_counts_by_condition" in path_lower:
        return "Cell counts for each cell type split by condition."

    if "celltype_proportions_by_condition" in key_lower or "celltype_proportions_by_condition" in path_lower:
        return "Cell-type proportions across conditions."

    if "celltype_counts" in key_lower or "celltype_counts.csv" in path_lower:
        return "Cell counts per annotated cell type."

    if "condition_de" in key_lower or "/condition_de/" in path_lower:
        if "significant" in key_lower or "significant" in path_lower:
            return "Significant whole-dataset differential expression genes."
        return "Full whole-dataset differential expression results."

    if "celltype_specific" in key_lower or "/celltype_specific_de/" in path_lower:
        if "significant" in key_lower or "significant" in path_lower:
            return "Significant cell-type-specific differential expression genes."
        return "Full cell-type-specific differential expression results."

    if "marker_genes" in key_lower or "/marker_genes/" in path_lower:
        if "all_marker_genes" in path_lower:
            return "All exploratory marker genes across groups."
        return "Exploratory marker genes for one group or cluster."

    if "integration_qc" in key_lower or "/integration_qc/" in path_lower:
        if "sample_counts" in key_lower or "sample_counts" in path_lower:
            return "Cell counts per sample."
        if "batch_counts" in key_lower or "batch_counts" in path_lower:
            return "Cell counts per batch."
        if "condition_counts" in key_lower or "condition_counts" in path_lower:
            return "Cell counts per condition."
        if "sample_condition_counts" in key_lower or "sample_condition_counts" in path_lower:
            return "Sample-by-condition cell count matrix."
        if "batch_condition_counts" in key_lower or "batch_condition_counts" in path_lower:
            return "Batch-by-condition cell count matrix."
        if "sample_qc_summary" in key_lower or "sample_qc_summary" in path_lower:
            return "Sample-level QC summary statistics."
        if "batch_qc_summary" in key_lower or "batch_qc_summary" in path_lower:
            return "Batch-level QC summary statistics."
        return "Multi-sample or batch-aware QC table."


    if path_lower.endswith(".csv"):
        return "Exported CSV analysis table."

    if path_lower.endswith(".json"):
        return "Exported JSON metadata file."

    return "Exported analysis table."


def _ordered_table_items(table_paths: Dict[str, str]) -> List[tuple[str, str]]:
    """
    Return table paths in a stable, human-friendly order.
    """

    priority = [
        "run_summary",
        "celltype_counts",
        "celltype_counts_by_condition",
        "celltype_proportions_by_condition",
        "integration_qc",
        "condition_de",
        "celltype_specific",
        "marker_genes",
    ]

    items = list(table_paths.items())

    def sort_key(item):
        key, path = item
        text = f"{key} {path}".lower()

        for idx, token in enumerate(priority):
            if token in text:
                return (idx, text)

        return (len(priority), text)

    return sorted(items, key=sort_key)


def format_table_section(results: Dict[str, Any]) -> str:
    """
    Build a text-report section summarizing exported tables.
    """

    table_paths = results.get("table_paths") or {}

    if not table_paths:
        return ""

    lines = []
    lines.append("Tables Generated")
    lines.append("----------------")
    lines.append(
        "The following machine-readable result tables were exported for downstream analysis."
    )
    lines.append("")

    for key, path in _ordered_table_items(table_paths):
        display_path = _display_table_path(path)
        description = _describe_table(key, path)
        lines.append(f"- {display_path}: {description}")

    lines.append("")

    return "\n".join(lines)


def append_table_section(report_text: str, results: Dict[str, Any]) -> str:
    """
    Insert the table summary section into an existing report.

    Preferred position:
    - before Methods section if present
    - otherwise append to the end
    """

    table_section = format_table_section(results)

    if not table_section:
        return report_text

    if "Tables Generated" in report_text:
        return report_text

    insertion_markers = [
        "\nMethods\n",
        "\n## Methods\n",
        "\n# Methods\n",
    ]

    for marker in insertion_markers:
        if marker in report_text:
            return report_text.replace(
                marker,
                "\n" + table_section + marker,
                1,
            )

    return report_text.rstrip() + "\n\n" + table_section + "\n"
