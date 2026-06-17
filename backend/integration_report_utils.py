from __future__ import annotations

from typing import Any, Dict


def format_integration_section(results: Dict[str, Any]) -> str:
    """
    Build a text section summarizing multi-sample and batch-aware QC.
    """

    integration_qc = results.get("integration_qc") or {}

    if not integration_qc:
        return ""

    status = integration_qc.get("status")
    message = integration_qc.get("message")
    sample_col = integration_qc.get("sample_col")
    batch_col = integration_qc.get("batch_col")
    condition_col = integration_qc.get("condition_col")
    n_samples = integration_qc.get("n_samples")
    n_batches = integration_qc.get("n_batches")
    method = integration_qc.get("method")

    lines = []
    lines.append("Integration and Batch QC")
    lines.append("------------------------")
    lines.append(f"Status: {status}")
    lines.append(f"Message: {message}")
    lines.append(f"Sample column: {sample_col}")
    lines.append(f"Batch column: {batch_col}")
    lines.append(f"Condition column: {condition_col}")
    lines.append(f"Number of samples: {n_samples}")
    lines.append(f"Number of batches: {n_batches}")
    lines.append(f"Integration method: {method}")
    lines.append("")

    if status == "integrated":
        lines.append(
            "Interpretation: batch-aware integration was completed successfully. "
            "Integrated UMAP figures should be used to inspect whether batch effects "
            "were reduced while biological structure was preserved."
        )

    elif status == "qc_only":
        lines.append(
            "Interpretation: multi-sample or batch-aware QC tables were exported, "
            "but integration was not performed. This can happen when no usable batch "
            "contrast exists or when the optional integration dependency is unavailable."
        )

    elif status == "skipped":
        lines.append(
            "Interpretation: no multi-sample or batch structure was detected, so "
            "integration was not required for this dataset."
        )

    else:
        lines.append(
            "Interpretation: integration status could not be fully determined. "
            "Review the exported integration_qc tables and run log."
        )

    lines.append("")

    return "\n".join(lines)


def append_integration_section(report_text: str, results: Dict[str, Any]) -> str:
    """
    Insert integration/QC section into a text report.
    """

    section = format_integration_section(results)

    if not section:
        return report_text

    if "Integration and Batch QC" in report_text:
        return report_text

    insertion_markers = [
        "\nTables Generated\n",
        "\nMethods\n",
        "\n## Methods\n",
        "\n# Methods\n",
    ]

    for marker in insertion_markers:
        if marker in report_text:
            return report_text.replace(
                marker,
                "\n" + section + marker,
                1,
            )

    return report_text.rstrip() + "\n\n" + section + "\n"
