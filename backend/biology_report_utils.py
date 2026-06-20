from __future__ import annotations

from typing import Any, Dict


def format_biology_section(results: Dict[str, Any]) -> str:
    validation = results.get("biology_validation") or {}

    if not validation:
        return ""

    summary = validation.get("summary")
    interpretation = validation.get("interpretation", [])

    lines = []
    lines.append("Biology Validation")
    lines.append("------------------")

    if summary is None or getattr(summary, "empty", True):
        lines.append("No biology validation results were available.")
        lines.append("")
        return "\n".join(lines)

    detected = summary[summary["detected"] == True].copy()

    if detected.empty:
        lines.append("No predefined biological signature reached the detection threshold.")
        lines.append("")
    else:
        detected = detected.sort_values(
            ["n_case_up_hits", "case_up_hit_fraction"],
            ascending=[False, False],
        )

        lines.append("Detected biological programs:")
        lines.append("")

        for _, row in detected.head(8).iterrows():
            lines.append(
                f"- {row['signature']}: {row['description']} "
                f"({int(row['n_case_up_hits'])}/{int(row['n_expected_genes'])} expected genes detected; "
                f"scope={row['scope']}; cell_type={row['cell_type']}; comparison={row['comparison']})"
            )

        lines.append("")

    if interpretation:
        lines.append("Interpretation:")
        lines.append("")

        for item in interpretation:
            lines.append(f"- {item}")

        lines.append("")

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
            return report_text.replace(
                marker,
                "\n" + section + marker,
                1,
            )

    return report_text.rstrip() + "\n\n" + section + "\n"
