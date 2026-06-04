from __future__ import annotations

from pathlib import Path

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Preformatted,
)


class RouterPDFReport:
    def save(
        self,
        report_text: str,
        filename: str = "reports/router_report.pdf"
    ) -> None:
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        doc = SimpleDocTemplate(
            str(path),
            pagesize=letter,
            rightMargin=0.7 * inch,
            leftMargin=0.7 * inch,
            topMargin=0.7 * inch,
            bottomMargin=0.7 * inch,
        )

        styles = getSampleStyleSheet()
        title_style = styles["Title"]
        heading_style = styles["Heading2"]
        body_style = styles["BodyText"]
        mono_style = ParagraphStyle(
            "Mono",
            parent=body_style,
            fontName="Courier",
            fontSize=8.5,
            leading=10,
            spaceAfter=6,
        )

        elements = []

        # Split report into sections and render in a readable way.
        lines = report_text.splitlines()
        buffer = []

        def flush_buffer():
            nonlocal buffer
            if buffer:
                elements.append(Preformatted("\n".join(buffer), mono_style))
                elements.append(Spacer(1, 0.12 * inch))
                buffer = []

        first_title_done = False

        for line in lines:
            stripped = line.strip()

            if stripped.startswith("==========") and not first_title_done:
                flush_buffer()
                elements.append(Paragraph("OmniSeqAI Router Report", title_style))
                elements.append(Spacer(1, 0.2 * inch))
                first_title_done = True
                continue

            if stripped.startswith("=====") or stripped.startswith("-----") or stripped.startswith("======"):
                flush_buffer()
                text = stripped.strip("=").strip("-").strip()
                if text:
                    elements.append(Paragraph(text, heading_style))
                    elements.append(Spacer(1, 0.08 * inch))
                continue

            if stripped == "":
                buffer.append("")
            else:
                buffer.append(line)

        flush_buffer()

        doc.build(elements)
        print(f"\nRouter PDF saved: {path}")
