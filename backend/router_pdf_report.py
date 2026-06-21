from __future__ import annotations

from pathlib import Path
from datetime import datetime
from typing import Any, Dict, Optional, Sequence, Union
from xml.sax.saxutils import escape

import pandas as pd

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Image,
)

from backend.report_table_utils import format_table_section
from backend.integration_report_utils import format_integration_section
from backend.biology_report_utils import format_biology_section

class RouterPDFReport:

    def _safe(self, value) -> str:
        if value is None:
            return ""
        return str(value)

    def _page_num(self, canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.grey)
        canvas.drawRightString(
            doc.pagesize[0] - 0.7 * inch,
            0.45 * inch,
            f"Page {doc.page}",
        )
        canvas.restoreState()

    def _kv_table(self, items, body_style):
        rows = [
            [
                Paragraph(f"<b>{self._safe(k)}</b>", body_style),
                Paragraph(self._safe(v), body_style),
            ]
            for k, v in items
        ]

        table = Table(rows, colWidths=[2.1 * inch, 4.8 * inch])
        table.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("ROWBACKGROUNDS", (0, 0), (-1, -1), [colors.whitesmoke, colors.white]),
                    ("BOX", (0, 0), (-1, -1), 0.4, colors.lightgrey),
                    ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                    ("LEFTPADDING", (0, 0), (-1, -1), 5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ]
            )
        )
        return table

    def _df_table(self, df: Optional[pd.DataFrame], body_style, max_rows: int = 8):
        if df is None or df.empty:
            return Paragraph("No table available.", body_style)

        out = df.copy()

        # Keep only useful columns for readability
        preferred = ["names", "logfoldchanges", "pvals_adj", "scores"]
        if set(preferred).intersection(out.columns):
            cols = [c for c in preferred if c in out.columns]
            out = out[cols]

        out = out.head(max_rows)

        for col in out.columns:
            if pd.api.types.is_numeric_dtype(out[col]):
                out[col] = out[col].map(lambda x: f"{x:.3g}")

        headers = list(out.columns)
        rows = [[Paragraph(f"<b>{h}</b>", body_style) for h in headers]]

        for _, row in out.iterrows():
            rows.append([Paragraph(self._safe(v), body_style) for v in row.tolist()])

        col_widths = [6.7 * inch / len(headers)] * len(headers)

        table = Table(rows, colWidths=col_widths, repeatRows=1)
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9E2F3")),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("TOPPADDING", (0, 0), (-1, -1), 3),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ]
            )
        )
        return table

    def _figure(self, path, caption, max_width=6.8 * inch, max_height=4.8 * inch):
        path = Path(path)

        if not path.exists():
            return []

        img = Image(str(path))
        iw, ih = img.imageWidth, img.imageHeight

        scale = min(max_width / iw, max_height / ih)
        img.drawWidth = iw * scale
        img.drawHeight = ih * scale

        return [
            img,
            Spacer(1, 0.06 * inch),
            Paragraph(f"<i>{caption}</i>", self.body),
            Spacer(1, 0.2 * inch),
        ]

    def _bullets(self, items):
        flow = []
        for item in items:
            flow.append(Paragraph(f"• {self._safe(item)}", self.body))
            flow.append(Spacer(1, 0.035 * inch))
        return flow

    def _table_section(self, results: Dict[str, Any]):
        """
        Build the exported-tables section for the PDF report.
        """

        section_text = format_table_section(results)

        if not section_text:
            return []

        flow = []

        flow.append(Paragraph("Tables Generated", self.h1))
        flow.append(
            Paragraph(
                "The following machine-readable result tables were exported for downstream analysis.",
                self.body,
            )
        )
        flow.append(Spacer(1, 0.08 * inch))

        for raw_line in section_text.splitlines():
            line = raw_line.strip()

            if not line:
                continue

            if line == "Tables Generated":
                continue

            if set(line) == {"-"}:
                continue

            if line.startswith("The following machine-readable"):
                continue

            if line.startswith("- "):
                item = line[2:].strip()

                if ": " in item:
                    path, description = item.split(": ", 1)
                    paragraph_text = (
                        f"• <b>{escape(path)}</b>: {escape(description)}"
                    )
                else:
                    paragraph_text = f"• {escape(item)}"

                flow.append(Paragraph(paragraph_text, self.body))
                flow.append(Spacer(1, 0.035 * inch))

        if len(flow) <= 3:
            return []

        return flow

    def _biology_validation_section(self, results: Dict[str, Any]):
        """
        Build the biology validation section for the PDF report.
        """

        section_text = format_biology_section(results)

        if not section_text:
            return []

        flow = []

        lines = [line.strip() for line in section_text.splitlines() if line.strip()]

        flow.append(Paragraph("Biology Validation", self.h1))

        for line in lines:
            if line == "Biology Validation":
                continue

            if set(line) == {"-"}:
                continue

            if line.startswith("- "):
                flow.append(Paragraph("• " + escape(line[2:]), self.body))
                flow.append(Spacer(1, 0.035 * inch))

            else:
                flow.append(Paragraph(escape(line), self.body))
                flow.append(Spacer(1, 0.035 * inch))

        return flow

    def _integration_qc_section(self, results: Dict[str, Any]):
        """
        Build the integration/QC section for the PDF report.
        """

        section_text = format_integration_section(results)

        if not section_text:
            return []

        flow = []

        lines = [line.strip() for line in section_text.splitlines() if line.strip()]

        flow.append(Paragraph("Integration and Batch QC", self.h1))

        for line in lines:
            if line == "Integration and Batch QC":
                continue

            if set(line) == {"-"}:
                continue

            if ":" in line and not line.startswith("Interpretation:"):
                key, value = line.split(":", 1)
                paragraph_text = f"<b>{escape(key.strip())}</b>: {escape(value.strip())}"
                flow.append(Paragraph(paragraph_text, self.body))
                flow.append(Spacer(1, 0.035 * inch))

            elif line.startswith("Interpretation:"):
                flow.append(Spacer(1, 0.08 * inch))
                flow.append(Paragraph(escape(line), self.body))

            else:
                flow.append(Paragraph(escape(line), self.body))

        return flow

    def _build_from_results(self, results: Dict[str, Any]):
        styles = getSampleStyleSheet()

        self.title = ParagraphStyle(
            "OmniTitle",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=21,
            leading=25,
            spaceAfter=12,
        )

        self.h1 = ParagraphStyle(
            "OmniH1",
            parent=styles["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=15,
            leading=18,
            spaceBefore=10,
            spaceAfter=7,
        )

        self.h2 = ParagraphStyle(
            "OmniH2",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=12,
            leading=14,
            spaceBefore=7,
            spaceAfter=5,
        )

        self.body = ParagraphStyle(
            "OmniBody",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=9,
            leading=11.5,
            spaceAfter=4,
        )

        story = []

        profile = results.get("profile", {}) or {}
        plan = results.get("plan", None)
        adata = results.get("adata", None)
        decision = results.get("decision", None)
        condition_de_results = results.get("condition_de_results", {}) or {}
        celltype_specific = results.get("celltype_specific", {}) or {}
        celltype_comparison = results.get("celltype_comparison", None)
        figure_paths = results.get("figure_paths", {}) or {}

        story.append(Paragraph("OmniSeqAI Analysis Report", self.title))
        story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}", self.body))

        if decision is not None:
            story.append(Paragraph(f"Selected mode: <b>{decision.mode}</b>", self.body))
            story.append(Paragraph(self._safe(decision.reason), self.body))

        story.append(Spacer(1, 0.15 * inch))

        story.append(Paragraph("Executive Summary", self.h1))

        summary = []

        if adata is not None:
            summary.append(f"Dataset size: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        summary.append(
            "Detected metadata: "
            f"cell type = {profile.get('cell_type_column')}, "
            f"condition = {profile.get('condition_column')}, "
            f"sample = {profile.get('sample_column')}, "
            f"patient = {profile.get('patient_column')}"
        )

        if plan is not None:
            mode = "pseudobulk" if getattr(plan, "use_pseudobulk", False) else "cell-level"
            summary.append(f"Primary analysis mode: {mode}")
            summary.append(f"Baseline/control: {getattr(plan, 'baseline', None)}")
            summary.append(f"Comparisons: {getattr(plan, 'pairwise_comparisons', None)}")

        if condition_de_results:
            pair, info = next(iter(condition_de_results.items()))
            summary.append(
                f"Whole-dataset DE: {pair[0]} vs {pair[1]} "
                f"({info.get('mode')}, {info.get('n_sig_genes')} significant genes)"
            )
            for item in info.get("interpretation", []):
                summary.append(item)

        story.extend(self._bullets(summary))

        story.append(Paragraph("Dataset Overview", self.h1))

        overview = [
            ("Cells", f"{adata.n_obs:,}" if adata is not None else ""),
            ("Genes", f"{adata.n_vars:,}" if adata is not None else ""),
            ("Cell type column", profile.get("cell_type_column")),
            ("Condition column", profile.get("condition_column")),
            ("Sample column", profile.get("sample_column")),
            ("Patient column", profile.get("patient_column")),
            ("Batch column", profile.get("batch_column")),
            ("Gene symbol source", profile.get("gene_symbol_source")),
        ]
        story.append(self._kv_table(overview, self.body))
        story.append(Spacer(1, 0.18 * inch))

        story.append(PageBreak())

        story.append(Paragraph("Figures", self.h1))

        figure_order = [
            (
                "umap_celltype",
                "Figure 1. UMAP annotated by curated cell-type labels.",
            ),
            (
                "umap_condition",
                "Figure 2. UMAP colored by experimental condition.",
            ),
            (
                "celltype_proportions",
                "Figure 3. Cell-type proportions across conditions.",
            ),
            (
                "umap_sample",
                "UMAP colored by sample.",
            ),
            (
                "umap_batch",
                "UMAP colored by batch.",
            ),
            (
                "umap_integrated_sample",
                "Integrated UMAP colored by sample.",
            ),
            (
                "umap_integrated_batch",
                "Integrated UMAP colored by batch.",
            ),
            (
                "umap_integrated_condition",
                "Integrated UMAP colored by condition.",
            ),
            (
                "volcano",
                "Figure 4. Volcano plot showing condition-associated differential expression.",
            ),
            (
                "biology_signature_hits",
                "Biology validation summary showing detected signature evidence.",
            ),
            (
                "biology_celltype_signature_heatmap",
                "Cell-type localization of detected biological signatures.",
            ),
            (
                "pseudobulk_heatmap",
                "Figure 5. Clustered heatmap of top pseudobulk differential genes.",
            ),
        ]

        for key, caption in figure_order:
            path = figure_paths.get(key)
            if path:
                story.extend(self._figure(path, caption))

        story.append(PageBreak())

        story.append(Paragraph("Whole-Dataset Differential Expression", self.h1))

        if condition_de_results:
            for pair, info in condition_de_results.items():
                story.append(Paragraph(f"{pair[0]} vs {pair[1]}", self.h2))

                story.append(
                    self._kv_table(
                        [
                            ("Mode", info.get("mode")),
                            ("Status", info.get("status")),
                            ("Significant genes", info.get("n_sig_genes")),
                        ],
                        self.body,
                    )
                )
                story.append(Spacer(1, 0.08 * inch))
                story.append(Paragraph("Top DE genes", self.h2))
                story.append(self._df_table(info.get("de_results"), self.body, max_rows=10))

                interpretation = info.get("interpretation", [])
                if interpretation:
                    story.append(Paragraph("Interpretation", self.h2))
                    story.extend(self._bullets(interpretation))
                break
        else:
            story.append(Paragraph("No whole-dataset DE results available.", self.body))

        story.append(PageBreak())

        story.append(Paragraph("Cell-Type-Specific Differential Expression", self.h1))

        ranked = []

        for cell_type, pair_map in celltype_specific.items():
            for pair, info in pair_map.items():
                if not isinstance(info, dict):
                    continue
                ranked.append(
                    (
                        int(info.get("n_sig_genes", 0) or 0),
                        cell_type,
                        pair,
                        info,
                    )
                )

        ranked = sorted(ranked, reverse=True)[:5]

        for _, cell_type, pair, info in ranked:
            story.append(Paragraph(f"{cell_type}: {pair[0]} vs {pair[1]}", self.h2))
            story.append(
                self._kv_table(
                    [
                        ("Mode", info.get("mode")),
                        ("Status", info.get("status")),
                        ("Condition counts", info.get("condition_counts")),
                        ("Significant genes", info.get("n_sig_genes")),
                    ],
                    self.body,
                )
            )
            story.append(Spacer(1, 0.06 * inch))
            story.append(self._df_table(info.get("de_results"), self.body, max_rows=5))

            interpretation = info.get("interpretation", [])
            if interpretation:
                story.extend(self._bullets(interpretation))

            story.append(Spacer(1, 0.12 * inch))

        # --------------------------------------------------
        # These sections must be OUTSIDE the cell-type loop.
        # --------------------------------------------------
        biology_flow = self._biology_validation_section(results)

        if biology_flow:
            story.append(PageBreak())
            story.extend(biology_flow)

        integration_flow = self._integration_qc_section(results)

        if integration_flow:
            story.append(PageBreak())
            story.extend(integration_flow)

        table_flow = self._table_section(results)

        if table_flow:
            story.append(PageBreak())
            story.extend(table_flow)

        story.append(PageBreak())

        story.append(Paragraph("Methods", self.h1))
        methods = [
            "OmniSeqAI automatically detected metadata columns for cell type, condition, sample, patient and batch where available.",
            "Condition analysis used pseudobulk aggregation when sample-level metadata was available.",
            "Differential expression was performed on condition contrasts selected by the comparison policy.",
            "Reported gene tables exclude obvious unmapped taxonomy identifiers where possible.",
            "Exported CSV and JSON tables are saved in the run folder for downstream analysis.",
            "Pathway enrichment is optional and disabled by default to avoid blocking report generation on remote web services.",
        ]
        story.extend(self._bullets(methods))

        return story

    def save(
        self,
        content: Union[Dict[str, Any], str],
        filename: str = "reports/router_report.pdf",
    ) -> None:
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        doc = SimpleDocTemplate(
            str(path),
            pagesize=letter,
            rightMargin=0.65 * inch,
            leftMargin=0.65 * inch,
            topMargin=0.65 * inch,
            bottomMargin=0.65 * inch,
            title="OmniSeqAI Analysis Report",
            author="OmniSeqAI",
        )

        if isinstance(content, dict):
            story = self._build_from_results(content)
        else:
            styles = getSampleStyleSheet()
            story = [Paragraph("OmniSeqAI Analysis Report", styles["Title"])]
            story.append(Paragraph(str(content).replace("\n", "<br/>"), styles["BodyText"]))

        doc.build(
            story,
            onFirstPage=self._page_num,
            onLaterPages=self._page_num,
        )

        print(f"\nRouter PDF saved: {path}")
