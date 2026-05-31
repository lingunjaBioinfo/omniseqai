import os

from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    PageBreak
)

from reportlab.lib.styles import (
    getSampleStyleSheet
)


class ReportGenerator:

    def __init__(self):

        os.makedirs(
            "reports",
            exist_ok=True
        )

    def generate_report(
        self,
        cluster_summaries
    ):

        report = []

        report.append(
            "OMNISEQAI SINGLE-CELL ANALYSIS REPORT\n"
        )

        report.append(
            "=" * 50 + "\n"
        )

        for cluster, info in cluster_summaries.items():

            report.append(
                f"\nCluster {cluster}"
            )

            report.append(
                f"\nCell Type: {info['cell_type']}"
            )

            report.append(
                "\nMarkers:"
            )

            report.append(
                ", ".join(
                    info["markers"][:10]
                )
            )

            report.append(
                "\nInterpretation:"
            )

            report.append(
                info["interpretation"]
            )

            report.append(
                "\n" + "-" * 50
            )

        return "\n".join(report)

    def save_report(
        self,
        report,
        filename="analysis_report.txt"
    ):

        path = f"reports/{filename}"

        with open(
            path,
            "w"
        ) as f:

            f.write(report)

        print(
            f"\nReport saved: {path}"
        )

    def save_pdf_report(
        self,
        cluster_summaries,
        filename="analysis_report.pdf"
    ):

        path = f"reports/{filename}"

        doc = SimpleDocTemplate(path)

        styles = getSampleStyleSheet()

        content = []

        content.append(
            Paragraph(
                "OmniSeqAI Single-Cell Analysis Report",
                styles["Title"]
            )
        )

        content.append(
            Spacer(1, 20)
        )

        content.append(
            Paragraph(
                "Automatically generated analysis report.",
                styles["Normal"]
            )
        )

        content.append(
            PageBreak()
        )

        figures = [
            "outputs/violin_qc_violin.png",
            "outputs/umap_clusters.png",
            "outputs/umap_celltypes.png"
        ]

        for fig in figures:

            content.append(
                Paragraph(
                    fig,
                    styles["Heading2"]
                )
            )

            content.append(
                Image(
                    fig,
                    width=400,
                    height=300
                )
            )

            content.append(
                Spacer(1, 20)
            )

        content.append(
            PageBreak()
        )

        for cluster, info in cluster_summaries.items():

            content.append(
                Paragraph(
                    f"Cluster {cluster}",
                    styles["Heading1"]
                )
            )

            content.append(
                Paragraph(
                    f"Cell Type: {info['cell_type']}",
                    styles["Normal"]
                )
            )

            content.append(
                Paragraph(
                    f"Markers: {', '.join(info['markers'][:10])}",
                    styles["Normal"]
                )
            )

            content.append(
                Paragraph(
                    info["interpretation"],
                    styles["Normal"]
                )
            )

            content.append(
                Spacer(1, 12)
            )

        doc.build(content)

        print(
            f"PDF report saved: {path}"
        )
