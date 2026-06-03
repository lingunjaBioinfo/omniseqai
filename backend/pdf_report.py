from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer
)

from reportlab.lib.styles import (
    getSampleStyleSheet
)


class PDFReport:

    def generate(
        self,
        summary,
        filename="reports/omniseqai_report.pdf"
    ):

        doc = SimpleDocTemplate(
            filename
        )

        styles = (
            getSampleStyleSheet()
        )

        elements = []

        elements.append(
            Paragraph(
                "OMNISEQAI REPORT",
                styles["Title"]
            )
        )

        elements.append(
            Spacer(1, 12)
        )

        for cluster in summary:

            data = summary[cluster]

            elements.append(
                Paragraph(
                    f"Cluster {cluster}",
                    styles["Heading2"]
                )
            )

            elements.append(
                Paragraph(
                    f"Cell Type: "
                    f"{data['cell_type']}",
                    styles["BodyText"]
                )
            )

            elements.append(
                Paragraph(
                    f"Conclusion: "
                    f"{data['conclusion']}",
                    styles["BodyText"]
                )
            )

            elements.append(
                Spacer(1, 10)
            )

        doc.build(elements)

        print(
            f"\nPDF report saved: "
            f"{filename}"
        )
