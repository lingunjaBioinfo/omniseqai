from backend.pdf_report import (
    PDFReport
)

summary = {
    "0": {
        "cell_type": "T cells",
        "conclusion":
        "Activated T-cell population."
    },
    "1": {
        "cell_type": "Monocytes",
        "conclusion":
        "Inflammatory activation."
    }
}

report = PDFReport()

report.generate(summary)
