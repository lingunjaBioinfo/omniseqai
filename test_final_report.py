from backend.final_report import FinalReport

cluster_summary = {
    "0": {
        "cell_type": "T cells",
        "markers": [
            "IL32",
            "LTB",
            "LDHB"
        ],
        "conclusion":
        "Activated T-cell population.",
        "disease_interpretation": [
            "Activated adaptive immune response."
        ]
    }
}

communication = {
    "0": {
        "signals": [
            "Adaptive immune signaling",
            "T-cell activation"
        ]
    }
}

de_findings = [
    "Enhanced T-cell activation detected."
]

reporter = FinalReport()

report = reporter.generate(
    cluster_summary,
    communication,
    de_findings
)

print(report)

reporter.save(report)
