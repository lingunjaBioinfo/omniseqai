from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from backend.analysis_router import AnalysisRouter
from backend.router_report import RouterReport
from backend.router_pdf_report import RouterPDFReport


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run OmniSeqAI on a single-cell AnnData file (.h5ad)."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input .h5ad file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="reports/router_report.txt",
        help="Path to output text report file"
    )
    parser.add_argument(
        "--pdf",
        default="reports/router_report.pdf",
        help="Path to output PDF report file"
    )
    parser.add_argument(
        "--preview-lines",
        type=int,
        default=300,
        help="How many lines of the report to print to console"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"\nLoading dataset: {input_path}")
    adata = sc.read_h5ad(str(input_path))

    router = AnalysisRouter()
    profile = router.inspect(adata)

    print("\n===== DETECTED PROFILE =====")
    for k, v in profile.items():
        print(f"{k}: {v}")

    results = router.run(adata, profile=profile)

    reporter = RouterReport()
    report = reporter.build(results)

    print("\n" + report[: args.preview_lines * 80])

    reporter.save(report, args.output)
    RouterPDFReport().save(report, args.pdf)

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
