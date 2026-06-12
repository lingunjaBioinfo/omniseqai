from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import scanpy as sc

# Quiet down library noise
sc.settings.verbosity = 0
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*storing '.*' as categorical.*")

try:
    import anndata as ad
    if hasattr(ad, "settings") and hasattr(ad.settings, "verbosity"):
        ad.settings.verbosity = 0
except Exception:
    pass

from backend.pipeline_orchestrator import PipelineOrchestrator


def parse_args():
    parser = argparse.ArgumentParser(description="Run OmniSeqAI.")
    parser.add_argument("-i", "--input", required=True, help="Input .h5ad file")
    parser.add_argument(
        "--mode",
        choices=["auto", "exploratory", "condition"],
        default="auto",
        help="Analysis mode",
    )
    parser.add_argument(
        "--output",
        default="reports/router_report.txt",
        help="Output text report path",
    )
    parser.add_argument(
        "--pdf",
        default="reports/router_report.pdf",
        help="Output PDF report path",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"\nLoading dataset: {input_path}")
    adata = sc.read_h5ad(str(input_path))

    orchestrator = PipelineOrchestrator()
    results = orchestrator.run(
        adata,
        mode=args.mode,
        report_txt=args.output,
        report_pdf=args.pdf,
        generate_pdf=True,
    )

    decision = results.get("decision", None)
    if decision:
        print(f"\nChosen mode: {decision.mode}")
        print(f"Reason: {decision.reason}")

    print("\nDone.")


if __name__ == "__main__":
    main()
