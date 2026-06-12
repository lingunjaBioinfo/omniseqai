from __future__ import annotations

from typing import Any, Dict

from backend.metadata_detector import MetadataDetector
from backend.mode_selector import ModeSelector
from backend.analysis_router import AnalysisRouter
from backend.router_report import RouterReport
from backend.router_pdf_report import RouterPDFReport
from backend.report_figures import ReportFigures

from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer
from backend.pathway_analysis import PathwayAnalyzer
from backend.conclusion_engine import ConclusionEngine
from backend.cell_communication import CellCommunication


class PipelineOrchestrator:
    """
    Single entry point for OmniSeqAI.

    Modes:
    - exploratory: QC -> clustering -> annotation -> markers -> pathways
    - condition: metadata -> routing -> pseudobulk/DE -> figures -> report
    - auto: choose based on metadata
    """

    def __init__(self):
        self.detector = MetadataDetector()
        self.mode_selector = ModeSelector()
        self.router = AnalysisRouter()
        self.reporter = RouterReport()
        self.pdf_reporter = RouterPDFReport()
        self.figures = ReportFigures()

    def run(
        self,
        adata,
        mode: str = "auto",
        generate_pdf: bool = True,
        report_txt: str = "reports/router_report.txt",
        report_pdf: str = "reports/router_report.pdf",
    ) -> Dict[str, Any]:

        # Detect metadata once for mode selection.
        # Standardization is handled inside the selected branch.
        profile = self.detector.detect(adata)

        decision = self.mode_selector.choose(
            profile,
            adata=adata,
            mode=mode
        )

        results: Dict[str, Any] = {
            "profile": profile,
            "decision": decision,
            "adata": adata,
        }

        # ==================================================
        # CONDITION MODE
        # ==================================================
        if decision.mode == "condition":

            # Router performs standardized metadata mapping once.
            routed = self.router.run(
                adata,
                profile=profile
            )

            results.update(routed)

            figure_paths: Dict[str, Any] = {}

            adata_out = results.get("adata")

            if adata_out is not None:
                figure_paths["umap_condition"] = self.figures.umap_by_column(
                    adata_out,
                    column="condition",
                    filename="umap_condition_standard.png",
                    title="UMAP colored by condition"
                )

                if "cell_type" in adata_out.obs.columns:
                    figure_paths["umap_celltype"] = self.figures.annotated_celltype_umap(
                        adata_out,
                        celltype_col="cell_type"
                    )

                if (
                    "condition" in adata_out.obs.columns
                    and "cell_type" in adata_out.obs.columns
                ):
                    figure_paths["celltype_proportions"] = self.figures.celltype_proportions(
                        adata_out,
                        celltype_col="cell_type",
                        condition_col="condition"
                    )

            condition_de_results = results.get(
                "condition_de_results",
                {}
            )

            for pair, info in condition_de_results.items():
                de_results = info.get("de_results")

                pair_name = (
                    f"{pair[0]}_vs_{pair[1]}"
                    .replace(" ", "_")
                    .replace("/", "_")
                )

                figure_paths["volcano"] = self.figures.volcano(
                    de_results,
                    title=f"{pair[0]} vs {pair[1]}",
                    filename=f"volcano_{pair_name}_standard.png"
                )

                pb = info.get("pseudobulk_adata")
                if pb is not None:
                    figure_paths["pseudobulk_heatmap"] = self.figures.pseudobulk_heatmap(
                        pb,
                        de_results,
                        condition_col="condition",
                        filename=f"pseudobulk_heatmap_{pair_name}.png",
                        top_n=40
                    )

                # Only one primary whole-dataset comparison for now.
                break

            results["figure_paths"] = figure_paths

            report_text = self.reporter.build(results)
            self.reporter.save(
                report_text,
                report_txt
            )

            if generate_pdf:
                self.pdf_reporter.save(
                    results,
                    report_pdf
                )

            results["report_text"] = report_text
            return results

        # ==================================================
        # EXPLORATORY MODE
        # ==================================================
        pipeline = SingleCellPipeline()

        if hasattr(pipeline, "calculate_qc"):
            adata = pipeline.calculate_qc(adata)

        if hasattr(pipeline, "filter_data"):
            adata = pipeline.filter_data(adata)

        if hasattr(pipeline, "normalize_data"):
            adata = pipeline.normalize_data(adata)

        if hasattr(pipeline, "identify_hvg"):
            adata = pipeline.identify_hvg(adata)

        if hasattr(pipeline, "run_pca"):
            adata = pipeline.run_pca(adata)

        if hasattr(pipeline, "compute_neighbors"):
            adata = pipeline.compute_neighbors(adata)

        if hasattr(pipeline, "run_umap"):
            adata = pipeline.run_umap(adata)

        if hasattr(pipeline, "cluster_cells"):
            adata = pipeline.cluster_cells(adata)

        annotator = CellAnnotator()

        try:
            adata = annotator.annotate(adata)
            cluster_map = annotator.cluster_annotations(adata)
        except Exception:
            cluster_map = {}

        marker = MarkerAnalyzer()
        adata = marker.find_markers(adata)

        if cluster_map:
            summary = marker.summarize_clusters(
                adata,
                cluster_map
            )
        else:
            summary = {}

        pathways = PathwayAnalyzer()
        conclusion = ConclusionEngine()
        communication = CellCommunication()

        if summary:
            for cluster in summary:
                genes = summary[cluster].get(
                    "markers",
                    []
                )

                try:
                    p = pathways.top_pathways(genes)
                    summary[cluster]["pathways"] = p
                except Exception:
                    summary[cluster]["pathways"] = None

                summary[cluster]["conclusion"] = conclusion.generate_conclusion(
                    summary[cluster].get("cell_type", "Unknown"),
                    summary[cluster].get("markers", []),
                    ""
                )

            comm = communication.analyze(summary)
        else:
            comm = {}

        results.update(
            {
                "adata": adata,
                "cluster_map": cluster_map,
                "summary": summary,
                "communication": comm,
            }
        )

        figure_paths: Dict[str, Any] = {}

        if "X_umap" in adata.obsm:
            if "cell_type" in adata.obs.columns:
                figure_paths["umap_celltype"] = self.figures.annotated_celltype_umap(
                    adata,
                    celltype_col="cell_type"
                )

            if "condition" in adata.obs.columns:
                figure_paths["umap_condition"] = self.figures.umap_by_column(
                    adata,
                    column="condition",
                    filename="umap_condition_standard.png",
                    title="UMAP colored by condition"
                )

        results["figure_paths"] = figure_paths

        report_lines = []
        report_lines.append(
            "========== OMNISEQAI EXPLORATORY REPORT ==========\n"
        )
        report_lines.append("Mode: exploratory")
        report_lines.append(f"Decision: {decision.reason}")
        report_lines.append(f"Cells: {adata.n_obs:,}")
        report_lines.append(f"Genes: {adata.n_vars:,}\n")

        report_lines.append("Cluster summary:")

        for cluster, info in summary.items():
            report_lines.append("\n" + "=" * 70)
            report_lines.append(f"Cluster {cluster}")
            report_lines.append(
                f"Cell type: {info.get('cell_type')}"
            )
            report_lines.append(
                f"Markers: {', '.join(info.get('markers', [])[:5])}"
            )
            report_lines.append(
                f"Conclusion: {info.get('conclusion')}"
            )
            report_lines.append("Communication:")

            for signal in comm.get(
                cluster,
                {}
            ).get(
                "signals",
                ["Unknown signaling"]
            ):
                report_lines.append(f"- {signal}")

        report_text = "\n".join(report_lines)

        self.reporter.save(
            report_text,
            report_txt
        )

        if generate_pdf:
            self.pdf_reporter.save(
                results,
                report_pdf
            )

        results["report_text"] = report_text

        return results
