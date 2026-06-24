from __future__ import annotations

from typing import Any, Dict
from types import SimpleNamespace

from backend.metadata_detector import MetadataDetector
from backend.mode_selector import ModeSelector
from backend.analysis_router import AnalysisRouter
from backend.gene_symbol_utils import ensure_de_gene_symbols
from backend.router_report import RouterReport
from backend.router_pdf_report import RouterPDFReport
from backend.report_figures import ReportFigures
from backend.biology_validator import BiologyValidator
from backend.biology_report_utils import append_biology_section


from backend.integration_qc import IntegrationQC
from backend.table_exporter import TableExporter
from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer
from backend.pathway_analysis import PathwayAnalyzer
from backend.conclusion_engine import ConclusionEngine
from backend.cell_communication import CellCommunication
from backend.report_table_utils import append_table_section
from backend.integration_report_utils import append_integration_section

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
        self.table_exporter = TableExporter()
        self.integration_qc = IntegrationQC()
        self.biology_validator = BiologyValidator()

    def run(
        self,
        data,
        mode: str = "auto",
        generate_pdf: bool = True,
        report_txt: str = "reports/router_report.txt",
        report_pdf: str = "reports/router_report.pdf",
        user_signature_path=None,
    ) -> Dict[str, Any]:

        adata = data

        # Preserve original gene annotation before router modifies adata.
        # This is required for mapping Ensembl IDs back to gene symbols
        # in volcano plots, heatmaps, and reports.
        gene_symbol_reference = SimpleNamespace(
            var=adata.var.copy(),
            var_names=adata.var_names.astype(str).copy(),
        )

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

        if user_signature_path:
            self.biology_validator = BiologyValidator(
                user_signature_path=user_signature_path
            )
            results["user_signature_path"] = str(user_signature_path)

        integration_qc = self.integration_qc.run(
            adata,
            profile=profile,
            attempt_integration=True,
        )

        results["integration_qc"] = integration_qc

        if integration_qc.get("status") in {"integrated", "qc_only"}:
            print(f"Integration/QC status: {integration_qc.get('status')}")
            print(f"Integration/QC message: {integration_qc.get('message')}")

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
                pb = info.get("pseudobulk_adata")

                symbol_source = gene_symbol_reference

                if (
                    de_results is not None
                    and hasattr(de_results, "empty")
                    and not de_results.empty
                ):
                    de_results = ensure_de_gene_symbols(de_results, symbol_source)
                    info["de_results"] = de_results

                pair_name = (
                    f"{pair[0]}_vs_{pair[1]}"
                    .replace(" ", "_")
                    .replace("/", "_")
                )

                # --------------------------------------------------
                # Volcano plot
                # --------------------------------------------------
                if de_results is None:
                    print(f"Skipping volcano for {pair_name}: de_results is None")

                elif hasattr(de_results, "empty") and de_results.empty:
                    print(f"Skipping volcano for {pair_name}: de_results is empty")

                else:
                    print(
                        f"Generating volcano for {pair_name}: "
                        f"{len(de_results)} genes"
                    )

                    volcano_path = self.figures.volcano(
                        de_results=de_results,
                        group1=pair[0],
                        group2=pair[1],
                        title=f"{pair[0]} vs {pair[1]}",
                        filename=f"volcano_{pair_name}_standard.png",
                    )

                    if volcano_path:
                        figure_paths["volcano"] = volcano_path
                    else:
                        print(f"Volcano function returned None for {pair_name}")

                # --------------------------------------------------
                # Pseudobulk heatmap
                # --------------------------------------------------
                if pb is not None:
                    heatmap_path = self.figures.pseudobulk_heatmap(
                        pb,
                        de_results,
                        condition_col="condition",
                        filename=f"pseudobulk_heatmap_{pair_name}.png",
                        top_n=40,
                    )

                    if heatmap_path:
                        figure_paths["pseudobulk_heatmap"] = heatmap_path

                # Only one primary whole-dataset comparison for now.
                break

            # --------------------------------------------------
            # Integration / batch-aware UMAP figures
            # --------------------------------------------------
            integration_qc = results.get("integration_qc", {}) or {}

            sample_col = integration_qc.get("sample_col")
            batch_col = integration_qc.get("batch_col")
            condition_col = integration_qc.get("condition_col")

            # Original-space UMAPs colored by sample/batch where available.
            if adata_out is not None and "X_umap" in adata_out.obsm:
                if sample_col and sample_col in adata_out.obs.columns:
                    figure_paths["umap_sample"] = self.figures.umap_with_labels(
                        adata_out,
                        label_col=sample_col,
                        title="UMAP colored by sample",
                        filename="umap_sample.png",
                    )

                if batch_col and batch_col in adata_out.obs.columns:
                    figure_paths["umap_batch"] = self.figures.umap_with_labels(
                        adata_out,
                        label_col=batch_col,
                        title="UMAP colored by batch",
                        filename="umap_batch.png",
                    )

            # Integrated-space UMAPs when optional integration succeeds.
            integrated_adata = integration_qc.get("integrated_adata")

            if integrated_adata is not None and "X_umap" in integrated_adata.obsm:
                if sample_col and sample_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_sample"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=sample_col,
                        title="Integrated UMAP colored by sample",
                        filename="umap_integrated_sample.png",
                    )

                if batch_col and batch_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_batch"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=batch_col,
                        title="Integrated UMAP colored by batch",
                        filename="umap_integrated_batch.png",
                    )

                if condition_col and condition_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_condition"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=condition_col,
                        title="Integrated UMAP colored by condition",
                        filename="umap_integrated_condition.png",
                    )

            # --------------------------------------------------
            # Biology validation
            # --------------------------------------------------
            results["biology_validation"] = self.biology_validator.run(results)

            biology_validation = results.get("biology_validation")
            if biology_validation:
                bio_plot = self.figures.biology_signature_hits(
                    biology_validation,
                    filename="biology_signature_hits.png",
                )

                if bio_plot:
                    figure_paths["biology_signature_hits"] = bio_plot

                bio_celltype_plot = self.figures.biology_celltype_signature_heatmap(
                    biology_validation,
                    filename="biology_celltype_signature_heatmap.png",
                )

                if bio_celltype_plot:
                    figure_paths["biology_celltype_signature_heatmap"] = bio_celltype_plot

            results["figure_paths"] = figure_paths
            results["table_paths"] = self.table_exporter.export(results)

            report_text = self.reporter.build(results)
            report_text = append_biology_section(report_text, results)
            report_text = append_integration_section(report_text, results)
            report_text = append_table_section(report_text, results)

            self.reporter.save(
                report_text,
                report_txt,
            )

            if generate_pdf:
                self.pdf_reporter.save(
                    results,
                    report_pdf,
                )

            results["report_text"] = report_text
            return results

        # ==================================================
        # EXPLORATORY MODE
        # ==================================================
        elif decision.mode == "exploratory":

            import numpy as np
            import scanpy as sc

            figure_paths: Dict[str, Any] = {}
            exploratory_log = []

            # --------------------------------------------------
            # Initialize classic exploratory pipeline
            # --------------------------------------------------
            pipeline = SingleCellPipeline()
            pipeline.adata = adata

            def _record(message: str):
                print(message)
                exploratory_log.append(message)

            def _is_anndata(obj) -> bool:
                return hasattr(obj, "obs") and hasattr(obj, "var") and hasattr(obj, "X")

            def _run_pipeline_step(method_name: str, *args, **kwargs):
                """
                Run a SingleCellPipeline method defensively.

                Some older methods modify pipeline.adata in place and return None.
                Some may return AnnData.
                Some may not exist.
                Some may fail on preprocessed datasets.

                Exploratory mode should continue instead of crashing.
                """

                method = getattr(pipeline, method_name, None)

                if method is None:
                    _record(f"Exploratory step skipped: {method_name} not available.")
                    return None

                try:
                    out = method(*args, **kwargs)

                    if _is_anndata(out):
                        pipeline.adata = out

                    _record(f"Exploratory step completed: {method_name}")
                    return out

                except TypeError as e:
                    # Most common bug: passing adata into methods that expect only self.
                    try:
                        out = method()

                        if _is_anndata(out):
                            pipeline.adata = out

                        _record(
                            f"Exploratory step completed after retry without arguments: {method_name}"
                        )
                        return out

                    except Exception as e2:
                        _record(
                            f"Exploratory step failed: {method_name} | {type(e2).__name__}: {e2}"
                        )
                        return None

                except Exception as e:
                    _record(
                        f"Exploratory step failed: {method_name} | {type(e).__name__}: {e}"
                    )
                    return None

            # --------------------------------------------------
            # Run available pipeline steps
            # --------------------------------------------------
            # These are intentionally attempted defensively.
            # If a method does not exist or fails, the fallback Scanpy section below
            # still tries to create useful exploratory outputs.

            for step in [
                "calculate_qc",
                "filter_data",
                "normalize_data",
                "run_pca",
                "compute_neighbors",
                "run_umap",
                "cluster_cells",
            ]:
                _run_pipeline_step(step)

            adata = pipeline.adata

            # --------------------------------------------------
            # Fallback preprocessing if pipeline did not produce UMAP
            # --------------------------------------------------
            # This makes exploratory mode work on new external datasets
            # like Paul15 even if the older SingleCellPipeline branch is incomplete.

            if "X_umap" not in adata.obsm:
                try:
                    _record("Fallback: computing PCA/neighbors/UMAP with Scanpy.")

                    # Avoid modifying raw count data too aggressively if already processed.
                    # Paul15 is already normalized/log-like, so PCA can usually run directly.
                    if "X_pca" not in adata.obsm:
                        n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
                        n_comps = max(2, n_comps)

                        sc.tl.pca(
                            adata,
                            n_comps=n_comps,
                            svd_solver="arpack",
                        )

                    n_pcs = min(30, adata.obsm["X_pca"].shape[1])

                    sc.pp.neighbors(
                        adata,
                        n_neighbors=15,
                        n_pcs=n_pcs,
                    )

                    sc.tl.umap(adata)

                    _record("Fallback UMAP completed.")

                except Exception as e:
                    _record(
                        f"Fallback UMAP failed: {type(e).__name__}: {e}"
                    )

            # --------------------------------------------------
            # Fallback clustering if no cluster column exists
            # --------------------------------------------------
            existing_cluster_cols = [
                col for col in ["leiden", "louvain", "cluster", "clusters"]
                if col in adata.obs.columns
            ]

            if not existing_cluster_cols:
                try:
                    _record("Fallback: computing Leiden clusters.")

                    if "neighbors" not in adata.uns:
                        if "X_pca" not in adata.obsm:
                            n_comps = min(50, adata.n_obs - 1, adata.n_vars - 1)
                            n_comps = max(2, n_comps)

                            sc.tl.pca(
                                adata,
                                n_comps=n_comps,
                                svd_solver="arpack",
                            )

                        n_pcs = min(30, adata.obsm["X_pca"].shape[1])

                        sc.pp.neighbors(
                            adata,
                            n_neighbors=15,
                            n_pcs=n_pcs,
                        )

                    sc.tl.leiden(
                        adata,
                        resolution=0.6,
                        key_added="cluster",
                    )

                    _record("Fallback Leiden clustering completed.")

                except Exception as e:
                    _record(
                        f"Fallback clustering failed: {type(e).__name__}: {e}"
                    )

            # --------------------------------------------------
            # Pick best label column for exploratory plots
            # --------------------------------------------------
            label_col = None

            for candidate in [
                "cell_type",
                "celltype",
                "cell_type_major",
                "annotation",
                "paul15_clusters",
                "bulk_labels",
                "leiden",
                "louvain",
                "cluster",
                "clusters",
            ]:
                if candidate in adata.obs.columns:
                    label_col = candidate
                    break

            if label_col is None:
                _record("No cell-type or cluster label column found for exploratory UMAP.")
            else:
                _record(f"Exploratory label column selected: {label_col}")

            # --------------------------------------------------
            # Marker genes
            # --------------------------------------------------
            if label_col is not None and "rank_genes_groups" not in adata.uns:
                try:
                    # Avoid marker testing if only one group exists.
                    n_groups = adata.obs[label_col].astype(str).nunique()

                    if n_groups > 1:
                        _record(f"Computing marker genes using groupby={label_col}.")

                        sc.tl.rank_genes_groups(
                            adata,
                            groupby=label_col,
                            method="wilcoxon",
                        )

                        _record("Marker gene calculation completed.")
                    else:
                        _record("Marker gene calculation skipped: only one label group.")

                except Exception as e:
                    _record(
                        f"Marker gene calculation failed: {type(e).__name__}: {e}"
                    )

            # --------------------------------------------------
            # Figures
            # --------------------------------------------------
            if "X_umap" in adata.obsm and label_col is not None:
                try:
                    figure_paths["umap_exploratory"] = self.figures.umap_by_column(
                        adata,
                        column=label_col,
                        filename="umap_exploratory_labels.png",
                        title=f"UMAP annotated by {label_col}",
                    )
                except Exception as e:
                    _record(
                        f"Exploratory UMAP figure failed: {type(e).__name__}: {e}"
                    )

            # Also plot clusters separately if cell_type and cluster both exist.
            cluster_col = None

            for candidate in ["leiden", "louvain", "cluster", "clusters"]:
                if candidate in adata.obs.columns:
                    cluster_col = candidate
                    break

            if (
                "X_umap" in adata.obsm
                and cluster_col is not None
                and cluster_col != label_col
            ):
                try:
                    figure_paths["umap_clusters"] = self.figures.umap_by_column(
                        adata,
                        column=cluster_col,
                        filename="umap_exploratory_clusters.png",
                        title=f"UMAP colored by {cluster_col}",
                    )
                except Exception as e:
                    _record(
                        f"Exploratory cluster UMAP failed: {type(e).__name__}: {e}"
                    )

            # --------------------------------------------------
            # Integration / batch-aware UMAP figures
            # --------------------------------------------------
            integration_qc = results.get("integration_qc", {}) or {}

            sample_col = integration_qc.get("sample_col")
            batch_col = integration_qc.get("batch_col")
            condition_col = integration_qc.get("condition_col")

            if "X_umap" in adata.obsm:
                if sample_col and sample_col in adata.obs.columns:
                    figure_paths["umap_sample"] = self.figures.umap_with_labels(
                        adata,
                        label_col=sample_col,
                        title="UMAP colored by sample",
                        filename="umap_sample.png",
                    )

                if batch_col and batch_col in adata.obs.columns:
                    figure_paths["umap_batch"] = self.figures.umap_with_labels(
                        adata,
                        label_col=batch_col,
                        title="UMAP colored by batch",
                        filename="umap_batch.png",
                    )

            integrated_adata = integration_qc.get("integrated_adata")

            if integrated_adata is not None and "X_umap" in integrated_adata.obsm:
                if sample_col and sample_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_sample"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=sample_col,
                        title="Integrated UMAP colored by sample",
                        filename="umap_integrated_sample.png",
                    )

                if batch_col and batch_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_batch"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=batch_col,
                        title="Integrated UMAP colored by batch",
                        filename="umap_integrated_batch.png",
                    )

                if condition_col and condition_col in integrated_adata.obs.columns:
                    figure_paths["umap_integrated_condition"] = self.figures.umap_with_labels(
                        integrated_adata,
                        label_col=condition_col,
                        title="Integrated UMAP colored by condition",
                        filename="umap_integrated_condition.png",
                    )

            # --------------------------------------------------
            # Store results
            # --------------------------------------------------
            results.update(
                {
                    "adata": adata,
                    "figure_paths": figure_paths,
                    "exploratory_results": {
                        "n_cells": int(adata.n_obs),
                        "n_genes": int(adata.n_vars),
                        "label_column": label_col,
                        "cluster_column": cluster_col,
                        "steps": exploratory_log,
                        "has_umap": "X_umap" in adata.obsm,
                        "has_marker_genes": "rank_genes_groups" in adata.uns,
                    },
                }
            )

            results["table_paths"] = self.table_exporter.export(results)

            report_text = self.reporter.build(results)
            report_text = append_integration_section(report_text, results)
            report_text = append_table_section(report_text, results)
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
