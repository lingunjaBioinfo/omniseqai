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

from backend.table_exporter import TableExporter
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
        self.table_exporter = TableExporter()

    def run(
        self,
        adata,
        mode: str = "auto",
        generate_pdf: bool = True,
        report_txt: str = "reports/router_report.txt",
        report_pdf: str = "reports/router_report.pdf",
    ) -> Dict[str, Any]:

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

                figure_paths["volcano"] = self.figures.volcano(
                    de_results=de_results,
                    group1=pair[0],
                    group2=pair[1],
                    title=f"{pair[0]} vs {pair[1]}",
                    filename=f"volcano_{pair_name}_standard.png",
                )

                if pb is not None:
                    figure_paths["pseudobulk_heatmap"] = self.figures.pseudobulk_heatmap(
                        pb,
                        de_results,
                        condition_col="condition",
                        filename=f"pseudobulk_heatmap_{pair_name}.png",
                        top_n=40,
                    )

                # Only one primary whole-dataset comparison for now.
                break

            results["figure_paths"] = figure_paths
            results["table_paths"] = self.table_exporter.export(results)
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
