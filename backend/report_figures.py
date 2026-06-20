from __future__ import annotations

from pathlib import Path
import contextlib
import io
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc


class ReportFigures:
    """
    Figure generation utilities for OmniSeqAI reports.

    All methods return the saved figure path as a string, or None if skipped.
    """

    def __init__(self, output_dir: str = "outputs/report_figures"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _path(self, filename: str) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        return self.output_dir / filename

    def _save(self, fig, filename: str):
        path = self._path(filename)
        fig.tight_layout()
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved figure: {path}")
        return str(path)

    def _safe_col(self, adata, col: Optional[str]) -> Optional[str]:
        if adata is None:
            return None

        if col is None:
            return None

        if col in adata.obs.columns:
            return col

        return None

    def _first_existing_col(self, adata, candidates):
        if adata is None:
            return None

        for col in candidates:
            if col in adata.obs.columns:
                return col

        return None

    # ------------------------------------------------------------------
    # UMAP plots
    # ------------------------------------------------------------------
    def umap_with_labels(
        self,
        adata,
        label_col: str,
        title: str = None,
        filename: str = "umap_labels.png",
        point_size: float = 8,
    ):
        """
        Generic UMAP plot colored by any obs column.

        Used for:
        - condition UMAP
        - cell-type UMAP
        - sample UMAP
        - batch UMAP
        - integrated UMAPs
        """

        if adata is None:
            print("No AnnData object available for UMAP plotting.")
            return None

        if "X_umap" not in adata.obsm:
            print("No UMAP found.")
            return None

        if label_col is None or label_col not in adata.obs.columns:
            print(f"UMAP label column not found: {label_col}")
            return None

        output_path = self._path(filename)
        plot_title = title or f"UMAP colored by {label_col}"

        try:
            fig, ax = plt.subplots(figsize=(7.8, 5.8), dpi=150)

            sc.pl.umap(
                adata,
                color=label_col,
                ax=ax,
                show=False,
                size=point_size,
                frameon=False,
                title=plot_title,
            )

            legend = ax.get_legend()

            if legend is not None:
                legend.set_bbox_to_anchor((1.04, 1.0))
                legend._loc = 2

            fig.tight_layout()
            fig.savefig(output_path, bbox_inches="tight")
            plt.close(fig)

            print(f"Saved figure: {output_path}")

            return str(output_path)

        except Exception as e:
            plt.close("all")
            print(f"Failed to generate UMAP for {label_col}: {e}")
            return None

    def umap_by_column(
        self,
        adata,
        column: str = None,
        label_col: str = None,
        title: str = None,
        filename: str = None,
        point_size: float = 8,
        **kwargs,
    ):
        """
        Backward-compatible wrapper for older pipeline_orchestrator calls.

        Accepts:
        - column=
        - label_col=
        - color=
        """

        chosen_col = column or label_col or kwargs.get("color")

        if chosen_col is None:
            print("No column provided for UMAP plot.")
            return None

        if filename is None:
            safe = (
                str(chosen_col)
                .replace(" ", "_")
                .replace("/", "_")
            )
            filename = f"umap_{safe}.png"

        return self.umap_with_labels(
            adata,
            label_col=chosen_col,
            title=title or f"UMAP colored by {chosen_col}",
            filename=filename,
            point_size=point_size,
        )

    def umap_condition_standard(
        self,
        adata,
        condition_col: str = "condition",
        filename: str = "umap_condition_standard.png",
    ):
        return self.umap_with_labels(
            adata,
            label_col=condition_col,
            title="UMAP colored by condition",
            filename=filename,
        )

    def condition_umap(
        self,
        adata,
        condition_col: str = "condition",
        filename: str = "umap_condition_standard.png",
    ):
        return self.umap_condition_standard(
            adata,
            condition_col=condition_col,
            filename=filename,
        )

    def umap_by_condition(
        self,
        adata,
        condition_col: str = "condition",
        filename: str = "umap_condition_standard.png",
    ):
        return self.umap_condition_standard(
            adata,
            condition_col=condition_col,
            filename=filename,
        )

    def umap_celltype_annotated(
        self,
        adata,
        celltype_col: str = "cell_type",
        filename: str = "umap_celltype_annotated.png",
    ):
        return self.umap_with_labels(
            adata,
            label_col=celltype_col,
            title="UMAP annotated by cell type",
            filename=filename,
        )

    def annotated_celltype_umap(
        self,
        adata,
        celltype_col: str = None,
        cell_type_col: str = None,
        label_col: str = None,
        column: str = None,
        title: str = None,
        filename: str = "umap_celltype_annotated.png",
        point_size: float = 8,
        **kwargs,
    ):
        """
        Backward-compatible wrapper for pipeline_orchestrator calls.

        Supports:
        - celltype_col=
        - cell_type_col=
        - label_col=
        - column=
        """

        chosen_col = (
            celltype_col
            or cell_type_col
            or label_col
            or column
            or kwargs.get("color")
            or "cell_type"
        )

        return self.umap_with_labels(
            adata,
            label_col=chosen_col,
            title=title or "UMAP annotated by cell type",
            filename=filename,
            point_size=point_size,
        )

    def celltype_umap(
        self,
        adata,
        celltype_col: str = "cell_type",
        filename: str = "umap_celltype_annotated.png",
    ):
        return self.umap_celltype_annotated(
            adata,
            celltype_col=celltype_col,
            filename=filename,
        )

    def umap_exploratory_labels(
        self,
        adata,
        label_col: str = "cell_type",
        filename: str = "umap_exploratory_labels.png",
    ):
        return self.umap_with_labels(
            adata,
            label_col=label_col,
            title=f"Exploratory UMAP colored by {label_col}",
            filename=filename,
        )

    def umap_exploratory_clusters(
        self,
        adata,
        cluster_col: str = "leiden",
        filename: str = "umap_exploratory_clusters.png",
    ):
        return self.umap_with_labels(
            adata,
            label_col=cluster_col,
            title=f"Exploratory UMAP colored by {cluster_col}",
            filename=filename,
        )

    # ------------------------------------------------------------------
    # Cell-type proportions
    # ------------------------------------------------------------------
    def celltype_proportions_by_condition(
        self,
        adata,
        celltype_col: str = "cell_type",
        condition_col: str = "condition",
        filename: str = "celltype_proportions_by_condition.png",
    ):
        if adata is None:
            print("No AnnData object available for cell-type proportions.")
            return None

        if celltype_col not in adata.obs.columns:
            print(f"Missing cell type column: {celltype_col}")
            return None

        if condition_col not in adata.obs.columns:
            print(f"Missing condition column: {condition_col}")
            return None

        counts = pd.crosstab(
            adata.obs[celltype_col].astype(str),
            adata.obs[condition_col].astype(str),
        )

        if counts.empty:
            print("No cell-type proportion table available.")
            return None

        proportions = counts.div(counts.sum(axis=0), axis=1) * 100.0

        fig_height = max(4.5, 0.35 * proportions.shape[0])
        fig, ax = plt.subplots(figsize=(8.0, fig_height), dpi=150)

        proportions.plot(
            kind="barh",
            stacked=False,
            ax=ax,
        )

        ax.set_xlabel("Percent of cells")
        ax.set_ylabel("Cell type")
        ax.set_title("Cell-type proportions by condition")
        ax.legend(title="Condition", bbox_to_anchor=(1.04, 1), loc="upper left")
        ax.grid(axis="x", alpha=0.25)

        return self._save(fig, filename)

    def celltype_proportions(
        self,
        adata,
        celltype_col: str = "cell_type",
        condition_col: str = "condition",
        filename: str = "celltype_proportions_by_condition.png",
    ):
        return self.celltype_proportions_by_condition(
            adata,
            celltype_col=celltype_col,
            condition_col=condition_col,
            filename=filename,
        )

    # ------------------------------------------------------------------
    # Volcano plot
    # ------------------------------------------------------------------
    def volcano(
        self,
        de_results,
        group1: str,
        group2: str,
        title: str = None,
        filename: str = "volcano_plot.png",
        padj_cutoff: float = 0.05,
        logfc_cutoff: float = 0.5,
        max_labels: int = 18,
    ):
        """
        Robust volcano plot for OmniSeqAI DE results.

        Positive logfoldchanges mean group2/case higher than group1/reference.
        """

        if de_results is None:
            print("No DE results available for volcano plot.")
            return None

        if not hasattr(de_results, "copy"):
            print("DE results object is not a DataFrame-like object.")
            return None

        df = de_results.copy()

        if df.empty:
            print("DE results table is empty; volcano plot skipped.")
            return None

        # --------------------------------------------------
        # Column detection
        # --------------------------------------------------
        gene_col = None
        for col in ["gene_symbol", "names", "gene", "symbol", "gene_name", "gene_id"]:
            if col in df.columns:
                gene_col = col
                break

        logfc_col = None
        for col in ["logfoldchanges", "logFC", "log2FoldChange", "avg_log2FC"]:
            if col in df.columns:
                logfc_col = col
                break

        padj_col = None
        for col in ["pvals_adj", "padj", "p_val_adj", "FDR", "fdr", "qval"]:
            if col in df.columns:
                padj_col = col
                break

        if gene_col is None:
            print(f"Volcano skipped: no gene column found. Columns: {list(df.columns)}")
            return None

        if logfc_col is None:
            print(f"Volcano skipped: no logFC column found. Columns: {list(df.columns)}")
            return None

        if padj_col is None:
            print(f"Volcano skipped: no adjusted p-value column found. Columns: {list(df.columns)}")
            return None

        # --------------------------------------------------
        # Numeric cleanup
        # --------------------------------------------------
        df[gene_col] = df[gene_col].astype(str)
        df[logfc_col] = pd.to_numeric(df[logfc_col], errors="coerce")
        df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")

        df = df.dropna(subset=[gene_col, logfc_col, padj_col]).copy()

        if df.empty:
            print("Volcano skipped: no valid rows after numeric cleanup.")
            return None

        df = df[np.isfinite(df[logfc_col])].copy()
        df = df[np.isfinite(df[padj_col])].copy()

        if df.empty:
            print("Volcano skipped: no finite DE values.")
            return None

        positive_padj = df.loc[df[padj_col] > 0, padj_col]

        if positive_padj.empty:
            min_positive = 1e-300
        else:
            min_positive = max(float(positive_padj.min()), 1e-300)

        df[padj_col] = df[padj_col].clip(lower=min_positive)
        df["neg_log10_padj"] = -np.log10(df[padj_col])

        df["significant"] = (
            (df[padj_col] <= padj_cutoff)
            & (df[logfc_col].abs() >= logfc_cutoff)
        )

        # --------------------------------------------------
        # Label selection
        # --------------------------------------------------
        priority_genes = {
            "ISG15",
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "IFI6",
            "IFI44",
            "IFI44L",
            "MX1",
            "MX2",
            "OAS1",
            "OAS2",
            "OAS3",
            "OASL",
            "IRF7",
            "IRF9",
            "STAT1",
            "STAT2",
            "RSAD2",
            "DDX58",
            "HERC5",
            "HERC6",
            "UBE2L6",
            "XAF1",
            "GBP1",
            "GBP2",
            "IFITM1",
            "IFITM2",
            "IFITM3",
            "CXCL10",
            "S100A8",
            "S100A9",
            "S100A12",
            "HLA-DRA",
            "HLA-DRB1",
            "CD74",
            "B2M",
        }

        df["_gene_upper"] = df[gene_col].astype(str).str.upper()

        label_df = df[
            df["significant"] & df["_gene_upper"].isin(priority_genes)
        ].copy()

        if label_df.empty:
            label_df = df[df["significant"]].copy()

        if not label_df.empty:
            label_df["label_score"] = (
                label_df["neg_log10_padj"]
                * label_df[logfc_col].abs()
            )
            label_df = label_df.sort_values(
                "label_score",
                ascending=False,
            ).head(max_labels)

        # --------------------------------------------------
        # Plot
        # --------------------------------------------------
        output_path = self._path(filename)

        fig, ax = plt.subplots(figsize=(7.8, 6.0), dpi=150)

        nonsig = df[~df["significant"]]
        sig = df[df["significant"]]

        ax.scatter(
            nonsig[logfc_col],
            nonsig["neg_log10_padj"],
            s=8,
            alpha=0.35,
            linewidths=0,
            label="Not significant",
        )

        ax.scatter(
            sig[logfc_col],
            sig["neg_log10_padj"],
            s=10,
            alpha=0.75,
            linewidths=0,
            label="Significant",
        )

        ax.axvline(logfc_cutoff, linestyle="--", linewidth=1)
        ax.axvline(-logfc_cutoff, linestyle="--", linewidth=1)
        ax.axhline(-np.log10(padj_cutoff), linestyle="--", linewidth=1)

        ax.set_xlabel(f"log fold change ({group2} vs {group1})")
        ax.set_ylabel("-log10 adjusted p-value")
        ax.set_title(title or f"{group1} vs {group2}")

        if not label_df.empty:
            texts = []

            for _, row in label_df.iterrows():
                texts.append(
                    ax.text(
                        row[logfc_col],
                        row["neg_log10_padj"],
                        row[gene_col],
                        fontsize=7,
                    )
                )

            try:
                from adjustText import adjust_text

                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    adjust_text(
                        texts,
                        ax=ax,
                        arrowprops=dict(
                            arrowstyle="-",
                            lw=0.4,
                            alpha=0.6,
                        ),
                    )

            except Exception:
                pass

        ax.legend(frameon=False, loc="best")
        ax.grid(alpha=0.2)

        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved figure: {output_path}")

        return str(output_path)

    # ------------------------------------------------------------------
    # Pseudobulk heatmap
    # ------------------------------------------------------------------
    def pseudobulk_heatmap(
        self,
        pb_adata,
        de_results,
        condition_col: str = "condition",
        filename: str = "pseudobulk_heatmap.png",
        top_n: int = 40,
    ):
        if pb_adata is None:
            print("Missing pseudobulk object for heatmap.")
            return None

        if de_results is None or not hasattr(de_results, "empty") or de_results.empty:
            print("Missing DE results for pseudobulk heatmap.")
            return None

        if condition_col not in pb_adata.obs.columns:
            print(f"Missing condition column in pseudobulk object: {condition_col}")
            return None

        df = de_results.copy()

        gene_candidates = []

        for col in ["gene_id", "gene_symbol", "names", "gene"]:
            if col in df.columns:
                gene_candidates.extend(df[col].astype(str).tolist())

        if not gene_candidates:
            print("No gene column found for pseudobulk heatmap.")
            return None

        if "pvals_adj" in df.columns:
            df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")
            df = df.sort_values("pvals_adj", ascending=True)

        elif "padj" in df.columns:
            df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
            df = df.sort_values("padj", ascending=True)

        if "gene_id" in df.columns:
            ranked_genes = df["gene_id"].astype(str).tolist()
        elif "names" in df.columns:
            ranked_genes = df["names"].astype(str).tolist()
        elif "gene_symbol" in df.columns:
            ranked_genes = df["gene_symbol"].astype(str).tolist()
        else:
            ranked_genes = gene_candidates

        valid_genes = [g for g in ranked_genes if g in pb_adata.var_names]

        if len(valid_genes) < 2:
            # Try gene_symbol/name mapping if var has feature_name.
            if "feature_name" in pb_adata.var.columns:
                symbol_to_var = {
                    str(symbol): str(var)
                    for var, symbol in zip(
                        pb_adata.var_names.astype(str),
                        pb_adata.var["feature_name"].astype(str),
                    )
                }

                mapped = []

                for g in ranked_genes:
                    if g in symbol_to_var:
                        mapped.append(symbol_to_var[g])

                valid_genes = [g for g in mapped if g in pb_adata.var_names]

        valid_genes = valid_genes[:top_n]

        if len(valid_genes) < 2:
            print("Not enough valid genes for pseudobulk heatmap.")
            return None

        sub = pb_adata[:, valid_genes].copy()

        X = sub.X

        if hasattr(X, "toarray"):
            X = X.toarray()

        X = np.asarray(X, dtype=float)

        # Samples x genes -> z-score genes across samples
        gene_means = np.nanmean(X, axis=0)
        gene_stds = np.nanstd(X, axis=0)
        gene_stds[gene_stds == 0] = 1.0

        Z = (X - gene_means) / gene_stds

        # Plot genes x samples
        Z_plot = Z.T

        sample_labels = sub.obs_names.astype(str).tolist()
        gene_labels = valid_genes

        conditions = sub.obs[condition_col].astype(str).tolist()
        condition_levels = list(dict.fromkeys(conditions))
        condition_to_code = {c: i for i, c in enumerate(condition_levels)}
        condition_codes = np.array([condition_to_code[c] for c in conditions])

        condition_cmap = plt.get_cmap("tab10", max(len(condition_levels), 1))
        condition_to_color = {
            condition: condition_cmap(i)
            for i, condition in enumerate(condition_levels)
        }

        condition_color_strip = np.array(
            [[condition_to_color[c] for c in conditions]]
        )

        fig_height = max(6.0, 0.18 * len(gene_labels))
        fig_width = max(9.0, 0.28 * len(sample_labels))

        fig = plt.figure(figsize=(fig_width, fig_height), dpi=150)

        gs = fig.add_gridspec(
            nrows=2,
            ncols=2,
            height_ratios=[0.25, 5.0],
            width_ratios=[5.0, 0.25],
            hspace=0.05,
            wspace=0.05,
        )

        ax_top = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[1, 0])
        ax_cbar = fig.add_subplot(gs[1, 1])

        ax_top.imshow(
            condition_color_strip,
            aspect="auto",
            interpolation="nearest",
        )

        ax_top.set_yticks([])
        ax_top.set_xticks([])
        ax_top.set_ylabel("Condition", rotation=0, labelpad=35, va="center")

        im = ax_heat.imshow(
            Z_plot,
            aspect="auto",
            interpolation="nearest",
            cmap="RdBu_r",
            vmin=-2.5,
            vmax=2.5,
        )

        ax_heat.set_yticks(np.arange(len(gene_labels)))
        ax_heat.set_yticklabels(gene_labels, fontsize=7)
        ax_heat.set_xticks([])
        ax_heat.set_xlabel("Pseudobulk samples")
        ax_heat.set_ylabel("Genes")
        fig.suptitle(
            "Top pseudobulk DE genes",
            fontsize=14,
            weight="bold",
            y=0.98,
        )

        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label("Z-score")

        # Legend for condition colors
        handles = []

        for condition in condition_levels:
            handles.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="s",
                    linestyle="",
                    label=condition,
                    markersize=8,
                    markerfacecolor=condition_to_color[condition],
                    markeredgecolor=condition_to_color[condition],
                )
            )

        ax_heat.legend(
            handles=handles,
            title="Condition",
            bbox_to_anchor=(1.22, 1.0),
            loc="upper left",
            frameon=False,
        )

        output_path = self._path(filename)
        fig.subplots_adjust(top=0.90)
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved figure: {output_path}")

        return str(output_path)

    # ------------------------------------------------------------------
    # Biology validation plot
    # ------------------------------------------------------------------
    def biology_signature_hits(
        self,
        biology_validation,
        filename: str = "biology_signature_hits.png",
        top_n: int = 10,
    ):
        """
        Bar plot showing detected biological signature hits.
        """

        if not biology_validation:
            return None

        summary = biology_validation.get("summary")

        if summary is None or getattr(summary, "empty", True):
            print("No biology validation summary available for plotting.")
            return None

        if "n_case_up_hits" not in summary.columns:
            print("Biology validation summary missing n_case_up_hits.")
            return None

        plot_df = summary.copy()

        # Aggregate repeated signature hits across cell types by best hit count.
        agg = (
            plot_df
            .groupby("signature", as_index=False)
            .agg(
                n_case_up_hits=("n_case_up_hits", "max"),
                case_up_hit_fraction=("case_up_hit_fraction", "max"),
            )
        )

        agg = agg.sort_values(
            ["n_case_up_hits", "case_up_hit_fraction"],
            ascending=[False, False],
        ).head(top_n)

        if agg.empty:
            print("No biology signatures available for plotting.")
            return None

        output_path = self._path(filename)

        labels = (
            agg["signature"]
            .astype(str)
            .str.replace("_", " ", regex=False)
        )

        values = agg["n_case_up_hits"].astype(int)

        fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=150)
        ax.barh(labels, values)
        ax.invert_yaxis()
        ax.set_xlabel("Detected expected genes")
        ax.set_title("Biology validation: signature evidence")
        ax.grid(axis="x", alpha=0.25)

        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved figure: {output_path}")

        return str(output_path)
