from __future__ import annotations

from pathlib import Path
from typing import Optional, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.cluster.hierarchy as sch
from matplotlib.patches import Patch


class ReportFigures:
    """
    Figure generator used by PipelineOrchestrator.

    Required methods:
    - umap_by_column
    - annotated_celltype_umap
    - celltype_proportions
    - volcano
    - pseudobulk_heatmap
    """

    def __init__(self, output_dir: str = "outputs/report_figures"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _save(self, fig, filename: str) -> str:
        path = self.output_dir / filename
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved figure: {path}")
        return str(path)

    def _clean_label(self, label: Any) -> str:
        label = str(label)

        replacements = {
            "MO.Classical": "Classical monocyte",
            "MO.Non-classical": "Non-classical monocyte",
            "MO.Classical_NEAT1hi": "Intermediate monocyte",
            "MO.NKG7": "NKG7+ monocyte",
            "DC.Plasmacytoid": "Plasmacytoid DC",
            "Dendritic": "Dendritic cell",
            "normal": "Healthy",
            "COVID-19": "COVID",
            "cv_19": "COVID",
            "HC": "Healthy",
            "Vax": "Vaccinated",
            "booster": "Vaccinated",
        }

        return replacements.get(label, label)

    def _palette(self, categories):
        cmap_names = ["tab20", "tab20b", "tab20c", "Set3", "Dark2"]
        colors = []

        for cmap_name in cmap_names:
            cmap = plt.get_cmap(cmap_name)
            n = getattr(cmap, "N", 20)
            colors.extend([cmap(i / max(n - 1, 1)) for i in range(n)])

        return {cat: colors[i % len(colors)] for i, cat in enumerate(categories)}

    def umap_by_column(
        self,
        adata,
        column: str,
        filename: Optional[str] = None,
        title: Optional[str] = None,
        max_categories: int = 12,
    ):
        if "X_umap" not in adata.obsm:
            print("No UMAP found.")
            return None

        if column not in adata.obs.columns:
            print(f"{column} not found in obs.")
            return None

        filename = filename or f"umap_{column}.png"
        title = title or f"UMAP colored by {column}"

        umap = adata.obsm["X_umap"]
        labels = adata.obs[column].astype(str).map(self._clean_label)

        counts = labels.value_counts()
        categories = counts.index[:max_categories].tolist()
        palette = self._palette(categories)

        fig, ax = plt.subplots(figsize=(8.8, 6.8), dpi=300)

        minor_mask = ~labels.isin(categories)
        if minor_mask.sum() > 0:
            ax.scatter(
                umap[minor_mask, 0],
                umap[minor_mask, 1],
                s=3,
                alpha=0.12,
                color="#BDBDBD",
                linewidths=0,
                rasterized=True,
                label=f"Other (n={int(minor_mask.sum()):,})",
            )

        for cat in categories:
            mask = labels == cat
            ax.scatter(
                umap[mask, 0],
                umap[mask, 1],
                s=4,
                alpha=0.75,
                color=palette[cat],
                linewidths=0,
                rasterized=True,
                label=f"{cat} (n={int(mask.sum()):,})",
            )

        ax.set_title(title, fontsize=16, pad=12, weight="bold")
        ax.set_xlabel("UMAP1", fontsize=11)
        ax.set_ylabel("UMAP2", fontsize=11)

        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            fontsize=8,
            markerscale=3,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        fig.tight_layout()
        return self._save(fig, filename)

    def annotated_celltype_umap(
        self,
        adata,
        celltype_col: str = "cell_type",
        min_fraction_for_label: float = 0.015,
        max_categories: int = 14,
    ):
        if "X_umap" not in adata.obsm:
            print("No UMAP found.")
            return None

        if celltype_col not in adata.obs.columns:
            print(f"{celltype_col} not found.")
            return None

        umap = adata.obsm["X_umap"]
        labels = adata.obs[celltype_col].astype(str).map(self._clean_label)

        counts = labels.value_counts()
        categories = counts.index[:max_categories].tolist()
        palette = self._palette(categories)

        fig, ax = plt.subplots(figsize=(9.2, 7.0), dpi=300)

        minor_mask = ~labels.isin(categories)
        if minor_mask.sum() > 0:
            ax.scatter(
                umap[minor_mask, 0],
                umap[minor_mask, 1],
                s=3,
                alpha=0.10,
                color="#BDBDBD",
                linewidths=0,
                rasterized=True,
                label=f"Other (n={int(minor_mask.sum()):,})",
            )

        for cat in categories:
            mask = labels == cat
            ax.scatter(
                umap[mask, 0],
                umap[mask, 1],
                s=4,
                alpha=0.72,
                color=palette[cat],
                linewidths=0,
                rasterized=True,
                label=f"{cat} (n={int(mask.sum()):,})",
            )

        total = len(labels)

        for cat in categories:
            mask = labels == cat
            frac = mask.sum() / total

            if frac < min_fraction_for_label:
                continue

            x = float(np.median(umap[mask, 0]))
            y = float(np.median(umap[mask, 1]))

            ax.text(
                x,
                y,
                cat,
                fontsize=8,
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="round,pad=0.25",
                    facecolor="white",
                    edgecolor="lightgrey",
                    alpha=0.85,
                ),
            )

        ax.set_title("UMAP annotated by cell type", fontsize=16, pad=12, weight="bold")
        ax.set_xlabel("UMAP1", fontsize=11)
        ax.set_ylabel("UMAP2", fontsize=11)

        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            fontsize=8,
            markerscale=3,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        fig.tight_layout()
        return self._save(fig, "umap_celltype_annotated.png")

    def celltype_proportions(
        self,
        adata,
        celltype_col: str = "cell_type",
        condition_col: str = "condition",
    ):
        if celltype_col not in adata.obs.columns or condition_col not in adata.obs.columns:
            print("Missing cell type or condition column.")
            return None

        obs = adata.obs[[celltype_col, condition_col]].copy()
        obs[celltype_col] = obs[celltype_col].astype(str).map(self._clean_label)
        obs[condition_col] = obs[condition_col].astype(str).map(self._clean_label)

        tab = pd.crosstab(
            obs[celltype_col],
            obs[condition_col],
            normalize="columns",
        )

        if tab.empty:
            print("No cell-type proportion table could be generated.")
            return None

        if tab.shape[1] >= 2:
            first, second = tab.columns[:2]
            tab["_delta"] = (tab[first] - tab[second]).abs()
            tab = tab.sort_values("_delta", ascending=False).drop(columns="_delta")
        else:
            tab = tab.sort_index()

        fig, ax = plt.subplots(figsize=(9.8, 5.8), dpi=300)
        tab.plot(kind="bar", ax=ax, width=0.82)

        ax.set_title("Cell-type proportions by condition", fontsize=16, pad=12, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Proportion of cells", fontsize=11)
        ax.legend(title="Condition", frameon=False, fontsize=9)

        ax.tick_params(axis="x", labelrotation=35)
        for label in ax.get_xticklabels():
            label.set_ha("right")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.22)

        fig.tight_layout()
        return self._save(fig, "celltype_proportions_by_condition.png")

    def volcano(
        self,
        de_results: pd.DataFrame,
        title: str = "Volcano plot",
        filename: str = "volcano_plot_standard.png",
        sig_p: float = 0.05,
        sig_lfc: float = 0.5,
        annotate_top_n: int = 8,
    ):
        try:
            from adjustText import adjust_text
            has_adjust_text = True
        except Exception:
            has_adjust_text = False

        if de_results is None or de_results.empty:
            print("No DE results for volcano.")
            return None

        df = de_results.copy()

        required = {"names", "logfoldchanges", "pvals_adj"}
        if not required.issubset(df.columns):
            print("DE results missing required columns.")
            return None

        df["names"] = df["names"].astype(str).str.strip()
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

        df = df.dropna(subset=["names", "logfoldchanges", "pvals_adj"])

        df = df[
            ~df["names"].str.startswith("ENSG", na=False)
            & ~df["names"].str.startswith("NCBITaxon:", na=False)
        ].copy()

        if df.empty:
            print("No symbol-mapped genes for volcano.")
            return None

        df["pvals_adj"] = df["pvals_adj"].clip(lower=1e-300, upper=1.0)
        df["neglog10_padj"] = -np.log10(df["pvals_adj"])

        df["direction"] = "Not significant"

        df.loc[
            (df["pvals_adj"] < sig_p) & (df["logfoldchanges"] >= sig_lfc),
            "direction",
        ] = "Upregulated"

        df.loc[
            (df["pvals_adj"] < sig_p) & (df["logfoldchanges"] <= -sig_lfc),
            "direction",
        ] = "Downregulated"

        up_n = int((df["direction"] == "Upregulated").sum())
        down_n = int((df["direction"] == "Downregulated").sum())

        fig, ax = plt.subplots(figsize=(9.0, 6.8), dpi=300)

        colors = {
            "Not significant": "#BDBDBD",
            "Upregulated": "#D62728",
            "Downregulated": "#1F77B4",
        }

        for group in ["Not significant", "Downregulated", "Upregulated"]:
            sub = df[df["direction"] == group]

            if sub.empty:
                continue

            ax.scatter(
                sub["logfoldchanges"],
                sub["neglog10_padj"],
                s=15 if group == "Not significant" else 22,
                alpha=0.30 if group == "Not significant" else 0.85,
                color=colors[group],
                linewidths=0,
                rasterized=True,
                label=f"{group} (n={len(sub):,})",
            )

        ax.axhline(-np.log10(sig_p), linestyle="--", linewidth=1, color="black", alpha=0.55)
        ax.axvline(sig_lfc, linestyle="--", linewidth=1, color="black", alpha=0.55)
        ax.axvline(-sig_lfc, linestyle="--", linewidth=1, color="black", alpha=0.55)

        sig = df[df["direction"] != "Not significant"].copy()

        if not sig.empty:
            label_df = sig.sort_values(
                ["pvals_adj", "logfoldchanges"],
                ascending=[True, False],
                kind="mergesort",
            ).head(annotate_top_n)
        else:
            label_df = df.sort_values(
                ["pvals_adj", "logfoldchanges"],
                ascending=[True, False],
                kind="mergesort",
            ).head(annotate_top_n)

        texts = []
        for _, row in label_df.iterrows():
            texts.append(
                ax.text(
                    row["logfoldchanges"],
                    row["neglog10_padj"],
                    row["names"],
                    fontsize=8,
                    weight="bold",
                )
            )

        if has_adjust_text and texts:
            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(
                    arrowstyle="-",
                    color="black",
                    lw=0.5,
                    alpha=0.6,
                ),
                expand_points=(1.2, 1.3),
                expand_text=(1.2, 1.3),
            )

        ax.set_title(title, fontsize=17, pad=14, weight="bold")
        ax.set_xlabel("Log fold change", fontsize=12)
        ax.set_ylabel("-log10 adjusted p-value", fontsize=12)

        ax.text(
            0.02,
            0.98,
            (
                f"Upregulated: {up_n:,}\n"
                f"Downregulated: {down_n:,}\n"
                f"FDR < {sig_p}, |logFC| ≥ {sig_lfc}"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.35",
                fc="white",
                ec="lightgrey",
                alpha=0.9,
            ),
        )

        ax.legend(frameon=False, fontsize=9, loc="center right")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(alpha=0.18)

        fig.tight_layout()
        return self._save(fig, filename)

    def pseudobulk_heatmap(
        self,
        pb_adata,
        de_results: pd.DataFrame,
        condition_col: str = "condition",
        filename: str = "pseudobulk_heatmap_top_genes.png",
        top_n: int = 40,
    ):
        if pb_adata is None or de_results is None or de_results.empty:
            print("Missing pseudobulk object or DE results.")
            return None

        if condition_col not in pb_adata.obs.columns:
            print(f"{condition_col} not found in pseudobulk obs.")
            return None

        df = de_results.copy()
        df["names"] = df["names"].astype(str).str.strip()
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")

        genes = (
            df[
                (df["pvals_adj"] < 0.05)
                & ~df["names"].str.startswith("ENSG", na=False)
                & ~df["names"].str.startswith("NCBITaxon:", na=False)
            ]
            .sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False])
            ["names"]
            .tolist()
        )

        genes = [g for g in genes if g in pb_adata.var_names][:top_n]

        if len(genes) < 3:
            print("Not enough valid genes for pseudobulk heatmap.")
            return None

        X = pb_adata[:, genes].X
        if sp.issparse(X):
            X = X.toarray()

        mat = pd.DataFrame(
            X,
            index=pb_adata.obs_names,
            columns=genes,
        )

        z = (mat - mat.mean(axis=0)) / mat.std(axis=0).replace(0, 1)
        z = z.fillna(0)
        z = z.clip(-2.5, 2.5)

        gene_linkage = sch.linkage(z.T, method="average", metric="euclidean")
        gene_order = sch.leaves_list(gene_linkage)
        z = z.iloc[:, gene_order]

        conditions = pb_adata.obs.loc[z.index, condition_col].astype(str).map(self._clean_label)

        ordered_samples = []

        for cond in conditions.value_counts().index.tolist():
            idx = conditions[conditions == cond].index.tolist()

            if len(idx) > 2:
                sub = z.loc[idx]
                sample_linkage = sch.linkage(sub, method="average", metric="euclidean")
                ordered = [idx[i] for i in sch.leaves_list(sample_linkage)]
            else:
                ordered = idx

            ordered_samples.extend(ordered)

        z = z.loc[ordered_samples]
        conditions = conditions.loc[ordered_samples]

        condition_categories = conditions.unique().tolist()
        palette = self._palette(condition_categories)

        fig = plt.figure(figsize=(12.0, 9.0), dpi=300)

        gs = fig.add_gridspec(
            nrows=3,
            ncols=2,
            height_ratios=[0.35, 0.25, 7],
            width_ratios=[1.25, 8],
            hspace=0.03,
            wspace=0.03,
        )

        ax_title = fig.add_subplot(gs[0, :])
        ax_colbar = fig.add_subplot(gs[1, 1])
        ax_gene_dendro = fig.add_subplot(gs[2, 0])
        ax_heat = fig.add_subplot(gs[2, 1])

        ax_title.axis("off")
        ax_title.text(
            0.5,
            0.5,
            "Top pseudobulk DE genes",
            ha="center",
            va="center",
            fontsize=17,
            weight="bold",
        )

        sch.dendrogram(
            gene_linkage,
            orientation="left",
            ax=ax_gene_dendro,
            no_labels=True,
            color_threshold=None,
        )

        ax_gene_dendro.axis("off")

        cond_colors = np.array([palette[c] for c in conditions])
        ax_colbar.imshow(cond_colors[np.newaxis, :, :], aspect="auto")
        ax_colbar.set_xticks([])
        ax_colbar.set_yticks([])
        ax_colbar.set_ylabel("Condition", fontsize=9, rotation=0, labelpad=35)

        im = ax_heat.imshow(
            z.T.values,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-2.5,
            vmax=2.5,
            interpolation="nearest",
        )

        ax_heat.set_yticks(range(len(z.columns)))
        ax_heat.set_yticklabels(z.columns, fontsize=7)
        ax_heat.set_xticks([])
        ax_heat.set_xlabel("Pseudobulk samples", fontsize=11)
        ax_heat.set_ylabel("Genes", fontsize=11)

        cond_values = conditions.values
        changes = np.where(cond_values[:-1] != cond_values[1:])[0]

        for c in changes:
            ax_heat.axvline(c + 0.5, color="black", linewidth=1.1)
            ax_colbar.axvline(c + 0.5, color="black", linewidth=1.1)

        cbar = fig.colorbar(im, ax=ax_heat, shrink=0.75, pad=0.015)
        cbar.set_label("Z-score", fontsize=10)

        handles = [
            Patch(facecolor=palette[c], edgecolor="none", label=c)
            for c in condition_categories
        ]

        ax_heat.legend(
            handles=handles,
            title="Condition",
            bbox_to_anchor=(1.18, 1.0),
            loc="upper left",
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )

        return self._save(fig, filename)
