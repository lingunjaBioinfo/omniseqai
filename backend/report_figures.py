from __future__ import annotations

from pathlib import Path
from typing import Optional, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.cluster.hierarchy as sch
from matplotlib.patches import Patch
from backend.signature_genes import priority_label_genes

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
        """
        Polished UMAP plot for condition/cell-type metadata columns.

        Features:
        - readable title and axes
        - minor categories grouped as Other
        - legend outside plot
        - small point size for large datasets
        - high-resolution output
        """

        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        if "X_umap" not in adata.obsm:
            print("No UMAP found.")
            return None

        if column not in adata.obs.columns:
            print(f"{column} not found in obs.")
            return None

        filename = filename or f"umap_{column}.png"
        title = title or f"UMAP colored by {column}"

        output_dir = (
            getattr(self, "output_dir", None)
            or getattr(self, "outdir", None)
            or getattr(self, "figures_dir", None)
            or "outputs/report_figures"
        )

        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, filename)

        umap = adata.obsm["X_umap"]

        clean_label = getattr(self, "_clean_label", lambda x: str(x))
        labels = adata.obs[column].astype(str).map(clean_label)

        counts = labels.value_counts()
        categories = counts.index[:max_categories].tolist()

        if hasattr(self, "_palette"):
            palette = self._palette(categories)
        else:
            cmap = plt.get_cmap("tab20")
            palette = {
                cat: cmap(i / max(len(categories) - 1, 1))
                for i, cat in enumerate(categories)
            }

        n_cells = adata.n_obs

        if n_cells > 100_000:
            point_size = 1.2
            alpha = 0.45
        elif n_cells > 50_000:
            point_size = 1.8
            alpha = 0.55
        elif n_cells > 20_000:
            point_size = 2.5
            alpha = 0.65
        else:
            point_size = 4.0
            alpha = 0.75

        fig, ax = plt.subplots(figsize=(9.2, 7.0), dpi=300)

        # Plot minor categories first as background.
        minor_mask = ~labels.isin(categories)

        if minor_mask.sum() > 0:
            ax.scatter(
                umap[minor_mask, 0],
                umap[minor_mask, 1],
                s=max(point_size * 0.8, 0.8),
                alpha=0.12,
                color="#BDBDBD",
                linewidths=0,
                rasterized=True,
                label=f"Other (n={int(minor_mask.sum()):,})",
            )

        # Plot major categories.
        for cat in categories:
            mask = labels == cat

            ax.scatter(
                umap[mask, 0],
                umap[mask, 1],
                s=point_size,
                alpha=alpha,
                color=palette.get(cat, "#333333"),
                linewidths=0,
                rasterized=True,
                label=f"{cat} (n={int(mask.sum()):,})",
            )

        ax.set_title(title, fontsize=16, fontweight="bold", pad=14)
        ax.set_xlabel("UMAP 1", fontsize=12)
        ax.set_ylabel("UMAP 2", fontsize=12)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        ax.tick_params(axis="both", labelsize=10)

        # Put legend outside so it does not cover cells.
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=9,
            markerscale=2.0,
            borderaxespad=0.0,
        )

        fig.tight_layout()

        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved figure: {output_path}")
        return output_path

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
        de_results=None,
        group1=None,
        group2=None,
        output_path=None,
        title=None,
        padj_cutoff: float = 0.05,
        lfc_cutoff: float = 0.5,
        max_labels: int = 12,
        **kwargs,
    ):
        """
        Generate a publication-readable volcano plot.

        Fixes:
        - parses group labels from title when needed
        - positive logFC direction is explicit
        - caps extreme x/y display values without changing statistics
        - prioritizes biologically meaningful signature genes for labels
        - avoids labeling Ensembl IDs unless no better labels exist

        Convention:
            group1 = baseline/reference
            group2 = case/test
            positive logFC = group2 higher than group1
        """

        import os
        import re
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        try:
            from backend.signature_genes import DEFAULT_PRIORITY_GENES, normalize_gene_name
        except Exception:
            DEFAULT_PRIORITY_GENES = [
                "ISG15",
                "ISG20",
                "IFIT1",
                "IFIT2",
                "IFIT3",
                "IFI6",
                "IFI44L",
                "MX1",
                "OAS1",
                "OAS2",
                "STAT1",
                "IRF7",
                "CXCL10",
                "DDX58",
                "HERC6",
                "UBE2L6",
                "BST2",
                "RSAD2",
                "EPSTI1",
                "GBP1",
                "GBP5",
                "CCL2",
                "CCL3",
                "CCL4",
                "CCL8",
                "IL1B",
                "TNF",
                "NFKBIA",
                "HLA-A",
                "HLA-B",
                "HLA-C",
                "HLA-E",
                "TAP1",
                "TAP2",
            ]

            def normalize_gene_name(x):
                return str(x).strip().upper()

        # --------------------------------------------------
        # Resolve DE table
        # --------------------------------------------------
        if de_results is None:
            de_results = (
                kwargs.get("df")
                or kwargs.get("data")
                or kwargs.get("de_table")
                or kwargs.get("de")
            )

        # --------------------------------------------------
        # Resolve title and groups
        # --------------------------------------------------
        if title is None:
            title = kwargs.get("plot_title") or kwargs.get("comparison_title")

        if group1 is None:
            group1 = (
                kwargs.get("group1")
                or kwargs.get("reference")
                or kwargs.get("ref")
                or kwargs.get("baseline")
                or kwargs.get("control")
            )

        if group2 is None:
            group2 = (
                kwargs.get("group2")
                or kwargs.get("case")
                or kwargs.get("test")
                or kwargs.get("condition")
                or kwargs.get("treatment")
            )

        comparison = kwargs.get("comparison")

        if comparison is not None and (group1 is None or group2 is None):
            try:
                group1, group2 = comparison
            except Exception:
                pass

        if (group1 is None or group2 is None) and title is not None:
            title_str = str(title).strip()

            match = re.match(
                r"^\s*(.*?)\s+(?:vs|versus)\s+(.*?)\s*$",
                title_str,
                flags=re.IGNORECASE,
            )

            if match:
                if group1 is None:
                    group1 = match.group(1).strip()

                if group2 is None:
                    group2 = match.group(2).strip()

            elif "_vs_" in title_str:
                parts = title_str.split("_vs_", 1)

                if len(parts) == 2:
                    if group1 is None:
                        group1 = parts[0].strip()

                    if group2 is None:
                        group2 = parts[1].strip()

        group1 = str(group1) if group1 is not None else "Reference"
        group2 = str(group2) if group2 is not None else "Case"

        if title is None:
            title = f"{group1} vs {group2}"

        # --------------------------------------------------
        # Resolve output path
        # --------------------------------------------------
        if output_path is None:
            output_path = (
                kwargs.get("save_path")
                or kwargs.get("fig_path")
                or kwargs.get("path")
                or kwargs.get("outfile")
            )

        filename = kwargs.get("filename")

        output_dir = (
            getattr(self, "output_dir", None)
            or getattr(self, "outdir", None)
            or getattr(self, "figures_dir", None)
            or "outputs/report_figures"
        )

        os.makedirs(output_dir, exist_ok=True)

        if output_path is None and filename is not None:
            output_path = os.path.join(output_dir, filename)

        if output_path is None:
            output_path = os.path.join(
                output_dir,
                f"volcano_{group1}_vs_{group2}_standard.png",
            )

        # --------------------------------------------------
        # Validate DE table
        # --------------------------------------------------
        if de_results is None:
            print("No DE results for volcano.")
            return None

        if not hasattr(de_results, "empty") or de_results.empty:
            print("No DE results for volcano.")
            return None

        required = {"names", "logfoldchanges", "pvals_adj"}

        if not required.issubset(set(de_results.columns)):
            print("DE results missing required columns for volcano.")
            return None

        df = de_results.copy()

        df["names"] = df["names"].astype(str)
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

        df = df.dropna(subset=["names", "logfoldchanges", "pvals_adj"]).copy()

        if df.empty:
            print("No usable DE results for volcano.")
            return None

        # --------------------------------------------------
        # Better display labels
        # --------------------------------------------------
        symbol_candidates = [
            "gene_symbol",
            "symbol",
            "gene",
            "features",
            "feature_name",
        ]

        display_col = None

        for col in symbol_candidates:
            if col in df.columns:
                display_col = col
                break

        if display_col is not None:
            df["display_gene"] = df[display_col].astype(str)
            bad_symbol = (
                df["display_gene"].isna()
                | (df["display_gene"].astype(str).str.strip() == "")
                | (df["display_gene"].astype(str).str.lower() == "nan")
            )
            df.loc[bad_symbol, "display_gene"] = df.loc[bad_symbol, "names"].astype(str)
        else:
            df["display_gene"] = df["names"].astype(str)

        df["display_gene"] = df["display_gene"].astype(str)
        df["display_upper"] = df["display_gene"].map(normalize_gene_name)

        def is_ensembl_like(value):
            value = str(value).strip().upper()
            return value.startswith("ENSG") or value.startswith("ENSMUSG")

        df["is_ensembl_label"] = df["display_gene"].map(is_ensembl_like)

        # --------------------------------------------------
        # Compute transformed axes
        # --------------------------------------------------
        tiny = np.nextafter(0, 1)
        df["pvals_adj_safe"] = df["pvals_adj"].clip(lower=tiny)
        df["neg_log10_padj_raw"] = -np.log10(df["pvals_adj_safe"])

        # Display caps prevent one or two extreme values from ruining the plot.
        abs_lfc = df["logfoldchanges"].abs().replace([np.inf, -np.inf], np.nan).dropna()
        y_raw = df["neg_log10_padj_raw"].replace([np.inf, -np.inf], np.nan).dropna()

        if len(abs_lfc) > 0:
            x_cap = float(np.nanpercentile(abs_lfc, 99.0))
            x_cap = max(3.0, min(x_cap, 8.0))
        else:
            x_cap = 5.0

        if len(y_raw) > 0:
            y_cap = float(np.nanpercentile(y_raw, 99.0))
            y_cap = max(10.0, min(y_cap, 60.0))
        else:
            y_cap = 20.0

        df["x_plot"] = df["logfoldchanges"].clip(lower=-x_cap, upper=x_cap)
        df["y_plot"] = df["neg_log10_padj_raw"].clip(upper=y_cap)

        # --------------------------------------------------
        # Categorize genes
        # --------------------------------------------------
        higher_case_label = f"Higher in {group2}"
        higher_ref_label = f"Higher in {group1}"

        df["category"] = "Not significant"

        df.loc[
            (df["pvals_adj"] < padj_cutoff)
            & (df["logfoldchanges"] >= lfc_cutoff),
            "category",
        ] = higher_case_label

        df.loc[
            (df["pvals_adj"] < padj_cutoff)
            & (df["logfoldchanges"] <= -lfc_cutoff),
            "category",
        ] = higher_ref_label

        up_count = int((df["category"] == higher_case_label).sum())
        down_count = int((df["category"] == higher_ref_label).sum())
        nonsig_count = int((df["category"] == "Not significant").sum())

        # --------------------------------------------------
        # Select labels
        # --------------------------------------------------
        df["abs_logfoldchanges"] = df["logfoldchanges"].abs()

        significant = df[
            (df["pvals_adj"] < padj_cutoff)
            & (df["logfoldchanges"].abs() >= lfc_cutoff)
        ].copy()

        selected_indices = []
        selected_upper = set()

        priority_upper = [normalize_gene_name(g) for g in DEFAULT_PRIORITY_GENES]

        # First: known signature genes.
        for gene_upper in priority_upper:
            match = significant[significant["display_upper"] == gene_upper].copy()

            if match.empty:
                continue

            match = match.sort_values(
                ["pvals_adj", "abs_logfoldchanges"],
                ascending=[True, False],
            )

            idx = match.index[0]
            label = str(df.loc[idx, "display_gene"])
            label_upper = normalize_gene_name(label)

            if label_upper not in selected_upper:
                selected_indices.append(idx)
                selected_upper.add(label_upper)

            if len(selected_indices) >= max_labels:
                break

        # Second: strongest non-Ensembl genes.
        if len(selected_indices) < max_labels:
            fallback = significant[
                ~significant["display_upper"].isin(selected_upper)
                & (~significant["is_ensembl_label"])
            ].copy()

            fallback = fallback.sort_values(
                ["pvals_adj", "abs_logfoldchanges"],
                ascending=[True, False],
            )

            for idx, row in fallback.iterrows():
                label = str(row["display_gene"])
                label_upper = normalize_gene_name(label)

                if label_upper not in selected_upper:
                    selected_indices.append(idx)
                    selected_upper.add(label_upper)

                if len(selected_indices) >= max_labels:
                    break

        if selected_indices:
            label_df = df.loc[selected_indices].copy()
            label_df = label_df.drop_duplicates(subset=["display_upper"])
        else:
            label_df = pd.DataFrame(columns=df.columns)

        

        # --------------------------------------------------
        # Plot
        # --------------------------------------------------
        fig, ax = plt.subplots(figsize=(9.8, 7.2))

        colors = {
            "Not significant": "#C7C7C7",
            higher_ref_label: "#377EB8",
            higher_case_label: "#E41A1C",
        }

        plot_order = [
            "Not significant",
            higher_ref_label,
            higher_case_label,
        ]

        for category in plot_order:
            sub = df[df["category"] == category]

            if sub.empty:
                continue

            if category == higher_ref_label:
                label = f"{higher_ref_label} (n={down_count})"
            elif category == higher_case_label:
                label = f"{higher_case_label} (n={up_count})"
            else:
                label = f"Not significant (n={nonsig_count})"

            ax.scatter(
                sub["x_plot"],
                sub["y_plot"],
                s=13,
                alpha=0.72,
                c=colors.get(category, "#C7C7C7"),
                edgecolors="none",
                label=label,
            )

        ax.axvline(
            lfc_cutoff,
            linestyle="--",
            linewidth=1.0,
            color="black",
            alpha=0.5,
        )
        ax.axvline(
            -lfc_cutoff,
            linestyle="--",
            linewidth=1.0,
            color="black",
            alpha=0.5,
        )
        ax.axhline(
            -np.log10(padj_cutoff),
            linestyle="--",
            linewidth=1.0,
            color="black",
            alpha=0.5,
        )

        texts = []

        for _, row in label_df.iterrows():
            text = ax.text(
                row["x_plot"],
                row["y_plot"],
                row["display_gene"],
                fontsize=8.5,
                fontweight="bold",
                ha="center",
                va="bottom",
                bbox=dict(
                    boxstyle="round,pad=0.18",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.78,
                ),
            )
            texts.append(text)

        if texts:
            try:
                from adjustText import adjust_text

                import contextlib
                import io

                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    adjust_text(
                        texts,
                        ax=ax,
                        expand_text=(1.08, 1.20),
                        expand_points=(1.05, 1.15),
                        force_text=(0.20, 0.35),
                        force_points=(0.05, 0.10),
                        only_move={"points": "y", "text": "xy"},
                        lim=200,
                    )
            except Exception:
                pass

        # --------------------------------------------------
        # Titles and labels
        # --------------------------------------------------
        ax.set_title(
            str(title),
            fontsize=15,
            fontweight="bold",
            pad=16,
        )

        ax.text(
            0.5,
            1.01,
            f"Positive logFC = higher in {group2}; negative logFC = higher in {group1}",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
        )

        ax.set_xlabel("log2 fold change", fontsize=12)
        ax.set_ylabel("-log10 adjusted p-value", fontsize=12)

        ax.set_xlim(-x_cap * 1.08, x_cap * 1.08)
        ax.set_ylim(-0.5, y_cap * 1.08)

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=True,
            fontsize=10,
            borderaxespad=0.0,
        )

        ax.grid(True, alpha=0.18, linewidth=0.6)

        ax.text(
            0.015,
            0.985,
            (
                f"FDR < {padj_cutoff}, |logFC| ≥ {lfc_cutoff}\n"
                f"Display capped at |logFC| ≤ {x_cap:.1f}, -log10 FDR ≤ {y_cap:.1f}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor="#BDBDBD",
                alpha=0.9,
            ),
        )

        fig.tight_layout(rect=[0, 0, 0.82, 1])

        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved figure: {output_path}")
        return output_path

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
