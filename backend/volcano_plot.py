from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class VolcanoPlot:

    def create(
        self,
        de_results: pd.DataFrame,
        output_file: str = "outputs/volcano_plot.png",
        title: str = "Condition Differential Expression",
        pval_col: str = "pvals_adj",
        lfc_col: str = "logfoldchanges",
        gene_col: str = "names",
        sig_p: float = 0.05,
        sig_lfc: float = 0.5,
        annotate_top_n: int = 8,
    ) -> None:
        if de_results is None or de_results.empty:
            print("No DE results available for volcano plot.")
            return

        required = {pval_col, lfc_col, gene_col}
        missing = required - set(de_results.columns)
        if missing:
            print(f"Volcano plot skipped. Missing columns: {sorted(missing)}")
            return

        df = de_results.copy()
        df[gene_col] = df[gene_col].astype(str).str.strip()
        df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce")
        df[lfc_col] = pd.to_numeric(df[lfc_col], errors="coerce")
        df = df.dropna(subset=[pval_col, lfc_col, gene_col])

        if df.empty:
            print("No valid rows available for volcano plot.")
            return

        df["neglog10_padj"] = -np.log10(np.clip(df[pval_col].astype(float), 1e-300, 1.0))
        df["sig"] = (df[pval_col] < sig_p) & (df[lfc_col].abs() >= sig_lfc)

        # Clip extreme values only for display, not analysis
        abs_lfc = df[lfc_col].abs().replace([np.inf, -np.inf], np.nan).dropna()
        if len(abs_lfc) > 0:
            clip_x = float(np.nanpercentile(abs_lfc, 99))
            clip_x = max(2.0, min(clip_x, 8.0))
        else:
            clip_x = 4.0

        df["plot_lfc"] = df[lfc_col].clip(-clip_x, clip_x)

        fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=140)

        nonsig = df[~df["sig"]]
        sig = df[df["sig"]]

        ax.scatter(
            nonsig["plot_lfc"],
            nonsig["neglog10_padj"],
            s=16,
            alpha=0.35,
            linewidths=0,
            label="Not significant",
        )

        if not sig.empty:
            up = sig[sig[lfc_col] > 0]
            down = sig[sig[lfc_col] < 0]

            if not down.empty:
                ax.scatter(
                    down["plot_lfc"],
                    down["neglog10_padj"],
                    s=18,
                    alpha=0.8,
                    linewidths=0,
                    label="Significant down",
                )

            if not up.empty:
                ax.scatter(
                    up["plot_lfc"],
                    up["neglog10_padj"],
                    s=18,
                    alpha=0.8,
                    linewidths=0,
                    label="Significant up",
                )

        ax.axhline(-np.log10(sig_p), linestyle="--", linewidth=1, color="grey")
        ax.axvline(sig_lfc, linestyle="--", linewidth=1, color="grey")
        ax.axvline(-sig_lfc, linestyle="--", linewidth=1, color="grey")

        ax.set_title(title, fontsize=16, pad=12)
        ax.set_xlabel("Log fold change", fontsize=12)
        ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)

        ax.set_xlim(-clip_x * 1.05, clip_x * 1.05)
        ymax = float(df["neglog10_padj"].max())
        ax.set_ylim(0, ymax * 1.10 if ymax > 0 else 1.0)

        ax.grid(True, alpha=0.18, linewidth=0.6)
        ax.legend(frameon=False, loc="upper right")

        # Annotate top genes
        if not sig.empty:
            ann = sig.sort_values(
                [pval_col, lfc_col],
                ascending=[True, False],
                kind="mergesort"
            ).head(annotate_top_n)
        else:
            ann = df.sort_values(
                [pval_col, lfc_col],
                ascending=[True, False],
                kind="mergesort"
            ).head(annotate_top_n)

        for _, row in ann.iterrows():
            ax.annotate(
                row[gene_col],
                (row["plot_lfc"], row["neglog10_padj"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
                alpha=0.9,
            )

        fig.tight_layout()

        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        fig.savefig(output_file, bbox_inches="tight")
        plt.close(fig)

        print(f"\nVolcano plot saved: {output_file}")
