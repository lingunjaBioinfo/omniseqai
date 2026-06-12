from __future__ import annotations

import os
from typing import Iterable, List, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp


class DEHeatmap:

    def _resolve_gene_names(self, adata, genes: Sequence[str]) -> List[str]:
        requested = [str(g).strip() for g in genes if pd.notna(g) and str(g).strip()]

        resolved = []
        available = set(map(str, adata.var_names))

        # Direct matches first
        for g in requested:
            if g in available and g not in resolved:
                resolved.append(g)

        # Alias matches through metadata columns
        alias_cols = [c for c in ["gene_symbol", "feature_name"] if c in adata.var.columns]
        if alias_cols:
            var_df = adata.var.copy()
            for g in requested:
                if g in resolved:
                    continue
                for col in alias_cols:
                    hits = var_df.index[var_df[col].astype(str).str.strip() == g].tolist()
                    if hits:
                        if hits[0] not in resolved:
                            resolved.append(hits[0])
                        break

        return resolved

    def create(
        self,
        adata,
        genes,
        groupby="condition",
        output_file="outputs/de_heatmap.png",
        max_genes: int = 15,
        cmap: str = "RdBu_r",
    ):
        if adata is None or len(genes) == 0:
            print("No genes available for heatmap.")
            return

        resolved = self._resolve_gene_names(adata, genes)
        if len(resolved) == 0:
            print("No valid genes found in adata.var_names.")
            return

        # Limit to top genes for readability
        resolved = resolved[:max_genes]

        if groupby not in adata.obs.columns:
            print(f"Heatmap skipped: '{groupby}' not found in adata.obs")
            return

        X = adata[:, resolved].X
        if sp.issparse(X):
            X = X.toarray()

        expr = pd.DataFrame(X, columns=resolved, index=adata.obs_names)
        expr[groupby] = adata.obs[groupby].astype(str).values

        group_order = list(dict.fromkeys(expr[groupby].tolist()))
        preferred = ["Healthy", "COVID", "Vaccinated", "Disease", "Control"]
        ordered = [g for g in preferred if g in group_order] + [g for g in group_order if g not in preferred]

        mean_expr = expr.groupby(groupby)[resolved].mean().loc[ordered]

        # z-score by gene
        z = (mean_expr - mean_expr.mean(axis=0)) / mean_expr.std(axis=0).replace(0, 1)
        z = z.fillna(0.0)

        # Order genes by average absolute deviation across groups
        gene_order = z.abs().mean(axis=0).sort_values(ascending=False).index.tolist()
        z = z[gene_order]

        fig_w = max(8.5, 0.45 * len(gene_order))
        fig_h = max(3.8, 0.65 * len(ordered))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=140)

        im = ax.imshow(z.values, aspect="auto", cmap=cmap, interpolation="nearest")
        ax.set_title("Condition DE Heatmap", fontsize=16, pad=12)
        ax.set_yticks(range(len(ordered)))
        ax.set_yticklabels(ordered, fontsize=11)
        ax.set_xticks(range(len(gene_order)))
        ax.set_xticklabels(gene_order, rotation=90, fontsize=9)

        cbar = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
        cbar.set_label("Z-score", fontsize=11)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()

        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        fig.savefig(output_file, bbox_inches="tight")
        plt.close(fig)

        print(f"\nHeatmap saved: {output_file}")
