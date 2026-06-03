import os

import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp


class DEHeatmap:

    def _resolve_gene_names(self, adata, genes):
        """
        Resolve requested gene names against:
        - adata.var_names
        - adata.var['gene_symbol']
        - adata.var['feature_name']
        """
        requested = [
            str(g).strip()
            for g in genes
            if pd.notna(g) and str(g).strip() != ""
        ]

        resolved_var_names = []
        display_names = []

        # Direct matches first
        for g in requested:
            if g in adata.var_names and g not in display_names:
                resolved_var_names.append(g)
                display_names.append(g)

        # Alias mapping from metadata columns
        alias_cols = []
        if "gene_symbol" in adata.var.columns:
            alias_cols.append("gene_symbol")
        if "feature_name" in adata.var.columns:
            alias_cols.append("feature_name")

        if alias_cols:
            var_df = adata.var.copy()
            for g in requested:
                if g in display_names:
                    continue

                matched = False
                for col in alias_cols:
                    hits = var_df.index[
                        var_df[col].astype(str).str.strip() == g
                    ].tolist()
                    if hits:
                        resolved_var_names.append(hits[0])
                        display_names.append(g)
                        matched = True
                        break

                if matched:
                    continue

        # Remove duplicates while preserving order
        seen = set()
        unique_resolved = []
        unique_display = []

        for var_name, disp in zip(resolved_var_names, display_names):
            if var_name not in seen:
                seen.add(var_name)
                unique_resolved.append(var_name)
                unique_display.append(disp)

        return unique_resolved, unique_display

    def create(
        self,
        adata,
        genes,
        groupby="condition",
        output_file="outputs/de_heatmap.png"
    ):

        if adata is None or len(genes) == 0:
            print("No genes available for heatmap.")
            return

        resolved_var_names, display_names = self._resolve_gene_names(
            adata,
            genes
        )

        if len(resolved_var_names) == 0:
            print("No valid genes found in adata.var_names.")
            print("Requested genes:", [str(g) for g in genes[:10]])
            print("Available genes:", list(map(str, adata.var_names[:10])))
            return

        X = adata[:, resolved_var_names].X

        if sp.issparse(X):
            X = X.toarray()

        expr = pd.DataFrame(
            X,
            columns=display_names,
            index=adata.obs_names
        )

        expr[groupby] = adata.obs[groupby].astype(str).values

        group_order = []
        for g in ["Healthy", "Disease"]:
            if g in expr[groupby].unique():
                group_order.append(g)

        for g in expr[groupby].unique():
            if g not in group_order:
                group_order.append(g)

        mean_expr = expr.groupby(groupby)[display_names].mean().loc[group_order]

        z = (
            mean_expr - mean_expr.mean(axis=0)
        ) / mean_expr.std(axis=0).replace(0, 1)

        plt.figure(
            figsize=(
                max(8, 0.45 * len(display_names)),
                max(3.5, 0.7 * len(group_order))
            )
        )

        plt.imshow(z.values, aspect="auto", cmap="viridis")
        plt.yticks(range(len(group_order)), group_order)
        plt.xticks(range(len(display_names)), display_names, rotation=90)
        plt.colorbar(label="Z-score")
        plt.title("Condition DE Heatmap")
        plt.tight_layout()

        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"\nHeatmap saved: {output_file}")
