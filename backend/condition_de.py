import pandas as pd
import scanpy as sc


class ConditionDE:

    def _preprocess_for_de(self, adata):
        """
        Ensure the matrix is suitable for rank_genes_groups.
        We only log-normalize if the object does not already advertise log1p.
        """
        adata = adata.copy()

        if "log1p" not in adata.uns:
            sc.pp.normalize_total(
                adata,
                target_sum=10000
            )
            sc.pp.log1p(adata)

        return adata

    def run(
        self,
        adata,
        condition_key="condition",
        group1="Disease",
        group2="Healthy"
    ):

        adata = self._preprocess_for_de(adata)

        sc.tl.rank_genes_groups(
            adata,
            groupby=condition_key,
            groups=[group1],
            reference=group2,
            method="wilcoxon",
            use_raw=False
        )

        print(
            f"\nCondition DE completed: "
            f"{group1} vs {group2}"
        )

        return adata

    def top_genes(
        self,
        adata,
        group="Disease",
        n_genes=20
    ):

        genes = sc.get.rank_genes_groups_df(
            adata,
            group=group
        )

        if genes is None or genes.empty:
            return genes

        genes = genes.copy()

        genes["names"] = (
            genes["names"]
            .astype(str)
            .str.strip()
        )

        for col in ["logfoldchanges", "pvals_adj", "scores"]:
            if col in genes.columns:
                genes[col] = pd.to_numeric(
                    genes[col],
                    errors="coerce"
                )

        genes = genes.sort_values(
            ["pvals_adj", "logfoldchanges"],
            ascending=[True, False],
            kind="mergesort"
        )

        return genes.head(n_genes).reset_index(drop=True)
