import scanpy as sc


class DifferentialExpression:

    def run(
        self,
        adata,
        groupby="leiden"
    ):

        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon"
        )

        print(
            "\nDifferential expression completed."
        )

        return adata

    def get_top_genes(
        self,
        adata,
        cluster,
        n_genes=20
    ):

        genes = sc.get.rank_genes_groups_df(
            adata,
            group=str(cluster)
        )

        return genes.head(n_genes)
