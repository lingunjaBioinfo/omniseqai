import os
import scanpy as sc


class Visualizer:

    def __init__(self):

        os.makedirs(
            "outputs",
            exist_ok=True
        )

    def save_qc_plots(self, adata):

        sc.pl.violin(
            adata,
            [
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt"
            ],
            multi_panel=True,
            show=False,
            save="_qc_violin.png"
        )

    def save_umap_clusters(self, adata):

        sc.pl.umap(
            adata,
            color="leiden",
            show=False,
            save="_clusters.png"
        )

    def save_umap_celltypes(self, adata):

        sc.pl.umap(
            adata,
            color="cell_type",
            legend_loc="on data",
            show=False,
            save="_celltypes.png"
        )
