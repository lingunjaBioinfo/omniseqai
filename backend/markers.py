import scanpy as sc
import pandas as pd


class MarkerAnalyzer:

    def __init__(self):
        pass

    def find_markers(
        self,
        adata,
        groupby="leiden"
    ):

        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon"
        )

        return adata

    def top_markers(
        self,
        adata,
        cluster,
        n_genes=10
    ):

        result = adata.uns["rank_genes_groups"]

        genes = result["names"][str(cluster)]

        return list(genes[:n_genes])

    def summarize_clusters(
        self,
        adata,
        cluster_map,
        n_genes=10
    ):

        summary = {}

        for cluster in cluster_map:

            markers = self.top_markers(
                adata,
                cluster,
                n_genes
            )

            summary[cluster] = {
                "cell_type": cluster_map[cluster],
                "markers": markers
            }

        return summary
