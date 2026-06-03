import scanpy as sc
import pandas as pd


class MarkerAnalyzer:

    def find_markers(self, adata, groupby="leiden"):
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            use_raw=False
        )
        print("\nMarker analysis complete.")
        return adata

    def _clean_markers(self, genes):
        bad_prefixes = ("RPL", "RPS", "MT-")
        bad_genes = {"MALAT1", "NEAT1", "XIST"}

        genes = genes[
            ~genes["names"].astype(str).str.startswith(bad_prefixes)
        ]
        genes = genes[
            ~genes["names"].astype(str).isin(bad_genes)
        ]
        genes = genes[
            ~genes["names"].astype(str).str.startswith("ENSG")
        ]
        return genes

    def summarize_clusters(self, adata, cluster_map, n_genes=10):
        summary = {}

        for cluster in adata.obs["leiden"].cat.categories:
            genes = sc.get.rank_genes_groups_df(adata, group=cluster)
            genes = self._clean_markers(genes)

            markers = genes["names"].head(n_genes).tolist()

            summary[str(cluster)] = {
                "cell_type": cluster_map.get(str(cluster), "Unknown"),
                "markers": markers
            }

        return summary
