import gseapy as gp


class PathwayAnalyzer:

    def enrich_markers(
        self,
        marker_genes,
        gene_sets="GO_Biological_Process_2023"
    ):

        results = gp.enrichr(
            gene_list=marker_genes,
            gene_sets=gene_sets,
            organism="human",
            outdir=None
        )

        return results.results

    def top_pathways(
        self,
        marker_genes,
        n_terms=10
    ):

        pathways = self.enrich_markers(
            marker_genes
        )

        return pathways.head(n_terms)
