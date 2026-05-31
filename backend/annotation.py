import celltypist
import pandas as pd


class CellAnnotator:

    def __init__(self):

        self.model = "Immune_All_Low.pkl"

    def annotate(self, adata):

        predictions = celltypist.annotate(
            adata,
            model=self.model
        )

        adata.obs["cell_type"] = (
            predictions.predicted_labels
        )

        return adata

    def summary(self, adata):

        print("\n===== CELL TYPES =====\n")

        print(
            adata.obs["cell_type"]
            .value_counts()
        )

    def cluster_annotations(self, adata):

        cluster_map = {}

        for cluster in adata.obs["leiden"].unique():

            subset = adata.obs[
                adata.obs["leiden"] == cluster
            ]

            major_type = (
                subset["cell_type"]
                .value_counts()
                .idxmax()
            )

            cluster_map[cluster] = major_type

        return cluster_map
