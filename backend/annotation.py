import celltypist
import pandas as pd


class CellAnnotator:

    def annotate(self, adata, model="Immune_All_Low.pkl"):
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=True
        )

        labels = predictions.predicted_labels
        if isinstance(labels, pd.DataFrame) and "majority_voting" in labels.columns:
            adata.obs["cell_type"] = labels["majority_voting"].astype(str).values
        else:
            adata.obs["cell_type"] = labels.astype(str).values

        return adata

    def cluster_annotations(self, adata, cluster_key="leiden"):
        cluster_map = {}

        for cluster in adata.obs[cluster_key].unique():
            subset = adata.obs[adata.obs[cluster_key] == cluster]
            major_type = subset["cell_type"].value_counts().idxmax()
            cluster_map[str(cluster)] = major_type

        print("\nCluster annotations complete.")
        return cluster_map
