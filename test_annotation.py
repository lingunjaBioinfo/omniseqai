from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator

pipeline = SingleCellPipeline()

pipeline.load_pbmc3k()

pipeline.calculate_qc()

pipeline.filter_data()

pipeline.normalize_data()

pipeline.identify_hvg()

pipeline.run_pca()

pipeline.compute_neighbors()

pipeline.run_umap()

pipeline.cluster_cells()

annotator = CellAnnotator()

adata = annotator.annotate(
    pipeline.adata
)

annotator.summary(adata)

cluster_map = annotator.cluster_annotations(
    adata
)

print("\n===== CLUSTER LABELS =====\n")

for cluster, label in cluster_map.items():

    print(
        f"Cluster {cluster}: {label}"
    )
