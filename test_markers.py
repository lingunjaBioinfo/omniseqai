from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer

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

cluster_map = annotator.cluster_annotations(
    adata
)

marker_engine = MarkerAnalyzer()

adata = marker_engine.find_markers(
    adata
)

summary = marker_engine.summarize_clusters(
    adata,
    cluster_map
)

for cluster, info in summary.items():

    print(
        f"\nCluster {cluster}"
    )

    print(
        f"Cell Type: {info['cell_type']}"
    )

    print(
        "Markers:"
    )

    print(
        ", ".join(info["markers"])
    )
