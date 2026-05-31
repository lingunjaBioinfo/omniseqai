import scanpy as sc

from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer
from backend.interpretation import BiologicalInterpreter
from backend.report import ReportGenerator
from backend.visualization import Visualizer

sc.settings.figdir = "outputs"

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

interpreter = BiologicalInterpreter()

for cluster in summary:

    summary[cluster][
        "interpretation"
    ] = interpreter.generate_summary(
        cluster,
        summary[cluster]["cell_type"],
        summary[cluster]["markers"]
    )

reporter = ReportGenerator()

report = reporter.generate_report(
    summary
)

print(report)

reporter.save_report(report)

# Generate figures

visualizer = Visualizer()

visualizer.save_qc_plots(adata)

visualizer.save_umap_clusters(adata)

visualizer.save_umap_celltypes(adata)

print("\nFigures saved to outputs/")

visualizer = Visualizer()

visualizer.save_qc_plots(adata)

visualizer.save_umap_clusters(adata)

visualizer.save_umap_celltypes(adata)

print("\nFigures saved to outputs/")

reporter = ReportGenerator()

report = reporter.generate_report(
    summary
)

print(report)

reporter.save_report(report)

reporter.save_pdf_report(
    summary
)
