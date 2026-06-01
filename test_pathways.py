from backend.pipeline import SingleCellPipeline
from backend.differential_expression import DifferentialExpression
from backend.pathway_analysis import PathwayAnalyzer


pipeline = SingleCellPipeline()

adata = pipeline.load_pbmc3k()

pipeline.calculate_qc()

pipeline.filter_data()

pipeline.normalize_data()

pipeline.run_pca()

pipeline.compute_neighbors()

pipeline.run_umap()

pipeline.cluster_cells()

adata = pipeline.adata


de = DifferentialExpression()

adata = de.run(adata)

genes = de.get_top_genes(
    adata,
    cluster=1,
    n_genes=100
)

marker_list = genes[
    genes["logfoldchanges"] > 1
]["names"].tolist()

analyzer = PathwayAnalyzer()

results = analyzer.top_pathways(
    marker_list
)

print(results[
    ["Term", "Adjusted P-value"]
])
