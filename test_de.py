from backend.pipeline import SingleCellPipeline
from backend.differential_expression import DifferentialExpression


pipeline = SingleCellPipeline()

pipeline.load_pbmc3k()

pipeline.calculate_qc()

pipeline.filter_data()

pipeline.normalize_data()

pipeline.run_pca()

pipeline.compute_neighbors()

pipeline.run_umap()

pipeline.cluster_cells()

# Get the processed AnnData object
adata = pipeline.adata

print("\nOBS COLUMNS:")
print(adata.obs.columns)

de = DifferentialExpression()

adata = de.run(adata)

top = de.get_top_genes(
    adata,
    cluster=0
)

print(top.head(10))
