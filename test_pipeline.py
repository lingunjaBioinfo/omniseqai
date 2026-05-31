from backend.pipeline import SingleCellPipeline

pipeline = SingleCellPipeline()

# Load data
pipeline.load_pbmc3k()

# Dataset summary
pipeline.summary()

# QC
pipeline.calculate_qc()
pipeline.qc_summary()

# Filtering
pipeline.filter_data()
pipeline.filter_summary()

# Normalization
pipeline.normalize_data()

# HVG
pipeline.identify_hvg()

# PCA
pipeline.run_pca()

# Neighbor graph
pipeline.compute_neighbors()

# UMAP
pipeline.run_umap()

# Clustering
pipeline.cluster_cells()
