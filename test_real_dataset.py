import scanpy as sc

from backend.gene_mapper import GeneMapper
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer
from backend.disease_interpreter import DiseaseInterpreter
from backend.cell_communication import CellCommunication


print("\nLoading dataset...")

adata = sc.read_h5ad("data/covid_pbmc/covid_pbmc.h5ad")

print(f"\nCells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")

mapper = GeneMapper()
adata = mapper.fix_gene_names(adata)

print("\nPreprocessing dataset...")

# Only normalize if the matrix looks raw
try:
    max_val = adata.X.max()
except Exception:
    max_val = None

if max_val is not None and max_val > 50:
    print("Raw counts detected; normalizing and log-transforming.")
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
else:
    print("Data already normalized.")

# Prefer existing metadata labels from the dataset
if "celltype.final" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["celltype.final"].astype(str)
elif "cell_type" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
else:
    # fallback to CellTypist only if no cell annotations exist
    annotator = CellAnnotator()
    adata = annotator.annotate(adata)

# cluster on the existing embedding / normalized data
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata = adata[:, adata.var.highly_variable].copy()

sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

print(f"\nClusters: {adata.obs['leiden'].nunique()}")

# cluster map from whichever annotation is available
cluster_map = {}
for cluster in adata.obs["leiden"].unique():
    subset = adata.obs[adata.obs["leiden"] == cluster]
    major_type = subset["cell_type"].value_counts().idxmax()
    cluster_map[str(cluster)] = major_type

marker = MarkerAnalyzer()
adata = marker.find_markers(adata, groupby="leiden")
summary = marker.summarize_clusters(adata, cluster_map)

interpreter = DiseaseInterpreter()
for cluster in summary:
    summary[cluster]["disease_interpretation"] = interpreter.interpret(
        summary[cluster]["markers"]
    )

communication = CellCommunication().analyze(summary)

for cluster in summary:
    print("\n")
    print("=" * 70)
    print(f"Cluster {cluster}")
    print(f"Cell Type: {summary[cluster]['cell_type']}")
    print(f"Markers: {', '.join(summary[cluster]['markers'][:5])}")

    print("\nDisease:")
    for item in summary[cluster]["disease_interpretation"]:
        print(f"- {item}")

    print("\nCommunication:")
    for signal in communication[cluster]["signals"]:
        print(f"- {signal}")

print("\nAnalysis complete.")
