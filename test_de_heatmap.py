import scanpy as sc

from backend.condition_de import ConditionDE
from backend.de_heatmap import DEHeatmap


adata = sc.datasets.pbmc3k()

adata.obs["condition"] = (
    ["Healthy"] * 1300 +
    ["Disease"] * (adata.n_obs - 1300)
)

sc.pp.normalize_total(
    adata,
    target_sum=10000
)

sc.pp.log1p(adata)

de = ConditionDE()

adata = de.run(adata)

results = de.top_genes(
    adata,
    n_genes=10
)

genes = results["names"].tolist()

heatmap = DEHeatmap()

heatmap.create(
    adata,
    genes
)
