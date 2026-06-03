import scanpy as sc

from backend.condition_de import ConditionDE
from backend.condition_pathways import (
    ConditionPathways
)


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
    n_genes=100
)

pathways = ConditionPathways()

enrichment = pathways.analyze(
    results
)

print(
    enrichment[
        ["Term", "Adjusted P-value"]
    ].head(10)
)
