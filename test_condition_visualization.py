import pandas as pd
import scanpy as sc

from backend.condition_visualization import (
    ConditionVisualization
)

# --------------------------
# MOCK PROPORTIONS
# --------------------------

proportions = pd.DataFrame(
    {
        "Healthy": [0.50, 0.33, 0.17],
        "Disease": [0.17, 0.50, 0.33]
    },
    index=[
        "T cells",
        "Monocytes",
        "NK cells"
    ]
)

viz = ConditionVisualization()

viz.celltype_proportions(
    proportions
)

# --------------------------
# UMAP TEST
# --------------------------

adata = sc.datasets.pbmc3k()

sc.pp.normalize_total(
    adata,
    target_sum=10000
)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata
)

sc.tl.pca(adata)

sc.pp.neighbors(adata)

sc.tl.umap(adata)

adata.obs["condition"] = (
    ["Healthy"] * 1300 +
    ["Disease"] * (
        adata.n_obs - 1300
    )
)

viz.condition_umap(
    adata
)
