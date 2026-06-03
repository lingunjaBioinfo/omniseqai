import scanpy as sc
import numpy as np

from backend.condition_de import ConditionDE
from backend.de_interpreter import DEInterpreter


# --------------------------
# LOAD DATA
# --------------------------

adata = sc.datasets.pbmc3k()

# --------------------------
# NORMALIZE
# --------------------------

sc.pp.normalize_total(
    adata,
    target_sum=10000
)

sc.pp.log1p(adata)

adata.X = adata.X.toarray()

# --------------------------
# CONDITIONS
# --------------------------

adata.obs["condition"] = "Healthy"

disease_cells = np.random.choice(
    adata.obs_names,
    size=500,
    replace=False
)

adata.obs.loc[
    disease_cells,
    "condition"
] = "Disease"

# --------------------------
# SIGNAL
# --------------------------

genes_to_modify = [
    "IL32",
    "LTB",
    "LDHB"
]

disease_mask = (
    adata.obs["condition"] == "Disease"
).values

for gene in genes_to_modify:

    idx = adata.var_names.get_loc(gene)

    adata.X[disease_mask, idx] += 5

# --------------------------
# DE
# --------------------------

de = ConditionDE()

adata = de.run(adata)

results = de.top_genes(
    adata,
    n_genes=20
)

# --------------------------
# INTERPRETATION
# --------------------------

interpreter = DEInterpreter()

findings = interpreter.interpret(
    results
)

print("\n===== DE INTERPRETATION =====\n")

for finding in findings:

    print(f"- {finding}")
