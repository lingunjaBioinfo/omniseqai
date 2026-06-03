import scanpy as sc
import pandas as pd

from backend.gene_mapper import GeneMapper
from backend.condition_de import ConditionDE
from backend.condition_pathways import ConditionPathways
from backend.condition_interpreter import ConditionInterpreter

print("\nLoading dataset...")

adata = sc.read_h5ad(
    "data/covid_pbmc/covid_pbmc.h5ad"
)

mapper = GeneMapper()
adata = mapper.fix_gene_names(adata)

adata.obs["condition"] = (
    adata.obs["disease"]
    .astype(str)
)

adata.obs["condition"] = adata.obs[
    "condition"
].replace({
    "normal": "Healthy",
    "COVID-19": "COVID"
})

adata.obs["cell_type"] = (
    adata.obs["celltype.final"]
    .astype(str)
)

print("\nCondition counts:")
print(
    adata.obs["condition"]
    .value_counts()
)

de_engine = ConditionDE()
pathway_engine = ConditionPathways()
interpreter = ConditionInterpreter()

results = []

celltypes = (
    adata.obs["cell_type"]
    .value_counts()
)

celltypes = celltypes[
    celltypes >= 200
].index.tolist()

for celltype in celltypes:

    print(
        f"\nAnalyzing {celltype}"
    )

    subset = adata[
        adata.obs["cell_type"]
        == celltype
    ].copy()

    counts = (
        subset.obs["condition"]
        .value_counts()
    )

    if (
        "COVID" not in counts
        or
        "Healthy" not in counts
    ):
        continue

    if (
        counts["COVID"] < 50
        or
        counts["Healthy"] < 50
    ):
        continue

    try:

        subset = de_engine.run(
            subset,
            condition_key="condition",
            group1="COVID",
            group2="Healthy"
        )

        de = de_engine.top_genes(
            subset,
            group="COVID",
            n_genes=200
        )

        sig = de[
            (de["pvals_adj"] < 0.05)
            &
            (abs(de["logfoldchanges"]) > 0.5)
        ]

        score = len(sig)

        pathways = pathway_engine.analyze(
            sig,
            top_n=50
        )

        interpretation = (
            interpreter.interpret(
                pathways
            )
        )

        results.append({

            "cell_type":
                celltype,

            "covid_cells":
                counts["COVID"],

            "healthy_cells":
                counts["Healthy"],

            "de_genes":
                score,

            "interpretation":
                interpretation

        })

    except Exception as e:

        print(
            f"Failed: {celltype}"
        )

        print(e)

results = pd.DataFrame(
    results
)

results = results.sort_values(
    "de_genes",
    ascending=False
)

print("\n")
print("=" * 80)
print("TOP COVID-RESPONSIVE CELL TYPES")
print("=" * 80)

print(
    results[
        [
            "cell_type",
            "de_genes",
            "covid_cells",
            "healthy_cells"
        ]
    ]
)

results.to_csv(
    "reports/celltype_covid_response_ranking.csv",
    index=False
)

print(
    "\nSaved:"
    " reports/celltype_covid_response_ranking.csv"
)
