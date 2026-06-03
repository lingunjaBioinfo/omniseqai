import os
import scanpy as sc
import pandas as pd

from backend.gene_mapper import GeneMapper
from backend.condition_comparison import ConditionComparison
from backend.condition_visualization import ConditionVisualization
from backend.condition_de import ConditionDE
from backend.condition_pathways import ConditionPathways
from backend.condition_interpreter import ConditionInterpreter
from backend.volcano_plot import VolcanoPlot
from backend.de_heatmap import DEHeatmap


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"
OUT_DIR = "outputs"
REPORT_FILE = "reports/real_status_analysis.txt"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs("reports", exist_ok=True)

print("\nLoading real COVID dataset...")
adata = sc.read_h5ad(DATA_FILE)

print("\nDataset loaded.")
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")

# --------------------------------------------------
# Fix gene identifiers
# --------------------------------------------------
mapper = GeneMapper()
adata = mapper.fix_gene_names(adata)

# --------------------------------------------------
# Use dataset-native cell type labels
# --------------------------------------------------
if "celltype.final" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["celltype.final"].astype(str)
elif "cell_type" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
else:
    raise ValueError("No usable cell type annotation found in the dataset.")

# --------------------------------------------------
# Build status condition from cv19_vax_boost_or_HC_status
# --------------------------------------------------
status_col = "cv19_vax_boost_or_HC_status"
if status_col not in adata.obs.columns:
    raise ValueError(f"No '{status_col}' column found in adata.obs.")

adata.obs["condition"] = adata.obs[status_col].astype(str)

# Normalize category names for consistency
adata.obs["condition"] = adata.obs["condition"].replace(
    {
        "cv_19": "COVID",
        "HC": "Healthy",
        "Vax": "Vaccinated",
        "booster": "Vaccinated"
    }
)

# Keep only the three main groups
adata = adata[
    adata.obs["condition"].isin(["Healthy", "COVID", "Vaccinated"])
].copy()

print("\nCondition counts:")
print(adata.obs["condition"].value_counts())

print("\nCell type counts:")
print(adata.obs["cell_type"].value_counts().head(15))

# --------------------------------------------------
# Cell-type proportions by condition
# --------------------------------------------------
comparison = ConditionComparison()
proportions = comparison.compare_celltypes(
    adata,
    condition_column="condition",
    celltype_column="cell_type"
)

print("\n===== CELL TYPE PROPORTIONS =====\n")
print(proportions)

viz = ConditionVisualization()
viz.celltype_proportions(
    proportions,
    output_file=os.path.join(OUT_DIR, "celltype_proportions_status.png")
)

if "X_umap" in adata.obsm:
    viz.condition_umap(
        adata,
        condition_key="condition",
        output_file=os.path.join(OUT_DIR, "condition_umap_status.png")
    )
else:
    print("\nNo UMAP found in the dataset; skipping condition UMAP.")

# --------------------------------------------------
# Helper for pairwise DE
# --------------------------------------------------
de_engine = ConditionDE()
pathway_engine = ConditionPathways()
interpreter = ConditionInterpreter()
volcano = VolcanoPlot()
heatmap = DEHeatmap()

def run_pairwise_de(
    adata_in,
    condition_key,
    group1,
    group2,
    label
):
    """
    Run DE only on cells in two chosen groups.
    """
    subset = adata_in[
        adata_in.obs[condition_key].isin([group1, group2])
    ].copy()

    subset.obs[condition_key] = pd.Categorical(
        subset.obs[condition_key],
        categories=[group2, group1]
    )

    print(f"\nRunning DE: {group1} vs {group2} ({label})")
    print(subset.obs[condition_key].value_counts())

    subset = de_engine.run(
        subset,
        condition_key=condition_key,
        group1=group1,
        group2=group2
    )

    de_results = de_engine.top_genes(
        subset,
        group=group1,
        n_genes=100
    )

    if de_results is None or de_results.empty:
        sig = de_results
    else:
        sig = de_results[
            (de_results["pvals_adj"] < 0.05) &
            (de_results["logfoldchanges"].abs() > 0.5)
        ].copy()

    if sig is None or sig.empty:
        sig_for_pathways = de_results
    else:
        sig_for_pathways = sig

    pathways = pathway_engine.analyze(
        sig_for_pathways,
        top_n=50
    )

    if pathways is not None and not pathways.empty:
        interpretation = interpreter.interpret(pathways)
    else:
        interpretation = ["No strong condition-specific pathway signature identified."]

    return {
        "subset": subset,
        "de_results": de_results,
        "sig_genes": sig,
        "pathways": pathways,
        "interpretation": interpretation,
        "label": label,
        "group1": group1,
        "group2": group2
    }

# --------------------------------------------------
# Run pairwise comparisons
# --------------------------------------------------
comparisons = {}

comparisons["COVID_vs_Healthy"] = run_pairwise_de(
    adata,
    "condition",
    "COVID",
    "Healthy",
    "COVID vs Healthy"
)

comparisons["Vaccinated_vs_Healthy"] = run_pairwise_de(
    adata,
    "condition",
    "Vaccinated",
    "Healthy",
    "Vaccinated vs Healthy"
)

comparisons["COVID_vs_Vaccinated"] = run_pairwise_de(
    adata,
    "condition",
    "COVID",
    "Vaccinated",
    "COVID vs Vaccinated"
)

# --------------------------------------------------
# Choose the strongest comparison for plots
# --------------------------------------------------
best_key = None
best_count = -1

for key, result in comparisons.items():
    sig_count = 0 if result["sig_genes"] is None or result["sig_genes"].empty else len(result["sig_genes"])
    if sig_count > best_count:
        best_count = sig_count
        best_key = key

best_result = comparisons[best_key]
best_subset = best_result["subset"]
best_de_results = best_result["de_results"]
best_sig_genes = best_result["sig_genes"]

safe_name = best_key.lower().replace(" ", "_").replace("/", "_")

if best_de_results is not None and not best_de_results.empty:
    volcano.create(
        best_de_results,
        output_file=os.path.join(
            OUT_DIR,
            f"volcano_{safe_name}.png"
        ),
        title=f"Volcano plot: {best_result['label']}"
    )

if best_subset is not None and best_sig_genes is not None and not best_sig_genes.empty:
    top_heatmap_genes = best_sig_genes["names"].head(15).tolist()

    heatmap.create(
        best_subset,
        top_heatmap_genes,
        groupby="condition",
        output_file=os.path.join(
            OUT_DIR,
            f"heatmap_{safe_name}.png"
        )
    )

# --------------------------------------------------
# Save report
# --------------------------------------------------
report_lines = []
report_lines.append("========== OMNISEQAI REAL STATUS ANALYSIS ==========\n")
report_lines.append(f"Cells: {adata.n_obs:,}")
report_lines.append(f"Genes: {adata.n_vars:,}\n")

report_lines.append("Condition counts:")
report_lines.append(adata.obs["condition"].value_counts().to_string())
report_lines.append("\nCell type proportions by condition:")
report_lines.append(proportions.to_string())

for key, result in comparisons.items():
    report_lines.append("\n" + "=" * 70)
    report_lines.append(f"Comparison: {result['label']}")

    report_lines.append("\nInterpretation:")
    for item in result["interpretation"]:
        report_lines.append(f"- {item}")

    report_lines.append("\nTop DE genes:")
    if result["sig_genes"] is not None and not result["sig_genes"].empty:
        report_lines.append(
            result["sig_genes"][["names", "logfoldchanges", "pvals_adj"]]
            .head(10)
            .to_string(index=False)
        )
    else:
        report_lines.append("No significant DE genes detected.")

with open(REPORT_FILE, "w") as f:
    f.write("\n".join(report_lines))

print(f"\nReport saved: {REPORT_FILE}")

# --------------------------------------------------
# Console summary
# --------------------------------------------------
print("\n===== STATUS-AWARE SUMMARY =====")
for key, result in comparisons.items():
    print("\n" + "-" * 60)
    print(f"Comparison: {result['label']}")
    print("Interpretation:")
    for item in result["interpretation"]:
        print(f"- {item}")

print("\nAnalysis complete.")
