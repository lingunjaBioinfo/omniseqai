import os

import scanpy as sc

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
REPORT_FILE = "reports/real_condition_analysis.txt"

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
# Use dataset-native labels
# --------------------------------------------------
if "celltype.final" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["celltype.final"].astype(str)
elif "cell_type" in adata.obs.columns:
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
else:
    raise ValueError("No usable cell type annotation found in the dataset.")

if "disease" not in adata.obs.columns:
    raise ValueError("No 'disease' column found in adata.obs.")

adata.obs["condition"] = adata.obs["disease"].astype(str)

adata = adata[adata.obs["condition"].isin(["normal", "COVID-19"])].copy()
adata.obs["condition"] = adata.obs["condition"].replace(
    {"normal": "Healthy", "COVID-19": "Disease"}
)

print("\nCondition counts:")
print(adata.obs["condition"].value_counts())

print("\nCell type counts:")
print(adata.obs["cell_type"].value_counts().head(15))

# --------------------------------------------------
# Condition comparison across cell types
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
    output_file=os.path.join(OUT_DIR, "celltype_proportions_real.png")
)

if "X_umap" in adata.obsm:
    viz.condition_umap(
        adata,
        condition_key="condition",
        output_file=os.path.join(OUT_DIR, "condition_umap_real.png")
    )
else:
    print("\nNo UMAP found in the dataset; skipping condition UMAP.")

# --------------------------------------------------
# Cell-type-specific DE: Disease vs Healthy
# --------------------------------------------------
de_engine = ConditionDE()
pathway_engine = ConditionPathways()
interpreter = ConditionInterpreter()
volcano = VolcanoPlot()
heatmap = DEHeatmap()

celltype_order = adata.obs["cell_type"].value_counts().index.tolist()
results_by_celltype = {}

representative_celltype = None
representative_sig_count = -1
representative_de_results = None
representative_subset = None
representative_top_genes = []

for cell_type in celltype_order:
    subset = adata[adata.obs["cell_type"] == cell_type].copy()

    condition_counts = subset.obs["condition"].value_counts()
    if condition_counts.get("Healthy", 0) < 20 or condition_counts.get("Disease", 0) < 20:
        continue

    print(f"\nRunning DE for cell type: {cell_type}")
    print(condition_counts)

    subset = de_engine.run(
        subset,
        condition_key="condition",
        group1="Disease",
        group2="Healthy"
    )

    de_results = de_engine.top_genes(
        subset,
        group="Disease",
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
        pathway_terms = pathways["Term"].head(10).tolist()
        interpretation = interpreter.interpret(pathways)
    else:
        pathway_terms = []
        interpretation = ["No strong condition-specific pathway signature identified."]

    results_by_celltype[cell_type] = {
        "condition_counts": condition_counts.to_dict(),
        "de_results": de_results,
        "sig_genes": sig,
        "pathways": pathways,
        "interpretation": interpretation
    }

    sig_count = 0 if sig is None or sig.empty else len(sig)
    if sig_count > representative_sig_count:
        representative_sig_count = sig_count
        representative_celltype = cell_type
        representative_de_results = de_results
        representative_subset = subset
        if sig is not None and not sig.empty:
            representative_top_genes = sig["names"].head(15).tolist()
        else:
            representative_top_genes = de_results["names"].head(15).tolist()

# --------------------------------------------------
# Representative plots
# --------------------------------------------------
if representative_celltype is not None and representative_de_results is not None:
    safe_name = (
        representative_celltype
        .replace("/", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )

    volcano.create(
        representative_de_results,
        output_file=os.path.join(
            OUT_DIR,
            f"volcano_{safe_name}_covid_vs_healthy.png"
        ),
        title=f"Volcano plot: {representative_celltype} (COVID-19 vs Healthy)"
    )

    heatmap.create(
        representative_subset,
        representative_top_genes,
        groupby="condition",
        output_file=os.path.join(
            OUT_DIR,
            f"heatmap_{safe_name}_covid_vs_healthy.png"
        )
    )

# --------------------------------------------------
# Save report
# --------------------------------------------------
report_lines = []
report_lines.append("========== OMNISEQAI REAL CONDITION ANALYSIS ==========\n")
report_lines.append(f"Cells: {adata.n_obs:,}")
report_lines.append(f"Genes: {adata.n_vars:,}\n")

report_lines.append("Condition counts:")
report_lines.append(adata.obs["condition"].value_counts().to_string())
report_lines.append("\nCell type proportions by condition:")
report_lines.append(proportions.to_string())

for cell_type, result in results_by_celltype.items():
    report_lines.append("\n" + "=" * 70)
    report_lines.append(f"Cell type: {cell_type}")
    report_lines.append(f"Condition counts: {result['condition_counts']}")

    report_lines.append("Top DE genes:")
    if result["sig_genes"] is not None and not result["sig_genes"].empty:
        report_lines.append(
            result["sig_genes"][
                ["names", "logfoldchanges", "pvals_adj"]
            ].head(10).to_string(index=False)
        )
    else:
        report_lines.append("No significant DE genes detected.")

    report_lines.append("Pathway interpretation:")
    for item in result["interpretation"]:
        report_lines.append(f"- {item}")

with open(REPORT_FILE, "w") as f:
    f.write("\n".join(report_lines))

print(f"\nReport saved: {REPORT_FILE}")

# --------------------------------------------------
# Console summary
# --------------------------------------------------
print("\n===== CONDITION-AWARE SUMMARY =====")
for cell_type, result in list(results_by_celltype.items())[:10]:
    print("\n" + "-" * 60)
    print(f"Cell type: {cell_type}")
    print(f"Counts: {result['condition_counts']}")
    print("Interpretation:")
    for item in result["interpretation"]:
        print(f"- {item}")

print("\nAnalysis complete.")
