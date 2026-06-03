import scanpy as sc

from backend.analysis_router import AnalysisRouter


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"

print("\nLoading dataset...")
adata = sc.read_h5ad(DATA_FILE)

router = AnalysisRouter()
profile = router.inspect(adata)

print("\n===== DETECTED PROFILE =====")
for k, v in profile.items():
    print(f"{k}: {v}")

results = router.run(adata, profile=profile)

print("\n===== ANALYSIS PLAN =====")
plan = results["plan"]
print(f"has_conditions: {plan.has_conditions}")
print(f"condition_column: {plan.condition_column}")
print(f"cell_type_column: {plan.cell_type_column}")
print(f"condition_groups: {plan.condition_groups}")
print(f"pairwise_comparisons: {plan.pairwise_comparisons}")

print("\n===== CELL-TYPE COMPARISON =====")
if results["celltype_comparison"] is not None:
    print(results["celltype_comparison"])
else:
    print("No comparison produced.")

print("\n===== CONDITION DE INTERPRETATION =====")
if "condition_de_interpretation" in results:
    for item in results["condition_de_interpretation"]:
        print(f"- {item}")
else:
    print("No condition DE interpretation available.")

print("\n===== TOP CELLTYPE-SPECIFIC RESULTS =====")
for cell_type, info in list(results["celltype_specific"].items())[:6]:
    print("\n" + "-" * 60)
    print(f"Cell type: {cell_type}")
    if "error" in info:
        print(f"Error: {info['error']}")
        continue

    print(f"Condition pair: {info['condition_pair']}")
    print(f"Condition counts: {info['condition_counts']}")
    print("Interpretation:")
    for item in info["interpretation"]:
        print(f"- {item}")

print("\nDone.")
