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

plan = results["plan"]

print("\n===== ANALYSIS PLAN =====")
print(f"has_conditions: {plan.has_conditions}")
print(f"condition_column: {plan.condition_column}")
print(f"cell_type_column: {plan.cell_type_column}")
print(f"condition_groups: {plan.condition_groups}")
print(f"baseline: {plan.baseline}")
print(f"pairwise_comparisons: {plan.pairwise_comparisons}")
print(f"use_pseudobulk: {plan.use_pseudobulk}")
print(f"pseudobulk_key: {plan.pseudobulk_key}")

print("\n===== CELL-TYPE COMPARISON =====")
if results["celltype_comparison"] is not None:
    print(results["celltype_comparison"])
else:
    print("No comparison produced.")

print("\n===== WHOLE-DATASET CONDITION DE =====")
for pair, info in results["condition_de_results"].items():
    print("\n" + "-" * 60)
    print(f"Comparison: {pair}")
    print(f"Mode: {info.get('mode')}")
    print(f"Significant genes: {info.get('n_sig_genes')}")
    print("Interpretation:")
    for item in info.get("interpretation", []):
        print(f"- {item}")

print("\n===== TOP CELL-TYPE-SPECIFIC RESULTS =====")
for cell_type, pair_map in list(results["celltype_specific"].items())[:6]:
    print("\n" + "-" * 60)
    print(f"Cell type: {cell_type}")

    for pair, info in pair_map.items():
        print(f"  Pair: {pair}")
        print(f"  Mode: {info.get('mode')}")
        if "error" in info and info["error"]:
            print(f"  Error: {info['error']}")
            continue

        print(f"  Significant genes: {info.get('n_sig_genes')}")
        print("  Interpretation:")
        for item in info.get("interpretation", []):
            print(f"  - {item}")

print("\nDone.")
