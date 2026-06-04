import scanpy as sc

from backend.metadata_detector import MetadataDetector
from backend.pseudobulk_de import PseudobulkDE
from backend.condition_pathways import ConditionPathways
from backend.condition_interpreter import ConditionInterpreter


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"

print("\nLoading dataset...")
adata = sc.read_h5ad(DATA_FILE)

detector = MetadataDetector()
profile = detector.detect(adata)
adata = detector.apply(adata, profile)

print("\n===== PROFILE =====")
for k, v in profile.items():
    print(f"{k}: {v}")

if "condition" not in adata.obs.columns:
    raise ValueError("No standardized 'condition' column found.")

if "sample_id" not in adata.obs.columns:
    raise ValueError("No standardized 'sample_id' column found.")

print("\nCondition counts:")
print(adata.obs["condition"].value_counts())

print("\nSample counts:")
print(adata.obs["sample_id"].value_counts().head(10))

pb_de = PseudobulkDE()

# Use Healthy vs COVID if those labels exist after normalization.
# The metadata detector usually standardizes normal/COVID-19 to Healthy/COVID.
group1 = "COVID"
group2 = "Healthy"

result = pb_de.run(
    adata,
    sample_key="sample_id",
    condition_key="condition",
    group1=group1,
    group2=group2,
    min_cells_per_sample=10,
)

pb = result.pseudobulk_adata

print("\n===== PSEUDOBULK SUMMARY =====")
print(pb)

print("\nPseudobulk condition counts:")
print(pb.obs["condition"].value_counts())

de = pb_de.top_genes(
    pb,
    group=group1,
    n_genes=30
)

print("\n===== TOP PSEUDOBULK DE GENES =====")
print(
    de[
        ["names", "logfoldchanges", "pvals_adj"]
    ]
)

pathways = ConditionPathways().analyze(
    de,
    top_n=50
)

print("\n===== PSEUDOBULK PATHWAYS =====")
if pathways is not None and not pathways.empty:
    print(pathways[["Term", "Adjusted P-value"]].head(10))
    interpretation = ConditionInterpreter().interpret(pathways)
    print("\n===== INTERPRETATION =====")
    for item in interpretation:
        print(f"- {item}")
else:
    print("No pathways returned.")

print("\nDone.")
