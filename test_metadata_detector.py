import scanpy as sc

from backend.metadata_detector import MetadataDetector


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"

print("\nLoading dataset...")
adata = sc.read_h5ad(DATA_FILE)

detector = MetadataDetector()
profile = detector.detect(adata)

print("\n===== DETECTED PROFILE =====")
for k, v in profile.items():
    print(f"{k}: {v}")

adata = detector.apply(adata, profile)

print("\n===== STANDARDIZED OBS COLUMNS =====")
for col in ["cell_type", "condition", "sample_id", "patient_id", "batch_id"]:
    if col in adata.obs.columns:
        print(f"{col}: present")
    else:
        print(f"{col}: missing")

print("\n===== CONDITION COUNTS =====")
if "condition" in adata.obs.columns:
    print(adata.obs["condition"].value_counts())

print("\n===== CELL TYPE COUNTS =====")
if "cell_type" in adata.obs.columns:
    print(adata.obs["cell_type"].value_counts().head(10))

print("\nDone.")
