from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


OUTDIR = Path("data/external/new_dataset_test")
FULL_PATH = OUTDIR / "pbmc68k_or_synthetic_test.h5ad"
MULTI_DIR = OUTDIR / "multisample"


def load_base_dataset():
    """
    Load a genuinely new dataset for testing.

    Preferred:
    - pbmc68k_reduced from Scanpy

    Fallback:
    - synthetic blobs dataset if pbmc68k_reduced is unavailable
    """

    try:
        print("Trying scanpy.datasets.pbmc68k_reduced()...")
        adata = sc.datasets.pbmc68k_reduced()
        dataset_name = "pbmc68k_reduced"
        print("Loaded pbmc68k_reduced.")
        return adata.copy(), dataset_name

    except Exception as e:
        print(f"Could not load pbmc68k_reduced: {e}")
        print("Falling back to scanpy.datasets.blobs().")

        adata = sc.datasets.blobs(
            n_variables=1200,
            n_observations=3000,
            n_centers=8,
            cluster_std=1.2,
            random_state=7,
        )

        dataset_name = "synthetic_blobs"
        print("Loaded synthetic blobs dataset.")
        return adata.copy(), dataset_name


def choose_cell_type_column(adata):
    candidates = [
        "cell_type",
        "bulk_labels",
        "louvain",
        "leiden",
        "clusters",
        "cluster",
        "blobs",
    ]

    for col in candidates:
        if col in adata.obs.columns:
            return col

    return None


def add_test_metadata(adata, dataset_name: str):
    """
    Add controlled sample/condition/batch metadata.

    This gives OmniSeqAI enough structure to test:
    - pseudobulk
    - multi-sample loading
    - batch-aware QC
    - condition analysis
    """

    rng = np.random.default_rng(7)

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if "feature_name" not in adata.var.columns:
        adata.var["feature_name"] = adata.var_names.astype(str)

    cell_type_col = choose_cell_type_column(adata)

    if cell_type_col is not None:
        adata.obs["cell_type"] = adata.obs[cell_type_col].astype(str)

    else:
        adata.obs["cell_type"] = "Unknown"

    samples = np.array(
        [
            "Healthy_1",
            "Healthy_2",
            "Healthy_3",
            "Disease_1",
            "Disease_2",
            "Disease_3",
        ]
    )

    sample_assignment = np.resize(samples, adata.n_obs)
    rng.shuffle(sample_assignment)

    adata.obs["sample_id"] = sample_assignment
    adata.obs["patient_id"] = sample_assignment

    adata.obs["condition"] = [
        "Healthy" if str(sample).startswith("Healthy") else "Disease"
        for sample in sample_assignment
    ]

    adata.obs["condition_raw"] = adata.obs["condition"].astype(str)

    batch_map = {
        "Healthy_1": "batch_1",
        "Disease_1": "batch_1",
        "Healthy_2": "batch_2",
        "Disease_2": "batch_2",
        "Healthy_3": "batch_3",
        "Disease_3": "batch_3",
    }

    adata.obs["batch_id"] = [
        batch_map.get(str(sample), "batch_unknown")
        for sample in sample_assignment
    ]

    adata.obs["source_dataset"] = dataset_name

    return adata


def write_outputs(adata):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    MULTI_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Writing full dataset: {FULL_PATH}")
    adata.write_h5ad(FULL_PATH)

    print(f"Writing multi-sample folder: {MULTI_DIR}")

    for sample_id in sorted(adata.obs["sample_id"].astype(str).unique()):
        sample_dir = MULTI_DIR / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        sample_adata = adata[adata.obs["sample_id"].astype(str) == sample_id].copy()

        sample_path = sample_dir / f"{sample_id}.h5ad"
        sample_adata.write_h5ad(sample_path)

        condition = sample_adata.obs["condition"].astype(str).iloc[0]
        batch = sample_adata.obs["batch_id"].astype(str).iloc[0]

        print(
            f"{sample_id}: {sample_adata.n_obs:,} cells, "
            f"condition={condition}, batch={batch} -> {sample_path}"
        )


def print_summary(adata):
    print("\nDataset summary")
    print("---------------")
    print(adata)
    print("\nCell types:")
    print(adata.obs["cell_type"].value_counts().head(20))

    print("\nSamples:")
    print(adata.obs["sample_id"].value_counts())

    print("\nConditions:")
    print(adata.obs["condition"].value_counts())

    print("\nBatches:")
    print(adata.obs["batch_id"].value_counts())

    print("\nSamples per condition:")
    print(
        adata.obs
        .assign(
            sample_id=adata.obs["sample_id"].astype(str),
            condition=adata.obs["condition"].astype(str),
        )
        .groupby("condition")["sample_id"]
        .nunique()
    )


def main():
    adata, dataset_name = load_base_dataset()
    adata = add_test_metadata(adata, dataset_name)
    write_outputs(adata)
    print_summary(adata)

    print("\nDone.")
    print(f"Full file: {FULL_PATH}")
    print(f"Multi-sample folder: {MULTI_DIR}")


if __name__ == "__main__":
    main()
