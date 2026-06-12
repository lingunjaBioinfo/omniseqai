from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Standardize an h5ad file for OmniSeqAI."
    )

    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output standardized .h5ad file")

    parser.add_argument("--condition-col", default=None)
    parser.add_argument("--celltype-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--patient-col", default=None)
    parser.add_argument("--batch-col", default=None)
    parser.add_argument("--gene-symbol-col", default=None)

    parser.add_argument(
        "--compute-umap",
        action="store_true",
        help="Compute UMAP if missing"
    )

    return parser.parse_args()


def pick_column(adata, explicit, candidates):
    if explicit is not None:
        if explicit not in adata.obs.columns and explicit not in adata.var.columns:
            raise ValueError(f"Column not found: {explicit}")
        return explicit

    for col in candidates:
        if col in adata.obs.columns:
            return col

    return None


def clean_condition(x):
    x = str(x).strip()

    mapping = {
        "ctrl": "Control",
        "CTRL": "Control",
        "control": "Control",
        "Control": "Control",
        "normal": "Control",
        "Normal": "Control",
        "healthy": "Control",
        "Healthy": "Control",

        "stim": "Stimulated",
        "STIM": "Stimulated",
        "stimulated": "Stimulated",
        "Stimulated": "Stimulated",

        "tumor": "Tumor",
        "Tumor": "Tumor",
        "cancer": "Tumor",
        "Cancer": "Tumor",

        "covid": "COVID",
        "COVID": "COVID",
        "COVID-19": "COVID",

        "IFN": "IFN_beta",
        "IFNB": "IFN_beta",
        "IFN-beta": "IFN_beta",
        "IFN_beta": "IFN_beta",
    }

    return mapping.get(x, x)


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"\nLoading: {input_path}")
    adata = sc.read_h5ad(input_path)

    print("\nOriginal dataset:")
    print(adata)

    print("\nOriginal obs columns:")
    print(list(adata.obs.columns))

    print("\nOriginal var columns:")
    print(list(adata.var.columns))

    # -----------------------------
    # Condition
    # -----------------------------
    condition_col = pick_column(
        adata,
        args.condition_col,
        [
            "condition",
            "label",
            "stim",
            "treatment",
            "group",
            "status",
            "disease",
            "phenotype",
            "case_control",
        ],
    )

    if condition_col is not None:
        print(f"\nUsing condition column: {condition_col}")
        adata.obs["condition"] = (
            adata.obs[condition_col]
            .astype(str)
            .map(clean_condition)
        )
    else:
        print("\nNo condition column found. Dataset will run in exploratory mode.")

    # -----------------------------
    # Cell type
    # -----------------------------
    celltype_col = pick_column(
        adata,
        args.celltype_col,
        [
            "cell_type",
            "celltype",
            "cell_type_label",
            "seurat_annotations",
            "seurat_annotation",
            "annotation",
            "cluster",
            "cell",
            "louvain",
            "leiden",
            "seurat_clusters",
        ],
    )

    if celltype_col is not None:
        print(f"Using cell type column: {celltype_col}")
        adata.obs["cell_type"] = adata.obs[celltype_col].astype(str)
    else:
        print("No cell type column found. Creating placeholder.")
        adata.obs["cell_type"] = "Unknown"

    # -----------------------------
    # Patient / donor
    # -----------------------------
    patient_col = pick_column(
        adata,
        args.patient_col,
        [
            "patient_id",
            "donor_id",
            "donor",
            "replicate",
            "ind",
            "individual",
            "subject",
        ],
    )

    if patient_col is not None:
        print(f"Using patient column: {patient_col}")
        adata.obs["patient_id"] = adata.obs[patient_col].astype(str)
    else:
        print("No patient column found. Creating patient_1.")
        adata.obs["patient_id"] = "patient_1"

    # -----------------------------
    # Sample ID
    # -----------------------------
    sample_col = pick_column(
        adata,
        args.sample_col,
        [
            "sample_id",
            "sample",
            "orig.ident",
            "library",
            "batch",
        ],
    )

    if sample_col is not None:
        print(f"Using sample column: {sample_col}")
        adata.obs["sample_id"] = adata.obs[sample_col].astype(str)
    else:
        print("No sample column found. Creating sample_id from condition + patient.")
        if "condition" in adata.obs.columns:
            adata.obs["sample_id"] = (
                adata.obs["condition"].astype(str)
                + "_"
                + adata.obs["patient_id"].astype(str)
            )
        else:
            adata.obs["sample_id"] = adata.obs["patient_id"].astype(str)

    # Critical safety:
    # If sample_id does not include condition, pseudobulk can mix groups.
    if "condition" in adata.obs.columns:
        mixed = (
            adata.obs.groupby("sample_id")["condition"]
            .nunique()
            .max()
        )

        if mixed > 1:
            print(
                "\nWARNING: Some sample_id values contain multiple conditions."
            )
            print(
                "Fixing sample_id by combining condition + sample_id."
            )
            adata.obs["sample_id"] = (
                adata.obs["condition"].astype(str)
                + "_"
                + adata.obs["sample_id"].astype(str)
            )

    # -----------------------------
    # Batch
    # -----------------------------
    batch_col = pick_column(
        adata,
        args.batch_col,
        [
            "batch_id",
            "batch",
            "library",
            "library_prep_batch",
        ],
    )

    if batch_col is not None:
        print(f"Using batch column: {batch_col}")
        adata.obs["batch_id"] = adata.obs[batch_col].astype(str)
    else:
        adata.obs["batch_id"] = "batch_1"

    # -----------------------------
    # Gene symbols
    # -----------------------------
    if args.gene_symbol_col is not None:
        if args.gene_symbol_col not in adata.var.columns:
            raise ValueError(f"Gene symbol column not found: {args.gene_symbol_col}")

        print(f"Using gene symbol column: {args.gene_symbol_col}")
        adata.var["gene_symbol"] = adata.var[args.gene_symbol_col].astype(str)

    elif "gene_symbol" in adata.var.columns:
        print("Using existing gene_symbol column.")

    elif "feature_name" in adata.var.columns:
        print("Using feature_name as gene symbols.")
        adata.var["gene_symbol"] = adata.var["feature_name"].astype(str)

    elif "name" in adata.var.columns:
        print("Using name as gene symbols.")
        adata.var["gene_symbol"] = adata.var["name"].astype(str)

    else:
        print("Using var_names as gene symbols.")
        adata.var["gene_symbol"] = adata.var_names.astype(str)

    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()

    # -----------------------------
    # Preserve raw/counts
    # -----------------------------
    if "counts" not in adata.layers:
        try:
            adata.layers["counts"] = adata.X.copy()
        except Exception:
            pass

    if adata.raw is None:
        try:
            adata.raw = adata.copy()
        except Exception:
            pass

    # -----------------------------
    # UMAP
    # -----------------------------
    if args.compute_umap and "X_umap" not in adata.obsm:
        print("\nComputing UMAP...")

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=3000,
            subset=False,
            flavor="seurat",
        )

        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        sc.tl.umap(adata)

    # -----------------------------
    # Summary
    # -----------------------------
    print("\nStandardized dataset:")
    print(adata)

    if "condition" in adata.obs.columns:
        print("\nCondition counts:")
        print(adata.obs["condition"].value_counts())

    print("\nCell type counts:")
    print(adata.obs["cell_type"].value_counts().head(30))

    print("\nSample counts:")
    print(adata.obs["sample_id"].value_counts().head(30))

    print("\nPatient counts:")
    print(adata.obs["patient_id"].value_counts().head(30))

    adata.write_h5ad(output_path)
    print(f"\nSaved standardized dataset: {output_path}")


if __name__ == "__main__":
    main()
