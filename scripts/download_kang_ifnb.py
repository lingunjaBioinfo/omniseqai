from pathlib import Path

import scanpy as sc


def pick_column(adata, candidates):
    for col in candidates:
        if col in adata.obs.columns:
            return col
    return None


def clean_condition(value):
    value = str(value).strip()

    mapping = {
        "ctrl": "Control",
        "CTRL": "Control",
        "control": "Control",
        "Control": "Control",
        "CTRL_": "Control",
        "stim": "IFN_beta",
        "STIM": "IFN_beta",
        "stimulated": "IFN_beta",
        "Stimulated": "IFN_beta",
        "stim_": "IFN_beta",
        "IFN": "IFN_beta",
        "IFNB": "IFN_beta",
        "IFN-beta": "IFN_beta",
        "IFN_beta": "IFN_beta",
    }

    return mapping.get(value, value)


def main():
    import pertpy as pt

    outdir = Path("data/kang_ifnb")
    outdir.mkdir(parents=True, exist_ok=True)

    print("\nDownloading Kang 2018 IFN-beta PBMC dataset...")
    adata = pt.dt.kang_2018()

    print("\nOriginal dataset:")
    print(adata)

    print("\nOriginal obs columns:")
    print(list(adata.obs.columns))

    # --------------------------------------------------
    # Preserve original data where possible
    # --------------------------------------------------
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

    # --------------------------------------------------
    # Standardize condition metadata
    # --------------------------------------------------
    condition_col = pick_column(
        adata,
        [
            "condition",
            "label",
            "stim",
            "treatment",
            "group",
            "status",
            "disease",
        ],
    )

    if condition_col is None:
        raise ValueError(
            "No condition column found. Available columns: "
            f"{list(adata.obs.columns)}"
        )

    print(f"\nUsing condition column: {condition_col}")
    print("\nRaw condition values:")
    print(adata.obs[condition_col].value_counts())

    adata.obs["condition"] = (
        adata.obs[condition_col]
        .astype(str)
        .map(clean_condition)
    )

    # --------------------------------------------------
    # Standardize cell type metadata
    # --------------------------------------------------
    celltype_col = pick_column(
        adata,
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
        print(f"\nUsing cell type column: {celltype_col}")
        adata.obs["cell_type"] = adata.obs[celltype_col].astype(str)
    else:
        print("\nNo cell type column found. Creating placeholder labels.")
        adata.obs["cell_type"] = "Unknown"

    # --------------------------------------------------
    # Standardize donor / replicate metadata
    # --------------------------------------------------
    donor_col = pick_column(
        adata,
        [
            "donor",
            "donor_id",
            "replicate",
            "ind",
            "individual",
            "patient_id",
        ],
    )

    if donor_col is not None:
        print(f"\nUsing donor/patient column: {donor_col}")
        adata.obs["patient_id"] = adata.obs[donor_col].astype(str)
    else:
        print("\nNo donor/patient column found. Creating placeholder patient IDs.")
        adata.obs["patient_id"] = "patient_1"

    # Important:
    # For Kang, each donor/replicate can have both control and stimulated cells.
    # Therefore sample_id must include BOTH donor and condition.
    # Otherwise pseudobulk would accidentally mix control and IFN_beta cells.
    adata.obs["sample_id"] = (
        adata.obs["condition"].astype(str)
        + "_"
        + adata.obs["patient_id"].astype(str)
    )

    # --------------------------------------------------
    # Standardize batch metadata
    # --------------------------------------------------
    batch_col = pick_column(
        adata,
        [
            "batch_id",
            "batch",
            "library",
            "library_prep_batch",
        ],
    )

    if batch_col is not None:
        print(f"\nUsing batch column: {batch_col}")
        adata.obs["batch_id"] = adata.obs[batch_col].astype(str)
    else:
        adata.obs["batch_id"] = "batch_1"

    # --------------------------------------------------
    # Standardize gene symbols
    # --------------------------------------------------
    if "gene_symbol" not in adata.var.columns:
        if "name" in adata.var.columns:
            print("\nUsing adata.var['name'] as gene symbols.")
            adata.var["gene_symbol"] = adata.var["name"].astype(str)
        else:
            print("\nUsing adata.var_names as gene symbols.")
            adata.var["gene_symbol"] = adata.var_names.astype(str)

    # Keep var_names usable for OmniSeqAI
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()

    # --------------------------------------------------
    # Use existing UMAP if present; compute only if missing
    # --------------------------------------------------
    if "X_umap" not in adata.obsm:
        print("\nNo UMAP found. Computing PCA, neighbors, and UMAP...")

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
    else:
        print("\nExisting UMAP detected.")

    # --------------------------------------------------
    # Final report
    # --------------------------------------------------
    print("\nStandardized dataset:")
    print(adata)

    print("\nCondition counts:")
    print(adata.obs["condition"].value_counts())

    print("\nCell type counts:")
    print(adata.obs["cell_type"].value_counts().head(30))

    print("\nPatient counts:")
    print(adata.obs["patient_id"].value_counts().head(20))

    print("\nSample counts:")
    print(adata.obs["sample_id"].value_counts().head(20))

    out = outdir / "kang_ifnb.h5ad"
    adata.write_h5ad(out)

    print(f"\nSaved: {out}")
    print("\nDone.")


if __name__ == "__main__":
    main()
