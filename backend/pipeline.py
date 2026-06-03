import scanpy as sc


class SingleCellPipeline:

    def preprocess(
        self,
        adata
    ):

        print("\nPreprocessing dataset...")

        if "total_counts" not in adata.obs:

            sc.pp.calculate_qc_metrics(
                adata,
                inplace=True
            )

        max_val = adata.X.max()

        if max_val > 50:

            print(
                "Raw counts detected."
            )

            sc.pp.normalize_total(
                adata,
                target_sum=10000
            )

            sc.pp.log1p(adata)

        else:

            print(
                "Data already normalized."
            )

        return adata

    # --------------------------
    # DATA LOADING
    # --------------------------

    def load_pbmc3k(self):

        self.adata = sc.datasets.pbmc3k()

        return self.adata

    def summary(self):

        if self.adata is None:
            raise ValueError("No dataset loaded")

        print(self.adata)

        print(f"\nCells: {self.adata.n_obs}")
        print(f"Genes: {self.adata.n_vars}")

    # --------------------------
    # QUALITY CONTROL
    # --------------------------

    def calculate_qc(self):

        self.adata.var["mt"] = (
            self.adata.var_names.str.startswith("MT-")
        )

        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=["mt"],
            inplace=True
        )

        return self.adata

    def qc_summary(self):

        print("\n===== QC SUMMARY =====")

        print(
            f"Mean genes per cell: "
            f"{self.adata.obs['n_genes_by_counts'].mean():.2f}"
        )

        print(
            f"Mean counts per cell: "
            f"{self.adata.obs['total_counts'].mean():.2f}"
        )

        print(
            f"Mean mitochondrial %: "
            f"{self.adata.obs['pct_counts_mt'].mean():.2f}"
        )

    # --------------------------
    # FILTERING
    # --------------------------

    def filter_data(
        self,
        min_genes=200,
        max_genes=2500,
        max_mt=5
    ):

        self.raw_cells = self.adata.n_obs

        self.adata = self.adata[
            self.adata.obs["n_genes_by_counts"] > min_genes
        ]

        self.adata = self.adata[
            self.adata.obs["n_genes_by_counts"] < max_genes
        ]

        self.adata = self.adata[
            self.adata.obs["pct_counts_mt"] < max_mt
        ].copy()

        self.filtered_cells = self.adata.n_obs

        self.removed_cells = (
            self.raw_cells -
            self.filtered_cells
        )

        return self.adata

    def filter_summary(self):

        print("\n===== FILTER SUMMARY =====")

        print(
            f"Input cells: "
            f"{self.raw_cells}"
        )

        print(
            f"Retained cells: "
            f"{self.filtered_cells}"
        )

        print(
            f"Removed cells: "
            f"{self.removed_cells}"
        )

    # --------------------------
    # NORMALIZATION
    # --------------------------

    def normalize_data(self):

        sc.pp.normalize_total(
            self.adata,
            target_sum=10000
        )

        sc.pp.log1p(self.adata)

        return self.adata

    # --------------------------
    # HIGHLY VARIABLE GENES
    # --------------------------

    def identify_hvg(
        self,
        n_top_genes=2000
    ):

        sc.pp.highly_variable_genes(
            self.adata,
            n_top_genes=n_top_genes
        )

        hvg_count = (
            self.adata.var["highly_variable"]
            .sum()
        )

        print(
            f"\nHighly Variable Genes: "
            f"{hvg_count}"
        )

        return self.adata

    # --------------------------
    # PCA
    # --------------------------

    def run_pca(self):

        sc.tl.pca(
            self.adata,
            svd_solver="arpack"
        )

        print("\nPCA completed.")

        return self.adata

    # --------------------------
    # NEIGHBORS
    # --------------------------

    def compute_neighbors(self):

        sc.pp.neighbors(
            self.adata,
            n_neighbors=10,
            n_pcs=40
        )

        print("\nNeighbor graph computed.")

        return self.adata

    # --------------------------
    # UMAP
    # --------------------------

    def run_umap(self):

        sc.tl.umap(self.adata)

        print("\nUMAP completed.")

        return self.adata

    # --------------------------
    # CLUSTERING
    # --------------------------

    def cluster_cells(
        self,
        resolution=0.5
    ):

        sc.tl.leiden(
            self.adata,
            resolution=resolution
        )

        n_clusters = (
            self.adata.obs["leiden"]
            .nunique()
        )

        print(
            f"\nClusters identified: "
            f"{n_clusters}"
        )

        return self.adata
