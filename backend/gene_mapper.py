import scanpy as sc


class GeneMapper:

    def fix_gene_names(self, adata):
        print("\nChecking gene identifiers...")

        if "feature_name" in adata.var.columns:
            print("Using feature_name as gene symbols.")

            symbols = adata.var["feature_name"].astype(str)

            # fall back to the original name when feature_name is missing/blank
            fallback = adata.var_names.astype(str)
            symbols = symbols.where(
                symbols.notna() & (symbols != "") & (symbols != "nan"),
                fallback
            )

            adata.var["gene_symbol"] = symbols
            adata.var_names = adata.var["gene_symbol"]
            adata.var_names_make_unique()
            return adata

        if len(adata.var_names) > 0 and not str(adata.var_names[0]).startswith("ENSG"):
            print("Gene symbols already detected.")
            return adata

        print("WARNING: Ensembl IDs detected but no feature_name column found.")
        return adata
