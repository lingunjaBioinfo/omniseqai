from __future__ import annotations


class GeneMapper:

    def fix_gene_names(self, adata, verbose: bool = False):
        if verbose:
            print("\nChecking gene identifiers...")

        if len(adata.var_names) > 0:
            first = str(adata.var_names[0])
            if (
                not first.startswith("ENSG")
                and "feature_name" not in adata.var.columns
                and "gene_symbol" not in adata.var.columns
            ):
                if verbose:
                    print("Gene symbols already detected.")
                return adata

        if "feature_name" in adata.var.columns:
            if verbose:
                print("Using feature_name as gene symbols.")

            symbols = adata.var["feature_name"].astype(str)

            if "gene_symbol" in adata.var.columns:
                fallback = adata.var["gene_symbol"].astype(str)
                mask = (
                    symbols.isna()
                    | (symbols == "")
                    | (symbols == "nan")
                    | symbols.str.startswith("ENSG", na=False)
                    | symbols.str.startswith("NCBITaxon:", na=False)
                )
                symbols = symbols.where(~mask, fallback)

            fallback = adata.var_names.astype(str)
            mask = (
                symbols.isna()
                | (symbols == "")
                | (symbols == "nan")
                | symbols.str.startswith("ENSG", na=False)
                | symbols.str.startswith("NCBITaxon:", na=False)
            )
            symbols = symbols.where(~mask, fallback)

            adata.var["gene_symbol"] = symbols.astype(str)
            adata.var_names = adata.var["gene_symbol"]
            adata.var_names_make_unique()
            return adata

        if "gene_symbol" in adata.var.columns:
            if verbose:
                print("Using gene_symbol as gene symbols.")
            symbols = adata.var["gene_symbol"].astype(str)

            mask = (
                symbols.isna()
                | (symbols == "")
                | (symbols == "nan")
                | symbols.str.startswith("ENSG", na=False)
                | symbols.str.startswith("NCBITaxon:", na=False)
            )

            symbols = symbols.where(~mask, adata.var_names.astype(str))

            adata.var["gene_symbol"] = symbols.astype(str)
            adata.var_names = adata.var["gene_symbol"]
            adata.var_names_make_unique()
            return adata

        if verbose:
            print("WARNING: No gene symbol column found. Keeping original var_names.")
        return adata
