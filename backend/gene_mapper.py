from __future__ import annotations


class GeneMapper:

    def fix_gene_names(self, adata):
        print("\nChecking gene identifiers...")

        # If the object already has sensible gene symbols in var_names,
        # do not touch them.
        if len(adata.var_names) > 0:
            first = str(adata.var_names[0])
            if (
                not first.startswith("ENSG")
                and "feature_name" not in adata.var.columns
                and "gene_symbol" not in adata.var.columns
            ):
                print("Gene symbols already detected.")
                return adata

        # Preferred source: feature_name
        if "feature_name" in adata.var.columns:
            print("Using feature_name as gene symbols.")

            symbols = adata.var["feature_name"].astype(str)

            # Replace blanks / missing / ENSG-like entries with gene_symbol if present
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

            # Final cleanup: keep original var_names only if still missing
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

        # Secondary source: gene_symbol
        if "gene_symbol" in adata.var.columns:
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

        print("WARNING: No gene symbol column found. Keeping original var_names.")
        return adata
