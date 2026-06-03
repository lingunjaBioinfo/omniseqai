import re
import pandas as pd

from backend.pathway_analysis import PathwayAnalyzer


class ConditionPathways:

    IMMUNE_HINTS = [
        "interferon",
        "antiviral",
        "cytokine",
        "chemokine",
        "inflammatory",
        "t cell",
        "b cell",
        "antigen presentation",
        "mhc",
        "innate immune",
        "adaptive immune",
        "nk cell",
        "leukocyte activation",
        "response to virus",
        "viral",
    ]

    HOUSEKEEPING_HINTS = [
        "translation",
        "ribosome",
        "spliceosome",
        "protein-containing complex assembly",
        "cellular respiration",
        "oxidative phosphorylation",
        "mrna splicing",
        "gene expression",
    ]

    def _is_housekeeping(self, term: str) -> bool:
        t = str(term).lower()
        return any(h in t for h in self.HOUSEKEEPING_HINTS)

    def _is_immune(self, term: str) -> bool:
        t = str(term).lower()
        return any(h in t for h in self.IMMUNE_HINTS)

    def analyze(
        self,
        de_results,
        top_n=100,
        min_logfc=0.25,
        max_terms=10
    ):

        if de_results is None or de_results.empty:
            return pd.DataFrame()

        df = de_results.copy()
        df["names"] = df["names"].astype(str).str.strip()

        for col in ["logfoldchanges", "pvals_adj", "scores"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # First pass: strongly significant, positive, non-housekeeping genes.
        filtered = df[
            (df["pvals_adj"] < 0.05) &
            (df["logfoldchanges"] > min_logfc) &
            ~df["names"].str.startswith(("RPL", "RPS", "MT-")) &
            ~df["names"].isin({"MALAT1", "NEAT1", "XIST"})
        ].copy()

        genes = filtered["names"].head(top_n).tolist()

        # Fallback if strict filter is too small.
        if len(genes) < 5:
            genes = (
                df.sort_values(
                    ["pvals_adj", "logfoldchanges"],
                    ascending=[True, False],
                    kind="mergesort"
                )["names"]
                .head(top_n)
                .tolist()
            )

        if len(genes) == 0:
            return pd.DataFrame()

        pathways = PathwayAnalyzer().top_pathways(genes)

        if pathways is None or pathways.empty:
            return pd.DataFrame()

        pathways = pathways.copy()
        pathways["Term"] = pathways["Term"].astype(str)

        # Rank immune terms above housekeeping terms.
        pathways["_immune"] = pathways["Term"].map(self._is_immune)
        pathways["_housekeeping"] = pathways["Term"].map(self._is_housekeeping)

        pathways = pathways.sort_values(
            by=["_immune", "_housekeeping", "Adjusted P-value"],
            ascending=[False, True, True],
            kind="mergesort"
        )

        # If we have immune pathways, keep them first.
        immune_only = pathways[pathways["_immune"] == True].drop(columns=["_immune", "_housekeeping"])
        if not immune_only.empty:
            return immune_only.head(max_terms).reset_index(drop=True)

        # Otherwise return the best available terms, but remove very generic housekeeping terms if possible.
        non_housekeeping = pathways[pathways["_housekeeping"] == False].drop(columns=["_immune", "_housekeeping"])
        if not non_housekeeping.empty:
            return non_housekeeping.head(max_terms).reset_index(drop=True)

        return pathways.drop(columns=["_immune", "_housekeeping"]).head(max_terms).reset_index(drop=True)
