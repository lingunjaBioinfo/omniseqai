import os
import pandas as pd

from backend.pathway_analysis import PathwayAnalyzer


class ConditionPathways:

    def analyze(
        self,
        de_results,
        top_n=100,
        min_logfc=0.25,
        max_terms=10
    ):
        if os.getenv("OMNISEQAI_ENABLE_PATHWAYS", "0") != "1":
            return pd.DataFrame()

        if de_results is None or de_results.empty:
            return pd.DataFrame()

        df = de_results.copy()
        df["names"] = df["names"].astype(str).str.strip()

        for col in ["logfoldchanges", "pvals_adj", "scores"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        filtered = df[
            (df["pvals_adj"] < 0.05) &
            (df["logfoldchanges"] > min_logfc) &
            ~df["names"].str.startswith(("RPL", "RPS", "MT-")) &
            ~df["names"].isin({"MALAT1", "NEAT1", "XIST"})
        ].copy()

        genes = filtered["names"].head(top_n).tolist()

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

        pathways = PathwayAnalyzer().top_pathways(genes, organism="human")

        if pathways is None or pathways.empty:
            return pd.DataFrame()

        if "Term" not in pathways.columns:
            return pd.DataFrame()

        return pathways.head(max_terms).reset_index(drop=True)
