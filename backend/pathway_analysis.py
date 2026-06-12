from __future__ import annotations

import os
from typing import Iterable, List, Optional

import pandas as pd


class PathwayAnalyzer:
    """
    Safe pathway enrichment wrapper.

    By default, pathway enrichment is disabled so OmniSeqAI can finish
    without waiting on remote Enrichr requests.
    Enable with:
        OMNISEQAI_ENABLE_PATHWAYS=1
    """

    def __init__(
        self,
        gene_sets: Optional[List[str]] = None,
        organism: str = "human",
    ):
        self.gene_sets = gene_sets or [
            "GO_Biological_Process_2021",
            "KEGG_2021_Human",
        ]
        self.organism = organism

    def _clean_genes(self, genes: Iterable[str]) -> List[str]:
        cleaned: List[str] = []

        for g in genes:
            if g is None:
                continue
            s = str(g).strip()
            if not s:
                continue
            if s.startswith("ENSG"):
                continue
            if s.startswith("NCBITaxon:"):
                continue
            if s not in cleaned:
                cleaned.append(s)

        return cleaned

    def enrich_markers(
        self,
        genes: Iterable[str],
        gene_sets: Optional[List[str]] = None,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        if os.getenv("OMNISEQAI_ENABLE_PATHWAYS", "0") != "1":
            return pd.DataFrame()

        genes = self._clean_genes(genes)
        if len(genes) == 0:
            return pd.DataFrame()

        try:
            import gseapy as gp
        except Exception as e:
            print(f"Pathway enrichment unavailable (gseapy import failed): {e}")
            return pd.DataFrame()

        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_sets or self.gene_sets,
                organism=organism or self.organism,
                outdir=None,
                cutoff=0.5,
            )

            if enr is None or not hasattr(enr, "results") or enr.results is None:
                return pd.DataFrame()

            df = enr.results.copy()
            if df.empty:
                return df

            return df

        except Exception as e:
            print(f"Pathway enrichment failed: {e}")
            return pd.DataFrame()

    def top_pathways(
        self,
        genes: Iterable[str],
        top_n: int = 10,
        gene_sets: Optional[List[str]] = None,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        df = self.enrich_markers(
            genes=genes,
            gene_sets=gene_sets,
            organism=organism,
        )

        if df is None or df.empty:
            return pd.DataFrame()

        if "Term" not in df.columns:
            return pd.DataFrame()

        df = df.copy()

        if "Adjusted P-value" in df.columns:
            df["Adjusted P-value"] = pd.to_numeric(df["Adjusted P-value"], errors="coerce")
            df = df.sort_values("Adjusted P-value", ascending=True, kind="mergesort")
        elif "P-value" in df.columns:
            df["P-value"] = pd.to_numeric(df["P-value"], errors="coerce")
            df = df.sort_values("P-value", ascending=True, kind="mergesort")

        return df.head(top_n).reset_index(drop=True)
