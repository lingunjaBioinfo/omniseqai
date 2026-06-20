from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from backend.biology_signatures import BIOLOGY_SIGNATURES, normalize_gene_symbol


class BiologyValidator:
    """
    Validate whether differential-expression results match known biological programs.
    """

    def __init__(
        self,
        padj_cutoff: float = 0.05,
        min_abs_logfc: float = 0.25,
        min_hits: int = 2,
    ):
        self.padj_cutoff = padj_cutoff
        self.min_abs_logfc = min_abs_logfc
        self.min_hits = min_hits

    def run(self, results: Dict[str, Any]) -> Dict[str, Any]:
        condition_de_results = results.get("condition_de_results", {}) or {}
        celltype_specific = results.get("celltype_specific", {}) or {}

        summary_rows: List[Dict[str, Any]] = []
        hit_rows: List[Dict[str, Any]] = []

        # Whole-dataset DE
        for pair, info in condition_de_results.items():
            if not isinstance(info, dict):
                continue

            de_results = info.get("de_results")

            rows_summary, rows_hits = self._score_de_table(
                de_results=de_results,
                comparison=pair,
                scope="whole_dataset",
                cell_type=None,
            )

            summary_rows.extend(rows_summary)
            hit_rows.extend(rows_hits)

        # Cell-type-specific DE
        for cell_type, pair_map in celltype_specific.items():
            if not isinstance(pair_map, dict):
                continue

            for pair, info in pair_map.items():
                if not isinstance(info, dict):
                    continue

                de_results = info.get("de_results")

                rows_summary, rows_hits = self._score_de_table(
                    de_results=de_results,
                    comparison=pair,
                    scope="celltype_specific",
                    cell_type=cell_type,
                )

                summary_rows.extend(rows_summary)
                hit_rows.extend(rows_hits)

        summary_df = pd.DataFrame(summary_rows)
        hits_df = pd.DataFrame(hit_rows)

        interpretation = self._build_interpretation(summary_df)

        status = "ok" if not summary_df.empty else "no_de_results"

        return {
            "status": status,
            "summary": summary_df,
            "hits": hits_df,
            "interpretation": interpretation,
        }

    # --------------------------------------------------
    # Scoring
    # --------------------------------------------------
    def _score_de_table(
        self,
        de_results: Optional[pd.DataFrame],
        comparison,
        scope: str,
        cell_type: Optional[str],
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        if de_results is None:
            return [], []

        if not hasattr(de_results, "empty") or de_results.empty:
            return [], []

        df = de_results.copy()

        gene_col = self._gene_column(df)
        lfc_col = self._lfc_column(df)
        padj_col = self._padj_column(df)

        if gene_col is None or lfc_col is None:
            return [], []

        if padj_col is None:
            df["_padj_for_validation"] = 1.0
            padj_col = "_padj_for_validation"

        df["_gene_norm"] = df[gene_col].map(normalize_gene_symbol)

        case_up = df[
            (df[lfc_col] >= self.min_abs_logfc)
            & (df[padj_col] <= self.padj_cutoff)
        ].copy()

        reference_up = df[
            (df[lfc_col] <= -self.min_abs_logfc)
            & (df[padj_col] <= self.padj_cutoff)
        ].copy()

        comparison_name = self._comparison_name(comparison)

        summary_rows = []
        hit_rows = []

        for signature_name, signature in BIOLOGY_SIGNATURES.items():
            expected_genes = {
                normalize_gene_symbol(g)
                for g in signature.get("genes", [])
            }

            case_hits = case_up[case_up["_gene_norm"].isin(expected_genes)]
            reference_hits = reference_up[reference_up["_gene_norm"].isin(expected_genes)]

            n_expected = len(expected_genes)
            n_case_hits = int(case_hits["_gene_norm"].nunique())
            n_reference_hits = int(reference_hits["_gene_norm"].nunique())

            hit_fraction_case = n_case_hits / n_expected if n_expected else 0.0
            hit_fraction_reference = n_reference_hits / n_expected if n_expected else 0.0

            detected = n_case_hits >= self.min_hits

            summary_rows.append(
                {
                    "scope": scope,
                    "cell_type": cell_type or "whole_dataset",
                    "comparison": comparison_name,
                    "signature": signature_name,
                    "description": signature.get("description"),
                    "expected_direction": signature.get("expected_direction"),
                    "n_expected_genes": n_expected,
                    "n_case_up_hits": n_case_hits,
                    "n_reference_up_hits": n_reference_hits,
                    "case_up_hit_fraction": hit_fraction_case,
                    "reference_up_hit_fraction": hit_fraction_reference,
                    "detected": detected,
                }
            )

            for _, row in case_hits.iterrows():
                hit_rows.append(
                    {
                        "scope": scope,
                        "cell_type": cell_type or "whole_dataset",
                        "comparison": comparison_name,
                        "signature": signature_name,
                        "direction": "case_up",
                        "gene": row.get(gene_col),
                        "logfoldchange": row.get(lfc_col),
                        "pvals_adj": row.get(padj_col),
                    }
                )

            for _, row in reference_hits.iterrows():
                hit_rows.append(
                    {
                        "scope": scope,
                        "cell_type": cell_type or "whole_dataset",
                        "comparison": comparison_name,
                        "signature": signature_name,
                        "direction": "reference_up",
                        "gene": row.get(gene_col),
                        "logfoldchange": row.get(lfc_col),
                        "pvals_adj": row.get(padj_col),
                    }
                )

        return summary_rows, hit_rows

    def _gene_column(self, df: pd.DataFrame) -> Optional[str]:
        for col in ["gene_symbol", "names", "gene", "symbol", "gene_name"]:
            if col in df.columns:
                return col

        return None

    def _lfc_column(self, df: pd.DataFrame) -> Optional[str]:
        for col in ["logfoldchanges", "logFC", "log2FoldChange"]:
            if col in df.columns:
                return col

        return None

    def _padj_column(self, df: pd.DataFrame) -> Optional[str]:
        for col in ["pvals_adj", "padj", "fdr", "qval"]:
            if col in df.columns:
                return col

        return None

    def _comparison_name(self, comparison) -> str:
        if isinstance(comparison, tuple) and len(comparison) == 2:
            return f"{comparison[0]} vs {comparison[1]}"

        return str(comparison)

    # --------------------------------------------------
    # Interpretation
    # --------------------------------------------------
    def _build_interpretation(self, summary_df: pd.DataFrame) -> List[str]:
        if summary_df is None or summary_df.empty:
            return [
                "No differential-expression results were available for biology validation."
            ]

        detected = summary_df[summary_df["detected"] == True].copy()

        if detected.empty:
            return [
                "No predefined biological signature reached the detection threshold."
            ]

        detected = detected.sort_values(
            ["n_case_up_hits", "case_up_hit_fraction"],
            ascending=[False, False],
        )

        top = detected.iloc[0]

        lines = [
            (
                f"Top detected biological program: {top['signature']} "
                f"({top['description']})"
            ),
            (
                f"Evidence: {int(top['n_case_up_hits'])} expected genes were "
                f"significantly upregulated in the case group for "
                f"{top['comparison']}."
            ),
        ]

        # Add secondary programs if present.
        top_signature = str(top["signature"])

        secondary = detected[
            detected["signature"].astype(str) != top_signature
        ].copy()

        if not secondary.empty:
            unique_programs = []

            for sig in secondary["signature"].astype(str).tolist():
                if sig not in unique_programs:
                    unique_programs.append(sig)

            unique_programs = unique_programs[:3]

            if unique_programs:
                programs = ", ".join(unique_programs)
                lines.append(f"Additional detected programs: {programs}.")

        return lines
