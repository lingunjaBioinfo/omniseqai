from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

import math
import numpy as np
import pandas as pd

from backend.biology_signatures import normalize_gene_symbol
from backend.signature_loader import merge_builtin_and_user_signatures


class BiologyValidator:
    """
    Validate differential-expression results against biological gene signatures.

    Supports:
    - built-in curated signatures
    - optional user-defined signatures loaded from CSV/TSV
    - whole-dataset condition DE
    - cell-type-specific DE
    """

    def __init__(
        self,
        signatures: Optional[Dict[str, Dict[str, Any]]] = None,
        user_signature_path: Optional[str] = None,
        padj_cutoff: float = 0.05,
        logfc_cutoff: float = 0.25,
        min_hits: int = 2,
        min_fraction: float = 0.20,
    ):
        if signatures is not None:
            self.signatures = signatures
        else:
            self.signatures = merge_builtin_and_user_signatures(
                user_signature_path=user_signature_path
            )

        self.padj_cutoff = padj_cutoff
        self.logfc_cutoff = logfc_cutoff
        self.min_hits = min_hits
        self.min_fraction = min_fraction

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def run(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run biology validation on an OmniSeqAI result dictionary.
        """

        summary_records: List[Dict[str, Any]] = []
        hit_records: List[Dict[str, Any]] = []

        # --------------------------------------------------------------
        # Whole-dataset condition DE
        # --------------------------------------------------------------
        for comparison, de_results in self._iter_condition_de(results):
            summary, hits = self._score_de_table(
                de_results=de_results,
                scope="whole_dataset",
                cell_type="whole_dataset",
                comparison=comparison,
            )

            summary_records.extend(summary)
            hit_records.extend(hits)

        # --------------------------------------------------------------
        # Cell-type-specific DE
        # --------------------------------------------------------------
        for cell_type, comparison, de_results in self._iter_celltype_de(results):
            summary, hits = self._score_de_table(
                de_results=de_results,
                scope="celltype_specific",
                cell_type=cell_type,
                comparison=comparison,
            )

            summary_records.extend(summary)
            hit_records.extend(hits)

        if not summary_records:
            return {
                "status": "no_de_results",
                "summary": pd.DataFrame(),
                "hits": pd.DataFrame(),
                "interpretation": [
                    "No differential-expression result tables were available for biology validation."
                ],
                "n_signatures": len(self.signatures),
            }

        summary_df = pd.DataFrame(summary_records)
        hits_df = pd.DataFrame(hit_records)

        if not summary_df.empty:
            summary_df = summary_df.sort_values(
                [
                    "detected",
                    "n_case_up_hits",
                    "case_up_hit_fraction",
                    "n_present_genes",
                ],
                ascending=[False, False, False, False],
            ).reset_index(drop=True)

        if not hits_df.empty:
            hits_df = hits_df.sort_values(
                [
                    "signature",
                    "scope",
                    "cell_type",
                    "comparison",
                    "pvals_adj",
                    "logfoldchanges",
                ],
                ascending=[True, True, True, True, True, False],
            ).reset_index(drop=True)

        interpretation = self._build_interpretation(summary_df)

        return {
            "status": "ok",
            "summary": summary_df,
            "hits": hits_df,
            "interpretation": interpretation,
            "n_signatures": len(self.signatures),
        }

    # ------------------------------------------------------------------
    # DE iterators
    # ------------------------------------------------------------------
    def _iter_condition_de(
        self,
        results: Dict[str, Any],
    ) -> Iterable[Tuple[str, pd.DataFrame]]:
        """
        Yield whole-dataset condition DE tables.

        Expected main structure:
        results["condition_de_results"][(group1, group2)]["de_results"]
        """

        containers = [
            results.get("condition_de_results"),
            results.get("condition_de"),
            results.get("de_results_by_condition"),
        ]

        for container in containers:
            if not isinstance(container, dict):
                continue

            for key, value in container.items():
                de_results = None

                if self._is_dataframe(value):
                    de_results = value

                elif isinstance(value, dict):
                    de_results = self._first_dataframe_from_dict(
                        value,
                        ["de_results", "results", "table"],
                    )

                if not self._is_dataframe(de_results):
                    continue

                comparison = self._format_comparison(key, value)

                yield comparison, de_results

    def _iter_celltype_de(
        self,
        results: Dict[str, Any],
    ) -> Iterable[Tuple[str, str, pd.DataFrame]]:
        """
        Yield cell-type-specific DE tables.

        Handles several possible nested structures because OmniSeqAI has
        evolved over time.
        """

        candidate_keys = [
            "celltype_specific",
            "celltype_specific_de",
            "celltype_specific_de_results",
            "celltype_results",
            "cell_type_de_results",
        ]

        for key in candidate_keys:
            container = results.get(key)

            if isinstance(container, dict):
                yield from self._walk_celltype_container(container, path=[])

    def _walk_celltype_container(
        self,
        obj: Any,
        path: List[Any],
    ) -> Iterable[Tuple[str, str, pd.DataFrame]]:
        """
        Recursively walk nested cell-type DE containers.
        """

        if isinstance(obj, dict):
            de_results = self._first_dataframe_from_dict(
                obj,
                ["de_results", "results", "table"],
            )

            if self._is_dataframe(de_results):
                cell_type = self._infer_cell_type(path, obj)
                comparison = self._infer_comparison(path, obj)

                yield cell_type, comparison, de_results
                return

            for key, value in obj.items():
                yield from self._walk_celltype_container(
                    value,
                    path=path + [key],
                )

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------
    def _score_de_table(
        self,
        de_results: pd.DataFrame,
        scope: str,
        cell_type: str,
        comparison: str,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """
        Score one DE table against all signatures.
        """

        summary_records: List[Dict[str, Any]] = []
        hit_records: List[Dict[str, Any]] = []

        if not self._is_dataframe(de_results) or de_results.empty:
            return summary_records, hit_records

        df = de_results.copy()

        gene_col = self._pick_col(
            df,
            [
                "gene_symbol",
                "names",
                "gene",
                "symbol",
                "gene_name",
                "gene_id",
            ],
        )

        logfc_col = self._pick_col(
            df,
            [
                "logfoldchanges",
                "logFC",
                "log2FoldChange",
                "avg_log2FC",
            ],
        )

        padj_col = self._pick_col(
            df,
            [
                "pvals_adj",
                "padj",
                "p_val_adj",
                "FDR",
                "fdr",
                "qval",
            ],
        )

        if gene_col is None or logfc_col is None or padj_col is None:
            return summary_records, hit_records

        df["_gene_symbol_norm"] = df[gene_col].map(normalize_gene_symbol)
        df["_logfc"] = pd.to_numeric(df[logfc_col], errors="coerce")
        df["_padj"] = pd.to_numeric(df[padj_col], errors="coerce")

        df = df.dropna(subset=["_gene_symbol_norm", "_logfc", "_padj"]).copy()

        if df.empty:
            return summary_records, hit_records

        df = df[np.isfinite(df["_logfc"]) & np.isfinite(df["_padj"])].copy()

        if df.empty:
            return summary_records, hit_records

        # Keep the best row per gene by adjusted p-value, then absolute LFC.
        df["_abs_lfc"] = df["_logfc"].abs()
        df = df.sort_values(
            ["_gene_symbol_norm", "_padj", "_abs_lfc"],
            ascending=[True, True, False],
        )
        df = df.drop_duplicates("_gene_symbol_norm", keep="first").copy()

        present_genes = set(df["_gene_symbol_norm"].tolist())

        for signature, payload in self.signatures.items():
            description = str(payload.get("description", signature))
            expected_direction = str(
                payload.get("expected_direction", "case_up")
            ).strip()
            source = str(payload.get("source", "built_in"))

            expected_genes = [
                normalize_gene_symbol(gene)
                for gene in payload.get("genes", [])
            ]

            expected_genes = [
                gene
                for gene in expected_genes
                if gene and gene.lower() != "nan"
            ]

            # Remove duplicates while preserving order.
            expected_genes = list(dict.fromkeys(expected_genes))

            n_expected = len(expected_genes)

            if n_expected == 0:
                continue

            expected_set = set(expected_genes)
            n_present = len(expected_set & present_genes)

            sig_df = df[df["_gene_symbol_norm"].isin(expected_set)].copy()

            if sig_df.empty:
                n_hits = 0
                hit_fraction = 0.0
                detected = False
            else:
                sig_df["_direction_match"] = sig_df["_logfc"].map(
                    lambda value: self._direction_match(
                        value,
                        expected_direction,
                    )
                )

                hit_df = sig_df[
                    (sig_df["_padj"] <= self.padj_cutoff)
                    & (sig_df["_direction_match"])
                ].copy()

                n_hits = int(hit_df["_gene_symbol_norm"].nunique())
                hit_fraction = n_hits / n_expected if n_expected else 0.0

                detected = self._is_detected(
                    n_hits=n_hits,
                    n_expected=n_expected,
                    hit_fraction=hit_fraction,
                )

                for _, row in hit_df.iterrows():
                    hit_records.append(
                        {
                            "signature": signature,
                            "description": description,
                            "source": source,
                            "expected_direction": expected_direction,
                            "scope": scope,
                            "cell_type": cell_type,
                            "comparison": comparison,
                            "gene_symbol": row["_gene_symbol_norm"],
                            "original_gene": str(row[gene_col]),
                            "logfoldchanges": float(row["_logfc"]),
                            "pvals_adj": float(row["_padj"]),
                            "padj_cutoff": self.padj_cutoff,
                            "logfc_cutoff": self.logfc_cutoff,
                        }
                    )

            summary_records.append(
                {
                    "signature": signature,
                    "description": description,
                    "source": source,
                    "expected_direction": expected_direction,
                    "scope": scope,
                    "cell_type": cell_type,
                    "comparison": comparison,
                    "n_expected_genes": int(n_expected),
                    "n_present_genes": int(n_present),
                    "n_case_up_hits": int(n_hits),
                    "case_up_hit_fraction": float(hit_fraction),
                    "detected": bool(detected),
                    "padj_cutoff": self.padj_cutoff,
                    "logfc_cutoff": self.logfc_cutoff,
                }
            )

        return summary_records, hit_records

    # ------------------------------------------------------------------
    # Detection logic
    # ------------------------------------------------------------------
    def _direction_match(
        self,
        logfc: float,
        expected_direction: str,
    ) -> bool:
        direction = str(expected_direction).strip().lower()

        if direction in {
            "case_up",
            "up",
            "upregulated",
            "positive",
            "higher",
        }:
            return logfc >= self.logfc_cutoff

        if direction in {
            "case_down",
            "down",
            "downregulated",
            "negative",
            "lower",
        }:
            return logfc <= -self.logfc_cutoff

        if direction in {
            "either",
            "any",
            "changed",
            "both",
        }:
            return abs(logfc) >= self.logfc_cutoff

        # Conservative default.
        return logfc >= self.logfc_cutoff

    def _is_detected(
        self,
        n_hits: int,
        n_expected: int,
        hit_fraction: float,
    ) -> bool:
        """
        Decide whether a signature is detected.

        Uses both absolute hit count and fraction so that small user-defined
        signatures can still be detected without making large panels too lax.
        """

        if n_expected <= 0:
            return False

        if n_expected <= 3:
            required_hits = 1
        else:
            required_hits = max(
                self.min_hits,
                int(math.ceil(self.min_fraction * n_expected)),
            )

        return n_hits >= required_hits and hit_fraction >= self.min_fraction

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------
    def _build_interpretation(
        self,
        summary_df: pd.DataFrame,
    ) -> List[str]:
        if summary_df is None or summary_df.empty:
            return []

        detected = summary_df[summary_df["detected"] == True].copy()

        if detected.empty:
            return [
                "No predefined or user-defined biological signature reached the detection threshold."
            ]

        detected = detected.sort_values(
            [
                "n_case_up_hits",
                "case_up_hit_fraction",
                "n_present_genes",
            ],
            ascending=[False, False, False],
        )

        top = detected.iloc[0]

        lines = [
            (
                "Top detected biological program: "
                f"{top['signature']} ({top['description']})"
            ),
            (
                "Evidence: "
                f"{int(top['n_case_up_hits'])} expected genes were significantly "
                f"upregulated in the case group for {top['comparison']}."
            ),
        ]

        top_signature = str(top["signature"])

        secondary = detected[
            detected["signature"].astype(str) != top_signature
        ].copy()

        if not secondary.empty:
            unique_programs = []

            for sig in secondary["signature"].astype(str).tolist():
                if sig not in unique_programs:
                    unique_programs.append(sig)

                if len(unique_programs) >= 3:
                    break

            if unique_programs:
                lines.append(
                    "Additional detected programs: "
                    + ", ".join(unique_programs)
                    + "."
                )

        user_detected = detected[
            detected.get("source", pd.Series(index=detected.index, dtype=str))
            .astype(str)
            .eq("user")
        ]

        if not user_detected.empty:
            user_names = list(dict.fromkeys(user_detected["signature"].astype(str)))
            lines.append(
                "User-defined signatures detected: "
                + ", ".join(user_names[:5])
                + "."
            )

        return lines

    # ------------------------------------------------------------------
    # Utility functions
    # ------------------------------------------------------------------
    @classmethod
    def _first_dataframe_from_dict(
        cls,
        obj: Dict[str, Any],
        keys: List[str],
    ):
        for key in keys:
            if key in obj and cls._is_dataframe(obj.get(key)):
                return obj.get(key)

        return None

    @staticmethod
    def _is_dataframe(obj: Any) -> bool:
        return hasattr(obj, "empty") and hasattr(obj, "columns")

    @staticmethod
    def _pick_col(
        df: pd.DataFrame,
        candidates: List[str],
    ) -> Optional[str]:
        lower_to_original = {
            str(col).lower(): col
            for col in df.columns
        }

        for candidate in candidates:
            key = candidate.lower()

            if key in lower_to_original:
                return lower_to_original[key]

        return None

    @staticmethod
    def _format_comparison(
        key: Any,
        value: Any = None,
    ) -> str:
        if isinstance(value, dict):
            if value.get("comparison"):
                return str(value.get("comparison"))

            group1 = value.get("group1") or value.get("reference")
            group2 = value.get("group2") or value.get("case")

            if group1 is not None and group2 is not None:
                return f"{group1} vs {group2}"

        if isinstance(key, tuple) and len(key) == 2:
            return f"{key[0]} vs {key[1]}"

        if isinstance(key, list) and len(key) == 2:
            return f"{key[0]} vs {key[1]}"

        return str(key)

    @staticmethod
    def _infer_cell_type(
        path: List[Any],
        info: Dict[str, Any],
    ) -> str:
        for key in [
            "cell_type",
            "celltype",
            "cell_type_label",
            "cluster",
        ]:
            if isinstance(info, dict) and info.get(key) is not None:
                return str(info.get(key))

        for item in path:
            if isinstance(item, tuple):
                continue

            text = str(item)

            if " vs " in text:
                continue

            if text.startswith("(") and "," in text:
                continue

            return text

        return "unknown_cell_type"

    @classmethod
    def _infer_comparison(
        cls,
        path: List[Any],
        info: Dict[str, Any],
    ) -> str:
        if isinstance(info, dict):
            if info.get("comparison"):
                return str(info.get("comparison"))

            group1 = info.get("group1") or info.get("reference")
            group2 = info.get("group2") or info.get("case")

            if group1 is not None and group2 is not None:
                return f"{group1} vs {group2}"

        for item in path:
            if isinstance(item, tuple) and len(item) == 2:
                return f"{item[0]} vs {item[1]}"

            text = str(item)

            if " vs " in text:
                return text

        return "unknown_comparison"
