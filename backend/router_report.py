from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd


class RouterReport:

    def _clean_gene_display(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        For DE tables, hide rows whose gene names are clearly unmapped IDs.
        """
        if df is None or df.empty or "names" not in df.columns:
            return df

        out = df.copy()
        names = out["names"].astype(str).str.strip()

        bad = (
            names.str.startswith("ENSG", na=False)
            | names.str.startswith("NCBITaxon:", na=False)
            | names.str.startswith("MIR", na=False)  # optional; comment out if you want these kept
        )

        # Keep rows that look like real symbols
        cleaned = out.loc[~bad].copy()

        # If cleaning would remove everything, fall back to the original table
        if cleaned.empty:
            return out

        return cleaned

    def _format_df_head(
        self,
        df: Optional[pd.DataFrame],
        cols=None,
        n: int = 10,
        clean_gene_names: bool = False
    ) -> str:
        if df is None:
            return "No table available."

        if df.empty:
            return "No results."

        out = df.copy()

        if clean_gene_names:
            out = self._clean_gene_display(out)

        if cols is not None:
            cols = [c for c in cols if c in out.columns]
            if cols:
                out = out[cols]

        if out.empty:
            return "No symbol-mapped results."

        if isinstance(out.index, pd.RangeIndex):
            return out.head(n).to_string(index=False)

        return out.head(n).to_string()

    def _format_profile(self, profile: Dict[str, Any]) -> str:
        lines = []
        lines.append("===== DETECTED PROFILE =====")
        for k in [
            "cell_type_column",
            "condition_column",
            "sample_column",
            "patient_column",
            "batch_column",
            "gene_symbol_source",
            "has_umap",
            "has_pca",
            "has_raw",
            "n_cells",
            "n_genes",
        ]:
            lines.append(f"{k}: {profile.get(k)}")
        return "\n".join(lines)

    def _format_plan(self, plan) -> str:
        lines = []
        lines.append("===== ANALYSIS PLAN =====")
        lines.append(f"has_conditions: {getattr(plan, 'has_conditions', None)}")
        lines.append(f"condition_column: {getattr(plan, 'condition_column', None)}")
        lines.append(f"cell_type_column: {getattr(plan, 'cell_type_column', None)}")
        lines.append(f"sample_column: {getattr(plan, 'sample_column', None)}")
        lines.append(f"patient_column: {getattr(plan, 'patient_column', None)}")
        lines.append(f"batch_column: {getattr(plan, 'batch_column', None)}")
        lines.append(f"condition_groups: {getattr(plan, 'condition_groups', None)}")
        lines.append(f"baseline: {getattr(plan, 'baseline', None)}")
        lines.append(f"pairwise_comparisons: {getattr(plan, 'pairwise_comparisons', None)}")
        lines.append(f"use_pseudobulk: {getattr(plan, 'use_pseudobulk', None)}")
        lines.append(f"pseudobulk_key: {getattr(plan, 'pseudobulk_key', None)}")
        return "\n".join(lines)

    def _format_interpretation(self, items) -> str:
        if not items:
            return "No interpretation available."
        return "\n".join([f"- {x}" for x in items])

    def _top_celltype_summary(self, results: Dict[str, Any], top_n: int = 3) -> str:
        celltype_specific = results.get("celltype_specific", {}) or {}

        ranked = []
        for cell_type, pair_map in celltype_specific.items():
            best = None

            for pair, info in pair_map.items():
                if not isinstance(info, dict):
                    continue
                if info.get("status") not in {None, "ok"}:
                    continue

                n_sig = info.get("n_sig_genes", 0) or 0
                if best is None or n_sig > best["n_sig_genes"]:
                    best = {
                        "cell_type": cell_type,
                        "pair": pair,
                        "n_sig_genes": int(n_sig),
                        "mode": info.get("mode"),
                        "interpretation": info.get("interpretation", []),
                    }

            if best is not None:
                ranked.append(best)

        if not ranked:
            return "No cell-type-specific results available."

        ranked = sorted(ranked, key=lambda x: x["n_sig_genes"], reverse=True)[:top_n]

        lines = []
        for item in ranked:
            interp = item["interpretation"] or []
            interp_text = "; ".join(interp[:2]) if interp else "No interpretation available."
            lines.append(
                f"- {item['cell_type']}: {item['n_sig_genes']} significant genes "
                f"({item['mode']}, {item['pair'][0]} vs {item['pair'][1]}). "
                f"{interp_text}"
            )

        return "\n".join(lines)

    def _top_dataset_summary(self, results: Dict[str, Any]) -> str:
        plan = results.get("plan", None)
        adata = results.get("adata", None)
        profile = results.get("profile", {}) or {}
        condition_de_results = results.get("condition_de_results", {}) or {}

        lines = []
        lines.append("===== EXECUTIVE SUMMARY =====")

        if adata is not None:
            lines.append(f"Dataset size: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        if profile:
            lines.append(
                "Detected metadata: "
                f"cell_type={profile.get('cell_type_column')}, "
                f"condition={profile.get('condition_column')}, "
                f"sample={profile.get('sample_column')}, "
                f"patient={profile.get('patient_column')}"
            )

        if plan is not None:
            mode = "pseudobulk" if getattr(plan, "use_pseudobulk", False) else "cell-level"
            lines.append(f"Primary analysis mode: {mode}")
            lines.append(f"Baseline/control: {getattr(plan, 'baseline', None)}")
            lines.append(f"Pairwise comparisons: {getattr(plan, 'pairwise_comparisons', None)}")

        if condition_de_results:
            for pair, info in condition_de_results.items():
                lines.append(
                    f"Whole-dataset DE: {pair[0]} vs {pair[1]} "
                    f"({info.get('mode')}, {info.get('n_sig_genes')} significant genes)"
                )
                interp = info.get("interpretation", [])
                if interp:
                    lines.append("Whole-dataset interpretation:")
                    for item in interp:
                        lines.append(f"- {item}")
                break
        else:
            lines.append("Whole-dataset DE: not available.")

        lines.append("Top responsive cell types:")
        lines.append(self._top_celltype_summary(results, top_n=3))

        if profile.get("gene_symbol_source") in {"feature_name", "gene_symbol"}:
            lines.append(
                "Note: a small number of features may still appear as Ensembl IDs "
                "if the source dataset did not provide a symbol for them."
            )

        return "\n".join(lines)

    def _build_exploratory_report(self, results: Dict[str, Any]) -> str:
        """
        Build a clean report for exploratory mode.

        Exploratory datasets may not have condition, sample, patient,
        or replicate metadata. This report avoids condition-DE and
        pseudobulk sections.
        """

        from datetime import datetime

        adata = results.get("adata")
        profile = results.get("profile")
        decision = results.get("decision")
        exploratory = results.get("exploratory_results", {})
        figure_paths = results.get("figure_paths", {})

        label_col = exploratory.get("label_column")
        cluster_col = exploratory.get("cluster_column")

        lines = []

        lines.append("OmniSeqAI Exploratory Analysis Report")
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
        lines.append("Selected mode: exploratory")

        if decision is not None:
            reason = getattr(decision, "reason", None)
            if reason:
                lines.append(reason)

        lines.append("")
        lines.append("Executive Summary")

        if adata is not None:
            lines.append(f"- Dataset size: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        if label_col:
            lines.append(f"- Main annotation column: {label_col}")

        if cluster_col:
            lines.append(f"- Cluster column: {cluster_col}")

        if exploratory.get("has_umap"):
            lines.append("- UMAP embedding available/generated")

        if exploratory.get("has_marker_genes"):
            lines.append("- Marker genes calculated")

        lines.append("")
        lines.append("Dataset Overview")

        if adata is not None:
            lines.append(f"Cells: {adata.n_obs:,}")
            lines.append(f"Genes: {adata.n_vars:,}")

        if profile is not None:
            profile_cell_type_col = getattr(profile, "cell_type_col", None)
            profile_condition_col = getattr(profile, "condition_col", None)
            profile_sample_col = getattr(profile, "sample_col", None)
            profile_patient_col = getattr(profile, "patient_col", None)
            profile_batch_col = getattr(profile, "batch_col", None)
        else:
            profile_cell_type_col = None
            profile_condition_col = None
            profile_sample_col = None
            profile_patient_col = None
            profile_batch_col = None

        # Use exploratory label column as fallback if detector did not classify it.
        display_cell_type_col = profile_cell_type_col or label_col

        lines.append(f"Cell type / annotation column: {display_cell_type_col or 'None detected'}")
        lines.append(f"Condition column: {profile_condition_col or 'None detected'}")
        lines.append(f"Sample column: {profile_sample_col or 'None detected'}")
        lines.append(f"Patient column: {profile_patient_col or 'None detected'}")
        lines.append(f"Batch column: {profile_batch_col or 'None detected'}")

        lines.append("")
        lines.append("Exploratory Workflow")

        steps = exploratory.get("steps", [])

        if steps:
            for step in steps:
                lines.append(f"- {step}")
        else:
            lines.append("- No exploratory workflow log available.")

        lines.append("")
        lines.append("Figures")

        if figure_paths:
            for name, path in figure_paths.items():
                lines.append(f"- {name}: {path}")
        else:
            lines.append("- No figures generated.")

        lines.append("")
        lines.append("Marker Genes")

        if adata is not None and "rank_genes_groups" in adata.uns:
            try:
                rgg = adata.uns["rank_genes_groups"]
                names = rgg["names"]

                if hasattr(names, "dtype") and names.dtype.names is not None:
                    groups = list(names.dtype.names)

                    for group in groups[:12]:
                        top_genes = [str(g) for g in names[group][:10]]
                        lines.append(f"{group}: {', '.join(top_genes)}")
                else:
                    lines.append("Marker genes were calculated but could not be summarized cleanly.")

            except Exception as e:
                lines.append(f"Marker gene summary unavailable: {type(e).__name__}: {e}")
        else:
            lines.append("No marker gene results available.")

        lines.append("")
        lines.append("Interpretation")

        lines.append(
            "This dataset was processed in exploratory mode because no valid "
            "condition/sample/patient structure was detected for pseudobulk "
            "condition analysis."
        )

        if label_col:
            lines.append(
                f"Cells were visualized and summarized using the annotation column '{label_col}'."
            )

        if cluster_col:
            lines.append(
                f"Unsupervised clustering results are available in '{cluster_col}'."
            )

        lines.append("")
        lines.append("Methods")

        lines.append("- OmniSeqAI detected available metadata columns before selecting the analysis mode.")
        lines.append("- Exploratory mode was used because no valid condition contrast was available.")
        lines.append("- QC, filtering, normalization, PCA, neighbor graph construction, UMAP, and clustering were run when supported.")
        lines.append("- Marker genes were calculated using Scanpy rank_genes_groups where possible.")
        lines.append("- No pseudobulk DE was run because biological condition and replicate metadata were unavailable.")
        lines.append("- No condition-specific interpretation was generated for this dataset.")

        # Remove accidental empty trailing whitespace lines.
        cleaned_lines = [str(line).rstrip() for line in lines]

        return "\n".join(cleaned_lines)

    def build(self, results: Dict[str, Any]) -> str:
        decision = results.get("decision")
        mode = getattr(decision, "mode", None)

        if mode == "exploratory":
            return self._build_exploratory_report(results)
        profile = results.get("profile", {})
        plan = results.get("plan", None)
        adata = results.get("adata", None)

        lines = []
        lines.append("========== OMNISEQAI ROUTER REPORT ==========\n")
        lines.append(self._top_dataset_summary(results))
        lines.append("")

        if adata is not None:
            lines.append(f"Cells: {adata.n_obs:,}")
            lines.append(f"Genes: {adata.n_vars:,}\n")

        lines.append(self._format_profile(profile))
        lines.append("")

        if plan is not None:
            lines.append(self._format_plan(plan))
            lines.append("")

        celltype_comparison = results.get("celltype_comparison", None)
        lines.append("===== CELL-TYPE COMPARISON =====")
        if celltype_comparison is not None:
            lines.append(self._format_df_head(celltype_comparison, n=20))
        else:
            lines.append("No comparison produced.")
        lines.append("")

        lines.append("===== WHOLE-DATASET CONDITION DE =====")
        condition_de_results = results.get("condition_de_results", {})
        if condition_de_results:
            for pair, info in condition_de_results.items():
                lines.append("\n" + "-" * 60)
                lines.append(f"Comparison: {pair}")
                lines.append(f"Mode: {info.get('mode')}")
                lines.append(f"Status: {info.get('status')}")
                lines.append(f"Significant genes: {info.get('n_sig_genes')}")
                if info.get("error"):
                    lines.append(f"Error: {info.get('error')}")
                lines.append("Top DE genes:")
                lines.append(
                    self._format_df_head(
                        info.get("de_results"),
                        cols=["names", "logfoldchanges", "pvals_adj"],
                        n=10,
                        clean_gene_names=True
                    )
                )
                lines.append("Interpretation:")
                lines.append(self._format_interpretation(info.get("interpretation", [])))
        else:
            lines.append("No whole-dataset DE results available.")

        lines.append("\n")
        lines.append("===== CELL-TYPE-SPECIFIC RESULTS =====")
        celltype_specific = results.get("celltype_specific", {})
        if celltype_specific:
            for cell_type, pair_map in celltype_specific.items():
                lines.append("\n" + "=" * 70)
                lines.append(f"Cell type: {cell_type}")

                for pair, info in pair_map.items():
                    lines.append(f"Pair: {pair}")
                    lines.append(f"Mode: {info.get('mode')}")
                    lines.append(f"Status: {info.get('status')}")
                    lines.append(f"Condition counts: {info.get('condition_counts')}")
                    lines.append(f"Significant genes: {info.get('n_sig_genes')}")

                    if info.get("error"):
                        lines.append(f"Error: {info.get('error')}")

                    lines.append("Top DE genes:")
                    lines.append(
                        self._format_df_head(
                            info.get("de_results"),
                            cols=["names", "logfoldchanges", "pvals_adj"],
                            n=5,
                            clean_gene_names=True
                        )
                    )
                    lines.append("Interpretation:")
                    lines.append(self._format_interpretation(info.get("interpretation", [])))
                    lines.append("")
        else:
            lines.append("No cell-type-specific results available.")

        return "\n".join(lines)

    def save(self, report_text: str, filename: str = "reports/router_report.txt") -> None:
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(report_text)
        print(f"\nRouter report saved: {path}")
