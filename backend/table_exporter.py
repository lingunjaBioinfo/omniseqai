from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd


class TableExporter:
    """
    Export OmniSeqAI result tables to CSV/JSON.

    Exports:
    - run summary
    - whole-dataset condition DE
    - significant genes
    - cell-type-specific DE where available
    - exploratory marker genes
    - cell-type counts/proportions
    """

    def __init__(self, output_dir: str = "outputs/tables"):
        self.output_dir = Path(output_dir)

    # --------------------------------------------------
    # Public API
    # --------------------------------------------------
    def export(self, results: Dict[str, Any], output_dir: Optional[str] = None) -> Dict[str, str]:
        outdir = Path(output_dir) if output_dir else self.output_dir
        outdir.mkdir(parents=True, exist_ok=True)

        table_paths: Dict[str, str] = {}

        summary_path = self._export_run_summary(results, outdir)
        table_paths["run_summary"] = str(summary_path)

        self._export_condition_de(results, outdir, table_paths)
        self._export_celltype_specific_de(results, outdir, table_paths)
        self._export_marker_genes(results, outdir, table_paths)
        self._export_celltype_counts(results, outdir, table_paths)
        self._export_integration_qc(results, outdir, table_paths)

        return table_paths

    # --------------------------------------------------
    # Helpers
    # --------------------------------------------------
    def _safe_name(self, value: Any) -> str:
        if isinstance(value, tuple):
            value = "_vs_".join(str(v) for v in value)

        value = str(value)
        value = value.replace(" ", "_")
        value = value.replace("/", "_")
        value = value.replace("\\", "_")
        value = re.sub(r"[^A-Za-z0-9_.+-]+", "_", value)
        value = value.strip("_")

        return value or "unnamed"

    def _is_dataframe(self, obj: Any) -> bool:
        return hasattr(obj, "to_csv") and hasattr(obj, "columns")

    def _significant_subset(self, df: pd.DataFrame) -> pd.DataFrame:
        if df is None or df.empty:
            return pd.DataFrame()

        if "pvals_adj" in df.columns:
            padj_col = "pvals_adj"
        elif "padj" in df.columns:
            padj_col = "padj"
        else:
            padj_col = None

        if "logfoldchanges" in df.columns:
            lfc_col = "logfoldchanges"
        elif "logFC" in df.columns:
            lfc_col = "logFC"
        elif "log2FoldChange" in df.columns:
            lfc_col = "log2FoldChange"
        else:
            lfc_col = None

        if padj_col is not None and lfc_col is not None:
            mask = (df[padj_col] < 0.05) & (df[lfc_col].abs() >= 0.5)
            return df.loc[mask].copy()

        if padj_col is not None:
            mask = df[padj_col] < 0.05
            return df.loc[mask].copy()

        return pd.DataFrame()

    # --------------------------------------------------
    # Run summary
    # --------------------------------------------------
    def _export_run_summary(self, results: Dict[str, Any], outdir: Path) -> Path:
        adata = results.get("adata")
        profile = results.get("profile")
        decision = results.get("decision")

        summary = {
            "selected_mode": getattr(decision, "mode", None),
            "reason": getattr(decision, "reason", None),
            "n_cells": int(adata.n_obs) if adata is not None else None,
            "n_genes": int(adata.n_vars) if adata is not None else None,
            "has_condition_de": bool(results.get("condition_de_results")),
            "has_exploratory_results": bool(results.get("exploratory_results")),
            "has_figures": bool(results.get("figure_paths")),
            "cell_type_col": getattr(profile, "cell_type_col", None) if profile else None,
            "condition_col": getattr(profile, "condition_col", None) if profile else None,
            "sample_col": getattr(profile, "sample_col", None) if profile else None,
            "patient_col": getattr(profile, "patient_col", None) if profile else None,
            "batch_col": getattr(profile, "batch_col", None) if profile else None,
        }

        path = outdir / "run_summary.json"

        with open(path, "w") as f:
            json.dump(summary, f, indent=2, default=str)

        return path

    # --------------------------------------------------
    # Condition DE
    # --------------------------------------------------
    def _export_condition_de(
        self,
        results: Dict[str, Any],
        outdir: Path,
        table_paths: Dict[str, str],
    ) -> None:
        condition_de_results = results.get("condition_de_results", {})

        if not condition_de_results:
            return

        cond_dir = outdir / "condition_de"
        cond_dir.mkdir(parents=True, exist_ok=True)

        for pair, info in condition_de_results.items():
            if not isinstance(info, dict):
                continue

            de_results = info.get("de_results")

            if not self._is_dataframe(de_results) or de_results.empty:
                continue

            comp_name = self._safe_name(pair)

            all_path = cond_dir / f"{comp_name}_all_genes.csv"
            de_results.to_csv(all_path, index=False)
            table_paths[f"condition_de.{comp_name}.all_genes"] = str(all_path)

            sig_df = self._significant_subset(de_results)

            if not sig_df.empty:
                sig_path = cond_dir / f"{comp_name}_significant_genes.csv"
                sig_df.to_csv(sig_path, index=False)
                table_paths[f"condition_de.{comp_name}.significant_genes"] = str(sig_path)

    # --------------------------------------------------
    # Cell-type-specific DE
    # --------------------------------------------------
    def _export_celltype_specific_de(
        self,
        results: Dict[str, Any],
        outdir: Path,
        table_paths: Dict[str, str],
    ) -> None:
        celltype_specific = results.get("celltype_specific", {})

        if not celltype_specific:
            return

        ct_dir = outdir / "celltype_specific_de"
        ct_dir.mkdir(parents=True, exist_ok=True)

        self._walk_de_objects(
            obj=celltype_specific,
            outdir=ct_dir,
            table_paths=table_paths,
            prefix="celltype_specific",
            name="celltype",
        )

    def _walk_de_objects(
        self,
        obj: Any,
        outdir: Path,
        table_paths: Dict[str, str],
        prefix: str,
        name: str,
    ) -> None:
        if obj is None:
            return

        if isinstance(obj, dict):
            de_results = obj.get("de_results")

            if self._is_dataframe(de_results) and not de_results.empty:
                safe = self._safe_name(name)

                all_path = outdir / f"{safe}_all_genes.csv"
                de_results.to_csv(all_path, index=False)
                table_paths[f"{prefix}.{safe}.all_genes"] = str(all_path)

                sig_df = self._significant_subset(de_results)

                if not sig_df.empty:
                    sig_path = outdir / f"{safe}_significant_genes.csv"
                    sig_df.to_csv(sig_path, index=False)
                    table_paths[f"{prefix}.{safe}.significant_genes"] = str(sig_path)

                return

            for key, value in obj.items():
                child_name = f"{name}_{self._safe_name(key)}"
                self._walk_de_objects(
                    value,
                    outdir=outdir,
                    table_paths=table_paths,
                    prefix=prefix,
                    name=child_name,
                )

    # --------------------------------------------------
    # Marker genes
    # --------------------------------------------------
    def _export_marker_genes(
        self,
        results: Dict[str, Any],
        outdir: Path,
        table_paths: Dict[str, str],
    ) -> None:
        adata = results.get("adata")

        if adata is None:
            return

        if "rank_genes_groups" not in adata.uns:
            return

        rgg = adata.uns["rank_genes_groups"]

        if "names" not in rgg:
            return

        names = rgg["names"]

        if not hasattr(names, "dtype") or names.dtype.names is None:
            return

        marker_dir = outdir / "marker_genes"
        marker_dir.mkdir(parents=True, exist_ok=True)

        groups = list(names.dtype.names)
        rows = []

        for group in groups:
            group_rows = []

            for rank, gene in enumerate(names[group], start=1):
                row = {
                    "group": group,
                    "rank": rank,
                    "gene": str(gene),
                }

                for key in ["scores", "logfoldchanges", "pvals", "pvals_adj"]:
                    try:
                        if key in rgg:
                            row[key] = rgg[key][group][rank - 1]
                    except Exception:
                        pass

                rows.append(row)
                group_rows.append(row)

            group_df = pd.DataFrame(group_rows)

            if not group_df.empty:
                group_path = marker_dir / f"{self._safe_name(group)}_markers.csv"
                group_df.to_csv(group_path, index=False)
                table_paths[f"marker_genes.{self._safe_name(group)}"] = str(group_path)

        all_df = pd.DataFrame(rows)

        if not all_df.empty:
            all_path = marker_dir / "all_marker_genes.csv"
            all_df.to_csv(all_path, index=False)
            table_paths["marker_genes.all"] = str(all_path)

    # --------------------------------------------------
    # Cell-type counts and proportions
    # --------------------------------------------------
    def _export_celltype_counts(
        self,
        results: Dict[str, Any],
        outdir: Path,
        table_paths: Dict[str, str],
    ) -> None:
        adata = results.get("adata")

        if adata is None:
            return

        if "cell_type" not in adata.obs.columns:
            return

        counts = (
            adata.obs["cell_type"]
            .astype(str)
            .value_counts()
            .rename_axis("cell_type")
            .reset_index(name="n_cells")
        )

        counts_path = outdir / "celltype_counts.csv"
        counts.to_csv(counts_path, index=False)
        table_paths["celltype_counts"] = str(counts_path)

        if "condition" in adata.obs.columns:
            ct = pd.crosstab(
                adata.obs["cell_type"].astype(str),
                adata.obs["condition"].astype(str),
            )

            ct_path = outdir / "celltype_counts_by_condition.csv"
            ct.to_csv(ct_path)
            table_paths["celltype_counts_by_condition"] = str(ct_path)

            proportions = ct.div(ct.sum(axis=0), axis=1) * 100.0

            prop_path = outdir / "celltype_proportions_by_condition.csv"
            proportions.to_csv(prop_path)
            table_paths["celltype_proportions_by_condition"] = str(prop_path)

    # --------------------------------------------------
    # Integration / batch-aware QC
    # --------------------------------------------------
    def _export_integration_qc(
        self,
        results: Dict[str, Any],
        outdir: Path,
        table_paths: Dict[str, str],
    ) -> None:
        integration_qc = results.get("integration_qc", {}) or {}

        tables = integration_qc.get("tables", {}) or {}

        if not tables:
            return

        int_dir = outdir / "integration_qc"
        int_dir.mkdir(parents=True, exist_ok=True)

        for name, table in tables.items():
            if not self._is_dataframe(table):
                continue

            if table.empty:
                continue

            safe_name = self._safe_name(name)
            path = int_dir / f"{safe_name}.csv"

            table.to_csv(path, index=True)

            table_paths[f"integration_qc.{safe_name}"] = str(path)
