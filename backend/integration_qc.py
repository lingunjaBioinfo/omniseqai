from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import pandas as pd
import scanpy as sc


@dataclass
class IntegrationQCResult:
    status: str
    message: str
    sample_col: Optional[str]
    batch_col: Optional[str]
    condition_col: Optional[str]
    n_samples: int
    n_batches: int
    method: Optional[str]
    integrated_adata: Any
    tables: Dict[str, pd.DataFrame]


class IntegrationQC:
    """
    Multi-sample and batch-aware QC helper.

    This module is intentionally conservative:
    - It summarizes sample/batch structure.
    - It optionally attempts Harmony integration when batch structure exists.
    - It does not modify the original AnnData object used for DE.
    """

    SAMPLE_CANDIDATES = [
        "sample_id",
        "sample",
        "Sample",
        "orig.ident",
        "orig_ident",
        "donor",
        "donor_id",
        "patient_id",
    ]

    BATCH_CANDIDATES = [
        "batch_id",
        "batch",
        "Batch",
        "library",
        "library_id",
        "lane",
        "sequencing_batch",
    ]

    CONDITION_CANDIDATES = [
        "condition",
        "Condition",
        "disease",
        "status",
        "group",
        "treatment",
        "stim",
    ]

    def run(self, adata, profile=None, attempt_integration: bool = True) -> Dict[str, Any]:
        sample_col = self._resolve_column(
            adata,
            profile,
            profile_attrs=["sample_col", "sample_column"],
            candidates=self.SAMPLE_CANDIDATES,
        )

        batch_col = self._resolve_column(
            adata,
            profile,
            profile_attrs=["batch_col", "batch_column"],
            candidates=self.BATCH_CANDIDATES,
        )

        condition_col = self._resolve_column(
            adata,
            profile,
            profile_attrs=["condition_col", "condition_column"],
            candidates=self.CONDITION_CANDIDATES,
        )

        tables = self._build_tables(
            adata=adata,
            sample_col=sample_col,
            batch_col=batch_col,
            condition_col=condition_col,
        )

        n_samples = self._n_unique(adata, sample_col)
        n_batches = self._n_unique(adata, batch_col)

        result = IntegrationQCResult(
            status="ok",
            message="Multi-sample QC completed.",
            sample_col=sample_col,
            batch_col=batch_col,
            condition_col=condition_col,
            n_samples=n_samples,
            n_batches=n_batches,
            method=None,
            integrated_adata=None,
            tables=tables,
        )

        if n_samples < 2 and n_batches < 2:
            result.status = "skipped"
            result.message = "No multi-sample or batch structure detected."
            return self._to_dict(result)

        if attempt_integration and batch_col is not None and n_batches >= 2:
            integration = self._try_harmony_integration(
                adata=adata,
                batch_col=batch_col,
            )

            result.method = integration.get("method")
            result.integrated_adata = integration.get("adata")

            if integration.get("status") == "ok":
                result.status = "integrated"
                result.message = integration.get("message")
            else:
                result.status = "qc_only"
                result.message = (
                    "Batch-aware QC completed, but integration was not performed. "
                    f"Reason: {integration.get('message')}"
                )

        else:
            result.status = "qc_only"
            result.message = "Multi-sample QC completed without integration."

        return self._to_dict(result)

    # --------------------------------------------------
    # Column detection
    # --------------------------------------------------
    def _resolve_column(self, adata, profile, profile_attrs, candidates):
        if profile is not None:
            for attr in profile_attrs:
                value = None

                if isinstance(profile, dict):
                    value = profile.get(attr)

                else:
                    value = getattr(profile, attr, None)

                if value in adata.obs.columns:
                    return value

        for col in candidates:
            if col in adata.obs.columns:
                return col

        return None

    def _n_unique(self, adata, col: Optional[str]) -> int:
        if col is None or col not in adata.obs.columns:
            return 0

        return int(adata.obs[col].astype(str).nunique())

    # --------------------------------------------------
    # QC tables
    # --------------------------------------------------
    def _build_tables(
        self,
        adata,
        sample_col: Optional[str],
        batch_col: Optional[str],
        condition_col: Optional[str],
    ) -> Dict[str, pd.DataFrame]:
        tables: Dict[str, pd.DataFrame] = {}

        if sample_col is not None:
            tables["sample_counts"] = self._count_table(
                adata,
                sample_col,
                output_col="sample",
            )

        if batch_col is not None:
            tables["batch_counts"] = self._count_table(
                adata,
                batch_col,
                output_col="batch",
            )

        if condition_col is not None:
            tables["condition_counts"] = self._count_table(
                adata,
                condition_col,
                output_col="condition",
            )

        if sample_col is not None and condition_col is not None:
            tables["sample_condition_counts"] = pd.crosstab(
                adata.obs[sample_col].astype(str),
                adata.obs[condition_col].astype(str),
            )

        if batch_col is not None and condition_col is not None:
            tables["batch_condition_counts"] = pd.crosstab(
                adata.obs[batch_col].astype(str),
                adata.obs[condition_col].astype(str),
            )

        if sample_col is not None:
            sample_qc = self._qc_summary_by_group(adata, sample_col)

            if sample_qc is not None and not sample_qc.empty:
                tables["sample_qc_summary"] = sample_qc

        if batch_col is not None:
            batch_qc = self._qc_summary_by_group(adata, batch_col)

            if batch_qc is not None and not batch_qc.empty:
                tables["batch_qc_summary"] = batch_qc

        return tables

    def _count_table(self, adata, col: str, output_col: str) -> pd.DataFrame:
        return (
            adata.obs[col]
            .astype(str)
            .value_counts()
            .rename_axis(output_col)
            .reset_index(name="n_cells")
        )

    def _qc_summary_by_group(self, adata, group_col: str) -> Optional[pd.DataFrame]:
        qc_cols = [
            col
            for col in [
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt",
                "pct_counts_ribo",
                "pct_counts_hb",
            ]
            if col in adata.obs.columns
        ]

        if not qc_cols:
            return None

        grouped = adata.obs.groupby(adata.obs[group_col].astype(str))[qc_cols]

        summary = grouped.agg(["mean", "median", "min", "max"])

        summary.columns = [
            f"{metric}_{stat}"
            for metric, stat in summary.columns.to_flat_index()
        ]

        summary = summary.reset_index().rename(columns={group_col: "group"})

        return summary

    # --------------------------------------------------
    # Optional Harmony integration
    # --------------------------------------------------
    def _try_harmony_integration(self, adata, batch_col: str) -> Dict[str, Any]:
        try:
            import scanpy.external as sce  # noqa: F401

        except Exception as e:
            return {
                "status": "skipped",
                "method": None,
                "message": f"scanpy.external not available: {e}",
                "adata": None,
            }

        try:
            adata_int = adata.copy()

            # Work on a copy only. Do not alter original counts used for DE.
            if "X_pca" not in adata_int.obsm:
                sc.pp.normalize_total(adata_int, target_sum=1e4)
                sc.pp.log1p(adata_int)

                n_top_genes = min(2000, adata_int.n_vars)

                try:
                    sc.pp.highly_variable_genes(
                        adata_int,
                        n_top_genes=n_top_genes,
                        flavor="seurat",
                    )
                    use_hvg = True
                except Exception:
                    use_hvg = False

                n_comps = min(50, adata_int.n_obs - 1, adata_int.n_vars - 1)
                n_comps = max(2, n_comps)

                sc.tl.pca(
                    adata_int,
                    n_comps=n_comps,
                    use_highly_variable=use_hvg,
                    svd_solver="arpack",
                )

            import scanpy.external as sce

            sce.pp.harmony_integrate(
                adata_int,
                key=batch_col,
                basis="X_pca",
                adjusted_basis="X_pca_harmony",
            )

            sc.pp.neighbors(
                adata_int,
                use_rep="X_pca_harmony",
            )

            sc.tl.umap(adata_int)

            return {
                "status": "ok",
                "method": "harmony",
                "message": "Harmony integration completed.",
                "adata": adata_int,
            }

        except Exception as e:
            return {
                "status": "failed",
                "method": "harmony",
                "message": str(e),
                "adata": None,
            }

    def _to_dict(self, result: IntegrationQCResult) -> Dict[str, Any]:
        return {
            "status": result.status,
            "message": result.message,
            "sample_col": result.sample_col,
            "batch_col": result.batch_col,
            "condition_col": result.condition_col,
            "n_samples": result.n_samples,
            "n_batches": result.n_batches,
            "method": result.method,
            "integrated_adata": result.integrated_adata,
            "tables": result.tables,
        }
