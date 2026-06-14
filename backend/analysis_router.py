from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

from backend.gene_symbol_utils import ensure_de_gene_symbols
from backend.signature_genes import find_signature_hits, format_signature_hits
from backend.pseudobulk_de import PseudobulkDE


@dataclass
class AnalysisPlan:
    has_conditions: bool = False
    condition_column: Optional[str] = None
    cell_type_column: Optional[str] = None
    sample_column: Optional[str] = None
    patient_column: Optional[str] = None
    batch_column: Optional[str] = None

    condition_groups: List[str] = field(default_factory=list)
    baseline: Optional[str] = None
    pairwise_comparisons: List[Tuple[str, str]] = field(default_factory=list)

    run_condition_de: bool = False
    use_pseudobulk: bool = False
    pseudobulk_key: Optional[str] = None


class AnalysisRouter:
    """
    OmniSeqAI condition-analysis router.

    Core DE convention:
        pair = (baseline, case)

    Therefore:
        ('Healthy', 'IFN_beta')

    means:
        positive logfoldchanges = IFN_beta higher than Healthy
    """

    def __init__(self):
        self.pseudobulk_de = PseudobulkDE()

    # ------------------------------------------------------------------
    # Metadata helpers
    # ------------------------------------------------------------------
    def _profile_get(self, profile, key: str, default=None):
        if profile is None:
            return default

        if isinstance(profile, dict):
            return profile.get(key, default)

        return getattr(profile, key, default)

    def _clean_condition(self, value: Any) -> str:
        value = str(value).strip()

        mapping = {
            "normal": "Healthy",
            "Normal": "Healthy",
            "healthy": "Healthy",
            "Healthy": "Healthy",
            "control": "Healthy",
            "Control": "Healthy",
            "ctrl": "Healthy",
            "CTRL": "Healthy",
            "untreated": "Healthy",
            "Untreated": "Healthy",

            "covid": "COVID",
            "COVID": "COVID",
            "covid-19": "COVID",
            "COVID-19": "COVID",
            "sars-cov-2": "COVID",
            "SARS-CoV-2": "COVID",

            "stim": "IFN_beta",
            "STIM": "IFN_beta",
            "stimulated": "IFN_beta",
            "Stimulated": "IFN_beta",
            "ifn": "IFN_beta",
            "IFN": "IFN_beta",
            "ifnb": "IFN_beta",
            "IFNB": "IFN_beta",
            "IFN-beta": "IFN_beta",
            "IFN_beta": "IFN_beta",
            "interferon": "IFN_beta",

            "tumor": "Tumor",
            "Tumor": "Tumor",
            "cancer": "Tumor",
            "Cancer": "Tumor",
            "disease": "Disease",
            "Disease": "Disease",
        }

        return mapping.get(value, value)

    def _standardize_metadata(self, adata, profile=None):
        """
        Standardize common dataset metadata columns.

        Creates/normalizes:
        - condition
        - condition_raw
        - cell_type
        - patient_id
        - sample_id
        - batch_id
        - gene_symbol
        """

        # -----------------------------
        # Condition column
        # -----------------------------
        condition_col = self._profile_get(profile, "condition_column", None)

        if condition_col is None:
            for col in [
                "condition",
                "condition_raw",
                "disease",
                "status",
                "group",
                "label",
                "treatment",
                "stim",
                "case_control",
                "phenotype",
            ]:
                if col in adata.obs.columns:
                    condition_col = col
                    break

        if condition_col is not None and condition_col in adata.obs.columns:
            adata.obs["condition_raw"] = adata.obs[condition_col].astype(str)
            adata.obs["condition"] = (
                adata.obs[condition_col]
                .astype(str)
                .map(self._clean_condition)
                .astype(str)
            )

        # -----------------------------
        # Cell type column
        # -----------------------------
        cell_type_col = self._profile_get(profile, "cell_type_column", None)

        if cell_type_col is None:
            for col in [
                "cell_type",
                "celltype",
                "cell_type_label",
                "seurat_annotations",
                "seurat_annotation",
                "annotation",
                "cluster",
                "louvain",
                "leiden",
                "seurat_clusters",
            ]:
                if col in adata.obs.columns:
                    cell_type_col = col
                    break

        if cell_type_col is not None and cell_type_col in adata.obs.columns:
            adata.obs["cell_type"] = adata.obs[cell_type_col].astype(str)
        elif "cell_type" not in adata.obs.columns:
            adata.obs["cell_type"] = "Unknown"

        # -----------------------------
        # Patient / donor column
        # -----------------------------
        patient_col = self._profile_get(profile, "patient_column", None)

        if patient_col is None:
            for col in [
                "patient_id",
                "donor_id",
                "donor",
                "replicate",
                "ind",
                "individual",
                "subject",
                "patient",
            ]:
                if col in adata.obs.columns:
                    patient_col = col
                    break

        if patient_col is not None and patient_col in adata.obs.columns:
            adata.obs["patient_id"] = adata.obs[patient_col].astype(str)
        elif "patient_id" not in adata.obs.columns:
            adata.obs["patient_id"] = "patient_1"

        # -----------------------------
        # Sample column
        # -----------------------------
        sample_col = self._profile_get(profile, "sample_column", None)

        if sample_col is None:
            for col in [
                "sample_id",
                "sample",
                "sampleName_with_day",
                "orig.ident",
                "library",
                "batch",
            ]:
                if col in adata.obs.columns:
                    sample_col = col
                    break

        if sample_col is not None and sample_col in adata.obs.columns:
            adata.obs["sample_id"] = adata.obs[sample_col].astype(str)
        elif "sample_id" not in adata.obs.columns:
            if "condition" in adata.obs.columns:
                adata.obs["sample_id"] = (
                    adata.obs["condition"].astype(str)
                    + "_"
                    + adata.obs["patient_id"].astype(str)
                )
            else:
                adata.obs["sample_id"] = adata.obs["patient_id"].astype(str)

        # Safety: sample_id must not mix conditions.
        if "condition" in adata.obs.columns and "sample_id" in adata.obs.columns:
            mixed_per_sample = (
                adata.obs.groupby("sample_id")["condition"]
                .nunique()
            )

            if len(mixed_per_sample) > 0 and mixed_per_sample.max() > 1:
                adata.obs["sample_id"] = (
                    adata.obs["condition"].astype(str)
                    + "_"
                    + adata.obs["sample_id"].astype(str)
                )

        # -----------------------------
        # Batch column
        # -----------------------------
        batch_col = self._profile_get(profile, "batch_column", None)

        if batch_col is None:
            for col in [
                "batch_id",
                "batch",
                "library",
                "library_prep_batch",
            ]:
                if col in adata.obs.columns:
                    batch_col = col
                    break

        if batch_col is not None and batch_col in adata.obs.columns:
            adata.obs["batch_id"] = adata.obs[batch_col].astype(str)
        elif "batch_id" not in adata.obs.columns:
            adata.obs["batch_id"] = "batch_1"

        # -----------------------------
        # Gene symbols
        # -----------------------------
        if "gene_symbol" not in adata.var.columns:
            if "feature_name" in adata.var.columns:
                adata.var["gene_symbol"] = adata.var["feature_name"].astype(str)
            elif "name" in adata.var.columns:
                adata.var["gene_symbol"] = adata.var["name"].astype(str)
            else:
                adata.var["gene_symbol"] = adata.var_names.astype(str)

        try:
            adata.var_names = adata.var["gene_symbol"].astype(str)
            adata.var_names_make_unique()
        except Exception:
            pass

        return adata

    # ------------------------------------------------------------------
    # Planning
    # ------------------------------------------------------------------
    def _choose_baseline(self, groups: List[str]) -> Optional[str]:
        if not groups:
            return None

        preferred = [
            "Healthy",
            "Control",
            "Normal",
            "Untreated",
            "Vehicle",
            "WT",
        ]

        for p in preferred:
            if p in groups:
                return p

        lowered = {g.lower(): g for g in groups}
        for p in preferred:
            if p.lower() in lowered:
                return lowered[p.lower()]

        return groups[0]

    def _make_plan(self, adata, profile=None) -> AnalysisPlan:
        has_conditions = (
            "condition" in adata.obs.columns
            and adata.obs["condition"].nunique() >= 2
        )

        condition_groups = []

        if has_conditions:
            condition_groups = (
                adata.obs["condition"]
                .astype(str)
                .value_counts()
                .index
                .tolist()
            )

        baseline = self._choose_baseline(condition_groups)

        pairwise_comparisons = []

        if baseline is not None:
            for group in condition_groups:
                if group != baseline:
                    pairwise_comparisons.append((baseline, group))

        use_pseudobulk = (
            has_conditions
            and "sample_id" in adata.obs.columns
            and adata.obs["sample_id"].nunique() >= 2
        )

        return AnalysisPlan(
            has_conditions=has_conditions,
            condition_column="condition" if has_conditions else None,
            cell_type_column="cell_type" if "cell_type" in adata.obs.columns else None,
            sample_column="sample_id" if "sample_id" in adata.obs.columns else None,
            patient_column="patient_id" if "patient_id" in adata.obs.columns else None,
            batch_column="batch_id" if "batch_id" in adata.obs.columns else None,
            condition_groups=condition_groups,
            baseline=baseline,
            pairwise_comparisons=pairwise_comparisons,
            run_condition_de=has_conditions and len(pairwise_comparisons) > 0,
            use_pseudobulk=use_pseudobulk,
            pseudobulk_key="sample_id" if use_pseudobulk else None,
        )

    # ------------------------------------------------------------------
    # Summary helpers
    # ------------------------------------------------------------------
    def _celltype_comparison(self, adata):
        if "cell_type" not in adata.obs.columns:
            return None

        if "condition" not in adata.obs.columns:
            return None

        try:
            return pd.crosstab(
                adata.obs["cell_type"].astype(str),
                adata.obs["condition"].astype(str),
                normalize="columns",
            )
        except Exception:
            return None

    def _count_sig_genes(self, de_results: pd.DataFrame) -> int:
        if de_results is None or de_results.empty:
            return 0

        if "pvals_adj" not in de_results.columns:
            return 0

        if "logfoldchanges" not in de_results.columns:
            return 0

        df = de_results.copy()
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")

        sig = df[
            (df["pvals_adj"] < 0.05)
            & (df["logfoldchanges"].abs() >= 0.5)
        ]

        return int(sig.shape[0])

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------
    def _basic_de_interpretation(self, de_results: pd.DataFrame) -> List[str]:
        """
        Interpret DE results using the full DE table.

        Positive logfoldchanges mean the case condition is higher than
        the reference condition.
        """

        if de_results is None or de_results.empty:
            return ["No differential expression results available."]

        hits = find_signature_hits(
            de_results,
            gene_col="names",
            lfc_col="logfoldchanges",
            padj_col="pvals_adj",
            padj_cutoff=0.05,
            lfc_cutoff=0.5,
        )

        lines = format_signature_hits(hits, max_genes_per_signature=15)

        if lines:
            return lines

        df = de_results.copy()

        required = {"names", "logfoldchanges", "pvals_adj"}

        if not required.issubset(set(df.columns)):
            return ["Differential expression table is missing required columns."]

        df["names"] = df["names"].astype(str)
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

        sig_up = df[
            (df["pvals_adj"] < 0.05)
            & (df["logfoldchanges"] > 0.5)
        ].copy()

        if sig_up.empty:
            return ["No strong condition-specific pathway signature identified."]

        top_up = (
            sig_up.sort_values(
                ["pvals_adj", "logfoldchanges"],
                ascending=[True, False],
            )
            ["names"]
            .astype(str)
            .head(8)
            .tolist()
        )

        return [
            "Differential expression detected, but no predefined pathway "
            f"signature dominated the ranked genes. Top increased genes include: {', '.join(top_up)}."
        ]

    # ------------------------------------------------------------------
    # Fallback DE
    # ------------------------------------------------------------------
    def _cell_level_fallback_de(
        self,
        adata,
        group1: str,
        group2: str,
        condition_key: str = "condition",
        n_genes: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Cell-level DE fallback.

        This is not preferred for real condition analysis, but it prevents
        complete failure for tiny cell types with too few pseudobulk samples.

        Direction:
            group1 = reference
            group2 = case

        Positive logfoldchanges = group2 higher than group1.
        """

        tmp = adata.copy()

        if condition_key not in tmp.obs.columns:
            raise ValueError(f"Missing condition column: {condition_key}")

        tmp.obs[condition_key] = tmp.obs[condition_key].astype(str)

        pair_mask = tmp.obs[condition_key].isin([group1, group2])
        tmp = tmp[pair_mask].copy()

        if tmp.n_obs < 2:
            raise ValueError("Too few cells for fallback DE.")

        counts = tmp.obs[condition_key].value_counts().to_dict()

        if counts.get(group1, 0) < 2 or counts.get(group2, 0) < 2:
            raise ValueError("Too few cells in one group for fallback DE.")

        tmp.obs[condition_key] = tmp.obs[condition_key].astype("category")

        sc.tl.rank_genes_groups(
            tmp,
            groupby=condition_key,
            groups=[group2],
            reference=group1,
            method="wilcoxon",
            n_genes=n_genes if n_genes is not None else tmp.n_vars,
        )

        df = sc.get.rank_genes_groups_df(tmp, group=group2)

        if df is None or df.empty:
            return pd.DataFrame(
                columns=[
                    "names",
                    "scores",
                    "logfoldchanges",
                    "pvals",
                    "pvals_adj",
                ]
            )

        for col in ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]:
            if col not in df.columns:
                df[col] = np.nan

        df["names"] = df["names"].astype(str)
        df["scores"] = pd.to_numeric(df["scores"], errors="coerce")
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        df["pvals"] = pd.to_numeric(df["pvals"], errors="coerce")
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

        df = df.dropna(subset=["names", "logfoldchanges", "pvals_adj"])

        df["abs_logfoldchanges"] = df["logfoldchanges"].abs()
        df["case"] = group2
        df["reference"] = group1

        df = df.sort_values(
            ["pvals_adj", "abs_logfoldchanges"],
            ascending=[True, False],
            kind="mergesort",
        ).reset_index(drop=True)

        return df

    # ------------------------------------------------------------------
    # Pairwise analysis
    # ------------------------------------------------------------------
    def _run_pairwise_analysis(
        self,
        adata,
        group1: str,
        group2: str,
        profile=None,
        cell_type: Optional[str] = None,
        n_genes: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Run condition DE for one comparison.

        group1 = baseline/reference
        group2 = case/test

        Positive logfoldchanges mean group2 higher than group1.
        """

        condition_key = "condition"
        sample_key = "sample_id"

        result = {
            "comparison": (group1, group2),
            "reference": group1,
            "case": group2,
            "mode": None,
            "status": "error",
            "error": None,
            "de_results": None,
            "pseudobulk_adata": None,
            "n_sig_genes": 0,
            "condition_counts": {},
            "interpretation": [],
        }

        if condition_key not in adata.obs.columns:
            result["error"] = "Missing condition column."
            result["interpretation"] = ["Condition metadata not available."]
            return result

        # Filter to comparison pair.
        pair_mask = adata.obs[condition_key].astype(str).isin([group1, group2])
        adata_pair = adata[pair_mask].copy()

        # Optional cell-type-specific subset.
        if cell_type is not None and "cell_type" in adata_pair.obs.columns:
            ct_mask = adata_pair.obs["cell_type"].astype(str) == str(cell_type)
            adata_pair = adata_pair[ct_mask].copy()

        if adata_pair.n_obs == 0:
            result["error"] = "No cells available for this comparison."
            result["interpretation"] = ["No cells available for this comparison."]
            return result

        result["condition_counts"] = (
            adata_pair.obs[condition_key]
            .astype(str)
            .value_counts()
            .to_dict()
        )

        if group1 not in result["condition_counts"] or group2 not in result["condition_counts"]:
            result["error"] = "One comparison group is absent after filtering."
            result["interpretation"] = [
                "One comparison group is absent after filtering."
            ]
            return result

        # --------------------------------------------------
        # Preferred path: pseudobulk DE
        # --------------------------------------------------
        try:
            pb_result = self.pseudobulk_de.pseudobulk(
                adata_pair,
                sample_key=sample_key,
                condition_key=condition_key,
                layer="counts" if "counts" in adata_pair.layers else None,
                min_cells_per_sample=10,
            )

            if pb_result.status != "ok":
                raise RuntimeError(pb_result.message)

            pb = pb_result.pseudobulk_adata

            de_results = self.pseudobulk_de.top_genes(
                pb_adata=pb,
                group=group2,          # case
                reference=group1,      # baseline/reference
                condition_key=condition_key,
                n_genes=None,          # keep full table
            )

            de_results = ensure_de_gene_symbols(de_results, adata_pair)

            result["mode"] = "pseudobulk"
            result["status"] = "ok"
            result["error"] = None
            result["de_results"] = de_results
            result["pseudobulk_adata"] = pb

            print(f"\nCondition DE completed: {group1} vs {group2}")

            return self._finalize_pairwise_result_dict(result)

        except Exception as e:
            pseudobulk_error = str(e)

            # If pseudobulk fails because there are too few biological
            # replicates, do NOT fall back to cell-level DE. Cell-level DE
            # would treat cells as independent replicates and can produce
            # misleading effect sizes.
            low_replicate_errors = [
                "Need at least 2 pseudobulk samples per group",
                "Too few pseudobulk samples",
            ]

            if any(msg in pseudobulk_error for msg in low_replicate_errors):
                result["mode"] = "skipped_low_replicates"
                result["status"] = "skipped"
                result["error"] = f"Pseudobulk skipped: {pseudobulk_error}"
                result["de_results"] = None
                result["pseudobulk_adata"] = None
                result["n_sig_genes"] = 0
                result["interpretation"] = [
                    "Differential expression skipped because this comparison has too few biological replicates for pseudobulk testing."
                ]
                return result

        # --------------------------------------------------
        # Fallback path: cell-level DE
        # --------------------------------------------------
        try:
            de_results = self._cell_level_fallback_de(
                adata_pair,
                group1=group1,
                group2=group2,
                condition_key=condition_key,
                n_genes=n_genes,
            )
            de_results = ensure_de_gene_symbols(de_results, adata_pair)
            result["mode"] = "cell_level_fallback"
            result["status"] = "ok"
            result["error"] = f"Pseudobulk failed: {pseudobulk_error}"
            result["de_results"] = de_results
            result["pseudobulk_adata"] = None

            print(f"\nCondition DE completed: {group1} vs {group2}")

            return self._finalize_pairwise_result_dict(result)

        except Exception as e:
            result["mode"] = None
            result["status"] = "error"
            result["error"] = (
                f"Pseudobulk failed: {pseudobulk_error}; "
                f"cell-level fallback failed: {e}"
            )
            result["interpretation"] = ["Differential expression failed."]
            return result

    # ------------------------------------------------------------------
    # Main router entry point
    # ------------------------------------------------------------------
    def run(self, adata, profile=None) -> Dict[str, Any]:
        adata = self._standardize_metadata(adata, profile=profile)
        plan = self._make_plan(adata, profile=profile)

        results = {
            "adata": adata,
            "profile": profile,
            "plan": plan,
            "condition_de_results": {},
            "condition_de_interpretation": {},
            "celltype_comparison": None,
            "celltype_specific": {},
        }

        results["celltype_comparison"] = self._celltype_comparison(adata)

        if not plan.run_condition_de:
            return results

        # --------------------------------------------------
        # Whole-dataset condition DE
        # --------------------------------------------------
        for group1, group2 in plan.pairwise_comparisons:
            result = self._run_pairwise_analysis(
                adata,
                group1=group1,
                group2=group2,
                profile=profile,
                cell_type=None,
                n_genes=None,
            )

            results["condition_de_results"][(group1, group2)] = result
            results["condition_de_interpretation"][(group1, group2)] = result.get(
                "interpretation",
                [],
            )

        # --------------------------------------------------
        # Cell-type-specific DE
        # --------------------------------------------------
        if "cell_type" in adata.obs.columns:
            cell_types = (
                adata.obs["cell_type"]
                .astype(str)
                .value_counts()
                .index
                .tolist()
            )

            for cell_type in cell_types:
                results["celltype_specific"][cell_type] = {}

                cell_mask = adata.obs["cell_type"].astype(str) == str(cell_type)
                adata_celltype = adata[cell_mask].copy()

                for group1, group2 in plan.pairwise_comparisons:
                    counts = (
                        adata_celltype.obs["condition"]
                        .astype(str)
                        .value_counts()
                        .to_dict()
                    )

                    # Skip extremely tiny comparisons.
                    if counts.get(group1, 0) < 20 or counts.get(group2, 0) < 20:
                        continue

                    result = self._run_pairwise_analysis(
                        adata,
                        group1=group1,
                        group2=group2,
                        profile=profile,
                        cell_type=cell_type,
                        n_genes=None,
                    )

                    results["celltype_specific"][cell_type][(group1, group2)] = result

        return results
