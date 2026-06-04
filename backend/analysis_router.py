from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from backend.metadata_detector import MetadataDetector
from backend.gene_mapper import GeneMapper
from backend.comparison_policy import ComparisonPolicy
from backend.condition_comparison import ConditionComparison
from backend.condition_de import ConditionDE
from backend.pseudobulk_de import PseudobulkDE
from backend.condition_pathways import ConditionPathways
from backend.condition_interpreter import ConditionInterpreter


@dataclass
class AnalysisPlan:
    has_conditions: bool = False
    condition_column: Optional[str] = None
    cell_type_column: Optional[str] = None
    sample_column: Optional[str] = None
    patient_column: Optional[str] = None
    batch_column: Optional[str] = None
    condition_groups: List[str] = field(default_factory=list)
    pairwise_comparisons: List[Tuple[str, str]] = field(default_factory=list)
    baseline: str = ""
    run_celltype_comparison: bool = False
    run_condition_de: bool = False
    run_celltype_specific_de: bool = False
    use_pseudobulk: bool = False
    pseudobulk_key: Optional[str] = None


class AnalysisRouter:

    def __init__(self):
        self.detector = MetadataDetector()
        self.gene_mapper = GeneMapper()
        self.comparison_policy = ComparisonPolicy()
        self.condition_comparison = ConditionComparison()
        self.condition_de = ConditionDE()
        self.pseudobulk_de = PseudobulkDE()
        self.condition_pathways = ConditionPathways()
        self.condition_interpreter = ConditionInterpreter()

    def inspect(self, adata) -> Dict[str, Any]:
        return self.detector.detect(adata)

    def standardize(
        self,
        adata,
        profile: Optional[Dict[str, Any]] = None
    ):
        if profile is None:
            profile = self.detector.detect(adata)

        adata = self.detector.apply(adata, profile)
        adata = self.gene_mapper.fix_gene_names(adata)

        return adata, profile

    def _choose_sample_key(self, adata) -> Optional[str]:
        for key in ("sample_id", "patient_id"):
            if key in adata.obs.columns:
                return key
        return None

    def build_plan(
        self,
        adata,
        profile: Optional[Dict[str, Any]] = None
    ) -> AnalysisPlan:
        if profile is None:
            profile = self.detector.detect(adata)

        plan = AnalysisPlan(
            cell_type_column=profile.get("cell_type_column"),
            condition_column=profile.get("condition_column"),
            sample_column=profile.get("sample_column"),
            patient_column=profile.get("patient_column"),
            batch_column=profile.get("batch_column"),
        )

        if "condition" in adata.obs.columns:
            try:
                comp = self.comparison_policy.summarize(
                    adata.obs["condition"]
                )
                plan.condition_groups = comp["all_groups"]
                plan.pairwise_comparisons = comp["comparisons"]
                plan.baseline = comp["baseline"]
                plan.has_conditions = len(plan.condition_groups) >= 2
                plan.run_celltype_comparison = plan.has_conditions
                plan.run_condition_de = plan.has_conditions
                plan.run_celltype_specific_de = plan.has_conditions
            except Exception:
                pass

        plan.use_pseudobulk = self._choose_sample_key(adata) is not None
        plan.pseudobulk_key = self._choose_sample_key(adata)

        return plan

    def _run_pairwise_analysis(
        self,
        adata,
        group1: str,
        group2: str,
        sample_key: Optional[str] = None,
        min_cells_per_group: int = 20,
        min_cells_per_sample: int = 10,
        n_genes: int = 100
    ) -> Dict[str, Any]:

        counts = adata.obs["condition"].value_counts()

        if counts.get(group1, 0) < min_cells_per_group or counts.get(group2, 0) < min_cells_per_group:
            return {
                "status": "skipped",
                "error": (
                    f"Not enough cells for {group1} vs {group2}. "
                    f"Counts: {counts.to_dict()}"
                ),
                "mode": None,
                "de_results": None,
                "pathways": None,
                "interpretation": [],
                "n_sig_genes": 0,
            }

        # Prefer pseudobulk when sample-like metadata exists.
        if sample_key is not None and sample_key in adata.obs.columns:
            try:
                pb_result = self.pseudobulk_de.run(
                    adata,
                    sample_key=sample_key,
                    condition_key="condition",
                    group1=group1,
                    group2=group2,
                    min_cells_per_sample=min_cells_per_sample,
                )

                pb = pb_result.pseudobulk_adata
                de_results = self.pseudobulk_de.top_genes(
                    pb,
                    group=group1,
                    n_genes=n_genes
                )

                mode = "pseudobulk"
            except Exception as e:
                # Fall back to cell-level DE if pseudobulk fails.
                try:
                    de_adata = self.condition_de.run(
                        adata.copy(),
                        condition_key="condition",
                        group1=group1,
                        group2=group2
                    )
                    de_results = self.condition_de.top_genes(
                        de_adata,
                        group=group1,
                        n_genes=n_genes
                    )
                    mode = "cell_level_fallback"
                    return self._finalize_pairwise_result(
                        group1,
                        group2,
                        mode,
                        de_results,
                        error=f"Pseudobulk failed: {e}"
                    )
                except Exception as e2:
                    return {
                        "status": "error",
                        "error": f"Pseudobulk failed: {e}; cell-level fallback failed: {e2}",
                        "mode": None,
                        "de_results": None,
                        "pathways": None,
                        "interpretation": [],
                        "n_sig_genes": 0,
                    }

            return self._finalize_pairwise_result(
                group1,
                group2,
                mode,
                de_results
            )

        # Otherwise use cell-level DE.
        try:
            de_adata = self.condition_de.run(
                adata.copy(),
                condition_key="condition",
                group1=group1,
                group2=group2
            )
            de_results = self.condition_de.top_genes(
                de_adata,
                group=group1,
                n_genes=n_genes
            )

            return self._finalize_pairwise_result(
                group1,
                group2,
                "cell_level",
                de_results
            )

        except Exception as e:
            return {
                "status": "error",
                "error": str(e),
                "mode": "cell_level",
                "de_results": None,
                "pathways": None,
                "interpretation": [],
                "n_sig_genes": 0,
            }

    def _finalize_pairwise_result(
        self,
        group1: str,
        group2: str,
        mode: str,
        de_results,
        error: Optional[str] = None
    ) -> Dict[str, Any]:

        if de_results is None or de_results.empty:
            return {
                "status": "ok",
                "mode": mode,
                "error": error,
                "de_results": de_results,
                "pathways": None,
                "interpretation": ["No DE results available."],
                "n_sig_genes": 0,
            }

        sig = de_results[
            (de_results["pvals_adj"] < 0.05) &
            (de_results["logfoldchanges"].abs() > 0.5)
        ].copy()

        pathways = self.condition_pathways.analyze(
            sig if not sig.empty else de_results,
            top_n=50
        )

        if pathways is not None and not pathways.empty:
            interpretation = self.condition_interpreter.interpret(pathways)
        else:
            interpretation = [
                "No strong condition-specific pathway signature identified."
            ]

        return {
            "status": "ok",
            "mode": mode,
            "error": error,
            "de_results": de_results,
            "significant_genes": sig,
            "pathways": pathways,
            "interpretation": interpretation,
            "n_sig_genes": 0 if sig.empty else int(len(sig)),
        }

    def run(
        self,
        adata,
        profile: Optional[Dict[str, Any]] = None,
        run_celltype_specific_de: bool = True,
        min_cells_per_group: int = 20,
        min_cells_per_celltype_group: int = 50,
        min_cells_per_sample: int = 10,
    ) -> Dict[str, Any]:

        adata, profile = self.standardize(adata, profile)
        plan = self.build_plan(adata, profile)

        results: Dict[str, Any] = {
            "profile": profile,
            "plan": plan,
            "adata": adata,
            "celltype_comparison": None,
            "condition_de_results": {},
            "condition_de_interpretation": {},
            "celltype_specific": {},
        }

        # Cell-type proportions across condition
        if plan.run_celltype_comparison:
            try:
                results["celltype_comparison"] = self.condition_comparison.compare_celltypes(
                    adata,
                    condition_column="condition",
                    celltype_column="cell_type"
                )
            except Exception as e:
                results["celltype_comparison_error"] = str(e)

        sample_key = plan.pseudobulk_key if plan.use_pseudobulk else None

        # Whole-dataset pairwise comparisons
        if plan.run_condition_de and len(plan.pairwise_comparisons) > 0:
            for group1, group2 in plan.pairwise_comparisons:
                result = self._run_pairwise_analysis(
                    adata,
                    group1=group1,
                    group2=group2,
                    sample_key=sample_key,
                    min_cells_per_group=min_cells_per_group,
                    min_cells_per_sample=min_cells_per_sample,
                    n_genes=100
                )

                results["condition_de_results"][(group1, group2)] = result
                results["condition_de_interpretation"][(group1, group2)] = result.get("interpretation", [])

        # Cell-type-specific comparisons
        if run_celltype_specific_de and plan.run_celltype_specific_de:
            celltype_results: Dict[str, Any] = {}

            for cell_type in adata.obs["cell_type"].astype(str).value_counts().index:
                subset = adata[
                    adata.obs["cell_type"].astype(str) == cell_type
                ].copy()

                cond_counts = subset.obs["condition"].value_counts()
                if len(cond_counts) < 2:
                    continue

                celltype_results[cell_type] = {}

                # Run each planned comparison inside this cell type.
                for group1, group2 in plan.pairwise_comparisons:
                    if cond_counts.get(group1, 0) < min_cells_per_celltype_group or cond_counts.get(group2, 0) < min_cells_per_celltype_group:
                        celltype_results[cell_type][(group1, group2)] = {
                            "status": "skipped",
                            "error": f"Not enough cells in {cell_type} for {group1} vs {group2}.",
                            "condition_counts": cond_counts.to_dict(),
                        }
                        continue

                    result = self._run_pairwise_analysis(
                        subset,
                        group1=group1,
                        group2=group2,
                        sample_key=sample_key,
                        min_cells_per_group=min_cells_per_celltype_group,
                        min_cells_per_sample=min_cells_per_sample,
                        n_genes=100
                    )

                    result["condition_counts"] = cond_counts.to_dict()
                    celltype_results[cell_type][(group1, group2)] = result

            results["celltype_specific"] = celltype_results

        return results
