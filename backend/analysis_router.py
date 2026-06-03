from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from backend.metadata_detector import MetadataDetector
from backend.gene_mapper import GeneMapper
from backend.condition_comparison import ConditionComparison
from backend.condition_de import ConditionDE
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
    run_celltype_comparison: bool = False
    run_condition_de: bool = False
    run_celltype_specific_de: bool = False


class AnalysisRouter:

    def __init__(self):
        self.detector = MetadataDetector()
        self.gene_mapper = GeneMapper()
        self.condition_comparison = ConditionComparison()
        self.condition_de = ConditionDE()
        self.condition_pathways = ConditionPathways()
        self.condition_interpreter = ConditionInterpreter()

    def inspect(self, adata) -> Dict[str, Any]:
        return self.detector.detect(adata)

    def standardize(self, adata, profile: Optional[Dict[str, Any]] = None):
        if profile is None:
            profile = self.detector.detect(adata)

        adata = self.detector.apply(adata, profile)
        adata = self.gene_mapper.fix_gene_names(adata)

        return adata, profile

    def build_plan(self, adata, profile: Optional[Dict[str, Any]] = None) -> AnalysisPlan:
        if profile is None:
            profile = self.detector.detect(adata)

        plan = AnalysisPlan(
            cell_type_column=profile.get("cell_type_column"),
            condition_column=profile.get("condition_column"),
            sample_column=profile.get("sample_column"),
            patient_column=profile.get("patient_column"),
            batch_column=profile.get("batch_column"),
        )

        # IMPORTANT:
        # Build groups from the standardized condition column if it exists.
        if "condition" in adata.obs.columns:
            values = (
                adata.obs["condition"]
                .astype(str)
                .replace({"nan": pd.NA})
                .dropna()
                .value_counts()
                .index
                .tolist()
            )
        else:
            values = []

        values = [v for v in values if str(v).strip() != ""]
        plan.has_conditions = len(values) >= 2
        plan.condition_groups = values

        if len(values) >= 2:
            plan.run_celltype_comparison = True
            plan.run_condition_de = True
            plan.run_celltype_specific_de = True

            # Use the two most common groups for whole-dataset DE.
            if len(values) == 2:
                plan.pairwise_comparisons = [(values[0], values[1])]
            else:
                plan.pairwise_comparisons = [(values[0], values[1])]

        return plan

    def run(
        self,
        adata,
        profile: Optional[Dict[str, Any]] = None,
        run_celltype_specific_de: bool = True,
        min_cells_per_group: int = 20,
        min_cells_per_celltype_group: int = 50,
    ) -> Dict[str, Any]:
        adata, profile = self.standardize(adata, profile)
        plan = self.build_plan(adata, profile)

        results: Dict[str, Any] = {
            "profile": profile,
            "plan": plan,
            "adata": adata,
            "celltype_comparison": None,
            "condition_de_results": {},
            "celltype_specific": {},
        }

        # Cell-type comparison across condition
        if plan.run_celltype_comparison:
            try:
                results["celltype_comparison"] = self.condition_comparison.compare_celltypes(
                    adata,
                    condition_column="condition",
                    celltype_column="cell_type"
                )
            except Exception as e:
                results["celltype_comparison_error"] = str(e)

        # Whole-dataset condition DE
        if plan.run_condition_de and len(plan.condition_groups) >= 2:
            group1, group2 = plan.condition_groups[:2]

            counts = adata.obs["condition"].value_counts()
            if counts.get(group1, 0) >= min_cells_per_group and counts.get(group2, 0) >= min_cells_per_group:
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
                        n_genes=100
                    )

                    results["condition_de_results"][(group1, group2)] = de_results

                    if de_results is not None and not de_results.empty:
                        pathways = self.condition_pathways.analyze(de_results, top_n=50)
                        results["condition_de_pathways"] = pathways

                        if pathways is not None and not pathways.empty:
                            results["condition_de_interpretation"] = self.condition_interpreter.interpret(pathways)
                        else:
                            results["condition_de_interpretation"] = [
                                "No strong condition-specific pathway signature identified."
                            ]
                    else:
                        results["condition_de_interpretation"] = [
                            "No DE results available."
                        ]

                except Exception as e:
                    results["condition_de_error"] = str(e)

        # Cell-type-specific DE
        if run_celltype_specific_de and plan.run_celltype_specific_de:
            celltype_results: Dict[str, Any] = {}

            for cell_type in adata.obs["cell_type"].astype(str).value_counts().index:
                subset = adata[adata.obs["cell_type"].astype(str) == cell_type].copy()

                if "condition" not in subset.obs.columns:
                    continue

                cond_counts = subset.obs["condition"].value_counts()
                if len(cond_counts) < 2:
                    continue

                groups = cond_counts.index.tolist()
                g1, g2 = groups[0], groups[1]

                if cond_counts.get(g1, 0) < min_cells_per_celltype_group or cond_counts.get(g2, 0) < min_cells_per_celltype_group:
                    continue

                try:
                    de_adata = self.condition_de.run(
                        subset,
                        condition_key="condition",
                        group1=g1,
                        group2=g2
                    )

                    de_results = self.condition_de.top_genes(
                        de_adata,
                        group=g1,
                        n_genes=100
                    )

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

                    celltype_results[cell_type] = {
                        "condition_pair": (g1, g2),
                        "condition_counts": cond_counts.to_dict(),
                        "de_results": de_results,
                        "significant_genes": sig,
                        "pathways": pathways,
                        "interpretation": interpretation,
                    }

                except Exception as e:
                    celltype_results[cell_type] = {
                        "error": str(e),
                        "condition_counts": cond_counts.to_dict(),
                    }

            results["celltype_specific"] = celltype_results

        return results
