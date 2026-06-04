from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class ComparisonPlan:
    all_groups: List[str] = field(default_factory=list)
    baseline: str = ""
    comparisons: List[Tuple[str, str]] = field(default_factory=list)
    control_like_groups: List[str] = field(default_factory=list)


class ComparisonPolicy:
    """
    Decide which condition should be used as baseline/control
    and which pairwise comparisons should be run.

    Supports:
    - automatic control-like label detection
    - explicit user overrides
    """

    CONTROL_HINTS = [
        "healthy",
        "normal",
        "control",
        "hc",
        "untreated",
        "baseline",
        "mock",
        "vehicle",
        "sham",
    ]

    DISEASE_HINTS = [
        "covid",
        "infect",
        "disease",
        "tumor",
        "cancer",
        "case",
        "treated",
        "treatment",
        "vaccin",
        "vax",
        "booster",
        "drug",
        "therapy",
    ]

    def _norm(self, value) -> str:
        return str(value).strip()

    def _score_control_like(self, label: str) -> int:
        s = label.lower()
        score = 0

        for hint in self.CONTROL_HINTS:
            if hint in s:
                score += 2

        if s in {"hc", "ctrl", "control", "healthy", "normal"}:
            score += 2

        return score

    def choose_baseline(self, labels: List[str]) -> str:
        unique_labels = [self._norm(x) for x in labels if self._norm(x) != ""]
        if not unique_labels:
            raise ValueError("No labels provided to ComparisonPolicy.")

        scored = []
        for idx, label in enumerate(unique_labels):
            scored.append((label, self._score_control_like(label), idx))

        control_like = [x for x in scored if x[1] > 0]
        if control_like:
            control_like.sort(key=lambda x: (-x[1], x[2]))
            return control_like[0][0]

        return unique_labels[0]

    def build_plan(
        self,
        condition_series: pd.Series,
        control: Optional[str] = None,
        case: Optional[str] = None,
        pairs: Optional[List[Tuple[str, str]]] = None,
    ) -> ComparisonPlan:
        """
        Build a comparison plan from a condition column.

        Parameters
        ----------
        control:
            Optional explicit baseline label.
        case:
            Optional explicit case label. If provided with control,
            one pair is created: (control, case).
        pairs:
            Optional explicit list of pairwise comparisons.
            Each pair must be (baseline, group).
        """
        values = (
            condition_series.astype(str)
            .replace({"nan": pd.NA, "None": pd.NA})
            .dropna()
            .map(self._norm)
        )
        values = [v for v in values if v != ""]
        groups_in_order = list(dict.fromkeys(values))

        if len(groups_in_order) < 2:
            raise ValueError("Need at least two condition groups to build a comparison plan.")

        control_like_groups = [
            g for g in groups_in_order if self._score_control_like(g) > 0
        ]

        if pairs is not None and len(pairs) > 0:
            baseline = pairs[0][0]
            return ComparisonPlan(
                all_groups=groups_in_order,
                baseline=baseline,
                comparisons=pairs,
                control_like_groups=control_like_groups,
            )

        if control is not None and case is not None:
            if control not in groups_in_order:
                raise ValueError(f"Control '{control}' not found in condition labels.")
            if case not in groups_in_order:
                raise ValueError(f"Case '{case}' not found in condition labels.")

            return ComparisonPlan(
                all_groups=groups_in_order,
                baseline=control,
                comparisons=[(control, case)],
                control_like_groups=control_like_groups,
            )

        baseline = self.choose_baseline(groups_in_order)
        comparisons = [(baseline, g) for g in groups_in_order if g != baseline]

        return ComparisonPlan(
            all_groups=groups_in_order,
            baseline=baseline,
            comparisons=comparisons,
            control_like_groups=control_like_groups,
        )

    def summarize(
        self,
        condition_series: pd.Series,
        control: Optional[str] = None,
        case: Optional[str] = None,
        pairs: Optional[List[Tuple[str, str]]] = None,
    ) -> Dict[str, object]:
        plan = self.build_plan(
            condition_series=condition_series,
            control=control,
            case=case,
            pairs=pairs,
        )
        return {
            "all_groups": plan.all_groups,
            "baseline": plan.baseline,
            "comparisons": plan.comparisons,
            "control_like_groups": plan.control_like_groups,
        }
