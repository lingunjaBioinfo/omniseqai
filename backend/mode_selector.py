from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any


@dataclass
class ModeDecision:
    mode: str
    reason: str


class ModeSelector:
    """
    Decide whether OmniSeqAI should run:
    - exploratory mode
    - condition mode
    - auto (choose based on metadata)
    """

    def choose(self, profile: Dict[str, Any], adata=None, mode: str = "auto") -> ModeDecision:
        if mode in {"exploratory", "condition"}:
            return ModeDecision(
                mode=mode,
                reason=f"User selected '{mode}' mode."
            )

        # auto mode
        condition_column = profile.get("condition_column")
        has_cell_type = bool(profile.get("cell_type_column"))
        has_sample = bool(profile.get("sample_column")) or bool(profile.get("patient_column"))

        if adata is not None and condition_column and condition_column in adata.obs.columns:
            groups = (
                adata.obs[condition_column]
                .astype(str)
                .replace({"nan": None, "None": None})
                .dropna()
                .value_counts()
                .index
                .tolist()
            )
            if len(groups) >= 2:
                return ModeDecision(
                    mode="condition",
                    reason="Detected a usable condition column with at least two groups."
                )

        if has_cell_type:
            return ModeDecision(
                mode="exploratory",
                reason="No usable condition comparison detected; defaulting to exploratory mode."
            )

        return ModeDecision(
            mode="exploratory",
            reason="Insufficient metadata for condition mode; defaulting to exploratory mode."
        )
