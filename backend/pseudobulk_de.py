from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import anndata as ad


@dataclass
class PseudobulkResult:
    pseudobulk_adata: Optional[ad.AnnData]
    status: str = "ok"
    message: str = ""


class PseudobulkDE:
    """
    Corrected pseudobulk differential expression engine for OmniSeqAI.

    Direction convention used everywhere:

        group1 = baseline / reference
        group2 = case / test condition

    Therefore:

        run(adata, group1="Healthy", group2="IFN_beta")

    means:

        positive logfoldchanges = IFN_beta higher than Healthy

    This fixes the previous DE direction/ranking bug.
    """

    def __init__(
        self,
        default_layer: Optional[str] = "counts",
        pseudocount: float = 1.0,
    ):
        self.default_layer = default_layer
        self.pseudocount = pseudocount

    # ------------------------------------------------------------------
    # Matrix helpers
    # ------------------------------------------------------------------
    def _get_matrix(self, adata, layer: Optional[str] = None):
        layer = layer if layer is not None else self.default_layer

        if layer is not None and layer in adata.layers:
            X = adata.layers[layer]
        else:
            X = adata.X

        if sp.issparse(X):
            return X.tocsr()

        return np.asarray(X)

    def _to_dense(self, X):
        if sp.issparse(X):
            return X.toarray()
        return np.asarray(X)

    def _clean_string_series(self, s: pd.Series) -> pd.Series:
        return s.astype(str).str.strip()

    def _benjamini_hochberg(self, pvals: np.ndarray) -> np.ndarray:
        pvals = np.asarray(pvals, dtype=float)
        pvals = np.where(np.isfinite(pvals), pvals, 1.0)

        n = len(pvals)
        order = np.argsort(pvals)
        ranked = pvals[order]

        adjusted = ranked * n / (np.arange(n) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0.0, 1.0)

        out = np.empty_like(adjusted)
        out[order] = adjusted

        return out

    # ------------------------------------------------------------------
    # Pseudobulk aggregation
    # ------------------------------------------------------------------
    def pseudobulk(
        self,
        adata,
        sample_key: str = "sample_id",
        condition_key: str = "condition",
        layer: Optional[str] = None,
        min_cells_per_sample: int = 10,
    ) -> PseudobulkResult:
        """
        Aggregate cells into sample-level pseudobulk profiles.

        Returned AnnData:
        - rows = pseudobulk samples
        - columns = genes
        - X = summed counts / expression values
        """

        if sample_key not in adata.obs.columns:
            return PseudobulkResult(
                pseudobulk_adata=None,
                status="failed",
                message=f"Missing sample column: {sample_key}",
            )

        if condition_key not in adata.obs.columns:
            return PseudobulkResult(
                pseudobulk_adata=None,
                status="failed",
                message=f"Missing condition column: {condition_key}",
            )

        obs = adata.obs.copy()
        obs[sample_key] = self._clean_string_series(obs[sample_key])
        obs[condition_key] = self._clean_string_series(obs[condition_key])

        # Safety check:
        # A pseudobulk sample must not mix multiple biological conditions.
        mixed = obs.groupby(sample_key)[condition_key].nunique()
        mixed_samples = mixed[mixed > 1].index.tolist()

        if mixed_samples:
            safe_sample_key = "__omniseqai_safe_sample_id"
            obs[safe_sample_key] = (
                obs[condition_key].astype(str)
                + "_"
                + obs[sample_key].astype(str)
            )
            sample_key_used = safe_sample_key
        else:
            sample_key_used = sample_key

        X = self._get_matrix(adata, layer=layer)

        groups = obs.groupby(sample_key_used, sort=False)

        bulk_rows = []
        bulk_obs = []
        bulk_names = []

        for sample_id, idx in groups.indices.items():
            idx = np.asarray(idx)

            if len(idx) < min_cells_per_sample:
                continue

            sample_conditions = obs.iloc[idx][condition_key].unique()

            if len(sample_conditions) != 1:
                continue

            condition = str(sample_conditions[0])

            summed = X[idx].sum(axis=0)
            summed = np.asarray(summed).ravel()

            bulk_rows.append(summed)
            bulk_names.append(str(sample_id))
            bulk_obs.append(
                {
                    "sample_id": str(sample_id),
                    "condition": condition,
                    "n_cells": int(len(idx)),
                }
            )

        if len(bulk_rows) < 2:
            return PseudobulkResult(
                pseudobulk_adata=None,
                status="failed",
                message="Too few pseudobulk samples after aggregation.",
            )

        bulk_X = np.vstack(bulk_rows)

        pb = ad.AnnData(
            X=bulk_X,
            obs=pd.DataFrame(bulk_obs, index=bulk_names),
            var=adata.var.copy(),
        )

        pb.var_names = adata.var_names.copy()

        pb.uns["pseudobulk_source"] = {
            "sample_key": sample_key,
            "sample_key_used": sample_key_used,
            "condition_key": condition_key,
            "layer": layer if layer is not None else self.default_layer,
            "min_cells_per_sample": min_cells_per_sample,
            "n_input_cells": int(adata.n_obs),
            "n_pseudobulk_samples": int(pb.n_obs),
        }

        return PseudobulkResult(
            pseudobulk_adata=pb,
            status="ok",
            message="Pseudobulk aggregation completed.",
        )

    # Compatibility aliases
    def create_pseudobulk(self, *args, **kwargs) -> PseudobulkResult:
        return self.pseudobulk(*args, **kwargs)

    def aggregate(self, *args, **kwargs) -> PseudobulkResult:
        return self.pseudobulk(*args, **kwargs)

    # ------------------------------------------------------------------
    # Comparison direction helpers
    # ------------------------------------------------------------------
    def _resolve_case_reference(
        self,
        pb_adata,
        group: str,
        reference: Optional[str],
        condition_key: str,
    ) -> Tuple[str, str]:
        conditions = (
            pb_adata.obs[condition_key]
            .astype(str)
            .str.strip()
            .unique()
            .tolist()
        )

        group = str(group).strip()

        if reference is not None:
            case = group
            reference = str(reference).strip()

            if case not in conditions:
                raise ValueError(
                    f"Case group not found: {case}. Available: {conditions}"
                )

            if reference not in conditions:
                raise ValueError(
                    f"Reference group not found: {reference}. Available: {conditions}"
                )

            return case, reference

        if group not in conditions:
            raise ValueError(
                f"Group not found: {group}. Available: {conditions}"
            )

        if len(conditions) != 2:
            raise ValueError(
                "reference must be supplied when more than two conditions are present."
            )

        # Backward-compatible OmniSeqAI behavior:
        # if only `group` is supplied, treat it as baseline/reference.
        reference = group
        case = [c for c in conditions if c != reference][0]

        return case, reference

    def _logcpm(self, counts: np.ndarray) -> np.ndarray:
        counts = np.asarray(counts, dtype=float)

        library_sizes = counts.sum(axis=1)
        library_sizes = np.where(library_sizes <= 0, 1.0, library_sizes)

        cpm = counts / library_sizes[:, None] * 1e6
        logcpm = np.log2(cpm + self.pseudocount)

        return logcpm

    # ------------------------------------------------------------------
    # Differential expression
    # ------------------------------------------------------------------
    def differential_expression(
        self,
        pb_adata,
        case: str,
        reference: str,
        condition_key: str = "condition",
        min_total_counts: float = 10.0,
        min_samples_expressed: int = 2,
    ) -> pd.DataFrame:
        """
        Run sample-level pseudobulk DE.

        Positive logfoldchanges mean higher expression in `case`
        compared with `reference`.
        """

        if pb_adata is None:
            raise ValueError("pb_adata is None.")

        if condition_key not in pb_adata.obs.columns:
            raise ValueError(
                f"Missing condition column in pseudobulk object: {condition_key}"
            )

        obs = pb_adata.obs.copy()
        obs[condition_key] = obs[condition_key].astype(str).str.strip()

        case = str(case).strip()
        reference = str(reference).strip()

        case_mask = obs[condition_key] == case
        ref_mask = obs[condition_key] == reference

        n_case = int(case_mask.sum())
        n_ref = int(ref_mask.sum())

        if n_case < 2 or n_ref < 2:
            raise ValueError(
                f"Need at least 2 pseudobulk samples per group. "
                f"Got {case}: {n_case}, {reference}: {n_ref}."
            )

        counts = self._to_dense(pb_adata.X).astype(float)

        expressed = (
            (counts.sum(axis=0) >= min_total_counts)
            & ((counts > 0).sum(axis=0) >= min_samples_expressed)
        )

        if int(expressed.sum()) == 0:
            raise ValueError("No genes passed pseudobulk expression filtering.")

        genes_all = np.asarray(pb_adata.var_names.astype(str))
        genes = genes_all[expressed]

        counts_f = counts[:, expressed]
        logcpm = self._logcpm(counts_f)

        case_values = logcpm[case_mask.values, :]
        ref_values = logcpm[ref_mask.values, :]

        mean_case = case_values.mean(axis=0)
        mean_ref = ref_values.mean(axis=0)

        logfc = mean_case - mean_ref

        test = stats.ttest_ind(
            case_values,
            ref_values,
            axis=0,
            equal_var=False,
            nan_policy="omit",
        )

        scores = np.asarray(test.statistic, dtype=float)
        pvals = np.asarray(test.pvalue, dtype=float)

        pvals = np.where(np.isfinite(pvals), pvals, 1.0)
        scores = np.where(np.isfinite(scores), scores, 0.0)

        pvals_adj = self._benjamini_hochberg(pvals)

        pct_case = (counts_f[case_mask.values, :] > 0).mean(axis=0)
        pct_ref = (counts_f[ref_mask.values, :] > 0).mean(axis=0)

        df = pd.DataFrame(
            {
                "names": genes,
                "scores": scores,
                "logfoldchanges": logfc,
                "pvals": pvals,
                "pvals_adj": pvals_adj,
                "mean_logcpm_case": mean_case,
                "mean_logcpm_reference": mean_ref,
                "pct_case": pct_case,
                "pct_reference": pct_ref,
                "case": case,
                "reference": reference,
                "n_case_samples": n_case,
                "n_reference_samples": n_ref,
            }
        )

        df["abs_logfoldchanges"] = df["logfoldchanges"].abs()

        df = df.sort_values(
            ["pvals_adj", "abs_logfoldchanges"],
            ascending=[True, False],
            kind="mergesort",
        ).reset_index(drop=True)

        return df

    def top_genes(
        self,
        pb_adata,
        group: str,
        reference: Optional[str] = None,
        condition_key: str = "condition",
        n_genes: Optional[int] = None,
        min_total_counts: float = 10.0,
        min_samples_expressed: int = 2,
    ) -> pd.DataFrame:
        """
        Backward-compatible DE entry point.

        If reference is provided:
            group = case
            reference = baseline

        If reference is None:
            group = baseline/reference
            other condition = case

        By default, this returns the full ranked DE table.
        """

        case, ref = self._resolve_case_reference(
            pb_adata=pb_adata,
            group=group,
            reference=reference,
            condition_key=condition_key,
        )

        df = self.differential_expression(
            pb_adata=pb_adata,
            case=case,
            reference=ref,
            condition_key=condition_key,
            min_total_counts=min_total_counts,
            min_samples_expressed=min_samples_expressed,
        )

        df.attrs["case"] = case
        df.attrs["reference"] = ref
        df.attrs["comparison"] = f"{case}_vs_{ref}"
        df.attrs["positive_logfc_means"] = f"{case} higher than {ref}"

                # Always return the full DE table.
        # Reporting code can decide how many rows to print.
        # This prevents biologically important genes ranked below top 100
        # from disappearing before interpretation/volcano/reporting.
        return df

    def run(
        self,
        adata,
        group: Optional[str] = None,
        reference: Optional[str] = None,
        group1: Optional[str] = None,
        group2: Optional[str] = None,
        sample_key: str = "sample_id",
        condition_key: str = "condition",
        layer: Optional[str] = None,
        min_cells_per_sample: int = 10,
        n_genes: Optional[int] = None,
        min_total_counts: float = 10.0,
        min_samples_expressed: int = 2,
        return_de: bool = False,
        **kwargs,
    ):
        """
        Router-compatible pseudobulk entry point.

        Default behavior:
            return PseudobulkResult

        This matches AnalysisRouter code that does:

            pb_result = self.pseudobulk_de.run(...)
            pb = pb_result.pseudobulk_adata

        Optional behavior:
            return_de=True returns the DE DataFrame directly.

        Direction convention:
            group1 = baseline/reference
            group2 = case/test
            positive logFC = group2 higher than group1
        """

        pb_result = self.pseudobulk(
            adata=adata,
            sample_key=sample_key,
            condition_key=condition_key,
            layer=layer,
            min_cells_per_sample=min_cells_per_sample,
        )

        if pb_result.status != "ok":
            raise RuntimeError(pb_result.message)

        # Default: return pseudobulk object for router compatibility
        if not return_de:
            return pb_result

        # Optional: return DE table directly
        if group1 is not None and group2 is not None:
            baseline = str(group1).strip()
            case = str(group2).strip()

            group = case
            reference = baseline

        elif group is not None and reference is not None:
            group = str(group).strip()
            reference = str(reference).strip()

        elif group is not None and reference is None:
            group = str(group).strip()

        else:
            raise ValueError(
                "For return_de=True, provide either group/reference or group1/group2."
            )

        de_results = self.top_genes(
            pb_adata=pb_result.pseudobulk_adata,
            group=group,
            reference=reference,
            condition_key=condition_key,
            n_genes=n_genes,
            min_total_counts=min_total_counts,
            min_samples_expressed=min_samples_expressed,
        )

        de_results.attrs["pseudobulk_adata"] = pb_result.pseudobulk_adata
        de_results.attrs["pseudobulk_status"] = pb_result.status
        de_results.attrs["pseudobulk_message"] = pb_result.message

        return de_results
