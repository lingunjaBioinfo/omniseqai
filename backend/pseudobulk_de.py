from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


@dataclass
class PseudobulkResult:
    pseudobulk_adata: sc.AnnData
    sample_key: str
    condition_key: str


class PseudobulkDE:
    """
    Create pseudobulk profiles by summing counts across cells from the same sample.
    Then run standard DE on the sample-level matrix.

    This is much more statistically appropriate than treating cells as independent
    replicates when sample/donor information exists.
    """

    def _get_matrix(self, adata, layer: Optional[str] = None):
        if layer is not None and layer in adata.layers:
            X = adata.layers[layer]
        else:
            X = adata.X

        if sp.issparse(X):
            return X.tocsr()
        return sp.csr_matrix(X)

    def _normalize_label(self, value) -> str:
        return str(value).strip()

    def build_pseudobulk(
        self,
        adata,
        sample_key: str = "sample_id",
        condition_key: str = "condition",
        layer: Optional[str] = None,
        min_cells_per_sample: int = 1,
    ) -> sc.AnnData:
        """
        Sum expression per sample. The sample-level condition is assigned by majority
        vote / mode across cells from that sample.
        """
        if sample_key not in adata.obs.columns:
            raise ValueError(f"'{sample_key}' not found in adata.obs")

        if condition_key not in adata.obs.columns:
            raise ValueError(f"'{condition_key}' not found in adata.obs")

        obs = adata.obs.copy()
        obs[sample_key] = obs[sample_key].astype(str).map(self._normalize_label)
        obs[condition_key] = obs[condition_key].astype(str).map(self._normalize_label)

        # Keep only samples with enough cells
        sample_sizes = obs[sample_key].value_counts()
        keep_samples = sample_sizes[sample_sizes >= min_cells_per_sample].index.tolist()
        obs = obs[obs[sample_key].isin(keep_samples)].copy()

        if obs.empty:
            raise ValueError("No samples passed the minimum cell count filter.")

        X = self._get_matrix(adata[obs.index], layer=layer)
        sample_ids = obs[sample_key].astype(str).tolist()
        conditions = obs[condition_key].astype(str).tolist()

        unique_samples = list(dict.fromkeys(sample_ids))
        rows = []
        out_obs = []

        for sample in unique_samples:
            idx = np.where(np.array(sample_ids) == sample)[0]
            if len(idx) == 0:
                continue

            sample_matrix = X[idx].sum(axis=0)
            if sp.issparse(sample_matrix):
                sample_matrix = sample_matrix.tocsr()
            else:
                sample_matrix = sp.csr_matrix(sample_matrix)

            rows.append(sample_matrix)

            sample_conditions = [conditions[i] for i in idx]
            condition_mode = pd.Series(sample_conditions).mode()
            if len(condition_mode) == 0:
                condition_value = sample_conditions[0]
            else:
                condition_value = str(condition_mode.iloc[0])

            out_obs.append(
                {
                    sample_key: sample,
                    condition_key: condition_value,
                    "n_cells": int(len(idx)),
                }
            )

        pbX = sp.vstack(rows).tocsr()
        pb_obs = pd.DataFrame(out_obs).set_index(sample_key)
        pb_var = adata.var.copy()

        pb = sc.AnnData(X=pbX, obs=pb_obs, var=pb_var)

        # Preserve useful metadata
        pb.uns["pseudobulk_source"] = {
            "sample_key": sample_key,
            "condition_key": condition_key,
            "layer": layer,
            "n_input_cells": int(adata.n_obs),
            "n_pseudobulk_samples": int(pb.n_obs),
        }

        return pb

    def run(
        self,
        adata,
        sample_key: str = "sample_id",
        condition_key: str = "condition",
        group1: str = "Disease",
        group2: str = "Healthy",
        layer: Optional[str] = None,
        min_cells_per_sample: int = 10,
    ) -> PseudobulkResult:
        """
        Build pseudobulk and run DE between group1 and group2 at sample level.
        """
        pb = self.build_pseudobulk(
            adata,
            sample_key=sample_key,
            condition_key=condition_key,
            layer=layer,
            min_cells_per_sample=min_cells_per_sample,
        )

        # Normalize pseudobulk sample-level counts
        sc.pp.normalize_total(pb, target_sum=10000)
        sc.pp.log1p(pb)

        # Make sure the requested groups exist
        groups = pb.obs[condition_key].astype(str).unique().tolist()
        if group1 not in groups or group2 not in groups:
            raise ValueError(
                f"Requested groups not found in pseudobulk data. "
                f"Available: {groups}"
            )

        sc.tl.rank_genes_groups(
            pb,
            groupby=condition_key,
            groups=[group1],
            reference=group2,
            method="wilcoxon",
            use_raw=False,
        )

        return PseudobulkResult(
            pseudobulk_adata=pb,
            sample_key=sample_key,
            condition_key=condition_key,
        )

    def _pretty_gene_names(self, adata, genes: pd.Series) -> pd.Series:
        if genes is None:
            return genes

        out = genes.astype(str).copy()

        for col in ("gene_symbol", "feature_name"):
            if col not in adata.var.columns:
                continue

            mapping = adata.var[col].astype(str).to_dict()
            mapped = out.map(mapping)

            mask = (
                mapped.notna()
                & (mapped != "")
                & (mapped != "nan")
                & ~mapped.str.startswith("ENSG", na=False)
                & ~mapped.str.startswith("NCBITaxon:", na=False)
            )

            out.loc[mask] = mapped.loc[mask]

        return out

    def top_genes(
        self,
        pseudobulk_adata: sc.AnnData,
        group: str = "Disease",
        n_genes: int = 20,
    ) -> pd.DataFrame:
        df = sc.get.rank_genes_groups_df(
            pseudobulk_adata,
            group=group
        )

        if df is None or df.empty:
            return df

        df = df.copy()
        df["names"] = self._pretty_gene_names(
            pseudobulk_adata,
            df["names"]
        )

        for col in ["logfoldchanges", "pvals_adj", "scores"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df = df.sort_values(
            ["pvals_adj", "logfoldchanges"],
            ascending=[True, False],
            kind="mergesort"
        )

        return df.head(n_genes).reset_index(drop=True)

