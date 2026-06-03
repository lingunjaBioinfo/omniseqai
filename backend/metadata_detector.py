from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import pandas as pd

from backend.gene_mapper import GeneMapper


@dataclass
class MetadataProfile:
    cell_type_column: Optional[str] = None
    condition_column: Optional[str] = None
    sample_column: Optional[str] = None
    patient_column: Optional[str] = None
    batch_column: Optional[str] = None
    gene_symbol_source: Optional[str] = None
    has_umap: bool = False
    has_pca: bool = False
    has_raw: bool = False
    n_cells: int = 0
    n_genes: int = 0


class MetadataDetector:
    """
    Detects the most likely metadata fields in an AnnData object
    and standardizes them into common OmniSeqAI columns:

    - obs['cell_type']
    - obs['condition']
    - obs['sample_id']
    - obs['patient_id']
    - obs['batch_id']
    """

    CELLTYPE_CANDIDATES = [
        "cell_type",
        "celltype",
        "celltype.final",
        "celltype_broad",
        "celltype.broad_clustering_annot",
        "cell_ontology_class",
        "annotation",
        "annot",
        "cell_annotation",
        "broad_cell_type",
        "major_cell_type",
    ]

    CONDITION_CANDIDATES = [
        "condition",
        "disease",
        "status",
        "group",
        "treatment",
        "phenotype",
        "diagnosis",
        "sample_group",
        "cv19_vax_boost_or_hc_status",
        "cv19_vax_boost_or_hc_status".lower(),
    ]

    SAMPLE_CANDIDATES = [
        "sample",
        "sample_id",
        "sample_name",
        "sampleName_with_day",
        "library_prep_batch",
        "donor_sample",
        "donor_sample_id",
    ]

    PATIENT_CANDIDATES = [
        "patient",
        "patient_id",
        "donor",
        "donor_id",
        "subject",
        "subject_id",
        "individual",
        "person_id",
    ]

    BATCH_CANDIDATES = [
        "batch",
        "batch_id",
        "library_prep_batch",
        "run",
        "sequencing_batch",
        "prep_batch",
    ]

    CONDITION_VALUE_HINTS = [
        "healthy",
        "normal",
        "covid",
        "covid-19",
        "covid19",
        "infected",
        "infection",
        "disease",
        "control",
        "case",
        "treated",
        "treatment",
        "vaccin",
        "vax",
        "booster",
        "hc",
    ]

    CONDITION_NORMALIZATION_MAP = {
        "normal": "Healthy",
        "healthy": "Healthy",
        "hc": "Healthy",
        "control": "Healthy",
        "covid": "COVID",
        "covid-19": "COVID",
        "covid19": "COVID",
        "cv_19": "COVID",
        "infected": "COVID",
        "infection": "COVID",
        "vax": "Vaccinated",
        "vaccinated": "Vaccinated",
        "booster": "Vaccinated",
        "treated": "Treated",
        "treatment": "Treated",
        "disease": "Disease",
    }

    def detect(self, adata) -> Dict[str, Any]:
        profile = MetadataProfile(
            cell_type_column=self._detect_by_name(adata.obs.columns, self.CELLTYPE_CANDIDATES),
            condition_column=self._detect_condition_column(adata),
            sample_column=self._detect_by_name(adata.obs.columns, self.SAMPLE_CANDIDATES),
            patient_column=self._detect_by_name(adata.obs.columns, self.PATIENT_CANDIDATES),
            batch_column=self._detect_by_name(adata.obs.columns, self.BATCH_CANDIDATES),
            gene_symbol_source=self._detect_gene_symbol_source(adata),
            has_umap="X_umap" in adata.obsm,
            has_pca="X_pca" in adata.obsm,
            has_raw=adata.raw is not None,
            n_cells=int(adata.n_obs),
            n_genes=int(adata.n_vars),
        )

        return profile.__dict__

    def apply(self, adata, profile: Optional[Dict[str, Any]] = None):
        """
        Standardize metadata columns in-place on a copy of adata:
        - gene symbols
        - cell_type
        - condition
        - sample_id
        - patient_id
        - batch_id
        """
        adata = adata.copy()

        if profile is None:
            profile = self.detect(adata)

        # ---- gene axis ----
        adata = GeneMapper().fix_gene_names(adata)

        # keep a record of the original gene identifiers
        if "gene_id" not in adata.var.columns:
            adata.var["gene_id"] = adata.var_names.astype(str)

        # ---- cell type ----
        ct_col = profile.get("cell_type_column")
        if ct_col and "cell_type" not in adata.obs.columns:
            adata.obs["cell_type"] = adata.obs[ct_col].astype(str)

        # ---- condition ----
        cond_col = profile.get("condition_column")
        if cond_col:
            if "condition_raw" not in adata.obs.columns:
                adata.obs["condition_raw"] = adata.obs[cond_col].astype(str)

            normalized = self._normalize_condition_series(adata.obs[cond_col])
            adata.obs["condition"] = normalized.astype(str)

        # ---- sample ----
        sample_col = profile.get("sample_column")
        if sample_col and "sample_id" not in adata.obs.columns:
            adata.obs["sample_id"] = adata.obs[sample_col].astype(str)

        # ---- patient ----
        patient_col = profile.get("patient_column")
        if patient_col and "patient_id" not in adata.obs.columns:
            adata.obs["patient_id"] = adata.obs[patient_col].astype(str)

        # ---- batch ----
        batch_col = profile.get("batch_column")
        if batch_col and "batch_id" not in adata.obs.columns:
            adata.obs["batch_id"] = adata.obs[batch_col].astype(str)

        return adata

    def _detect_by_name(self, columns, candidates):
        lower_to_original = {str(c).lower(): c for c in columns}

        for candidate in candidates:
            c = candidate.lower()
            if c in lower_to_original:
                return lower_to_original[c]

        for col in columns:
            col_lower = str(col).lower()
            for candidate in candidates:
                if candidate.lower() in col_lower:
                    return col

        return None

    def _is_textual(self, series: pd.Series) -> bool:
        return (
            pd.api.types.is_object_dtype(series)
            or pd.api.types.is_string_dtype(series)
            or isinstance(series.dtype, pd.CategoricalDtype)
        )

    def _detect_condition_column(self, adata):
        # First try name-based detection.
        by_name = self._detect_by_name(adata.obs.columns, self.CONDITION_CANDIDATES)
        if by_name is not None:
            return by_name

        # Then fall back to value-based detection on text-like columns.
        best_col = None
        best_score = 0

        for col in adata.obs.columns:
            series = adata.obs[col]

            if not self._is_textual(series):
                continue

            values = series.astype(str).str.lower().dropna().unique().tolist()

            # Skip high-cardinality columns that are unlikely to be condition labels.
            if len(values) < 2 or len(values) > 12:
                continue

            score = 0
            for v in values:
                if any(hint in v for hint in self.CONDITION_VALUE_HINTS):
                    score += 1

            if score > best_score and score >= 2:
                best_score = score
                best_col = col

        return best_col

    def _detect_gene_symbol_source(self, adata):
        if "feature_name" in adata.var.columns:
            return "feature_name"
        if "gene_symbol" in adata.var.columns:
            return "gene_symbol"
        if len(adata.var_names) > 0 and not str(adata.var_names[0]).startswith("ENSG"):
            return "var_names"
        return None

    def _normalize_condition_series(self, series: pd.Series) -> pd.Series:
        def normalize_one(v):
            s = str(v).strip().lower()

            # direct map first
            for key, mapped in self.CONDITION_NORMALIZATION_MAP.items():
                if key in s:
                    return mapped

            # keep original text if unknown
            return str(v)

        return series.astype(str).map(normalize_one)
