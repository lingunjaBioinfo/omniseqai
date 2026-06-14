from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Set

import numpy as np
import pandas as pd


SIGNATURE_GENE_SETS: Dict[str, Set[str]] = {
    "Interferon-stimulated antiviral response": {
        "ISG15",
        "ISG20",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "IFI6",
        "IFI27",
        "IFI35",
        "IFI44",
        "IFI44L",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "OAS3",
        "OASL",
        "STAT1",
        "STAT2",
        "IRF7",
        "IRF9",
        "DDX58",
        "DDX60",
        "IFIH1",
        "HERC5",
        "HERC6",
        "UBE2L6",
        "BST2",
        "RSAD2",
        "EPSTI1",
        "TRIM14",
        "BATF2",
        "CXCL10",
        "ZBP1",
        "GBP1",
        "GBP2",
        "GBP4",
        "GBP5",
        "PARP14",
        "SP110",
        "WARS",
    },
    "Inflammatory cytokine and chemokine signaling": {
        "IL1B",
        "IL1A",
        "IL6",
        "TNF",
        "CXCL8",
        "IL8",
        "CXCL10",
        "CXCL11",
        "CCL2",
        "CCL3",
        "CCL4",
        "CCL7",
        "CCL8",
        "NFKBIA",
        "NFKBIZ",
        "JUN",
        "FOS",
        "FOSB",
        "JUNB",
        "IER3",
    },
    "Antigen presentation and MHC response": {
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-E",
        "B2M",
        "TAP1",
        "TAP2",
        "PSMB8",
        "PSMB9",
        "PSME1",
        "PSME2",
        "CIITA",
    },
    "Cytotoxic lymphocyte activity": {
        "NKG7",
        "GNLY",
        "GZMB",
        "GZMA",
        "GZMK",
        "PRF1",
        "KLRD1",
        "KLRB1",
        "KLRF1",
        "KLRG1",
        "FGFBP2",
        "CTSW",
        "CD8A",
        "CD8B",
    },
    "Monocyte activation": {
        "LYZ",
        "S100A8",
        "S100A9",
        "S100A12",
        "FCN1",
        "VCAN",
        "LST1",
        "CTSS",
        "FCGR3A",
        "MS4A7",
        "CD14",
        "CD68",
        "LGALS3",
        "AIF1",
    },
}


DEFAULT_PRIORITY_GENES: List[str] = [
    "ISG15",
    "ISG20",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFI6",
    "IFI44L",
    "MX1",
    "OAS1",
    "OAS2",
    "STAT1",
    "IRF7",
    "CXCL10",
    "DDX58",
    "HERC6",
    "UBE2L6",
    "BST2",
    "RSAD2",
    "EPSTI1",
    "GBP1",
    "GBP5",
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL8",
    "IL1B",
    "TNF",
    "NFKBIA",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "TAP1",
    "TAP2",
]


def normalize_gene_name(gene: object) -> str:
    return str(gene).strip().upper()


def _prepare_de_table(
    de_results: pd.DataFrame,
    gene_col: str = "names",
    lfc_col: str = "logfoldchanges",
    padj_col: str = "pvals_adj",
) -> pd.DataFrame:
    if de_results is None or de_results.empty:
        return pd.DataFrame()

    df = de_results.copy()

    required = {gene_col, lfc_col, padj_col}
    missing = required - set(df.columns)

    if missing:
        return pd.DataFrame()

    df[gene_col] = df[gene_col].astype(str)
    df["gene_upper"] = df[gene_col].map(normalize_gene_name)
    df[lfc_col] = pd.to_numeric(df[lfc_col], errors="coerce")
    df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")

    df = df.dropna(subset=[gene_col, lfc_col, padj_col]).copy()
    df["abs_logfoldchanges"] = df[lfc_col].abs()

    return df


def find_signature_hits(
    de_results: pd.DataFrame,
    gene_col: str = "names",
    lfc_col: str = "logfoldchanges",
    padj_col: str = "pvals_adj",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 0.5,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Return significant up/down signature genes.

    Positive logFC is interpreted as case higher than reference.
    """

    df = _prepare_de_table(
        de_results=de_results,
        gene_col=gene_col,
        lfc_col=lfc_col,
        padj_col=padj_col,
    )

    output: Dict[str, Dict[str, List[str]]] = {}

    if df.empty:
        return output

    sig_up = df[(df[padj_col] < padj_cutoff) & (df[lfc_col] >= lfc_cutoff)].copy()
    sig_down = df[(df[padj_col] < padj_cutoff) & (df[lfc_col] <= -lfc_cutoff)].copy()

    for signature_name, genes in SIGNATURE_GENE_SETS.items():
        genes_upper = {normalize_gene_name(g) for g in genes}

        up_hits = (
            sig_up[sig_up["gene_upper"].isin(genes_upper)]
            .sort_values([padj_col, "abs_logfoldchanges"], ascending=[True, False])
            [gene_col]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )

        down_hits = (
            sig_down[sig_down["gene_upper"].isin(genes_upper)]
            .sort_values([padj_col, "abs_logfoldchanges"], ascending=[True, False])
            [gene_col]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )

        if up_hits or down_hits:
            output[signature_name] = {
                "up": up_hits,
                "down": down_hits,
            }

    return output


def format_signature_hits(
    hits: Dict[str, Dict[str, List[str]]],
    max_genes_per_signature: int = 15,
) -> List[str]:
    lines: List[str] = []

    for signature_name, directions in hits.items():
        up = directions.get("up", [])
        down = directions.get("down", [])

        if up:
            genes = ", ".join(up[:max_genes_per_signature])
            lines.append(f"{signature_name} is increased ({genes}).")

        if down:
            genes = ", ".join(down[:max_genes_per_signature])
            lines.append(f"{signature_name} is decreased ({genes}).")

    return lines


def priority_label_genes(
    de_results: pd.DataFrame,
    gene_col: str = "names",
    lfc_col: str = "logfoldchanges",
    padj_col: str = "pvals_adj",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 0.5,
    max_labels: int = 14,
    extra_priority_genes: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    Pick volcano labels.

    Priority:
    1. Significant known signature genes.
    2. Strongest remaining genes by adjusted p-value and effect size.
    """

    df = _prepare_de_table(
        de_results=de_results,
        gene_col=gene_col,
        lfc_col=lfc_col,
        padj_col=padj_col,
    )

    if df.empty:
        return []

    priority = list(DEFAULT_PRIORITY_GENES)

    if extra_priority_genes:
        priority.extend([str(g) for g in extra_priority_genes])

    priority_upper = [normalize_gene_name(g) for g in priority]

    sig = df[
        (df[padj_col] < padj_cutoff)
        & (df[lfc_col].abs() >= lfc_cutoff)
    ].copy()

    selected: List[str] = []
    selected_upper: Set[str] = set()

    for gene_upper in priority_upper:
        match = sig[sig["gene_upper"] == gene_upper]

        if match.empty:
            continue

        gene = str(match.sort_values(padj_col).iloc[0][gene_col])

        if normalize_gene_name(gene) not in selected_upper:
            selected.append(gene)
            selected_upper.add(normalize_gene_name(gene))

        if len(selected) >= max_labels:
            return selected

    remaining = (
        sig[~sig["gene_upper"].isin(selected_upper)]
        .sort_values([padj_col, "abs_logfoldchanges"], ascending=[True, False])
        [gene_col]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    for gene in remaining:
        gene_upper = normalize_gene_name(gene)

        if gene_upper not in selected_upper:
            selected.append(gene)
            selected_upper.add(gene_upper)

        if len(selected) >= max_labels:
            break

    return selected
