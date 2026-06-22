# OmniSeqAI Validation

OmniSeqAI is currently validated using the Kang IFN-beta PBMC single-cell RNA-seq dataset.

The validation run checks whether the pipeline correctly detects the expected interferon-stimulated antiviral response in IFN-beta-treated immune cells.

## Regression Command

```bash
./scripts/check_omniseqai_kang.sh
```

## Expected Result

The regression check should finish with:

```text
OmniSeqAI regression check PASSED
```

## Expected Biological Finding

The Kang IFN-beta validation run should detect:

```text
Dominant program: interferon_antiviral_response
Evidence strength: 28/28 expected genes detected in Healthy vs IFN_beta
```

## Expected Driver Cell Types

The dominant interferon response should localize mainly to immune cell populations including:

* CD4 T cells
* CD14+ Monocytes
* B cells
* FCGR3A+ Monocytes
* Dendritic cells
* NK cells

## Expected Figures

The validation run should generate:

```text
biology_celltype_signature_heatmap.png
biology_signature_hits.png
celltype_proportions_by_condition.png
pseudobulk_heatmap_Healthy_vs_IFN_beta.png
umap_batch.png
umap_celltype_annotated.png
umap_condition_standard.png
umap_sample.png
volcano_Healthy_vs_IFN_beta_standard.png
```

## Expected Tables

The validation run should generate:

```text
tables/run_summary.json
tables/biology_validation/signature_summary.csv
tables/biology_validation/signature_hits.csv
tables/celltype_counts.csv
tables/celltype_counts_by_condition.csv
tables/celltype_proportions_by_condition.csv
```

## Validation Purpose

This validation confirms that OmniSeqAI can:

1. Load a real scRNA-seq dataset.
2. Detect metadata columns.
3. Run condition-based differential expression.
4. Generate report figures.
5. Export machine-readable result tables.
6. Generate text and PDF reports.
7. Detect a known biological response using curated gene-signature panels.
8. Localize the detected biological program to specific cell types.

