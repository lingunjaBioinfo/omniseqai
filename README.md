# OmniSeqAI

**OmniSeqAI** is a single-cell RNA-seq analysis pipeline that turns raw or processed scRNA-seq datasets into interpretable biological reports.

It supports automated mode selection, condition-based differential expression, exploratory clustering, pseudobulk analysis, report figures, machine-readable result tables, PDF reporting, and biological signature validation.

The current validated demo uses the Kang IFN-beta PBMC dataset and correctly detects a strong interferon-stimulated antiviral response in IFN-beta-treated immune cells.

---
## Installation

Recommended setup uses Conda:

```bash
conda env create -f environment.yml
conda activate omniseqai
```

To update an existing environment:

```bash
conda env update -f environment.yml --prune
conda activate omniseqai
```

Pip-only fallback:

```bash
pip install -r requirements.txt
```

For PDF text validation in the regression script, `pdftotext` is required. The Conda environment installs it through `poppler`.

On Ubuntu, it can also be installed manually:

```bash
sudo apt install poppler-utils
```

## Core Capabilities

OmniSeqAI currently supports:

* Flexible input loading for `.h5ad`, 10x `.h5`, `.loom`, `.csv`, `.tsv`, `.txt`, 10x MTX folders, and multi-sample directories.
* Automatic or user-selected analysis mode.
* Condition-based differential expression.
* Cell-type-specific pseudobulk differential expression.
* Exploratory clustering and marker-gene discovery.
* UMAP visualization by condition, cell type, sample, and batch.
* Volcano plots with biologically meaningful gene labels.
* Pseudobulk heatmaps.
* Cell-type proportion plots.
* Machine-readable CSV/JSON output tables.
* PDF and text reports.
* Biology validation using built-in curated gene-signature panels.
* Cell-type localization of detected biological programs.
* Regression testing using a validated Kang IFN-beta run.

---

## Analysis Modes

OmniSeqAI supports three main modes:

| Mode          | Purpose                                                                                    |
| ------------- | ------------------------------------------------------------------------------------------ |
| `auto`        | Automatically choose the most appropriate analysis mode.                                   |
| `condition`   | Compare biological conditions such as Healthy vs Disease or Control vs Treatment.          |
| `exploratory` | Run clustering, UMAP, and marker-gene discovery when no condition comparison is available. |

Example:

```bash
python run_omniseqai.py \
  --input data/kang_ifnb/kang_ifnb.h5ad \
  --mode condition \
  --output report.txt \
  --pdf report.pdf
```

For organized run folders, use:

```bash
./scripts/run_omniseqai_runfolder.sh \
  --input data/kang_ifnb/kang_ifnb.h5ad \
  --mode condition \
  --run-name kang_biology_validation
```

---

## Validated Kang IFN-beta Demo

The Kang IFN-beta PBMC dataset is used as the main biological validation test.

Run:

```bash
./scripts/run_omniseqai_runfolder.sh \
  --input data/kang_ifnb/kang_ifnb.h5ad \
  --mode condition \
  --run-name kang_biology_validation
```

Outputs are written to:

```text
runs/kang_biology_validation/
```

Expected core outputs:

```text
report.txt
report.pdf
figures/
tables/
```

Expected figures:

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

Expected biology result:

```text
Dominant program: interferon_antiviral_response
Evidence strength: 28/28 expected genes detected in Healthy vs IFN_beta
Main driver cell types: CD4 T cells, CD14+ Monocytes, B cells, FCGR3A+ Monocytes, Dendritic cells, NK cells
Representative supporting genes: HERC6, DDX58, IRF7, STAT1, OAS2, GBP1, IFITM2, ISG15, UBE2L6, HERC5, IFIT2, MX1
```

This confirms that OmniSeqAI detects the expected interferon-stimulated antiviral response after IFN-beta stimulation.

---

## Biology Validation Layer

OmniSeqAI includes a biology validation engine that checks whether differential expression results match known biological programs.

Current supported signatures include:

* Interferon antiviral response
* Inflammatory myeloid activation
* Antigen presentation and MHC activity
* Cytotoxic T/NK activation
* T cell exhaustion/checkpoint state
* B cell/plasma activation
* Cell-cycle proliferation
* Complement macrophage/monocyte program

The biology report layer produces:

1. A dominant biological conclusion.
2. Driver cell types for the detected program.
3. Representative supporting genes.
4. Additional detected biological programs.
5. Signature evidence tables.
6. Cell-type signature localization heatmaps.

---

## Regression Check

Before committing future changes, run:

```bash
./scripts/check_omniseqai_kang.sh
```

The regression check verifies:

* Python compilation
* Kang IFN-beta condition-mode execution
* text report generation
* PDF report generation
* expected figures
* expected result tables
* biological conclusion section
* interferon signature detection
* PDF biology content

Expected final output:

```text
==========================================
 OmniSeqAI regression check PASSED
 Run directory: runs/kang_biology_validation_regression
==========================================
```

---

## Output Tables

OmniSeqAI exports machine-readable results under `tables/`.

Common outputs include:

```text
run_summary.json
celltype_counts.csv
celltype_counts_by_condition.csv
celltype_proportions_by_condition.csv
biology_validation/signature_summary.csv
biology_validation/signature_hits.csv
condition_de/
celltype_specific_de/
integration_qc/
```

These tables are intended for downstream review, external plotting, manuscript preparation, and validation.

---

## Project Structure

```text
backend/
  analysis_router.py
  biology_report_utils.py
  biology_signatures.py
  biology_validator.py
  gene_symbol_utils.py
  input_loader.py
  integration_qc.py
  integration_report_utils.py
  pipeline_orchestrator.py
  pseudobulk_de.py
  report_figures.py
  report_table_utils.py
  router_pdf_report.py
  router_report.py
  table_exporter.py

scripts/
  run_omniseqai_runfolder.sh
  check_omniseqai_kang.sh

run_omniseqai.py
```

---

## Current Status

OmniSeqAI currently performs validated end-to-end analysis on the Kang IFN-beta dataset.

The pipeline can produce:

* UMAP figures
* volcano plots
* pseudobulk heatmaps
* biology signature evidence plots
* cell-type biology localization heatmaps
* text reports
* PDF reports
* result tables
* regression-tested analysis outputs

---

## Development Rule

Before pushing future changes, run:

```bash
./scripts/check_omniseqai_kang.sh
```

Only commit if the regression check passes.

