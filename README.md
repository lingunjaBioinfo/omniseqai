# OmniSeqAI

AI-powered single-cell RNA-seq analysis platform.

## Features

### Data Processing
- Quality control
- Cell filtering
- Normalization
- Highly variable gene selection
- PCA
- UMAP
- Leiden clustering

### Cell Annotation
- Automated cell type annotation using CellTypist
- Cluster-level cell type assignment

### Marker Detection
- Marker gene identification
- Cluster summaries

### Biological Interpretation
- Automated interpretation of cluster biology
- Marker-based explanations

### Reporting
- Text report generation
- PDF report generation
- Embedded visualizations

### Visualization
- QC plots
- UMAP clusters
- UMAP cell-type maps

### Differential Expression
- Cluster-vs-rest differential expression
- Top marker identification

---

## Project Structure

```text
omniseqai/
│
├── backend/
│   ├── pipeline.py
│   ├── annotation.py
│   ├── markers.py
│   ├── interpretation.py
│   ├── differential_expression.py
│   ├── report.py
│   └── visualization.py
│
├── frontend/
│   └── app.py
│
├── data/
├── outputs/
├── reports/
│
├── test_pipeline.py
├── test_annotation.py
├── test_markers.py
├── test_interpretation.py
├── test_de.py
├── test_report.py
│
└── README.md
