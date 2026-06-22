# OmniSeqAI Roadmap

## Completed

* Flexible scRNA-seq input loading
* Automatic analysis mode selection
* Condition-based differential expression
* Cell-type-specific pseudobulk differential expression
* Exploratory clustering and marker detection
* UMAP report figures
* Volcano plots
* Pseudobulk heatmaps
* Cell-type proportion plots
* Text reports
* PDF reports
* Machine-readable table export
* Integration and batch QC summaries
* Built-in curated gene-signature biology validation
* Cell-type biology localization heatmap
* Biological conclusion report layer
* Kang IFN-beta regression check
* Reproducible Conda and pip installation files

## Near-Term Priorities

### 1. User-defined biology signatures

Allow users to provide their own gene-signature files.

Example:

```text
signature_name,gene
interferon_response,ISG15
interferon_response,IFIT1
```

### 2. Disease-specific signature panels

Add curated panels for:

* cancer immune microenvironment
* viral infection
* autoimmune inflammation
* neuroinflammation
* fibrosis
* cell-cycle/proliferation
* exhaustion/checkpoint states

### 3. Stronger evidence grading

Separate detected biology into:

* strong evidence
* moderate evidence
* weak supportive signal

### 4. Better multi-batch integration support

Improve Harmony and batch-QC reporting for true multi-batch datasets.

### 5. Better demo packaging

Add small public demo instructions and example output screenshots.

### 6. CLI improvements

Add commands such as:

```bash
omniseqai run
omniseqai validate
omniseqai demo
```

### 7. Automated tests

Add unit and integration tests around:

* input loading
* analysis routing
* table export
* biology validation
* report generation

## Long-Term Goals

* web dashboard
* interactive report viewer
* plugin-style disease signature modules
* pathway enrichment integration
* multi-omics extension
* cloud execution support

