# User-defined biology signatures

OmniSeqAI supports built-in curated biology signatures and optional user-defined gene-signature files.

User-defined signatures allow researchers to test their own biological programs, pathways, marker panels, or disease-specific gene sets against OmniSeqAI differential-expression results.

## CLI usage

```bash
python run_omniseqai.py \
  --input data/kang_ifnb/kang_ifnb.h5ad \
  --mode condition \
  --output report.txt \
  --pdf report.pdf \
  --signatures examples/custom_ifn_signature.csv
```

The run-folder wrapper also supports the same option:

```bash
./scripts/run_omniseqai_runfolder.sh \
  --input data/kang_ifnb/kang_ifnb.h5ad \
  --mode condition \
  --run-name kang_custom_signature_test \
  --signatures examples/custom_ifn_signature.csv
```

## Required columns

A signature file must contain at least:

```csv
signature_name,gene
custom_ifn_response,ISG15
custom_ifn_response,IFIT1
custom_ifn_response,MX1
```

Accepted signature-name columns:

```text
signature_name
signature
pathway
program
gene_set
geneset
```

Accepted gene columns:

```text
gene
gene_symbol
symbol
gene_name
```

## Optional columns

```text
description
expected_direction
```

Example:

```csv
signature_name,gene,description,expected_direction
custom_ifn_response,ISG15,User-defined interferon response,case_up
custom_ifn_response,IFIT1,User-defined interferon response,case_up
custom_ifn_response,MX1,User-defined interferon response,case_up
```

## Expected direction

The `expected_direction` column tells OmniSeqAI how to interpret the signature.

Accepted values include:

```text
case_up
up
upregulated
case_down
down
downregulated
either
changed
```

For a condition comparison such as:

```text
Healthy vs IFN_beta
```

`case_up` means the signature is expected to be increased in the second condition, here `IFN_beta`.

## Output files

User-defined signatures are included in the same biology-validation outputs as built-in signatures:

```text
tables/biology_validation/signature_summary.csv
tables/biology_validation/signature_hits.csv
report.txt
report.pdf
figures/biology_signature_hits.png
figures/biology_celltype_signature_heatmap.png
```

In `signature_summary.csv`, user signatures are marked with:

```text
source=user
```

## Example

The repository includes:

```text
examples/custom_ifn_signature.csv
```

This example tests a small interferon-response panel against the Kang IFN-beta dataset.
