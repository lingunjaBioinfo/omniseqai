#!/usr/bin/env bash

set -euo pipefail

INPUT="data/kang_ifnb/kang_ifnb.h5ad"
RUN_NAME="kang_biology_validation_regression"
RUN_DIR="runs/${RUN_NAME}"

echo "=========================================="
echo " OmniSeqAI Kang IFN-beta regression check "
echo "=========================================="
echo ""

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input file not found: ${INPUT}"
  exit 1
fi

echo "[1/7] Compiling Python files..."
find backend -name "*.py" -print0 | xargs -0 -n1 python -m py_compile
python -m py_compile run_omniseqai.py
echo "Python compile check passed."
echo ""

echo "[2/7] Removing old regression run..."
rm -rf "${RUN_DIR}"
echo "Removed: ${RUN_DIR}"
echo ""

echo "[3/7] Running OmniSeqAI condition mode..."
./scripts/run_omniseqai_runfolder.sh \
  --input "${INPUT}" \
  --mode condition \
  --run-name "${RUN_NAME}"
echo ""

assert_file() {
  local path="$1"

  if [[ ! -f "${path}" ]]; then
    echo "ERROR: expected file missing: ${path}"
    exit 1
  fi

  echo "OK: ${path}"
}

assert_dir() {
  local path="$1"

  if [[ ! -d "${path}" ]]; then
    echo "ERROR: expected directory missing: ${path}"
    exit 1
  fi

  echo "OK: ${path}"
}

echo "[4/7] Checking main reports..."
assert_file "${RUN_DIR}/report.txt"
assert_file "${RUN_DIR}/report.pdf"
assert_dir "${RUN_DIR}/figures"
assert_dir "${RUN_DIR}/tables"
echo ""

echo "[5/7] Checking expected figures..."
assert_file "${RUN_DIR}/figures/umap_condition_standard.png"
assert_file "${RUN_DIR}/figures/umap_celltype_annotated.png"
assert_file "${RUN_DIR}/figures/celltype_proportions_by_condition.png"
assert_file "${RUN_DIR}/figures/volcano_Healthy_vs_IFN_beta_standard.png"
assert_file "${RUN_DIR}/figures/pseudobulk_heatmap_Healthy_vs_IFN_beta.png"
assert_file "${RUN_DIR}/figures/umap_sample.png"
assert_file "${RUN_DIR}/figures/umap_batch.png"
assert_file "${RUN_DIR}/figures/biology_signature_hits.png"
assert_file "${RUN_DIR}/figures/biology_celltype_signature_heatmap.png"
echo ""

echo "[6/7] Checking expected tables..."
assert_file "${RUN_DIR}/tables/run_summary.json"
assert_file "${RUN_DIR}/tables/biology_validation/signature_summary.csv"
assert_file "${RUN_DIR}/tables/biology_validation/signature_hits.csv"
assert_file "${RUN_DIR}/tables/celltype_counts.csv"
assert_file "${RUN_DIR}/tables/celltype_counts_by_condition.csv"
assert_file "${RUN_DIR}/tables/celltype_proportions_by_condition.csv"
echo ""

echo "[7/7] Checking biological report content..."

if ! grep -q "Biological conclusion" "${RUN_DIR}/report.txt"; then
  echo "ERROR: Biological conclusion section missing from report.txt"
  exit 1
fi

echo "OK: Biological conclusion section found."

if ! grep -q "interferon_antiviral_response" "${RUN_DIR}/report.txt"; then
  echo "ERROR: interferon_antiviral_response missing from report.txt"
  exit 1
fi

echo "OK: interferon_antiviral_response found."

if grep -q "Additional interpretation" "${RUN_DIR}/report.txt"; then
  echo "ERROR: redundant Additional interpretation section still present."
  exit 1
fi

echo "OK: redundant Additional interpretation section absent."

if command -v pdftotext >/dev/null 2>&1; then
  if ! pdftotext "${RUN_DIR}/report.pdf" - | grep -q "Biological conclusion"; then
    echo "ERROR: Biological conclusion missing from PDF."
    exit 1
  fi

  echo "OK: Biological conclusion found in PDF."
else
  echo "WARNING: pdftotext not installed; skipping PDF text validation."
fi

echo ""
echo "=========================================="
echo " OmniSeqAI regression check PASSED "
echo " Run directory: ${RUN_DIR}"
echo "=========================================="
