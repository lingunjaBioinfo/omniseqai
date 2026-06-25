#!/usr/bin/env bash
set -euo pipefail

INPUT="data/kang_ifnb/kang_ifnb.h5ad"
SIG_FILE="examples/custom_ifn_signature.csv"
RUN_NAME="kang_custom_signature_regression"
RUN_DIR="runs/${RUN_NAME}"
REPORT_TXT="${RUN_DIR}/report.txt"
SUMMARY_TABLE="${RUN_DIR}/tables/biology_validation/signature_summary.csv"
HITS_TABLE="${RUN_DIR}/tables/biology_validation/signature_hits.csv"

echo "OmniSeqAI custom-signature regression"
echo "Input: ${INPUT}"
echo "Signature file: ${SIG_FILE}"
echo "Run directory: ${RUN_DIR}"

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: Missing input dataset: ${INPUT}" >&2
  exit 1
fi

if [[ ! -f "${SIG_FILE}" ]]; then
  echo "ERROR: Missing signature file: ${SIG_FILE}" >&2
  exit 1
fi

rm -rf "${RUN_DIR}"

./scripts/run_omniseqai_runfolder.sh \
  --input "${INPUT}" \
  --mode condition \
  --run-name "${RUN_NAME}" \
  --signatures "${SIG_FILE}"

if [[ ! -f "${REPORT_TXT}" ]]; then
  echo "ERROR: Missing report: ${REPORT_TXT}" >&2
  exit 1
fi

if [[ ! -f "${SUMMARY_TABLE}" ]]; then
  echo "ERROR: Missing biology summary table: ${SUMMARY_TABLE}" >&2
  exit 1
fi

if [[ ! -f "${HITS_TABLE}" ]]; then
  echo "ERROR: Missing biology hits table: ${HITS_TABLE}" >&2
  exit 1
fi

grep -q "custom_ifn_response" "${REPORT_TXT}"
grep -q "custom_ifn_response" "${SUMMARY_TABLE}"
grep -q "custom_ifn_response" "${HITS_TABLE}"
grep -q "user" "${SUMMARY_TABLE}"
grep -q "True" "${SUMMARY_TABLE}"

echo "Custom-signature regression: PASSED"
