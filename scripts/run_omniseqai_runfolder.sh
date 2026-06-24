#!/usr/bin/env bash
set -euo pipefail

INPUT=""
MODE="auto"
RUN_NAME=""
SIGNATURES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT="$2"
      shift 2
      ;;
    --signatures)
      SIGNATURES="$2"
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    --run-name)
      RUN_NAME="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT" ]]; then
  echo "ERROR: --input is required"
  exit 1
fi

if [[ ! -f "$INPUT" && ! -d "$INPUT" ]]; then
  echo "ERROR: input does not exist: $INPUT"
  exit 1
fi

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
INPUT_BASE=$(basename "$INPUT")
INPUT_BASE="${INPUT_BASE%.*}"

if [[ -z "$RUN_NAME" ]]; then
  RUN_NAME="${TIMESTAMP}_${INPUT_BASE}_${MODE}"
fi

RUN_DIR="runs/${RUN_NAME}"

mkdir -p "$RUN_DIR"
mkdir -p "$RUN_DIR/figures"
mkdir -p "$RUN_DIR/logs"
mkdir -p "$RUN_DIR/metadata"
mkdir -p "$RUN_DIR/tables"

rm -rf "$RUN_DIR/figures"/*
rm -rf "$RUN_DIR/tables"/*

REPORT_TXT="$RUN_DIR/report.txt"
REPORT_PDF="$RUN_DIR/report.pdf"
LOG_FILE="$RUN_DIR/logs/run.log"

echo "Run directory: $RUN_DIR"
echo "Input: $INPUT"
echo "Mode: $MODE"
echo "Text report: $REPORT_TXT"
echo "PDF report: $REPORT_PDF"
echo ""

# ------------------------------------------------------------
# Clear global figure directory before this run.
# This prevents stale figures from older runs being copied into
# the current run folder.
# ------------------------------------------------------------
mkdir -p outputs/report_figures
mkdir -p outputs/tables

rm -f outputs/report_figures/*.png
rm -rf outputs/tables/*

# ------------------------------------------------------------
# Run OmniSeqAI.
# Capture both stdout and stderr so warnings are stored in run.log.
# ------------------------------------------------------------
SIGNATURE_ARGS=()

if [[ -n "${SIGNATURES}" ]]; then
  SIGNATURE_ARGS+=(--signatures "${SIGNATURES}")
fi
python run_omniseqai.py \
  --input "$INPUT" \
  --mode "$MODE" \
  --output "$REPORT_TXT" \
  --pdf "$REPORT_PDF" \
  "${SIGNATURE_ARGS[@]}" 2>&1 | tee "$LOG_FILE"

# ------------------------------------------------------------
# Copy only figures generated during this run.
# Since outputs/report_figures was cleared before the run,
# these should belong only to the current run.
# ------------------------------------------------------------
if [[ -d "outputs/report_figures" ]]; then
  find outputs/report_figures -maxdepth 1 -type f -name "*.png" -exec cp -f {} "$RUN_DIR/figures/" \;
fi
if [[ -d "outputs/tables" ]]; then
  rm -rf "$RUN_DIR/tables"/*
  cp -a outputs/tables/. "$RUN_DIR/tables/"
fi

# ------------------------------------------------------------
# Write run manifest.
# ------------------------------------------------------------
cat > "$RUN_DIR/metadata/run_manifest.json" <<EOF
{
  "run_name": "$RUN_NAME",
  "run_dir": "$RUN_DIR",
  "input": "$INPUT",
  "mode": "$MODE",
  "report_txt": "$REPORT_TXT",
  "report_pdf": "$REPORT_PDF",
  "log_file": "$LOG_FILE",
  "timestamp": "$TIMESTAMP"
}
EOF

echo ""
echo "Run complete."
echo "Outputs saved in: $RUN_DIR"
