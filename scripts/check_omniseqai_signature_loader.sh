#!/usr/bin/env bash
set -euo pipefail

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

VALID_ONE_GENE_PER_ROW="${TMP_DIR}/valid_one_gene_per_row.csv"
VALID_GENE_LIST="${TMP_DIR}/valid_gene_list.csv"
INVALID_MISSING_GENE="${TMP_DIR}/invalid_missing_gene.csv"
INVALID_DIRECTION="${TMP_DIR}/invalid_direction.csv"

cat > "${VALID_ONE_GENE_PER_ROW}" <<'CSV'
signature_name,gene,description,expected_direction
test_ifn,ISG15,Test interferon panel,case_up
test_ifn,IFIT1,Test interferon panel,case_up
test_ifn,MX1,Test interferon panel,case_up
CSV

cat > "${VALID_GENE_LIST}" <<'CSV'
signature_name,genes,description,expected_direction
test_multi,"ISG15, IFIT1; MX1|OAS1",Test multi-gene panel,either
CSV

cat > "${INVALID_MISSING_GENE}" <<'CSV'
signature_name,description
bad_signature,Missing gene column
CSV

cat > "${INVALID_DIRECTION}" <<'CSV'
signature_name,gene,description,expected_direction
bad_direction,ISG15,Bad direction,wrong_direction
CSV

python - <<PY
from backend.signature_loader import load_user_signatures

one = load_user_signatures("${VALID_ONE_GENE_PER_ROW}")
assert "test_ifn" in one
assert one["test_ifn"]["genes"] == ["ISG15", "IFIT1", "MX1"]
assert one["test_ifn"]["source"] == "user"

multi = load_user_signatures("${VALID_GENE_LIST}")
assert "test_multi" in multi
assert multi["test_multi"]["genes"] == ["ISG15", "IFIT1", "MX1", "OAS1"]
assert multi["test_multi"]["expected_direction"] == "either"

try:
    load_user_signatures("${INVALID_MISSING_GENE}")
except ValueError as exc:
    assert "gene" in str(exc).lower()
else:
    raise AssertionError("Missing gene column did not raise ValueError")

try:
    load_user_signatures("${INVALID_DIRECTION}")
except ValueError as exc:
    assert "invalid expected_direction" in str(exc).lower()
else:
    raise AssertionError("Invalid expected_direction did not raise ValueError")

print("Signature-loader regression: PASSED")
PY
