#!/usr/bin/env bash
# Verify the published snapshot without R or the large input data.
# Run from the repository root: bash tools/verify_release.sh
set -u

cd "$(dirname "$0")/.." || exit 1
failures=0

# Check the manifest entries whose files are present; report absent ones.
check_present() {
  local label="$1" manifest="$2" exclude="${3:-}"
  local present absent ok
  present=$(awk -v ex="$exclude" '$2 != ex' "$manifest" | while read -r hash path; do
    [ -f "$path" ] && printf '%s  %s\n' "$hash" "$path"
  done)
  absent=$(awk -v ex="$exclude" '$2 != ex' "$manifest" | while read -r hash path; do
    [ -f "$path" ] || printf '%s\n' "$path"
  done | grep -c . || true)
  if [ -z "$present" ]; then
    echo "FAIL  $label: no listed files are present"
    failures=$((failures + 1))
    return
  fi
  if ok=$(printf '%s\n' "$present" | sha256sum -c --quiet - 2>&1); then
    echo "PASS  $label: $(printf '%s\n' "$present" | grep -c .) present files match; $absent not distributed"
  else
    echo "FAIL  $label:"
    printf '%s\n' "$ok" | sed 's/^/        /'
    failures=$((failures + 1))
  fi
}

echo "== Provenance =="
check_present "validated-run source code (source_checksums.sha256)" \
  analysis_v2/provenance/source_checksums.sha256
check_present "validated-run committed outputs (output_manifest.sha256)" \
  analysis_v2/provenance/output_manifest.sha256
check_present "frozen inputs (input_checksums.sha256)" \
  analysis_v2/provenance/input_checksums.sha256
check_present "legacy 2025 baseline (baseline_2025-09-22.sha256)" \
  provenance/baseline_2025-09-22.sha256 "DepMapCrispr/DepMapCrispr.R"

echo
echo "== Sanitization =="
leaks=$(grep -rIlE '(/home/[A-Za-z0-9._-]+/|/Users/[A-Za-z0-9._-]+/)' \
  --exclude-dir=.git --exclude=verify_release.sh . || true)
if [ -n "$leaks" ]; then
  echo "FAIL  absolute home-directory paths found in:"
  printf '%s\n' "$leaks" | sed 's/^/        /'
  failures=$((failures + 1))
else
  echo "PASS  no absolute home-directory paths in text files"
fi

echo
echo "== Inventory (static counts) =="
tests_dir=analysis_v2/tests/testthat
echo "Test files:          $(find "$tests_dir" -name 'test-*.R' | wc -l)"
echo "test_that blocks:    $(grep -ro 'test_that(' "$tests_dir" | wc -l)"
echo "expect_* calls:      $(grep -roE 'expect_[a-z_]+\(' "$tests_dir" | wc -l)"
echo "Recorded warnings:   $(grep -c '^W[0-9]' analysis_v2/provenance/warnings.csv)"
echo "Audit findings:      $(grep -cE '^\| L[0-9]+ ' docs/AUDIT_REMEDIATION.md) remediated, $(grep -cE '^\| O[0-9]+ ' docs/AUDIT_REMEDIATION.md) open"

echo
if [ "$failures" -eq 0 ]; then
  echo "All checks passed."
else
  echo "$failures check(s) failed."
fi
exit "$failures"
