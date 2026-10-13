#!/bin/sh
# Run the r1 tests; outputs go to outputs/. Exit status is nonzero if any check fails.
cd "$(dirname "$0")"
mkdir -p outputs
st=0
# Failure detection (2026-10-12).  Primary: the exit status of each command (every test script exits 1 when one of
# its checks fails or a falsifier does not fail).  Secondary: a scan of the output for failure lines in any of the
# formats used in tests_cm: [FAIL], RESULT|..|FAIL, a bare FAIL at the start of a line, SOME FAILED / SOME FAIL,
# DID NOT FAIL, EXIT STATUS 1, "<n> FAIL" / "<n> failures" with n >= 1.  PASS lines and "fails as required" are ignored,
# so "0 FAIL" and "0 failures" do not count.
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
nfail_lines() { grep -E "$FAIL_PAT" "$1" | grep -v -E "$OK_PAT" | wc -l | tr -d ' '; }

for t in test_classgroup test_extraction_rows test_step45_bookkeeping; do
  python3 -I $t.py > outputs/$t.txt 2>&1; rc=$?
  nf=$(nfail_lines outputs/$t.txt)
  if [ $rc -ne 0 ] || [ $nf -ne 0 ]; then st=1; echo "[FAIL] $t: exit $rc, failure lines $nf"; fi
  tail -1 outputs/$t.txt
done
exit $st
