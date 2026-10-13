#!/bin/sh
# Runs the s6 checks (about 10 s). Optional first argument: path to tests_cm/part_B_lines_0550_0997.
cd "$(dirname "$0")" || exit 1
mkdir -p outputs
export PYTHON="${PYTHON:-python3}"
status=0
LIB="${1:-../../part_B_lines_0550_0997}"
# Failure detection (2026-10-12).  Primary: the exit status of each command (every test script exits 1 when one of
# its checks fails or a falsifier does not fail).  Secondary: a scan of the output for failure lines in any of the
# formats used in tests_cm: [FAIL], RESULT|..|FAIL, a bare FAIL at the start of a line, SOME FAILED / SOME FAIL,
# DID NOT FAIL, EXIT STATUS 1, "<n> FAIL" / "<n> failures" with n >= 1.  PASS lines and "fails as required" are ignored,
# so "0 FAIL" and "0 failures" do not count.
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
nfail_lines() { grep -E "$FAIL_PAT" "$1" | grep -v -E "$OK_PAT" | wc -l | tr -d ' '; }
run() {
  name=$1; shift
  "$PYTHON" "$@" > "outputs/$name.txt" 2>&1
  rc=$?
  nf=$(nfail_lines "outputs/$name.txt")
  if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi
  echo "$name: exit $rc, failure lines $nf"
  if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then echo "[FAIL] $name: exit $rc, failure lines $nf"; fi
}
run s6_twisted_paired s6_twisted_paired.py "$LIB"
run s6_twisted_reduction_120_170 s6_twisted_reduction.py "$LIB" 120 170
run s6_twisted_reduction_120_40 s6_twisted_reduction.py "$LIB" 120 40
run s6_twisted_reduction_200_300 s6_twisted_reduction.py "$LIB" 200 300
run s6_margins_and_rows s6_margins_and_rows.py "$LIB"
cat outputs/*.txt | grep -E 'RESULT|falsifier'
exit "$status"
