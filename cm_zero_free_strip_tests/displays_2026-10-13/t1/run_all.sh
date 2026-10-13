#!/bin/bash
# Agent t1: every test behind rows.json (displays 0.1 and 2.1-6.5 of the CM guide).
#   1. the new tests of this folder (t1_*.py), outputs in outputs/<name>.txt;
#   2. rerun_existing.sh: the existing tests that cover the other displays, re-run read-only, outputs in
#      outputs/existing_<name>.txt.
# Exit status 1 if any script exits nonzero or prints a failure line.  CM_TESTS_FORCE_FAIL=1 makes the first check of
# each new script fail (self-test of the wiring).  Needs python3 (numpy, sympy, mpmath) and PARI/GP.
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/outputs"; mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
PY="${PYTHON:-python3}"
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
status=0
T0=$(date +%s)
for t in t1_pointcount_steps t1_euler_chain t1_mellin t1_sizes; do
  t1=$(date +%s)
  ( cd "$HERE" && "$PY" "$t.py" ) < /dev/null > "$OUT/$t.txt" 2>&1; rc=$?
  nf=$(grep -E "$FAIL_PAT" "$OUT/$t.txt" | grep -v -E "$OK_PAT" | wc -l | tr -d ' ')
  echo "== $t: exit $rc, failure lines $nf, $(( $(date +%s) - t1 )) s"
  grep -E '^\[(PASS|FAIL|SANITY)\]' "$OUT/$t.txt" | cut -c1-160
  if [ $rc -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi
done
if [ "${CM_TESTS_FORCE_FAIL:-0}" != "1" ]; then
  # PIPESTATUS[0] is the status of rerun_existing.sh itself, not of tee (Codex review, 2026-10-13)
  "$HERE/rerun_existing.sh" | tee "$OUT/rerun_existing_summary.txt"
  [ "${PIPESTATUS[0]}" -eq 0 ] || status=1
fi
echo "agent t1: total $(( $(date +%s) - T0 )) s, overall status $status"
exit $status
