#!/bin/bash
# Run every check of the proofs added to the CM guide on 2026-10-12 (blocks b1-b4). Outputs go to outputs/.
# Exit status 1 if any script exits nonzero or prints a failure line.
cd "$(dirname "$0")" || exit 1
mkdir -p outputs
export PYTHON="${PYTHON:-python3}"
LIBB=../part_B_lines_0550_0997
status=0
# Failure detection (2026-10-12).  Primary: the exit status of each command (every test script exits 1 when one of
# its checks fails or a falsifier does not fail).  Secondary: a scan of the output for failure lines in any of the
# formats used in tests_cm: [FAIL], RESULT|..|FAIL, a bare FAIL at the start of a line, SOME FAILED / SOME FAIL,
# DID NOT FAIL, EXIT STATUS 1, "<n> FAIL" / "<n> failures" with n >= 1.  PASS lines and "fails as required" are ignored,
# so "0 FAIL" and "0 failures" do not count.
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
nfail_lines() { grep -E "$FAIL_PAT" "$1" | grep -v -E "$OK_PAT" | wc -l | tr -d ' '; }
run() { name=$1; shift; echo "== $name"; ( "$@" ) < /dev/null > outputs/$name.txt 2>&1; rc=$?; tail -n 25 outputs/$name.txt
        nf=$(nfail_lines outputs/$name.txt); echo "-- exit $rc, failure lines $nf"
        if [ $rc -ne 0 ] || [ $nf -ne 0 ]; then status=1; fi; }
run b1_points      "$PYTHON" b1/t_points.py
run b1_extract     "$PYTHON" b1/t_extract.py
run b1_bd          gp -q b1/t_bd.gp
run b2_reduction_170 "$PYTHON" b2/b2_reduction_check.py $LIBB 120 170
run b2_reduction_40  "$PYTHON" b2/b2_reduction_check.py $LIBB 120 40
run b2_normalization "$PYTHON" b2/b2_normalization_check.py $LIBB
run b3_symbolic    "$PYTHON" b3/b3_symbolic_checks.py
run b3_endtoend    "$PYTHON" b3/t_b3_prop_endtoend.py
run b4_checks      "$PYTHON" b4/b4_checks.py
run b4_weight_decay "$PYTHON" b4/b4_weight_decay.py
echo "overall status $status"
exit $status
