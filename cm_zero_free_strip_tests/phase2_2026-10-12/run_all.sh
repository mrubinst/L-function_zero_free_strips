#!/bin/bash
# Phase 2 (2026-10-12): referee checks (r1-r3) and checks of the six Part II appendices (s1-s6).
# Needs python3 (numpy, sympy, mpmath, scipy) and PARI/GP. About 15-20 minutes. Outputs: outputs/<suite>.txt
# Exit status 1 if any script or sub-runner exits nonzero, or prints a failure line.
cd "$(dirname "$0")" || exit 1
mkdir -p outputs
export PYTHON="${PYTHON:-python3}"
status=0
# Failure detection (2026-10-12).  Primary: the exit status of each command (every test script exits 1 when one of
# its checks fails or a falsifier does not fail).  Secondary: a scan of the output for failure lines in any of the
# formats used in tests_cm: [FAIL], RESULT|..|FAIL, a bare FAIL at the start of a line, SOME FAILED / SOME FAIL,
# DID NOT FAIL, EXIT STATUS 1, "<n> FAIL" / "<n> failures" with n >= 1.  PASS lines and "fails as required" are ignored,
# so "0 FAIL" and "0 failures" do not count.
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
nfail_lines() { grep -E "$FAIL_PAT" "$1" | grep -v -E "$OK_PAT" | wc -l | tr -d ' '; }
run() { name=$1; dir=$2; shift 2; echo "== $name"; ( cd "$dir" && "$@" ) < /dev/null > "outputs/$name.txt" 2>&1; rc=$?
        nf=$(nfail_lines "outputs/$name.txt")
        echo "   exit $rc, failure lines $nf"; if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi; }
run r1_classgroup   r1 "$PYTHON" -I test_classgroup.py
run r1_extraction   r1 "$PYTHON" -I test_extraction_rows.py
run r1_bookkeeping  r1 "$PYTHON" -I test_step45_bookkeeping.py
run r2_t2           r2 "$PYTHON" -B t2_grouping_cusp_data.py
run r2_t1           r2 "$PYTHON" -B t1_prop_grouped_endtoend.py
run r2_t3           r2 "$PYTHON" -B t3_gl_family.py
run r3_t1           r3 "$PYTHON" t1_gauss_and_curves.py
run r3_t2           r3 "$PYTHON" t2_higher_weight.py
run r3_t3           r3 "$PYTHON" t3_final_inequality.py
run s1_local        s1 "$PYTHON" t_s1_local_phases.py
run s1_endtoend_g0  s1 "$PYTHON" t_s1_endtoend_active.py gauss 0
run s1_endtoend_g1  s1 "$PYTHON" t_s1_endtoend_active.py gauss 1
run s1_endtoend_b0  s1 "$PYTHON" t_s1_endtoend_active.py bessel 0
run s1_endtoend_b1  s1 "$PYTHON" t_s1_endtoend_active.py bessel 1
run s1_kernel       s1 "$PYTHON" t_s1_kernel.py
run s2_t1           s2 "$PYTHON" t1_theta_fe.py
run s2_t2           s2 "$PYTHON" t2_growth.py
run s2_t3           s2 "$PYTHON" t3_reflection.py
run s2_t4           s2 "$PYTHON" t4_conductors.py
run s3_all          s3 ./run_all.sh
run s4_checks       s4 "$PYTHON" s4_checks.py
run s5_all          s5 ./run_all.sh
run s6_all          s6 ./run_all.sh
echo "overall status $status"; exit $status
