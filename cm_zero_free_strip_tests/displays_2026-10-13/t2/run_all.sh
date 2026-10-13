#!/bin/bash
# t2: tests of the numbered displays (7.1)-(9.4) of cm-eleven-twelfths-guide.tex (Sections 7-9), 2026-10-13.
#   ./run_all.sh           new tests t2a-t2d, then reruns of the existing tests that cover the other displays (~8 min)
#   ./run_all.sh --quick   skip the 4-minute end-to-end reflection test (part B t5, eq:weight)
# New-test outputs go to outputs/, reruns to outputs/reruns/ (nothing is written into the other suites: python3 -B).
# Exit status 1 if any script exits nonzero or prints a failure line (same patterns as the other tests_cm runners).
cd "$(dirname "$0")" || exit 1
HERE=$PWD; TC=$(cd ../.. && pwd); OUT=$HERE/outputs; RR=$OUT/reruns
mkdir -p "$OUT" "$RR"
PY="${PYTHON:-python3}"
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
nfail_lines() { grep -E "$FAIL_PAT" "$1" | grep -v -E "$OK_PAT" | wc -l | tr -d ' '; }
status=0; T0=$(date +%s)
run() { # run <outfile> <dir> <cmd...>
  out=$1; dir=$2; shift 2; t=$(date +%s)
  ( cd "$dir" && "$@" ) < /dev/null > "$out" 2>&1; rc=$?
  nf=$(nfail_lines "$out")
  printf '%-45s exit %d, failure lines %s, %ds\n' "$(basename "$out")" "$rc" "$nf" "$(( $(date +%s) - t ))"
  if [ $rc -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi; }

echo "== new tests"
run "$OUT/t2a_gauss_identities.txt"   "$HERE" "$PY" -B t2a_gauss_identities.py
run "$OUT/t2b_sxi_endtoend.txt"       "$HERE" "$PY" -B t2b_sxi_endtoend.py
run "$OUT/t2c_change_of_variables.txt" "$HERE" "$PY" -B t2c_change_of_variables.py
run "$OUT/t2d_sizecheck_hyp.txt"      "$HERE" "$PY" -B t2d_sizecheck_hyp.py

echo "== reruns of existing tests (unchanged scripts, read-only)"
PB=$TC/part_B_lines_0550_0997
run "$RR/partB_t2_poisson.txt"  "$PB" "$PY" -B t2_poisson.py
run "$RR/partB_t3_gauss.txt"    "$PB" "$PY" -B t3_gauss.py
run "$RR/partB_t4_theta.txt"    "$PB" "$PY" -B t4_theta.py
run "$RR/partB_t6_sieve.txt"    "$PB" "$PY" -B t6_sieve.py
if [ "$1" != "--quick" ]; then
  run "$RR/partB_t5_reflection_j2.txt" "$PB" "$PY" -B t5_reflection_j2.py
fi
PR=$TC/proofs_2026-10-12
run "$RR/b2_reduction_170.txt"      "$TC" "$PY" -B proofs_2026-10-12/b2/b2_reduction_check.py part_B_lines_0550_0997 120 170
run "$RR/b2_reduction_40.txt"       "$TC" "$PY" -B proofs_2026-10-12/b2/b2_reduction_check.py part_B_lines_0550_0997 120 40
run "$RR/b2_reduction_D300_H360.txt" "$TC" "$PY" -B proofs_2026-10-12/b2/b2_reduction_check.py part_B_lines_0550_0997 300 360
run "$RR/b2_normalization.txt"      "$TC" "$PY" -B proofs_2026-10-12/b2/b2_normalization_check.py part_B_lines_0550_0997
run "$RR/b3_symbolic.txt"           "$PR" "$PY" -B b3/b3_symbolic_checks.py
run "$RR/b4_checks.txt"             "$PR" "$PY" -B b4/b4_checks.py
run "$RR/r1_test_classgroup.txt"    "$TC/phase2_2026-10-12/r1" "$PY" -I -B test_classgroup.py
run "$RR/r1_test_step45_bookkeeping.txt" "$TC/phase2_2026-10-12/r1" "$PY" -I -B test_step45_bookkeeping.py
run "$RR/r2_t2_grouping_cusp_data.txt" "$TC/phase2_2026-10-12/r2" "$PY" -B t2_grouping_cusp_data.py
run "$RR/r2_t3_gl_family.txt"       "$TC/phase2_2026-10-12/r2" "$PY" -B t3_gl_family.py
run "$RR/sec10_descent.txt"         "$TC/section10_descent" "$PY" -B test_sec10_descent.py
echo "total $(( $(date +%s) - T0 )) s; overall status $status"
exit $status
