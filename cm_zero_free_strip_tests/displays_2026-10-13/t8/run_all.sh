#!/bin/sh
# Agent t8 (displays of Appendices E and G of cm-eleven-twelfths-guide.tex): runs the new t8 tests and re-runs the
# existing suites that test displays of these appendices (tests_cm/phase2_2026-10-12/s4 and s6), writing every output
# to t8/outputs/ (the other suites' folders are not written to).  Exit status 1 if any script exits nonzero or
# prints a failure line.  About 70 s.
cd "$(dirname "$0")" || exit 1
mkdir -p outputs
PY="${PYTHON:-python3}"
S4=../../phase2_2026-10-12/s4
S6=../../phase2_2026-10-12/s6
LIB=../../part_B_lines_0550_0997
status=0
FAIL_PAT='^[[:space:]]*\[(FAIL|SANITY-FAIL)\]|^RESULT[|][^|]*[|]FAIL|DID NOT FAIL|EXIT STATUS [1-9]|<-- PROBLEM'
run() {
  name=$1; shift
  t0=$(date +%s)
  "$PY" "$@" > "outputs/$name.txt" 2>&1
  rc=$?
  t1=$(date +%s)
  nf=$(grep -E "$FAIL_PAT" "outputs/$name.txt" | grep -v 'fails as required' | wc -l | tr -d ' ')
  if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; echo "[FAIL] $name: exit $rc, failure lines $nf ($((t1-t0)) s)";
  else echo "[PASS] $name: exit 0 ($((t1-t0)) s)"; fi
}
T0=$(date +%s)
# new tests
run t8_coefficients   t8_coefficients.py
run t8_signal_init    t8_signal_init.py
run t8_beta_block     t8_beta_block.py
run t8_euler_mellin   t8_euler_mellin.py
run t8_exponents      t8_exponents.py
run t8_tex_compare    t8_tex_compare.py
run t8_size_checks    t8_size_checks.py
# re-runs of existing suites (read-only; outputs here)
run rerun_s4_checks   "$S4/s4_checks.py"
run rerun_s6_twisted_paired              "$S6/s6_twisted_paired.py" "$LIB"
run rerun_s6_twisted_reduction_120_170   "$S6/s6_twisted_reduction.py" "$LIB" 120 170
run rerun_s6_twisted_reduction_120_40    "$S6/s6_twisted_reduction.py" "$LIB" 120 40
run rerun_s6_twisted_reduction_200_300   "$S6/s6_twisted_reduction.py" "$LIB" 200 300
run rerun_s6_margins_and_rows            "$S6/s6_margins_and_rows.py" "$LIB"
T1=$(date +%s)
echo "total wall time $((T1-T0)) s; status $status"
exit "$status"
