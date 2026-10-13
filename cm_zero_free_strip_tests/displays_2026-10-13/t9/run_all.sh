#!/bin/sh
# t9 display tests (Appendix F of the CM guide).  Runs the six new scripts d1-d6 and re-runs, read-only, the six
# scripts of tests_cm/phase2_2026-10-12/s5 (their output goes to outputs/rerun_s5_*.txt here; nothing is written
# into the s5 folder).  Exit status 1 if any script exits nonzero or prints a failure line.
cd "$(dirname "$0")" || exit 1
HERE=$(pwd)
S5="$HERE/../../phase2_2026-10-12/s5"
mkdir -p outputs
export PYTHON="${PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1
status=0
FAIL_PAT='^RESULT[|][^|]*[|]FAIL|DID NOT FAIL|EXIT STATUS [1-9]|Traceback'
start=$(date +%s)
run() {   # $1 = label, $2 = directory, $3 = script, $4 = output file
  t0=$(date +%s)
  (cd "$2" && "$PYTHON" "$3") > "$4" 2>&1
  rc=$?
  t1=$(date +%s)
  nf=$(grep -E -c "$FAIL_PAT" "$4")
  echo "$1: rc=$rc failure-lines=$nf time=$((t1 - t0))s  $(grep SUMMARY "$4" | tail -1)"
  if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; echo "[FAIL] $1"; fi
}
for t in d1_alpha_mult d2_rowchar d3_poisson d4_mellin d5_lattice d6_fourth_moment; do
  run "$t" "$HERE" "$t.py" "outputs/$t.txt"
done
for t in t1_dual_profile t2_hankel t3_reflection t4_lattice t5_row_reflection t6_kappa; do
  run "s5/$t (rerun)" "$S5" "$t.py" "$HERE/outputs/rerun_s5_$t.txt"
done
echo "total time $(( $(date +%s) - start ))s; status $status"
exit "$status"
