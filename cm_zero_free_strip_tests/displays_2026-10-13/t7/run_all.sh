#!/bin/sh
# Agent t7 display tests (CM guide, Appendices C and D).  Runs the new tests of this folder, then (unless
# T7_SKIP_RERUN=1) re-runs the existing phase-2 suites s2 (t1-t4) and s3 (t1-t5) through rerun_existing.sh.
# Needs python3 with numpy, scipy, sympy, mpmath, and PARI/GP (gp).  Outputs: outputs/*.txt.
# Exit status 1 if any script exits nonzero or prints a [FAIL] line.
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE" || exit 1
mkdir -p outputs
export PYTHONDONTWRITEBYTECODE=1
PY="${PYTHON:-python3}"
status=0
/bin/rm -f outputs/run_all_summary.txt
T0=$(date +%s)
for t in c1_pres_gauss c2_fourier c3_theta_mellin c4_fe_pari d1_members; do
  start=$(date +%s)
  "$PY" "$t.py" > "outputs/$t.txt" 2>&1
  rc=$?
  end=$(date +%s)
  nf=$(grep -c '^\[FAIL\]' "outputs/$t.txt")
  echo "$t: exit $rc, [FAIL] lines $nf, $((end - start)) s" | tee -a outputs/run_all_summary.txt
  grep '^\[' "outputs/$t.txt" | cut -c1-150
  if [ "$rc" -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi
done
if [ "${T7_SKIP_RERUN:-0}" != "1" ]; then
  /bin/rm -f outputs/rerun_summary.txt
  sh ./rerun_existing.sh || status=1
  cat outputs/rerun_summary.txt >> outputs/run_all_summary.txt
fi
T1=$(date +%s)
echo "run_all: overall exit $status, total $((T1 - T0)) s" | tee -a outputs/run_all_summary.txt
exit $status
