#!/bin/bash
# t6 (2026-10-13): all checks for the displays of Appendix app:s1-general (B.1-B.38).
#  1. re-runs of the existing suites (rerun_existing.sh: phase2 s1 x6, proofs b3 x2), in parallel;
#  2. the new scripts of this folder, in parallel.
# Outputs in out/.  Exit status 1 if any script exits nonzero.  Wall time about 4 minutes on 8+ cores.
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE" || exit 1
PY="${PYTHON:-python3}"; mkdir -p out
T0=$(date +%s)
./rerun_existing.sh > out/rerun_existing.log 2>&1 &
PIDR=$!
run() { name=$1; shift; t0=$(date +%s); "$PY" -B "$@" < /dev/null > "out/$name.txt" 2>&1; rc=$?
        echo $rc > "out/$name.rc"; echo "$name: exit $rc, $(( $(date +%s)-t0 )) s"; }
run t6_theta_translates_pi0 t6_theta_translates.py 0 &
run t6_theta_translates_pi1 t6_theta_translates.py 1 &
run t6_mellin_rep           t6_mellin_rep.py &
run t6_kernel_bookkeeping   t6_kernel_bookkeeping.py &
run t6_moment_sanity        t6_moment_sanity.py &
wait
cat out/rerun_existing.log
status=0
for f in out/*.rc; do [ "$(cat "$f")" = "0" ] || { echo "NONZERO EXIT: $f"; status=1; }; done
echo "t6 overall status $status; wall time $(( $(date +%s)-T0 )) s"; exit $status
