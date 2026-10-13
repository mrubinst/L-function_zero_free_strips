#!/bin/bash
# t6 (2026-10-13): re-run the existing tests of Appendix app:s1-general (phase2 s1) and of prop:b3-reflection (proofs b3).
# The scripts are run read-only in their own folders (python -B: no bytecode written there); outputs go to
# out/rerun_<name>.txt, exit codes to out/rerun_<name>.rc.  Exit status 1 if any script exits nonzero.
HERE="$(cd "$(dirname "$0")" && pwd)"; TC="$(cd "$HERE/../.." && pwd)"
PY="${PYTHON:-python3}"; mkdir -p "$HERE/out"
S1="$TC/phase2_2026-10-12/s1"; B3="$TC/proofs_2026-10-12/b3"
run() { name=$1; dir=$2; shift 2; t0=$(date +%s)
        ( cd "$dir" && "$PY" -B "$@" ) < /dev/null > "$HERE/out/rerun_$name.txt" 2>&1; rc=$?
        echo $rc > "$HERE/out/rerun_$name.rc"
        echo "rerun_$name: exit $rc, $(( $(date +%s)-t0 )) s"; }
run s1_local       "$S1" t_s1_local_phases.py &
run s1_kernel      "$S1" t_s1_kernel.py &
run s1_endtoend_g0 "$S1" t_s1_endtoend_active.py gauss 0 &
run s1_endtoend_g1 "$S1" t_s1_endtoend_active.py gauss 1 &
run s1_endtoend_b0 "$S1" t_s1_endtoend_active.py bessel 0 &
run s1_endtoend_b1 "$S1" t_s1_endtoend_active.py bessel 1 &
run b3_symbolic    "$B3" b3_symbolic_checks.py &
run b3_endtoend    "$B3" t_b3_prop_endtoend.py &
wait
status=0
for f in "$HERE"/out/rerun_*.rc; do [ "$(cat "$f")" = "0" ] || { echo "NONZERO EXIT: $f"; status=1; }; done
echo "rerun overall status $status"; exit $status
