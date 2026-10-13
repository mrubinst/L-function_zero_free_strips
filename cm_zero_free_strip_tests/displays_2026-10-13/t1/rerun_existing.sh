#!/bin/bash
# Agent t1 (displays of the front matter and Sections 2-6 of the CM guide): re-run, unmodified and read-only, the
# existing tests that cover these displays. Output of each run goes to outputs/existing_<name>.txt in THIS folder
# (the other suites' own outputs/ folders are not touched; PYTHONDONTWRITEBYTECODE=1 avoids writing caches there).
# Exit status 1 if any script exits nonzero or prints a failure line (same failure patterns as the other suites).
HERE="$(cd "$(dirname "$0")" && pwd)"
TC="$(cd "$HERE/../.." && pwd)"            # tests_cm/
OUT="$HERE/outputs"; mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
PY="${PYTHON:-python3}"
FAIL_PAT='^[[:space:]]*\[(FAIL|CITED-X)\]|^RESULT[|][^|]*[|]FAIL|^[[:space:]]*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)'
OK_PAT='^[[:space:]]*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT[|][^|]*[|]PASS[|]|PASS([^A-Za-z_]|$))|fails as required'
status=0
run() { name=$1; dir=$2; shift 2; t0=$(date +%s)
        ( cd "$dir" && "$@" ) < /dev/null > "$OUT/existing_$name.txt" 2>&1; rc=$?
        nf=$(grep -E "$FAIL_PAT" "$OUT/existing_$name.txt" | grep -v -E "$OK_PAT" | wc -l | tr -d ' ')
        t1=$(date +%s); echo "== $name: exit $rc, failure lines $nf, $((t1-t0)) s"
        if [ $rc -ne 0 ] || [ "$nf" -ne 0 ]; then status=1; fi; }
T0=$(date +%s)
A="$TC/part_A_lines_0001_0549"
run A_t0_intro    "$A" "$PY" t0_intro.py
run A_t2_integer  "$A" "$PY" t2_integer.py
run A_t3_lattice  "$A" "$PY" t3_lattice.py
run A_t4_average  "$A" "$PY" t4_average.py
run A_t5_sixth    "$A" "$PY" t5_sixth.py
run B_t1_zeros    "$TC/part_B_lines_0550_0997" "$PY" t1_zeros.py
P="$TC/proofs_2026-10-12"
run b1_points     "$P" "$PY" b1/t_points.py
run b1_extract    "$P" "$PY" b1/t_extract.py
run b1_bd         "$P" gp -q b1/t_bd.gp
run r1_extraction_rows "$TC/phase2_2026-10-12/r1" "$PY" -I test_extraction_rows.py
run deuring_identity "$TC/deuring_identity" "$PY" hecke_identity.py
echo "re-run of existing tests: total $(( $(date +%s) - T0 )) s, overall status $status"
exit $status
