#!/bin/sh
# Re-run (read-only) the existing phase-2 suites that test displays of Appendices C and D of the CM guide.
# The scripts are not modified; their outputs go to outputs/rerun_<suite>_<test>.txt here.
# Exit status 1 if any script exits nonzero.
HERE="$(cd "$(dirname "$0")" && pwd)"
P2="$HERE/../../phase2_2026-10-12"
OUT="$HERE/outputs"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
PY="${PYTHON:-python3}"
status=0
for t in s2/t1_theta_fe s2/t2_growth s2/t3_reflection s2/t4_conductors \
         s3/t1_labels s3/t2_conductors s3/t3_fe s3/t4_gamma s3/t5_citations; do
  name=$(echo "$t" | tr '/' '_')
  start=$(date +%s)
  "$PY" "$P2/$t.py" > "$OUT/rerun_$name.txt" 2>&1
  rc=$?
  end=$(date +%s)
  echo "$t: exit $rc, $((end - start)) s" | tee -a "$OUT/rerun_summary.txt"
  if [ "$rc" -ne 0 ]; then status=1; fi
done
echo "rerun_existing: overall exit $status" | tee -a "$OUT/rerun_summary.txt"
exit $status
