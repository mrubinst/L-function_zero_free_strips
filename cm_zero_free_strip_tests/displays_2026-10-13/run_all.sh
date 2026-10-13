#!/bin/bash
# Run the display tests of all ten groups t0-t9. Exit status 1 if any group fails.
# Each group's README names the script behind every display (rows.json has one row per display).
cd "$(dirname "$0")" || exit 2
status=0
run() { echo "===== $1: $2"; ( cd "$1" && eval "$2" ); rc=$?; echo "===== $1 exit $rc"; [ $rc -ne 0 ] && status=1; }
run t0 ./run_all.sh
run t1 ./run_all.sh
run t2 ./run_all.sh
run t3 "python3 run_t3.py"
run t4 "python3 run_all.py"
run t5 "python3 run_all.py --reruns"   # --reruns: also rerun the eleven existing scripts t5 relies on
run t6 ./run_all.sh
run t7 ./run_all.sh
run t8 ./run_all.sh
run t9 ./run_all.sh
echo "===== overall exit $status"
exit $status
