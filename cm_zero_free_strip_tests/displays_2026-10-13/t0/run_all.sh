#!/bin/bash
# t0: tests of the ten numbered displays of the Introduction, (0.1)-(0.10). Exit status 1 if the script fails.
cd "$(dirname "$0")" || exit 2
PY="${PYTHON:-python3}"; export PYTHONDONTWRITEBYTECODE=1
mkdir -p outputs
"$PY" t0_intro.py < /dev/null > outputs/t0_intro.txt 2>&1; rc=$?
grep -E '^\[(PASS|FAIL)\]|^SUMMARY' outputs/t0_intro.txt | cut -c1-170
nf=$(grep -cE '^\[FAIL\]' outputs/t0_intro.txt)
echo "t0: exit $rc, failure lines $nf"
[ $rc -eq 0 ] && [ "$nf" -eq 0 ]
