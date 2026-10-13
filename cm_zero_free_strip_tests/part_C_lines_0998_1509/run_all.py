"""Run every test of the Part C audit (guide lines 998-1509) and print one line per claim.

    python3 run_all.py            # all scripts (about 4 minutes; t04 takes ~90 s, uses PARI/GP)
    python3 run_all.py --quick    # skip t04 (brute-force conductors) and t08 (argument principle)

Statuses: PASS, FAIL, CITED (restatement matches OAI78/G78), CITED-X (does not), INFO (sanity numbers, no criterion),
UNTESTED (reason given).  Outputs of each script are saved in outputs/<script>.txt.  Exit status 1 on FAIL or CITED-X,
or if any script exits nonzero (every script exits 1 when one of its checks fails).
"""
import os, sys, subprocess, time, re
from collections import Counter
# Secondary check (added 2026-10-12): the exit status of each script is primary; in addition the output is scanned for
# failure lines in any of the formats used in tests_cm ([FAIL], RESULT|..|FAIL, bare FAIL, SOME FAILED, DID NOT FAIL,
# "<n> FAIL"/"<n> failures" with n >= 1, EXIT STATUS 1), ignoring PASS lines and "fails as required".
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
def failure_lines(text):
    return [l for l in text.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["t01_curves.py", "t02_gauss.py", "t03_reflection.py", "t04_growth_fe.py", "t05_rows_moments.py",
           "t06_final_inequality.py", "t07_source_crossrefs.py", "t08_theorem_sanity.py", "t09_untested.py"]
quick = "--quick" in sys.argv
os.makedirs(os.path.join(HERE, "outputs"), exist_ok=True)
lines = []; crashed = []
for s in SCRIPTS:
    if quick and s in ("t04_growth_fe.py", "t08_theorem_sanity.py"):
        continue
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, s)], capture_output=True, text=True, cwd=HERE)
    with open(os.path.join(HERE, "outputs", s.replace(".py", ".txt")), "w") as f:
        f.write(r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr.strip() else ""))
    got = [l for l in r.stdout.split("\n") if l.startswith("[")]
    if r.returncode != 0:
        crashed.append(s)
        got.append(f"[FAIL] EXIT {s} :: exit status {r.returncode}: {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}")
    elif failure_lines(r.stdout):
        crashed.append(s)
        got.append(f"[FAIL] OUTPUT {s} :: failure lines: {failure_lines(r.stdout)[0][:100]}")
    lines += got
    print(f"# {s}: {len(got)} claims, {time.time() - t0:.1f} s", flush=True)
    for l in got:
        print(l, flush=True)

st = Counter(re.match(r"\[([A-Z-]+)\]", l).group(1) for l in lines)
cls = Counter((re.search(r"\{([A-Z])\}", l) or re.search(r"(x)", "x")).group(1) for l in lines)
print("\n# summary by status:", dict(st))
print("# summary by class (E exact, I inequality/exponent, A asymptotic/analytic, C cited):", dict(cls))
sys.exit(1 if st.get("FAIL", 0) or st.get("CITED-X", 0) or crashed else 0)
