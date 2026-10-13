"""t4: tests of the numbered displays (11.1)-(15.2) and (17.1)-(17.9) of cm-eleven-twelfths-guide.tex, 2026-10-13.

    python3 run_all.py          # new tests t4_*.py, then the reruns of the existing tests (about 3.5 minutes)
    python3 run_all.py --new    # only the new tests (about 25 s)

Outputs: outputs/<script>.txt; the reruns also write outputs/rerun_<suite>_<script>.txt (nothing is written into the
other suites: PYTHONDONTWRITEBYTECODE=1, their stdout is captured here).  Exit status 1 if any script exits nonzero or
prints a failure line (same patterns as the other tests_cm runners).  The row table is rows.json.
"""
import os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs"); os.makedirs(OUT, exist_ok=True)
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
SCRIPTS = ["t4_sources.py", "t4_local_euler.py", "t4_lfunctions.py", "t4_moments.py"]
if "--new" not in sys.argv:
    SCRIPTS.append("rerun_existing.py")
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
bad = 0; T0 = time.time()
for s in SCRIPTS:
    t = time.time()
    r = subprocess.run([sys.executable, s], cwd=HERE, env=env, capture_output=True, text=True)
    txt = r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr.strip() else "")
    open(os.path.join(OUT, s.replace(".py", ".txt")), "w").write(txt)
    fl = [l for l in txt.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]
    bad += (r.returncode != 0) or bool(fl)
    print(f"{s:24s} exit {r.returncode}, failure lines {len(fl)}, {time.time() - t:.1f} s", flush=True)
    for l in r.stdout.splitlines():
        if l.startswith("["):
            print("   " + l[:200])
print(f"# total {time.time() - T0:.1f} s; {'ALL PASS' if not bad else f'{bad} script(s) FAILED'}")
sys.exit(1 if bad else 0)
