"""Run all Part D tests (lines 1510-end of cm-eleven-twelfths-guide.tex) and print one PASS/FAIL line per claim.

Usage:  python3 run_all.py            (all tests; t08 takes about 3-4 minutes)
        python3 run_all.py --fast     (skip t08, the end-to-end reflection with B = 80000)
Outputs of each script are also saved in outputs/<script>.txt.
Exit status 1 if any claim fails or any script exits nonzero (every script exits 1 when one of its checks fails).
"""
import os, sys, subprocess, time, re
# Secondary check (added 2026-10-12): the exit status of each script is primary; in addition the output is scanned for
# failure lines in any of the formats used in tests_cm ([FAIL], RESULT|..|FAIL, bare FAIL, SOME FAILED, DID NOT FAIL,
# "<n> FAIL"/"<n> failures" with n >= 1, EXIT STATUS 1), ignoring PASS lines and "fails as required".
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
def failure_lines(text):
    return [l for l in text.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["t01_eta_quotients.py", "t02_chain_rule.py", "t03_mellin.py", "t04_kernels.py", "t05_misc.py",
           "t06_sympowers.py", "t07_theta_appendix.py", "t08_reflection_j.py"]
if "--fast" in sys.argv:
    SCRIPTS = SCRIPTS[:-1]
os.makedirs(os.path.join(HERE, "outputs"), exist_ok=True)
npass = nfail = 0; failed_scripts = []
for s in SCRIPTS:
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, s)], capture_output=True, text=True, cwd=HERE)
    out = r.stdout + (("\n[stderr]\n" + r.stderr) if r.returncode else "")
    open(os.path.join(HERE, "outputs", s.replace(".py", ".txt")), "w").write(out)
    for line in r.stdout.splitlines():
        if line.startswith("[PASS]") or line.startswith("[FAIL]"):
            print(line, flush=True)
            if line.startswith("[PASS]"): npass += 1
            else: nfail += 1
    if r.returncode:
        failed_scripts.append(s); print(f"[FAIL] {s} exited with code {r.returncode}: {r.stderr.strip().splitlines()[-1:]}")
        nfail += 1
    elif failure_lines(r.stdout):
        failed_scripts.append(s); print(f"[FAIL] {s} printed failure lines: {failure_lines(r.stdout)[:2]}")
        nfail += 1
    print(f"[INFO] {s}: {time.time() - t0:.0f}s", flush=True)
print(f"\nSUMMARY: {npass} PASS, {nfail} FAIL" + (f"; failed scripts: {failed_scripts}" if failed_scripts else ""))
sys.exit(1 if nfail else 0)
