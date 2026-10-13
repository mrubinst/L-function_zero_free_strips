"""t4 (guide Sections 11-15 and 17): re-run the existing tests that check the displays of this inventory.

The scripts of the other suites are run read-only (cwd = their folder, PYTHONDONTWRITEBYTECODE=1); their stdout is
saved in outputs/rerun_<suite>_<script>.txt here (never in the other suites' outputs/).  A display counts as
re-tested only if (a) the script exits with status 0 and (b) every test id listed for it is reported with status
PASS (INFO for the sanity lines marked as such).  Exit status 1 otherwise.

    python3 rerun_existing.py            # all (about 4 minutes wall clock, scripts run in parallel)
"""
import os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, "..", ".."))           # tests_cm/
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

C = "part_C_lines_0998_1509"
D = "part_D_lines_1510_end"
R3 = "phase2_2026-10-12/r3"

# (suite, script) -> list of (test id, required status)
NEED = {
    (C, "t01_curves.py"): [("C02", "PASS"), ("C03", "PASS"), ("C04", "PASS"), ("C05", "PASS"), ("C07", "PASS"),
                           ("C08", "PASS"), ("C09", "PASS")],
    (C, "t02_gauss.py"): [("C13", "PASS"), ("C14", "PASS"), ("C15", "PASS"), ("C17", "PASS"), ("C19", "PASS"),
                          ("C23", "PASS")],
    (C, "t04_growth_fe.py"): [("C37", "PASS"), ("C38", "PASS"), ("C39", "PASS"), ("C43", "PASS")],
    (C, "t05_rows_moments.py"): [("C45", "PASS"), ("C46", "PASS"), ("C103", "INFO")],
    (C, "t06_final_inequality.py"): [("C58", "PASS"), ("C59", "PASS"), ("C60", "PASS")],
    (C, "t08_theorem_sanity.py"): [("C98", "PASS")],
    (D, "t01_eta_quotients.py"): [("D01.11", "PASS"), ("D01.12", "PASS"), ("D01.13", "PASS")],
    (D, "t02_chain_rule.py"): [("D02.1", "PASS"), ("D02.3", "PASS")],
    (D, "t03_mellin.py"): [("D03.1", "PASS"), ("D03.3", "PASS")],
    (D, "t04_kernels.py"): [("D04.1", "PASS"), ("D04.2", "PASS"), ("D04.5", "PASS")],
    (D, "t05_misc.py"): [("D05.2", "PASS"), ("D05.7", "PASS")],
    (D, "t06_sympowers.py"): [("D06.1", "PASS"), ("D06.2", "PASS")],
    (R3, "t1_gauss_and_curves.py"): [("G2", "PASS"), ("G3", "PASS"), ("G4", "PASS"), ("C1", "PASS"), ("C2", "PASS")],
    (R3, "t2_higher_weight.py"): [("H1", "PASS"), ("H2", "PASS"), ("H3", "PASS"), ("H6", "PASS"), ("H7", "PASS")],
    (R3, "t3_final_inequality.py"): [("E4", "PASS"), ("E4f", "PASS"), ("E5", "PASS"), ("E6", "PASS")],
}

def run(key):
    suite, script = key
    d = os.path.join(TESTS, suite)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    t0 = time.time()
    r = subprocess.run([sys.executable, script], cwd=d, env=env, capture_output=True, text=True)
    dt = time.time() - t0
    tag = suite.replace("/", "_").replace("phase2_2026-10-12_", "")
    with open(os.path.join(OUT, f"rerun_{tag}_{script.replace('.py', '.txt')}"), "w") as f:
        f.write(r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr.strip() else "") + f"\n--- exit {r.returncode}, {dt:.1f} s\n")
    return key, r.returncode, r.stdout, dt

T0 = time.time()
with ThreadPoolExecutor(max_workers=6) as ex:
    results = list(ex.map(run, NEED))
bad = 0
for (suite, script), rc, out, dt in results:
    lines = out.splitlines()
    for tid, st in NEED[(suite, script)]:
        pat = re.compile(r"^\[(\w[\w-]*)\]\s+" + re.escape(tid) + r"\s")
        hits = [m.group(1) for l in lines for m in [pat.match(l)] if m]
        ok = rc == 0 and hits and all(h == st for h in hits)
        bad += not ok
        print(f"[{'PASS' if ok else 'FAIL'}] {suite}/{script} {tid}: reported {hits or 'nothing'} (need {st}), "
              f"exit {rc}, {dt:.1f} s")
print(f"# wall clock {time.time() - T0:.1f} s; CPU sum {sum(r[3] for r in results):.1f} s; {bad} failure(s)")
sys.exit(1 if bad else 0)
