"""Runner for the Section 10 display checks of agent t3 (cm-eleven-twelfths-guide.tex, displays (10.1)-(10.41)).

Re-runs, read-only, the existing tests that cover Section 10, plus the new checks in t3_sec10_new.py:
  section10_descent/test_sec10_descent.py           S10.*  (old display numbers (10.1)-(10.23))
  proofs_2026-10-12/b4/b4_checks.py                 C1-C11 (Poisson identities of the transfer, regrouping, long cubes)
  part_B_lines_0550_0997/t7_descent.py              B68-B87
  part_B_lines_0550_0997/t2,t3,t4,t5                B17, B27, B46, B50, B53, B55 (the table (10.27) = B88)
  displays_2026-10-13/t3/t3_sec10_new.py            D10.*, S10.32, S10.34 (new)
The scripts run in parallel with PYTHONDONTWRITEBYTECODE=1; nothing is written into the other suites.  Each output goes
to outputs/<script>.txt here.  Exit status 1 if any script exits nonzero, prints a failure line, or does not report
PASS for one of the test ids that rows.json relies on.

Run: python3 run_t3.py            (about 3-4 minutes wall clock, dominated by t5_reflection_j2.py)
     python3 run_t3.py --no-table (skips t2-t5, i.e. the B88 constituents; about 1 minute)
"""
import os, sys, subprocess, time, re
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "outputs"); os.makedirs(OUT, exist_ok=True)
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
B = "part_B_lines_0550_0997"
JOBS = {   # name: (script relative to tests_cm, ids that must PASS)
    "test_sec10_descent": ("section10_descent/test_sec10_descent.py",
        ["S10.1", "S10.2-3", "S10.4", "S10.5-7", "S10.8", "S10.9", "S10.10", "S10.11a", "S10.11", "S10.11b",
         "S10.12-13", "S10.15", "S10.16-20", "S10.18", "S10.23", "S10.cite"]),
    "b4_checks": ("proofs_2026-10-12/b4/b4_checks.py", ["C1", "C2", "C3", "C4", "C5", "C6", "C10", "C10b", "C11"]),
    "t7_descent": (f"{B}/t7_descent.py", ["B68", "B69", "B70", "B72", "B74", "B75", "B76", "B77", "B78", "B79", "B80",
                                         "B81", "B82", "B83", "B84", "B85", "B86", "B87"]),
    "t3_sec10_new": ("displays_2026-10-13/t3/t3_sec10_new.py",
        ["D10.10", "D10.30", "D10.11", "D10.16", "D10.21", "D10.38", "D10.39", "D10.35", "D10.41", "D10.27",
         "S10.32", "S10.34"]),
    "t2_poisson": (f"{B}/t2_poisson.py", ["B17"]),
    "t3_gauss": (f"{B}/t3_gauss.py", ["B27"]),
    "t4_theta": (f"{B}/t4_theta.py", ["B46", "B50"]),
    "t5_reflection_j2": (f"{B}/t5_reflection_j2.py", ["B53", "B55"]),
}
if "--no-table" in sys.argv:
    for k in ("t2_poisson", "t3_gauss", "t4_theta", "t5_reflection_j2"): JOBS.pop(k)

def passed_ids(text):
    ids = set()
    for l in text.splitlines():
        m = re.match(r"^RESULT\|([^|]+)\|PASS\|", l)
        if m: ids.add(m.group(1))
        m = re.match(r"^\s*\[PASS\]\s+(\S+)", l)
        if m: ids.add(m.group(1))
    return ids

env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
t0 = time.time(); procs = {}
for name, (rel, _) in JOBS.items():
    path = os.path.join(TESTS, rel)
    fh = open(os.path.join(OUT, name + ".txt"), "w")
    procs[name] = (subprocess.Popen([sys.executable, path], cwd=os.path.dirname(path), stdout=fh, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, env=env), fh, time.time())
status = 0; lines = []
while procs:
    for name in list(procs):
        p, fh, ts = procs[name]
        if p.poll() is None: continue
        fh.close(); dt = time.time() - ts; del procs[name]
        text = open(os.path.join(OUT, name + ".txt")).read()
        fails = [l for l in text.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]
        ids = passed_ids(text); need = JOBS[name][1]; missing = [i for i in need if i not in ids]
        nfals = text.count("fails as required")
        ok = p.returncode == 0 and not fails and not missing
        status |= not ok
        lines.append(f"[{'PASS' if ok else 'FAIL'}] {name}: exit {p.returncode}, {len(need)} required ids "
                     f"({'all PASS' if not missing else 'missing/failed: ' + ', '.join(missing)}), {nfals} falsifiers failed "
                     f"as required, {len(fails)} failure lines, {dt:.0f} s")
        print(lines[-1], flush=True)
    time.sleep(0.5)
tot = time.time() - t0
summary = "\n".join(sorted(lines)) + f"\nwall clock {tot:.0f} s\n" + ("ALL PASS" if status == 0 else "SOME FAILURES") + "\n"
open(os.path.join(OUT, "run_t3_summary.txt"), "w").write(summary)
print(f"wall clock {tot:.0f} s"); print("ALL PASS" if status == 0 else "SOME FAILURES")
sys.exit(1 if status else 0)
