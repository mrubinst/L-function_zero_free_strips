"""Run every Part A test (lines 1-549 of cm-eleven-twelfths-guide.tex) and print one line per claim.

Status words:
  [PASS]      exact identity / explicit inequality verified on the stated range AND its falsifier failed
              (for analytic toy checks marked [A], the explicit finite inequality held);
  [FAIL]      a mismatch was found, or the falsifier was not detected;
  [SANITY]    analytic (A) claim: sizes printed, no pass/fail criterion;
  [UNTESTED]  not finitely testable here (main theorems, remarks, forward references, citations).
Usage:  python3 run_all.py            (about 1-2 minutes; outputs also written to outputs/*.txt)
Exit status 1 if any script exits nonzero (every script exits 1 when one of its checks fails) or prints a failure line.
"""
import os, sys, subprocess, time, re
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['t0_intro.py', 't2_integer.py', 't3_lattice.py', 't4_average.py', 't5_sixth.py']
# Secondary check (added 2026-10-12): the exit status of each script is primary; in addition the output is scanned for
# failure lines in any of the formats used in tests_cm ([FAIL], RESULT|..|FAIL, bare FAIL, SOME FAILED, DID NOT FAIL,
# "<n> FAIL"/"<n> failures" with n >= 1, EXIT STATUS 1), ignoring PASS lines and "fails as required".
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
def failure_lines(text):
    return [l for l in text.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]

def main():
    os.makedirs(os.path.join(HERE, 'outputs'), exist_ok=True)
    counts = {}; failed = []
    t0 = time.time()
    for s in SCRIPTS:
        t = time.time()
        r = subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=HERE, capture_output=True, text=True)
        out = r.stdout
        if r.returncode != 0:
            failed.append(f"{s} (exit {r.returncode})")
            out += f"\n[FAIL] {s} exited with status {r.returncode}:\n{r.stderr}"
        elif failure_lines(out):
            failed.append(f"{s} (failure lines)")
        with open(os.path.join(HERE, 'outputs', s.replace('.py', '.txt')), 'w') as f:
            f.write(out)
        print(f"# {s} ({time.time() - t:.0f}s)")
        for line in out.splitlines():
            print(line)
            m = re.match(r'\[(PASS|FAIL|SANITY|UNTESTED)\]', line)
            if m:
                counts[m.group(1)] = counts.get(m.group(1), 0) + 1
    print(f"\nSUMMARY ({time.time() - t0:.0f}s): " + ", ".join(f"{k} {counts.get(k, 0)}"
                                                         for k in ['PASS', 'FAIL', 'SANITY', 'UNTESTED']))
    if failed:
        print("FAILED SCRIPTS:", ", ".join(failed))
    return 1 if (counts.get('FAIL', 0) or failed) else 0

if __name__ == '__main__':
    sys.exit(main())
