"""Run every part-B script and print one line per claim of lines 550-997 of cm-eleven-twelfths-guide.tex.

    python3 run_all.py            # all scripts, about 3-4 minutes; outputs saved to outputs/<script>.txt
    python3 run_all.py --quick    # skip t5_reflection_j2.py (the 2-minute end-to-end reflection test)

Each claim line is  [STATUS] id  lines  (label)  [class]  claim  ::  test detail
STATUS is PASS or FAIL for tested claims, UNTESTED for cited claims that are not tested (reason given), SKIPPED in
--quick mode.  A claim also FAILs if one of its falsifiers (deliberately wrong variants) did not fail.
B88 (the table of zeta -> CM changes) is assembled from B17, B27, B46, B50, B53, B55, B68 and their zeta-form controls.
Exit status 1 if any claim fails or any script exits nonzero (every script exits 1 when one of its checks fails).
"""
import os, sys, subprocess, time, re
HERE = os.path.dirname(os.path.abspath(__file__))
# Secondary check (added 2026-10-12): the exit status of each script is primary; in addition the output is scanned for
# failure lines in any of the formats used in tests_cm ([FAIL], RESULT|..|FAIL, bare FAIL, SOME FAILED, DID NOT FAIL,
# "<n> FAIL"/"<n> failures" with n >= 1, EXIT STATUS 1), ignoring PASS lines and "fails as required".
FAIL_RE = re.compile(r"^\s*\[(FAIL|CITED-X)\]|^RESULT\|[^|]*\|FAIL|^\s*FAIL([^A-Za-z_]|$)|SOME FAIL|DID NOT FAIL"
                     r"|FALSIFIER NOT DETECTED|EXIT STATUS [1-9]|<-- PROBLEM|(^|[^0-9])[1-9][0-9]* (FAIL|failures)([^A-Za-z_]|$)")
OK_RE = re.compile(r"^\s*(\[(PASS|INFO|SANITY|UNTESTED|CITED)\]|RESULT\|[^|]*\|PASS\||PASS([^A-Za-z_]|$))|fails as required")
def failure_lines(text):
    return [l for l in text.splitlines() if FAIL_RE.search(l) and not OK_RE.search(l)]
sys.path.insert(0, HERE)
from claims import CLAIMS, WHY_UNTESTED

SCRIPTS = ["t1_zeros.py", "t2_poisson.py", "t3_gauss.py", "t4_theta.py", "t5_reflection_j2.py", "t6_sieve.py", "t7_descent.py"]
quick = "--quick" in sys.argv
os.makedirs(os.path.join(HERE, "outputs"), exist_ok=True)
results = {}; falsif = {}; crashed = []; ran = set()
for s in SCRIPTS:
    if quick and s == "t5_reflection_j2.py":
        continue
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, s)], capture_output=True, text=True, cwd=HERE)
    out = p.stdout + (("\n--- stderr ---\n" + p.stderr) if p.stderr.strip() else "")
    open(os.path.join(HERE, "outputs", s.replace(".py", ".txt")), "w").write(out)
    ran.add(s)
    if p.returncode != 0:
        crashed.append(f"{s} (exit {p.returncode})")
    elif failure_lines(p.stdout):
        crashed.append(f"{s} (failure lines)")
    for line in p.stdout.splitlines():
        if line.startswith("RESULT|"):
            _, cid, st, det = line.split("|", 3)
            results[cid] = (st, det)
        m = re.match(r"\s*falsifier for (B\d+) \[(.*)\]: err=(\S+) -> (.*)", line)
        if m:
            falsif.setdefault(m.group(1), []).append((m.group(2), m.group(4).startswith("fails as required")))
    print(f"# {s}: {('exit status %d' % p.returncode) if p.returncode else 'ok'} ({time.time() - t0:.0f} s)", flush=True)

# B88: the table of changes, from its constituent claims and their zeta-form controls
parts = ["B17", "B27", "B46", "B50", "B53", "B55", "B68"]
if all(c in results for c in parts):
    ok = all(results[c][0] == "PASS" for c in parts) and all(all(w for _, w in falsif.get(c, [])) for c in parts)
    results["B88"] = ("PASS" if ok else "FAIL",
                      "rows: mu alpha (B17 end to end, control alpha(n1)conj(alpha(n2)) fails); conj(alpha)^2 gamma_2 (B27, zeta "
                      "coefficients fail); d_zbar^2 (B50, B53: one-derivative kernel fails); conj(alpha(b))^6 (B46, B68: conj(alpha(h))^3 "
                      "fails); alpha(l)^2 (B55: alpha(l)^1 fails)")

counts = {}; nfail = 0
print()
for cid, lines, label, cls, claim, script in CLAIMS:
    if script is None:
        st, det = "UNTESTED", WHY_UNTESTED.get(cid, "")
    elif cid in results:
        st, det = results[cid]
        bad = [n for n, w in falsif.get(cid, []) if not w]
        if bad:
            st = "FAIL"; det += f"  [falsifier did not fail: {bad}]"
        nf = len(falsif.get(cid, []))
        if nf: det += f"  [{nf} falsifier(s) failed as required]"
    elif quick and script == "t5_reflection_j2.py" or (quick and cid == "B88"):
        st, det = "SKIPPED", "--quick"
    else:
        st, det = "FAIL", f"no result reported by {script}"
    nfail += st == "FAIL"
    counts[(cls, st)] = counts.get((cls, st), 0) + 1
    print(f"[{st}] {cid}  L{lines}  ({label})  [{cls}]  {claim}  ::  {det}")
print()
classes = sorted(set(c for c, _ in counts))
for c in classes:
    print(f"class {c}: " + ", ".join(f"{st} {n}" for (cc, st), n in sorted(counts.items()) if cc == c))
tot = {}
for (c, st), n in counts.items(): tot[st] = tot.get(st, 0) + n
print("total: " + ", ".join(f"{st} {n}" for st, n in sorted(tot.items())) + f"  ({len(CLAIMS)} claims)")
if crashed: print("FAILED SCRIPTS (nonzero exit status or failure lines):", crashed)
sys.exit(1 if (nfail or crashed) else 0)
