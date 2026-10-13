"""Run the t5 display tests (Appendix A of the CM guide) and, with --reruns, the existing tests they rely on.

    python3 -B run_all.py            # t5a-t5d, about 1.5 minutes; outputs in outputs/<script>.txt
    python3 -B run_all.py --reruns   # also rerun the existing scripts (about 8 minutes wall clock, run in parallel);
                                     # outputs in outputs/reruns/ (nothing is written into the other suites)
Exit status 1 if any script exits nonzero.
"""
import os, sys, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, '..', '..'))
OWN = ["t5a_finite_identities.py", "t5b_theta_pointwise.py", "t5c_mellin.py", "t5d_decay.py"]
RERUNS = [("proofs_2026-10-12/b3", "b3_symbolic_checks.py"), ("proofs_2026-10-12/b3", "t_b3_prop_endtoend.py"),
          ("phase2_2026-10-12/r2", "t1_prop_grouped_endtoend.py"), ("phase2_2026-10-12/r2", "t2_grouping_cusp_data.py"),
          ("part_B_lines_0550_0997", "t4_theta.py"), ("part_B_lines_0550_0997", "t5_reflection_j2.py"),
          ("part_D_lines_1510_end", "t02_chain_rule.py"), ("part_D_lines_1510_end", "t03_mellin.py"),
          ("part_D_lines_1510_end", "t04_kernels.py"), ("part_D_lines_1510_end", "t07_theta_appendix.py"),
          ("part_D_lines_1510_end", "t08_reflection_j.py")]
os.makedirs(os.path.join(HERE, "outputs", "reruns"), exist_ok=True)
bad = []; T0 = time.time()
for s in OWN:
    t0 = time.time()
    p = subprocess.run([sys.executable, "-B", os.path.join(HERE, s)], capture_output=True, text=True, cwd=HERE)
    open(os.path.join(HERE, "outputs", s.replace(".py", ".txt")), "w").write(p.stdout + p.stderr)
    print(f"{s}: exit {p.returncode} ({time.time() - t0:.0f} s)", flush=True)
    if p.returncode: bad.append(s)
if "--reruns" in sys.argv:
    PREFIX = {"proofs_2026-10-12/b3": "", "phase2_2026-10-12/r2": "r2_", "part_B_lines_0550_0997": "partB_",
              "part_D_lines_1510_end": "partD_"}
    procs = []
    for d, s in RERUNS:
        out = open(os.path.join(HERE, "outputs", "reruns", PREFIX[d] + s.replace(".py", ".txt")), "w")
        procs.append((d, s, time.time(), subprocess.Popen([sys.executable, "-B", s], cwd=os.path.join(TESTS, d), stdout=out,
                                                          stderr=subprocess.STDOUT)))
    summ = []
    for d, s, t0, p in procs:
        rc = p.wait()
        summ.append(f"{d}/{s} exit={rc} wall={time.time() - t0:.0f}s (all started together)")
        print(summ[-1], flush=True)
        if rc: bad.append(f"{d}/{s}")
    open(os.path.join(HERE, "outputs", "reruns", "SUMMARY.txt"), "w").write("\n".join(summ) + "\n")
print(f"total wall time {time.time() - T0:.0f} s")
print("ALL PASS" if not bad else f"FAILED: {bad}")
sys.exit(1 if bad else 0)
