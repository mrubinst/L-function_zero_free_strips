"""Run every numerical check of the 11/12 zeta guide and summarize the results.

Usage:  python3 run_all.py            (all scripts, about 2 minutes)
        python3 run_all.py --quick    (skip the 75-second theta reflection test)

Each script's output is written to outputs/<script name>.txt.  The summary counts the [PASS] and [FAIL]
lines; lines without a tag are sanity checks of bounds and are not counted.  The exit status is 1 if any
check fails or any script crashes, and 0 otherwise.
"""
import os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["test_sec1_2.py", "test_sec3.py", "test_sec4_5.py", "test_sec4_8_extra.py", "test_sec6.py",
           "test_sec7_8.py", "test_sec7_completed_sanity.py", "test_sec7_reflection_j1.py"]

def main():
    quick = "--quick" in sys.argv[1:]
    os.makedirs(os.path.join(HERE, "outputs"), exist_ok=True)
    total_pass = total_fail = 0; bad_scripts = []
    for name in SCRIPTS:
        if quick and name == "test_sec7_reflection_j1.py":
            print(f"{name:32s} skipped (--quick)")
            continue
        t0 = time.time()
        proc = subprocess.run([sys.executable, "-u", os.path.join(HERE, name)], cwd=HERE,
                              capture_output=True, text=True)
        out = proc.stdout + (("\n" + proc.stderr) if proc.stderr else "")
        with open(os.path.join(HERE, "outputs", name.replace(".py", ".txt")), "w") as fh:
            fh.write(out)
        npass = len(re.findall(r"^\[PASS\]", out, re.M)); nfail = len(re.findall(r"^\[FAIL\]", out, re.M))
        total_pass += npass; total_fail += nfail
        status = "ok" if proc.returncode == 0 and nfail == 0 else "PROBLEM"
        if status != "ok":
            bad_scripts.append(name)
        print(f"{name:32s} {npass:3d} pass  {nfail:2d} fail  exit {proc.returncode}  {time.time()-t0:6.1f}s  {status}")
    print(f"\nTotal: {total_pass} checks passed, {total_fail} failed.")
    if bad_scripts:
        print("Problems in:", ", ".join(bad_scripts))
    sys.exit(1 if bad_scripts else 0)

if __name__ == "__main__":
    main()
