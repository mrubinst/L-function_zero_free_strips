"""Check that every display-test runner reports a failure of any script it runs, however deeply nested.

For each group t0-t9, and for the top-level run_all.sh, the runner is executed on a scratch copy of the test folder
(so no saved output is overwritten), with a sitecustomize.py that intercepts every Python script started:
  * control run: every script prints "[PASS] skipped" and exits 0 -- the runner must exit 0, and the scripts it
    starts are recorded;
  * one run per recorded script: that script prints "[FAIL] forced" and exits 1, all others are skipped -- the
    runner must exit nonzero.
Some runners (t3, t4) check the content of their scripts' output, so the skipped control fails. For these the real
run is made once with every script's output recorded; the control and the forced failures then replay it.
Runner scripts (run_all.py, run_t3.py, rerun_existing.py) run normally, so their own bookkeeping is exercised.
Scripts started with python3 -I, and PARI/GP scripts, cannot be intercepted; they run normally and are listed.

    python3 check_runner_propagation.py            # all groups and the top level, about 10 minutes
    python3 check_runner_propagation.py t1 t5      # selected groups
Exit status 1 if any forced failure is not reported, or any control run fails.
"""
import os, sys, shutil, subprocess, tempfile, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)                      # tests_cm
GROUPS = {"t0": ["./run_all.sh"], "t1": ["./run_all.sh"], "t2": ["./run_all.sh"], "t3": ["python3", "run_t3.py"],
          "t4": ["python3", "run_all.py"], "t5": ["python3", "run_all.py", "--reruns"], "t6": ["./run_all.sh"],
          "t7": ["./run_all.sh"], "t8": ["./run_all.sh"], "t9": ["./run_all.sh"]}
RUNNERS = {"run_all.py", "run_t3.py", "rerun_existing.py"}

SITECUSTOMIZE = r'''
import os, sys, hashlib
_mode = os.environ.get("CMRP_MODE")
if _mode and sys.argv and sys.argv[0] not in ("", "-c", "-m"):
    _p = os.path.realpath(sys.argv[0])
    if os.path.basename(_p) not in %r:
        _key = hashlib.sha1(repr((_p, sys.argv[1:])).encode()).hexdigest()
        _cache = os.path.join(os.environ["CMRP_CACHE"], _key)
        if _mode == "record":
            class _Tee:
                def __init__(self, a, b): self.a, self.b = a, b
                def write(self, x): self.b.write(x); self.b.flush(); return self.a.write(x)
                def flush(self): self.a.flush(); self.b.flush()
                def __getattr__(self, n): return getattr(self.a, n)
            _f = open(_cache, "a")
            sys.stdout = _Tee(sys.stdout, _f); sys.stderr = _Tee(sys.stderr, _f)
        else:
            with open(os.environ["CMRP_LOG"], "a") as _f:
                _f.write(_p + "\n")
            if _p == os.environ.get("CMRP_TARGET"):
                print("[FAIL] forced failure (runner propagation check)"); sys.stdout.flush(); os._exit(1)
            if _mode == "replay" and os.path.exists(_cache):
                sys.stdout.write(open(_cache).read())
            else:
                print("[PASS] skipped (runner propagation check)")
            sys.stdout.flush(); os._exit(0)
''' % (RUNNERS,)


def make_copy():
    """A scratch copy of the whole test folder (which is self-contained: its sources are in sources/)."""
    tmp = os.path.realpath(tempfile.mkdtemp(prefix="cm_runner_check_"))
    dst = os.path.join(tmp, os.path.basename(TESTS))
    shutil.copytree(TESTS, dst, symlinks=True)
    sc = os.path.join(tmp, "sitecustomize")
    os.makedirs(sc); os.makedirs(os.path.join(tmp, "cache"))
    open(os.path.join(sc, "sitecustomize.py"), "w").write(SITECUSTOMIZE)
    return tmp, dst, sc


def run(cwd, cmd, sc, log, target=None, mode="skip"):
    env = dict(os.environ, CMRP_MODE=mode, CMRP_LOG=log, CMRP_CACHE=os.path.join(os.path.dirname(sc), "cache"),
               PYTHONPATH=sc, PYTHONDONTWRITEBYTECODE="1")
    env.pop("CMRP_TARGET", None)
    if target:
        env["CMRP_TARGET"] = target
    open(log, "w").close()
    p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                       stdin=subprocess.DEVNULL)
    started = [l for l in open(log).read().splitlines() if l]
    return p.returncode, started, p.stdout


def main():
    sel = [a for a in sys.argv[1:] if a in GROUPS or a == "top"] or list(GROUPS) + ["top"]
    tmp, dst, sc = make_copy()
    disp = os.path.join(dst, "displays_2026-10-13")
    log = os.path.join(tmp, "started.txt")
    bad, report = [], {}
    for g in sel:
        cwd, cmd = (disp, ["./run_all.sh"]) if g == "top" else (os.path.join(disp, g), GROUPS[g])
        t0 = time.time()
        mode = "skip"
        if g == "top":
            # reuse the outputs recorded for the content-checking groups (same copy, same scripts and arguments)
            mode = "replay"
        rc, started, out = run(cwd, cmd, sc, log, mode=mode)
        if rc != 0 and g != "top":
            # this runner checks the content of its scripts' output: record the real output once, then replay it
            mode = "replay"
            rec_rc, _, rec_out = run(cwd, cmd, sc, log, mode="record")
            if rec_rc != 0:
                bad.append(f"{g}: the real (recording) run exited {rec_rc}")
                os.makedirs(os.path.join(HERE, "runner_check"), exist_ok=True)
                open(os.path.join(HERE, "runner_check", f"record_{g}.txt"), "w").write(rec_out)
            rc, started, out = run(cwd, cmd, sc, log, mode="replay")
        scripts = list(dict.fromkeys(started))
        ok_control = rc == 0
        if not ok_control:
            bad.append(f"{g}: control run exited {rc}")
            os.makedirs(os.path.join(HERE, "runner_check"), exist_ok=True)
            open(os.path.join(HERE, "runner_check", f"control_{g}.txt"), "w").write(out)
        targets = scripts if g != "top" else [s for i, s in enumerate(scripts) if i in (0, len(scripts) // 2, len(scripts) - 1)]
        missed = []
        for s in targets:
            rc2, _, out2 = run(cwd, cmd, sc, log, target=s, mode=mode)
            if rc2 == 0:
                missed.append(os.path.relpath(s, dst))
        bad += [f"{g}: forced failure of {m} not reported" for m in missed]
        report[g] = dict(mode=mode, control_exit=rc, scripts_started=len(scripts), forced=len(targets), missed=missed,
                         seconds=round(time.time() - t0))
        print(f"{g} ({mode}): control exit {rc}; {len(scripts)} Python scripts started; {len(targets)} forced failures, "
              f"{len(targets) - len(missed)} reported, {len(missed)} missed ({time.time() - t0:.0f} s)", flush=True)
        for m in missed:
            print(f"   [FAIL] not reported: {m}")
    json.dump(report, open(os.path.join(HERE, "runner_propagation.json"), "w"), indent=1)
    shutil.rmtree(tmp, ignore_errors=True)
    print("ALL FORCED FAILURES REPORTED" if not bad else "PROBLEMS:\n  " + "\n  ".join(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
