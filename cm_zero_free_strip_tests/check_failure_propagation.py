"""Check that failed checks and nested suites produce a nonzero exit status.

Run with an interpreter containing the math-suite dependencies. Temporary runner
copies and fake commands isolate these checks from the mathematical test logs.
"""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
ENV = os.environ.copy()
ENV.pop("CM_TESTS_FORCE_FAIL", None)
N_CHECKS = 0


def expect(command, expected, *, env, cwd=None):
    global N_CHECKS
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)
    if result.returncode != expected:
        raise AssertionError(
            f"{command}: expected exit {expected}, got {result.returncode}\n"
            + result.stdout + result.stderr
        )
    N_CHECKS += 1


def check_helpers():
    helpers = [
        (ROOT / "part_B_lines_0550_0997", "cmlib",
         "m.report('exit-check', ERROR, 0)",
         "m.falsifier('exit-check', 'ineffective control', 0, 1)"),
        (ROOT / "phase2_2026-10-12/s3", "s3lib",
         "m.report('exit-check', 'self-test', 'check', not ERROR)", None),
        (ROOT / "phase2_2026-10-12/s5", "s5lib",
         "m.report('exit-check', ERROR, 0)",
         "m.falsifier('exit-check', 'ineffective control', 0, 1)"),
    ]
    for directory, module, report, falsifier in helpers:
        prefix = "import sys; sys.path.insert(0, sys.argv[1]); import " + module + " as m; "
        for error in (0, 1):
            code = prefix + report.replace("ERROR", str(error))
            expect([sys.executable, "-B", "-c", code, str(directory)], error, env=ENV)
        forced = ENV | {"CM_TESTS_FORCE_FAIL": "1"}
        code = prefix + report.replace("ERROR", "0")
        expect([sys.executable, "-B", "-c", code, str(directory)], 1, env=forced)
        if falsifier:
            expect([sys.executable, "-B", "-c", prefix + falsifier, str(directory)], 1, env=ENV)


def check_runners():
    runners = [
        ("proofs_2026-10-12/run_all.sh", "t_points.py"),
        ("phase2_2026-10-12/run_all.sh", "t1_labels.py"),
        ("phase2_2026-10-12/s3/run_all.sh", "t1_labels.py"),
        ("phase2_2026-10-12/s5/run_all.sh", "t1_dual_profile.py"),
        ("phase2_2026-10-12/s6/run_all.sh", "s6_twisted_paired.py"),
    ]
    with tempfile.TemporaryDirectory(prefix="cm-exit-check-") as temporary:
        scratch = Path(temporary)
        for relative, _ in runners:
            target = scratch / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, target)
        for folder in ("r1", "r2", "r3", "s1", "s2", "s3", "s4", "s5", "s6"):
            (scratch / "phase2_2026-10-12" / folder).mkdir(exist_ok=True)
        fake_bin = scratch / "bin"
        fake_bin.mkdir()
        fake = fake_bin / "python"
        fake.write_text(
            '#!/bin/sh\n'
            'for arg do\n'
            '  if [ "${arg##*/}" = "$CM_FAIL_SCRIPT" ]; then\n'
            '    echo "simulated unsuccessful check"; exit 2\n'
            '  fi\n'
            'done\n'
            'echo "0 FAIL; 0 failures"\n'
            'exit 0\n'
        )
        fake.chmod(0o755)
        shutil.copy2(fake, fake_bin / "gp")
        env = ENV | {"PYTHON": str(fake), "PATH": str(fake_bin) + os.pathsep + ENV["PATH"]}
        for relative, failed_script in runners:
            command = ["bash", str(scratch / relative)]
            # A zero-failure summary must not be a false positive.
            expect(command, 0, env=env | {"CM_FAIL_SCRIPT": "never.py"})
            # One early failure followed by successful tests must survive.
            expect(command, 1, env=env | {"CM_FAIL_SCRIPT": failed_script})


if __name__ == "__main__":
    check_helpers()
    check_runners()
    print(f"PASS: {N_CHECKS} failure-propagation checks")
