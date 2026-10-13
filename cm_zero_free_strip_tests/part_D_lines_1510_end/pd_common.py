"""Shared helpers for the Part D tests (lines 1510-end of cm-eleven-twelfths-guide.tex).

Conventions are the paper's: Z[omega] elements are pairs (a, b) = a + b*omega, primary = 1 mod 3,
lambda = 1 + 2 omega = i sqrt 3, chi_n(x) = (x/n)_6, alpha(n) = n/|n|, e(z) = exp(4 pi i Im z/sqrt 3),
ce(z) = exp(2 pi i (z + conj z)).  zw.py and fastgauss.py are unmodified copies of the helpers in
tests_guide/.

Every test prints exactly one line per claim:
    [PASS] <id> L<line> <label>: <description> | <details>
or [FAIL] ...; a falsifier (a deliberately wrong variant) is part of each test, and the test passes only
if the true statement holds AND the falsifier fails.  Lines starting with [INFO] carry no verdict.
"""
import os, sys, math, cmath, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import *          # noqa: F401,F403  (ring arithmetic, sextic symbols)

RESULTS = []

# Exit status (added 2026-10-12): a script that reports a failed check (or a falsifier that did not fail) through
# the helpers below exits with status 1 (atexit hook).  Self-test: with CM_TESTS_FORCE_FAIL=1 in the environment the
# first check reported here counts as failed.
import atexit
_FORCE_FAIL = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]
_NFAILED = [0]

def _forced():
    """True exactly once (for the first check) when CM_TESTS_FORCE_FAIL=1."""
    if _FORCE_FAIL[0]:
        _FORCE_FAIL[0] = False
        return True
    return False

def _exit_status_hook():
    if _NFAILED[0]:
        try:
            print(f"EXIT STATUS 1: {_NFAILED[0]} check(s) failed", file=sys.stderr)
            sys.stdout.flush(); sys.stderr.flush()
        finally:
            os._exit(1)

atexit.register(_exit_status_hook)

def report(tid, line, label, desc, ok, details=""):
    """Print one PASS/FAIL line and record it."""
    if _forced():
        ok = False; details = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + details
    if not ok:
        _NFAILED[0] += 1
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {tid} L{line} {label}: {desc} | {details}", flush=True)
    RESULTS.append((tid, ok))
    return ok

def info(msg):
    print(f"[INFO] {msg}", flush=True)

def gp(script, timeout=600):
    """Run a PARI/GP script (string) and return stdout."""
    r = subprocess.run(["gp", "-q", "-D", "parisizemax=2000000000"], input=script, capture_output=True,
                       text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return r.stdout

def eis_elements_upto(B):
    """All nonzero (a, b) with N(a + b omega) <= B, as a dict norm -> list of pairs."""
    out = {}
    R = int(2*math.sqrt(B/3)) + 2
    for a in range(-R - 2, R + 3):
        for b in range(-R - 2, R + 3):
            n = a*a - a*b + b*b
            if 0 < n <= B:
                out.setdefault(n, []).append((a, b))
    return out

def rational_primes(B):
    s = bytearray([1])*(B + 1); s[0] = s[1] = 0
    for i in range(2, int(B**.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(B + 1) if s[i]]

def split_primary_pair(p):
    """For p = 1 mod 3, the two primary primes pi, conj(pi) of norm p."""
    b = 1
    while 3*b*b <= 4*p:
        r = 4*p - 3*b*b; s = math.isqrt(r)
        if s*s == r and (s + b) % 2 == 0:
            a = (s + b)//2
            pi = prim_assoc((a, b))
            return pi, conj(pi)
        b += 1
    raise ValueError(p)

def mod2_class(g):
    """For g coprime to 2, the cube root of unity (as a pair) congruent to g mod 2."""
    for u in [(1, 0), (0, 1), (-1, -1)]:
        if divides((2, 0), sub(g, u)):
            return u
    raise ValueError(g)

def cubic_sym_prime(x, pi):
    """(x/pi)_3 for a prime pi coprime to 3 (including -2, norm 4), as a complex number; 0 if pi | x."""
    if divides(pi, x):
        return 0
    r = powmod(x, (N(pi) - 1)//3, pi)
    for u, val in [((1, 0), 1), ((0, 1), W), ((-1, -1), W*W)]:
        if divides(pi, sub(r, u)):
            return val
    raise ValueError((x, pi))
