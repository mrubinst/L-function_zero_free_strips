"""Shared reporting for the new display tests of agent t1 (tests_cm/displays_2026-10-13/t1).

Each check prints one line
    [PASS] <display label> (<printed number>) <test id>: message | falsifier detected
    [FAIL] ...
A check passes only if the identity/inequality holds on the tested range AND, when a falsifier (a deliberately wrong
variant) is supplied, that variant FAILS.  The script exits with status 1 if any check failed (atexit hook).
CM_TESTS_FORCE_FAIL=1 makes the first check fail (self-test of the wiring).

Helpers of tests_cm/part_A_lines_0001_0549 (cmcore.py, zw.py) are imported read-only; no bytecode is written there.
"""
import os, sys, atexit
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
PART_A = os.path.normpath(os.path.join(HERE, '..', '..', 'part_A_lines_0001_0549'))
sys.path.insert(0, PART_A)

_FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']
_NF = [0]
_N = [0]


def _hook():
    print(f"SUMMARY: {_N[0]} checks, {_NF[0]} failed")
    sys.stdout.flush()
    if _NF[0]:
        print(f"EXIT STATUS 1: {_NF[0]} check(s) failed", file=sys.stderr)
        sys.stderr.flush()
        os._exit(1)


atexit.register(_hook)


def check(tid, label, number, ok, msg, falsifier=None):
    """falsifier: None (no falsifier) or a bool that is True when the wrong variant was detected (failed)."""
    if _FORCE[0]:
        _FORCE[0] = False
        ok = False
        msg = 'CM_TESTS_FORCE_FAIL=1 (forced) ' + msg
    good = bool(ok) and (falsifier is None or bool(falsifier))
    _N[0] += 1
    if not good:
        _NF[0] += 1
    fs = '' if falsifier is None else (' | falsifier detected' if falsifier else ' | FALSIFIER NOT DETECTED')
    print(f"[{'PASS' if good else 'FAIL'}] {label} ({number}) {tid}: {msg}{fs}", flush=True)
    return good


def sanity(tid, label, number, msg):
    print(f"[SANITY] {label} ({number}) {tid}: {msg}", flush=True)


# ------------------------------------------------------------------ q_d(m) from PARI a_p
def gp_ap(d, X):
    """{p: a_p(E_d)} for primes p <= X (PARI ellap on y^2 = x^3 + d)."""
    from cmcore import gp
    out = gp(f"E=ellinit([0,0,0,0,{d}]);forprime(p=2,{X},print(p,\" \",ellap(E,p)))", timeout=1200)
    return {int(a): int(b) for a, b in (l.split() for l in out)}


def qd_array(d, M, ap):
    """numpy array q[0..M] with q[m] = q_d(m) of eq:qd (q[0] = 0): multiplicative, q(p) = -a_p/sqrt p,
    q(p^2) = 1, q(p^j) = 0 for j >= 3, and q(m) = 0 if (m, 6d) > 1."""
    import math
    import numpy as np
    q = np.zeros(M + 1)
    spf = list(range(M + 1))
    for i in range(2, int(M**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, M + 1, i):
                if spf[j] == j:
                    spf[j] = i
    for m in range(1, M + 1):
        if math.gcd(m, 6*d) > 1:
            continue
        r = 1.0; mm = m
        while mm > 1:
            p = spf[mm]; e = 0
            while mm % p == 0:
                mm //= p; e += 1
            if e == 1:
                r *= -ap[p]/math.sqrt(p)
            elif e >= 3:
                r = 0.0; break
        q[m] = r
    return q


def gl_panels(t0, t1, h, k=16):
    """Gauss-Legendre nodes and weights on [t0, t1] split into panels of width about h, k nodes per panel."""
    import numpy as np
    x, w = np.polynomial.legendre.leggauss(k)
    n = max(1, int(np.ceil((t1 - t0)/h)))
    edges = np.linspace(t0, t1, n + 1)
    a = edges[:-1, None]; b = edges[1:, None]
    T = (a + b)/2 + (b - a)/2*x[None, :]
    Wt = (b - a)/2*w[None, :]
    return T.ravel(), Wt.ravel()
