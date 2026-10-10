"""Shared helpers for the numerical checks of the 11/12 zeta guide (zeta_11_12_guide.pdf).

Conventions are the guide's: primary = 1 mod 3, lambda = 1 + 2 omega = i sqrt(3), chi_n(x) = (x/n)_6,
e(z) = exp(4 pi i Im z / sqrt(3)).  Elements of Z[omega] are pairs (a, b) = a + b omega (see zw.py).

Every test prints one line
    [PASS] <test id> [guide <display>] <description>: max error <err> (tol <tol>) <details>
where <display> is the equation number (or the text near it) in the guide that the test checks.
The map from test ids to guide displays is the dictionary GUIDE below.  Lines without [PASS]/[FAIL]
are sanity checks of the size of a bound; they print numbers but have no pass/fail criterion.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import *
from functools import lru_cache
import numpy as np
from fastgauss import primary_primes

PR = primary_primes(200000)                 # primary primes of norm <= 200000, including -2
PR_ODD = [p for p in PR if p != (-2, 0)]    # ... coprime to 6

@lru_cache(maxsize=None)
def factor(n):
    """Factor a primary n coprime to 6 into primary primes (tuple, with multiplicity), N(n) <= 200000^2."""
    out = []; m = n
    for p in PR_ODD:
        if N(m) == 1:
            break
        if N(p) > N(m):
            break
        while divides(p, m):
            out.append(p); m = div(m, p)
    assert N(m) == 1, (n, m)
    return tuple(out)

@lru_cache(maxsize=None)
def sym_index(x, p):
    """(x/p)_6 as the exponent j of zeta^j (zeta = e^{i pi/3}), or None if p | x; x is reduced mod p first."""
    return sym6_prime(canon_mod(x, p), p)

ZETA = [cmath.exp(1j*math.pi*j/3) for j in range(6)]

def chi(n, x):
    """chi_n(x) = (x/n)_6 (guide (3.1)), n primary coprime to 6; 0 if x and n are not coprime."""
    r = 1 + 0j
    for p in factor(n):
        j = sym_index(x, p)
        if j is None:
            return 0
        r *= ZETA[j]
    return r

def mu(n):
    """Moebius function mu_K(n) of a primary n coprime to 6 (guide Section 2)."""
    f = factor(n)
    return 0 if len(set(f)) < len(f) else (-1)**len(f)

def odd(x):
    """True if x is coprime to 2."""
    return not divides((2, 0), x)

def primaries(B, sqfree=False):
    """All primary n coprime to 6 with 1 <= N(n) <= B (only squarefree n if sqfree), sorted by norm."""
    out = []
    R = int(2*math.sqrt(B)) + 3
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            n = (a, b)
            if primary(n) and 0 < N(n) <= B and odd(n):
                if (not sqfree) or mu(n) != 0:
                    out.append(n)
    return sorted(out, key=lambda x: (N(x), x))

def gamma_j(j, n):
    """gamma_j(n) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n) (guide (6.6)), by brute force."""
    nn = cx(n); s = 0
    for x in residues(n):
        c = chi(n, x)
        if c != 0:
            s += c**j*e_add(cx(x)/nn)
    return s/math.sqrt(N(n))

def G(n):
    """G(n) = conj(chi_n(4)) gamma_3(n) (guide (6.6))."""
    return chi(n, (4, 0)).conjugate()*gamma_j(3, n)

# Test id -> the display of the guide it checks (equation numbers of zeta_11_12_guide.pdf).
GUIDE = {
    # Section 1 (test_sec1_2.py)
    '1.1': '(1.1)', '1.2a': '(1.2)', '1.2b': '(1.2)', '1.3a': '(1.3)', '1.3b': 'text after (1.3)',
    '1.3c': 'text after (1.3)', '1.5a': '(1.5)', '1.5b': 'text after (1.5)', '1.5c': '(1.5)',
    # Section 2 (test_sec1_2.py)
    '2.0': 'Sec. 2, omega', '2.1': '(2.1)', '2.1b': 'text after (2.1)', '2.2': '(2.2) and text after',
    '2.3a': 'text before (2.3)', '2.3b': 'helper for (2.3)', '2.3c': 'text before (2.3)', '2.3': '(2.3)',
    '2.3d': 'text after (2.3)', '2.4': '(1.4) = (2.4)',
    # Section 3 (test_sec3.py)
    '3.0a': 'text before (3.1)', '3.0b': 'text after (3.1)', '3.0c': 'text after (3.1)',
    '3.1': '(3.4) and text after', '3.2a': '(3.4)', '3.2b': '(3.5), (3.6)', '3.2c': '(3.7)',
    # Section 4 (test_sec4_5.py, test_sec4_8_extra.py)
    '4.1': '(4.1), (4.2)', '4.2': '(4.3)', '4.3': '(4.5)', '4.4': 'text after (4.5)', '4.5': '(4.7)',
    '4.6': '(4.8)', '4.7a': '(2.3), (3.7)', '4.7b': '(4.9)', '4.8': '(4.10)', '4.10': '(4.6)', '4.11': '(4.3)',
    # Section 5 (test_sec4_5.py)
    '5.1': '(5.1), (5.2)', '5.2': '(1.5), (5.4)', '5.3': '(5.5) and text after',
    # Section 6 (test_sec6.py, test_sec4_8_extra.py)
    '6.1': '(6.2)', '6.2': 'text after (6.2)', '6.3a': '(6.1)', '6.3b': '(6.4)', '6.5': '(6.7)',
    '6.6': 'text after (6.7)', '6.7': '(6.8)', '6.8': 'text before (6.10)', '6.9': '(6.10)',
    '6.10': '(3.5) + (6.4)', '6.11': '(6.5)', '6.12': '(7.10)', '6.13': '(7.11)',
    '6.14': 'Sec. 7.3, quadratic twists',
    # Section 7 (test_sec7_8.py, test_sec4_8_extra.py)
    '7.1': '(7.1)', '7.2': 'text after (7.1)', '7.3': 'text before (7.1)', '7.4a': '(7.3), (7.4)',
    '7.4b': '(7.3), (7.5)', '7.5': 'theta automorphy (input)', '7.9': 'text after (7.6)',
    '7.10': 'text after (7.8)', '7.11': '(7.6)', '7.12a': '(7.4)', '7.12b': '(7.5)',
    # Section 8 (test_sec7_8.py)
    '8.0': '(7.1), helper', '8.1': '(8.1)', '8.2': '(8.2)-(8.11)', '8.3': '(6.11), (8.13)', '8.4': '(8.7)',
}

def report(name, err, tol, extra=""):
    """Print one PASS/FAIL line.  name starts with the test id; err is the largest error found
    (or a count of failures), and the test passes if err <= tol.  Returns True on PASS."""
    ok = err <= tol
    tid = name.split()[0]
    ref = f" [guide {GUIDE[tid]}]" if tid in GUIDE else ""
    print(f"[{'PASS' if ok else 'FAIL'}] {tid}{ref}{name[len(tid):]}: max error {err:.3g} (tol {tol:g}) {extra}")
    return ok
