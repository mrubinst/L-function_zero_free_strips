"""Shared helpers for the part-B checks (lines 550-997 of cm-eleven-twelfths-guide.tex).

Conventions are the paper's: omega = e^{2 pi i/3}; a + b omega is the integer pair (a, b); primary = 1 mod 3;
lambda = 1 + 2 omega = i sqrt 3; chi_n(x) = (x/n)_6 (sixth root of unity = x^{(N pi - 1)/6} mod pi at primes,
zero on multiples, multiplicative in n); e(z) = exp(4 pi i Im z / sqrt 3); ce(z) = exp(2 pi i (z + conj z));
alpha(n) = n/|n|; gamma_j(n) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n); G(n) = conj(chi_n(4)) gamma_3(n).

The scalar helpers (chi, mu, factor, primaries, gamma_j) come from common.py / zw.py (copied unchanged from
tests_guide).  This module adds numpy-vectorized symbols and Gauss sums (checked against brute force in t3),
and the report() function used by every script.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import *            # zw (mul, N, conj, divides, ...), chi, mu, factor, primaries, gamma_j, G, PR_ODD

SQ3 = math.sqrt(3.0)

# ------------------------------------------------------------------------------------------- reporting
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

def report(cid, err, tol, detail="", note=""):
    """One line per claim: 'RESULT|id|PASS/FAIL|detail'.  Passes iff err <= tol."""
    ok = (err <= tol)
    if _forced():
        ok = False; detail = "CM_TESTS_FORCE_FAIL=1 (forced failure); " + detail
    if not ok:
        _NFAILED[0] += 1
    print(f"RESULT|{cid}|{'PASS' if ok else 'FAIL'}|err={err:.3g} tol={tol:g}; {detail}" + (f" NOTE: {note}" if note else ""),
          flush=True)
    return ok

def falsifier(cid, name, err, must_exceed):
    """A deliberately wrong variant: it must FAIL, i.e. its error must exceed must_exceed.  Returns True if it failed
    as required (so the falsifier 'works')."""
    worked = err > must_exceed
    if _forced():
        worked = False; name = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + name
    if not worked:
        _NFAILED[0] += 1
    print(f"      falsifier for {cid} [{name}]: err={err:.3g} -> {'fails as required' if worked else 'DID NOT FAIL (bad falsifier)'}",
          flush=True)
    return worked

# ------------------------------------------------------------------------------------------- vectorized symbols
UNIT_PAIRS = UNITS                  # zeta^j, zeta = 1 + omega, j = 0..5 (from zw.py)
ZETA6 = np.exp(1j*np.pi*np.arange(6)/3)

def _split_data(pi):
    a, b = pi; p = N(pi)
    w0 = (-a*pow(b, -1, p)) % p     # image of omega in F_p
    units = [ (u[0] + u[1]*w0) % p for u in UNIT_PAIRS ]
    return p, w0, units

def _powmod_vec(x, e, p):
    r = np.ones_like(x); b = x % p
    while e:
        if e & 1: r = (r*b) % p
        b = (b*b) % p; e >>= 1
    return r

def sym6_index_vec(A, B, pi):
    """Index j (0..5) of chi_pi(A + B omega) = zeta^j, or -1 if pi | x; A, B int64 arrays; pi a primary prime != -2."""
    A = np.asarray(A, dtype=np.int64); B = np.asarray(B, dtype=np.int64)
    out = np.full(A.shape, -1, dtype=np.int64)
    if pi[1] != 0:                                   # split prime
        p, w0, units = _split_data(pi)
        x = (A % p + (B % p)*w0) % p
        r = _powmod_vec(x, (p - 1)//6, p)
        for j, u in enumerate(units):
            out[(r == u) & (x != 0)] = j
        assert np.all((out >= 0) | (x == 0))
        return out
    q = -pi[0]; assert q > 0                          # inert prime -q, residue field F_q[omega]
    a = A % q; b = B % q
    e = (q*q - 1)//6
    ra = np.ones_like(a); rb = np.zeros_like(a); ba = a.copy(); bb = b.copy()
    while e:
        if e & 1:
            ra, rb = (ra*ba - rb*bb) % q, (ra*bb + rb*ba - rb*bb) % q
        ba, bb = (ba*ba - bb*bb) % q, (2*ba*bb - bb*bb) % q
        e >>= 1
    nz = (a != 0) | (b != 0)
    for j, u in enumerate(UNIT_PAIRS):
        out[(ra == u[0] % q) & (rb == u[1] % q) & nz] = j
    assert np.all((out >= 0) | (~nz))
    return out

def chi_index_vec(n, A, B):
    """Index of chi_n(A + B omega) for primary n coprime to 6 (any n, product over its prime factors); -1 for 0."""
    A = np.asarray(A, dtype=np.int64); B = np.asarray(B, dtype=np.int64)
    tot = np.zeros(A.shape, dtype=np.int64); zero = np.zeros(A.shape, dtype=bool)
    for p in factor(n):
        j = sym6_index_vec(A, B, p)
        zero |= (j < 0); tot += np.where(j < 0, 0, j)
    return np.where(zero, -1, tot % 6)

def chi_vec(n, A, B, power=1):
    """chi_n(x)^power as a complex array (zero where chi_n(x)=0, for every power, as in the paper)."""
    j = chi_index_vec(n, A, B)
    v = ZETA6[(power*np.where(j < 0, 0, j)) % 6]
    return np.where(j < 0, 0, v)

def hnf_reps(n):
    """A complete residue system of O/n as arrays (A, B): 0 <= A < h1, 0 <= B < h2 with h1 h2 = N(n)."""
    nn = N(n)
    w = mul(n, (0, 1))
    h2 = math.gcd(n[1], w[1]); h1 = nn//h2
    assert h1*h2 == nn
    A, B = np.meshgrid(np.arange(h1, dtype=np.int64), np.arange(h2, dtype=np.int64), indexing='ij')
    return A.ravel(), B.ravel()

def e_phase_vec(A, B, n, mult=(1, 0)):
    """e(mult * x / n) for x = A + B omega: x mult conj(n) = P + Q omega and e = exp(2 pi i Q / N(n))."""
    nn = N(n); c, d = mul(mult, conj(n))
    Q = (A*d + B*c - B*d) % nn
    return np.exp(2j*np.pi*Q/nn)

def ce_phase_vec(A, B, n, mult=(1, 0)):
    """ce(mult * x / n) = exp(2 pi i Tr(x mult conj n)/N(n))."""
    nn = N(n); c, d = mul(mult, conj(n))
    P = A*c - B*d; Q = A*d + B*c - B*d
    return np.exp(2j*np.pi*((2*P - Q) % nn)/nn)

def gamma_j_fast(j, n):
    """gamma_j(n) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n) (vectorized; checked against brute force)."""
    if n == (1, 0): return 1.0 + 0j
    A, B = hnf_reps(n)
    return complex(np.sum(chi_vec(n, A, B, j)*e_phase_vec(A, B, n))/math.sqrt(N(n)))

def gauss_psi_fast(z1, z2):
    """gamma(conj chi_{z1} chi_{z2}) modulo q = z1 z2 (vectorized)."""
    q = mul(z1, z2)
    if q == (1, 0): return 1.0 + 0j
    A, B = hnf_reps(q)
    v = np.conj(chi_vec(z1, A, B))*chi_vec(z2, A, B)
    return complex(np.sum(v*e_phase_vec(A, B, q))/math.sqrt(N(q)))

def gauss_psi_general(z1, z2, conj_first=True):
    """gamma of conj(chi_{z1})chi_{z2} (conj_first) or chi_{z1}conj(chi_{z2})."""
    q = mul(z1, z2)
    if q == (1, 0): return 1.0 + 0j
    A, B = hnf_reps(q)
    if conj_first: v = np.conj(chi_vec(z1, A, B))*chi_vec(z2, A, B)
    else: v = chi_vec(z1, A, B)*np.conj(chi_vec(z2, A, B))
    return complex(np.sum(v*e_phase_vec(A, B, q))/math.sqrt(N(q)))

def G_fast(n):
    return complex(chi(n, (4, 0))).conjugate()*gamma_j_fast(3, n)

def cubic3(n, x):
    """(x/n)_3 computed independently by Euler's criterion with exponent (N(pi)-1)/3 at each prime."""
    from fastgauss import cubic_prime
    r = 1 + 0j
    for p in factor(n):
        r *= cubic_prime(x, p)
    return r

def lattice(Rn):
    """All (a, b) with N(a + b omega) <= Rn (including 0)."""
    R = int(2*math.sqrt(Rn)) + 2
    return [(a, b) for a in range(-R, R+1) for b in range(-R, R+1) if N((a, b)) <= Rn]

def gcd_primary(a, b):
    fa = set(factor(a)); fb = set(factor(b)); g = (1, 0)
    for p in fa & fb: g = mul(g, p)
    return g

def divisors_sqfree(g):
    ps = sorted(set(factor(g))); out = []
    for mask in range(1 << len(ps)):
        d = (1, 0)
        for i, p in enumerate(ps):
            if mask >> i & 1: d = mul(d, p)
        out.append(d)
    return out

def alpha_c(n):
    z = cx(n); return z/abs(z)

def bump(y, a=1.0, b=2.0):
    return math.exp(-1/((y - a)*(b - y))) if a < y < b else 0.0
