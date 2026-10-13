"""Core helpers for the Part A audit (lines 1-549 of cm-eleven-twelfths-guide.tex).

Conventions are the paper's:
  * omega = e^{2 pi i/3}; an Eisenstein integer A + B*omega is the pair (A, B)   (zw.py, copied from tests_guide);
  * primary means = 1 mod 3, i.e. A = 1, B = 0 mod 3;
  * chi_pi(u) = (u/pi)_6 is the sixth root of unity = u^{(N pi - 1)/6} mod pi, 0 if pi | u, extended
    multiplicatively in the denominator; sixth roots of unity are stored as an index j (zeta^j, zeta = 1 + omega
    = e^{i pi/3}), with None (or -1 in numpy arrays) for the value 0;
  * mu_K is the Moebius function on primary elements; alpha(n) = n/|n|.

All results are printed by report(): one line per claim, '[PASS]' or '[FAIL]'.  For an exact identity or explicit
inequality the test passes only if (i) the claim holds on the whole tested range and (ii) a deliberately wrong
variant (the falsifier) FAILS on the same range.
"""
import os, sys, math, cmath, subprocess, functools, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from zw import (mul, add, sub, neg, conj, N, cx, divides, div, pw, primary, UNITS, LAM,
                canon_mod, powmod, sym6_prime)
import numpy as np

ZETA = [cmath.exp(1j*math.pi*j/3) for j in range(6)]
SQ3 = math.sqrt(3.0)

# ------------------------------------------------------------------ reporting
RESULTS = []

# Exit status (added 2026-10-12): a script that reports a FAIL through report() exits with status 1 (atexit hook
# below).  Self-test: with CM_TESTS_FORCE_FAIL=1 in the environment the first check reported here counts as failed.
import atexit
_FORCE_FAIL = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']
_NFAILED = [0]

def _exit_status_hook():
    if _NFAILED[0]:
        try:
            print(f"EXIT STATUS 1: {_NFAILED[0]} check(s) failed", file=sys.stderr)
            sys.stdout.flush(); sys.stderr.flush()
        finally:
            os._exit(1)

atexit.register(_exit_status_hook)

def report(cid, line, label, ok, msg, falsifier=None, kind='E'):
    """Print one result line.  falsifier: None if no falsifier applies, else True when the wrong variant was
    detected (i.e. it failed).  A claim PASSes only if ok and (falsifier is None or falsifier)."""
    if _FORCE_FAIL[0]:
        _FORCE_FAIL[0] = False; ok = False; msg = 'CM_TESTS_FORCE_FAIL=1 (forced failure) ' + msg
    status = 'PASS' if ok and (falsifier is None or falsifier) else 'FAIL'
    if status == 'FAIL':
        _NFAILED[0] += 1
    lab = f", {label}" if label else ""
    fs = "" if falsifier is None else (" | falsifier detected" if falsifier else " | FALSIFIER NOT DETECTED")
    print(f"[{status}] {cid} (l.{line}{lab}) [{kind}] {msg}{fs}", flush=True)
    RESULTS.append((cid, status))
    return status == 'PASS'

def info(cid, line, label, msg, kind='A'):
    """A sanity line for an analytic (A) claim: prints numbers, no pass/fail criterion."""
    lab = f", {label}" if label else ""
    print(f"[SANITY] {cid} (l.{line}{lab}) [{kind}] {msg}", flush=True)
    RESULTS.append((cid, 'SANITY'))

def untested(cid, line, label, msg, kind='A'):
    lab = f", {label}" if label else ""
    print(f"[UNTESTED] {cid} (l.{line}{lab}) [{kind}] {msg}", flush=True)
    RESULTS.append((cid, 'UNTESTED'))

# ------------------------------------------------------------------ rational integers
SPF_LIMIT = 2_000_000
_spf = None

def spf_table():
    global _spf
    if _spf is None:
        s = np.zeros(SPF_LIMIT + 1, dtype=np.int64)
        s[1] = 1
        for i in range(2, SPF_LIMIT + 1):
            if s[i] == 0:
                s[i] = i
                if i*i <= SPF_LIMIT:
                    blk = s[i*i::i]
                    blk[blk == 0] = i
        _spf = s
    return _spf

def factorint(m):
    """Factorisation of a positive integer m as a dict {p: e}."""
    m = abs(int(m)); out = {}
    if m <= SPF_LIMIT:
        s = spf_table()
        while m > 1:
            p = int(s[m]); out[p] = out.get(p, 0) + 1; m //= p
        return out
    p = 2
    while p*p <= m:
        while m % p == 0:
            out[p] = out.get(p, 0) + 1; m //= p
        p += 1
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out

def is_prime(n):
    return n > 1 and factorint(n) == {n: 1}

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]

def mobius(m):
    f = factorint(m)
    return 0 if any(e > 1 for e in f.values()) else (-1)**len(f)

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1)//2, p) == 1 else -1

def jacobi(a, m):
    """Jacobi symbol (a/m), m odd positive."""
    assert m > 0 and m % 2 == 1
    r = 1
    for p, e in factorint(m).items():
        r *= legendre(a, p)**e
    return r

def chi_m3(t):
    """chi_{-3} as defined at l.226-228."""
    return 0 if t % 3 == 0 else (1 if t % 3 == 1 else -1)

# ------------------------------------------------------------------ Eisenstein primes
def prim_assoc(x):
    c = [mul(u, x) for u in UNITS if primary(mul(u, x))]
    assert len(c) == 1, x
    return c[0]

@functools.lru_cache(maxsize=None)
def primes_above(p):
    """Primary primes above the rational prime p (p != 3): two conjugates if p = 1 mod 3, [-p] otherwise."""
    if p == 3:
        return ()
    if p % 3 == 2:
        return ((-p, 0),)
    b = 1
    while 3*b*b <= 4*p:
        r = 4*p - 3*b*b; s = math.isqrt(r)
        if s*s == r and (s + b) % 2 == 0:
            pi = prim_assoc(((s + b)//2, b))
            assert N(pi) == p
            return tuple(sorted({pi, conj(pi)}))
        b += 1
    raise ValueError(p)

@functools.lru_cache(maxsize=None)
def factor_primary(n):
    """Primary prime factorisation (with multiplicity) of a primary n coprime to 3."""
    out = []; m = n
    for p, e in sorted(factorint(N(n)).items()):
        for pi in primes_above(p):
            while N(m) > 1 and divides(pi, m):
                out.append(pi); m = div(m, pi)
    assert m == (1, 0), (n, m)
    return tuple(out)

def muK(n):
    f = factor_primary(n)
    return 0 if len(set(f)) < len(f) else (-1)**len(f)

def alpha(n):
    z = cx(n); return z/abs(z)

@functools.lru_cache(maxsize=None)
def _sym(xc, pi):
    return sym6_prime(xc, pi)

def sym6(x, pi):
    """Index j of (x/pi)_6 = zeta^j for a primary prime pi coprime to 6; None if pi | x."""
    return _sym(canon_mod(x, pi), pi)

def chi_idx(n, x):
    """Index j of chi_n(x) = (x/n)_6 (n primary, coprime to 6); None if not coprime."""
    j = 0
    for pi in factor_primary(n):
        s = sym6(x, pi)
        if s is None:
            return None
        j += s
    return j % 6

def chi(n, x):
    j = chi_idx(n, x)
    return 0 if j is None else ZETA[j]

def coprime6(n):
    """n coprime to 6 in Z[omega]: not divisible by 2 (A, B not both even) nor by lambda (A + B = 0 mod 3)."""
    return not divides((2, 0), n) and not divides(LAM, n)

def columns(X):
    """All primary n coprime to 6 with N(n) <= X (paper's column set), via the coordinate description
    A = 1 (3), B = 0 (3), A, B not both even (l.383-384)."""
    out = []
    R = int(2*math.sqrt(X/3)) + 3
    for A in range(-R, R + 1):
        if A % 3 != 1:
            continue
        for B in range(-R, R + 1):
            if B % 3 != 0 or (A % 2 == 0 and B % 2 == 0):
                continue
            if 0 < N((A, B)) <= X:
                out.append((A, B))
    return sorted(out, key=lambda z: (N(z), z))

def by_norm(cols):
    d = {}
    for n in cols:
        d.setdefault(N(n), []).append(n)
    return d

# ------------------------------------------------------------------ vectorised sextic symbols
@functools.lru_cache(maxsize=None)
def sym_table(pi):
    """Lookup table for (x + y*omega / pi)_6 on numpy arrays.  Returns (mod, w0, T, kind):
    split pi of norm p: residue r = x + y*w0 mod p (w0 = image of omega), T[r] = index or -1;
    inert pi = -q: residue r = (x mod q)*q + (y mod q), T[r] = index or -1."""
    if pi[1] != 0:
        p = N(pi); a, b = pi
        if b % p == 0:
            raise ValueError(pi)
        w0 = (-a*pow(b, -1, p)) % p
        r = np.arange(p, dtype=np.int64)
        e = (p - 1)//6
        res = np.ones_like(r); base = r.copy()
        while e:
            if e & 1:
                res = (res*base) % p
            base = (base*base) % p; e >>= 1
        T = -np.ones(p, dtype=np.int64)
        for j in range(6):
            u = UNITS[j]; img = (u[0] + u[1]*w0) % p
            T[(res == img) & (r != 0)] = j
        assert np.all(T[1:] >= 0)
        return (p, w0, T, 'split')
    q = -pi[0]
    A, B = np.meshgrid(np.arange(q, dtype=np.int64), np.arange(q, dtype=np.int64), indexing='ij')
    A = A.ravel(); B = B.ravel()
    e = (q*q - 1)//6
    ra = np.ones_like(A); rb = np.zeros_like(A); ba = A.copy(); bb = B.copy()
    while e:
        if e & 1:
            ra, rb = (ra*ba - rb*bb) % q, (ra*bb + rb*ba - rb*bb) % q
        ba, bb = (ba*ba - bb*bb) % q, (2*ba*bb - bb*bb) % q
        e >>= 1
    T = -np.ones(q*q, dtype=np.int64)
    for j in range(6):
        u = UNITS[j]
        T[(ra == u[0] % q) & (rb == u[1] % q) & ((A != 0) | (B != 0))] = j
    assert np.all(T[1:] >= 0)
    return (q, None, T, 'inert')

def sym_vec(pi, x, y):
    """Indices of ((x + y omega)/pi)_6 for integer arrays x, y; -1 where pi divides."""
    m, w0, T, kind = sym_table(pi)
    if kind == 'split':
        return T[((x % m) + (y % m)*w0) % m]
    return T[(x % m)*m + (y % m)]

def row_matrix(xs, ys, cols):
    """Complex matrix C[i, k] = chi_{n_k}(x_i + y_i omega) for squarefree columns n_k (zero where not coprime)."""
    xs = np.asarray(xs, dtype=np.int64); ys = np.asarray(ys, dtype=np.int64)
    C = np.empty((len(xs), len(cols)), dtype=complex)
    cache = {}
    Z = np.array(ZETA + [0j])
    for k, n in enumerate(cols):
        idx = np.zeros(len(xs), dtype=np.int64); zero = np.zeros(len(xs), dtype=bool)
        for pi in factor_primary(n):
            if pi not in cache:
                cache[pi] = sym_vec(pi, xs, ys)
            s = cache[pi]; zero |= (s < 0); idx += np.where(s < 0, 0, s)
        idx %= 6; idx[zero] = 6
        C[:, k] = Z[idx]
    return C

def lattice_points(H, include_zero=False):
    """All (x, y) with 0 < Q(x, y) <= H (Q = x^2 - xy + y^2)."""
    R = int(2*math.sqrt(H/3)) + 2
    x, y = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing='ij')
    x = x.ravel(); y = y.ravel(); q = x*x - x*y + y*y
    keep = (q <= H) & ((q > 0) | include_zero)
    return x[keep], y[keep]

# ------------------------------------------------------------------ weights
def W_test(y):
    """A complex smooth compactly supported weight on [1, 2] (bump times a phase)."""
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y, dtype=complex)
    m = (y > 1) & (y < 2)
    t = y[m]
    out[m] = np.exp(-1.0/((t - 1)*(2 - t)) + 4)*np.exp(1j*t)
    return out

# ------------------------------------------------------------------ PARI/GP
def gp(script, timeout=600):
    """Run a gp script (quiet) and return stdout lines."""
    r = subprocess.run(['gp', '-q', '-f', '-s', '400000000'], input=script, capture_output=True, text=True,
                       timeout=timeout)
    if r.returncode != 0 or r.stderr.strip():
        raise RuntimeError(r.stderr)
    return r.stdout.strip().splitlines()

@functools.lru_cache(maxsize=None)
def ellan(d, M):
    """Dirichlet coefficients a_1..a_M of L(E_d, s), E_d: y^2 = x^3 + d (PARI ellan)."""
    out = gp(f"v=ellan(ellinit([0,0,0,0,{d}]),{M});for(i=1,#v,print(v[i]))")
    return [None] + [int(t) for t in out]

@functools.lru_cache(maxsize=None)
def conductor(d):
    return int(gp(f"print(ellglobalred(ellinit([0,0,0,0,{d}]))[1])")[0])

def ap_count(d, p):
    """a_p(d) = p + 1 - #E_d(F_p) by brute-force point count (affine points + 1 at infinity)."""
    sq = [0]*p
    for y in range(p):
        sq[y*y % p] += 1
    cnt = 1 + sum(sq[(x*x*x + d) % p] for x in range(p))
    return p + 1 - cnt

def ap_legendre(d, p):
    """-sum_x ((x^3 + d)/p)."""
    return -sum(legendre(x*x*x + d, p) for x in range(p))

# ------------------------------------------------------------------ the integer coefficients q_d
def Qd_int(d, m, ap):
    """Q(m) = sqrt(m) * q_d(m), an integer: multiplicative with Q(p) = -a_p, Q(p^2) = p, Q(p^j) = 0 (j >= 3);
    Q(m) = 0 if (m, 6d) > 1.  ap: function p -> a_p(d)."""
    if math.gcd(m, 6*d) > 1:
        return 0
    r = 1
    for p, e in factorint(m).items():
        if e == 1:
            r *= -ap(p)
        elif e == 2:
            r *= p
        else:
            return 0
    return r

def qd(d, m, ap):
    return Qd_int(d, m, ap)/math.sqrt(m)
