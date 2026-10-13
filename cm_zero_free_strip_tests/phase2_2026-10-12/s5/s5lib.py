"""Helpers for the s5 tests (plain fourth moment with an angular factor).

Conventions (those of the CM guide and of OAI78, Section sec:arithmetic, L557-605):
  omega = exp(2 pi i/3); an Eisenstein integer a + b*omega is stored as a pair (a, b) of Python ints;
  primary = congruent to 1 mod 3; lambda = 1 + 2 omega = i sqrt 3;
  e(z) = exp(4 pi i Im z / sqrt 3)  (the pairing for which Z[omega] is self-dual, [QF, Prop. 16]);
  alpha(n) = n/|n|;  chi_p(x) = (x/p)_6, zero on multiples of p.
Poisson ([QF, Thm. 19]):  sum_{x in O} f(x) = (2/sqrt3) sum_{y in O} fhat(y),  fhat(y) = int f(z) e(-zy) d^2z.
"""
import os, sys, math, cmath
import numpy as np

SQ3 = math.sqrt(3.0)
OM = complex(-0.5, SQ3 / 2)                       # omega
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]   # zeta^j, zeta = 1 + omega = -omega^2

def cx(x):
    return x[0] + x[1] * OM

def mul(x, y):
    a, b = x; c, d = y                            # (a + b w)(c + d w), w^2 = -1 - w
    return (a * c - b * d, a * d + b * c - b * d)

def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def neg(x): return (-x[0], -x[1])
def conj(x):                                      # conj(a + b w) = a + b w^2 = (a - b) - b w
    return (x[0] - x[1], -x[1])
def N(x):
    a, b = x
    return a * a - a * b + b * b

def divmod_e(x, m):
    """Euclidean division: x = q m + r with N(r) < N(m)."""
    z = cx(x) / cx(m)
    # nearest lattice point to z: solve z = u + v w
    v = z.imag / (SQ3 / 2); u = z.real + v / 2
    best = None
    for du in (math.floor(u), math.ceil(u)):
        for dv in (math.floor(v), math.ceil(v)):
            q = (du, dv); r = sub(x, mul(q, m))
            if best is None or N(r) < N(best[1]):
                best = (q, r)
    return best

def mod(x, m):
    return divmod_e(x, m)[1]

def divides(d, x):
    return N(mod(x, d)) == 0

def exact_div(x, d):
    q, r = divmod_e(x, d); assert N(r) == 0; return q

def is_primary(x):
    return (x[0] - 1) % 3 == 0 and x[1] % 3 == 0

def primary_assoc(x):
    """(u, x0) with x = u * x0, u a unit, x0 primary.  Requires (x, 3) = 1."""
    for u in UNITS:
        y = mul(x, conj(u))   # conj(u) = u^{-1} for units
        if is_primary(y):
            return u, y
    raise ValueError("not coprime to 3")

def unit_index(u):
    return UNITS.index(u)

def ucx(u):
    return cx(u)

def e_add(z):
    return cmath.exp(4j * math.pi * z.imag / SQ3)

# ---------------------------------------------------------------- residues and residue symbols
def residues(m):
    """A complete residue system mod m (m an Eisenstein integer), as a list of pairs."""
    n = N(m); seen = {}; out = []
    R = int(math.isqrt(n)) + 2
    for a in range(-2 * R, 2 * R + 1):
        for b in range(-2 * R, 2 * R + 1):
            r = mod((a, b), m)
            if r not in seen:
                seen[r] = 1; out.append(r)
            if len(out) == n:
                return out
    assert len(out) == n, (len(out), n)
    return out

def canon(x, m):
    """Canonical representative of x mod m (deterministic)."""
    return mod(x, m)

def powmod(x, e, m):
    r = (1, 0); b = mod(x, m)
    while e:
        if e & 1: r = mod(mul(r, b), m)
        b = mod(mul(b, b), m); e >>= 1
    return r

def chi_p(x, p):
    """Sextic symbol (x/p)_6 for a prime p of O not dividing 6; 0 if p | x."""
    if divides(p, x):
        return 0
    r = powmod(x, (N(p) - 1) // 6, p)
    for j, u in enumerate(UNITS):
        if N(sub(r, u)) % N(p) == 0 and divides(p, sub(r, u)):
            return cmath.exp(1j * math.pi * j / 3)
    raise ValueError("symbol not a sixth root of unity")

def primes_upto(B):
    """Primary primes of O coprime to 6 with norm <= B (one generator per prime ideal)."""
    out = []
    for p in range(5, B + 1):
        if any(p % q == 0 for q in range(2, int(math.isqrt(p)) + 1)):
            continue
        if p % 3 == 1:
            # find a + b w of norm p
            found = None
            R = int(math.isqrt(4 * p)) + 2
            for a in range(-R, R + 1):
                for b in range(-R, R + 1):
                    if a * a - a * b + b * b == p:
                        found = (a, b); break
                if found: break
            u, x0 = primary_assoc(found)
            out.append(x0); out.append(conj(x0))
        elif p % 3 == 2 and p * p <= B:
            out.append((-p, 0))
    return out

# ---------------------------------------------------------------- reporting
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

def report(cid, err, tol, detail=""):
    ok = err <= tol
    if _forced():
        ok = False; detail = "CM_TESTS_FORCE_FAIL=1 (forced failure); " + detail
    if not ok:
        _NFAILED[0] += 1
    RESULTS.append((cid, ok))
    print(f"RESULT|{cid}|{'PASS' if ok else 'FAIL'}|err={err:.3g} tol={tol:g}; {detail}", flush=True)
    return ok

def falsifier(cid, name, err, must_exceed):
    ok = err > must_exceed
    if _forced():
        ok = False; name = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + name
    if not ok:
        _NFAILED[0] += 1
    RESULTS.append((cid + ":falsifier:" + name, ok))
    print(f"      falsifier for {cid} [{name}]: err={err:.3g} -> "
          f"{'fails as required' if ok else 'DID NOT FAIL (bad falsifier)'}", flush=True)
    return ok

def summary():
    bad = [c for c, ok in RESULTS if not ok]
    print(f"SUMMARY: {len(RESULTS)} checks, {len(bad)} failures" + (f": {bad}" if bad else ""))
    return len(bad) == 0
