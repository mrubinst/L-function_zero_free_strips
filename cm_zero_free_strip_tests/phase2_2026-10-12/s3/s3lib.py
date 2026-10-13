"""Shared arithmetic for the s3 tests (the zero detector with the label set X_u; CM guide, item 3 of Section 16).

Conventions (OAI78 / the guide): omega = e^{2 pi i/3}, O = Z[omega], lambda = 1 + 2 omega, primary = 1 mod 3,
chi_n(u) = (u/n)_6 multiplicative in n with the zero convention (0 when (n,u) != 1), alpha(n) = n/|n|.
Elements are pairs (a, b) = a + b*omega (see lib/zw.py, copied unchanged from tests_cm/part_C_lines_0998_1509/lib).

Factorisation of elements of large norm uses sympy.factorint on the rational norm, so composite moduli of norm up to
~1e12 can be handled (the copied common.py only factors elements whose prime factors have norm <= 2e5).
"""
import os, sys, math, cmath, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import zw
from zw import mul, add, sub, conj, N, cx, divides, div, pw, primary, prim_assoc, UNITS, LAM, reduce_mod
from sympy import factorint, sqrt_mod
from functools import lru_cache

ZETA6 = [cmath.exp(1j*math.pi*j/3) for j in range(6)]


def alpha(x):
    z = cx(x)
    return z/abs(z)


def egcd(x, y):
    """A gcd of x, y in Z[omega] (Euclidean with nearest-point division)."""
    while y != (0, 0):
        x, y = y, reduce_mod(x, y)
    return x


@lru_cache(maxsize=None)
def split_prime(p):
    """The primary prime pi of norm p (p = 1 mod 3) with the smaller (a, b) representation; its conjugate is the other."""
    r = int(sqrt_mod(-3, p))                  # r^2 = -3 mod p; omega = (-1 + sqrt(-3))/2
    w0 = ((-1 + r)*pow(2, -1, p)) % p        # a cube root of unity mod p
    g = egcd((p, 0), sub((0, 1), (w0, 0)))  # gcd(p, omega - w0)
    assert N(g) == p, (p, g)
    return prim_assoc(g)


def factor_elt(x):
    """Factor a nonzero x in Z[omega]: list of (prime, exponent) with primes LAM (above 3), (-2, 0) (above 2),
    (-q, 0) for inert q = 2 mod 3, and primary split primes; the unit left over is ignored."""
    out = []
    m = x
    for p, e in sorted((int(a), int(b)) for a, b in factorint(int(N(x))).items()):
        if p == 3:
            v = 0
            while divides(LAM, m):
                m = div(m, LAM); v += 1
            out.append((LAM, v))
        elif p % 3 == 2:
            q = (-p, 0); v = 0
            while divides(q, m):
                m = div(m, q); v += 1
            out.append((q, v))
        else:
            pi = split_prime(p)
            for P in (pi, prim_assoc(conj(pi))):
                v = 0
                while divides(P, m):
                    m = div(m, P); v += 1
                if v:
                    out.append((P, v))
    assert N(m) == 1, (x, m)
    return out


@lru_cache(maxsize=None)
def sym6(u, pi):
    """(u/pi)_6 as an index j (zeta6^j), None if pi | u; pi primary prime coprime to 6."""
    return zw.sym6_prime(u, pi)


def chi(n, u):
    """chi_n(u) = (u/n)_6 for n coprime to 6 (any associate; the symbol depends only on the ideal), zero convention."""
    r = 1 + 0j
    for P, v in factor_elt(n):
        j = sym6(u, P)
        if j is None:
            return 0
        r *= ZETA6[(j*v) % 6]
    return r


def chipow(c, s):
    """c^s with the zero convention (0 stays 0)."""
    return 0 if c == 0 else c**s


def coprime_to(x, primes):
    return all(not divides(P, x) for P in primes)


class Member:
    """A presentation psi_{u,nu,sigma}(n) = alpha(n)^{k sigma} theta(n) chi_n(u)^sigma, n good primary.
    theta(n) = prod_b chi_n(b)^{e_b} over fixed b supported on {2, 3} (a finite-order ray character with conductor in S);
    k = 0 (finite-order member) or k = 1 (CM member, co-oriented: alpha^sigma with chi^sigma)."""

    def __init__(self, u_factors, unit=(1, 0), theta=(), sigma=-1, k=1, alpha_exp=None):
        self.u_factors = list(u_factors)          # [(prime, j)], primes not dividing 6, 1 <= j <= 5
        u = unit
        for P, j in self.u_factors:
            u = mul(u, pw(P, j))
        self.u = u
        self.theta = tuple(theta)                 # ((b, e_b), ...)
        self.sigma = sigma
        self.k = k
        self.aexp = k*sigma if alpha_exp is None else alpha_exp   # exponent of alpha (alpha_exp used for falsifiers)
        self.bad = [LAM, (-2, 0)] + [P for P, _ in self.u_factors]

    def theta_val(self, n):
        r = 1 + 0j
        for b, e in self.theta:
            r *= chi(n, b)**e
        return r

    def __call__(self, n):
        """Value at a good primary n (coprime to 6)."""
        c = chi(n, self.u)
        if c == 0:
            return 0
        return alpha(n)**self.aexp*self.theta_val(n)*c**self.sigma

    def conj_member(self):
        """The label (conj nu, -sigma): conj theta = theta^{-1}, alpha^{k sigma} -> alpha^{-k sigma}."""
        th = tuple((b, -e) for b, e in self.theta)
        return _with_u(self, th, -self.sigma)


def _with_u(m, theta, sigma):
    new = Member.__new__(Member)
    new.u_factors = m.u_factors; new.u = m.u; new.theta = theta; new.sigma = sigma; new.k = m.k
    new.aexp = m.k*sigma; new.bad = m.bad
    return new


def sample_valid(g, member, target, rng, samples=50, R=30, maxtries=4000):
    """True if member(prim_assoc(x)) == target(x) for `samples` random x = 1 + g*y coprime to 6u."""
    cnt = 0; tries = 0
    while cnt < samples and tries < maxtries:
        tries += 1
        y = (rng.randint(-R, R), rng.randint(-R, R))
        x = add((1, 0), mul(g, y))
        if N(x) == 0 or not coprime_to(x, member.bad):
            continue
        if abs(member(prim_assoc(x)) - target(x)) > 1e-8:
            return False
        cnt += 1
    if cnt < samples:
        raise RuntimeError("not enough samples for modulus %s" % (g,))
    return True


def modulus(exps, primes):
    g = (1, 0)
    for P, a in zip(primes, exps):
        g = mul(g, pw(P, a))
    return g


def conductor(member, target, emax=(4, 6), rng=None, samples=50):
    """Conductor exponents of the Hecke character inducing `member`, assuming f | 2^emax[0] lam^emax[1] prod pi^2.
    `target(x)` is the infinity-type value psi((x)) for x = 1 mod f (1 for finite order, alpha(x)^sigma for CM).
    Returns (dict prime -> exponent, norm) or None if the big modulus is not valid."""
    rng = rng or random.Random(1)
    primes = [(-2, 0), LAM] + [P for P, _ in member.u_factors]
    big = [emax[0], emax[1]] + [2]*len(member.u_factors)
    if not sample_valid(modulus(big, primes), member, target, rng, samples):
        return None
    ex = list(big)
    for i in range(len(primes)):
        for a in range(0, big[i] + 1):
            trial = list(big); trial[i] = a
            if sample_valid(modulus(trial, primes), member, target, rng, samples):
                ex[i] = a
                break
    f = modulus(ex, primes)
    assert sample_valid(f, member, target, rng, samples)
    return dict(zip([str(P) for P in primes], ex)), N(f), f


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


def report(cid, ref, desc, ok, details="", falsifier_failed=None, status=None):
    """One line per claim: [STATUS] <id> <ref> <desc> :: <details>.  PASS needs ok and (if a falsifier is given)
    falsifier_failed; a falsifier that passes turns the line into FAIL (the test would have no power)."""
    if status is None:
        if falsifier_failed is None:
            status = "PASS" if ok else "FAIL"
        else:
            status = "PASS" if (ok and falsifier_failed) else "FAIL"
            details += " [falsifier failed as required]" if falsifier_failed else " [FALSIFIER DID NOT FAIL]"
    if _forced():
        status = "FAIL"; details = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + details
    if status == "FAIL":
        _NFAILED[0] += 1
    line = f"[{status}] {cid} {ref} {desc} :: {details}"
    print(line, flush=True)
    RESULTS.append((status, cid))
    return status


def summary():
    from collections import Counter
    c = Counter(s for s, _ in RESULTS)
    print("SUMMARY:", dict(c))
    return c.get("FAIL", 0) == 0
