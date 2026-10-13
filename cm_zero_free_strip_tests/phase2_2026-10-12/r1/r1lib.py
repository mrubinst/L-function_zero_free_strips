"""r1lib: a small, independent implementation of the arithmetic of Z[omega] used by the r1 tests.

Written from scratch (it does not import the tests_cm helpers), so that the r1 tests are an
independent check of the guide's conventions:
  * an element a + b*omega is the pair (a, b), omega = e^{2 pi i/3};
  * primary means = 1 (mod 3), i.e. a = 1, b = 0 (mod 3);
  * chi_n(x) = (x/n)_6, the sixth root of unity = x^{(N pi - 1)/6} (mod pi) at each prime pi | n,
    extended multiplicatively in n, with value 0 when (x, n) != 1 (zero values kept for every power);
  * alpha(n) = n/|n|;  e(z) = exp(4 pi i Im z / sqrt 3);
  * gamma_j(n) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n), with the primary generator n;
  * G(n) = conj(chi_n(4)) gamma_3(n).
Exponents of sixth roots of unity are stored as integers j mod 6, meaning zeta6^j, zeta6 = e^{i pi/3};
the value 0 of a symbol is represented by None.
"""
import cmath
import math
from functools import lru_cache

SQ3 = math.sqrt(3.0)
OMEGA = complex(-0.5, SQ3 / 2)


# ---------------------------------------------------------------- ring arithmetic
def mul(x, y):
    a, b = x
    c, d = y
    # (a + b w)(c + d w) = ac + (ad + bc) w + bd w^2,  w^2 = -1 - w
    return (a * c - b * d, a * d + b * c - b * d)


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def neg(x):
    return (-x[0], -x[1])


def conj(x):
    # conj(a + b w) = a + b w^2 = (a - b) - b w
    return (x[0] - x[1], -x[1])


def N(x):
    a, b = x
    return a * a - a * b + b * b


def cx(x):
    return x[0] + x[1] * OMEGA


def pw(x, e):
    r = (1, 0)
    for _ in range(e):
        r = mul(r, x)
    return r


ONE = (1, 0)
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]  # zeta6^j for j = 0..5, zeta6 = 1 + omega


def divides(d, x):
    """d | x in Z[omega]."""
    n = N(d)
    y = mul(x, conj(d))
    return y[0] % n == 0 and y[1] % n == 0


def exact_div(x, d):
    n = N(d)
    y = mul(x, conj(d))
    assert y[0] % n == 0 and y[1] % n == 0
    return (y[0] // n, y[1] // n)


def reduce_mod(x, m):
    """A representative of x modulo m (nearest-lattice-point rounding)."""
    n = N(m)
    y = mul(x, conj(m))  # x/m = y/n
    qa = round(y[0] / n)
    qb = round(y[1] / n)
    return sub(x, mul(m, (qa, qb)))


def congruent(x, y, m):
    return divides(m, sub(x, y))


def is_primary(x):
    return x[0] % 3 == 1 and x[1] % 3 == 0


def primary_assoc(x):
    for u in UNITS:
        y = mul(u, x)
        if is_primary(y):
            return y
    raise ValueError("not coprime to 3: %r" % (x,))


def residues(m):
    """A complete residue system modulo m (N(m) elements), via the Hermite form of the lattice mZ[omega]."""
    v1 = m
    v2 = mul(m, (0, 1))
    d2 = math.gcd(v1[1], v2[1])
    n = N(m)
    d1 = n // d2
    return [(x, y) for y in range(d2) for x in range(d1)]


def powmod(x, e, m):
    r = (1, 0)
    base = reduce_mod(x, m)
    while e > 0:
        if e & 1:
            r = reduce_mod(mul(r, base), m)
        base = reduce_mod(mul(base, base), m)
        e >>= 1
    return r


# ---------------------------------------------------------------- primes and factorisation
def _isprime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = int(math.isqrt(n))
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


@lru_cache(maxsize=None)
def prime_above(p):
    """For a rational prime p = 1 (mod 3): the two primary primes of norm p (pi, conj pi)."""
    r = int(math.isqrt(4 * p // 3)) + 1
    for b in range(1, r + 1):
        disc = 4 * p - 3 * b * b
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s == disc and (b + s) % 2 == 0:
            a = (b + s) // 2
            pi = primary_assoc((a, b))
            assert N(pi) == p
            return pi, primary_assoc(conj(pi))
    raise ValueError(p)


@lru_cache(maxsize=None)
def primary_primes(bound):
    """All primary primes pi of Z[omega] with pi coprime to 6 and N(pi) <= bound, sorted by norm."""
    out = []
    for p in range(5, bound + 1):
        if not _isprime(p):
            continue
        if p % 3 == 1:
            out.extend(prime_above(p))
        elif p * p <= bound:  # p = 2 (mod 3), inert, primary generator -p
            out.append((-p, 0))
    out.sort(key=N)
    return tuple(out)


def _rational_factor(n):
    res = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            res[d] = res.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        res[n] = res.get(n, 0) + 1
    return res


@lru_cache(maxsize=None)
def _factor_cached(n):
    m = n
    res = {}
    for p, k in sorted(_rational_factor(N(n)).items()):
        if p in (2, 3):
            raise ValueError("not coprime to 6: %r" % (n,))
        if p % 3 == 2:
            assert k % 2 == 0
            res[(-p, 0)] = k // 2
            for _ in range(k // 2):
                m = exact_div(m, (-p, 0))
        else:
            for pi in prime_above(p):
                while divides(pi, m):
                    m = exact_div(m, pi)
                    res[pi] = res.get(pi, 0) + 1
    assert m == (1, 0), (n, m)
    return tuple(sorted(res.items()))


def factor(n):
    """Factor a primary n coprime to 6 into primary primes: returns dict prime -> exponent."""
    assert is_primary(n), n
    return dict(_factor_cached(n))


def mu(n):
    f = factor(n)
    if any(e > 1 for e in f.values()):
        return 0
    return (-1) ** len(f)


def primaries(bound, squarefree=False):
    """Primary n coprime to 6 with N(n) <= bound (n = 1 included)."""
    out = []
    r = int(math.isqrt(4 * bound // 3)) + 2
    for b in range(-r, r + 1):
        if b % 3:
            continue
        for a in range(-r, r + 1):
            if a % 3 != 1:
                continue
            x = (a, b)
            nx = N(x)
            if nx > bound:
                continue
            if a % 2 == 0 and b % 2 == 0:
                continue
            if squarefree and mu(x) == 0:
                continue
            out.append(x)
    out.sort(key=lambda x: (N(x), x))
    return out


def gcd_primary(a, b):
    """Primary gcd of primary a, b (via factorisation)."""
    fa = factor(a)
    fb = factor(b)
    g = (1, 0)
    for p, e in fa.items():
        if p in fb:
            g = mul(g, pw(p, min(e, fb[p])))
    return g


# ---------------------------------------------------------------- sextic symbol
def _unit_index_mod(y, pi):
    for j, u in enumerate(UNITS):
        if congruent(y, u, pi):
            return j
    raise ValueError("Euler criterion did not give a sixth root of unity")


@lru_cache(maxsize=None)
def chi_prime_idx(pi, x):
    """Index j with (x/pi)_6 = zeta6^j, or None if pi | x."""
    if divides(pi, x):
        return None
    y = powmod(x, (N(pi) - 1) // 6, pi)
    return _unit_index_mod(y, pi)


def chi_idx(n, x):
    """Index of chi_n(x) (sum over primes with multiplicity), or None if (n, x) != 1."""
    tot = 0
    for p, e in factor(n).items():
        j = chi_prime_idx(p, x)
        if j is None:
            return None
        tot += e * j
    return tot % 6


ZETA6 = [cmath.exp(1j * math.pi * j / 3) for j in range(6)]


def chi(n, x, power=1):
    j = chi_idx(n, x)
    if j is None:
        return 0j
    return ZETA6[(power * j) % 6]


def alpha(n):
    z = cx(n)
    return z / abs(z)


def e_add(z):
    return cmath.exp(4j * math.pi * z.imag / SQ3)


@lru_cache(maxsize=None)
def gamma(j, n):
    """gamma_j(n), brute force over a complete residue system modulo the primary generator n."""
    if n == (1, 0):
        return 1.0 + 0j
    s = 0j
    zn = cx(n)
    fac = factor(n)
    for x in residues(n):
        tot = 0
        ok = True
        for p, e in fac.items():
            k = chi_prime_idx(p, x)
            if k is None:
                ok = False
                break
            tot += e * k
        if not ok:
            continue
        s += ZETA6[(j * tot) % 6] * e_add(cx(x) / zn)
    return s / math.sqrt(N(n))


def gauss_char(z1, z2, conj_first=True):
    """gamma(conj(chi_z1) chi_z2) (conj_first) or gamma(chi_z1 conj(chi_z2)), modulus q = z1 z2 (primary)."""
    q = mul(z1, z2)
    if q == (1, 0):
        return 1.0 + 0j
    zq = cx(q)
    f1 = factor(z1) if z1 != (1, 0) else {}
    f2 = factor(z2) if z2 != (1, 0) else {}
    s1 = -1 if conj_first else 1
    s = 0j
    for x in residues(q):
        tot = 0
        ok = True
        for p, e in f1.items():
            k = chi_prime_idx(p, x)
            if k is None:
                ok = False
                break
            tot += s1 * e * k
        if not ok:
            continue
        for p, e in f2.items():
            k = chi_prime_idx(p, x)
            if k is None:
                ok = False
                break
            tot -= s1 * e * k
        if not ok:
            continue
        s += ZETA6[tot % 6] * e_add(cx(x) / zq)
    return s / math.sqrt(N(q))


def G(n):
    return chi(n, (4, 0)).conjugate() * gamma(3, n)


# ---------------------------------------------------------------- classes modulo 4 and modulo 72
def cls4(x):
    return (x[0] % 4, x[1] % 4)


def units_mod4():
    """The 12 classes of (Z[omega]/4)^x, as pairs (a mod 4, b mod 4) with a + b w odd (not both even)."""
    out = []
    for a in range(4):
        for b in range(4):
            if a % 2 == 0 and b % 2 == 0:
                continue
            out.append((a, b))
    return out


def mul4(x, y):
    return cls4(mul(x, y))


def inv4(x):
    for y in units_mod4():
        if mul4(x, y) == (1, 0):
            return y
    raise ValueError


_FORCE_FAIL = [__import__("os").environ.get("CM_TESTS_FORCE_FAIL") == "1"]   # self-test (added 2026-10-12)


def _forced():
    """True exactly once (for the first check) when CM_TESTS_FORCE_FAIL=1."""
    if _FORCE_FAIL[0]:
        _FORCE_FAIL[0] = False
        return True
    return False


class Reporter:
    def __init__(self, name):
        self.name = name
        self.rows = []

    def check(self, cid, ok, detail):
        if _forced():
            ok = False; detail = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + detail
        tag = "PASS" if ok else "FAIL"
        line = "[%s] %s %s" % (tag, cid, detail)
        print(line, flush=True)
        self.rows.append((cid, tag, detail))
        return ok

    def falsifier(self, cid, name, fails, detail):
        if _forced():
            fails = False; name = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + name
        tag = "fails as required" if fails else "DID NOT FAIL (bad falsifier)"
        line = "      falsifier %s [%s]: %s -> %s" % (cid, name, detail, tag)
        print(line, flush=True)
        self.rows.append((cid + "-F", "PASS" if fails else "FAIL", name + ": " + detail))
        return fails

    def summary(self):
        nf = sum(1 for r in self.rows if r[1] == "FAIL")
        print("SUMMARY %s: %d checks, %d FAIL" % (self.name, len(self.rows), nf))
        return nf
