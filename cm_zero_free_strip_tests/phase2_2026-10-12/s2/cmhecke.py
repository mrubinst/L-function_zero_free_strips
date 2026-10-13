"""Small exact library for Hecke characters of K=Q(sqrt(-3)) with infinity type zbar/|z|.

Elements of O=Z[omega] are pairs (a,b) meaning a+b*omega, omega=exp(2 pi i/3).
A character in F_CM is given by its finite part phi: (O/f)^x -> C^x with phi(zeta)=zeta
for the six units; then psi((x)) = phi(x) * conj(x)/|x|  (well defined on ideals).

Conventions follow the guide / [QF]: e(z)=exp(4 pi i Im z/sqrt3); Poisson
sum_{x in O} F(x) = (2/sqrt3) sum_{y in O} Fhat(y), Fhat(y)=int F(z) e(-zy) d^2z.
"""
import cmath
import math

OMEGA = complex(-0.5, math.sqrt(3) / 2)
UNITS = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]  # 1,-1,w,-w,w^2=-1-w,-w^2=1+w


def cval(x):
    return x[0] + x[1] * OMEGA


def mul(x, y):
    a, b = x
    c, d = y
    # (a+bw)(c+dw) = ac + (ad+bc) w + bd w^2, w^2 = -1-w
    return (a * c - b * d, a * d + b * c - b * d)


def conj(x):
    a, b = x
    # conj(a+bw) = a + b w^2 = a - b - b w
    return (a - b, -b)


def norm(x):
    a, b = x
    return a * a - a * b + b * b


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def divisible(x, f):
    """x in f*O ?"""
    n = norm(f)
    y = mul(x, conj(f))
    return y[0] % n == 0 and y[1] % n == 0


def key_mod(x, f):
    """canonical key of the class of x modulo f (injective on O/f)."""
    n = norm(f)
    y = mul(x, conj(f))
    return (y[0] % n, y[1] % n)


def residue_reps(f):
    """one representative of each class of O/f."""
    n = norm(f)
    reps = {}
    for a in range(n):
        for b in range(n):
            k = key_mod((a, b), f)
            if k not in reps:
                reps[k] = (a, b)
        if len(reps) == n:
            break
    assert len(reps) == n
    return list(reps.values())


def coprime(x, f):
    """gcd((x),(f)) = 1 ?  (test: x invertible mod f)"""
    # brute force: x invertible mod f iff some y with x*y = 1 mod f
    for y in residue_reps(f):
        if divisible(sub(mul(x, y), (1, 0)), f):
            return True
    return False


def unit_mod3(x):
    """the unit zeta with zeta = x mod 3 (x prime to 3); returns its index in UNITS."""
    for i, u in enumerate(UNITS):
        if divisible(sub(x, u), (3, 0)):
            return i
    return None


def unit_power_value(i, k):
    """value (as complex) of UNITS[i]**k"""
    return cval(UNITS[i]) ** k


class SexticSymbol:
    """chi_pi(x) = (x/pi)_6 for a primary prime pi of prime norm p = 1 mod 6."""

    def __init__(self, pi):
        self.pi = pi
        self.p = norm(pi)
        a, b = pi
        p = self.p
        assert b % p != 0
        # omega = -a/b mod pi
        self.w = (-a * pow(b, -1, p)) % p
        self.unit_res = [((u[0] + u[1] * self.w) % p) for u in UNITS]

    def __call__(self, x):
        p = self.p
        r = (x[0] + x[1] * self.w) % p
        if r == 0:
            return 0
        t = pow(r, (p - 1) // 6, p)
        i = self.unit_res.index(t)
        return cval(UNITS[i])


class CMChar:
    """psi((x)) = phi(x) conj(x)/|x| with phi a primitive character mod f, phi(unit)=unit."""

    def __init__(self, name, f, phi):
        self.name = name
        self.f = f
        self.Q = norm(f)
        self.phi = phi  # function of element, 0 if not coprime to f
        self.reps = residue_reps(f)
        # check phi(zeta) = zeta
        for u in UNITS:
            assert abs(phi(u) - cval(u)) < 1e-12, (name, u, phi(u))

    def gauss(self):
        """gamma(phi) = Q^{-1/2} sum_{x mod f} phi(x) e(x/f), e(z)=exp(4 pi i Im z/sqrt3)."""
        fc = cval(self.f)
        s = 0
        for x in self.reps:
            v = self.phi(x)
            if v != 0:
                z = cval(x) / fc
                s += v * cmath.exp(4j * math.pi * z.imag / math.sqrt(3))
        return s / math.sqrt(self.Q)

    def root_number(self):
        """W(psi) = - gamma(phi) sqrt(Q) / f."""
        return -self.gauss() * math.sqrt(self.Q) / cval(self.f)

    def psi_elem(self, x):
        v = self.phi(x)
        if v == 0:
            return 0
        z = cval(x)
        return v * z.conjugate() / abs(z)


def lattice(nmax):
    """all nonzero x in O with N(x) <= nmax."""
    out = []
    r = int(math.isqrt(4 * nmax // 3 + 4)) + 2
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a * a - a * b + b * b
            if 0 < n <= nmax:
                out.append((a, b))
    return out


def primary_primes_of_norm(p):
    out = []
    r = int(math.isqrt(4 * p // 3 + 4)) + 2
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            if a * a - a * b + b * b == p and divisible(sub((a, b), (1, 0)), (3, 0)):
                out.append((a, b))
    return out


def make_test_chars():
    chars = []
    # (A) psi = alphabar: phi(x) = unit = x mod 3, f = 3, Q = 9 (L(s,psi)=L(27a,s+1/2))
    def phiA(x):
        i = unit_mod3(x)
        return 0 if i is None else cval(UNITS[i])
    chars.append(CMChar("A: alphabar, f=3", (3, 0), phiA))
    # (B) phi = chi_pi, N pi = 7, f = pi, Q = 7 (unramified at sqrt(-3))
    pi7 = primary_primes_of_norm(7)[0]
    chi7 = SexticSymbol(pi7)
    chars.append(CMChar("B: chi_pi, N(pi)=7, f=pi", pi7, lambda x: chi7(x)))
    # (C) phi = zeta(x)^5 chi_pi(x), N pi = 13, f = 3 pi, Q = 117
    pi13 = primary_primes_of_norm(13)[0]
    chi13 = SexticSymbol(pi13)

    def phiC(x):
        i = unit_mod3(x)
        if i is None:
            return 0
        c = chi13(x)
        return 0 if c == 0 else cval(UNITS[i]) ** 5 * c
    chars.append(CMChar("C: zeta^5 chi_pi, N(pi)=13, f=3pi", mul((3, 0), pi13), phiC))
    # (D) phi = zeta(x)^3 chi_pi(x)^4, N pi = 7, f = lambda pi, Q = 21 (exponent 1 at lambda)
    lam = (1, 2)

    def phiD(x):
        i = unit_mod3(x)
        if i is None:
            return 0
        c = chi7(x)
        return 0 if c == 0 else cval(UNITS[i]) ** 3 * c ** 4
    chars.append(CMChar("D: zeta^3 chi_pi^4, N(pi)=7, f=lambda pi", mul(lam, pi7), phiD))
    return chars


def is_primitive(ch):
    """phi primitive mod f: for every proper divisor g of f (g = f/p, p prime | f),
    phi is nontrivial on {x = 1 mod g, x prime to f}."""
    f = ch.f
    # prime divisors of f: test the candidate primes lambda, and primes of norm | Q
    Q = ch.Q
    cands = []
    for p in range(2, Q + 1):
        if Q % p == 0 and all(p % q for q in range(2, int(p ** 0.5) + 1)):
            for x in lattice(p):
                if norm(x) == p and divisible(f, x):
                    cands.append(x)
            if divisible(f, (p, 0)) and norm((p, 0)) == p * p:
                cands.append((p, 0))
    ok = True
    for pr in cands:
        g = None
        # g = f / pr
        fc = cval(f) / cval(pr)
        a = round(fc.real + fc.imag / math.sqrt(3))
        b = round(2 * fc.imag / math.sqrt(3))
        g = (a, b)
        assert mul(g, pr) == f, (g, pr, f)
        nontriv = False
        for x in ch.reps:
            if ch.phi(x) != 0 and divisible(sub(x, (1, 0)), g) and abs(ch.phi(x) - 1) > 1e-9:
                nontriv = True
                break
        ok = ok and nontriv
    return ok
