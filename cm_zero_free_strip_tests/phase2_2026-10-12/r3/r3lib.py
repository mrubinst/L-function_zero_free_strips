"""Self-contained arithmetic in Z[omega] for the r3 referee tests (independent of tests_cm helpers).
Elements are pairs (a, b) = a + b*omega, omega = e^{2 pi i/3}, omega^2 = -1 - omega.
Conventions of the guide and of OAI78: primary = 1 mod 3; chi_n(x) = (x/n)_6 with zero values
retained; alpha(n) = n/|n|; e(z) = exp(4 pi i Im z / sqrt 3);
gamma_j(c) = N(c)^{-1/2} sum_{v mod c} chi_c(v)^j e(v/c); G(c) = conj(chi_c(4)) gamma_3(c) (squarefree c)."""
import cmath, math, subprocess
from math import gcd

W = complex(-0.5, math.sqrt(3) / 2)

def mul(x, y):
    a, b = x; c, d = y
    # (a+bw)(c+dw) = ac + (ad+bc)w + bd w^2 = (ac-bd) + (ad+bc-bd) w
    return (a * c - b * d, a * d + b * c - b * d)

def conj(x):
    a, b = x
    return (a - b, -b)

def norm(x):
    a, b = x
    return a * a - a * b + b * b

def cx(x):
    return x[0] + x[1] * W

def alpha(x):
    z = cx(x)
    return z / abs(z)

def is_primary(x):
    return x[0] % 3 == 1 and x[1] % 3 == 0

UNITS = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]  # 1, w, w^2, -1, -w, -w^2

def primary_assoc(x):
    for u in UNITS:
        y = mul(u, x)
        if is_primary(y):
            return y
    return None

def is_prime_int(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = int(math.isqrt(n))
    for p in range(3, r + 1, 2):
        if n % p == 0:
            return False
    return True

def primary_primes(B):
    """Primary primes of norm <= B, coprime to 6 (both primes of each split pair, and inert -p)."""
    out = []
    for p in range(5, B + 1):
        if not is_prime_int(p):
            continue
        if p % 3 == 1:
            found = set()
            for b in range(1, int(math.isqrt(4 * p // 3)) + 1):
                t = 4 * p - 3 * b * b
                r = math.isqrt(t)
                if r * r == t:
                    for s_ in (r, -r):
                        if (s_ + b) % 2 == 0:
                            a = (s_ + b) // 2
                            for bb in (b, -b):
                                for aa in (a, a - b + bb if bb != b else a):
                                    if aa * aa - aa * bb + bb * bb == p:
                                        q = primary_assoc((aa, bb))
                                        if q is not None:
                                            found.add(q)
                    # also the conjugates
            found |= {conj(q) for q in found}
            out.extend(sorted(found))
        elif p * p <= B:
            out.append((-p, 0))
    return out

# ---- residue symbols -----------------------------------------------------------------------
def _split_root(pi):
    """For a split prime pi = a+bw of norm p, return (p, r) with pi | (w - r)."""
    a, b = pi
    p = norm(pi)
    r = (-a * pow(b, -1, p)) % p
    assert (r * r + r + 1) % p == 0
    return p, r

def _fp2_mul(x, y, p):
    a, b = x; c, d = y
    return ((a * c - b * d) % p, (a * d + b * c - b * d) % p)

def _fp2_pow(x, e, p):
    r = (1, 0)
    while e:
        if e & 1:
            r = _fp2_mul(r, x, p)
        x = _fp2_mul(x, x, p)
        e >>= 1
    return r

def chi_prime(pi, x):
    """(x/pi)_6 for a primary prime pi coprime to 6, as a unit (a,b) or 0."""
    n = norm(pi)
    if pi[1] != 0:  # split
        p, r = _split_root(pi)
        v = (x[0] + x[1] * r) % p
        if v == 0:
            return 0
        t = pow(v, (p - 1) // 6, p)
        for u in UNITS:
            if (u[0] + u[1] * r - t) % p == 0:
                return u
        raise ValueError
    else:           # inert, pi = -p
        p = -pi[0]
        xx = (x[0] % p, x[1] % p)
        if xx == (0, 0):
            return 0
        t = _fp2_pow(xx, (p * p - 1) // 6, p)
        for u in UNITS:
            if ((u[0] - t[0]) % p, (u[1] - t[1]) % p) == (0, 0):
                return u
        raise ValueError

def chi(primes, x):
    """chi_c(x) for c = product of the listed primary primes (with multiplicity), complex or 0."""
    val = 1 + 0j
    for p in primes:
        c = chi_prime(p, x)
        if c == 0:
            return 0
        val *= cx(c)
    return val

def prod(primes):
    r = (1, 0)
    for p in primes:
        r = mul(r, p)
    return r

def residues(c):
    """A complete residue system mod c: {x + y w : 0<=x<n1, 0<=y<g2}."""
    a, b = c
    g2 = gcd(abs(b), abs(a - b))
    n1 = norm(c) // g2
    return [(x, y) for x in range(n1) for y in range(g2)]

def e_add(z):
    return cmath.exp(4j * math.pi * z.imag / math.sqrt(3))

def gauss(primes, j):
    """gamma_j(c) for squarefree c = prod(primes)."""
    c = prod(primes)
    cz = cx(c)
    s = 0j
    for v in residues(c):
        ch = chi(primes, v)
        if ch == 0:
            continue
        s += ch ** j * e_add(cx(v) / cz)
    return s / math.sqrt(norm(c))

def mu(primes):
    return (-1) ** len(primes) if len(set(primes)) == len(primes) else 0

def G(primes):
    c4 = chi(primes, (4, 0))
    return c4.conjugate() * gauss(primes, 3)

def gp(script, prec=38):
    out = subprocess.run(["gp", "-q", "-f", "-D", f"realprecision={prec}"], input=script,
                         capture_output=True, text=True, timeout=1200)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    return [l for l in out.stdout.split("\n") if l.strip()]

def gauss_all(primes, js=(1, 2, 3, -1)):
    """{j: gamma_j(c)} for squarefree c = prod(primes), computing chi_c(v) once per residue."""
    c = prod(primes)
    cz = cx(c)
    s = {j: 0j for j in js}
    for v in residues(c):
        ch = chi(primes, v)
        if ch == 0:
            continue
        ev = e_add(cx(v) / cz)
        for j in js:
            s[j] += ch ** j * ev
    r = math.sqrt(norm(c))
    return {j: s[j] / r for j in js}
