"""Hecke characters of K = Q(sqrt-3) of infinity type (a/|a|)^e, given through their finite part.

A character xi of infinity type (a/|a|)^e and modulus m is xi((a)) = phi(a) (a/|a|)^e for a coprime to m, where phi is a
character of (O/m)^x with phi(u) = u^{-e} on the six units.  We build phi from a function on primary elements
(xi_prim(a0), a0 = primary generator), reduce it to its conductor, and evaluate sums over elements with numpy.
"""
import math, cmath
import numpy as np
from s5lib import (mul, sub, add, conj, N, cx, primary_assoc, divides, exact_div, e_add, UNITS, OM, SQ3)

def lattice_hnf(m):
    """Lattice m*O in (a,b)-coordinates has a basis (n1, 0), (k, n2).  Returns (n1, k, n2)."""
    c, d = m
    g = math.gcd(c, d)
    # find x, y with x*d + y*(c-d) = g  (b-coordinate of x*m + y*m*omega)
    def egcd(a, b):
        if b == 0: return (a, 1, 0)
        q, r = divmod(a, b); gg, s, t = egcd(b, r); return (gg, t, s - q * t)
    gg, x, y = egcd(d, c - d)
    if gg < 0: gg, x, y = -gg, -x, -y
    assert gg == g
    v = add(mul((x, 0), m), mul((y, 0), mul(m, (0, 1))))
    assert v[1] == g
    n1 = N(m) // g
    return n1, v[0] % n1, g

def canon_index(A, B, m, hnf=None):
    """Canonical index in [0, N(m)) of the class of A + B omega mod m (numpy int arrays)."""
    n1, k, n2 = hnf if hnf is not None else lattice_hnf(m)
    A = np.asarray(A, dtype=np.int64); B = np.asarray(B, dtype=np.int64)
    t = np.floor_divide(B, n2)
    A = A - t * k; B = B - t * n2
    A = np.mod(A, n1)
    return A + n1 * B

def residue_reps(m):
    n1, k, n2 = lattice_hnf(m)
    reps = [(a, b) for b in range(n2) for a in range(n1)]   # index a + n1*b
    return reps

def coprime(x, primes):
    return all(not divides(p, x) for p in primes)

class Char:
    """phi on (O/m)^x stored as a table over canonical indices; e is the infinity exponent."""
    def __init__(self, m, primes, e, table):
        self.m = m; self.primes = list(primes); self.e = e; self.table = np.asarray(table, dtype=complex)
        self.hnf = lattice_hnf(m)

    def __call__(self, A, B):
        return self.table[canon_index(A, B, self.m, self.hnf)]

def from_primary_function(xi_prim, m, primes, e, check=True):
    """phi(x) = xi((x)) (x/|x|)^{-e}, xi((x)) = xi_prim(primary associate of x)."""
    reps = residue_reps(m)
    tab = np.zeros(len(reps), dtype=complex)
    for i, x in enumerate(reps):
        if not coprime(x, primes):
            continue
        u, x0 = primary_assoc(x)
        z = cx(x)
        tab[i] = xi_prim(x0) * (z / abs(z)) ** (-e)
    ch = Char(m, primes, e, tab)
    if check:     # well-definedness modulo m and the unit rule phi(u) = u^{-e}
        import random
        rnd = random.Random(1)
        for _ in range(60):
            x = (rnd.randint(-40, 40), rnd.randint(-40, 40))
            if not coprime(x, primes): continue
            y = add(x, mul(m, (rnd.randint(-5, 5), rnd.randint(-5, 5))))
            u, x0 = primary_assoc(y); z = cx(y)
            val = xi_prim(x0) * (z / abs(z)) ** (-e)
            assert abs(val - ch(np.array([x[0]]), np.array([x[1]]))[0]) < 1e-9, "not a character mod m"
        for u in UNITS:
            assert abs(ch(np.array([u[0]]), np.array([u[1]]))[0] - cx(u) ** (-e)) < 1e-12
    return ch

def is_trivial_on_kernel(ch, mprime, pi):
    """Is ch trivial on {c coprime to m: c = 1 mod mprime}?  (c = 1 + mprime*t, t over residues mod pi)."""
    for t in residue_reps(pi):
        c = add((1, 0), mul(mprime, t))
        if not coprime(c, ch.primes):
            continue
        if abs(ch(np.array([c[0]]), np.array([c[1]]))[0] - 1) > 1e-9:
            return False
    return True

def primitive(ch):
    """Reduce to the conductor.  Returns a new Char modulo the conductor f (with prime list of f)."""
    m = ch.m; plist = list(ch.primes)          # plist: primes of m with multiplicity
    changed = True
    cur_m, cur_list = m, plist
    while changed:
        changed = False
        for i, p in enumerate(cur_list):
            mp_ = exact_div(cur_m, p)
            if is_trivial_on_kernel(ch, mp_, p):
                cur_m = mp_; cur_list = cur_list[:i] + cur_list[i + 1:]
                changed = True; break
    f = cur_m; fprimes = sorted(set(cur_list))
    reps = residue_reps(f)
    tab = np.zeros(len(reps), dtype=complex)
    # cofactor for lifting
    for i, x in enumerate(reps):
        if not coprime(x, fprimes):
            continue
        # lift x to a class mod m coprime to all primes of m
        val = None
        for t in residue_reps(exact_div(m, f)) if N(f) < N(m) else [(0, 0)]:
            y = add(x, mul(f, t))
            if coprime(y, ch.primes):
                val = ch(np.array([y[0]]), np.array([y[1]]))[0]; break
        assert val is not None
        tab[i] = val
    out = Char(f, fprimes, ch.e, tab)
    out.prime_list = cur_list
    return out

def gauss(ch):
    """gamma(phi) = N(f)^{-1/2} sum_{x mod f} phi(x) e(x/f)."""
    f = ch.m; s = 0
    for i, x in enumerate(residue_reps(f)):
        if ch.table[i] != 0:
            s += ch.table[i] * e_add(cx(x) / cx(f))
    return s / math.sqrt(N(f))

def elements(Nmax):
    """All nonzero a + b omega with N <= Nmax (numpy arrays A, B, and complex values)."""
    R = int(math.ceil(2 * math.sqrt(Nmax / 3.0))) + 2
    a = np.arange(-R, R + 1); A, B = np.meshgrid(a, a, indexing='ij'); A = A.ravel(); B = B.ravel()
    n = A * A - A * B + B * B
    keep = (n > 0) & (n <= Nmax)
    A = A[keep]; B = B[keep]
    Z = A + B * OM
    return A, B, Z, (A * A - A * B + B * B)
