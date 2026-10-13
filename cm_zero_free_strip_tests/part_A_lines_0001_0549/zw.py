"""Exact arithmetic in the Eisenstein integers Z[omega] and sextic residue symbols.

An element a + b*omega of Z[omega] (omega = e^{2 pi i/3}, so omega^2 = -1 - omega) is stored as the
pair of Python integers (a, b).  All ring operations are exact; complex numbers appear only when a
residue symbol or an additive character is evaluated.

Conventions (as in the guide, zeta_11_12_guide.pdf):
  * primary            n = 1 mod 3, i.e. a = 1 and b = 0 mod 3 (guide (2.2));
  * lambda             1 + 2*omega = i*sqrt(3), stored as LAM = (1, 2);
  * chi_n(x)           the sextic residue symbol (x/n)_6, extended multiplicatively in n (guide (3.1));
                       its value at a prime pi is the sixth root of unity congruent to
                       x^{(N(pi)-1)/6} mod pi, and 0 if pi | x;
  * e(z)               exp(4 pi i Im(z)/sqrt(3)), the additive character used in the guide (6.1);
  * alpha(n)           n/|n|.

The sixth roots of unity are indexed by j = 0..5 as zeta^j with zeta = 1 + omega = e^{i pi/3}; the
list UNITS gives the corresponding elements of Z[omega].
"""
import cmath, math

W = cmath.exp(2j*math.pi/3)                         # omega as a complex number


# ---------------------------------------------------------------- ring operations on pairs (a, b)
def mul(x, y):
    """Product of x = a + b*omega and y = c + d*omega (uses omega^2 = -1 - omega)."""
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)

def add(x, y):
    """Sum x + y."""
    return (x[0] + y[0], x[1] + y[1])

def sub(x, y):
    """Difference x - y."""
    return (x[0] - y[0], x[1] - y[1])

def neg(x):
    """Negative -x."""
    return (-x[0], -x[1])

def conj(x):
    """Complex conjugate: conj(a + b*omega) = (a - b) - b*omega."""
    a, b = x
    return (a - b, -b)

def N(x):
    """Norm N(a + b*omega) = |a + b*omega|^2 = a^2 - ab + b^2 (guide (2.1))."""
    a, b = x
    return a*a - a*b + b*b

def tr(x):
    """Trace x + conj(x) = 2a - b."""
    return 2*x[0] - x[1]

def cx(x):
    """The element x as a Python complex number."""
    return x[0] + x[1]*W

def divides(d, x):
    """True if d divides x in Z[omega] (d nonzero): x*conj(d)/N(d) must have integer coordinates."""
    n = N(d); q = mul(x, conj(d))
    return q[0] % n == 0 and q[1] % n == 0

def div(x, d):
    """Exact quotient x/d; raises AssertionError if d does not divide x."""
    n = N(d); q = mul(x, conj(d))
    assert q[0] % n == 0 and q[1] % n == 0, (x, d)
    return (q[0]//n, q[1]//n)

def pw(x, e):
    """Power x^e for an integer e >= 0."""
    r = (1, 0)
    for _ in range(e):
        r = mul(r, x)
    return r


# ---------------------------------------------------------------- units, lambda, primary elements
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]   # zeta^j, j = 0..5, zeta = 1 + omega
LAM = (1, 2)                                                     # lambda = 1 + 2 omega = i sqrt(3)

def primary(x):
    """True if x is primary, i.e. x = 1 mod 3 (coordinates a = 1, b = 0 mod 3)."""
    return x[0] % 3 == 1 and x[1] % 3 == 0

def prim_assoc(x):
    """The unique primary associate u*x (x coprime to 3)."""
    c = [mul(u, x) for u in UNITS if primary(mul(u, x))]
    assert len(c) == 1, x
    return c[0]

def isprime_int(n):
    """Trial-division primality test for a rational integer."""
    return n > 1 and all(n % k for k in range(2, int(n**.5) + 1))


# ---------------------------------------------------------------- reduction and powers modulo m
def reduce_mod(y, m):
    """A small representative of y modulo m (subtract the nearest multiple of m)."""
    n = N(m); q = mul(y, conj(m))
    qa = round(q[0]/n); qb = round(q[1]/n)
    return sub(y, mul((qa, qb), m))

def powmod(x, e, m):
    """x^e modulo m by repeated squaring (result reduced with reduce_mod)."""
    r = (1, 0); b = reduce_mod(x, m)
    while e:
        if e & 1:
            r = reduce_mod(mul(r, b), m)
        b = reduce_mod(mul(b, b), m); e >>= 1
    return r

def canon_mod(y, m):
    """Canonical representative of y modulo m: y - floor(y/m)*m with a coordinatewise floor.

    Two elements are congruent modulo m exactly when they have the same canonical representative,
    so the result can be used as a dictionary key for residue classes."""
    n = N(m); q = mul(y, conj(m))
    qa = q[0]//n; qb = q[1]//n
    return sub(y, mul((qa, qb), m))

def residues(m):
    """A complete residue system of Z[omega]/m (canonical representatives, N(m) of them)."""
    n = N(m); seen = set()
    for a in range(n):
        for b in range(n):
            seen.add(canon_mod((a, b), m))
            if len(seen) == n:
                return sorted(seen)
    assert len(seen) == n, (m, len(seen))
    return sorted(seen)


# ---------------------------------------------------------------- primes, symbols, Moebius
def primes_upto(B):
    """Primary primes pi coprime to 6 with N(pi) <= B, sorted by norm (slow; see fastgauss.primary_primes)."""
    out = []
    for p in range(5, B + 1):
        if not isprime_int(p):
            continue
        if p % 3 == 1:                       # split: two conjugate primes of norm p
            r = int(2*p**.5) + 2
            for a in range(-r, r + 1):
                for b in range(-r, r + 1):
                    if N((a, b)) == p and primary((a, b)):
                        out.append((a, b))
        elif p*p <= B:                       # inert: -p is the primary generator, norm p^2
            out.append((-p, 0))
    return sorted(set(out), key=lambda x: (N(x), x))

def sym6_prime(x, pi):
    """(x/pi)_6 for a primary prime pi coprime to 6, as the index j in 0..5 of zeta^j; None if pi | x.

    Euler's criterion: the symbol is the sixth root of unity congruent to x^{(N(pi)-1)/6} mod pi."""
    if divides(pi, x):
        return None
    v = powmod(x, (N(pi) - 1)//6, pi)
    for j, u in enumerate(UNITS):
        if divides(pi, sub(v, u)):
            return j
    raise ValueError((x, pi))

_pc = {}
def PRIMES_CACHE(B):
    """Primary primes coprime to 6 of norm <= B, from a cached list."""
    B0 = max([k for k in _pc if k >= B], default=None)
    if B0 is None:
        B0 = max(B, 1000); _pc[B0] = primes_upto(B0)
    return [p for p in _pc[B0] if N(p) <= B]

def factor_primary(n, plist=None):
    """Factor a primary n coprime to 6 into primary primes (list with multiplicity), by trial division."""
    out = []; m = n
    for pi in PRIMES_CACHE(N(n)):
        while N(m) > 1 and divides(pi, m):
            out.append(pi); m = div(m, pi)
        if N(m) == 1:
            break
    assert m == (1, 0), (n, m)
    return out

def chi(n, x):
    """chi_n(x) = (x/n)_6 as a complex number, for any primary n coprime to 6 (product over the primes of n)."""
    r = 1 + 0j
    for pi in factor_primary(n):
        j = sym6_prime(x, pi)
        if j is None:
            return 0
        r *= cmath.exp(1j*math.pi*j/3)
    return r

def mu(n):
    """Moebius function mu_K of a primary n coprime to 6."""
    f = factor_primary(n)
    return 0 if len(set(f)) < len(f) else (-1)**len(f)

def alpha(n):
    """alpha(n) = n/|n| (guide (6.6))."""
    return cx(n)/abs(cx(n))

def e_add(z):
    """Additive character e(z) = exp(4 pi i Im(z)/sqrt(3)) of the guide (6.1), for a complex number z."""
    return cmath.exp(4j*math.pi*z.imag/math.sqrt(3))
