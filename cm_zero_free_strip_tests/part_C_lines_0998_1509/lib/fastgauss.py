"""Fast lists of primary primes and fast cubic Gauss sums.

The cubic Gauss sum of a squarefree primary c (possibly divisible by -2) is

    g(c) = sum_{x mod c} (x/c)_3 ce(x/c),        ce(z) = exp(2 pi i (z + conj z)),

the additive character used by Kubota, Patterson and Dunn-Radziwill.  It is computed from the prime
Gauss sums g(pi) (vectorised with numpy) and the twisted multiplicativity

    g(c1 c2) = (c2/c1)_3 (c1/c2)_3 g(c1) g(c2)        (c1, c2 coprime).

These are used for the theta coefficients (guide (7.1)) in test_sec7_8.py, test_sec7_completed_sanity.py
and test_sec7_reflection_j1.py.  Test 8.0 checks them against the brute-force Gauss sums gamma_2.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import mul, add, sub, conj, N, cx, divides, div, powmod, primary, UNITS, reduce_mod, pw
import numpy as np

W = cmath.exp(2j*math.pi/3)
CUBE = {(1, 0): 1, (0, 1): W, (-1, -1): W*W}      # cube roots of unity 1, omega, omega^2 in Z[omega]

def sieve(n):
    """Rational primes <= n (numpy array)."""
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0]

def prim_assoc(x):
    """The primary associate of x (x coprime to 3)."""
    for u in UNITS:
        y = mul(u, x)
        if primary(y):
            return y
    raise ValueError(x)

def primary_primes(B):
    """All primary primes of norm <= B, including -2 (norm 4), sorted by norm.

    p = 1 mod 3 contributes two conjugate primes of norm p (found from 4p = (2a-b)^2 + 3b^2);
    p = 2 mod 3 contributes -p, of norm p^2; the ramified prime lambda is excluded."""
    out = []
    for p in sieve(B):
        p = int(p)
        if p == 3:
            continue
        if p == 2:
            if 4 <= B:
                out.append((-2, 0))
            continue
        if p % 3 == 2:
            if p*p <= B:
                out.append((-p, 0))
            continue
        b = 1
        while 3*b*b <= 4*p:
            r = 4*p - 3*b*b; s = math.isqrt(r)
            if s*s == r and (s + b) % 2 == 0:
                a = (s + b)//2; pi = prim_assoc((a, b))
                out.append(pi); out.append(conj(pi)); break
            b += 1
    return sorted(set(out), key=lambda x: (N(x), x))

def cubic_prime(x, pi):
    """Cubic residue symbol (x/pi)_3 for a prime pi coprime to 3 (including -2); 0 if pi | x."""
    if divides(pi, x):
        return 0
    r = powmod(x, (N(pi) - 1)//3, pi)
    for u, val in CUBE.items():
        if divides(pi, sub(r, u)):
            return val
    raise ValueError((x, pi, r))

def powmod_vec(x, e, p):
    """Elementwise x^e mod p for a numpy integer array x and a rational prime p."""
    r = np.ones_like(x); b = x % p
    while e:
        if e & 1:
            r = (r*b) % p
        b = (b*b) % p; e >>= 1
    return r

def g_prime(pi):
    """The cubic Gauss sum g(pi) = sum_{x mod pi} (x/pi)_3 ce(x/pi) for a primary prime pi.

    For -2 the value is 2.  For an inert prime -q the residues mod q are enumerated as a + b*omega
    (0 <= a, b < q) and the symbol is computed by vectorised exponentiation in F_{q^2}.  For a split
    prime pi = a + b*omega of norm p, Z[omega]/pi = F_p via omega -> -a/b mod p."""
    if pi == (-2, 0):
        return 2.0 + 0j
    n = N(pi)
    if pi[1] == 0:                                   # inert prime -q
        q = -pi[0]
        A, Bb = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
        A = A.ravel().astype(np.int64); Bb = Bb.ravel().astype(np.int64)
        e = (q*q - 1)//3
        ra = np.ones_like(A); rb = np.zeros_like(A); ba = A.copy(); bb = Bb.copy()
        while e:
            if e & 1:
                ra, rb = (ra*ba - rb*bb) % q, (ra*bb + rb*ba - rb*bb) % q
            ba, bb = (ba*ba - bb*bb) % q, (2*ba*bb - bb*bb) % q
            e >>= 1
        chi = np.zeros(len(A), complex)
        chi[(ra == 1) & (rb == 0)] = 1; chi[(ra == 0) & (rb == 1)] = W; chi[(ra == q-1) & (rb == q-1)] = W*W
        nz = (A != 0) | (Bb != 0)
        assert np.all(np.abs(chi[nz]) > 0)
        ph = np.exp(-2j*np.pi*(2*A - Bb)/q)          # ce((a + b omega)/(-q)) = exp(-2 pi i (2a - b)/q)
        return np.sum(chi*ph)
    p = n; a, b = pi                                 # split prime of norm p
    x = np.arange(1, p, dtype=np.int64)
    r = powmod_vec(x, (p - 1)//3, p)
    w0 = (-a*pow(b, -1, p)) % p                      # image of omega in F_p
    chi = np.zeros(p - 1, complex); chi[r == 1] = 1; chi[r == w0] = W; chi[r == (w0*w0) % p] = W*W
    assert np.all(np.abs(chi) > 0)
    tr = 2*a - b                                     # ce(x/pi) = exp(2 pi i x tr(pi)/p) for rational x
    ph = np.exp(2j*np.pi*((x*tr) % p)/p)
    return np.sum(chi*ph)

_pair = {}
def pairsym(p1, p2):
    """(p1/p2)_3, cached."""
    k = (p1, p2)
    if k not in _pair:
        _pair[k] = cubic_prime(p1, p2)
    return _pair[k]

def squarefree_list(B, primes):
    """All squarefree primary c with N(c) <= B as pairs (c, [prime factors]); c = 1 is included."""
    out = [((1, 0), [])]
    def dfs(start, cur, fac, nrm):
        for i in range(start, len(primes)):
            p = primes[i]; nn = nrm*N(p)
            if nn > B:
                break
            c = mul(cur, p); f = fac + [p]; out.append((c, f)); dfs(i + 1, c, f, nn)
    dfs(0, (1, 0), [], 1)
    return out
