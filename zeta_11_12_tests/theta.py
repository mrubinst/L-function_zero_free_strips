"""Kubota's cubic theta function from Patterson's Fourier coefficients.

Normalization of Dunn-Radziwill [DR, Section 5.1, (5.6)], also used in OpenAI's 11/12 paper:

    theta(z, v) = (3^{5/2}/2) v^{2/3}
                  + sum_{0 != l in lambda^{-3} O} tau(l) v K_{1/3}(4 pi |l| v) ce(l z),

with ce(z) = exp(2 pi i (z + conj z)), (z, v) in C x R_{>0} (upper half-space), and tau(l) given by
Patterson's formulas [DR, (5.7)-(5.8)].  The function tau(m) below returns tau(lambda^{-3} m) for
m in O = Z[omega]: writing l = u lambda^k c d^3 with u a unit, c squarefree primary and d primary,
the coefficient is nonzero only for k = 3 nu - 4 with nu >= 1 (any unit u) or k = 3 nu - 3 with
nu >= 0 (u = +-1), and is then a conjugated cubic Gauss sum of c times |d|/|c| times a power of 3.

This module is used by test_sec7_8.py:
  * test 7.1 compares d_0(lambda^{-3} n b^3) = conj(tau(-lambda^{-3} n b^3)) with guide (7.1);
  * test 7.3 checks which indices lambda^{-3} x (x primary) carry nonzero coefficients;
  * test 7.5 checks the automorphy theta(gamma w) = theta(w) under z -> z + 1 and under the inversion
    (z, v) -> (-conj z, v)/(|z|^2 + v^2), which is an independent check of the transcribed coefficients.

[DR] A. Dunn and M. Radziwill, Bias in cubic Gauss sums: Patterson's conjecture, Ann. of Math. 200 (2024).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import *
import numpy as np
from scipy.special import kv
from collections import Counter

LAMc = cx(LAM)                                       # lambda = i sqrt(3) as a complex number

def ce(z):
    """ce(z) = exp(2 pi i (z + conj z)) = exp(4 pi i Re z)."""
    return cmath.exp(4j*math.pi*z.real)

def cubic_sym(x, c):
    """(x/c)_3 for c primary and squarefree, possibly divisible by -2."""
    r = 1 + 0j; m = c
    if divides((2, 0), m):
        m = div(m, (-2, 0))
        xr = reduce_mod(x, (2, 0))
        if xr == (0, 0):
            return 0
        for val, rep in [(1, (1, 0)), (W, (0, 1)), (W*W, (-1, -1))]:   # residues mod 2: 1, omega, omega^2
            if divides((2, 0), sub(x, rep)):
                r *= val; break
        else:
            raise ValueError
    if N(m) > 1:
        for pi in factor_primary(m):
            j = sym6_prime(x, pi)
            if j is None:
                return 0
            r *= cmath.exp(2j*math.pi*j/3)           # (x/pi)_3 = (x/pi)_6^2
    return r

_g = {}
def _cubic_factory(c):
    """Return a function x -> (x/c)_3, with the factorization of c computed once."""
    m = c; two = False
    if divides((2, 0), m):
        m = div(m, (-2, 0)); two = True
    fac = factor_primary(m) if N(m) > 1 else []
    def f(x):
        r = 1 + 0j
        if two:
            if divides((2, 0), x):
                return 0
            for val, rep in [(1, (1, 0)), (W, (0, 1)), (W*W, (-1, -1))]:
                if divides((2, 0), sub(x, rep)):
                    r *= val; break
        for pi in fac:
            j = sym6_prime(x, pi)
            if j is None:
                return 0
            r *= cmath.exp(2j*math.pi*j/3)
        return r
    return f

def gsum(r, c):
    """Twisted cubic Gauss sum g(r, c) = sum_{x mod c} (x/c)_3 ce(r x / c), brute force, cached."""
    key = (r, c)
    if key in _g:
        return _g[key]
    cc = cx(c); s = 0; f = _cubic_factory(c)
    for x in residues(c):
        k = f(x)
        if k:
            s += k*ce(cx(mul(r, x))/cc)
    _g[key] = s
    return s

def lam_val(m):
    """(k, m') with m = lambda^k m' and lambda not dividing m'."""
    k = 0
    while divides(LAM, m):
        m = div(m, LAM); k += 1
    return k, m

def is_sqfree_primary_cd3(m):
    """For m coprime to lambda, return (u, c, d) with m = u c d^3, u a unit, c squarefree primary,
    d primary; return None if some prime occurs to an exponent = 2 mod 3 (then the coefficient is 0)."""
    p = prim_assoc(m); u = div(m, p)
    c = (1, 0); d = (1, 0); rest = p
    if divides((2, 0), rest):
        e = 0
        while divides((2, 0), rest):
            rest = div(rest, (-2, 0)); e += 1
        q, r = divmod(e, 3)
        if r == 2:
            return None
        c = mul(c, pw((-2, 0), r)); d = mul(d, pw((-2, 0), q))
    if N(rest) > 1:
        for pi, e in Counter(factor_primary(rest)).items():
            q, r = divmod(e, 3)
            if r == 2:
                return None
            c = mul(c, pw(pi, r)); d = mul(d, pw(pi, q))
    return u, c, d

def tau(m):
    """Patterson's coefficient tau(l) at l = lambda^{-3} m [DR, (5.7)-(5.8)]."""
    if m == (0, 0):
        return 0
    k, mp = lam_val(m); kk = k - 3                   # exponent of lambda in l
    dec = is_sqfree_primary_cd3(mp)
    if dec is None:
        return 0
    u, c, d = dec
    absdc = abs(cx(d))/abs(cx(c))
    L2 = mul(LAM, LAM)
    if kk % 3 == 2:                                  # kk = 3 nu - 4, nu >= 1; all six units
        nu = (kk + 4)//3
        if nu < 1:
            return 0
        if u in [(1, 0), (-1, 0)]:
            return gsum(L2, c).conjugate()*absdc*3**(nu/2 + 2)
        if u in [(0, 1), (0, -1)]:
            return cmath.exp(-2j*math.pi/9)*gsum(mul((0, 1), L2), c).conjugate()*absdc*3**(nu/2 + 2)
        if u in [(-1, -1), (1, 1)]:
            return cmath.exp(2j*math.pi/9)*gsum(mul((-1, -1), L2), c).conjugate()*absdc*3**(nu/2 + 2)
    if kk % 3 == 0:                                  # kk = 3 nu - 3, nu >= 0; units +-1 only
        nu = (kk + 3)//3
        if nu < 0:
            return 0
        if u in [(1, 0), (-1, 0)]:
            return gsum((1, 0), c).conjugate()*absdc*3**(nu/2 + 5/2)
    return 0

class Theta:
    """theta(z, v) truncated to Fourier indices |l| <= R (the K-Bessel factor makes the tail negligible
    for the heights used in the tests)."""
    def __init__(self, R=7.0):
        L3 = cx(LAM)**-3
        B = int(R*abs(cx(LAM))**3) + 3
        self.mu = []; self.coef = []
        for a in range(-2*B, 2*B + 1):
            for b in range(-2*B, 2*B + 1):
                m = (a, b)
                if m == (0, 0):
                    continue
                z = L3*cx(m)
                if abs(z) > R:
                    continue
                t = tau(m)
                if t != 0:
                    self.mu.append(z); self.coef.append(t)
        self.mu = np.array(self.mu); self.coef = np.array(self.coef)
        self.sig0 = 3**2.5/2

    def __call__(self, z, v, dz=0, dzbar=0, const=True):
        """theta(z, v), or a derivative: d/dz ce(l z) = 2 pi i l ce(l z), d/dzbar gives 2 pi i conj(l)."""
        mu = self.mu
        ph = np.exp(4j*np.pi*(mu*z).real)
        fac = (2j*np.pi*mu)**dz*(2j*np.pi*np.conj(mu))**dzbar
        val = np.sum(self.coef*fac*v*kv(1/3, 4*np.pi*np.abs(mu)*v)*ph)
        if const and dz == 0 and dzbar == 0:
            val += self.sig0*v**(2/3)
        return val

def act(g, z, v):
    """Action of g = (a, b, c, d) in SL_2(Z[omega]) on upper half-space:
    (z, v) -> ( ((az+b) conj(cz+d) + a conj(c) v^2)/D , v/D ),  D = |cz+d|^2 + |c|^2 v^2."""
    a, b, c, d = [cx(x) if isinstance(x, tuple) else x for x in g]
    den = abs(c*z + d)**2 + abs(c)**2*v*v
    return (((a*z + b)*np.conj(c*z + d) + a*np.conj(c)*v*v)/den, v/den)
