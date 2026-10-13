"""The theta coefficients of the appendix (lines 1700-1748), implemented literally from the paper's tables.

  g(r, c) = sum_{x mod c} (x/c)_3 ce(r x/c),   c primary squarefree (possibly even), (x/(-2))_3 = x mod 2.
  tau(l):   +-lambda^{3nu-4} c b^3 (nu>=1):  3^{nu/2+2} conj g(lambda^2, c) |b/c|
            +-omega lambda^{3nu-4} c b^3:    3^{nu/2+2} zeta_9^{-1} conj g(omega lambda^2, c) |b/c|
            +-omega^2 lambda^{3nu-4} c b^3:  3^{nu/2+2} zeta_9 conj g(omega^2 lambda^2, c) |b/c|
            +-lambda^{3nu-3} c b^3 (nu>=0):  3^{nu/2+5/2} conj g(1, c) |b/c|;   0 otherwise.
  tau_1(omega^j lambda^-4 c b^3), tau_2(-omega^j lambda^-4 c b^3): the 3x2 table of line 1724; 0 elsewhere.
  t_0 = tau, t_+(l) = omega tau_2(omega l) ce(l), t_-(l) = omega^2 tau_1(omega^2 l) ce(l), d_s(l) = conj(t_s(-l)).

An index l in lambda^{-4} O is passed as m = lambda^4 l = 9 l in O (a pair).
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import mul, add, sub, neg, conj, N, cx, divides, div, pw, primary, UNITS, LAM, residues, powmod, W
from fastgauss import primary_primes
from collections import Counter
import numpy as np

Z9 = cmath.exp(2j*math.pi/9)
OMEGA = (0, 1); OMEGA2 = (-1, -1)
LAMc = complex(cx(LAM))
_PR = primary_primes(20000)            # includes -2

def ce(z):
    return cmath.exp(4j*math.pi*z.real)

def cubic_sym_prime(x, pi):
    """(x/pi)_3 for a primary prime pi (pi = -2 allowed: then x^{(4-1)/3} = x mod 2)."""
    if divides(pi, x):
        return 0
    r = powmod(x, (N(pi) - 1)//3, pi)
    for u, val in [((1, 0), 1), ((0, 1), W), ((-1, -1), W*W)]:
        if divides(pi, sub(r, u)):
            return val
    raise ValueError((x, pi))

def factor_any(n):
    """Factor a primary n (coprime to 3, possibly even) into primary primes (including -2)."""
    out = []; m = n
    for p in _PR:
        if N(m) == 1: break
        while divides(p, m):
            out.append(p); m = div(m, p)
    if N(m) != 1:
        raise ValueError(("factor", n, m))
    return out

_gcache = {}; _chartab = {}
def cubic_char_table(c):
    if c not in _chartab:
        fac = factor_any(c) if N(c) > 1 else []
        tab = {}
        for x in residues(c):
            v = 1 + 0j
            for p in fac:
                s = cubic_sym_prime(x, p)
                if s == 0: v = 0; break
                v *= s
            tab[x] = v
        _chartab[c] = tab
    return _chartab[c]

def g(r, c):
    """Brute-force twisted cubic Gauss sum g(r, c)."""
    key = (r, c)
    if key not in _gcache:
        cc = complex(cx(c)); s = 0
        for x, v in cubic_char_table(c).items():
            if v != 0:
                s += v*ce(complex(cx(mul(r, x)))/cc)
        _gcache[key] = s
    return _gcache[key]

def prim_assoc(x):
    for u in UNITS:
        y = mul(u, x)
        if primary(y):
            return u, y
    raise ValueError(x)

def decompose(m):
    """m = lambda^k * u * c * b^3 with u a unit, c primary squarefree, b primary (both possibly even).
    Returns (k, u, c, b), or None if some prime exponent is 2 mod 3."""
    k = 0
    while divides(LAM, m):
        m = div(m, LAM); k += 1
    uinv, p = prim_assoc(m)               # p = uinv * m primary
    u = None
    for uu in UNITS:
        if mul(uu, uinv) == (1, 0): u = uu
    c = (1, 0); b = (1, 0)
    if N(p) > 1:
        for pi, e in Counter(factor_any(p)).items():
            q, r = divmod(e, 3)
            if r == 2: return None
            c = mul(c, pw(pi, r)); b = mul(b, pw(pi, q))
    return k, u, c, b

def _absbc(b, c):
    return abs(complex(cx(b)))/abs(complex(cx(c)))

L2 = mul(LAM, LAM)
def tau(m):
    """tau(l) at l = m/9 (m = lambda^4 l in O), from the table of line 1711."""
    if m == (0, 0): return 0
    dec = decompose(m)
    if dec is None: return 0
    k, u, c, b = dec
    kk = k - 4                                     # exponent of lambda in l
    f = _absbc(b, c)
    if kk % 3 == 2:                                # kk = 3 nu - 4
        nu = (kk + 4)//3
        if nu < 1: return 0
        A = 3**(nu/2 + 2)*f
        if u in [(1, 0), (-1, 0)]:   return A*np.conj(g(L2, c))
        if u in [OMEGA, (0, -1)]:    return A*Z9**-1*np.conj(g(mul(OMEGA, L2), c))
        if u in [OMEGA2, (1, 1)]:    return A*Z9*np.conj(g(mul(OMEGA2, L2), c))
    if kk % 3 == 0:                                # kk = 3 nu - 3
        nu = (kk + 3)//3
        if nu < 0: return 0
        if u in [(1, 0), (-1, 0)]:   return 3**(nu/2 + 5/2)*f*np.conj(g((1, 0), c))
    return 0

def tau1(m, variant=None):
    """tau_1(omega^j lambda^-4 c b^3), line 1726 (first column)."""
    if m == (0, 0): return 0
    dec = decompose(m)
    if dec is None: return 0
    k, u, c, b = dec
    if k != 0: return 0                           # l = m/9 has lambda-exponent -4 iff lambda does not divide m
    f = 9*_absbc(b, c)
    z9 = Z9 if variant != 'z9swap' else Z9**-1
    if u == (1, 0):  return f*W*np.conj(g(L2, c))
    if u == OMEGA:   return f*z9**-1*W*W*np.conj(g(mul(OMEGA, L2), c))
    if u == OMEGA2:  return f*z9*np.conj(g(mul(OMEGA2, L2), c))
    return 0

def tau2(m):
    """tau_2(-omega^j lambda^-4 c b^3), line 1726 (second column)."""
    if m == (0, 0): return 0
    dec = decompose(m)
    if dec is None: return 0
    k, u, c, b = dec
    if k != 0: return 0
    f = 9*_absbc(b, c)
    if u == (-1, 0):  return f*W*W*np.conj(g(L2, c))
    if u == (0, -1):  return f*Z9**-1*np.conj(g(mul(OMEGA, L2), c))
    if u == (1, 1):   return f*W*Z9*np.conj(g(mul(OMEGA2, L2), c))   # -omega^2 = 1 + omega
    return 0

def t_coeff(sig, m, variant=None):
    """t_sigma(l) at l = m/9."""
    l = complex(cx(m))/9
    if sig == '0': return tau(m)
    if sig == '+': return W*tau2(mul(OMEGA, m))*ce(l)
    if sig == '-': return W*W*tau1(mul(OMEGA2, m), variant)*ce(l)
    raise ValueError(sig)

def d_coeff(sig, m, variant=None):
    """d_sigma(l) = conj(t_sigma(-l))."""
    return np.conj(t_coeff(sig, neg(m), variant))
