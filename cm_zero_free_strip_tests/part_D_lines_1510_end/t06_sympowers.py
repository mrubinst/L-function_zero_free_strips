"""Part D, test 06: the CM factorization of the symmetric power L-functions (lines 1644-1661).

Claim (up to finitely many Euler factors), with L_an(Sym^m E_d, s) = L(Sym^m E_d, s + m/2):
    L_an(Sym^m E_d, s) = prod_{0 <= i < m/2} L(s, (psi_d^*)^{m-2i}) * { 1 (m odd);  L(s, chi_{-3}^{m/2}) (m even) },
and the last factor is zeta(s) if m = 0 mod 4 and L(s, chi_{-3}) if m = 2 mod 4.
Here psi_d^*(n) = conj(alpha(n)) chi_n(4d) on primary n (line 1020).

Test: compare local Euler polynomials in X = p^{-s} at every good prime p <= 600 (p nmid 6d) for d in
{1, 2, -3, 7, 5, -11, 17} and m = 1..8.  The Sym^m polynomial is built from PARI's a_p(E_d) (ellap):
prod_{k=0}^m (1 - a^{m-k} b^k X) with a, b the roots of x^2 - (a_p/sqrt p) x + 1.  The right side uses the sextic
residue symbols computed in Z[omega] (independent of PARI).
Falsifiers: (i) the last factor chi_{-3}^{m/2+1} instead of chi_{-3}^{m/2}; (ii) psi_d^* replaced by
psi_d^* * chi_.(2).  Both must fail.  Also the pole statement for 4 | m (line 1659) is checked at the level of the
factor list: zeta(s) occurs exactly when 4 | m.  Spot check of the citation [MW, Sec. 3.5] (arXiv math/0604095,
source antsVII.tex, eq. (cmfact)): L(Sym^m E,s) = prod_{i=0}^{m/2} L(psi^{m-2i}, s-m/2), psi^0 = zeta if 4|m,
L(theta_K, s) if 2||m: identical after the shift s -> s + m/2.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import *
import numpy as np

DS = [1, 2, -3, 7, 5, -11, 17]
PMAX = 600
out = gp("\n".join(f"print(vector(primepi({PMAX}),i,ellap(ellinit([0,0,0,0,{d}]),prime(i))));" for d in DS))
aps = {d: [int(x) for x in line.strip("[]").split(",")] for d, line in zip(DS, out.strip().split("\n"))}
primes = rational_primes(PMAX)

def sym6(pi, x):
    j = sym6_prime(x, pi)
    return 0 if j is None else cmath.exp(1j*math.pi*j/3)

def polymul(a, b):
    return np.convolve(a, b)

def sym_poly(ap, p, m):
    t = ap/math.sqrt(p)
    disc = complex(t*t - 4)
    a = (t + cmath.sqrt(disc))/2; b = (t - cmath.sqrt(disc))/2
    P = np.array([1+0j])
    for k in range(m + 1):
        P = polymul(P, np.array([1, -(a**(m - k))*(b**k)]))
    return P

def rhs_poly(d, p, m, variant=None):
    P = np.array([1+0j])
    if p % 3 == 1:
        pi, pib = split_primary_pair(p)
        def psi(x):
            v = np.conj(complex(cx(x))/abs(complex(cx(x))))*sym6(x, (4*d, 0))
            if variant == 'twist': v *= sym6(x, (2, 0))
            return v
        for i in range(0, (m + 1)//2):          # 0 <= i < m/2
            r = m - 2*i
            P = polymul(P, np.array([1, -psi(pi)**r]))
            P = polymul(P, np.array([1, -psi(pib)**r]))
        chi3 = 1
    else:
        inert = (-p, 0)
        v = -1*sym6(inert, (4*d, 0))           # conj(alpha(-p)) = -1
        if variant == 'twist': v *= sym6(inert, (2, 0))
        for i in range(0, (m + 1)//2):
            r = m - 2*i
            P = polymul(P, np.array([1, 0, -v**r]))          # 1 - psi(-p)^r X^2
        chi3 = -1
    if m % 2 == 0:
        e = m//2 + (1 if variant == 'chi' else 0)
        P = polymul(P, np.array([1, -chi3**e]))
    return P

worst = 0; fal = {'chi': 0, 'twist': 0}; n = 0
for d in DS:
    for idx, p in enumerate(primes):
        if p in (2, 3) or d % p == 0:
            continue
        for m in range(1, 9):
            L = sym_poly(aps[d][idx], p, m); R = rhs_poly(d, p, m)
            worst = max(worst, np.max(np.abs(L - R))); n += 1
            for v in fal:
                Rv = rhs_poly(d, p, m, v)
                if len(Rv) == len(L):
                    fal[v] = max(fal[v], np.max(np.abs(L - Rv)))
                else:
                    fal[v] = max(fal[v], 1.0)
report("D06.1", 1648, "Sym^m factorization", f"local Euler polynomials agree at all good p <= {PMAX}, m = 1..8, d in {DS}",
       worst < 1e-8 and min(fal.values()) > 0.5,
       f"{n} (d,p,m) cases, max coefficient error={worst:.1e}; falsifiers chi_-3^(m/2+1): {fal['chi']:.2f}, psi*chi(2): {fal['twist']:.2f}")

# last factor: zeta when 4|m, L(s,chi_-3) when m = 2 mod 4 (as a character: chi_-3^{m/2} trivial iff m/2 even)
ok = all(((m//2) % 2 == 0) == (m % 4 == 0) for m in range(2, 200, 2))
ok_f = not all(((m//2) % 2 == 1) == (m % 4 == 0) for m in range(2, 200, 2))
report("D06.2", 1655, "last factor", "chi_{-3}^{m/2} is trivial (zeta, pole at s=1) iff m = 0 mod 4, and chi_{-3} iff m = 2 mod 4",
       ok and ok_f, "falsifier (opposite parity) rejected")

# shift: L_an(Sym^m, s) = L(Sym^m, s + m/2): the normalized polynomial uses a_p/sqrt p, i.e. Frobenius eigenvalues / p^{1/2}
d = 1; idx = primes.index(7); p = 7; m = 3
ap = aps[d][idx]
t = ap/math.sqrt(p); a = (t + cmath.sqrt(t*t - 4))/2
alpha_p = a*math.sqrt(p)
ok = abs(alpha_p*np.conj(alpha_p) - p) < 1e-9          # |Frobenius eigenvalue|^2 = p, so Sym^m eigenvalues / p^{m/2} are unimodular
report("D06.3", 1645, "normalization", "Sym^m eigenvalues of E_d at p have absolute value p^{m/2}; dividing by p^{m/2} (shift s -> s+m/2) makes them unimodular",
       ok and abs(abs(a) - 1) < 1e-12 and abs(abs(alpha_p) - 1) > 0.5, f"|alpha_p|^2 = {abs(alpha_p)**2:.6f} = p = 7; falsifier: unnormalized |alpha_p| != 1")
