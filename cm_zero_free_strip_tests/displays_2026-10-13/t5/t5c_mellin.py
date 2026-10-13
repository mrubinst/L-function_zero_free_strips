"""t5c: the Mellin steps of Appendix A (Steps 3 and 5), checked directly.

Displays tested (labels of cm-eleven-twelfths-guide.tex, numbers of the 2026-10-13 build):
  (A.10) eq:b3-bessel     int_0^oo v^{2s} K_{1/3}(4 pi|l|v) dv = 2^{2s-1} G(s+1/3) G(s+2/3) (4 pi|l|)^{-2s-1}, Re s > -1/3
                          (mpmath, 30 digits, including s close to -1/3 and complex s)
  (A.19)/(A.20) eq:d106, eq:b3-inversion  T(X;Psi) from the definition (eq:T, with W = V/sqrt(y)) equals
                          (1/2 pi i) int_(sigma_1) Vhat(w-1/2) T(w) X^{w-1/2} dw, T(w) the product of the n- and b-series
  (A.24) eq:b3-J          J(s) = int_0^oo f(v) v^{2s-1} dv, f = d_zbar^2 Theta_Psi at z = 0 computed from the twisted theta
                          series (QF coefficients), equals (3 pi/2)(27/(2 pi)^2)^s G(s+1/3) G(s+2/3) T(s-1/2), with T built
                          from alpha, gamma_2 (brute force) and Psi as in (A.19) -- same finite index set on both sides
  (A.32) eq:b3-dual       the dual Mellin transform: single mode (mpmath) and the truncated series G at three cusps
  (A.36) eq:d116          T(w) = J(w+1/2)/(...) is (A.24) at s = w+1/2 (sympy)
Falsifiers: 2^{2s} in (A.10); X^w in place of X^{w-1/2} and N(b)^{-3w} in (A.20); one derivative, T(s) in place of
T(s-1/2) and 3 pi in place of 3 pi/2 in (A.24); conj(alpha(l))^2 in (A.32).
Usage: python3 -B t5c_mellin.py      (about 1 minute)
"""
import os, sys, time, math, cmath
t0 = time.time()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t5lib import *
import cmlib
import mpmath as mp
from collections import Counter
import sympy as sp
from scipy.special import kv, gamma as cgamma

# ------------------------------------------------------------------------------------------------ (A.10)
mp.mp.dps = 30
worst = 0; fal = 0
for al in [mp.mpf('0.2'), mp.mpf('1.3')]:
    for s in [mp.mpf('-0.3'), mp.mpc('0.25', '1.5'), mp.mpf(2), mp.mpc('3.7', '-2')]:
        # substitution v = u^15 removes the integrable endpoint singularity v^(2s-1/3) (Re s close to -1/3)
        f = lambda v: v**(2*s)*mp.besselk(mp.mpf(1)/3, 4*mp.pi*al*v)
        lhs = mp.quad(lambda u: f(u**15)*15*u**14, [0, mp.mpf('0.5'), 1, mp.mpf('1.1'), mp.inf])
        rhs = 2**(2*s - 1)*mp.gamma(s + mp.mpf(1)/3)*mp.gamma(s + mp.mpf(2)/3)*(4*mp.pi*al)**(-2*s - 1)
        worst = max(worst, abs(lhs - rhs)/abs(rhs))
        fal = max(fal, abs(lhs - 2*rhs)/abs(rhs))
report("t5c.A10", "eq:b3-bessel", worst < 1e-20,
       f"8 cases (|l| = 0.2, 1.3; s = -0.3, 0.25+1.5i, 2, 3.7-2i), mpmath 30 digits: max rel err {mp.nstr(worst, 3)}; "
       f"same as [QF (5.21)] and [DLMF 10.43.19]")
falsifier("t5c.A10", "2^{2s} in place of 2^{2s-1}", float(fal), 1e-3)

# ------------------------------------------------------------------------------------------------ (A.32) single mode
worst = 0; fal = 0
for l in [mp.mpc('0.3', '0.7'), mp.mpc('-0.11', '0.02'), mp.mpc('1.5', '-0.4')]:
    al = l/abs(l); Nl = abs(l)**2
    for s in [mp.mpf('1.7'), mp.mpc('2.2', '0.7'), mp.mpc('4', '-3')]:
        lhs = mp.quad(lambda v: -4*mp.pi**2*Nl*al**2*v*mp.besselk(mp.mpf(1)/3, 4*mp.pi*abs(l)*v)*v**(2*s - 1), [0, 1, mp.inf])
        rhs = -(mp.pi/2)*(2*mp.pi)**(-2*s)*mp.gamma(s + mp.mpf(1)/3)*mp.gamma(s + mp.mpf(2)/3)*al**2*Nl**(mp.mpf(1)/2 - s)
        worst = max(worst, abs(lhs - rhs)/abs(rhs)); fal = max(fal, abs(lhs - rhs*mp.conj(al)**4)/abs(rhs))
w32a = float(worst); f32a = float(fal)

# ------------------------------------------------------------------------------------------------ twist and arithmetic data
PI7 = [p for p in fg.primary_primes(7) if N(p) == 7][0]
J7 = 2
def xi0(x): return psi4(x)*om3(x)
def Psi(n):
    if not prime_coprime6(n): return 0
    ch = sext_prime(n, PI7)
    return xi0(n)*(ch**J7 if ch != 0 else 0)
NMAX = 4000
xs = [(a, b) for a in range(-75, 76) for b in range(-75, 76) if prime_coprime6((a, b)) and N((a, b)) <= NMAX]
pairs = []
g2 = {}
for x in xs:
    _, ps = factor_prim(x)
    cnt = Counter(ps)
    if any(e % 3 == 2 for e in cnt.values()): continue
    n = (1, 0); b = (1, 0)
    for p, e in cnt.items():
        q, r = divmod(e, 3); n = mul(n, pw(p, r)); b = mul(b, pw(p, q))
    if n not in g2: g2[n] = cmlib.gamma_j_fast(2, n)              # gamma_2(n) = N(n)^-1/2 sum chi_n(x)^2 e(x/n)
    pairs.append((x, n, b))
print(f"   twist Psi = psi4 (omega/.)_3 chi_pi7^2; {len(pairs)} pairs (n, b) with N(n b^3) <= {NMAX}; "
      f"{len(g2)} gamma_2(n) by brute force, {time.time()-t0:.1f}s", flush=True)
def acoef(n, b):
    """the coefficient of (A.19): conj(alpha(n))^2 gamma_2(n) Psi(n) conj(alpha(b))^6 Psi(b)^3."""
    an = complex(cx(n))/abs(cx(n)); ab = complex(cx(b))/abs(cx(b))
    return np.conj(an)**2*g2[n]*Psi(n)*np.conj(ab)**6*Psi(b)**3
CO = np.array([acoef(n, b) for x, n, b in pairs])
NN = np.array([N(n) for x, n, b in pairs], float); NB = np.array([N(b) for x, n, b in pairs], float)
def Tw(w, variant=None):
    e_b = 3*w if variant == 'Nb3w' else 3*w - 0.5
    return np.sum(CO*NN**(-w)*NB**(-e_b))

# ------------------------------------------------------------------------------------------------ (A.24)
TT = ThetaTables(3*NMAX + 10, verbose=False)
tb = TT.table('0')
dmap = dict(zip(tb['m'], tb['d']))
# Theta_Psi coefficients at l = lambda^-3 x: d_0(l) phi(x) chi_pi(x)^2, phi(x) = (lambda/x)_3 xi_0(x)  (A.11)
L = []; Cth = []
for x, n, b in pairs:
    d0 = dmap.get(mul(LAM, x), 0)
    ch = sext_prime(x, PI7)
    Cth.append(d0*lam3_euler(x)*xi0(x)*(ch**J7 if ch != 0 else 0)); L.append(LAMc**-3*complex(cx(x)))
L = np.array(L); Cth = np.array(Cth); AL = np.abs(L)
def f_trunc(v, k=2):
    return np.sum(Cth*(2j*np.pi*np.conj(L))**k*v*kv(1/3, 4*np.pi*AL*v))
def J_quad(s, k=2, step=0.01):
    ts = np.arange(-12.0, 3.0, step)
    vals = np.array([f_trunc(math.exp(t), k)*math.exp(2*s.real*t)*np.exp(2j*s.imag*t) for t in ts])
    return np.sum(vals)*step
e24 = 0; q24 = 0; f24 = {'one derivative': 0, 'T(s)': 0, '3 pi': 0}
for s in [2.5 + 0j, 3 + 1.2j]:
    Jq = J_quad(s); Jq2 = J_quad(s, step=0.0125)
    q24 = max(q24, abs(Jq - Jq2)/abs(Jq))
    pref = (3*np.pi/2)*(27/(2*np.pi)**2)**s*cgamma(s + 1/3)*cgamma(s + 2/3)
    rhs = pref*Tw(s - 0.5)
    e24 = max(e24, abs(rhs - Jq)/abs(Jq))
    f24['one derivative'] = max(f24['one derivative'], abs(rhs - J_quad(s, k=1))/abs(rhs))
    f24['T(s)'] = max(f24['T(s)'], abs(pref*Tw(s) - Jq)/abs(Jq))
    f24['3 pi'] = max(f24['3 pi'], abs(2*rhs - Jq)/abs(Jq))
report("t5c.A24", "eq:b3-J", e24 < 1e-9 and q24 < 1e-11,
       f"{len(pairs)} modes l = lambda^-3 n b^3 (N(nb^3) <= {NMAX}); J(s) by trapezoid quadrature in log v of f(v) v^(2s-1) "
       f"(steps 0.01/0.0125 agree to {q24:.1e}) vs (3pi/2)(27/(2pi)^2)^s G(s+1/3)G(s+2/3)T(s-1/2), s = 2.5, 3+1.2i: "
       f"max rel err {e24:.1e}")
for k, v in f24.items():
    falsifier("t5c.A24", k, v, 1e-3)

# ------------------------------------------------------------------------------------------------ (A.19), (A.20)
s_w = 0.5
def V(y): return np.exp(-np.log(y)**2/(2*s_w*s_w))
def Vhat(sv): return math.sqrt(2*math.pi)*s_w*np.exp(s_w*s_w*sv*sv/2)       # int V(y) y^{s-1} dy
# product index set: n squarefree with N(n) <= 1000, b with N(b) <= 40 (both primary, coprime to 6): T(w) = (sum_n)(sum_b)
ns = sorted(set(n for x, n, b in pairs if N(n) <= 1000)); bs = sorted(set(b for x, n, b in pairs if N(b) <= 40))
an = np.array([np.conj(complex(cx(n))/abs(cx(n)))**2*g2[n]*Psi(n) for n in ns]); Nn = np.array([N(n) for n in ns], float)
ab = np.array([np.conj(complex(cx(b))/abs(cx(b)))**6*Psi(b)**3 for b in bs]); Nb = np.array([N(b) for b in bs], float)
def Tprod(w, variant=None):
    e_b = 3*w if variant == 'Nb3w' else 3*w - 0.5
    return np.sum(an*Nn**(-w))*np.sum(ab*Nb**(-e_b))
e20 = 0; f20 = {'X^w': 0, 'N(b)^-3w': 0}
for X in [60.0, 300.0]:
    # eq:T literally: X^{-1/2} sum_n sum_b conj(alpha(n))^2 gamma_2 Psi(n) conj(alpha(b))^6 Psi(b)^3 N(b)^{1/2} W(N(n)N(b)^3/X)
    Y = np.outer(Nn, Nb**3)/X
    Wv = V(Y)/np.sqrt(Y)
    Tdef = np.sum(np.outer(an, ab*np.sqrt(Nb))*Wv)/math.sqrt(X)
    for s1 in [1.5, 2.0]:
        us = np.arange(-30, 30, 0.002); w = s1 + 1j*us
        integ = np.array([Vhat(wi - 0.5)*Tprod(wi)*X**(wi - 0.5) for wi in w])
        Tm = np.sum(integ).real*0.002/(2*np.pi) + 1j*np.sum(integ).imag*0.002/(2*np.pi)
        e20 = max(e20, abs(Tm - Tdef)/abs(Tdef))
        if s1 == 2.0:
            f20['X^w'] = max(f20['X^w'], abs(np.sum(integ*X**0.5)*0.002/(2*np.pi) - Tdef)/abs(Tdef))
            iw = np.array([Vhat(wi - 0.5)*Tprod(wi, 'Nb3w')*X**(wi - 0.5) for wi in w])
            f20['N(b)^-3w'] = max(f20['N(b)^-3w'], abs(np.sum(iw)*0.002/(2*np.pi) - Tdef)/abs(Tdef))
report("t5c.A20", "eq:b3-inversion (with eq:d106)", e20 < 1e-9,
       f"{len(ns)} n x {len(bs)} b, log-normal V; T(X;Psi) from eq:T (W = V/sqrt(y)) vs the w-integral on Re w = 1.5, 2 with "
       f"T(w) = (sum_n)(sum_b) as in (A.19); X = 60, 300: max rel err {e20:.1e}")
for k, v in f20.items():
    falsifier("t5c.A20", k, v, 1e-3)

# ------------------------------------------------------------------------------------------------ (A.32) series at three cusps
TT2 = ThetaTables(20000, verbose=False)
e32 = 0; f32 = 0
for sig, (a, c) in [('0', ((1, 0), (-3, 0))), ('+', ((1, 0), (-1, -2))), ('-', ((0, 1), (1, 0)))]:
    cd = cusp_data(a, c); assert cd['sigma'] == sig
    t = TT2.table(sig); l = t['l']; d = t['d']
    keep = np.abs(l) <= 6.0
    l = l[keep]; d = d[keep]
    cc = complex(cx(cd['c'])); dl = complex(cx(cd['delta']))
    ph = ce(-dl*l/cc)
    def Gt(vp, conjang=False):
        ang = np.conj(l/np.abs(l))**2 if conjang else (l/np.abs(l))**2
        return -4*np.pi**2*np.sum(d*np.abs(l)**2*ang*vp*kv(1/3, 4*np.pi*np.abs(l)*vp)*ph)
    for s1 in [2.2 + 0.7j, 3.0 + 0j]:
        ts = np.arange(-12.0, 4.6, 0.01)
        Jq = np.sum([Gt(math.exp(tt))*np.exp(2*s1*tt) for tt in ts])*0.01
        Jf = np.sum([Gt(math.exp(tt), True)*np.exp(2*s1*tt) for tt in ts])*0.01
        rhs = -(np.pi/2)*(2*np.pi)**(-2*s1)*cgamma(s1 + 1/3)*cgamma(s1 + 2/3)*np.sum(d*(l/np.abs(l))**2*np.abs(l)**(-2*(s1 - 0.5))*ph)
        e32 = max(e32, abs(Jq - rhs)/abs(rhs)); f32 = max(f32, abs(Jf - rhs)/abs(rhs))
report("t5c.A32", "eq:b3-dual", w32a < 1e-20 and e32 < 1e-9,
       f"single mode (mpmath, 9 cases): max rel err {w32a:.1e}; truncated G (|l| <= 6) at cusps 1/(-3) (sigma 0), 1/(-lambda) (+), "
       f"omega/1 (-), s' = 2.2+0.7i, 3: quadrature vs series, max rel err {e32:.1e}")
falsifier("t5c.A32", "conj(alpha(l))^2 (single mode)", f32a, 1e-3)
falsifier("t5c.A32", "conj(alpha(l))^2 (series)", f32, 1e-3)

# ------------------------------------------------------------------------------------------------ (A.36)
w, s = sp.symbols('w s')
Tsym = sp.Function('T'); Jsym = sp.Function('J')
J_of_s = sp.Rational(3, 2)*sp.pi*(27/(2*sp.pi)**2)**s*sp.gamma(s + sp.Rational(1, 3))*sp.gamma(s + sp.Rational(2, 3))*Tsym(s - sp.Rational(1, 2))
expr = J_of_s.subs(s, w + sp.Rational(1, 2))/(sp.Rational(3, 2)*sp.pi*(27/(2*sp.pi)**2)**(w + sp.Rational(1, 2))
                                               *sp.gamma(w + sp.Rational(5, 6))*sp.gamma(w + sp.Rational(7, 6)))
ok36 = sp.simplify(expr - Tsym(w)) == 0
report("t5c.A36", "eq:d116", ok36, "sympy: (A.24) at s = w + 1/2 divided by (3pi/2)(27/(2pi)^2)^(w+1/2) G(w+5/6) G(w+7/6) is T(w)")
finish(t0)
