"""Part D, test 03: the direct and reflected Mellin transforms with j zbar-derivatives (lines 1576-1586, 1600).

Claims tested
  L1579  int_0^oo (2 pi i lbar)^j v K_{1/3}(4 pi |l| v) v^{2s+j-2} dv
           = (i^j/4) (2 pi)^{-2s} Gamma(s+j/2-1/6) Gamma(s+j/2+1/6) conj(alpha(l))^j N(l)^{-s}
  L1576  the integral converges absolutely exactly for Re s > 1/6 - j/2
  L1583  j = 1, 2 give the table entries Gamma(s+1/3)Gamma(s+2/3) and Gamma(s+5/6)Gamma(s+7/6) (line 1223)
  L1585  reflected side: the same substitution gives prod_pm Gamma((j+2)/2 - s +- 1/6), the cusp factor
         conj(alpha(c))^{2j} and the dual angular factor alpha(mu)^j
  L1600  the reflected argument N(mu) X / N(c)^2 does not depend on j
The reflected-side integral is the Mellin transform (same weight v^{2s+j-2}) of the image of one Fourier
mode under the chain rule of L1571:
  (-(c v)^{-2})^j (2 pi i mu)^j v' K_{1/3}(4 pi |mu| v') ce(mu z'),   v' = 1/(N(c) v),
whose closed form (derived by the substitution w = 1/(N(c) v)) is
  (-1)^j i^j/4 * conj(alpha(c))^{2j} alpha(mu)^j (2 pi)^{2s-2} N(c)^{1-2s} N(mu)^{s-1}
     * Gamma((j+2)/2 - s - 1/6) Gamma((j+2)/2 - s + 1/6);
only the Gamma factors, conj(alpha(c))^{2j} and alpha(mu)^j are claimed by the paper; the remaining
factors are checked here as well.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import report
import mpmath as mp
from fractions import Fraction as Fr
import sympy as sp

mp.mp.dps = 40
third = mp.mpf(1)/3

def direct_lhs(j, s, l):
    f = lambda v: (2j*mp.pi*mp.conj(l))**j*v*mp.besselk(third, 4*mp.pi*abs(l)*v)*v**(2*s + j - 2)
    m = 30   # v = t^m near 0 removes the weak endpoint singularity v^{2 Re s + j - 4/3} (> -1)
    head = mp.quad(lambda t: f(t**m)*m*t**(m - 1), [0, mp.mpf(0.05)**(mp.mpf(1)/m)])
    return head + mp.quad(f, [0.05, 0.3, 1, 3, mp.inf])

def direct_rhs(j, s, l, variant=None):
    a = l/abs(l); q = abs(l)**2
    jj = j if variant != 'gshift' else j - 1
    ang = mp.conj(a)**j if variant != 'angle' else a**j
    pref = (1j**j/4) if variant != 'phase' else (1j**(j + 1)/4)
    return pref*(2*mp.pi)**(-2*s)*mp.gamma(s + mp.mpf(jj)/2 - mp.mpf(1)/6)*mp.gamma(s + mp.mpf(jj)/2 + mp.mpf(1)/6)*ang*q**(-s)

tests = [(mp.mpc(0.7, -1.3), mp.mpc(1.4, 0.8)), (mp.mpc(-0.4, 0.25), mp.mpc(0.35, -1.1)),
         (mp.mpc(2.1, 0.6), mp.mpc(-0.3, 0.5))]
worst = 0; fal = {'gshift': [], 'angle': [], 'phase': []}
for j in range(1, 7):
    for l, s in tests:
        if mp.re(s) <= mp.mpf(1)/6 - mp.mpf(j)/2:
            continue
        L = direct_lhs(j, s, l); R = direct_rhs(j, s, l)
        worst = max(worst, abs(L - R)/abs(R))
        for v in fal:
            fal[v].append(abs(L - direct_rhs(j, s, l, v))/abs(R))
falmin = {v: min(x) for v, x in fal.items()}
report("D03.1", 1579, "direct Mellin", "formula for j = 1..6, three (l, s) pairs each, 40-digit quadrature (paper: j <= 5, err < 1e-25)",
       worst < 1e-25 and all(x > 1e-3 for x in falmin.values()),
       f"max rel err={mp.nstr(worst, 3)}; falsifiers (Gamma(s+(j-1)/2+-1/6), alpha^j, i^(j+1)) min rel err "
       + ", ".join(f"{k}={mp.nstr(v, 3)}" for k, v in falmin.items()))

# ---- L1576: convergence exactly for Re s > 1/6 - j/2.  Near v = 0, K_{1/3}(x) ~ (Gamma(1/3)/2)(2/x)^{1/3},
#      so the integrand is ~ C v^{2s+j-4/3}; integrable at 0 iff Re(2s+j-4/3) > -1.  At infinity K decays exponentially.
xs = sp.symbols('x', positive=True)
lead = sp.limit(sp.besselk(sp.Rational(1, 3), xs)*xs**sp.Rational(1, 3), xs, 0)
ok_lead = sp.simplify(lead - sp.gamma(sp.Rational(1, 3))*2**sp.Rational(-2, 3)) == 0
sj, jj = sp.symbols('s j', real=True)
bound = sp.solve(sp.Eq(2*sj + jj - sp.Rational(4, 3), -1), sj)[0]
ok_bound = sp.simplify(bound - (sp.Rational(1, 6) - jj/2)) == 0
# numerical illustration: int_eps^1 grows like eps^{-delta} just left of the boundary, converges just right of it
def piece(j, s, eps):
    return mp.quad(lambda v: v*mp.besselk(third, 4*mp.pi*v)*v**(2*s + j - 2), [eps, 1])
grow = []; conv = []
for j in (1, 2, 3):
    b = mp.mpf(1)/6 - mp.mpf(j)/2
    left = [abs(piece(j, b - 0.05, e)) for e in (1e-4, 1e-8)]
    right = [abs(piece(j, b + 0.25, e)) for e in (1e-8, 1e-12)]
    grow.append(left[1]/left[0]); conv.append(abs(right[1] - right[0])/right[1])
report("D03.2", 1576, "convergence", "integrand ~ v^{2s+j-4/3} at 0: absolute convergence iff Re s > 1/6 - j/2",
       ok_lead and ok_bound and min(grow) > 2 and max(conv) < 1e-2,
       f"K_1/3(x) x^1/3 -> Gamma(1/3)2^(-2/3): {ok_lead}; at Re s = bound-0.05 the partial integral grows x{mp.nstr(min(grow), 3)} "
       f"(falsifier side); at bound+0.25, eps 1e-8 -> 1e-12 changes it by {mp.nstr(max(conv), 2)}")

# ---- L1583: j = 1, 2 give the table entries
def shifts(j): return sorted([Fr(j, 2) - Fr(1, 6), Fr(j, 2) + Fr(1, 6)])
ok = shifts(1) == [Fr(1, 3), Fr(2, 3)] and shifts(2) == [Fr(5, 6), Fr(7, 6)] and shifts(3) != [Fr(5, 6), Fr(7, 6)]
report("D03.3", 1583, "table entries", "j=1: Gamma(s+1/3)Gamma(s+2/3); j=2: Gamma(s+5/6)Gamma(s+7/6) (table line 1223)", ok,
       f"shifts j=1 {shifts(1)}, j=2 {shifts(2)}; falsifier j=3 gives {shifts(3)}")

# ---- L1585: reflected side
def refl_lhs(j, s, mu, c):
    Nc = abs(c)**2
    f = lambda v: v**(2*s + j - 2)*(-(c*v)**-2)**j*(2j*mp.pi*mu)**j*(1/(Nc*v))*mp.besselk(third, 4*mp.pi*abs(mu)/(Nc*v))
    return mp.quad(f, [0, 0.01, 0.1, 1, 10, 100, mp.inf])

def refl_rhs(j, s, mu, c, variant=None):
    Nc = abs(c)**2; Nm = abs(mu)**2
    ac = c/abs(c); am = mu/abs(mu)
    g0 = mp.mpf(j + 2)/2 if variant != 'gamma' else mp.mpf(j + 1)/2
    cuspp = mp.conj(ac)**(2*j) if variant != 'cusp' else mp.conj(ac)**j
    dual = am**j if variant != 'dual' else mp.conj(am)**j
    return ((-1)**j*1j**j/4*cuspp*dual*(2*mp.pi)**(2*s - 2)*Nc**(1 - 2*s)*Nm**(s - 1)
            *mp.gamma(g0 - s - mp.mpf(1)/6)*mp.gamma(g0 - s + mp.mpf(1)/6))

cases = [(mp.mpc(0.6, -0.9), mp.mpc(1.3, 0.4), mp.mpc(0.2, 0.7)), (mp.mpc(-1.1, 0.3), mp.mpc(0.5, -1.7), mp.mpc(0.45, -1.3)),
         (mp.mpc(0.25, 1.5), mp.mpc(-2.0, 0.9), mp.mpc(-0.4, 0.3))]
worst = 0; fl = {'gamma': [], 'cusp': [], 'dual': []}
for j in range(1, 5):
    for mu, c, s in cases:
        if mp.re(s) >= mp.mpf(j)/2 + mp.mpf(5)/6:
            continue
        L = refl_lhs(j, s, mu, c); R = refl_rhs(j, s, mu, c)
        worst = max(worst, abs(L - R)/abs(R))
        for v in fl:
            fl[v].append(abs(L - refl_rhs(j, s, mu, c, v))/abs(R))
flmin = {v: min(x) for v, x in fl.items()}
report("D03.4", 1585, "reflected side", "prod Gamma((j+2)/2 - s +- 1/6), conj(alpha(c))^{2j}, alpha(mu)^j (closed form incl. all constants), j=1..4",
       worst < 1e-20 and all(x > 1e-3 for x in flmin.values()),
       f"max rel err={mp.nstr(worst, 3)}; falsifiers (Gamma((j+1)/2-s..), conj(alpha(c))^j, conj(alpha(mu))^j) min rel err "
       + ", ".join(f"{k}={mp.nstr(v, 3)}" for k, v in flmin.items()))

# ---- L1600: j-independence of the reflected argument.  The ratio reflected/direct (per mode) is
#      (-1)^j conj(alpha(c))^{2j} alpha(mu)^j / alpha(l)bar^j * [(2pi)^4 N(mu) N(l)/N(c)^2]^s * (j-free s-free) * Gamma ratio.
s_, j_ = sp.symbols('s j'); Nc, Nm, Nl = sp.symbols('N_c N_mu N_l', positive=True)
refl_s = (2*sp.pi)**(2*s_ - 2)*Nc**(1 - 2*s_)*Nm**(s_ - 1)        # s-dependent prefactor of refl_rhs (j-free)
dir_s = (2*sp.pi)**(-2*s_)*Nl**(-s_)
dlog = sp.simplify(sp.diff(sp.log(refl_s/dir_s), s_))
ok = sp.simplify(dlog - sp.log((2*sp.pi)**4*Nm*Nl/Nc**2)) == 0 and not dlog.has(j_)
# numerical confirmation: with N(mu)/N(c)^2 fixed, refl_lhs * N(c) N(mu)/gamma-angular is the same for (mu, c) and (3 mu, sqrt3 c)
j = 2; s = mp.mpc(0.3, 0.6)
mu1, c1 = mp.mpc(0.6, -0.9), mp.mpc(1.3, 0.4)
mu2, c2 = 3*mu1, mp.sqrt(3)*c1
def normalized(mu, c):
    return refl_lhs(j, s, mu, c)/refl_rhs(j, s, mu, c)*(2*mp.pi)**(2*s - 2)*(abs(mu)**2/abs(c)**4)**s
r1, r2 = normalized(mu1, c1), normalized(mu2, c2)
num_ok = abs(r1 - r2)/abs(r1) < 1e-20
# falsifier: changing N(mu)/N(c)^2 changes the value
r3 = normalized(2*mu1, c1)
report("D03.5", 1600, "dual length", "reflected s-dependence is [(2pi)^4 N(mu) X/N(c)^2]^s, base independent of j",
       ok and num_ok and abs(r3 - r1)/abs(r1) > 1e-3,
       f"d/ds log(prefactor) = {dlog}; invariance under (mu,c)->(3mu,sqrt3 c): {num_ok}; falsifier (2mu, c) changes value")
