"""Part D, test 04: the reflected kernels R_j and the gamma quotients (lines 1587-1598, 1605-1607, 1615).

Claims tested
  L1589-91  R_j(t) = prod_pm Gamma((j+1)/2+t+-1/6)/Gamma((j+1)/2-t+-1/6); R_1, R_2 are the table kernels
            (line 1224-1225; also eq. (weight), line 773, and OAI78's R(t), build/paper.tex)
  L1591-92  R_j is holomorphic for Re t > 1/6-(j+1)/2, which is left of -1/4 for every j >= 1
  L1593     |R_j(sigma+iu)| << (1+|u|)^{4 sigma} on -1/4 <= sigma <= A   (paper: checked j <= 5)
  L1594     the reflected weight obeys << min(x^{1/4}, x^{-A})            (partial numerical test)
  L1597     R_j(1) = (3j-4)(3j-2)(3j+2)(3j+4)/1296 ~ j^4/16
  L1606     exponent of |Gamma(1-s+k/2)/Gamma(s+k/2)| on Re s = sigma is 1-2 sigma, independent of k
  L1615     Gamma(s+k/2)/Gamma(1-s+k/2) is a Hankel kernel of order k (Mellin transform of J_k)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import report
import mpmath as mp
import sympy as sp
from fractions import Fraction as Fr

mp.mp.dps = 30
SIX = mp.mpf(1)/6

def R(j, t):
    a = mp.mpf(j + 1)/2
    return mp.gamma(a + t + SIX)*mp.gamma(a + t - SIX)/(mp.gamma(a - t + SIX)*mp.gamma(a - t - SIX))

# ---- D04.1 R_1, R_2 = table kernels (exact parameter sets)
def params(j): return sorted([Fr(j + 1, 2) + Fr(1, 6), Fr(j + 1, 2) - Fr(1, 6)])
ok = params(1) == [Fr(5, 6), Fr(7, 6)] and params(2) == [Fr(4, 3), Fr(5, 3)] and params(3) != [Fr(4, 3), Fr(5, 3)]
t0 = mp.mpc(0.3, 2.7)
num1 = mp.gamma(1 + t0 + SIX)*mp.gamma(1 + t0 - SIX)/(mp.gamma(1 - t0 + SIX)*mp.gamma(1 - t0 - SIX))
num2 = mp.gamma(mp.mpf(4)/3 + t0)*mp.gamma(mp.mpf(5)/3 + t0)/(mp.gamma(mp.mpf(4)/3 - t0)*mp.gamma(mp.mpf(5)/3 - t0))
ok &= abs(R(1, t0) - num1) < 1e-25 and abs(R(2, t0) - num2) < 1e-25 and abs(R(3, t0) - num2) > 1e-3
report("D04.1", 1591, "R_1, R_2", "R_1 = prod Gamma(1+t+-1/6)/Gamma(1-t+-1/6), R_2 = Gamma(4/3+t)Gamma(5/3+t)/(Gamma(4/3-t)Gamma(5/3-t))",
       ok, f"parameters j=1 {params(1)}, j=2 {params(2)}; falsifier R_3 != table CM kernel")

# ---- D04.2 holomorphy half-plane and its position
ok = True; details = []
for j in range(1, 41):
    a = Fr(j + 1, 2)
    poles = sorted([-n - a - Fr(1, 6) for n in range(3)] + [-n - a + Fr(1, 6) for n in range(3)])
    rightmost = poles[-1]
    ok &= rightmost == Fr(1, 6) - Fr(j + 1, 2) and rightmost < Fr(-1, 4)
    if j <= 6:
        # it is a genuine pole: denominator finite and nonzero there, |R| ~ C/eps
        e1, e2 = mp.mpf('1e-8'), mp.mpf('1e-10')
        tt = mp.mpf(rightmost.numerator)/rightmost.denominator
        r1, r2 = abs(R(j, tt + e1)), abs(R(j, tt + e2))
        ok &= 90 < r2/r1 < 110
        details.append(f"j={j}: pole {rightmost}")
report("D04.2", 1592, "holomorphy", "rightmost pole of R_j is t = 1/6-(j+1)/2 (a genuine simple pole), < -1/4 for j=1..40",
       ok, "; ".join(details[:3]) + "; falsifier: just left of the boundary R_j blows up like 1/eps")

# ---- D04.3 growth (1+|u|)^{4 sigma}.  The claim is a bound with an implied constant depending on (j, A); we test
#      that the ratio is finite on the grid and tends to 1 as u -> oo (Stirling), so its sup is a finite constant.
A = 3
sigmas = [-0.25, -0.1, 0, 0.25, 0.5, 1, 2, 3]
us = [0, 0.5, 1, 3, 10, 30, 100, 1e3, 1e4, 1e5]
sups = {}; stir = 0; fal_growth = []
for j in range(1, 7):
    for sg in sigmas:
        vals = [abs(R(j, mp.mpc(sg, u)))/(1 + u)**(4*sg) for u in us]
        sups[j] = max(sups.get(j, 0), max(vals))
        stir = max(stir, abs(abs(R(j, mp.mpc(sg, 1e7)))/mp.mpf(1e7)**(4*sg) - 1))
        f = [abs(R(j, mp.mpc(sg, u)))/(1 + u)**(4*sg - 0.2) for u in (1e2, 1e5)]
        fal_growth.append(f[1]/f[0])
report("D04.3", 1593, "kernel growth", f"|R_j(sigma+iu)|/(1+|u|)^(4 sigma) finite on j=1..6, sigma in [-1/4,{A}], u<=1e5, and -> 1 as u -> oo",
       all(mp.isfinite(v) for v in sups.values()) and stir < 1e-5 and min(fal_growth) > 3,
       "sup ratio per j: " + ", ".join(f"{j}:{mp.nstr(v, 3)}" for j, v in sups.items())
       + f"; |ratio(u=1e7)-1| <= {mp.nstr(stir, 2)}; falsifier exponent 4sigma-0.2: ratio grows x{mp.nstr(min(fal_growth), 3)}")

# ---- D04.4 V#(x) << min(x^{1/4}, x^{-A}) for a test V (partial)
s_w = mp.mpf('0.5'); K = (2*mp.pi)**4/27
def Vhat(t): return mp.sqrt(2*mp.pi)*s_w*mp.exp(s_w**2*t**2/2)     # Mellin transform of exp(-(log y)^2/(2 s_w^2))
def Vsharp(j, x):
    # contour Re t = sg; V# is independent of sg in the holomorphic strip (R_j holomorphic for Re t > 1/6-(j+1)/2,
    # V^ entire).  sg = -1/4 for x < 1 and sg = 3 for x >= 1 avoids cancellation at the two ends.
    sg = mp.mpf(-0.25) if x < 1 else mp.mpf(3)
    f = lambda u: Vhat(-(sg + 1j*u))*R(j, sg + 1j*u)*(K*x)**(-(sg + 1j*u))
    return mp.quad(f, [-30, -10, 0, 10, 30])/(2*mp.pi)
# contour independence (consistency of the shifted contours)
_c = [abs(mp.quad(lambda u: Vhat(-(sg + 1j*u))*R(2, sg + 1j*u)*(K*mp.mpf(0.7))**(-(sg + 1j*u)), [-30, 0, 30]))
      for sg in (-0.25, 0, 3)]
assert max(_c) - min(_c) < 1e-15*max(_c), _c
xs = [mp.mpf(10)**e for e in range(-8, 9)]
ok = True; ratios = []; leads = []; ends = []
for j in (1, 2, 3):
    vals = [abs(Vsharp(j, x)) for x in xs]
    r = max(v/min(x**0.25, x**-2) for v, x in zip(vals, xs))
    ratios.append(r)
    ends += [vals[0]/min(xs[0]**0.25, xs[0]**-2), vals[-1]/min(xs[-1]**0.25, xs[-1]**-2)]
    # sharpness / falsifier: V#(x) ~ C x^{(j+1)/2-1/6} as x -> 0, so the bound x^{(j+1)/2-1/6+0.1} fails
    e = mp.mpf(j + 1)/2 - SIX
    leads.append(vals[0]/xs[0]**(e + 0.1))
report("D04.4", 1594, "reflected weight", "|V#_j(x)|/min(x^{1/4}, x^{-2}) bounded on x in [1e-8,1e8] and decreasing to 0 at both ends, j=1,2,3, Gaussian test V (partial: one V, one A)",
       max(ratios) < 1e3 and max(ends) < max(ratios)/1000 and min(leads) > 1,
       f"max ratio={mp.nstr(max(ratios), 3)}, ratio at x=1e-8 and 1e8 <= {mp.nstr(max(ends), 2)}; falsifier x^((j+1)/2-1/6+0.1) violated at x=1e-8 (ratio {mp.nstr(min(leads), 3)})")

# ---- D04.5 R_j(1)
ok = True
for j in range(1, 61):
    exact = sp.gammasimp(sp.prod([sp.gamma(sp.Rational(j + 1, 2) + 1 + e)/sp.gamma(sp.Rational(j + 1, 2) - 1 + e)
                                  for e in (sp.Rational(1, 6), -sp.Rational(1, 6))]))
    claim = sp.Rational((3*j - 4)*(3*j - 2)*(3*j + 2)*(3*j + 4), 1296)
    ok &= sp.simplify(exact - claim) == 0 and abs(R(j, 1) - mp.mpf(claim.p)/claim.q) < 1e-20
fal = sp.Rational((3*5 - 4)*(3*5 - 2)*(3*5 + 2)*(3*5 + 4), 648) != sp.gammasimp(
    sp.gamma(sp.Rational(3) + 1 + sp.Rational(1, 6))*sp.gamma(sp.Rational(3) + 1 - sp.Rational(1, 6))
    /(sp.gamma(sp.Rational(3) - 1 + sp.Rational(1, 6))*sp.gamma(sp.Rational(3) - 1 - sp.Rational(1, 6))))
jj = sp.symbols('j', positive=True)
asym = sp.limit((3*jj - 4)*(3*jj - 2)*(3*jj + 2)*(3*jj + 4)/1296/(jj**4/16), jj, sp.oo)
report("D04.5", 1597, "R_j(1)", "R_j(1) = (3j-4)(3j-2)(3j+2)(3j+4)/1296 exactly for j=1..60; ratio to j^4/16 -> 1",
       ok and fal and asym == 1, f"limit={asym}; e.g. R_1(1)=-35/1296; falsifier /648 rejected")

# ---- D04.6 |Gamma(1-s+k/2)/Gamma(s+k/2)| ~ |t|^{1-2 sigma}, independent of k
worst = 0; fal = []
for k in range(0, 9):
    for sg in (-0.5, -0.1, 0.25, 0.5, 0.75, 1.1):
        s = mp.mpc(sg, 1e6)
        r = abs(mp.gamma(1 - s + mp.mpf(k)/2)/mp.gamma(s + mp.mpf(k)/2))/mp.mpf(1e6)**(1 - 2*sg)
        worst = max(worst, abs(r - 1))
        if k >= 1:
            fal.append(abs(abs(mp.gamma(1 - s + mp.mpf(k)/2)/mp.gamma(s + mp.mpf(k)/2))/mp.mpf(1e6)**(1 - 2*sg + k/2) - 1))
report("D04.6", 1606, "Hecke gamma quotient", "|Gamma(1-s+k/2)/Gamma(s+k/2)| / |t|^{1-2sigma} -> 1 at t=1e6 for k=0..8, six sigma",
       worst < 1e-4 and min(fal) > 0.9, f"max |ratio-1|={mp.nstr(worst, 3)}; falsifier exponent 1-2sigma+k/2: min |ratio-1|={mp.nstr(min(fal), 3)}")

# ---- D04.7 Hankel kernel of order k: int_0^oo J_k(x) x^{2s-1} dx = 2^{2s-1} Gamma(s+k/2)/Gamma(1-s+k/2)
mp.mp.dps = 20
worst = 0; fal = []
for k in range(0, 5):
    for s in (mp.mpc(0.1, 0.3), mp.mpc(0.2, -0.7)):
        if mp.re(s) <= -mp.mpf(k)/2: continue
        # v = t^20 on [0,1] removes the integrable endpoint singularity; quadosc on [1, oo)
        lhs = mp.quad(lambda t: mp.besselj(k, t**20)*(t**20)**(2*s - 1)*20*t**19, [0, 1]) \
              + mp.quadosc(lambda x: mp.besselj(k, x)*x**(2*s - 1), [1, mp.inf], omega=1)
        rhs = 2**(2*s - 1)*mp.gamma(s + mp.mpf(k)/2)/mp.gamma(1 - s + mp.mpf(k)/2)
        worst = max(worst, abs(lhs - rhs)/abs(rhs))
        lhs_f = mp.quad(lambda t: mp.besselj(k + 1, t**20)*(t**20)**(2*s - 1)*20*t**19, [0, 1]) \
                + mp.quadosc(lambda x: mp.besselj(k + 1, x)*x**(2*s - 1), [1, mp.inf], omega=1)
        fal.append(abs(lhs_f - rhs)/abs(rhs))
report("D04.7", 1615, "Hankel kernel", "Mellin transform of J_k at 2s equals 2^{2s-1}Gamma(s+k/2)/Gamma(1-s+k/2), k=0..4",
       worst < 1e-15 and min(fal) > 1e-2, f"max rel err={mp.nstr(worst, 3)}; falsifier J_(k+1): min rel err={mp.nstr(min(fal), 3)}")
