"""d4: Lemma s5-sharp (CM guide, Appendix F): (F.17) eq:d182 and a size check of (F.16) eq:d181.

F.17 eq:d182  int_0^oo J_1(2 sqrt(xy)) y^{s-1} dy = x^{-s} G(s),  G(s) = Gamma(s+1/2)/Gamma(3/2-s),  -1/2 < Re s < 3/4.
     Tested by direct integration in y (mpmath quadosc over the zeros j_{1,n}^2/(4x) of the kernel, 30 digits) for
     x = 0.3, 2.5 and s = -0.45+0.2i, 0.1, 0.5+1.1i, 0.7-0.4i (near both ends of the strip); also the classical
     Mellin transform int_0^oo t^{mu-1} J_1(t) dt = 2^{mu-1} Gamma((1+mu)/2)/Gamma((3-mu)/2), -1 < Re mu < 3/2.
     Falsifiers: Gamma(s)/Gamma(1-s) (OAI78's finite-order quotient); x^{+s}.
F.16 eq:d181  |(y d/dy)^j W#(y)| << min(y^{1/2}, y^{-A}) for an annular C_c^oo bump W on [1,2], j = 0, 1, 2:
     (i) the explicit small-y bound of the proof, |(y d/dy)^j W#(y)| <= y^{1/2} int |W_j(x)| x^{1/2} dx with
         W_j = (-d/dx x)^j W, on 81 points y in [1e-6, 1e4] (an inequality with explicit constant; checked exactly);
     (ii) the identity (y d/dy)^j W# = int W_j(x) J_1(2 sqrt(xy)) dx used in the proof;
     (iii) size check of the decay: sup of y^A |(y d/dy)^j W#| per decade, A = 1, 2, for y up to 1e6; the last decade
         is below the maximum (not a proof: the constants are not computable this way).
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
import sympy as sp
from scipy import special
from t9lib import report, falsifier, summary

mp.mp.dps = 30
def G(s): return mp.gamma(s + mp.mpf(1) / 2) / mp.gamma(mp.mpf(3) / 2 - s)
def Gfin(s): return mp.gamma(s) / mp.gamma(1 - s)

def lhs182(x, s):
    """Direct integral in y.  On [0, 1/x] subtract the leading term sqrt(x) y^{1/2} of J_1(2 sqrt(xy)) and integrate
    it exactly; on [1/x, oo) substitute y = v^2, where the kernel oscillates with period pi/sqrt(x) in v."""
    x = mp.mpf(x)
    f = lambda y: mp.besselj(1, 2 * mp.sqrt(x * y)) * y ** (s - 1)
    y0 = 1 / x
    half = mp.mpf(1) / 2
    a = mp.quad(lambda y: f(y) - mp.sqrt(x) * y ** (s - half), [0, y0]) + mp.sqrt(x) * y0 ** (s + half) / (s + half)
    g = lambda v: 2 * v * mp.besselj(1, 2 * mp.sqrt(x) * v) * v ** (2 * s - 2)
    b = mp.quadosc(g, [mp.sqrt(y0), mp.inf], period=mp.pi / mp.sqrt(x))
    return a + b

err = 0; f1 = 1e9; f2 = 1e9
for x in ["0.3", "2.5"]:
    for s in [mp.mpc("-0.45", "0.2"), mp.mpf("0.1"), mp.mpc("0.5", "1.1"), mp.mpc("0.7", "-0.4")]:
        v = lhs182(x, s)
        target = mp.mpf(x) ** (-s) * G(s)
        err = max(err, abs(v - target) / abs(target))
        f1 = min(f1, abs(v - mp.mpf(x) ** (-s) * Gfin(s)) / abs(target))
        f2 = min(f2, abs(v - mp.mpf(x) ** (s) * G(s)) / abs(target))
report("F17", float(err), 1e-15, "direct y-integral = x^{-s} Gamma(s+1/2)/Gamma(3/2-s), x = 0.3, 2.5, Re s from -0.45 to 0.7")
falsifier("F17", "Gamma(s)/Gamma(1-s)", float(f1), 1e-3)
falsifier("F17", "x^{+s}", float(f2), 1e-3)

err = 0
for mu in [mp.mpf("-0.9"), mp.mpc("0.3", "0.5"), mp.mpf("1.0"), mp.mpf("1.45")]:
    f = lambda t: t ** (mu - 1) * mp.besselj(1, t)
    v = (mp.quad(lambda t: (mp.besselj(1, t) - t / 2) * t ** (mu - 1), [0, 1]) + mp.mpf(1) / 2 / (mu + 1)
         + mp.quadosc(f, [1, mp.inf], period=2 * mp.pi))
    target = 2 ** (mu - 1) * mp.gamma((1 + mu) / 2) / mp.gamma((3 - mu) / 2)
    err = max(err, abs(v - target) / abs(target))
report("F17-DLMF-10.22.43", float(err), 1e-15, "int t^{mu-1} J_1(t) dt = 2^{mu-1} Gamma((1+mu)/2)/Gamma((3-mu)/2), mu = -0.9, 0.3+0.5i, 1, 1.45")

# ------------------------------------------------------------------ F.16 size check
# y d/dy = (t/2) d/dt for t = 2 sqrt(xy):  D^0 J_1 = J_1, D J_1 = (t/2) J_1', D^2 J_1 = (t/4) J_1' + (t^2/4) J_1''
fD = [lambda t: special.j1(t),
      lambda t: t / 2 * special.jvp(1, t, 1),
      lambda t: t / 4 * special.jvp(1, t, 1) + t ** 2 / 4 * special.jvp(1, t, 2)]
xs = sp.symbols('x', positive=True)
u = 2 * xs - 3
Wsym = sp.exp(-1 / (1 - u ** 2))
Wj = [Wsym]
for j in range(2):
    Wj.append(-sp.diff(xs * Wj[-1], xs))
fW = [sp.lambdify(xs, w, modules='numpy') for w in Wj]
xg, xw = special.roots_legendre(6000)
xx = 1.5 + 0.5 * xg; wx = 0.5 * xw
inner = np.abs(2 * xx - 3) < 0.995          # outside, exp(-1/(1-u^2)) < e^{-100}: set W_j = 0 (avoids 0/0)
Wvals = []
for j in range(3):
    v = np.zeros_like(xx)
    with np.errstate(all='ignore'):
        v[inner] = fW[j](xx[inner])
    assert np.all(np.isfinite(v)), "nonfinite W_j"
    Wvals.append(v)
const = [np.sum(np.abs(Wvals[j]) * np.sqrt(xx) * wx) for j in range(3)]
ys = np.logspace(-6, 4, 81)
viol = 0; ident = 0
for j in range(3):
    direct = np.array([np.sum(Wvals[0] * fD[j](2 * np.sqrt(xx * y)) * wx) for y in ys])
    via = np.array([np.sum(Wvals[j] * special.j1(2 * np.sqrt(xx * y)) * wx) for y in ys])
    ident = max(ident, np.abs(direct - via).max())
    assert np.all(np.isfinite(direct)) and np.isfinite(const[j])
    viol += int(np.sum(~(np.abs(direct) <= np.sqrt(ys) * const[j] * (1 + 1e-12))))
    print(f"   j={j}: sup_y |(y d/dy)^j W#| / (y^1/2 int|W_j| x^1/2) = {np.max(np.abs(direct) / (np.sqrt(ys) * const[j])):.4f}")
report("F16-ibp-identity", ident, 1e-11, "(y d/dy)^j W#(y) = int W_j(x) J_1(2 sqrt(xy)) dx, W_j = (-d/dx x)^j W, j = 0, 1, 2")
report("F16-small-y", viol, 0, "|(y d/dy)^j W#(y)| <= y^{1/2} int |W_j| x^{1/2} dx at 81 points y in [1e-6, 1e4], j = 0, 1, 2")
ok = True
for j in range(3):
    for A in (1, 2):
        sups = []
        for k in range(0, 6):
            yy = np.logspace(k, k + 1, 40)
            v = np.array([np.sum(Wvals[0] * fD[j](2 * np.sqrt(xx * y)) * wx) for y in yy])
            sups.append((yy ** A * np.abs(v)).max())
        ok &= sups[-1] < max(sups)
        print(f"   j={j}, A={A}: sup y^A |(y d/dy)^j W#| per decade [1,10],...,[1e5,1e6]: " + ", ".join(f"{s:.2e}" for s in sups))
report("F16-large-y", 0.0 if ok else 1.0, 0.0, "size check: y^A |(y d/dy)^j W#| (A = 1, 2) peaks at moderate y and decreases by y = 1e6")
summary()
