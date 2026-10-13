"""t5d: the reflected weight V# of a compactly supported W (Lemma lem:b3-decay).

Displays (labels of cm-eleven-twelfths-guide.tex, numbers of the 2026-10-13 build):
  (A.7)  eq:d101   |(x d/dx)^j V#(x)| << ||W||_{C^J} min(x^{1/4}, x^{-A})        -- class A: SIZE CHECK only
  (A.9)  eq:d102   |(x d/dx)^j V#(x)| <= Y^{-eta}/(2 pi) int |eta+iu|^j |Vhat(-eta-iu)| |R(eta+iu)| du, Y = (2pi)^4 x/27,
                   for eta = -1/4 and eta = A                                      -- explicit inequality, checked
Test function: V(y) = sqrt(y) W(y) with W(y) = y^{-1/2} p(log y), p(s) = (s(L-s))^K/(L/2)^{2K}, L = log 2, K = 18, i.e.
W in C^{17} with support I = [1, 2].  The proof of Lemma lem:b3-decay uses only ||W||_{C^J}, J = ceil(4A)+j+2 <= 12
here, so this W is admissible for the (A, j) tested (A = 1, 2; j = 0, 1, 2).  Vhat(-t) = int_0^L p(s) e^{-ts} ds is
computed by Gauss-Legendre quadrature for |t| <= 70 and by the exact finite expansion
sum_{i=K}^{2K} p^(i)(0)(1 - (-1)^i e^{-tL})/t^{i+1} for |t| > 70 (the two agree at |t| = 70).
(x d/dx)^j V#(x) = (1/2 pi i) int_(eta) (-t)^j Vhat(-t) R(t) ((2pi)^4 x/27)^{-t} dt with R the CM gamma quotient of eq:weight,
computed on the four lines Re t = -1/4, 0, 1, 2 (trapezoid in u).
Checks: (i) the four lines give the same values at every x (to within roundoff of each integral): R is holomorphic in
Re t > -4/3, so no pole is crossed; (ii) (A.7) size check: |(x d/dx)^j V#(x)| / min(x^{1/4}, x^{-A}) for A = 1, 2, 3 is
bounded on the measurable range (|V#| above roundoff) and tends to 0 as x -> 0; the fitted small-x slope is 4/3 (first pole
of R at t = -4/3), so the exponent 1/4 is valid but not sharp; at large x the log-log slope over the last measurable decade
is steeper than -3, and beyond it |V#| is below 1e-16; (iii) (A.9): the inequality holds at every x, j, for eta = -1/4, 1, 2,
with the left side computed on the line Re t = 0.
Falsifier: the small-x exponent 2 in place of 1/4 (the ratio grows by a factor > 100 from x = 1e-4 to 1e-8).
Usage: python3 -B t5d_decay.py      (about 30 s)
"""
import os, sys, time, math
t0 = time.time()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t5lib import report, falsifier, finish
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.special import loggamma

K = 18; L = math.log(2.0); T1 = 70.0
# p(s) = (s (L - s))^K / (L/2)^(2K) as a power series in s
base = np.array([0.0, L, -1.0])                     # s(L - s) = L s - s^2
pc = np.array([1.0])
for _ in range(K): pc = P.polymul(pc, base)
pc = pc/(L/2)**(2*K)
fact = np.array([math.factorial(i) for i in range(2*K + 1)], dtype=float)
p0der = pc*fact                                      # p^(i)(0) = i! * coefficient
xg, wg = np.polynomial.legendre.leggauss(400)
sg = (xg + 1)*L/2; wg = wg*L/2
psg = ((sg*(L - sg))/(L/2)**2)**K                  # p at the nodes, evaluated in product form (no cancellation)
def Vhat_minus(t):
    """Vhat(-t) = int_0^L p(s) e^{-t s} ds (vector t)."""
    t = np.asarray(t, complex); out = np.empty(t.shape, complex)
    small = np.abs(t) <= T1
    ts = t[small]
    for i0 in range(0, len(ts), 2000):
        out_idx = np.nonzero(small)[0][i0:i0 + 2000]
        out[out_idx] = np.exp(-np.outer(ts[i0:i0 + 2000], sg)) @ (psg*wg)
    tl = t[~small]
    if len(tl):
        acc = np.zeros(len(tl), complex); eL = np.exp(-tl*L)
        for i in range(K, 2*K + 1):
            acc += p0der[i]*(1 - (-1)**i*eL)/tl**(i + 1)
        out[~small] = acc
    return out
# consistency of the two evaluations at |t| = 300
tt = np.array([0.3 + 70j, -0.25 + 70j, 2 - 70j, 70 + 0j])
quad = np.exp(-np.outer(tt, sg)) @ (psg*wg)
acc = np.zeros(len(tt), complex)
for i in range(K, 2*K + 1): acc += p0der[i]*(1 - (-1)**i*np.exp(-tt*L))/tt**(i + 1)
vh_check = np.max(np.abs(quad - acc))/np.max(np.abs(acc))

def R(t):
    return np.exp(loggamma(4/3 + t) + loggamma(5/3 + t) - loggamma(4/3 - t) - loggamma(5/3 - t))
CX = (2*np.pi)**4/27
ETAS = [-0.25, 0.0, 1.0, 2.0]
H = 0.01
LINES = {}
for eta in ETAS:
    U = 3000.0
    u = np.arange(-U, U + H/2, H); t = eta + 1j*u
    LINES[eta] = (t, Vhat_minus(t)*R(t))
print(f"   Vhat: quadrature vs expansion at |t| = 70: rel diff {vh_check:.1e}; lines ready {time.time()-t0:.1f}s", flush=True)
XS = 10.0**np.arange(-8, 8.01, 0.5)
def Vsharp_j(x, j, eta):
    t, F = LINES[eta]
    return np.sum((-t)**j*F*np.exp(-t*np.log(CX*x)))*H/(2*np.pi)
def bound(x, j, eta):
    t, F = LINES[eta]
    return (CX*x)**(-eta)/(2*np.pi)*np.sum(np.abs(t)**j*np.abs(F))*H
vals = {}
for j in (0, 1, 2):
    for x in XS:
        vals[(j, x)] = {eta: Vsharp_j(x, j, eta) for eta in ETAS}
# roundoff scale of a contour integral: 1e-15 times the integral of the absolute integrand
def noise(x, j, eta): return 1e-15*bound(x, j, eta)
# (i) contour independence (R holomorphic in Re t > -4/3): every line agrees with Re t = 0 to within roundoff
cbad = 0; cmax = 0
for (j, x), d in vals.items():
    for e in ETAS:
        tol = 100*(noise(x, j, e) + noise(x, j, 0.0))
        diff = abs(d[e] - d[0.0])
        cmax = max(cmax, diff/tol)
        cbad += diff > tol
# values: Re t = -1/4 for x <= 1 (the proof's line), Re t = 0 for x > 1 (smallest roundoff)
def best(j, x): return vals[(j, x)][-0.25] if x <= 1 else vals[(j, x)][0.0]
def meas(j, x): return abs(best(j, x)) > 1e3*noise(x, j, -0.25 if x <= 1 else 0.0)
rows = []; okA = True
for j in (0, 1, 2):
    xm = [x for x in XS if meas(j, x)]
    for A in (1, 2, 3):
        r = np.array([abs(best(j, x))/min(x**0.25, x**-A) for x in xm])
        okA &= r[0] < 1e-3*r.max()
        okA &= r[-1] < r.max()                               # the maximum is attained inside the measurable range
        rows.append(f"j={j},A={A}: max {r.max():.2e} at x={xm[int(np.argmax(r))]:.0e} (measurable x <= {max(xm):.0e}; "
                    f"ratio {r[0]:.1e} at x=1e-8, {r[-1]:.1e} at the last x)")
# small-x slope (expected 4/3: residue of the pole at t = -4/3)
sl = [(math.log(abs(best(j, 1e-7))) - math.log(abs(best(j, 1e-8))))/math.log(10) for j in (0, 1, 2)]
# large x: log-log slope over the last measurable decade is steeper than -3 (faster than the powers tested)
sl_inf = []
for j in (0, 1, 2):
    xm = [x for x in XS if meas(j, x) and x > 1]
    i2 = list(XS).index(xm[-1]); x2 = XS[i2]; x1 = XS[i2 - 2]          # one decade (grid step 1/2 decade)
    sl_inf.append((math.log(abs(best(j, x2))) - math.log(abs(best(j, x1))))/math.log(10))
okslope = all(abs(s_ - 4/3) < 0.01 for s_ in sl) and all(s_ < -3 for s_ in sl_inf)
report("t5d.A7", "eq:d101 (size check)", okA and okslope and cbad == 0,
       f"W = y^-1/2 (log y log(2/y))^18 (normalized, C^17, support [1,2]); x in [1e-8, 1e8]: lines Re t = -1/4, 0, 1, 2 agree "
       f"within roundoff ({cbad} disagreements, max diff/(100*roundoff) {cmax:.1e}); |(x d/dx)^j V#|/min(x^1/4, x^-A): "
       + "; ".join(rows) + f"; slope at x -> 0: {', '.join(f'{s_:.4f}' for s_ in sl)} (= 4/3 >= 1/4); log-log slope over the "
       f"last measurable decade: {', '.join(f'{s_:.1f}' for s_ in sl_inf)} (< -3); beyond that |V#| is below roundoff 1e-16")
fal = (abs(best(0, 1e-8))/1e-8**2)/(abs(best(0, 1e-4))/1e-4**2)
falsifier("t5d.A7", "exponent 2 at x -> 0: ratio(1e-8)/ratio(1e-4)", fal, 100)
# (iii) (A.9): |(x d/dx)^j V#(x)| <= bound(eta) for eta = -1/4, 1, 2, left side from Re t = 0 (or -1/4 when eta = 0 is not used)
viol = 0; tight = {}
for j in (0, 1, 2):
    for x in XS:
        for eta in (-0.25, 1.0, 2.0):
            lhs = abs(vals[(j, x)][0.0]) if eta != 0.0 else abs(vals[(j, x)][-0.25])
            bd = bound(x, j, eta)
            if lhs > bd*(1 + 1e-9) + noise(x, j, 0.0): viol += 1
            key = (eta, 'x<=1' if x <= 1 else 'x>1')
            tight[key] = max(tight.get(key, 0), lhs/bd)
report("t5d.A9", "eq:d102", viol == 0,
       f"{3*len(XS)*3} cases (j = 0, 1, 2; 33 values of x in [1e-8, 1e8]; eta = -1/4, 1, 2), left side from the line Re t = 0: "
       f"{viol} violations; largest LHS/bound: " + ", ".join(f"eta={k[0]} {k[1]}: {v:.2e}" for k, v in sorted(tight.items())))
finish(t0)
