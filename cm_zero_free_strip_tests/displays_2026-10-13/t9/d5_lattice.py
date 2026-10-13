"""d5: the proof of Lemma s5-lattice (CM guide, Appendix F): (F.20) eq:d184, (F.21) eq:d185, size checks for
(F.22) eq:d186.  S = {lambda, 2}; theta in {1, (z/2)_3, (-1)^{(N z - 1)/2}} (finite-order, conductor primes in S).

F.20 eq:d184  L^0(X,t) = sum_{l outside S} alpha(l)^e theta(l) N(l)^{it} W(N(l)/X)   [ideals enumerated one by one,
              primary generator found by primary_assoc]
              = sum_{z in O} w(z) g_{X,t}(z),  g_{X,t}(z) = X^{it} h_t(z/sqrt X),  h_t(z) = (z/|z|)^e |z|^{2it} W(|z|^2),
              with w(z) = 1_{z = 1 mod 3} 1_{(z,2)=1} theta(z) stored as a TABLE on O/N_0, N_0 = 12 (|w| <= 1 checked).
F.21 eq:d185  L^0(X,t) = (2/sqrt3) N(N_0)^{-1} sum_w what(w) X^{1+it} hhat_t(sqrt X w/N_0),
              what(w) = sum_{t0 mod N_0} w(t0) e(t0 w/N_0).  Profile W(x) = x^{11/2} e^{-x} (h_t Schwartz enough: the
              dual terms decay like |w|^{-13}); hhat_t(xi) = pi e (conj xi/|xi|)^e (b/2) Gamma(mu+3/2) 1F1(mu+3/2;2;-b^2/4),
              mu = 11/2 + it, b = 4 pi |xi|/sqrt3, a closed form checked here against the 2-d DEFINITION of hhat_t
              by polar quadrature; also hhat_t(0) = 0 (the vanishing zero frequency).
F.22 eq:d186  size checks: (i) X sum_{w != 0} (1 + sqrt X |w|/|N_0|)^{-J0} is bounded for J0 = 3, X = 1 .. 1e6 (it is
              about 2 pi |N_0|^2/sqrt3 for X << |N_0|^2 and decays like X^{-1/2} beyond), and diverges for J0 = 2
              (falsifier: the hypothesis J0 >= 3 is needed);
              (ii) for an annular bump, sup_xi |hhat_t(xi)| (1+|xi|)^3 / (1+|t|)^3 stays bounded for t = 0 .. 64.
Falsifiers: w without the primary restriction; X^{-it} in g; X^{1-it} in (F.21); e(-t0 w/N_0) in what.
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from scipy import special
from t9lib import (report, falsifier, summary, cx, mul, add, N, OM, SQ3, UNITS, primary_assoc, is_primary, divides,
                   elements, residue_reps, canon_index, e_add, LAM, TWO, alpha_el, cubic_mod2, sextant_generator)

N0 = (12, 0)
reps = residue_reps(N0)                      # 144 classes, index = canon_index
def theta_fun(name, z):
    if name == "1": return 1.0
    if name == "(z/2)_3": return cubic_mod2(z)
    if name == "(-1)^((Nz-1)/2)": return 1.0 if N(z) % 4 == 1 else -1.0
THETAS = ["1", "(z/2)_3", "(-1)^((Nz-1)/2)"]

def w_value(name, z):
    if not is_primary(z) or divides(TWO, z): return 0.0
    return theta_fun(name, z)
def w_table(name):
    return np.array([w_value(name, r) for r in reps], dtype=complex)
def w_nonprimary(name, z):     # falsifier: drop the primary restriction (theta of the primary associate)
    if divides(LAM, z) or divides(TWO, z): return 0.0
    return theta_fun(name, primary_assoc(z)[1])

# well-definedness of w mod N_0 and |w| <= 1
rng = np.random.default_rng(3)
bad = 0; maxw = 0
for name in THETAS:
    tab = w_table(name); maxw = max(maxw, np.abs(tab).max())
    for _ in range(400):
        z = (int(rng.integers(-60, 61)), int(rng.integers(-60, 61)))
        y = (int(rng.integers(-9, 10)), int(rng.integers(-9, 10)))
        z2 = add(z, mul(N0, y))
        idx = canon_index(np.array([z2[0]]), np.array([z2[1]]), N0)[0]
        if abs(w_value(name, z) - tab[idx]) > 1e-15: bad += 1
report("F20-w-table", bad, 0, "w(z) = 1_{z=1 mod 3} 1_{(z,2)=1} theta(z) depends only on z mod N_0 = 12 (400 random z, z + 12y per theta)")
report("F20-|w|", max(0.0, maxw - 1), 1e-15, "|w| <= 1")

def bump(x):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x); m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out
K = 5
def Wsm(x): return np.asarray(x, dtype=float) ** (K + 0.5) * np.exp(-np.asarray(x, dtype=float))

# ================================================================== F.20
BI = 6000
RR = int(math.ceil(2 * math.sqrt(BI / 3.0))) + 2
ideals = set()
for a in range(-RR, RR + 1):
    for b in range(-RR, RR + 1):
        if 0 < N((a, b)) <= BI: ideals.add(sextant_generator((a, b)))
ideals = list(ideals)
err = 0; f1 = 0; f2 = 0
for name in THETAS:
    tab = w_table(name)
    for e in (1, -1):
        for t in (0.0, 2.3):
            for X, W in [(40.0, bump), (2500.0, bump), (60.0, Wsm)]:
                cut = min(BI, 2 * X if W is bump else 80 * X)
                lhs = 0
                for g in ideals:
                    if N(g) > cut or divides(LAM, g) or divides(TWO, g): continue        # outside S only
                    u, z0 = primary_assoc(g)
                    nn = N(z0)
                    lhs += alpha_el(z0) ** e * theta_fun(name, z0) * nn ** (1j * t) * W(np.array([nn / X]))[0]
                A, B, Z, n = elements(int(cut) + 1)
                k = n <= cut; A, B, Z, n = A[k], B[k], Z[k], n[k]
                wz = tab[canon_index(A, B, N0)]
                zs = Z / math.sqrt(X)
                h = (zs / np.abs(zs)) ** e * np.abs(zs) ** (2j * t) * W(np.abs(zs) ** 2)
                rhs = np.sum(wz * X ** (1j * t) * h)
                err = max(err, abs(lhs - rhs) / max(1.0, abs(lhs)))
                wnp = np.array([w_nonprimary(name, (int(a), int(b))) for a, b in zip(A, B)])
                f1 = max(f1, abs(lhs - np.sum(wnp * X ** (1j * t) * h)))
                if t != 0: f2 = max(f2, abs(lhs - np.sum(wz * X ** (-1j * t) * h)))
report("F20", err, 1e-11, "ideal-by-ideal L^0(X,t) = sum_z w(z) g_{X,t}(z): 3 theta, e = +-1, t = 0, 2.3, bump X = 40, 2500, x^{11/2}e^{-x} X = 60")
falsifier("F20", "w without the primary restriction", f1, 1e-3)
falsifier("F20", "X^{-it} in g_{X,t}", f2, 1e-3)

# ================================================================== F.21
mp.mp.dps = 20
def H_closed(rho, t):
    """int_0^oo x^{it} W(x) J_1(4 pi rho sqrt(x/3)) dx for W = x^{K+1/2} e^{-x}."""
    b = 4 * mp.pi * rho / mp.sqrt(3); mu = K + mp.mpf(1) / 2 + 1j * t
    return complex((b / 2) * mp.gamma(mu + mp.mpf(3) / 2) * mp.hyp1f1(mu + mp.mpf(3) / 2, 2, -b * b / 4))
def hhat_closed(xi, t, e):
    if xi == 0: return 0j
    return np.pi * e * (np.conj(xi) / abs(xi)) ** e * H_closed(abs(xi), t)
rg, rw = special.roots_legendre(400)
NT = 512; th = np.arange(NT) * 2 * np.pi / NT
def hhat_def(xi, t, e, rmax=10.0):
    r = 0.5 * rmax * (rg + 1); wr = 0.5 * rmax * rw
    Z = r[:, None] * np.exp(1j * th)[None, :]
    I = np.exp(1j * e * th)[None, :] * (r ** (2j * t) * Wsm(r ** 2))[:, None] * np.exp(-4j * np.pi * (Z * xi).imag / SQ3)
    return np.sum(I * (r * wr)[:, None]) * (2 * np.pi / NT)
errc = 0; zero = 0
for t in (0.0, 1.3):
    for e in (1, -1):
        zero = max(zero, abs(hhat_def(0.0, t, e)))
        for xi in [0.07 + 0.02j, -0.3 + 0.4j, 0.9 - 0.2j, 1.7j]:
            errc = max(errc, abs(hhat_def(xi, t, e) - hhat_closed(xi, t, e)))
report("F21-closed-form", errc, 1e-10, "closed form of hhat_t (Hankel + Kummer) = 2-d definition by polar quadrature, t = 0, 1.3, e = +-1")
report("F21-zero-frequency", zero, 1e-12, "hhat_t(0) = int (z/|z|)^e |z|^{2it} W(|z|^2) d^2z = 0")

err = 0; fa = 0; fb = 0
for name in THETAS[:2] + [THETAS[2]]:
    tab = w_table(name)
    for e in (1, -1):
        for t in (0.0, 1.3):
            for X in (20.0, 50.0):
                cut = 90 * X
                A, B, Z, n = elements(int(cut) + 1)
                k = n <= cut; A, B, Z, n = A[k], B[k], Z[k], n[k]
                wz = tab[canon_index(A, B, N0)]
                zs = Z / math.sqrt(X)
                lhs = np.sum(wz * X ** (1j * t) * (zs / np.abs(zs)) ** e * np.abs(zs) ** (2j * t) * Wsm(np.abs(zs) ** 2))
                rho_c = 14.0
                Wmax = (rho_c * 12 / math.sqrt(X)) ** 2
                Aw, Bw, Zw, nw = elements(int(Wmax) + 1)
                k = nw <= Wmax; Aw, Bw, Zw, nw = Aw[k], Bw[k], Zw[k], nw[k]
                # what(w) for each w
                tr = [cx(r) for r in reps]
                what = np.array([sum(tab[i] * e_add(tr[i] * complex(z) / 12) for i in range(len(reps)) if tab[i] != 0) for z in Zw])
                whatm = np.array([sum(tab[i] * e_add(-tr[i] * complex(z) / 12) for i in range(len(reps)) if tab[i] != 0) for z in Zw])
                # radial part once per norm
                cache = {}
                hh = np.empty(len(Zw), dtype=complex)
                for i, (z, nn) in enumerate(zip(Zw, nw)):
                    if nn not in cache: cache[nn] = H_closed(math.sqrt(X * nn) / 12, t)
                    xi = math.sqrt(X) * complex(z) / 12
                    hh[i] = np.pi * e * (np.conj(xi) / abs(xi)) ** e * cache[nn]
                pref = 2 / SQ3 / 144
                rhs = pref * np.sum(what * X ** (1 + 1j * t) * hh)           # w = 0 term vanishes
                err = max(err, abs(lhs - rhs) / max(1.0, abs(lhs)))
                if t != 0: fa = max(fa, abs(lhs - pref * np.sum(what * X ** (1 - 1j * t) * hh)))
                fb = max(fb, abs(lhs - pref * np.sum(whatm * X ** (1 + 1j * t) * hh)))
report("F21", err, 1e-10, "Poisson: L^0 = (2/sqrt3)/N(N_0) sum_w what(w) X^{1+it} hhat_t(sqrt X w/N_0); 3 theta, e = +-1, t = 0, 1.3, X = 20, 50")
falsifier("F21", "X^{1-it}", fa, 1e-4)
falsifier("F21", "e(-t0 w/N_0) in what", fb, 1e-4)

# ================================================================== F.22 size checks
c = 12.0
def latsum(X, J0, Rw):
    A, B, Z, n = elements(int(Rw ** 2) + 1)
    k = n <= Rw ** 2
    return X * np.sum((1 + math.sqrt(X) * np.abs(Z[k]) / c) ** (-J0))
vals = []
for X in [1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6]:
    Rw = 400.0
    s = latsum(X, 3, Rw)
    tail = (2 / SQ3) * 2 * math.pi * c ** 3 / (math.sqrt(X) * Rw)      # int_{|w|>R} (sqrt X |w|/c)^{-3} density 2/sqrt3
    vals.append(s + tail)
limit = 2 * math.pi * c ** 2 / SQ3
print("   X sum_{w!=0} (1+sqrt X|w|/12)^-3 for X = 1..1e6:", ", ".join(f"{v:.1f}" for v in vals), f"(2 pi |N_0|^2/sqrt3 = {limit:.1f})")
report("F22-latsum-J0=3", max(0.0, max(vals) / limit - 1.05), 0.0, "bounded by 1.05 x 2 pi |N_0|^2/sqrt3 uniformly for X = 1..1e6")
g2 = [latsum(1e4, 2, R) for R in (100.0, 400.0, 1600.0)]
print("   J0 = 2, X = 1e4, truncation |w| <= 100, 400, 1600:", ", ".join(f"{v:.1f}" for v in g2))
falsifier("F22-latsum", "J0 = 2 (the sum diverges logarithmically)", g2[-1] - g2[0], 100.0)

xg, xwts = special.roots_legendre(3000)
xx = 1.5 + 0.5 * xg; wx = 0.5 * xwts
rhos = np.linspace(0.0, 120.0, 4001)[1:]
ratios = []; ratios0 = []
for t in [0.0, 4.0, 16.0, 64.0]:
    Wt = bump(xx) * xx ** (1j * t)
    H = np.array([np.sum(Wt * special.j1(4 * np.pi * r * np.sqrt(xx / 3)) * wx) for r in rhos])
    mag = np.pi * np.abs(H)
    ratios.append((mag * (1 + rhos) ** 3).max() / (1 + t) ** 3)
    ratios0.append((mag * (1 + rhos) ** 3).max())
print("   sup |hhat_t| (1+|xi|)^3 / (1+|t|)^3 for t = 0, 4, 16, 64:", ", ".join(f"{v:.3g}" for v in ratios),
      "; without the (1+|t|)^3:", ", ".join(f"{v:.3g}" for v in ratios0))
report("F22-hhat-t", max(ratios) / ratios[0] - 1 if max(ratios) > ratios[0] else 0.0, 1.0,
       "size check: sup |hhat_t(xi)|(1+|xi|)^3/(1+|t|)^3 does not grow with t (bump, t <= 64)")
summary()
