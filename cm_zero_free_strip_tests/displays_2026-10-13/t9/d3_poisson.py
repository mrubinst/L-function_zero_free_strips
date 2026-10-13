"""d3: the proof of Lemma s5-reflection (CM guide, Appendix F), display by display.

Notation: Xi((a)) = phi(a)(a/|a|)^e, phi primitive mod f, phi(u) = u^{-e}; c(a) = phi(a)(a/|a|)^e;
g_X(z) = (z/|z|)^e W(|z|^2/X); hat h(w) = int_C h(z) e(-zw) d^2z, e(z) = exp(4 pi i Im z/sqrt3).
F.9  eq:d176      sum_ideals Xi(a) W(N a/X) = (1/6) sum_{a in O} phi(a) g_X(a); also c(ua) = c(a).
F.10 eq:d177      sum_{y in O} g_X(t + f y) = (2/sqrt3) N(f)^{-1} sum_w e(tw/f) ghat_X(w/f)   (t mod f, incl. (t,f) != 1)
F.11 eq:s5-step1  sum_a phi(a) g_X(a) = (2/sqrt3) gamma(phi) N(f)^{-1/2} sum_w conj(phi(w)) ghat_X(w/f);
                  and [QF, Lemma 17]: sum_{t mod f} phi(t) e(tw/f) = conj(phi(w)) sqrt(N f) gamma(phi) for every w mod f.
F.12 eq:d178      int_0^{2pi} e^{i e th} exp(-4 pi i r rho sin(th+beta)/sqrt3) dth = 2pi J_e(.) e^{-i e beta}
                  = 2 pi e J_1(.) (conj w/|w|)^e.
F.13 eq:s5-hankel ghat_X(w) = pi e X (conj w/|w|)^e int W(x) J_1(4 pi |w| sqrt(Xx/3)) dx, ghat_X(0) = 0, for X != 1.
F.14 eq:d179      the three identities (sympy, symbolic) and the second one on character values.
F.15 eq:d180      sum_a phi(a) g_X(a) = (2 pi e X/sqrt3) gamma N(f)^{-1/2} (f/|f|)^e sum_{w != 0} conj c(w) W#(N(w)/Y).
In F.10 and F.11 ghat_X is computed from its DEFINITION by 2-d polar quadrature (independent of F.13).
"""
import os, sys, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import special
import sympy as sp
from t9lib import (report, falsifier, summary, cx, mul, N, OM, SQ3, UNITS, chi_p, primes_upto, from_primary_function,
                   primitive, gauss, elements, residue_reps, e_add, LAM, TWO, alpha_el, cubic_mod2, sextant_generator)

P = primes_upto(40)
p7 = [p for p in P if N(p) == 7][0]; p13 = [p for p in P if N(p) == 13][0]; p19 = [p for p in P if N(p) == 19][0]
cases = [
    ("alpha", lambda a: alpha_el(a), (3, 0), [LAM, LAM], 1),
    ("conj alpha", lambda a: 1 / alpha_el(a), (3, 0), [LAM, LAM], -1),
    ("alpha chi_p7^5", lambda a: alpha_el(a) * chi_p(a, p7) ** 5, mul((3, 0), p7), [LAM, LAM, p7], 1),
    ("alpha^-1 chi_p13^2", lambda a: chi_p(a, p13) ** 2 / alpha_el(a), mul((3, 0), p13), [LAM, LAM, p13], -1),
    ("conj alpha (a/2)_3", lambda a: cubic_mod2(a) / alpha_el(a), (6, 0), [LAM, LAM, TWO], -1),
    ("alpha chi_p7 chi_p19^3", lambda a: alpha_el(a) * chi_p(a, p7) * chi_p(a, p19) ** 3, mul(mul((3, 0), p7), p19),
     [LAM, LAM, p7, p19], 1),
]
chars = []
for name, xi, m, pr, e in cases:
    ch = primitive(from_primary_function(xi, m, pr, e))
    chars.append((name, ch, e))
    print(f"{name}: e={e}, f={ch.m}, N(f)={N(ch.m)}")

def W12(x): return np.sqrt(x) * np.exp(-x)
def W12s(y): return np.sqrt(y) * np.exp(-y)
def W32(x): return x ** 1.5 * np.exp(-x)
def W32s(y): return np.sqrt(y) * np.exp(-y) * (2 - y)

def gX(Z, e, W, X):
    out = np.zeros(Z.shape, dtype=complex); nz = np.abs(Z) > 0
    out[nz] = (Z[nz] / np.abs(Z[nz])) ** e * W(np.abs(Z[nz]) ** 2 / X)
    return out

# ------------------------------------------------------------------ 2-d polar quadrature of the definition
rg, rw = np.polynomial.legendre.leggauss(220)
NT = 384
th = np.arange(NT) * 2 * np.pi / NT
def ghat_def(W, e, X, w, rmax):
    """int_C (z/|z|)^e W(|z|^2/X) e(-z w) d^2 z by Gauss-Legendre in r on [0, rmax] and the trapezoid rule in theta."""
    r = 0.5 * rmax * (rg + 1); wr = 0.5 * rmax * rw
    Z = r[:, None] * np.exp(1j * th)[None, :]
    integrand = np.exp(1j * e * th)[None, :] * W(r ** 2 / X)[:, None] * np.exp(-4j * np.pi * (Z * w).imag / SQ3)
    return np.sum(integrand * (r * wr)[:, None]) * (2 * np.pi / NT)

def hankel_side(W, e, X, w, xr, kernel=special.j1, nx=600):
    x0, x1 = xr
    xg, xw = np.polynomial.legendre.leggauss(nx)
    x = 0.5 * (x1 - x0) * (xg + 1) + x0; wx = 0.5 * (x1 - x0) * xw
    return np.pi * e * X * (np.conj(w) / abs(w)) ** e * np.sum(W(x) * kernel(4 * np.pi * abs(w) * np.sqrt(X * x / 3)) * wx)

# ================================================================== F.12 (eq:d178)
rng = np.random.default_rng(1)
NT2 = 1024; th2 = np.arange(NT2) * 2 * np.pi / NT2
err = 0; f1 = 0; f2 = 0
for _ in range(60):
    r, rho, beta = rng.uniform(0.05, 4), rng.uniform(0.05, 4), rng.uniform(-np.pi, np.pi)
    for e in (1, -1):
        c = 4 * np.pi * r * rho / SQ3
        lhs = np.sum(np.exp(1j * e * th2) * np.exp(-1j * c * np.sin(th2 + beta))) * 2 * np.pi / NT2
        mid = 2 * np.pi * special.jv(e, c) * np.exp(-1j * e * beta)
        wbar = np.exp(-1j * beta)
        right = 2 * np.pi * e * special.j1(c) * wbar ** e
        err = max(err, abs(lhs - mid), abs(lhs - right))
        f1 = max(f1, abs(lhs - 2 * np.pi * special.j1(c) * wbar ** e))          # without the factor e
        f2 = max(f2, abs(lhs - 2 * np.pi * e * special.j1(c) * np.exp(1j * e * beta)))   # unconjugated angle
report("F12", err, 1e-12, "Jacobi-Anger: both equalities of (F.12), 60 random (r, rho, beta), e = +-1")
falsifier("F12", "2 pi J_1 (conj w/|w|)^e without the factor e", f1, 1e-3)
falsifier("F12", "angle e^{+i e beta}", f2, 1e-3)

# ================================================================== F.13 (eq:s5-hankel) with X != 1
def bump(x):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x); m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out
SIG = 0.6
def lognormal(x):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x); m = x > 0
    out[m] = np.exp(-np.log(x[m]) ** 2 / (2 * SIG ** 2)); return out
profiles13 = [("bump[1,2]", bump, (1.0, 2.0), lambda X: math.sqrt(2 * X)),
              ("x^1/2 e^-x", W12, (0.0, 45.0), lambda X: math.sqrt(45 * X)),
              ("lognormal", lognormal, (0.0, 80.0), lambda X: math.sqrt(80 * X))]
ws13 = [0.3 + 0.1j, -0.7 + 0.45j, 1.1 - 0.9j, 0.05j, 1.6]
for pname, W, xr, rmax in profiles13:
    for e in (1, -1):
        err = 0; zero = 0; fX = 0; fsq = 0
        for X in (0.37, 1.0, 4.2):
            if pname == "bump[1,2]":
                # annular: integrate r over [sqrt X, sqrt 2X] only
                def gd(w):
                    r0, r1 = math.sqrt(X), math.sqrt(2 * X)
                    r = 0.5 * (r1 - r0) * rg + 0.5 * (r1 + r0); wr = 0.5 * (r1 - r0) * rw
                    Z = r[:, None] * np.exp(1j * th)[None, :]
                    I = np.exp(1j * e * th)[None, :] * W(r ** 2 / X)[:, None] * np.exp(-4j * np.pi * (Z * w).imag / SQ3)
                    return np.sum(I * (r * wr)[:, None]) * (2 * np.pi / NT)
            else:
                gd = lambda w: ghat_def(W, e, X, w, rmax(X))
            zero = max(zero, abs(gd(0.0)))
            for w in ws13:
                d = gd(w)
                h = hankel_side(W, e, X, w, xr)
                err = max(err, abs(d - h))
                fX = max(fX, abs(d - h / X))                                       # X factor dropped
                fsq = max(fsq, abs(d - np.pi * e * X * (np.conj(w) / abs(w)) ** e *
                                   hankel_side(W, e, 1.0, w, xr) / (np.pi * e * (np.conj(w) / abs(w)) ** e)))  # no sqrt X in J_1
        report(f"F13-{pname}-e{e:+d}", err, 1e-10, "ghat_X(w) = pi e X (conj w/|w|)^e int W J_1(4 pi |w| sqrt(Xx/3)) dx, X = 0.37, 1, 4.2")
        report(f"F13-{pname}-e{e:+d}-zero", zero, 1e-12, "ghat_X(0) = 0")
        falsifier(f"F13-{pname}-e{e:+d}", "prefactor without X", fX, 1e-4)
        falsifier(f"F13-{pname}-e{e:+d}", "J_1 argument without sqrt X", fsq, 1e-4)

# ================================================================== F.9 (eq:d176)
err = 0; ecu = 0; fal9 = 0
for name, ch, e in chars:
    for X, W, cut in [(0.8, W12, 60.0), (5.0, W12, 60.0), (30.0, bump, 2.0)]:
        A, B, Z, n = elements(int(cut * X) + 1)
        k = n <= cut * X; A, B, Z, n = A[k], B[k], Z[k], n[k]
        phi = ch(A, B)
        rhs = np.sum(phi * gX(Z.astype(complex), e, W, X)) / 6
        lhs = 0
        seen = set()
        for a_, b_ in zip(A, B):
            g = sextant_generator((int(a_), int(b_)))
            if g in seen: continue
            seen.add(g)
            z = cx(g)
            lhs += ch(np.array([g[0]]), np.array([g[1]]))[0] * (z / abs(z)) ** e * W(np.array([abs(z) ** 2 / X]))[0]
        err = max(err, abs(lhs - rhs))
        fal9 = max(fal9, abs(lhs - np.sum(phi * W(n / X)) / 6))
        # c(ua) = c(a)
        c = phi * (Z / np.abs(Z)) ** e
        for u in UNITS:
            Au, Bu = mul((1, 0), u)
            uA = A * u[0] - B * u[1]; uB = A * u[1] + B * u[0] - B * u[1]
            cu = ch(uA, uB) * ((uA + uB * OM) / np.abs(uA + uB * OM)) ** e
            ecu = max(ecu, np.abs(cu - c).max())
report("F9", err, 1e-12, "ideal sum (one generator per ideal) = (1/6) element sum, 6 characters, X = 0.8, 5 (x^1/2 e^-x), 30 (bump)")
report("F9-units", ecu, 1e-12, "c(ua) = c(a) for the six units")
falsifier("F9", "element sum of phi(a) W without the angular factor", fal9, 1e-3)

# ================================================================== F.10, F.11 (eq:d177, eq:s5-step1) and [QF, Lemma 17]
err10 = 0; err11 = 0; err17 = 0; f10a = 0; f10b = 0; f11 = 0; f17 = 0
for name, ch, e in chars:
    f = ch.m; fz = cx(f); Nf = N(f)
    if Nf > 50: continue
    gam = gauss(ch)
    reps = residue_reps(f)
    # Lemma 17 for every w mod f
    for w in reps:
        lhs = sum(ch.table[i] * e_add(cx(t) * cx(w) / fz) for i, t in enumerate(reps))
        phw = ch(np.array([w[0]]), np.array([w[1]]))[0]
        err17 = max(err17, abs(lhs - np.conj(phw) * math.sqrt(Nf) * gam))
    for X in (0.6, 1.5, 4.0):
        rmax = math.sqrt(45 * X)
        # dual frequencies w with N(w/f) X <= 3.5
        Wn = 3.5 * Nf / X
        Aw, Bw, Zw, nw = elements(int(Wn) + 1)
        k = nw <= Wn; Aw, Bw, Zw = Aw[k], Bw[k], Zw[k]
        Aw = np.concatenate([[0], Aw]); Bw = np.concatenate([[0], Bw]); Zw = np.concatenate([[0j], Zw])
        gh = np.array([ghat_def(W12, e, X, complex(z) / fz, rmax) for z in Zw])
        # F.11
        A, B, Z, n = elements(int(60 * X) + 1)
        k = n <= 60 * X; A, B, Z = A[k], B[k], Z[k]
        lhs11 = np.sum(ch(A, B) * gX(Z.astype(complex), e, W12, X))
        phw = ch(Aw, Bw)
        rhs11 = 2 / SQ3 * gam / math.sqrt(Nf) * np.sum(np.conj(phw) * gh)
        err11 = max(err11, abs(lhs11 - rhs11))
        f11 = max(f11, abs(lhs11 - 2 / SQ3 * gam / math.sqrt(Nf) * np.sum(phw * gh)))
        # F.10 for every t mod f (including t not prime to f)
        tmax = max(abs(cx(t)) for t in reps)
        R = int((math.sqrt(60 * X) + tmax) / math.sqrt(Nf) * 2 / SQ3) + 3     # covers |t + f y|^2 <= 60 X
        ya = np.arange(-R, R + 1); YA, YB = np.meshgrid(ya, ya, indexing='ij'); YA = YA.ravel(); YB = YB.ravel()
        for t in reps:
            zt = cx(t) + fz * (YA + YB * OM)
            keep = np.abs(zt) ** 2 <= 60 * X
            lhs10 = np.sum(gX(zt[keep], e, W12, X))
            ph = np.array([e_add(cx(t) * complex(z) / fz) for z in Zw])
            rhs10 = 2 / SQ3 / Nf * np.sum(ph * gh)
            err10 = max(err10, abs(lhs10 - rhs10))
            f10a = max(f10a, abs(lhs10 - 2 / SQ3 / Nf * np.sum(np.conj(ph) * gh)))       # e(-tw/f)
            f10b = max(f10b, abs(lhs10 - 1 / Nf * np.sum(ph * gh)))                        # 1 in place of 2/sqrt3
report("F10", err10, 1e-12, "Poisson in the class t + fO for every t mod f, X = 0.6, 1.5, 4, five characters with N(f) <= 39")
falsifier("F10", "e(-tw/f)", f10a, 1e-4)
falsifier("F10", "1 in place of 2/sqrt3", f10b, 1e-4)
report("F11", err11, 1e-12, "sum_a phi(a) g_X(a) = (2/sqrt3) gamma N(f)^{-1/2} sum_w conj(phi(w)) ghat_X(w/f)")
falsifier("F11", "phi(w) unconjugated", f11, 1e-4)
report("F11-QF-Lemma17", err17, 1e-11, "sum_t phi(t) e(tw/f) = conj(phi(w)) sqrt(N f) gamma(phi) for every w mod f")
# imprimitive character: Lemma 17 fails at w not prime to f
name, ch, e = chars[0]
mimp = mul((3, 0), p7)
chimp = from_primary_function(lambda a: alpha_el(a) * (0 if chi_p(a, p7) == 0 else 1), mimp, [LAM, LAM, p7], 1)
reps = residue_reps(mimp); gi = sum(chimp.table[i] * e_add(cx(t) / cx(mimp)) for i, t in enumerate(reps)) / math.sqrt(N(mimp))
for w in reps:
    lhs = sum(chimp.table[i] * e_add(cx(t) * cx(w) / cx(mimp)) for i, t in enumerate(reps))
    f17 = max(f17, abs(lhs - np.conj(chimp(np.array([w[0]]), np.array([w[1]]))[0]) * math.sqrt(N(mimp)) * gi))
falsifier("F11-QF-Lemma17", "imprimitive character (alpha induced to modulus 3 p7)", f17, 1e-3)

# ================================================================== F.15 (eq:d180)
err15 = 0; f15a = 0; f15b = 0
for name, ch, e in chars:
    f = ch.m; fz = cx(f); Nf = N(f); gam = gauss(ch); C = 3 * Nf / (4 * math.pi ** 2)
    for W, Ws in [(W12, W12s), (W32, W32s)]:
        for X in (0.6, 1.5, 4.0, 11.0):
            Y = C / X
            A, B, Z, n = elements(int(70 * X) + 1)
            k = n <= 70 * X; A, B, Z = A[k], B[k], Z[k]
            lhs = np.sum(ch(A, B) * gX(Z.astype(complex), e, W, X))
            Aw, Bw, Zw, nw = elements(int(70 * Y) + 1)
            k = nw <= 70 * Y; Aw, Bw, Zw, nw = Aw[k], Bw[k], Zw[k], nw[k]
            cbar = np.conj(ch(Aw, Bw) * (Zw / np.abs(Zw)) ** e)
            S = np.sum(cbar * Ws(nw / Y))
            pref = 2 * np.pi * e * X / SQ3 * gam / math.sqrt(Nf)
            ang = (fz / abs(fz)) ** e
            err15 = max(err15, abs(lhs - pref * ang * S))
            if abs(ang - 1) > 1e-6: f15a = max(f15a, abs(lhs - pref * S))
            f15b = max(f15b, abs(lhs - pref * ang * S * e))     # e^2 = 1: catches the missing e for e = -1
report("F15", err15, 1e-11, "(F.15) for 6 characters (N(f) up to 1197), x^{1/2}e^{-x} and x^{3/2}e^{-x}, X = 0.6..11")
falsifier("F15", "without (f/|f|)^e", f15a, 1e-4)
falsifier("F15", "without the factor e", f15b, 1e-4)

# ================================================================== F.14 (eq:d179), symbolic
rho, sig, b_, g_, X_, x_, Nf_ = sp.symbols('rho sigma beta gamma X x N_f', positive=True)
ok = True
for e in (1, -1):
    w = rho * sp.exp(sp.I * b_); fz = sig * sp.exp(sp.I * g_)
    lhs = (sp.conjugate(w / fz) / sp.Abs(w / fz)) ** e
    rhs = (sp.conjugate(w) / sp.Abs(w)) ** e * (fz / sp.Abs(fz)) ** e
    ok &= sp.simplify(sp.expand_complex(lhs - rhs)) == 0
    wrong = (sp.conjugate(w) / sp.Abs(w)) ** e * (fz / sp.Abs(fz)) ** (-e)
    fal14 = sp.simplify(sp.expand_complex(lhs - wrong)) != 0
Y_ = 3 * sig ** 2 / (4 * sp.pi ** 2 * X_)
third = 4 * sp.pi * rho / sig * sp.sqrt(X_ * x_ / 3) - 2 * sp.sqrt(x_ * rho ** 2 / Y_)
ok3 = sp.simplify(third) == 0
third_wrong = 4 * sp.pi * rho / sig * sp.sqrt(X_ * x_ / 3) - 2 * sp.sqrt(x_ * rho ** 2 / (3 * sig ** 2 * X_ / (4 * sp.pi ** 2)))
fal3 = sp.simplify(third_wrong) != 0
# second identity on character values: conj(phi(w)) (conj w/|w|)^e = conj(c(w))
err2 = 0
for name, ch, e in chars:
    A, B, Z, n = elements(400)
    ph = ch(A, B); cc = ph * (Z / np.abs(Z)) ** e
    err2 = max(err2, np.abs(np.conj(ph) * (np.conj(Z) / np.abs(Z)) ** e - np.conj(cc)).max())
report("F14-first", 0.0 if ok else 1.0, 0.0, "sympy: (conj(w/f)/|w/f|)^e = (conj w/|w|)^e (f/|f|)^e, e = +-1, w, f in polar form")
report("F14-third", 0.0 if ok3 else 1.0, 0.0, "sympy: 4 pi |w/f| sqrt(Xx/3) = 2 sqrt(x N(w)/Y) with Y = 3N(f)/(4 pi^2 X)")
report("F14-second", err2, 1e-14, "conj(phi(w)) (conj w/|w|)^e = conj(c(w)) on 6 characters, N(w) <= 400")
falsifier("F14-first", "(f/|f|)^{-e}", 1.0 if fal14 else 0.0, 0.5)
falsifier("F14-third", "Y = C_Xi X instead of C_Xi / X", 1.0 if fal3 else 0.0, 0.5)
summary()
