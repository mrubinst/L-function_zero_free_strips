"""c3: theta series, Poisson step, Mellin transform and the split identity of Lemma s2-fe (Appendix C of the CM guide):
displays eq:d148 (C.4), eq:d149 (C.5), eq:d150 (C.6), eq:d153 (C.9), eq:d154 (C.10), eq:d155 (C.11), eq:d156 (C.12),
eq:s2-split (C.13).

Theta_phi(v) = sum_{a in O} phi(a) conj(a) exp(-2 pi v N(a)/sqrt3) is always computed DIRECTLY from its definition
(eq:d147), never through its transformation law, with all terms of exp(...) >= e^{-46}.  phi is tabulated on O/f and
looked up by the injective key (x conj f mod Q).  Integrals over (0, oo) are computed in x = log v by the trapezoid rule
(the integrand is analytic and negligible at the cut points v_lo = 0.0864/Q, v_hi = 20, so the rule converges
spectrally; checked against the step 2h); integrals over (0, v0) and (v0, oo) by Gauss-Legendre in log v (checked
against more nodes).  Characters: t7lib.test_chars (eight primitive phi, Q = 7 ... 225).

Claims
  T7-C9a eq:d153 (Poisson step)  for every residue t mod f and v in {0.5, 1, 2}/sqrt Q:
         sum_{z in O} G(t + f z) = (2/(sqrt3 Q)) sum_{y in O} e(t y/f) Ghat(y/f),  G(w) = conj(w) e^{-2 pi v|w|^2/sqrt3},
         Ghat from eq:d152.  Falsifiers: Jacobian Q omitted; e(-t y/f) in place of e(t y/f).
  T7-C9b eq:d153 / eq:d148  first line = second line = Theta_phi(v) = -gamma/(f sqrt Q) v^-2 conj(Theta_phi(1/(Qv))),
         at 5 values of v.  Falsifiers: no conjugation, conj(gamma), sign +, conj(f).
  T7-C5  eq:d149 / eq:d154  for Re s > 1: int_0^oo Theta v^{s+1/2} dv/v (quadrature)
             = (2 pi/sqrt3)^{-s-1/2} Gamma(s+1/2) sum_{a != 0} phi(a) conj(a) N(a)^{-s-1/2}   (element sum, eq:d154)
             = 6 (2 pi/sqrt3)^{-s-1/2} Gamma(s+1/2) L(s, psi)                                   (ideal sum, eq:d149),
         L(s,psi) summed over ideals with the canonical generator x + y zeta (x > 0, y >= 0); N <= 2e5 (tail < 1e-7);
         s = 2.5+3i, 3, 3.5-2i; the same for conj(Theta) and conj(psi).
         Falsifiers: factor 6 omitted; Gamma(s) for Gamma(s+1/2); (2 pi)^{-s-1/2} for (2 pi/sqrt3)^{-s-1/2}.
  T7-C11 eq:d155  B(w) and B^-(w) by quadrature over (v0, oo) equal the incomplete-gamma series
         sum phi(a) conj(a) (c N)^{-(w+1/2)} Gamma(w+1/2, c N v0) used in tests_cm/phase2_2026-10-12/s2/t1_theta_fe.py.
  T7-C12 eq:d156  for Re s > 1 (s = 1.5+2i, 2.5-1i): int_0^{v0} Theta v^{s+1/2} dv/v
             = -gamma/(f sqrt Q) int_0^{v0} conj(Theta(1/(Qv))) v^{s-3/2} dv/v = -gamma/(f sqrt Q) Q^{3/2-s} B^-(1-s),
         v0 = Q^{-1/2}, each term by its own quadrature.  Falsifiers: Q^{s-3/2}; B^-(s) for B^-(1-s).
  T7-C13 eq:s2-split  Q^{s/2} I_psi(s) = Q^{s/2} B(s) + W Q^{(1-s)/2} B^-(1-s), and the conjugate identity
         Q^{s/2} I_conjpsi(s) = Q^{s/2} B^-(s) + conj(W) Q^{(1-s)/2} B(1-s), at s = 1.5+2i, 2.5-1i, and (I entire) at
         s = 0.3+5i, -0.4-2i; plus the symbolic identities -gamma/(f sqrt Q) = W/Q and
         (2 pi/sqrt3)^{-s-1/2} = 3^{1/4}(2 pi)^{-1/2} 3^{s/2} (2 pi)^{-s}  (the constant kappa).
         Falsifiers: W -> conj(W) (for non-real W), W -> -W.
  T7-C6  eq:d150  independent of the theta transformation: with I computed by direct quadrature for s in and
         around the critical strip, Q^{s/2} I_psi(s) = W Q^{(1-s)/2} I_conjpsi(1-s) (this is Lambda(s,psi) =
         W Lambda(1-s, conj psi) times kappa), W = -gamma sqrt Q/f, |W| = 1; s = 0.5+3i, 0.2-4i, 0.8+0.3i, -0.3+1.5i,
         1.4-2i.  Falsifiers: W -> -W, W -> conj(W) (non-real W), N = Q in place of the scaling (Q^{s/2} -> (3Q)^{s/2}
         on the left only).
"""
import os, sys, math, cmath, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t7lib import *
import numpy as np
import mpmath as mp
import sympy as sp

t0 = time.time()
c2 = 2 * math.pi / SQ3
CUT = 46.0
chars = test_chars()


def table(ch):
    """phi on O/f as a Q x Q lookup array indexed by the key (x conj f mod Q)."""
    Q = ch.Q
    T = np.zeros((Q, Q), dtype=complex)
    fc = conj(ch.f)
    for x in ch.reps:
        k = mul(x, fc)
        T[k[0] % Q, k[1] % Q] = ch.phi(x)
    return T


def keys(ch, A, B):
    c, d = conj(ch.f)
    k0 = (A * c - B * d) % ch.Q
    k1 = (A * d + B * c - B * d) % ch.Q
    return k0, k1


def lattice_np(nmax):
    r = int(math.isqrt(int(4 * nmax / 3) + 4)) + 2
    a = np.arange(-r, r + 1)
    A, B = np.meshgrid(a, a, indexing="ij")
    A = A.ravel(); B = B.ravel()
    n = A * A - A * B + B * B
    m = (n > 0) & (n <= nmax)
    A, B, n = A[m], B[m], n[m]
    o = np.argsort(n, kind="stable")
    return A[o], B[o], n[o].astype(float)


OM = complex(-0.5, SQ3 / 2)


class ThetaData:
    def __init__(self, ch):
        self.ch = ch
        Q = ch.Q
        self.Q = Q
        self.T = table(ch)
        self.vlo = 0.0864 / Q
        self.vhi = 20.0
        nmax = CUT / (c2 * self.vlo) + 10
        A, B, n = lattice_np(nmax)
        k0, k1 = keys(ch, A, B)
        ph = self.T[k0, k1]
        m = ph != 0
        self.A, self.B, self.n = A[m], B[m], n[m]
        self.z = self.A + self.B * OM
        self.coef = ph[m] * np.conj(self.z)          # phi(a) conj(a)
        self.gamma = ch.gamma()
        self.fc = cx(ch.f)
        self.W = ch.root_number()

    def theta(self, v):
        k = np.searchsorted(self.n, CUT / (c2 * v), side="right")
        return np.sum(self.coef[:k] * np.exp(-c2 * v * self.n[:k]))


# ---------------------------------------------------------------- T7-C9a / C9b
pois_err = 0.0; pois_f1 = 1e9; pois_f2 = 1e9
lines_err = 0.0; fal_b = {"no conj": [], "conj gamma": [], "+ sign": [], "conj f": []}
for ch in chars:
    td = ThetaData(ch)
    Q = ch.Q
    fc = conj(ch.f)
    for cv in (0.5, 1.0, 2.0):
        v = cv / math.sqrt(Q)
        # LHS by residue class: sum over a = t mod f of conj(a) e^{-c2 v N(a)}, all a (not only coprime ones)
        A, B, n = lattice_np(CUT / (c2 * v) + 10)
        k0, k1 = keys(ch, A, B)
        z = A + B * OM
        val = np.conj(z) * np.exp(-c2 * v * n)
        lhs = {}
        for a0, a1, w in zip(k0, k1, val):
            lhs[(a0, a1)] = lhs.get((a0, a1), 0) + w
        # RHS: (2/(sqrt3 Q)) sum_y e(t y/f) Ghat(y/f), Ghat(eta) = -(sqrt3/(2v^2)) eta exp(-2 pi |eta|^2/(sqrt3 v))
        Ay, By, ny = lattice_np(CUT * Q * v / c2 + 10)
        Ay = np.concatenate([[0], Ay]); By = np.concatenate([[0], By]); ny = np.concatenate([[0.0], ny])
        eta = (Ay + By * OM) / cx(ch.f)
        gh = -(SQ3 / (2 * v * v)) * eta * np.exp(-c2 * ny / (Q * v))
        e_t = f1_t = f2_t = 0.0
        for t in ch.reps:
            p, q = mul(t, fc)
            bb = (Ay * q + By * (p - q)) % Q                       # omega-coordinate of t y conj(f)
            ev = np.exp(2j * math.pi * bb / Q)
            rhs = 2 / (SQ3 * Q) * np.sum(ev * gh)
            rhs_m = 2 / (SQ3 * Q) * np.sum(np.conj(ev) * gh)        # falsifier e(-ty/f)
            l = lhs.get((p % Q, q % Q), 0)
            sc = max(abs(l), 1e-3)
            e_t = max(e_t, abs(l - rhs) / sc)
            f1_t = max(f1_t, abs(l - rhs * Q) / sc)                 # falsifier: Jacobian omitted
            f2_t = max(f2_t, abs(l - rhs_m) / sc)
        pois_err = max(pois_err, e_t)
        pois_f1 = min(pois_f1, f1_t)
        pois_f2 = min(pois_f2, f2_t)
    # lines of eq:d153 at 5 values of v
    for cv in (0.35, 0.6, 1.25, 1.7, 2.6):   # (v0 = Q^{-1/2} avoided: Theta(v0) = 0 when Theta is real and W = -1, chars F, G)
        v = cv / math.sqrt(Q)
        th = td.theta(v)
        Ay, By, ny = lattice_np(CUT * Q * v / c2 + 10)
        k0, k1 = keys(ch, Ay, By)
        ph = td.T[k0, k1]
        y = Ay + By * OM
        ex = np.exp(-c2 * ny / (Q * v))
        line1 = 2 / (SQ3 * Q) * np.sum(np.conj(ph) * math.sqrt(Q) * td.gamma * (-SQ3 / (2 * v * v)) * (y / td.fc) * ex)
        S = np.sum(np.conj(ph) * y * ex)
        line2 = -td.gamma / (td.fc * math.sqrt(Q)) * v**-2 * S
        th2 = td.theta(1 / (Q * v))
        rhs148 = -td.gamma / (td.fc * math.sqrt(Q)) * v**-2 * np.conj(th2)
        sc = abs(th)
        lines_err = max(lines_err, abs(line1 - th) / sc, abs(line2 - th) / sc, abs(S - np.conj(th2)) / abs(S), abs(rhs148 - th) / sc)
        fal_b["no conj"].append(abs(-td.gamma / (td.fc * math.sqrt(Q)) * v**-2 * th2 - th) / sc)
        fal_b["conj gamma"].append(abs(-np.conj(td.gamma) / (td.fc * math.sqrt(Q)) * v**-2 * np.conj(th2) - th) / sc)
        fal_b["+ sign"].append(abs(td.gamma / (td.fc * math.sqrt(Q)) * v**-2 * np.conj(th2) - th) / sc)
        fal_b["conj f"].append(abs(-td.gamma / (np.conj(td.fc) * math.sqrt(Q)) * v**-2 * np.conj(th2) - th) / sc)
report("T7-C9a", "eq:d153", "Poisson step: sum_z G(t+fz) = (2/(sqrt3 Q)) sum_y e(ty/f) Ghat(y/f), every t mod f",
       pois_err < 1e-11, f"8 characters x 3 v x all residues t: max rel. error {pois_err:.1e}; falsifiers (max over t, min over "
       f"characters and v): Jacobian omitted {pois_f1:.2f}, e(-ty/f) {pois_f2:.2f}", falsifier_failed=(pois_f1 > 1e-2 and pois_f2 > 1e-2))
fb = {k: max(v) for k, v in fal_b.items()}
report("T7-C9b", "eq:d153, eq:d148", "line 1 = line 2 = Theta(v) = -gamma/(f sqrt Q) v^-2 conj Theta(1/(Qv))",
       lines_err < 1e-11, f"8 characters x 5 v: max rel. error {lines_err:.1e}; falsifiers (max over characters, must be large "
       f"for some character): " + ", ".join(f"{k} {v:.2f}" for k, v in fb.items()),
       falsifier_failed=all(v > 1e-3 for v in fb.values()))


# ---------------------------------------------------------------- Mellin transforms
def trap_nodes(td, h):
    xs = np.arange(math.log(td.vlo), math.log(td.vhi) + h / 2, h)
    th = np.array([td.theta(math.exp(x)) for x in xs])
    return xs, th


def I_trap(xs, th, s, h, conjugate=False):
    f = (np.conj(th) if conjugate else th) * np.exp(xs * (s + 0.5))
    return h * np.sum(f)


def gl_nodes(td, lo, hi, m):
    g, w = np.polynomial.legendre.leggauss(m)
    xs = 0.5 * (hi - lo) * g + 0.5 * (hi + lo)
    return xs, 0.5 * (hi - lo) * w


L_dir_cache = {}


def L_ideal(td, s, X=200000):
    """L(s, psi) = sum over ideals (canonical generator x + y zeta, x > 0, y >= 0) of psi((a)) N(a)^{-s}."""
    key = (td.ch.name, X)
    if key not in L_dir_cache:
        r = int(math.isqrt(X)) + 2
        xx, yy = np.meshgrid(np.arange(1, r + 1), np.arange(0, r + 1), indexing="ij")
        xx = xx.ravel(); yy = yy.ravel()
        n = xx * xx + xx * yy + yy * yy
        m = n <= X
        xx, yy, n = xx[m], yy[m], n[m]
        A = xx + yy; B = yy                                  # x + y zeta = (x + y) + y omega
        k0, k1 = keys(td.ch, A, B)
        ph = td.T[k0, k1]
        z = A + B * OM
        psi = ph * np.conj(z) / np.abs(z)
        mm = ph != 0
        L_dir_cache[key] = (psi[mm], n[mm].astype(float))
    psi, n = L_dir_cache[key]
    return np.sum(psi * n ** (-s))


def elem_sum(td, s, X=200000):
    A, B, n = lattice_np(X)
    k0, k1 = keys(td.ch, A, B)
    ph = td.T[k0, k1]
    z = A + B * OM
    return np.sum(ph * np.conj(z) * n ** (-(s + 0.5)))


def gfac(s, kind="ok"):
    s = complex(s)
    if kind == "ok":
        return complex(mp.power(c2, -(s + 0.5)) * mp.gamma(s + 0.5))
    if kind == "gamma(s)":
        return complex(mp.power(c2, -(s + 0.5)) * mp.gamma(s))
    if kind == "2pi":
        return complex(mp.power(2 * mp.pi, -(s + 0.5)) * mp.gamma(s + 0.5))


h = 0.02
c5_err = 0.0; c5_quad = 0.0; c5_f = {"no 6": 1e9, "Gamma(s)": 1e9, "(2pi)^": 1e9}
c11_err = 0.0; c12_err = 0.0; c12_f = {"Q^(s-3/2)": 1e9, "B-(s)": 1e9}
c13_err = 0.0; c13_f = {"conj W": 0.0, "-W": 1e9}
c6_err = 0.0; c6_quad = 0.0; c6_f = {"-W": 1e9, "conj W": 0.0, "scale": 1e9}
lines = []
for ch in chars:
    td = ThetaData(ch)
    Q = ch.Q
    W = td.W
    xs, th = trap_nodes(td, h)
    xs2, th2 = xs[::2], th[::2]
    # --- C5 / C10
    for s in (2.5 + 3j, 3.0, 3.5 - 2j):
        for conjugate in (False, True):
            I = I_trap(xs, th, s, h, conjugate)
            I2 = I_trap(xs2, th2, s, 2 * h, conjugate)
            c5_quad = max(c5_quad, abs(I - I2) / abs(I))
            if not conjugate:
                Es = elem_sum(td, s)
                Ld = L_ideal(td, s)
            else:
                Es = np.conj(elem_sum(td, np.conj(s)))
                Ld = np.conj(L_ideal(td, np.conj(s)))
            r154 = gfac(s) * Es
            r149 = 6 * gfac(s) * Ld
            c5_err = max(c5_err, abs(I - r154) / abs(I), abs(I - r149) / abs(I))
            c5_f["no 6"] = min(c5_f["no 6"], abs(I - gfac(s) * Ld) / abs(I))
            c5_f["Gamma(s)"] = min(c5_f["Gamma(s)"], abs(I - 6 * gfac(s, "gamma(s)") * Ld) / abs(I))
            c5_f["(2pi)^"] = min(c5_f["(2pi)^"], abs(I - 6 * gfac(s, "2pi") * Ld) / abs(I))
    # --- C11, C12, C13
    v0 = Q ** -0.5
    lo_x, lo_w = gl_nodes(td, math.log(td.vlo), math.log(v0), 160)
    hi_x, hi_w = gl_nodes(td, math.log(v0), math.log(td.vhi), 160)
    th_lo = np.array([td.theta(math.exp(x)) for x in lo_x])
    th_hi = np.array([td.theta(math.exp(x)) for x in hi_x])
    th_inv = np.array([td.theta(1 / (Q * math.exp(x))) for x in lo_x])

    def Bq(w, conjugate=False):
        f = (np.conj(th_hi) if conjugate else th_hi) * np.exp(hi_x * (w + 0.5))
        return np.sum(hi_w * f)

    def Bser(w, conjugate=False):
        w = mp.mpc(w)
        tot = mp.mpc(0)
        k = np.searchsorted(td.n, CUT / (c2 * v0) + 5, side="right")
        for cf, n in zip(td.coef[:k], td.n[:k]):
            cf = complex(cf)
            cc = cf.conjugate() if conjugate else cf
            tot += cc * mp.power(c2 * n, -(w + 0.5)) * mp.gammainc(w + 0.5, c2 * n * v0)
        return complex(tot)

    for w in (1.5 + 2j, 2.5 - 1j, -0.5 - 2j, -1.5 + 1j):
        for cj in (False, True):
            a1, a2 = Bq(w, cj), Bser(w, cj)
            c11_err = max(c11_err, abs(a1 - a2) / abs(a2))
    for s in (1.5 + 2j, 2.5 - 1j, 0.3 + 5j, -0.4 - 2j):
        Iful = I_trap(xs, th, s, h)
        Ifulc = I_trap(xs, th, s, h, True)
        Bs, Bms, B1s, Bm1s = Bq(s), Bq(s, True), Bq(1 - s), Bq(1 - s, True)
        if s.real > 1:
            left = np.sum(lo_w * th_lo * np.exp(lo_x * (s + 0.5)))
            mid = -td.gamma / (td.fc * math.sqrt(Q)) * np.sum(lo_w * np.conj(th_inv) * np.exp(lo_x * (s - 1.5)))
            right = -td.gamma / (td.fc * math.sqrt(Q)) * Q ** (1.5 - s) * Bm1s
            c12_err = max(c12_err, abs(left - mid) / abs(left), abs(left - right) / abs(left))
            c12_f["Q^(s-3/2)"] = min(c12_f["Q^(s-3/2)"], abs(left + td.gamma / (td.fc * math.sqrt(Q)) * Q ** (s - 1.5) * Bm1s) / abs(left))
            c12_f["B-(s)"] = min(c12_f["B-(s)"], abs(left + td.gamma / (td.fc * math.sqrt(Q)) * Q ** (1.5 - s) * Bms) / abs(left))
        lhs = Q ** (s / 2) * Iful
        rhs = Q ** (s / 2) * Bs + W * Q ** ((1 - s) / 2) * Bm1s
        lhsc = Q ** (s / 2) * Ifulc
        rhsc = Q ** (s / 2) * Bms + np.conj(W) * Q ** ((1 - s) / 2) * B1s
        c13_err = max(c13_err, abs(lhs - rhs) / abs(lhs), abs(lhsc - rhsc) / abs(lhsc))
        if abs(W.imag) > 1e-6:
            c13_f["conj W"] = max(c13_f["conj W"], abs(lhs - (Q ** (s / 2) * Bs + np.conj(W) * Q ** ((1 - s) / 2) * Bm1s)) / abs(lhs))
        if s.real < 1:      # (for Re s > 1 the reflected term is small, so -W is tested where both terms matter)
            c13_f["-W"] = min(c13_f["-W"], abs(lhs - (Q ** (s / 2) * Bs - W * Q ** ((1 - s) / 2) * Bm1s)) / abs(lhs))
    # --- C6: FE from direct quadratures
    for s in (0.5 + 3j, 0.2 - 4j, 0.8 + 0.3j, -0.3 + 1.5j, 1.4 - 2j):
        lhs = Q ** (s / 2) * I_trap(xs, th, s, h)
        Ic = I_trap(xs, th, 1 - s, h, True)
        rhs = W * Q ** ((1 - s) / 2) * Ic
        lhs2 = Q ** (s / 2) * I_trap(xs2, th2, s, 2 * h)
        c6_quad = max(c6_quad, abs(lhs - lhs2) / abs(lhs))
        c6_err = max(c6_err, abs(lhs - rhs) / abs(lhs))
        c6_f["-W"] = min(c6_f["-W"], abs(lhs + rhs) / abs(lhs))
        if abs(W.imag) > 1e-6:
            c6_f["conj W"] = max(c6_f["conj W"], abs(lhs - np.conj(W) * Q ** ((1 - s) / 2) * Ic) / abs(lhs))
        c6_f["scale"] = min(c6_f["scale"], abs((3 * Q) ** (s / 2) * I_trap(xs, th, s, h) - rhs) / abs(lhs))
    lines.append(abs(abs(W) - 1))
    print(f"   {ch.name}: done ({time.time() - t0:.0f}s)", flush=True)

report("T7-C5", "eq:d149, eq:d154", "Mellin transform of Theta = (2pi/sqrt3)^{-s-1/2} Gamma(s+1/2) sum_a = 6 (...) L(s,psi); also conj",
       c5_err < 1e-7 and c5_quad < 1e-10,
       f"8 characters x 3 s x (psi, conj psi): max rel. error {c5_err:.1e} (Dirichlet tail < 1e-7); trapezoid h vs 2h: "
       f"{c5_quad:.1e}; falsifiers (min rel. error): " + ", ".join(f"{k} {v:.2f}" for k, v in c5_f.items()),
       falsifier_failed=all(v > 1e-2 for v in c5_f.values()))
report("T7-C11", "eq:d155", "B(w), B^-(w) by quadrature over (v0, oo) = incomplete-gamma series (as implemented in s2/t1)",
       c11_err < 1e-10, f"8 characters x 4 w x (B, B^-): max rel. error {c11_err:.1e}")
report("T7-C12", "eq:d156", "int_0^{v0} Theta v^{s+1/2} dv/v = -gamma/(f sqrt Q) int_0^{v0} conj Theta(1/(Qv)) v^{s-3/2} dv/v = -gamma/(f sqrt Q) Q^{3/2-s} B^-(1-s)",
       c12_err < 1e-10, f"8 characters x s in {{1.5+2i, 2.5-i}}: max rel. error {c12_err:.1e}; falsifiers (min rel. error): " +
       ", ".join(f"{k} {v:.2f}" for k, v in c12_f.items()), falsifier_failed=all(v > 1e-2 for v in c12_f.values()))
# symbolic parts of C13
s_ = sp.symbols("s")
Qs, gam, f_ = sp.symbols("Q gamma f", positive=True)
Wsym = -gam * sp.sqrt(Qs) / f_
sym_a = sp.simplify(-gam / (f_ * sp.sqrt(Qs)) - Wsym / Qs) == 0
sym_b = sp.simplify(sp.expand_power_base((2 * sp.pi / sp.sqrt(3)) ** (-s_ - sp.Rational(1, 2)), force=True)
                    - 3 ** sp.Rational(1, 4) * (2 * sp.pi) ** sp.Rational(-1, 2) * 3 ** (s_ / 2) * (2 * sp.pi) ** (-s_)) == 0
sym_c = sp.simplify(sp.powsimp(Qs ** (s_ / 2) * (Wsym / Qs) * Qs ** (sp.Rational(3, 2) - s_) - Wsym * Qs ** ((1 - s_) / 2), force=True)) == 0
report("T7-C13", "eq:s2-split", "Q^{s/2} I(s) = Q^{s/2} B(s) + W Q^{(1-s)/2} B^-(1-s) and its conjugate; kappa; -gamma/(f sqrt Q) = W/Q",
       c13_err < 1e-9 and sym_a and sym_b and sym_c,
       f"8 characters x 4 s (two with Re s < 1, I by direct quadrature): max rel. error {c13_err:.1e}; symbolic: "
       f"{sym_a}, {sym_b}, {sym_c}; falsifiers: conj W max rel. error {c13_f['conj W']:.2f}, -W min rel. error (Re s < 1) {c13_f['-W']:.2f}",
       falsifier_failed=(c13_f["conj W"] > 1e-2 and c13_f["-W"] > 1e-2))
report("T7-C6", "eq:d150", "Q^{s/2} I_psi(s) = W Q^{(1-s)/2} I_conjpsi(1-s) (Lambda FE) from direct theta quadratures, W = -gamma sqrt Q/f",
       c6_err < 1e-8 and c6_quad < 1e-8 and max(lines) < 1e-12,
       f"8 characters x 5 s (Re s from -0.3 to 1.4): max rel. error {c6_err:.1e}; trapezoid h vs 2h {c6_quad:.1e}; max ||W| - 1| = {max(lines):.1e}"
       f"; falsifiers: -W min rel. error {c6_f['-W']:.2f}, conj W max {c6_f['conj W']:.2f}, (3Q)^(s/2) min {c6_f['scale']:.2f}",
       falsifier_failed=(c6_f["-W"] > 1e-2 and c6_f["conj W"] > 1e-2 and c6_f["scale"] > 1e-2))
print(f"time {time.time() - t0:.1f}s")
finish()
