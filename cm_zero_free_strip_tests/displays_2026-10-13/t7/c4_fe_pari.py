"""c4: functional equation in the L-normalization and its consequences, with PARI/GP L-functions built from the Euler
(Dirichlet) coefficients: eq:d150 (C.6, explicit root number), eq:s2-fe-unc (C.14), eq:d157 (C.15), the size check for
eq:d158 (C.16), and for members of X_u: eq:d168 (D.10), eq:s3-reflected (D.8).

PARI convention: lfuncreate([a, 1, [1/2, 3/2], 1, N, eps]) is Lambda(s) = N^{s/2} Gamma_R(s+1/2) Gamma_R(s+3/2) L(s)
= 2 (2 pi)^{-1/2} (3Q)^{s/2} (2 pi)^{-s} Gamma(s+1/2) L(s) for N = 3Q, with Lambda(s) = eps conj-Lambda(1 - s); so eps
is the guide's W(psi).  lfuncheckfeq returns log2 of the FE defect (about -100 to -120 = holds to working precision).

Claims
  T7-C6p  eq:d150  with eps = W(psi) = -gamma(phi) sqrt(Q)/f computed in Python (NOT by lfunrootres), the FE holds for
          the eight characters of t7lib (lfuncheckfeq <= -100 bits at 60 digits; coefficients and W built in gp from exact data), and for conj psi with eps = conj W.
          Falsifiers: eps = -W, eps = conj W (non-real W), N = Q, Vga = [0, 1] (Gamma(s)): >= -15 bits.
          Normalization: PARI's L(2.5+it), t = 1, 40, -60, equals the Dirichlet series (Python, N <= 2e5).
  T7-C14  eq:s2-fe-unc  L(s,psi) = W (3Q/(2pi)^2)^{1/2-s} Gamma(3/2-s)/Gamma(1/2+s) L(1-s, conj psi) at 6 points s
          (Re s from -1.3 to 1.1, |Im s| <= 12), with both values from PARI; symbolic: eq:d150 with
          Lambda(s) = (3Q)^{s/2}(2pi)^{-s}Gamma(s+1/2)L(s) is equivalent to eq:s2-fe-unc (sympy).
          Falsifiers: finite-order quotient Gamma(1-s)/Gamma(s); Q/(2pi)^2 without the 3; conj W.
  T7-C15  eq:d157  |L(-1/10+it)| = (3Q)^{3/5}(2pi)^{-6/5} G(-1/10,t) |L(11/10-it, conj psi)|, t in {0, 1, -3, 10, -25, 60};
          and the size of |L(-1/10+it)|/(Q^{3/5}(3+|t|)^{6/5}) on that grid (bounded, as the << claims).
          Falsifier: (3Q)^{1/2} in place of (3Q)^{3/5}.
  T7-C16  eq:d158 (size check only; the lemma is analytic): a = 0.6, e = 1e-4, t in {0, 10, 40}: the winding number of
          L around the circle |s - (2+it)| = 2-a-2e is 0 (no zero in the disk: the hypothesis holds there), and the
          maximum of |L| + |1/L| on the circle of radius 2-a-6e (hence on the closed disk, by the maximum principle)
          is reported together with log(|L|+|1/L|)/log(2Q(3+|t|)^2).
  T7-D10  eq:d168  for CM members phi of X_u (sigma = -1 and +1; s3lib.Member, conductor by s3lib.conductor, psi* on
          primes as in tests_cm/phase2_2026-10-12/s3/t3_fe.py, coefficients built in gp at 60 digits): eps = lfunrootres has
          |eps| = 1, lfuncheckfeq <= -100,
          and L(s,phi*) = eps (3Q)^{1/2-s}(2pi)^{2s-1} Gamma(3/2-s)/Gamma(1/2+s) L(1-s, conj phi*) at 4 points.
  T7-D8   eq:s3-reflected  for the same members, 51/100 <= a <= 1, e in {1e-4, 9.9e-4}, t in {0, +-2, +-10, +-40}:
          R = |L(1-a-6e+it,phi*)| / ((3Q)^{a-1/2+6e}(3+|t|)^{2a-1+12e}|L(a+6e-it, conj phi*)|) <= 1 (so C_0 = 1 works on
          this grid), and R equals (2pi)^{1-2a-12e}|Gamma(x1-it)/Gamma(x2+it)|/(3+|t|)^{2a-1+12e} (the proof's identity).
          Falsifiers: t-exponent 2a-2+12e or Q-exponent a-1 makes max R > 1.
"""
import os, sys, math, cmath, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t7lib import *
import numpy as np
import mpmath as mp
import sympy as sp

t0 = time.time()
mp.mp.dps = 30
chars = test_chars()
B = 1500           # lfuncost for N <= 675 on |Re s - 1/2| <= 3, |Im s| <= 80 asks for <= 919 coefficients


def gpf(x):
    return float(x.replace(" ", ""))


def gq(s):
    s = mp.mpc(s)
    return complex(mp.gamma(1.5 - s) / mp.gamma(0.5 + s))


def zetaK(sig):
    sig = mp.mpf(sig)
    return mp.zeta(sig) * (mp.zeta(sig, mp.mpf(1) / 3) - mp.zeta(sig, mp.mpf(2) / 3)) / mp.power(3, sig)


def L_dir(ch, s, X=200000):
    a = ideal_coeffs(ch, X)
    n = np.arange(1, X + 1, dtype=float)
    return np.sum(a[1:] * n ** (-complex(s)))


S14 = [complex(-0.1, 3), complex(0.25, -7), complex(0.5, 0.5), complex(0.8, 12), complex(1.1, -2), complex(-1.3, 1)]
T15 = [0, 1, -3, 10, -25, 60]

c6_bits = []; c6_fal = []; c6n = 0.0; c6_wdiff = 0.0
c14_err = 0.0; c14_f = {"Gamma(1-s)/Gamma(s)": 1e9, "no 3": 1e9, "conj W": 0.0}
c15_err = 0.0; c15_f = 1e9; c15_ratio = 0.0
c16_lines = []; c16_ok = True
for ch in chars:
    Q = ch.Q
    W = ch.root_number()
    sc = [gp_char_setup(ch, B, "a"), "ac=conj(a); W=Wg;"]
    sc.append(f"L1=lfuncreate([a,1,[1/2,3/2],1,{3*Q},W]); L2=lfuncreate([ac,1,[1/2,3/2],1,{3*Q},conj(W)]);")
    sc.append("print(lfuncheckfeq(L1)); print(lfuncheckfeq(L2)); print(abs(W-(" + gp_c(W) + ")));")
    sc.append(f"print(lfuncheckfeq(lfuncreate([a,1,[1/2,3/2],1,{3*Q},-W])));")
    sc.append(f"print(lfuncheckfeq(lfuncreate([a,1,[1/2,3/2],1,{3*Q},conj(W)])));")
    sc.append(f"print(lfuncheckfeq(lfuncreate([a,1,[1/2,3/2],1,{Q},W])));")
    sc.append(f"print(lfuncheckfeq(lfuncreate([a,1,[0,1],1,{3*Q},W])));")
    pts = [complex(2.5, 1), complex(2.5, 40), complex(2.5, -60)] + S14 + [1 - z for z in S14]
    pts15 = [complex(-0.1, t) for t in T15] + [complex(1.1, -t) for t in T15]
    # disk data for C16
    a16, e16 = 0.6, 1e-4
    circ = []
    for t in (0, 10, 40):
        for R, m in ((2 - a16 - 2 * e16, 256), (2 - a16 - 6 * e16, 64)):
            circ += [complex(2, t) + R * cmath.exp(2j * math.pi * k / m) for k in range(m)]
    sc.append("L1i=lfuninit(L1,[3.0,80]); L2i=lfuninit(L2,[3.0,80]);")
    allpts = pts + pts15 + circ
    sc.append("P=[" + ",".join(gp_c(z) for z in allpts) + "];")
    sc.append("for(k=1,#P, z=lfun(L1i,P[k]); print(real(z)); print(imag(z)); z=lfun(L2i,P[k]); print(real(z)); print(imag(z)));")
    out = gp_run("\n".join(sc) + "\n", prec=60)
    bits = [int(out[0]), int(out[1])] + [int(x) for x in out[3:7]]
    wdiff = gpf(out[2])
    vals = [gpf(x) for x in out[7:]]
    V1 = {}; V2 = {}
    for k, z in enumerate(allpts):
        V1[z] = complex(vals[4 * k], vals[4 * k + 1])
        V2[z] = complex(vals[4 * k + 2], vals[4 * k + 3])
    c6_bits += bits[:2]
    c6_wdiff = max(c6_wdiff, wdiff)
    fal = [bits[2], bits[4], bits[5]] + ([bits[3]] if abs(W.imag) > 1e-6 else [])
    c6_fal += fal
    # normalization
    for z in (complex(2.5, 1), complex(2.5, 40), complex(2.5, -60)):
        ld = L_dir(ch, z)
        c6n = max(c6n, abs(V1[z] - ld) / abs(ld))
    # C14
    for s in S14:
        lhs = V1[s]
        rhs = W * complex(mp.power(3 * Q / (2 * mp.pi) ** 2, 0.5 - s)) * gq(s) * V2[1 - s]
        c14_err = max(c14_err, abs(lhs - rhs) / abs(lhs))
        f1 = W * complex(mp.power(3 * Q / (2 * mp.pi) ** 2, 0.5 - s)) * complex(mp.gamma(1 - mp.mpc(s)) / mp.gamma(mp.mpc(s))) * V2[1 - s]
        f2 = W * complex(mp.power(Q / (2 * mp.pi) ** 2, 0.5 - s)) * gq(s) * V2[1 - s]
        c14_f["Gamma(1-s)/Gamma(s)"] = min(c14_f["Gamma(1-s)/Gamma(s)"], abs(lhs - f1) / abs(lhs))
        c14_f["no 3"] = min(c14_f["no 3"], abs(lhs - f2) / abs(lhs))
        if abs(W.imag) > 1e-6:
            c14_f["conj W"] = max(c14_f["conj W"], abs(lhs - rhs * np.conj(W) / W) / abs(lhs))
    # C15
    for t in T15:
        s = complex(-0.1, t)
        lhs = abs(V1[s])
        G = abs(gq(s))
        rhs = (3 * Q) ** 0.6 * (2 * math.pi) ** -1.2 * G * abs(V2[complex(1.1, -t)])
        c15_err = max(c15_err, abs(lhs - rhs) / rhs)
        c15_f = min(c15_f, abs(lhs - (3 * Q) ** 0.5 * (2 * math.pi) ** -1.2 * G * abs(V2[complex(1.1, -t)])) / rhs)
        c15_ratio = max(c15_ratio, lhs / (Q ** 0.6 * (3 + abs(t)) ** 1.2))
    # C16
    k = len(pts) + len(pts15)
    for t in (0, 10, 40):
        big = [V1[z] for z in allpts[k:k + 256]]; k += 256
        small = [V1[z] for z in allpts[k:k + 64]]; k += 64
        wind = sum(cmath.phase(big[(j + 1) % 256] / big[j]) for j in range(256)) / (2 * math.pi)
        mx = max(abs(z) + 1 / abs(z) for z in small)
        Xh = 2 * Q * (3 + abs(t)) ** 2
        c16_ok &= abs(wind) < 0.01
        c16_lines.append(f"{ch.name.split(':')[0]} t={t}: winding {wind:+.3f}, max|L|+|1/L| = {mx:.3f}, "
                         f"log ratio {math.log(mx) / math.log(Xh):.3f}")
    print(f"   {ch.name}: feq bits {bits[:2]}, falsifiers {fal} ({time.time() - t0:.0f}s)", flush=True)

report("T7-C6p", "eq:d150", "FE with the explicit root number W = -gamma sqrt Q/f (PARI lfuncheckfeq, coefficients from psi)",
       max(c6_bits) <= -100 and c6n < 1e-8 and c6_wdiff < 1e-15,
       f"8 characters, psi and conj psi (coefficients and W built in gp at 60 digits): feq bits max {max(c6_bits)}; "
       f"|W(gp) - W(python)| <= {c6_wdiff:.1e}; PARI L(2.5+it), t = 1, 40, -60, vs Dirichlet series (N <= 2e5): max rel. error {c6n:.1e}; "
       f"falsifiers (-W, N=Q, Gamma(s), conj W): min bits {min(c6_fal)}", falsifier_failed=min(c6_fal) >= -15)
# symbolic equivalence of eq:d150 and eq:s2-fe-unc
s_ = sp.symbols("s")
Qs = sp.symbols("Q", positive=True)
Lam = lambda z: (3 * Qs) ** (z / 2) * (2 * sp.pi) ** (-z) * sp.gamma(z + sp.Rational(1, 2))
ratio = sp.simplify(sp.powsimp(sp.expand_power_base(Lam(1 - s_) / Lam(s_)
                    / ((3 * Qs / (2 * sp.pi) ** 2) ** (sp.Rational(1, 2) - s_) * sp.gamma(sp.Rational(3, 2) - s_) / sp.gamma(sp.Rational(1, 2) + s_)), force=True), force=True))
sym14 = (ratio == 1)
report("T7-C14", "eq:s2-fe-unc", "L(s,psi) = W (3Q/(2pi)^2)^{1/2-s} Gamma(3/2-s)/Gamma(1/2+s) L(1-s, conj psi)",
       c14_err < 1e-12 and sym14,
       f"8 characters x 6 s (PARI values): max rel. error {c14_err:.1e}; symbolic equivalence with eq:d150: {sym14}; "
       f"falsifiers: Gamma(1-s)/Gamma(s) min {c14_f['Gamma(1-s)/Gamma(s)']:.2f}, Q/(2pi)^2 min {c14_f['no 3']:.2f}, "
       f"conj W max {c14_f['conj W']:.2f}", falsifier_failed=all(v > 1e-2 for v in c14_f.values()))
Cbound = float(mp.power(3, 0.6) * mp.power(2 * mp.pi, -1.2) * zetaK(1.1))
report("T7-C15", "eq:d157", "|L(-1/10+it)| = (3Q)^{3/5}(2pi)^{-6/5} G(-1/10,t)|L(11/10-it, conj psi)|; size << Q^{3/5}(3+|t|)^{6/5}",
       c15_err < 1e-12 and c15_ratio <= Cbound,
       f"8 characters x 6 t: max rel. error {c15_err:.1e}; max |L(-1/10+it)|/(Q^(3/5)(3+|t|)^(6/5)) = {c15_ratio:.3f} "
       f"(<= 3^(3/5)(2pi)^(-6/5) zeta_K(11/10) = {Cbound:.3f}, times sup G/(3+|t|)^(6/5) <= 1); falsifier (3Q)^(1/2): min "
       f"rel. error {c15_f:.2f}", falsifier_failed=c15_f > 1e-2)
for l in c16_lines:
    print("   " + l)
report("T7-C16", "eq:d158", "size check: no zero in the disk of radius 2-a-2e (winding 0); max |L|+|1/L| on radius 2-a-6e",
       c16_ok, f"a = 0.6, e = 1e-4, 8 characters x t in {{0, 10, 40}}; see lines above (not a proof: the lemma is analytic)")

# ---------------------------------------------------------------- CM members of X_u (D.10, D.8)
P2 = os.path.join(TESTS, "phase2_2026-10-12", "s3")
sys.path.insert(0, P2)
import s3lib
from s3lib import Member, conductor as s3_conductor, split_prime, alpha as s3_alpha
from sympy import primerange


def prime_ideals(Bn):
    out = [(LAM, 3)]
    if 4 <= Bn:
        out.append(((-2, 0), 4))
    for p in primerange(5, Bn + 1):
        p = int(p)
        if p % 3 == 2:
            if p * p <= Bn:
                out.append(((-p, 0), p * p))
        else:
            pi = split_prime(p)
            out.append((pi, p)); out.append((prim_assoc(conj(pi)), p))
    return out


def s_prime_value(member, f, P, sigma, k):
    """psi*((P)) for P in {lambda, -2} not dividing f, via x = eps P y = 1 mod f (as in s3/t3_fe.py)."""
    for a_ in range(-60, 61):
        for b_ in range(-60, 61):
            y = (a_, b_)
            if not primary(y) or not s3lib.coprime_to(y, member.bad):
                continue
            for eps in UNITS:
                x = mul(mul(eps, P), y)
                if divides(f, sub(x, (1, 0))):
                    return s3_alpha(x) ** (k * sigma) / member(y)
    raise RuntimeError("no auxiliary y")


def member_prime_data(member, f, sigma, Bn):
    """[(q, j, A, B)] for the prime ideals of norm q <= Bn prime to f, with phi*(P) = zeta6^j alpha(A + B omega)^sigma."""
    out = []
    for P, q in prime_ideals(Bn):
        if divides(P, f):
            continue
        if P in (LAM, (-2, 0)):
            for a_ in range(-60, 61):
                done = False
                for b_ in range(-60, 61):
                    y = (a_, b_)
                    if not primary(y) or not s3lib.coprime_to(y, member.bad):
                        continue
                    for eps in UNITS:
                        x = mul(mul(eps, P), y)
                        if divides(f, sub(x, (1, 0))):
                            v = s3_alpha(x) ** sigma / member(y)          # = psi*((P)), as in s3/t3_fe.py
                            elt = mul(x, conj(y))
                            done = True
                            break
                    if done:
                        break
                if done:
                    break
            assert done
        else:
            v = member(P)
            elt = P
        r = v / s3_alpha(elt) ** sigma
        j = round(cmath.phase(r) / (math.pi / 3)) % 6
        assert abs(r - ZETA6[j]) < 1e-9, (P, r)
        out.append((q, j, elt[0], elt[1]))
    return out


def gp_member_setup(data, Bn, sigma):
    return "\n".join([GP_HEAD, "PD=[" + ";".join(f"{q},{j},{A},{b}" for q, j, A, b in data) + "];",
                      f"a=vector({Bn}); a[1]=1; for(k=1,#PD~, my(v=z6(PD[k,2])*al(PD[k,3]+PD[k,4]*w)^({sigma}), q=PD[k,1]); "
                      f"forstep(m=q,{Bn},q, a[m] += v*a[m\\q]));"])


rng = random.Random(7)
CASES = [("abar chibar(pi7)", [(split_prime(7), 1)], (1, 0), (), -1),
         ("alpha chi(2) chi(-5)", [((-5, 0), 1)], (1, 0), (((2, 0), 1),), 1),
         ("alpha chi(omega) chi(omega pi19)", [(split_prime(19), 1)], (0, 1), (((0, 1), 1),), 1),
         ("abar chi(lam)^3 chi(-1) chibar(pi7)", [(split_prime(7), 1)], (1, 0), ((LAM, 3), ((-1, 0), 1)), -1)]
d10_err = 0.0; d10_ok = True; d8_max = 0.0; d8_id = 0.0; d8_f1 = 0.0; d8_f2 = 0.0; mlines = []
AGRID = [0.51, 0.75, 1.0]; EGRID = [1e-4, 9.9e-4]; TGRID = [0, 2, -2, 10, -10, 40, -40]
S10 = [complex(0.3, 4), complex(-0.2, -9), complex(0.6, 0.2), complex(1.2, 15)]
for name, uf, unit, th, sigma in CASES:
    m = Member(uf, unit=unit, theta=th, sigma=sigma, k=1)
    ex, Q, f = s3_conductor(m, lambda x, s=sigma: s3_alpha(x) ** s, rng=rng)
    cost = gp_run(f"L=lfuncreate([x->1,1,[1/2,3/2],1,{3*Q},1]); print(lfuncost(L,[1.5,45])[1]);", prec=60)
    Bn = int(int(cost[0]) * 1.3) + 100
    data = member_prime_data(m, f, sigma, Bn)
    sc = [gp_member_setup(data, Bn, sigma), "ac=conj(a);",
          f"L0=lfuncreate([a,1,[1/2,3/2],1,{3*Q},0]); w=lfunrootres(L0)[3];",
          f"L1=lfuncreate([a,1,[1/2,3/2],1,{3*Q},w]); L2=lfuncreate([ac,1,[1/2,3/2],1,{3*Q},conj(w)]);",
          "print(lfuncheckfeq(L1)); print(lfuncheckfeq(L2)); print(real(w)); print(imag(w));",
          "L1i=lfuninit(L1,[1.5,45]); L2i=lfuninit(L2,[1.5,45]);"]
    pts = []
    for s in S10:
        pts += [(1, s), (2, 1 - s)]
    for aa in AGRID:
        for ee in EGRID:
            for t in TGRID:
                pts += [(1, complex(1 - aa - 6 * ee, t)), (2, complex(aa + 6 * ee, -t))]
    sc.append("P=[" + ",".join(gp_c(z) for _, z in pts) + "]; K=[" + ",".join(str(k) for k, _ in pts) + "];")
    sc.append("for(j=1,#P, z=if(K[j]==1, lfun(L1i,P[j]), lfun(L2i,P[j])); print(real(z)); print(imag(z)));")
    out = gp_run("\n".join(sc) + "\n", prec=60)
    b1, b2 = int(out[0]), int(out[1])
    w = complex(gpf(out[2]), gpf(out[3]))
    vals = [complex(gpf(out[4 + 2 * j]), gpf(out[5 + 2 * j])) for j in range(len(pts))]
    V = {(k, z): v for (k, z), v in zip(pts, vals)}
    d10_ok &= (b1 <= -100 and b2 <= -100 and abs(abs(w) - 1) < 1e-15)
    for s in S10:
        lhs = V[(1, s)]
        rhs = w * complex(mp.power(3 * Q, 0.5 - s) * mp.power(2 * mp.pi, 2 * s - 1)) * gq(s) * V[(2, 1 - s)]
        d10_err = max(d10_err, abs(lhs - rhs) / abs(lhs))
    for aa in AGRID:
        for ee in EGRID:
            for t in TGRID:
                num = abs(V[(1, complex(1 - aa - 6 * ee, t))])
                den0 = abs(V[(2, complex(aa + 6 * ee, -t))])
                Rr = num / ((3 * Q) ** (aa - 0.5 + 6 * ee) * (3 + abs(t)) ** (2 * aa - 1 + 12 * ee) * den0)
                x1, x2 = 0.5 + aa + 6 * ee, 1.5 - aa - 6 * ee
                pred = (2 * math.pi) ** (1 - 2 * aa - 12 * ee) * abs(complex(mp.gamma(mp.mpc(x1, -t)) / mp.gamma(mp.mpc(x2, t)))) \
                    / (3 + abs(t)) ** (2 * aa - 1 + 12 * ee)
                d8_max = max(d8_max, Rr)
                d8_id = max(d8_id, abs(Rr - pred) / pred)
                d8_f1 = max(d8_f1, Rr * (3 + abs(t)))
                d8_f2 = max(d8_f2, Rr * (3 * Q) ** (0.5 + 6 * ee))
    mlines.append(f"{name} (sigma={sigma:+d}): Q={Q}, coefficients {Bn}, feq bits {b1}/{b2}, |eps|-1 = {abs(w)-1:.1e}")
    print("   " + mlines[-1] + f" ({time.time() - t0:.0f}s)", flush=True)
report("T7-D10", "eq:d168", "L(s,phi*) = eps (3Q)^{1/2-s}(2pi)^{2s-1} Gamma(3/2-s)/Gamma(1/2+s) L(1-s, conj phi*), CM members of X_u",
       d10_ok and d10_err < 1e-10, f"{len(CASES)} CM members (sigma = -1, +1): FE defect and |eps| ok: {d10_ok}; identity at 4 s: "
       f"max rel. error {d10_err:.1e}")
report("T7-D8", "eq:s3-reflected", "|L(1-a-6e+it,phi*)| <= C0 (3Q)^{a-1/2+6e}(3+|t|)^{2a-1+12e}|L(a+6e-it, conj phi*)| with C0 = 1 on the grid",
       d8_max <= 1 and d8_id < 1e-9,
       f"{len(CASES)} CM members x a in {AGRID} x e in {EGRID} x t in {TGRID}: max ratio {d8_max:.3f}; ratio = (2pi)^(1-2a-12e) x "
       f"normalized gamma quotient: max rel. error {d8_id:.1e}; falsifiers: t-exponent 2a-2+12e max {d8_f1:.2f}, Q-exponent "
       f"a-1 max {d8_f2:.2f}", falsifier_failed=(d8_f1 > 1 and d8_f2 > 1))
print(f"time {time.time() - t0:.1f}s")
finish()
