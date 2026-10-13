"""Section 10 of cm-eleven-twelfths-guide.tex, revised 2026-10-12: 'Remove the cubes, then reduce the number of rows'.

The revision ports the presentation of Section 8 of the zeta guide (precise transfer estimate, explicit induction,
retained support, exclusion removal) to the CM setting.  Every display (10.1)-(10.23) that is an identity or an explicit
inequality is tested here.  Test id S10.n refers to display (10.n).  The cited analytic estimates (10.4), (10.14), (10.17),
(10.21), (10.22) are not provable by computation; they get size checks only.

S10.1      (10.1)  cube inversion with mu_K(h) conj(alpha(h))^6 Psi(h)^3 / N(h): exact, X = 2000, two (k, f, xi).
                   Falsifier: the zeta factor conj(alpha(h))^3.
S10.2-3    (10.2), (10.3)  N(h)^3 <= H_c^3 = min(X, X^2/H^2) and H <= Sigma give H + H^2 F N(h)^3/X <= 2 Sigma, X/N(h)^3 >= 1
                   (exact rationals).  Falsifier: H_c^3 = min(X, X^2/H).
S10.4      text after (10.4): H^2 <= X empties the supremum; a nonempty supremum has N(b)^3 > X^2/H^2 (exact rationals).
S10.5-7    (10.5)-(10.7)  Y_{C,e0} with N(e0)/N(C)^2 and c_* suppressed is Y_enlarged; N(h) ~ L^2/H, y = h f^2 gives
                   Y_natural; Y_enlarged/Y_natural = Sigma/(LF) >= 1 when LF <= Sigma (sympy).
S10.8      (10.8)  H' = HL/(Sigma F R^2) equals H/(N(b)^3 F^2 R^2) at L = X/N(b)^3, Sigma = XF; without enlargement
                   (Sigma reset to LF) it is H/(F^2 R^2), equal to H at F = R = 1 (sympy).
S10.9      (10.9)  local cancellation, both cases.  Falsifier: the third term without its indicator.
S10.10     (10.10) retained support with actual Eisenstein elements f' | k' (k' = f' m): F' <= N(f') <= N(k') <= HL/(Sigma F N(r)^2)
                   <= H' whenever N(r) >= R.  Falsifier: N(r) < R.
S10.11a    twisted multiplicativity a(vn) = a(v) a(n) chi_n(v)^4 for coprime squarefree v, n, a = conj(alpha)^2 gamma_2 xi.
                   Falsifier: chi_n(v)^2.
S10.11     (10.11) the exclusion-removal identity, exact, with the CM coefficients, X' = 2000, two characters xi,
                   divisors r with one to three primes.  Falsifiers: C_1(X'/N(j); k', f') (j dropped from the twist);
                   no sign mu_K(j).
S10.11b    the consequence E_r(H, X, F) <= tau(r) sum_{j | r} E(H, X/N(j), F N(j)), both sides computed from (7.4).
S10.12-13  (10.12), (10.13)  the replacement preserves H', Sigma' = X'F' = L/R and H'/Sigma'; X'/N(j) >= Sigma F/H >= 1
                   (200 000 random retained sets).
S10.15     (10.15) the three scale inequalities, from H' = HL/(Sigma F R^2), Sigma' = L/R, R, F >= 1 (sympy + sampling).
S10.16-20  (10.16), (10.19), (10.20)  100 000 random chains through the induction (cube removal, transfer, exclusion removal):
                   standing hypotheses at every node, the transfer hypotheses max(H, LF) <= Sigma with Sigma = XF unchanged,
                   the contraction (10.19), H' <= D^{(j-1) kappa}, the surviving hypotheses (10.20), and chain length
                   <= ceil(C0/kappa).  Half the chains take the extremal (longest) child.
S10.18     (10.18) base case: exactly six k with 0 < N(k) <= 1; |a_xi(n)| <= 1; E <= 6 (#f/F)(sum|W|)^2/X <= C_I X ||W||^2
                   on computed CM averages.
S10.23     (10.23) (D^2/(H B^2))/(D/B) = D^{-theta}/B at H = D^{1+theta}; kappa = theta/2 for D >= (C/B)^{2/theta}.
S10.cite   the citations of OAI1112 in Section 10 (Lemma 4.4, Lemma 5.3, Proposition 5.4 with C^{4m+12} and [u/16, 4v],
                   Proposition 5.1, Section 5.5, Section 7, Lemma 7.3), recomputed from the LaTeX source of the manuscript.
S10.22     sanity (printed): E(H, X, F)/(Sigma ||W||_inf^2) with the CM coefficients stays bounded; control a(n) = 1.

Run time: a few minutes.  Run from this folder: python3 test_sec10_descent.py
"""
import os, sys, math, random, itertools, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "part_B_lines_0550_0997"))
from cmlib import *
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
random.seed(1012)
allok = True

def Wb(y): return math.exp(-1/((y - 1)*(2 - y))) if 1 < y < 2 else 0.0
WINF = math.exp(-4.0)                                   # max of Wb, at y = 3/2
def xi_om(n): return complex(cubic3(n, (0, 1)))         # (omega/n)_3, a class function mod 9
XIS = {'trivial': lambda n: 1.0, '(omega/n)_3': xi_om}
g2c = {}
def g2(n):
    if n not in g2c: g2c[n] = gamma_j_fast(2, n)
    return g2c[n]
def a(n, xi='trivial'):
    return alpha_c(n).conjugate()**2*g2(n)*XIS[xi](n)
P7, P13 = [p for p in PR_ODD if N(p) == 7][0], [p for p in PR_ODD if N(p) == 13][0]
P19, P31 = [p for p in PR_ODD if N(p) == 19][0], [p for p in PR_ODD if N(p) == 31][0]

# ---------------------------------------------------------------- S10.1 cube inversion
X = 2000.0
SQX = primaries(int(2*X), sqfree=True)
bs = primaries(int((2*X)**(1/3)) + 2)
def T(Y, k, f, xi, bexp=6):
    s = 0
    for b in bs:
        Nb3 = N(b)**3
        if Nb3 >= 2*Y: continue
        Pb = XIS[xi](b)*chi(b, k)*chi(b, f)**4
        if Pb == 0: continue
        wb = alpha_c(b).conjugate()**bexp*Pb**3*math.sqrt(N(b))
        for n in SQX:
            w = Wb(N(n)*Nb3/Y)
            if w: s += alpha_c(n).conjugate()**2*g2(n)*XIS[xi](n)*chi(n, k)*chi(n, f)**4*wb*w
    return s/math.sqrt(Y)
err = fa3 = 0
for (k, f, xi) in [((5, 3), (7, 9), '(omega/n)_3'), ((-1, 0), mul(P7, P13), 'trivial')]:
    Psi = lambda n: XIS[xi](n)*chi(n, k)*chi(n, f)**4
    lhs = sum(a(n, xi)*chi(n, k)*chi(n, f)**4*Wb(N(n)/X) for n in SQX)/math.sqrt(X)
    rhs = r3 = 0
    for h in bs:
        if h != (1, 0) and mu(h) == 0: continue
        Nh3 = N(h)**3
        if Nh3 >= 2*X: continue
        m = mu(h) if h != (1, 0) else 1
        Th = T(X/Nh3, k, f, xi)
        rhs += m*alpha_c(h).conjugate()**6*Psi(h)**3/N(h)*Th
        r3 += m*alpha_c(h).conjugate()**3*Psi(h)**3/N(h)*Th
    err = max(err, abs(lhs - rhs)/abs(lhs)); fa3 = max(fa3, abs(lhs - r3)/abs(lhs))
allok &= report("S10.1", err, 1e-12, f"(10.1) cube inversion, X={X:g}, two (k, f, xi)")
allok &= falsifier("S10.1", "conj(alpha(h))^3 (zeta factor)", fa3, 1e-6)

# ---------------------------------------------------------------- S10.2-3 and S10.4 (exact rationals)
def rnd_fr(lo, hi): return Fr(int(10**random.uniform(lo, hi)*1000) + 1, 1000)
viol = fal = cnt = 0
for _ in range(20000):
    F = rnd_fr(0, 3); Xv = rnd_fr(0, 7); Sig = Xv*F; Hv = Sig*rnd_fr(-5, 0)
    if Hv > Sig or Hv < 1: continue
    Hc3 = min(Xv, Xv*Xv/(Hv*Hv)); Nh = 1
    while Fr(Nh**3) <= Hc3:
        cnt += 1; viol += not (Hv + Hv*Hv*F*Nh**3/Xv <= 2*Sig and Xv/Nh**3 >= 1); Nh += max(1, Nh//7)
    Hc3b = min(Xv, Xv*Xv/Hv); Nh = 1
    while Fr(Nh**3) <= Hc3b:
        fal += not (Hv + Hv*Hv*F*Nh**3/Xv <= 2*Sig); Nh += max(1, Nh//7)
allok &= report("S10.2-3", viol, 0, f"(10.2), (10.3) {cnt} exact cases: H + H^2 F N(h)^3/X <= 2 Sigma and X/N(h)^3 >= 1")
allok &= falsifier("S10.2-3", "H_c^3 = min(X, X^2/H)", fal, 0)
viol = bad = cnt = 0
for _ in range(20000):
    Xv = rnd_fr(0, 6); Hv = rnd_fr(0, 4); Hc3 = min(Xv, Xv*Xv/(Hv*Hv))
    for Nb in range(1, int(Xv**Fr(1, 3)) + 3):
        insup = Fr(Nb**3) > Hc3 and Xv/Nb**3 > 1; cnt += 1
        viol += (Hv*Hv <= Xv and insup); bad += (insup and not Fr(Nb**3) > Xv*Xv/(Hv*Hv))
allok &= report("S10.4", viol + bad, 0, f"text after (10.4), {cnt} exact cases: H^2 <= X empties the supremum; "
                "a nonempty supremum has N(b)^3 > X^2/H^2")

# ---------------------------------------------------------------- S10.5-7, S10.8, S10.15, S10.23 (sympy)
Hs, Ls, Fs, Ss, Xs, Rs, Nb3, Ne0, NC, cst, Ds, th, Bs = sp.symbols('H L F Sigma X R Nb3 Ne0 NC c_* D theta B', positive=True)
Y_CE = cst*Ss*Ls*Fs*Ne0/(Hs*NC**2); Y_enl = Ss*Ls*Fs/Hs; Y_nat = Ls**2*Fs**2/Hs
ok = sp.simplify(Y_CE.subs({cst: 1, Ne0: 1, NC: 1}) - Y_enl) == 0
ok &= sp.simplify((Ls**2/Hs)*Fs**2 - Y_nat) == 0
ok &= sp.simplify(Y_enl/Y_nat - Ss/(Ls*Fs)) == 0
allok &= report("S10.5-7", 0 if ok else 1, 0, "(10.5)-(10.7): Y_{C,e0} -> Y_enlarged; N(h) F^2 ~ Y_natural; ratio Sigma/(LF) (>= 1 iff LF <= Sigma)")
Hp = Hs*Ls/(Ss*Fs*Rs**2)
ok = sp.simplify(Hp.subs({Ls: Xs/Nb3, Ss: Xs*Fs}) - Hs/(Nb3*Fs**2*Rs**2)) == 0
nat = sp.simplify(Hp.subs(Ss, Ls*Fs))
ok &= sp.simplify(nat - Hs/(Fs**2*Rs**2)) == 0 and sp.simplify(nat.subs({Fs: 1, Rs: 1}) - Hs) == 0
allok &= report("S10.8", 0 if ok else 1, 0, "(10.8): H' = H/(N(b)^3 F^2 R^2); without enlargement H/(F^2 R^2), = H at F = R = 1")
Sp_ = Ls/Rs
ok = sp.simplify(Hp/Sp_ - Hs/(Ss*Fs*Rs)) == 0                     # H'/Sigma' = H/(Sigma F R) <= H/Sigma
viol = 0
for _ in range(100000):
    H_ = 10**random.uniform(0, 8); F_ = 10**random.uniform(0, 3); S_ = max(H_, 1.0)*10**random.uniform(0, 4)
    L_ = S_/F_*10**random.uniform(-6, 0); R_ = 10**random.uniform(0, 3)
    h = H_*L_/(S_*F_*R_**2); s_ = L_/R_
    viol += not (h <= H_*L_/(S_*F_)*(1 + 1e-12) and h/s_ <= H_/S_*(1 + 1e-12) and s_ <= L_*(1 + 1e-12))
allok &= report("S10.15", viol + (0 if ok else 1), 0, "(10.15): H'/Sigma' = H/(Sigma F R); the three inequalities on 100000 samples with R, F >= 1")
Hinit = Ds**2/(Ds**(1 + th)*Bs**2); Sinit = Ds/Bs
ok = sp.simplify(Hinit/Sinit - Ds**(-th)/Bs) == 0
viol = 0
for _ in range(20000):
    t_ = random.uniform(0.001, 0.1); B_ = 10**random.uniform(0, 3); C_ = 10**random.uniform(-2, 4)
    lD_ = (2/t_)*math.log(C_/B_) + random.uniform(0, 7) if C_ > B_ else random.uniform(0, 70)   # log D, D >= (C/B)^{2/theta}
    lhs = math.log(C_) - t_*lD_ - math.log(B_); rhs = -t_/2*lD_
    viol += not (lhs <= rhs + 1e-9)
allok &= report("S10.23", viol + (0 if ok else 1), 0, "(10.23): H/Sigma = D^-theta/B at H = D^{1+theta}; "
                "C D^-theta/B <= D^{-theta/2} once D >= (C/B)^{2/theta} (20000 samples)")

# ---------------------------------------------------------------- S10.9 local cancellation
ok = all((1 + d) - d - (1 - d) == d for d in (0, 1))
fal = sum((1 + d) - d - 1 != d for d in (0, 1))
allok &= report("S10.9", 0 if ok else 1, 0, "(10.9) for pi | k' and pi does not divide k'")
allok &= falsifier("S10.9", "third term 1 instead of 1_{pi not| k'}", fal, 0)

# ---------------------------------------------------------------- S10.10 retained support with actual elements
SQf = primaries(3000, sqfree=True); LAT = [z for z in lattice(400) if z != (0, 0)]
viol = fal = cnt = 0
for _ in range(100000):
    fp = random.choice(SQf); m = random.choice(LAT); kp = mul(fp, m)       # f' | k' in Z[omega]
    Fp = 2**math.floor(math.log2(N(fp)))                                    # dyadic: F' <= N(f') < 2F'
    R_ = 2**random.randint(0, 6); Nr = R_*random.uniform(1, 2)              # R <= N(r) < 2R
    top = N(kp)*10**random.uniform(0, 3)                                    # H L/(Sigma F N(r)^2): r is retained
    F_ = 10**random.uniform(0, 2); L_ = top*Nr**2*F_*10**random.uniform(0, 2)
    S_ = L_*F_*10**random.uniform(0, 2); Hh = top*Nr**2*S_*F_/L_           # then H <= Sigma and LF <= Sigma
    cnt += 1; Hp_ = Hh*L_/(S_*F_*R_**2)
    viol += not (Fp <= N(fp) <= N(kp) <= Hh*L_/(S_*F_*Nr**2)*(1 + 1e-12) <= Hp_*(1 + 1e-9) and Hp_ >= Fp >= 1
                 and Hh <= S_*(1 + 1e-12) and L_*F_ <= S_*(1 + 1e-12))
    Nr_bad = R_*random.uniform(0.2, 0.9)
    fal += Hh*L_/(S_*F_*Nr_bad**2) > Hp_*(1 + 1e-12)                         # with N(r) < R the last step breaks
allok &= report("S10.10", viol, 0, f"(10.10) {cnt} retained (r, f', k') with f' | k' in Z[omega]: F' <= N(f') <= N(k') <= HL/(Sigma F N(r)^2) <= H', H' >= F' >= 1")
allok &= falsifier("S10.10", "N(r) < R (the last inequality fails)", fal, 0)

# ---------------------------------------------------------------- S10.11a twisted multiplicativity
SQm = primaries(1500, sqfree=True)
err = fal = 0; cnt = 0
while cnt < 300:
    v, n = random.choice(SQm), random.choice(SQm)
    if v == (1, 0) or n == (1, 0) or chi(n, v) == 0 or N(v)*N(n) > 60000: continue
    vn = mul(v, n)
    for xi in XIS:
        lhs = a(vn, xi)
        err = max(err, abs(lhs - a(v, xi)*a(n, xi)*chi(n, v)**4))
        fal = max(fal, abs(lhs - a(v, xi)*a(n, xi)*chi(n, v)**2))
    cnt += 1
allok &= report("S10.11a", err, 1e-12, "a(vn) = a(v) a(n) chi_n(v)^4 for 300 coprime squarefree pairs, two xi, a = conj(alpha)^2 gamma_2 xi")
allok &= falsifier("S10.11a", "chi_n(v)^2", fal, 1e-6)

# ---------------------------------------------------------------- S10.11 exclusion-removal identity
SQ = primaries(6000, sqfree=True)
def cols(Y): return [n for n in SQ if Y < N(n) < 2*Y]
def C1(Y, k, g, xi): return sum(a(n, xi)*chi(n, k)*chi(n, g)**4*Wb(N(n)/Y) for n in cols(Y))
def divs(ps):
    out = []
    for m_ in range(len(ps) + 1):
        for sub in itertools.combinations(ps, m_):
            d = (1, 0)
            for p in sub: d = mul(d, p)
            out.append((d, len(sub)))
    return out
X = 2000.0; err = f1 = f2 = 0; ncase = 0; mags = []
for xi in XIS:
    for (k, f, rps) in [((5, 3), (1, 0), [P7]), ((2, 0), (7, 9), [P13, P19]), ((-1, 1), mul(P7, (1, 3)), [P13, P19]),
                        ((4, -7), (1, 0), [P7, P13, P19])]:
        r = (1, 0)
        for p in rps: r = mul(r, p)
        lhs = sum(a(n, xi)*chi(n, k)*chi(n, f)**4*Wb(N(n)/X) for n in cols(X) if chi(n, r) != 0)
        rhs = b1 = b2 = 0
        for j, w in divs(rps):
            if j != (1, 0) and chi(j, f) == 0: continue
            ext = a(j, xi)*chi(j, k)*chi(j, f)**4 if j != (1, 0) else 1.0
            c = C1(X/N(j), k, mul(j, f), xi)
            rhs += (-1)**w*ext*c; b2 += ext*c
            b1 += (-1)**w*ext*C1(X/N(j), k, f, xi)
        err = max(err, abs(lhs - rhs)/abs(lhs)); f1 = max(f1, abs(lhs - b1)/abs(lhs)); f2 = max(f2, abs(lhs - b2)/abs(lhs))
        ncase += 1; mags.append(abs(lhs))
allok &= report("S10.11", err, 1e-12, f"(10.11) exact, {ncase} cases (X' = {X:g}, two xi, r with 1-3 primes); min |LHS| = {min(mags):.3g}")
allok &= falsifier("S10.11", "C_1(X'/N(j); k', f') without j in the twist", f1, 1e-6)
allok &= falsifier("S10.11", "no sign mu_K(j)", f2, 1e-6)

# ---------------------------------------------------------------- S10.11b and S10.18: computed averages (7.4)
SQe = primaries(2600, sqfree=True)
def energy(H, Xv, F, excl=None, coef=None):
    coef = coef or (lambda n: a(n))
    ns = [n for n in SQe if Xv < N(n) < 2*Xv and (excl is None or chi(n, excl) != 0)]
    ks = [k for k in lattice(H) if k != (0, 0)]; fs = [f for f in SQe if F <= N(f) < 2*F]
    A_ = np.array([k[0] for k in ks]); B_ = np.array([k[1] for k in ks])
    chik = {n: chi_vec(n, A_, B_) for n in ns}
    tot = 0.0
    for f in fs:
        v = np.zeros(len(ks), complex)
        for n in ns:
            w = coef(n)*complex(chi(n, f))**4*Wb(N(n)/Xv)
            if w != 0: v += w*chik[n]
        tot += float(np.sum(np.abs(v)**2))
    return tot/(Xv*F), len(fs), len(ns)
worst = 0; lines = []
for (H, Xv, F, rps) in [(12, 300.0, 1, [P7]), (25, 500.0, 1, [P7, P13]), (10, 400.0, 7, [P13])]:
    r = (1, 0)
    for p in rps: r = mul(r, p)
    Er = energy(H, Xv, F, excl=r)[0]
    rhs = 2**len(rps)*sum(energy(H, Xv/N(j), F*N(j))[0] for j, w in divs(rps))
    worst = max(worst, Er/rhs); lines.append(f"H={H},X={Xv:g},F={F},N(r)={N(r)}: {Er:.4g} <= {rhs:.4g}")
allok &= report("S10.11b", max(0.0, worst - 1), 0, f"E_r <= tau(r) sum_j E(H, X/N(j), F N(j)); largest ratio {worst:.3f}; " + "; ".join(lines))
units = [k for k in lattice(1) if k != (0, 0)]
amax = max(abs(a(n, xi)) for n in SQe for xi in XIS)
worst = 0; lines = []
for (Xv, F) in [(200.0, 1), (500.0, 1), (300.0, 7), (800.0, 13)]:
    E, nf, nn = energy(1, Xv, F)
    s1 = sum(Wb(N(n)/Xv) for n in SQe if Xv < N(n) < 2*Xv)
    tri = 6*nf*s1**2/(Xv*F); CI = 6*(nf/F)*(nn/Xv)**2
    worst = max(worst, E/tri, tri/(CI*Xv*WINF**2))
    lines.append(f"X={Xv:g},F={F}: E/(X||W||^2)={E/(Xv*WINF**2):.3g}, C_I={CI:.3g}")
allok &= report("S10.18", (len(units) != 6) + max(0.0, amax - 1 - 1e-12) + max(0.0, worst - 1), 1e-12,
                f"(10.18) {len(units)} units; max |a_xi(n)| = {amax:.15f} (N(n) <= 2600, two xi); " + "; ".join(lines))

# ---------------------------------------------------------------- S10.12-13 replacement scale algebra
bad = 0
for _ in range(200000):
    F = 10**random.uniform(0, 3); Xv = 10**random.uniform(0, 8); Sig = Xv*F; H = Sig*10**random.uniform(-6, 0)
    L = Xv/10**random.uniform(0, math.log10(Xv)) if Xv > 1 else 1.0
    top = H*L/(Sig*F)
    if top < 1: continue
    R = 10**random.uniform(0, 0.5*math.log10(top)); Nr = R*random.uniform(1, 2); cap = top/Nr**2
    if cap < 1: continue
    Hp_ = H*L/(Sig*F*R*R); Fp = 10**random.uniform(0, math.log10(cap)); Xp = L/(R*Fp); Nj = random.uniform(1, Nr)
    Xc, Fc = Xp/Nj, Fp*Nj
    ok = (abs(Xc*Fc - L/R) <= 1e-9*L/R and abs(Hp_/(Xc*Fc) - Hp_/(Xp*Fp)) <= 1e-9*Hp_/(Xp*Fp)
          and Xc >= Sig*F/H*(1 - 1e-9) and Sig*F/H >= 1 - 1e-12 and Fc >= 1)
    bad += not ok
allok &= report("S10.12-13", bad, 0, "(10.12), (10.13) on 200000 retained sets: H', Sigma' = L/R, H'/Sigma' preserved; X'/N(j) >= Sigma F/H >= 1")

# ---------------------------------------------------------------- S10.16-20 random induction chains
random.seed(17); viol = {}; nchains = 0; maxsteps = 0; leaves = {}
def chk(name, cond):
    if not cond: viol[name] = viol.get(name, 0) + 1
def level(H, D, kap): return 0 if H <= 1 else math.ceil(math.log(H)/(kap*math.log(D)) - 1e-12)
while nchains < 100000:
    D = 10**random.uniform(3, 30); C0 = random.uniform(1, 4); kap = random.uniform(0.02, 0.5); lD = math.log(D)
    Sig = D**(C0*random.random())
    if Sig < D**kap: continue
    F = Sig**random.random(); X = Sig/F
    H = math.exp(random.uniform(0, math.log(Sig) - kap*lD))
    nchains += 1; jmax = math.ceil(C0/kap); steps = 0
    chk("root level <= ceil(C0/kappa)", level(H, D, kap) <= jmax)
    while True:
        Sig = X*F
        chk("(10.16) H, X, F >= 1", H >= 1 - 1e-9 and X >= 1 - 1e-9 and F >= 1 - 1e-9)
        chk("(10.16) Sigma <= D^C0", Sig <= D**C0*(1 + 1e-9))
        chk("(10.16) H <= Sigma D^-kappa", H <= Sig*D**(-kap)*(1 + 1e-9))
        j = level(H, D, kap)
        if H <= 1: leaves['base case'] = leaves.get('base case', 0) + 1; break
        if H*H <= X: leaves['empty cube supremum'] = leaves.get('empty cube supremum', 0) + 1; break
        extremal = random.random() < 0.5
        lo, hi = max(0.0, math.log(X*X/(H*H))), math.log(X)
        Nb3 = math.exp(lo + 1e-9) if extremal else math.exp(random.uniform(lo, hi)); L = X/Nb3
        chk("transfer: max(H, LF) <= Sigma, Sigma = XF unchanged", max(H, L*F) <= Sig*(1 + 1e-9))
        chk("transfer: H, L, F, Sigma in [1, D^C0]", min(H, L, F, Sig) >= 1 - 1e-9 and max(H, L, F, Sig) <= D**C0*(1 + 1e-9))
        top = H*L/(Sig*F)
        if top < 1: leaves['no retained r'] = leaves.get('no retained r', 0) + 1; break
        R = 1.0 if extremal else math.exp(random.uniform(0, 0.5*math.log(top))); Nr = R if extremal else R*random.uniform(1, 2)
        if top/Nr**2 < 1: Nr = R
        cap = top/Nr**2; Hp_ = H*L/(Sig*F*R*R)
        Fp = 1.0 if extremal else math.exp(random.uniform(0, math.log(cap))); Xp = L/(R*Fp)
        Nj = 1.0 if extremal else random.uniform(1, Nr); Xc, Fc = Xp/Nj, Fp*Nj
        chk("(10.15) H' <= HL/(Sigma F)", Hp_ <= top*(1 + 1e-9))
        chk("(10.15)/(10.20) H'/Sigma' <= H/Sigma <= D^-kappa", Hp_/(Xc*Fc) <= H/Sig*(1 + 1e-9) and H/Sig <= D**(-kap)*(1 + 1e-9))
        chk("(10.20) Sigma' <= L <= Sigma <= D^C0", Xc*Fc <= L*(1 + 1e-9) and L <= Sig*(1 + 1e-9))
        chk("(10.20) H', X', F' >= 1", Hp_ >= 1 - 1e-9 and Xc >= 1 - 1e-9 and Fc >= 1 - 1e-9)
        chk("(10.19) H' <= H/(N(b)^3 F^2) < H (H/XF)^2 <= H D^-2kappa",
            Hp_ <= H/(Nb3*F*F)*(1 + 1e-9) and H/(Nb3*F*F) < H*(H/Sig)**2*(1 + 1e-9) and H*(H/Sig)**2 <= H*D**(-2*kap)*(1 + 1e-9))
        chk("H' <= D^{(j-1) kappa}", Hp_ <= D**((j - 1)*kap)*(1 + 1e-9))
        chk("(10.10) F' <= H'", Fp <= Hp_*(1 + 1e-9))
        H, X, F = Hp_, Xc, Fc; steps += 1
    maxsteps = max(maxsteps, steps)
    chk("chain length <= ceil(C0/kappa)", steps <= jmax)
allok &= report("S10.16-20", sum(viol.values()), 0, (str(viol) + " " if viol else "") +
                f"100000 chains; longest {maxsteps} steps; leaves {leaves}")

# ---------------------------------------------------------------- S10.cite: numbering of the cited results in OAI1112
_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
src = cmsources.OAI1112
nums = {}; sec = 0; sub = 0; thm = 0; appendix = False
envs = ("theorem", "proposition", "lemma", "corollary", "definition", "slogan", "remark")
cur = None
for line in open(src):
    if line.startswith("\\appendix"): appendix = True
    if re.match(r"\\section[\[{]", line): sec += 1; sub = 0; thm = 0; cur = ("sec", sec)
    elif re.match(r"\\subsection[\[{]", line): sub += 1; cur = ("sub", sec, sub)
    m = re.match(r"\\begin\{(" + "|".join(envs) + r")\}", line)
    if m: thm += 1; cur = ("thm", sec, thm)
    for lab in re.findall(r"\\label\{([^}]*)\}", line):
        if cur and not appendix: nums[lab] = ".".join(str(x) for x in cur[1:])
want = {"lem:remove-exclusions": "4.4", "lem:cube-reduction": "5.3", "prop:transfer": "5.4", "prop:canonical": "5.1",
        "sec:canonical-proof": "5.5", "sec:transfer-proof": "7", "lem:second-transfer": "7.3"}
text = open(src).read()
pt = text[text.index("\\begin{proposition}\\label{prop:transfer}"):]
pt = pt[:pt.index("\\end{proposition}")]
bad = [f"{k}: {nums.get(k)} != {v}" for k, v in want.items() if nums.get(k) != v]
bad += [] if "C^{4m+12}(I)" in pt else ["C^{4m+12} missing"]
bad += [] if "[u/16,4v]" in pt else ["[u/16,4v] missing"]
allok &= report("S10.cite", len(bad), 0, ("; ".join(bad) + " ") if bad else
                "Lemma 4.4, Lemma 5.3, Proposition 5.4 (C^{4m+12}, [u/16,4v]), Proposition 5.1, Section 5.5, Section 7, Lemma 7.3 all match the source")

# ---------------------------------------------------------------- S10.22 sanity
print("      S10.22 sanity [(10.22) is a bound]: E(H,X,F)/(Sigma ||W||_inf^2) with a = conj(alpha)^2 gamma_2; control a(n) = 1")
for (Xv, F) in [(300.0, 1), (600.0, 1), (1200.0, 1), (300.0, 7)]:
    row = []
    for H in [1, 10, 40, 160]:
        E = energy(H, Xv, F)[0]; Ec = energy(H, Xv, F, coef=lambda n: 1.0)[0]
        row.append(f"H={H}: {E/(Xv*F*WINF**2):.3g} (control {Ec/(Xv*F*WINF**2):.3g})")
    print(f"          X={Xv:g}, F={F}:  " + ";  ".join(row))
print("ALL PASS" if allok else "SOME FAILURES")
sys.exit(0 if allok else 1)
