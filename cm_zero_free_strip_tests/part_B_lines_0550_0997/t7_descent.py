"""Section 'Remove the cubes, then reduce the number of rows' (lines 851-994): claims B68-B70, B72, B74-B87.
(B71, eq:remove, is a cited estimate and is not tested; B88, the table, is assembled in run_all.py from the CM-vs-zeta
falsifiers of B17, B27, B46, B50, B53, B55, B68.)

B68  eq:inversion with the CM factors, X = 4000, three (k, f, xi); falsifiers: conj(alpha(h))^3 (zeta factor);
     N(h)^{1/2} in place of N(h).
B69  sum_{h | b'} mu_K(h) = 1_{b'=1} (primary divisors, N(b') <= 3000); alpha and Psi completely multiplicative.
B70  H <= Sigma, N(h)^3 <= H_c^3 = min(X, X^2/H^2): H + H^2 F N(h)^3/X <= 2 Sigma and X/N(h)^3 >= 1 (exact rationals).
     Falsifier: H_c^3 = min(X, X^2/H).
B72  H^2 <= X empties the supremum; a nonempty supremum has N(b)^3 > X^2/H^2 (exact rationals, integer norms).
B74  Y_{C,e0} -> Y_enlarged when N(e0)/N(C)^2 and c_* are suppressed; the natural y-range is <= Y_{C,e0} (sympy).
B75  N(h) ~ L^2/H and N(h f^2) ~ L^2 F^2/H (sympy).      B76 Y_enl/Y_nat = Sigma/(LF) >= 1.
B77  without enlargement H'_nat = H/(F^2 R^2) (= H at F = R = 1).      B78 local cancellation.
B79  retained-support chain: F' <= N(f') <= N(k') <= H L/(Sigma F N(r)^2) <= H', X'/N(j) >= Sigma F/H >= 1 (sampling).
B80  transfer ranges with H' = H L/(Sigma F R^2), Sigma' = L/R (sampling + sympy).
B81  eq:contract (sampling + sympy).  Falsifier: H D^{-3 kappa}.     B82 H'/Sigma' <= D^-kappa preserved.
B83  random chains reach a leaf within ceil(C0/(2 kappa)) steps.
B84  exactly six k with 0 < N(k) <= 1; base-case counting bound on computed CM averages.
B85  sanity: E(H, X, F)/(Sigma ||W||^2) with CM coefficients a = conj(alpha)^2 gamma_2 (bounded).
B86  (D^2/(H B^2))/(D/B) = D^-theta/B.      B87 kappa = theta/2 works for D >= (C/B)^{2/theta}.
Run time: about 1-2 minutes.
"""
import os, sys, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import *
import sympy as sp
from fractions import Fraction as Fr
random.seed(994)

def Wb(y): return math.exp(4 - 1/((y - 1)*(2 - y))) if 1 < y < 2 else 0.0      # ||W||_inf = 1
def xi_om(n): return complex(cubic3(n, (0, 1)))                                    # (omega/n)_3, a class function mod 9
XIS = {'trivial': lambda n: 1.0, '(omega/n)_3': xi_om}

# ---------------- B68: cube inversion
X = 4000.0
sqX = primaries(int(2*X), sqfree=True)
g2 = {n: gamma_j_fast(2, n) for n in sqX}
bs = primaries(int((2*X)**(1/3)) + 2)
def T(Y, k, f, xi, bexp=6):
    s = 0
    for b in bs:
        Nb3 = N(b)**3
        if Nb3 >= 2*Y: continue
        Pb = xi(b)*chi(b, k)*chi(b, f)**4
        if Pb == 0: continue
        wb = alpha_c(b).conjugate()**bexp*Pb**3*math.sqrt(N(b))
        for n in sqX:
            w = Wb(N(n)*Nb3/Y)
            if w: s += alpha_c(n).conjugate()**2*g2[n]*xi(n)*chi(n, k)*chi(n, f)**4*wb*w
    return s/math.sqrt(Y)
err = 0; fa3 = 0; fNh = 0; lhs_mag = []
P7 = [p for p in PR_ODD if N(p) == 7][0]; P13 = [p for p in PR_ODD if N(p) == 13][0]
for (k, f, xiname) in [((5, 3), (7, 9), '(omega/n)_3'), ((2, 0), (1, 0), 'trivial'), ((-1, 0), mul(P7, P13), '(omega/n)_3')]:
    xi = XIS[xiname]
    Psi = lambda n: xi(n)*chi(n, k)*chi(n, f)**4
    lhs = sum(alpha_c(n).conjugate()**2*g2[n]*Psi(n)*Wb(N(n)/X) for n in sqX)/math.sqrt(X)
    rhs = 0; r3 = 0; rN = 0
    for h in bs:
        if h != (1, 0) and mu(h) == 0: continue
        Nh3 = N(h)**3
        if Nh3 >= 2*X: continue
        m = mu(h) if h != (1, 0) else 1
        Th = T(X/Nh3, k, f, xi)
        rhs += m*alpha_c(h).conjugate()**6*Psi(h)**3/N(h)*Th
        r3 += m*alpha_c(h).conjugate()**3*Psi(h)**3/N(h)*Th
        rN += m*alpha_c(h).conjugate()**6*Psi(h)**3/math.sqrt(N(h))*Th
    err = max(err, abs(lhs - rhs)/abs(lhs)); fa3 = max(fa3, abs(lhs - r3)/abs(lhs)); fNh = max(fNh, abs(lhs - rN)/abs(lhs))
    lhs_mag.append(f"{abs(lhs):.3g}")
nh = sum(1 for h in bs if h != (1, 0) and N(h)**3 < 2*X and mu(h) != 0)
report("B68", err, 1e-12, f"X={X:g}, three (k, f, xi), {nh} nontrivial cube indices h; |LHS| = {', '.join(lhs_mag)}")
falsifier("B68", "conj(alpha(h))^3 (zeta factor) with the CM T", fa3, 1e-6)
falsifier("B68", "N(h)^{1/2} in place of N(h)", fNh, 1e-6)

# ---------------- B69
from collections import Counter
bad = 0; cnt = 0
for bp in primaries(3000):
    fac = Counter(factor(bp)); s = 0
    for exps in itertools.product(*[range(e + 1) for e in fac.values()]):
        h = (1, 0)
        for p, e in zip(fac.keys(), exps): h = mul(h, pw(p, e))
        s += (mu(h) if h != (1, 0) else 1)
    bad += s != (1 if bp == (1, 0) else 0); cnt += 1
ea = 0; ep = 0; sq = primaries(800)
for _ in range(400):
    m, n = random.choice(sq), random.choice(sq)
    ea = max(ea, abs(alpha_c(mul(m, n)) - alpha_c(m)*alpha_c(n)))
    for k, f in [((5, 3), (7, 9)), ((2, 1), P7)]:
        Psi = lambda x: xi_om(x)*chi(x, k)*chi(x, f)**4
        ep = max(ep, abs(Psi(mul(m, n)) - Psi(m)*Psi(n)))
report("B69", bad + ea + ep, 1e-12, f"{cnt} primary b' (N<=3000): sum over primary divisors of mu_K = 1_{{b'=1}}; alpha and Psi completely "
       f"multiplicative on 400 random pairs (not necessarily coprime): {ea:.1e}, {ep:.1e}")

# ---------------- B70, B72 (exact rationals)
def rnd_fr(lo, hi):            # random rational in [10^lo, 10^hi] (as Fraction)
    return Fr(int(10**random.uniform(lo, hi)*1000) + 1, 1000)
viol70 = 0; fal70 = 0; n70 = 0
for _ in range(20000):
    F = rnd_fr(0, 3); Xv = rnd_fr(0, 7); Sig = Xv*F; Hv = Sig*rnd_fr(-5, 0)
    if Hv > Sig or Hv < 1: continue
    Hc3 = min(Xv, Xv*Xv/(Hv*Hv))
    Nh = 1
    while Fr(Nh**3) <= Hc3:
        n70 += 1
        viol70 += not (Hv + Hv*Hv*F*Nh**3/Xv <= 2*Sig and Xv/Nh**3 >= 1)
        Nh += max(1, Nh//7)
    Hc3b = min(Xv, Xv*Xv/Hv); Nh = 1
    while Fr(Nh**3) <= Hc3b:
        fal70 += not (Hv + Hv*Hv*F*Nh**3/Xv <= 2*Sig)
        Nh += max(1, Nh//7)
report("B70", viol70, 0, f"{n70} exact rational cases (H <= Sigma, N(h)^3 <= H_c^3, integer N(h)): H + H^2 F N(h)^3/X <= 2 Sigma and X/N(h)^3 >= 1 "
       "(with N(f) < 2F the completed bound is <= 3 Sigma)")
falsifier("B70", "H_c^3 = min(X, X^2/H)", fal70, 0)
viol72 = 0; nonempty_bad = 0; fal72 = 0; n72 = 0
for _ in range(20000):
    Xv = rnd_fr(0, 6); Hv = rnd_fr(0, 4)
    Hc3 = min(Xv, Xv*Xv/(Hv*Hv))
    for Nb in range(1, int(Xv**Fr(1, 3)) + 3):
        Lb = Xv/Nb**3
        insup = (Fr(Nb**3) > Hc3) and Lb > 1
        n72 += 1
        if Hv*Hv <= Xv and insup: viol72 += 1
        if insup and not (Fr(Nb**3) > Xv*Xv/(Hv*Hv)): nonempty_bad += 1
        if Hv*Hv > Xv and insup: fal72 += 1
report("B72", viol72 + nonempty_bad, 0, f"{n72} exact (X, H, N(b)) cases: H^2 <= X never has N(b) > H_c with L_b > 1; every member of a nonempty supremum has N(b)^3 > X^2/H^2")
falsifier("B72", "the same emptiness for H^2 > X (supremum nonempty in some cases)", fal72, 0)

# ---------------- B74-B78, B80, B81, B86 (sympy)
Hs, Ls, Fs, Ss, Xs, Ds, kap, th, Bs, Rs, Nb3, Ne0, NC, cst, CPh, vv = sp.symbols('H L F Sigma X D kappa theta B R Nb3 Ne0 NC c_* C_Phi v', positive=True)
Y_CE = cst*Ss*Ls*Fs*Ne0/(Hs*NC**2)
Y_enl = Ss*Ls*Fs/Hs; Y_nat = Ls**2*Fs**2/Hs
ok74 = sp.simplify(Y_CE.subs({cst: 1, Ne0: 1, NC: 1}) - Y_enl) == 0
# natural y-range from OAI (source lines 2331-2334): N(y) <= 4 C_Phi v^2 L^2 F^2 N(e0)/(H N(C)^2); <= Y_{C,e0} iff c_* Sigma >= 4 C_Phi v^2 L F
natural = 4*CPh*vv**2*Ls**2*Fs**2*Ne0/(Hs*NC**2)
ratio74 = sp.simplify(Y_CE/natural)                         # = c_* Sigma/(4 C_Phi v^2 L F) >= 1 when c_* >= 4 C_Phi v^2 and Sigma >= LF
ok74b = sp.simplify(ratio74 - cst*Ss/(4*CPh*vv**2*Ls*Fs)) == 0
report("B74", 0 if ok74 and ok74b else 1, 0, "Y_{C,e0}|_{c_*=N(e0)=N(C)=1} = Sigma L F/H = Y_enlarged; Y_{C,e0}/(natural y-range) = c_* Sigma/(4 C_Phi v^2 L F) >= 1 "
       "when c_* >= 4 C_Phi v^2 and Sigma >= LF; weights H N(C)/(N(e0) L^2 F) and ell = L/(N(C)N(t)) match OAI (source lines 2236-2240)")
Nh_size = Ls**2/Hs                                           # N(h) <~ N(e0) N(z1 z2)/H with N(z_i) ~ L (suppressing divisors)
ok75 = sp.simplify(Nh_size*Fs**2 - Y_nat) == 0
report("B75", 0 if ok75 else 1, 0, "N(h) ~ L*L/H and N(y) = N(h)N(f)^2 ~ L^2 F^2/H = Y_natural (sympy)")
ok76 = sp.simplify(Y_enl/Y_nat - Ss/(Ls*Fs)) == 0
viol76 = 0
for _ in range(50000):
    Xv = 10**random.uniform(0, 8); Fv = 10**random.uniform(0, 3); Nb3v = 10**random.uniform(0, max(0.0, math.log10(Xv)))
    Lv = Xv/Nb3v; viol76 += (Xv*Fv)/(Lv*Fv) < 1 - 1e-12
report("B76", (0 if ok76 else 1) + viol76, 0, "Y_enl/Y_nat = Sigma/(LF) (sympy); >= 1 since L = X/N(b)^3 <= X (50000 samples)")
Hp = Hs*Ls/(Ss*Fs*Rs**2)
ok77 = sp.simplify(Hp*Ss/(Ls*Fs) - Hs/(Fs**2*Rs**2)) == 0 and sp.simplify((Hs/(Fs**2*Rs**2)).subs({Fs: 1, Rs: 1}) - Hs) == 0
report("B77", 0 if ok77 else 1, 0, "H'_nat = H' Sigma/(LF) = H/(F^2 R^2), equal to H at F = R = 1 (sympy)")
viol78 = sum(1 for a in (0, 1) if (1 + a) - a - (1 - a) != a)
report("B78", viol78, 0, "(1+1_{pi|k'}) - 1_{pi|k'} - 1_{pi nmid k'} = 1_{pi|k'} for both values")
falsifier("B78", "(1+1_{pi|k'}) - 1_{pi nmid k'} = 1_{pi|k'}", sum(1 for a in (0, 1) if (1 + a) - (1 - a) != a), 0)
Sp_ = Ls/Rs
ok80 = (sp.simplify(Hp/Sp_ - Hs/(Ss*Fs*Rs)) == 0)
viol80 = 0
for _ in range(100000):
    Fv = 10**random.uniform(0, 3); Xv = 10**random.uniform(0, 8); Sv = Xv*Fv; Hv = Sv*10**random.uniform(-6, 0)
    Lv = Xv/10**random.uniform(0, math.log10(Xv)) if Xv > 1 else 1.0; Rv = 10**random.uniform(0, 3)
    Hpv = Hv*Lv/(Sv*Fv*Rv*Rv); Spv = Lv/Rv
    viol80 += not (Hpv <= Hv*Lv/(Sv*Fv)*(1 + 1e-12) and Hpv/Spv <= Hv/Sv*(1 + 1e-12) and Spv <= Lv*(1 + 1e-12))
report("B80", (0 if ok80 else 1) + viol80, 0, "H'/Sigma' = H/(Sigma F R) (sympy); the three ranges hold on 100000 random sets with F, R >= 1",
       note="the transfer estimate itself (OAI Prop. 5.4) is cited and untested")
# B81, B82: contraction
viol81 = 0; viol82 = 0; fal81 = 0; n81 = 0
for _ in range(200000):
    Dv = 10**random.uniform(1, 30); kv = random.uniform(0.01, 0.5)
    Fv = 10**random.uniform(0, 3); Xv = 10**random.uniform(0, 10); Sv = Xv*Fv
    Hv = Sv*Dv**(-kv)*10**random.uniform(-3, 0)
    if Hv < 1: continue
    lo = max(0.0, 2*math.log10(Xv/Hv)); hi = math.log10(Xv)
    if lo >= hi: continue
    Nb3v = 10**random.uniform(lo, hi); Lv = Xv/Nb3v          # nonempty supremum: X^2/H^2 < N(b)^3 < X
    Hpv = Hv*Lv/(Sv*Fv); n81 += 1
    viol81 += not (abs(Hpv - Hv/(Nb3v*Fv*Fv)) <= 1e-9*Hpv and Hpv < Hv*(Hv/Sv)**2*(1 + 1e-12) and Hv*(Hv/Sv)**2 <= Hv*Dv**(-2*kv)*(1 + 1e-9))
    Rv = 10**random.uniform(0, 2); Spv = Lv/Rv; Hp2 = Hv*Lv/(Sv*Fv*Rv*Rv)
    viol82 += Hp2/Spv > Dv**(-kv)*(1 + 1e-9)
    fal81 += Hpv > Hv*Dv**(-3*kv)
ok81 = sp.simplify((Hs*Ls/(Ss*Fs)).subs({Ls: Xs/Nb3, Ss: Xs*Fs}) - Hs/(Nb3*Fs**2)) == 0
report("B81", viol81 + (0 if ok81 else 1), 0, f"{n81} random admissible sets with a nonempty supremum: H' = H/(N(b)^3F^2) < H(H/XF)^2 <= H D^-2kappa (and the first step symbolically)")
falsifier("B81", "H' <= H D^{-3 kappa}", fal81, 0)
report("B82", viol82, 0, f"H'/Sigma' <= H/Sigma <= D^-kappa on the same {n81} sets (R >= 1)")
ok86 = sp.simplify((Ds**2/(Ds**(1 + th)*Bs**2))/(Ds/Bs) - Ds**(-th)/Bs) == 0
report("B86", 0 if ok86 else 1, 0, "(D^2/(H B^2))/(D/B) = D^-theta/B with H = D^(1+theta) (sympy)")
falsifier("B86", "(D^2/(H B))/(D/B) = D^-theta/B", 0 if sp.simplify((Ds**2/(Ds**(1 + th)*Bs))/(Ds/Bs) - Ds**(-th)/Bs) == 0 else 1, 0)
# B87
viol87 = 0
for _ in range(20000):
    thv = random.uniform(0.001, 0.1); Cv = 10**random.uniform(0, 3); Bv = 10**random.uniform(0, 3)
    lD0 = (2/thv)*max(0.0, math.log(Cv/Bv))                      # log of the threshold D0 = (C/B)^(2/theta)
    for lD in [lD0 + 1e-7, lD0 + math.log(10), lD0 + math.log(1e6)]:  # compare logarithms (D is huge)
        viol87 += math.log(Cv) - thv*lD - math.log(Bv) > -thv/2*lD + 1e-9
report("B87", viol87, 0, "C D^-theta/B <= D^-theta/2 for D >= (C/B)^(2/theta) (20000 random theta in (0,0.1], C, B)")

# ---------------- B79: retained-support chain (zeta guide eq:retained, eq:childcolumn, which the CM text summarizes)
viol79 = 0; n79 = 0
while n79 < 100000:
    Fv = 10**random.uniform(0, 3); Xv = 10**random.uniform(0, 8); Sv = Xv*Fv; Hv = Sv*10**random.uniform(-6, 0)
    Lv = Xv/10**random.uniform(0, math.log10(Xv)) if Xv > 1 else 1.0
    top = Hv*Lv/(Sv*Fv)
    if top < 1: continue
    R = 10**random.uniform(0, 0.5*math.log10(top)); Nr = R*random.uniform(1, 2); cap = top/Nr**2
    if cap < 1: continue
    Hpv = top/R**2
    Nk = random.uniform(1, cap); Nf = random.uniform(1, Nk); Fp = Nf/random.uniform(1, 2)   # F' <= N(f') <= N(k') <= cap
    Fp = max(1.0, Fp)
    Xp = Lv/(R*Fp); Nj = random.uniform(1, Nr)
    viol79 += not (Hpv >= Fp*(1 - 1e-12) and Fp >= 1 and Xp/Nj >= Sv*Fv/Hv*(1 - 1e-9) and Sv*Fv/Hv >= 1 - 1e-12)
    n79 += 1
report("B79", viol79, 0, f"{n79} random sets with a retained r: H' >= F' >= 1 and child column scale X'/N(j) >= Sigma F/H >= 1",
       note="the inequality chain is the zeta guide's (eq:retained, eq:childcolumn); the CM text only summarizes it")

# ---------------- B83: random chains through the induction
random.seed(17); viol = 0; maxsteps = 0; nch = 0; leaves = Counter()
while nch < 50000:
    Dv = 10**random.uniform(3, 30); C0 = random.uniform(1, 4); kv = random.uniform(0.02, 0.5)
    Sv = Dv**(C0*random.random())
    if Sv < Dv**kv: continue
    Fv = Sv**random.random(); Xv = Sv/Fv
    Hv = math.exp(random.uniform(0, math.log(Sv) - kv*math.log(Dv)))
    nch += 1; steps = 0; bound = math.ceil(C0/(2*kv))
    while True:
        Sv = Xv*Fv
        if Hv <= 1: leaves['H<=1'] += 1; break
        if Hv*Hv <= Xv: leaves['empty cube supremum'] += 1; break
        extremal = random.random() < 0.5
        lo, hi = max(0.0, math.log(Xv*Xv/(Hv*Hv))), math.log(Xv)
        Nb3v = math.exp(lo + 1e-9) if extremal else math.exp(random.uniform(lo, hi)); Lv = Xv/Nb3v
        top = Hv*Lv/(Sv*Fv)
        if top < 1: leaves['no retained r'] += 1; break
        R = 1.0 if extremal else math.exp(random.uniform(0, 0.5*math.log(top)))
        Hp = top/R**2; Fp = 1.0 if extremal else math.exp(random.uniform(0, math.log(max(1.0, Hp))))
        Xp = Lv/(R*Fp)
        viol += Hp > Hv*Dv**(-2*kv)*(1 + 1e-9) or Hp/(Xp*Fp) > Dv**(-kv)*(1 + 1e-9)
        Hv, Xv, Fv = Hp, Xp, Fp; steps += 1
    maxsteps = max(maxsteps, steps); viol += steps > bound
report("B83", viol, 0, f"{nch} random chains: every step contracts by D^-2kappa and keeps H'/Sigma' <= D^-kappa; all end within ceil(C0/(2 kappa)) steps "
       f"(longest {maxsteps}); leaves {dict(leaves)}")

# ---------------- B84, B85: base case and sanity of E with CM coefficients
units = [k for k in lattice(1) if k != (0, 0)]
SQ2 = primaries(2600, sqfree=True)
def cols(Y): return [n for n in SQ2 if Y < N(n) < 2*Y]
def energy(H, Xv, F, coef):
    ns = cols(Xv); ks = [k for k in lattice(H) if k != (0, 0)]; fs = [f for f in SQ2 if F <= N(f) < 2*F]
    A_ = np.array([k[0] for k in ks]); B_ = np.array([k[1] for k in ks])
    chik = {n: chi_vec(n, A_, B_) for n in ns}
    tot = 0.0
    for f in fs:
        v = np.zeros(len(ks), complex)
        for n in ns:
            w = coef(n)*complex(chi(n, f))**4*Wb(N(n)/Xv)
            if w != 0: v += w*chik[n]
        tot += np.sum(np.abs(v)**2)
    return tot/(Xv*F), len(fs), len(ns)
g2b = {}
def a_cm(n):
    if n not in g2b: g2b[n] = alpha_c(n).conjugate()**2*gamma_j_fast(2, n)
    return g2b[n]
worst = 0; lines = []
for (Xv, F) in [(200.0, 1), (500.0, 1), (300.0, 7)]:
    E, nf, nn = energy(1, Xv, F, a_cm)
    s1 = sum(Wb(N(n)/Xv) for n in cols(Xv))
    tri = 6*nf*s1**2/(Xv*F)
    worst = max(worst, E/tri); lines.append(f"X={Xv:g},F={F}: E={E:.3g} <= 6(#f/F)(sum W)^2/X={tri:.3g}, E/X={E/Xv:.3g}")
report("B84", (len(units) != 6) + max(0.0, worst - 1), 1e-12, f"{len(units)} lattice points with 0<N(k)<=1 (the units); " + "; ".join(lines))
rowsE = []; worstE = 0
for (Xv, F) in [(300.0, 1), (600.0, 1), (1200.0, 1), (300.0, 7)]:
    rr = []
    for H in [10, 40, 160]:
        E = energy(H, Xv, F, a_cm)[0]; rr.append(f"H={H}: {E/(Xv*F):.3g}"); worstE = max(worstE, E/(Xv*F))
    rowsE.append(f"X={Xv:g},F={F}: " + ", ".join(rr))
Ec = energy(40, 600.0, 1, lambda n: 1.0)[0]/600.0
report("B85", worstE, 1.0, "E(H,X,F)/(Sigma ||W||^2) with a = conj(alpha)^2 gamma_2: " + "; ".join(rowsE) + f"; control a(n)=1 at X=600,H=40: {Ec:.3g}",
       note="analytic bound (D^eps, unspecified constant); sanity only")
