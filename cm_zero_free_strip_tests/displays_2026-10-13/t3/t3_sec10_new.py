"""New checks for Section 10 of cm-eleven-twelfths-guide.tex (displays (10.1)-(10.41) as numbered on 2026-10-13).

These cover the displays of Section 10 that the existing suites (section10_descent, part_B t7_descent,
proofs_2026-10-12/b4) do not test directly.  Test ids are D10.<display number>.

D10.10   (10.10) eq:d059  r = tg/e, f' = Cew, k' = e0 e h2: for 3000 admissible (C, t, e0, g, e, w) built from actual
                  Eisenstein primes, r and f' are squarefree and coprime (so they are an exclusion divisor and an
                  auxiliary divisor); the definitions agree with OAI1112, proof of Lemma 7.3 (source text).
D10.30   (10.30) eq:d062 / (10.8) and the text after (10.33): with c_* = max(1, 4 C_Phi v^2) and LF <= Sigma, every
                  y = h f^2 with N(h) <= C_Phi v^2 L^2 N(e0)/(H N(C)^2), F <= N(f) < 2F has N(y) <= Y_{C,e0}
                  (exact rationals, actual elements).  Falsifier: c_* = max(1, C_Phi v^2).
D10.11   (10.11)/(10.37) kernel support of the second Poisson summation: N(h2) <= C_Phi b_*^2 ell^2 N(e)/(Y N(g)^2)
                  with b_* = 2v, Y = Y_{C,e0} gives N(k') <= HL/(Sigma F N(r)^2) (exact rationals, actual elements).
                  Falsifier: b_* = 4v.
D10.16   (10.16) eq:childcolumn: X'/N(j) >= Sigma F N(r)^2/(H R N(j)) >= Sigma F/H >= 1 for retained (r, f', k'),
                  j | r, all three inequalities, exact rationals with actual elements.  Falsifier: N(r) in place of
                  N(r)^2 in the middle term.
D10.21   (10.21) eq:basecase: (F/(XF))(X ||W||)^2 = X ||W||^2 <= Sigma ||W||^2 (sympy, F >= 1).
D10.35   (10.35) eq:d065: sum_{f,h} |P_{C,e0,t}(hf^2;U)|^2 <= m_max sum_{N(y) <= Y_{C,e0}} |P_{C,e0,t}(y;U)|^2 <= m_max Q
                  (Phi >= 1 on [0,1]), m_max = max multiplicity of (f, h) -> hf^2 <= tau(y), computed block by block
                  with the CM coefficients mu conj(alpha) xi_1.  Falsifier: containment with c_* = max(1, C_Phi v^2).
D10.38   (10.38) eq:d067: the character combination, exact on 3000 random tuples; conj(chi_z(e)) = chi_z(e)^5.
                  Falsifiers: psi(e) read as chi_{z1}(e) conj(chi_{z2}(e)); chi_z(Ce)^2.
D10.39   (10.39) eq:d068: 1_{(n,gCtw)=1} chi_n(e0 e h2) chi_n(Cew)^4 = 1_{(n,tg/e)=1} (same), exact on 4000 random
                  tuples.  Falsifier: 1_{(n,t)=1} (g/e dropped).
D10.41   (10.41) eq:d070 and (10.17): J_1 = 2m+4 gives 2J_1 + 4 = 4m + 12 = J_t(m) (sympy); E(H,L,F) <= A(W) on
                  computed CM averages with a Phi >= 1 on [0,1]; derivative orders C^{2j+4} (Lemma 7.1),
                  C^{2m+4} (Lemma 7.3), C^{4m+12} and [u/16,4v] (Proposition 5.4) in the OAI1112 source.
D10.27   (10.27) eq:d060, zeta column: mu_K(n), conj(alpha) gamma_2, conj(alpha(b))^3, one zbar-derivative in the zeta
                  guide source (sources/zeta_11_12_guide.tex); alpha(ell) and d_zbar in OAI1112.  (The CM column is B88.)
S10.32   (10.32) size check: Z/(H Phihat(0) ||W||^2) = (1/(LF)) sum_f sum_C |W(N(C)/L)|^2 prod(1 - 1/N(pi)) stays
                  bounded (the O_I(1) constant), L up to 3000, F up to 100.
S10.34   (10.34) size check: sup over the Phihat argument c of ||K_c||_{C^q(I^2)}/||W||_{C^q(I)}^2, q <= 2.

Run: python3 t3_sec10_new.py   (about 1 minute).  Exit status 1 if any check fails or a falsifier does not fail.
"""
import os, sys, math, random, itertools, re
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(TESTS, "part_B_lines_0550_0997"))
from cmlib import *                 # report, falsifier, chi, mu, factor, primaries, mul, div, divides, N, ...
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
random.seed(20261013)
_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
OAI = cmsources.OAI1112
ZG = cmsources.ZG                                     # source of the zeta guide [ZG]

ONE = (1, 0)
SMALLP = [p for p in PR_ODD if N(p) <= 200]           # primary primes coprime to 6, norm <= 200
def prod(ps):
    r = ONE
    for p in ps: r = mul(r, p)
    return r
def coprime(a, b): return all(not divides(p, b) for p in set(factor(a))) if a != ONE else True
def sqfree(a): return a == ONE or mu(a) != 0
def muu(n): return 1 if n == ONE else mu(n)
def al(n): return 1.0 + 0j if n == ONE else alpha_c(n)
def chi1(n, x): return (1.0 + 0j) if n == ONE else complex(chi(n, x))
def tau(n): return 1 if n == ONE else math.prod(e + 1 for e in __import__("collections").Counter(factor(n)).values())
def divs_sqfree(n): return [ONE] if n == ONE else divisors_sqfree(n)
LOESCH = sorted({a*a - a*b + b*b for a in range(-200, 201) for b in range(-200, 201)})   # norms of Z[omega]; complete up to 30000
def max_norm_le(x):
    """Largest norm of a nonzero element of Z[omega] that is <= x (0 if none)."""
    best = 0
    for m in LOESCH:
        if m == 0: continue
        if m > x: break
        best = m
    return best

# ------------------------------------------------------------------ D10.10: r = tg/e, f' = Cew squarefree and coprime
def rand_sqfree(pool, kmax):
    return prod(random.sample(pool, random.randint(0, kmax)))
bad = 0; cnt = 0
while cnt < 3000:
    ps = random.sample(SMALLP, 9)
    C = prod(ps[0:random.randint(0, 2)]); t = prod(ps[2:2 + random.randint(0, 2)])
    gps = ps[4:4 + random.randint(0, 3)]; g = prod(gps)
    e = prod([p for p in gps if random.random() < 0.5]); w = prod(ps[7:7 + random.randint(0, 2)])
    e0 = prod([p for p in (factor(C) if C != ONE else []) if random.random() < 0.5])
    # admissibility (10.36): (g, Ct) = (w, gCt) = 1 (built in), e | g, e0 | C
    assert coprime(g, mul(C, t)) and coprime(w, mul(mul(g, C), t)) and divides(e, g) and divides(e0, C)
    r = mul(t, div(g, e)); fp = mul(mul(C, e), w)
    bad += not (sqfree(r) and sqfree(fp) and coprime(r, fp)); cnt += 1
src = open(OAI).read()
i0 = src.index("\\begin{lemma}\\label{lem:second-transfer}")
ok_src = "r=tg/e,\\quad f'=Cew,\\quad k'=deh." in src[i0:i0 + 3000]
report("D10.10", bad + (0 if ok_src else 1), 0, f"(10.10) {cnt} admissible (C,t,e0,g,e,w) from actual primes: r = tg/e and "
       f"f' = Cew squarefree and coprime; OAI1112 proof of Lemma 7.3 has r=tg/e, f'=Cew, k'=deh (d = e0): {ok_src}")

# ------------------------------------------------------------------ D10.30 and D10.11: support bookkeeping
def rnd(lo, hi): return Fr(int(10**random.uniform(lo, hi)*64) + 1, 64)
LAT_SMALL = [z for z in lattice(3000) if z != (0, 0)]
viol30 = fal30 = viol11 = fal11 = n30 = n11 = 0
for it in range(4000):
    v = random.choice([Fr(2), Fr(3, 2), Fr(4), Fr(5, 2)]); CPhi = rnd(-0.3, 0.7)
    cst = max(Fr(1), 4*CPhi*v*v); cbad = max(Fr(1), CPhi*v*v)
    ps = random.sample(SMALLP[:20], 8)
    C = prod(ps[0:random.randint(0, 2)]); t = prod(ps[2:2 + random.randint(0, 1)])
    e0 = prod([p for p in (factor(C) if C != ONE else []) if random.random() < 0.5])
    F = rnd(0, 1.2); L = rnd(0.5, 2.5) * N(C) * N(t)
    Sig = L*F*rnd(0, 1) if random.random() < 0.7 else L*F            # LF <= Sigma, often extremal
    H = Sig/rnd(0, 1.5)
    if H < 1: H = Fr(1)
    if H > Sig: continue
    Y = cst*Sig*L*F*N(e0)/(H*N(C)**2); Ybad = cbad*Sig*L*F*N(e0)/(H*N(C)**2)
    # (10.30): first Poisson frequencies and y = h f^2
    hmax = CPhi*v*v*L*L*N(e0)/(H*N(C)**2)
    nh = max_norm_le(hmax)
    if nh > 0:
        fs = [f for f in primaries(int(2*F) + 1, sqfree=True) if F <= N(f) < 2*F and coprime(f, mul(C, t))]
        if fs:
            f = max(fs, key=N); n30 += 1
            Ny = nh*N(f)**2
            viol30 += not (Ny <= 4*CPhi*v*v*L*L*F*F*N(e0)/(H*N(C)**2) <= Y)
            fal30 += Ny > Ybad
    # (10.11): second Poisson frequency h2 and k' = e0 e h2, r = tg/e
    ell = L/(N(C)*N(t)); bst = 2*v
    gps = ps[4:4 + random.randint(0, 2)]; g = prod(gps); e = prod([p for p in gps if random.random() < 0.5])
    if N(g) > bst*ell: continue
    h2max = CPhi*bst*bst*ell*ell*N(e)/(Y*N(g)**2)
    nh2 = max_norm_le(h2max)
    if nh2 == 0: continue
    n11 += 1
    Nk = N(e0)*N(e)*nh2; r = mul(t, div(g, e))
    viol11 += not (Nk <= H*L/(Sig*F*N(r)**2))
    h2bad = CPhi*(2*bst)**2*ell*ell*N(e)/(Y*N(g)**2)
    fal11 += N(e0)*N(e)*max_norm_le(h2bad) > H*L/(Sig*F*N(r)**2)
report("D10.30", viol30, 0, f"(10.30)/(10.8) {n30} cases (exact rationals, the largest actual N(h) and N(f)): "
       "N(hf^2) <= 4 C_Phi v^2 L^2F^2 N(e0)/(H N(C)^2) <= Y_{C,e0} with c_* = max(1, 4C_Phi v^2), LF <= Sigma")
falsifier("D10.30", "c_* = max(1, C_Phi v^2)", fal30, 0)
report("D10.11", viol11, 0, f"(10.11)/(10.37) {n11} cases: N(h2) <= C_Phi (2v)^2 ell^2 N(e)/(Y_{{C,e0}} N(g)^2) gives "
       "N(k') = N(e0 e h2) <= HL/(Sigma F N(r)^2), r = tg/e (kernel support of the second Poisson summation)")
falsifier("D10.11", "b_* = 4v (wrong support of U)", fal11, 0)

# ------------------------------------------------------------------ D10.16: child column chain
SQf = primaries(2000, sqfree=True)
viol = fal = cnt = 0
while cnt < 50000:
    ps = random.sample(SMALLP[:30], 6)
    r = prod(ps[:random.randint(0, 3)]); R = Fr(2)**int(math.floor(math.log2(N(r))))      # R <= N(r) < 2R
    fp = prod(ps[3:3 + random.randint(0, 2)]); Fp = Fr(2)**int(math.floor(math.log2(N(fp))))
    m = random.choice(LAT_SMALL[:400]); kp = mul(fp, m)                                   # f' | k'
    F = rnd(0, 2); top = Fr(N(kp))*rnd(0, 1)                                              # HL/(Sigma F N(r)^2) >= N(k')
    L = top*N(r)**2*F*rnd(0, 2); Sig = L*F*rnd(0, 1.5); H = top*N(r)**2*Sig*F/L
    if H > Sig or H < 1: continue
    cnt += 1
    Xp = L/(R*Fp)
    for j in divs_sqfree(r):
        a1 = Xp/N(j); a2 = Sig*F*N(r)**2/(H*R*N(j)); a3 = Sig*F/H
        viol += not (a1 >= a2 >= a3 >= 1)
        fal += not (Sig*F*N(r)/(H*R*N(j)) >= a3)
report("D10.16", viol, 0, f"(10.16) {cnt} retained (r, f', k') with f' | k' in Z[omega] and all j | r (exact rationals): "
       "X'/N(j) >= Sigma F N(r)^2/(H R N(j)) >= Sigma F/H >= 1")
falsifier("D10.16", "N(r) in place of N(r)^2 in the middle term", fal, 0)

# ------------------------------------------------------------------ D10.21: base case algebra
Xs, Fs, Ws, Ss = sp.symbols("X F W Sigma", positive=True)
ok = sp.simplify(Fs/(Xs*Fs)*(Xs*Ws)**2 - Xs*Ws**2) == 0
ok &= sp.simplify((Xs*Fs*Ws**2 - Xs*Ws**2) - Xs*Ws**2*(Fs - 1)) == 0      # Sigma W^2 - X W^2 = X W^2 (F - 1) >= 0
report("D10.21", 0 if ok else 1, 0, "(10.21) F/(XF) (X||W||)^2 = X||W||^2 and Sigma||W||^2 - X||W||^2 = X||W||^2 (F-1) >= 0 (sympy)")

# ------------------------------------------------------------------ D10.38: character combination
def chpow(z, x, k): return chi1(z, x)**k
viol = f1 = f2 = 0; cnt = 0; c5 = 0
while cnt < 3000:
    ps = random.sample(SMALLP[:25], 9)
    z1 = prod(ps[0:random.randint(0, 2)]); z2 = prod(ps[2:2 + random.randint(0, 2)])
    pool = ps[4:] + ps[:4]                 # C, e, e0 may share primes with z1, z2 (then both sides vanish)
    C = prod(random.sample(pool, random.randint(0, 2)))
    e = prod(random.sample(pool, random.randint(0, 2)))
    if not sqfree(C) or not sqfree(e): continue
    e0 = prod([p for p in (factor(C) if C != ONE else []) if random.random() < 0.5])
    h2 = random.choice(LAT_SMALL)
    psi = lambda x: chi1(z1, x).conjugate()*chi1(z2, x)
    C4e0 = mul(pw(C, 4), e0)
    lhs = psi(e)*psi(h2).conjugate()*chi1(z1, C4e0)*chi1(z2, C4e0).conjugate()
    d = mul(mul(e0, e), h2); Ce = mul(C, e)
    rhs = chi1(z1, d)*chpow(z1, Ce, 4)*(chi1(z2, d)*chpow(z2, Ce, 4)).conjugate()
    viol = max(viol, abs(lhs - rhs))
    psib = lambda x: chi1(z1, x)*chi1(z2, x).conjugate()
    lb = psib(e)*psi(h2).conjugate()*chi1(z1, C4e0)*chi1(z2, C4e0).conjugate()
    f1 = max(f1, abs(lb - rhs))
    rb = chi1(z1, d)*chpow(z1, Ce, 2)*(chi1(z2, d)*chpow(z2, Ce, 2)).conjugate()
    f2 = max(f2, abs(lhs - rb))
    for z in (z1, z2):
        c5 = max(c5, abs(chi1(z, e).conjugate() - chi1(z, e)**5)) if z != ONE else c5
    cnt += 1
report("D10.38", max(viol, c5), 1e-12, f"(10.38) {cnt} random (z1, z2, C, e0 | C, e, h2), shared primes allowed: "
       "psi(e)conj(psi(h2)) chi_{z1}(C^4e0) conj(chi_{z2}(C^4e0)) = chi_{z1}(e0eh2)chi_{z1}(Ce)^4 conj(...); "
       "conj(chi_z(e)) = chi_z(e)^5")
falsifier("D10.38", "psi(e) read as chi_{z1}(e) conj(chi_{z2}(e))", f1, 1e-6)
falsifier("D10.38", "chi_z(Ce)^2 in place of chi_z(Ce)^4", f2, 1e-6)

# ------------------------------------------------------------------ D10.39: exclusion indicator
viol = fal = 0; cnt = 0; nontriv = 0
while cnt < 4000:
    ps = random.sample(SMALLP[:22], 10)
    C = prod(ps[0:random.randint(0, 2)]); t = prod(ps[2:2 + random.randint(0, 1)])
    gps = ps[3:3 + random.randint(0, 3)]; g = prod(gps); e = prod([p for p in gps if random.random() < 0.5])
    w = prod(ps[6:6 + random.randint(0, 2)])
    e0 = prod([p for p in (factor(C) if C != ONE else []) if random.random() < 0.5])
    h2 = random.choice(LAT_SMALL)
    if not coprime(w, h2): continue                       # (w, gCth2) = 1
    # column n: random squarefree, often sharing primes with g, e, t, C, w
    npool = random.sample(SMALLP[:22], 3)
    n = prod(npool[:random.randint(1, 3)])
    gCtw = mul(mul(mul(g, C), t), w); tge = mul(t, div(g, e))
    core = chi1(n, mul(mul(e0, e), h2))*chi1(n, mul(mul(C, e), w))**4
    lhs = (1 if coprime(n, gCtw) else 0)*core
    rhs = (1 if coprime(n, tge) else 0)*core
    bad_ = (1 if coprime(n, t) else 0)*core
    viol = max(viol, abs(lhs - rhs)); fal = max(fal, abs(lhs - bad_)); nontriv += abs(lhs) > 0.5
    cnt += 1
report("D10.39", viol, 1e-12, f"(10.39) {cnt} random tuples with (g,Ct) = (w,gCth2) = 1, e | g, e0 | C, n squarefree "
       f"({nontriv} with nonzero left side): 1_(n,gCtw)=1 chi_n(e0eh2)chi_n(Cew)^4 = 1_(n,tg/e)=1 (same)")
falsifier("D10.39", "1_{(n,t)=1} (g/e dropped)", fal, 0.5)

# ------------------------------------------------------------------ D10.35: the multiplicity / containment step
def bump(x, a, b):
    x = np.asarray(x, float); out = np.zeros_like(x); m = (x > a) & (x < b)
    out[m] = np.exp(-1.0/((x[m] - a)*(b - x[m]))); return out
xi1 = lambda n: 1.0 + 0j if n == ONE else complex(chi(n, (-1, 0)))*complex(cubic3(n, (0, 1)))
u_, v_ = 1.0, 2.0; CPhi = 1.0; cst = max(1.0, 4*CPhi*v_**2)
Lp, Fp_, Sigp, Hp = 20.0, 7.0, 150.0, 10.0              # LF = 140 <= Sigma, H <= Sigma
Ustar = lambda x: bump(x, u_/2, 2*v_)*np.exp(1j*np.asarray(x, float))
SQ = primaries(400, sqfree=True)
def Nvec(A, B): return A*A - A*B + B*B
blocks = []
for C in SQ:
    for t in SQ:
        if N(C)*N(t) > v_*Lp or not coprime(C, t): continue
        for e0 in divs_sqfree(C): blocks.append((C, t, e0))
LATC = {}
worst = 0; maxmult = 0; ncont = 0; nbad_c = 0; tot_l = tot_r = 0.0; nterms = 0; taubad = 0
for (C, t, e0) in blocks:
    ell = Lp/(N(C)*N(t)); Y = cst*Sigp*Lp*Fp_*N(e0)/(Hp*N(C)**2); Ybad = Y/cst*max(1.0, CPhi*v_**2)
    w = Hp*N(C)/(N(e0)*Lp**2*Fp_)
    hmax = CPhi*v_**2*Lp**2*N(e0)/(Hp*N(C)**2)
    hs = [z for z in lattice(hmax) if z != (0, 0) and all(not divides(p, z) for p in (factor(t) if t != ONE else []))]
    fs = [f for f in SQ if Fp_ <= N(f) < 2*Fp_ and coprime(f, mul(C, t))]
    if not hs or not fs: continue
    if Y not in LATC:
        pts = lattice(Y); LATC[Y] = (pts, np.array([p[0] for p in pts], np.int64), np.array([p[1] for p in pts], np.int64),
                                     {p: i for i, p in enumerate(pts)})
    pts, A, B, pos = LATC[Y]
    cols = [n for n in SQ if u_/2*ell < N(n) < 2*v_*ell and coprime(n, t)]
    P = np.zeros(len(pts), complex)
    for n in cols:
        cf = muu(n)*al(n).conjugate()*xi1(n)*chi1(n, C)**4*chi1(n, e0)*complex(Ustar(N(n)/ell))
        if cf == 0: continue
        P += cf*np.conj(chi_vec(n, A, B)) if n != ONE else cf
    mult = {}
    for f in fs:
        f2 = mul(f, f)
        for h in hs:
            y = mul(h, f2); mult[y] = mult.get(y, 0) + 1
    lhs = 0.0
    for f in fs:
        f2 = mul(f, f)
        for h in hs:
            y = mul(h, f2); nterms += 1
            if y not in pos: ncont += 1; continue          # containment failure (must not happen)
            lhs += abs(P[pos[y]])**2
            nbad_c += N(y) > Ybad
    for y, c in mult.items():
        tau_y = sum(1 for f in SQ if N(f)**2 <= N(y) and divides(mul(f, f), y))
        taubad += c > tau_y
    mb = max(mult.values()); maxmult = max(maxmult, mb)
    rhs_in = float(np.sum(np.abs(P)**2))                  # sum over N(y) <= Y, where Phi(N(y)/Y) >= 1
    worst = max(worst, lhs/(mb*rhs_in) if rhs_in > 0 else (0 if lhs == 0 else 9e9))
    tot_l += w*lhs; tot_r += w*rhs_in
report("D10.35", ncont + taubad + max(0.0, worst - 1), 0, f"(10.35) {len(blocks)} blocks (C,t,e0), L={Lp:g}, F={Fp_:g}, "
       f"Sigma={Sigp:g}, H={Hp:g}, {nterms} pairs (f,h): every hf^2 lies in N(y) <= Y_{{C,e0}}; multiplicity <= "
       f"#{{f: f^2 | y}} <= tau(y) (max {maxmult}); blockwise sum_(f,h)|P(hf^2)|^2 <= m_max sum_(N(y)<=Y)|P(y)|^2, "
       f"worst ratio {worst:.3f}; weighted totals {tot_l:.4g} <= {tot_r:.4g}")
falsifier("D10.35", "containment with c_* = max(1, C_Phi v^2)", nbad_c, 0)

# ------------------------------------------------------------------ D10.41: exponent arithmetic, E <= A(W), source orders
m_ = sp.symbols("m", integer=True, nonnegative=True)
J1 = 2*m_ + 4
ok = sp.expand(2*J1 + 4 - (4*m_ + 12)) == 0
def energy_pair(H, Xv, F, Phi):
    """E(H, X, F) of (7.4) with 0 < N(k) <= H, and A(W) with weight Phi(N(k)/H) over all k (cut at Phi < 1e-18)."""
    ns = [n for n in primaries(int(2*Xv) + 2, sqfree=True) if Xv < N(n) < 2*Xv]
    fs = [f for f in primaries(int(2*F) + 1, sqfree=True) if F <= N(f) < 2*F]
    ks = lattice(H*45); A_ = np.array([k[0] for k in ks]); B_ = np.array([k[1] for k in ks])
    Nk = Nvec(A_, B_); wk = Phi(Nk/H); inE = (Nk > 0) & (Nk <= H)
    E = Aw = 0.0
    for f in fs:
        vv = np.zeros(len(ks), complex)
        for n in ns:
            cf = al(n).conjugate()**2*gamma_j_fast(2, n)*complex(cubic3(n, (0, 1)))*chi1(n, f)**4*bump(N(n)/Xv, 1.0, 2.0)
            if cf != 0: vv += cf*chi_vec(n, A_, B_)
        E += float(np.sum(np.abs(vv[inE])**2)); Aw += float(np.sum(wk*np.abs(vv)**2))
    return E/(Xv*F), Aw/(Xv*F)
PhiT = lambda t: np.exp(1.0 - np.asarray(t, float))          # nonnegative, >= 1 on [0, 1]
lines = []; bad = 0
for (H, Xv, F) in [(10, 150.0, 1.0), (30, 200.0, 7.0), (5, 120.0, 13.0)]:
    E, Aw = energy_pair(H, Xv, F, PhiT); bad += E > Aw*(1 + 1e-12); lines.append(f"H={H},L={Xv:g},F={F:g}: E={E:.4g} <= A={Aw:.4g}")
p = src.index("\\begin{proposition}\\label{prop:transfer}"); pt = src[p:src.index("\\end{proposition}", p)]
l1 = src.index("\\begin{lemma}\\label{lem:first-transfer}"); l1t = src[l1:src.index("\\end{lemma}", l1)]
l3 = src.index("\\begin{lemma}\\label{lem:second-transfer}"); l3t = src[l3:src.index("\\end{lemma}", l3)]
srcok = ("C^{4m+12}(I)" in pt and "[u/16,4v]" in pt and "C^{2j+4}(I)" in l1t and "C^{2m+4}(I_*)" in l3t)
report("D10.41", (0 if ok else 1) + bad + (0 if srcok else 1), 0, "(10.41) 2(2m+4)+4 = 4m+12 (sympy); E(H,L,F) <= A(W) "
       "with Phi(t) = e^{1-t} >= 1 on [0,1]: " + "; ".join(lines) + f"; OAI1112 Lemma 7.1 C^(2j+4), Lemma 7.3 C^(2m+4), "
       f"Prop. 5.4 C^(4m+12) and [u/16,4v] in the source: {srcok}")

# ------------------------------------------------------------------ D10.27: zeta column of the table (10.27)
zg = open(ZG).read()
checks = {
    "mu_K(n) (original coefficients, A_u)": "A_u(D)&=\\sum_n\\muk(n)\\chi_n(u)W(\\N(n)/D)" in zg,
    "conj(alpha) gamma_2 (cubic coefficients)": "$a_\\xi(n)=\\overline{\\alpha(n)}\\gamma_2(n)\\xi(n)$" in zg,
    "conj(alpha(b))^3 (cubes, eq:T)": "\\overline{\\alpha(b_1)}^3\\Psi(b_1)^3" in zg,
    "conj(alpha(h))^3 (cube inversion)": "\\muk(h)\\overline{\\alpha(h)}^3\\Psi(h)^3" in zg,
    "one zbar-derivative (zeta guide numerics)": "one $\\bar z$-derivative" in zg,
    "d_zbar Theta at z=0 (OAI1112)": "\\partial_{\\bar z}\\Theta_\\Psi(z,v)\\Bigr|_{z=0}" in src,
    "alpha(ell) in the dual sum (OAI1112)": "\\frac{d(\\ell)\\alpha(\\ell)}{\\sqrt{\\NK(\\ell)}}" in src,
}
miss = [k for k, v in checks.items() if not v]
report("D10.27", len(miss), 0, "(10.27) zeta column found in the sources: " + "; ".join(f"{k}: {v}" for k, v in checks.items()))

# ------------------------------------------------------------------ S10.32 size check: zero-frequency constant
print("      S10.32 size check [(10.32): Z << H ||W||^2]: c(L,F) = (1/(LF)) sum_f sum_(C,f)=1 |W(N(C)/L)|^2 prod(1-1/N(pi)), ||W||_inf = 1")
Wn = lambda x: math.exp(4 - 1/((x - 1)*(2 - x))) if 1 < x < 2 else 0.0
SQbig = primaries(6000, sqfree=True)
cmax = 0; rowsz = []
for L in (100.0, 300.0, 1000.0, 3000.0):
    cs = []
    for F in (1.0, 7.0, 30.0, 100.0):
        fs = [f for f in SQbig if F <= N(f) < 2*F]; Cs = [c for c in SQbig if L < N(c) < 2*L]
        s = 0.0
        for f in fs:
            for c in Cs:
                if not coprime(c, f): continue
                s += Wn(N(c)/L)**2*math.prod(1 - 1/N(p) for p in factor(c))
        cs.append(s/(L*F)); cmax = max(cmax, s/(L*F))
    rowsz.append(f"L={L:g}: " + ", ".join(f"{c:.3f}" for c in cs))
print("          F = 1, 7, 30, 100:  " + ";  ".join(rowsz))
report("S10.32", max(0.0, cmax - 1.0), 0, f"(10.32) size check: Z/(H Phihat(0)||W||^2) <= {cmax:.3f} for L <= 3000, F <= 100 "
       "(bounded, no growth); Z <= Sigma||W||^2 then uses H <= Sigma", note="class A: the constant is O_I(1); sanity only")

# ------------------------------------------------------------------ S10.34 size check: kernel derivatives uniform in h
z1, z2, cc = sp.symbols("z1 z2 c", positive=True)
Wsym = lambda z: sp.exp(-1/((z - 1)*(2 - z)))
Cphi_ = 1.0
Phh = lambda s: sp.exp(-1/(Cphi_ - s))                     # smooth, vanishes for s >= C_Phi (model of Phihat)
K = (z1*z2)**sp.Rational(-1, 2)*Wsym(z1)*Wsym(z2)*Phh(cc/(z1*z2))
grid = np.linspace(1.002, 1.998, 70); Z1, Z2 = np.meshgrid(grid, grid, indexing="ij")
Wnorm = {}
for q in range(3):
    Wnorm[q] = max(float(np.max(np.abs(sp.lambdify(z1, sp.diff(Wsym(z1), z1, i), "numpy")(grid)))) for i in range(q + 1))
ratios = {}
cgrid = np.linspace(0.0, Cphi_*4.0*0.999, 25)              # c/(z1 z2) < C_Phi requires c < 4 C_Phi on I^2
for q in range(3):
    worst = 0.0
    fns = [sp.lambdify((z1, z2, cc), sp.diff(K, z1, a, z2, b), "numpy") for a in range(q + 1) for b in range(q + 1 - a)]
    for c in cgrid:
        mask = (c/(Z1*Z2) < Cphi_*(1 - 1e-9))
        nm = 0.0
        for fn in fns:
            with np.errstate(all="ignore"):
                val = np.where(mask, fn(Z1, Z2, c), 0.0)
            nm = max(nm, float(np.nanmax(np.abs(val))))
        worst = max(worst, nm)
    ratios[q] = worst/Wnorm[q]**2
report("S10.34", 0 if all(np.isfinite(list(ratios.values()))) else 1, 0, "(10.34) size check: sup over c in [0, 4C_Phi) of "
       "||K_c||_{C^q(I^2)}/||W||_{C^q(I)}^2 = " + ", ".join(f"q={q}: {r:.3g}" for q, r in ratios.items()) +
       " (I = [1,2], W and Phihat smooth bumps); K_c = 0 for c >= 4 C_Phi", note="class A: uniformity in (C, e0, h); sanity only")
print("ALL DONE")
