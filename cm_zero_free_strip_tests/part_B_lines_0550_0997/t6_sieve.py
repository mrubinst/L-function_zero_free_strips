"""Section 'The quadratic average and the completed bound' (lines 810-849): claims B62-B67.

B62  eq:qls (cited, Goldmakher-Louvel): sanity ratio LHS/((K+U) sum|beta|^2), K=U=400, rows squarefree primary coprime
     to 6, columns squarefree primary (even columns included), six coefficient choices.  Pass criterion (arbitrary): < 2.
B63  'square columns would give repeated values 1': beta = 1 on squares b^2 gives (b^2/k)_2 = 1 for (b,k) = 1, and the
     ratio grows like sqrt(U) (no cancellation); the squarefree control stays bounded.
B64  eq:completed (analytic): sanity ratio sum_{0<N(k)<=H}|T(X;k,f)|^2 / (H + H^2 N(f)/X) with the CM T of eq:T
     (conj(alpha)^2 gamma_2, conj(alpha(b))^6, all lattice rows k).  Pass criterion (arbitrary): < 5.
B65  every k != 0 is u s v0^2 with u a unit and s squarefree (all k with N(k) <= 3000).
B66  chi_n(v0^2) chi_n(f)^4 = chi_n(f v0^2)^4, and hence T(X; u s v0^2, f) = T(X; u s, f v0^2) exactly.
     Falsifier: chi_n(f v0)^4.
B67  sum_{v0 != 0} N(v0)^{-2} = 6 zeta(2) L(2, chi_-3) < oo (partial sums), and the bookkeeping
     (H/N(v0)^2) + (H/N(v0)^2)^2 N(f v0^2)/X = N(v0)^{-2} (H + H^2 N(f)/X)  (sympy).
Run time: about 1 minute.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import *
import sympy as sp, mpmath as mp
random.seed(849)
LAMc = cx(LAM)

# ---------------- B62, B63
K = U = 400
rows = primaries(K, sqfree=True)                                   # squarefree primary, coprime to 6
def sqfree_even_cols(Ub):
    out = []
    R = int(2*math.sqrt(Ub)) + 3
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            n = (a, b)
            if not primary(n) or not (0 < N(n) <= Ub): continue
            m = n; tw = 0
            while divides((2, 0), m): m = div(m, (-2, 0)); tw += 1
            if tw > 1: continue
            if m != (1, 0) and mu(m) == 0: continue
            out.append(n)
    return out
cols = sqfree_even_cols(U)
A = np.array([n[0] for n in cols]); B = np.array([n[1] for n in cols])
QS = np.array([chi_vec(k, A, B, 3) for k in rows])                 # (n/k)_2 = chi_k(n)^3
worst = 0; rng = np.random.default_rng(1)
for trial in range(6):
    beta = np.ones(len(cols)) if trial == 0 else rng.normal(size=len(cols)) + 1j*rng.normal(size=len(cols))
    lhs = np.sum(np.abs(QS @ beta)**2)
    worst = max(worst, lhs/((K + U)*np.sum(np.abs(beta)**2)))
neven = sum(1 for n in cols if divides((2, 0), n))
report("B62", worst, 2.0, f"K=U=400: {len(rows)} rows, {len(cols)} squarefree columns ({neven} even); max LHS/((K+U)sum|beta|^2) = {worst:.3f} over 6 choices",
       note="cited bound with unspecified constant; sanity only")
# B63: beta = 1 on squares b^2 (b primary coprime to 6) versus beta = 1 on squarefree columns, with K = U growing
ratios_sq = []; ratios_sf = []; ones = True
for Ub in [400, 1600, 6400]:
    rws = primaries(Ub, sqfree=True)
    bsq = primaries(int(math.sqrt(Ub)))
    sqc = [pw(b, 2) for b in bsq]
    As = np.array([n[0] for n in sqc]); Bs = np.array([n[1] for n in sqc])
    vals = np.array([chi_vec(k, As, Bs, 3) for k in rws])
    nz = vals[np.abs(vals) > 1e-9]
    ones &= bool(np.all(np.abs(nz - 1) < 1e-9))
    ratios_sq.append(np.sum(np.abs(vals.sum(axis=1))**2)/((2*Ub)*len(sqc)))
    sf = sqfree_even_cols(Ub)
    Af = np.array([n[0] for n in sf]); Bf = np.array([n[1] for n in sf])
    tot = 0.0
    for k in rws:
        tot += abs(np.sum(chi_vec(k, Af, Bf, 3)))**2
    ratios_sf.append(tot/((2*Ub)*len(sf)))
grows = ratios_sq[2] > 1.5*ratios_sq[1] and ratios_sq[1] > 1.5*ratios_sq[0]
report("B63", 0 if (ones and grows) else 1, 0, f"(b^2/k)_2 = 1 whenever (b,k)=1: {ones}; K=U=400,1600,6400: ratio LHS/((K+U)sum|beta|^2) with beta=1 on squares "
       f"{', '.join(f'{r:.2f}' for r in ratios_sq)} (grows like sqrt U); with beta=1 on squarefree columns {', '.join(f'{r:.3f}' for r in ratios_sf)} (bounded)")
# ---------------- B64: completed estimate sanity with the CM T
def Wb(y): return math.exp(4 - 1/((y - 1)*(2 - y))) if 1 < y < 2 else 0.0      # ||W||_inf = 1
X = 1400.0
sqn = primaries(int(2*X), sqfree=True)
bset = primaries(int((2*X)**(1/3)) + 2)
def xi(n): return complex(cubic3(n, (0, 1)))
g2 = {n: gamma_j_fast(2, n) for n in sqn}
pairs = []
for b in bset:
    for n in sqn:
        y = N(n)*N(b)**3/X
        w = Wb(y)
        if w == 0: continue
        pairs.append((n, b, alpha_c(n).conjugate()**2*g2[n]*xi(n)*alpha_c(b).conjugate()**6*xi(b)**3*math.sqrt(N(b))*w/math.sqrt(X)))
Hmax = 600
krows = [k for k in lattice(Hmax) if k != (0, 0)]
KA = np.array([k[0] for k in krows]); KB = np.array([k[1] for k in krows]); KN = np.array([N(k) for k in krows])
primes_used = sorted(set(p for n, b, c in pairs for p in factor(n) + factor(b)))
J = {p: sym6_index_vec(KA, KB, p) for p in primes_used}
def Tvals(f):
    tot = np.zeros(len(krows), complex)
    for n, b, c in pairs:
        e = {}
        for p in factor(n): e[p] = e.get(p, 0) + 1
        for p in factor(b): e[p] = e.get(p, 0) + 3
        fac = complex(chi(n, f))**4*complex(chi(b, f))**12          # chi_{n b^3}(f)^4, zero convention
        if fac == 0: continue
        idx = np.zeros(len(krows), np.int64); zero = np.zeros(len(krows), bool)
        for p, m in e.items():
            jp = J[p]; zero |= jp < 0; idx += m*np.where(jp < 0, 0, jp)
        tot += np.where(zero, 0, c*fac*ZETA6[idx % 6])
    return tot
P7 = [p for p in PR_ODD if N(p) == 7][0]; P13 = [p for p in PR_ODD if N(p) == 13][0]
worst = 0; rowsout = []
for f in [(1, 0), P7, mul(P7, P13)]:
    T = Tvals(f)
    for H in [20, 150, 600]:
        S = np.sum(np.abs(T[KN <= H])**2)
        r = S/(H + H*H*N(f)/X); worst = max(worst, r); rowsout.append(f"N(f)={N(f)},H={H}: {r:.3f}")
report("B64", worst, 5.0, f"X={X:g}, {len(pairs)} (n,b) terms ({sum(1 for n, b, c in pairs if b != (1, 0))} with b!=1), all lattice rows: "
       "ratio sum|T|^2/(H+H^2N(f)/X): " + "; ".join(rowsout), note="analytic bound with D^eps and unspecified constant; sanity only")

# ---------------- B65
PRIMES_ALL = [LAM, (-2, 0)] + list(PR_ODD)
def usv(k):
    rest = k; s = (1, 0); v = (1, 0)
    for p in PRIMES_ALL:
        if N(rest) == 1: break
        if N(p) > N(rest): continue
        e = 0
        while divides(p, rest): rest = div(rest, p); e += 1
        s = mul(s, pw(p, e % 2)); v = mul(v, pw(p, e//2))
    return rest, s, v
bad = 0; cnt = 0; nonprim = 0
for k in lattice(3000):
    if k == (0, 0): continue
    u, s_, v = usv(k); cnt += 1
    bad += (N(u) != 1) or (mul(mul(u, s_), pw(v, 2)) != k)
    fs = factor(s_) if (N(s_) > 1 and not divides(LAM, s_) and not divides((2, 0), s_)) else []
    nonprim += divides(LAM, s_) or divides((2, 0), s_)
report("B65", bad, 0, f"all {cnt} nonzero k with N(k)<=3000 factor as u s v0^2 (s squarefree; {nonprim} have lambda or 2 in s, which the row sum allows)")

# ---------------- B66
err = 0; f66 = 0; sq = primaries(600, sqfree=True)
for _ in range(400):
    n = random.choice(sq); v0 = (random.randint(-15, 15), random.randint(-15, 15)); f = (random.randint(-15, 15), random.randint(-15, 15))
    if v0 == (0, 0) or f == (0, 0): continue
    lhs = chi(n, pw(v0, 2))*chi(n, f)**4
    err = max(err, abs(lhs - chi(n, mul(f, pw(v0, 2)))**4)); f66 = max(f66, abs(lhs - chi(n, mul(f, v0))**4))
# T(X; u s v0^2, f) = T(X; u s, f v0^2) exactly (with the CM T, all (n, b) terms above)
def T_single(k, f):
    return sum(c*chi(n, k)*chi(b, k)**3*chi(n, f)**4*chi(b, f)**12 for n, b, c in pairs)
eT = 0
for (u, s_, v0) in [((0, 1), P13, P7), ((-1, -1), (5, 3), (2, 1)), ((1, 0), (1, 0), LAM), ((1, 1), P7, mul(P13, (-2, 0)))]:
    for f in [(1, 0), P7]:
        k = mul(mul(u, s_), pw(v0, 2))
        eT = max(eT, abs(T_single(k, f) - T_single(mul(u, s_), mul(f, pw(v0, 2)))))
report("B66", max(err, eT), 1e-10, f"400 random (n, v0, f) incl. non-coprime cases; T(X;u s v0^2,f) = T(X;u s,f v0^2) for 8 cases (err {eT:.1e})")
falsifier("B66", "chi_n(f v0)^4", f66, 1e-3)

# ---------------- B67
zK2 = 6*mp.zeta(2)*mp.nsum(lambda k: 1/(3*k + 1)**2 - 1/(3*k + 2)**2, [0, mp.inf])
pts = []
for R in [100, 1000, 10000]:
    Sm = sum(1.0/N(v)**2 for v in lattice(R) if v != (0, 0)); pts.append((R, Sm, float(zK2) - Sm))
tail_ok = all(0 < d < 15.0/R for R, Sm, d in pts)
Hs, Nf, Xs, Nv = sp.symbols('H Nf X Nv', positive=True)
ok = sp.simplify((Hs/Nv**2 + (Hs/Nv**2)**2*Nf*Nv**2/Xs) - (Hs + Hs**2*Nf/Xs)/Nv**2) == 0
report("B67", 0 if (tail_ok and ok) else 1, 0, f"6 zeta(2) L(2,chi_-3) = {float(zK2):.6f}; partial sums over N(v0)<=R: " +
       ", ".join(f"R={R}: {Sm:.6f} (tail {d:.1e})" for R, Sm, d in pts) + "; bookkeeping identity (sympy)")
