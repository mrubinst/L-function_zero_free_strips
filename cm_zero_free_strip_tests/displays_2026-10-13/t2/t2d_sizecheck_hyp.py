"""t2d: size check (sanity only) for the hypothesis eq:b2-hyp (7.25) in the ranges eq:b2-ranges (7.26).

eq:b2-hyp says E(calH, X, F; xi, W) << ||W||^2 D^eps X F for X = D/(BF), 0 < calH <= C D^2/(H B^2), H = D^{1+theta}.
It is an analytic bound with unspecified constants (proved via eq:finalenergy of Section 10); no computation proves it.
Here we compute E literally from eq:energy (CM coefficients a_xi = conj(alpha)^2 gamma_2 xi, two characters xi) at
theta = 1/10, C = 1, for several (D, B, F), and print E/(XF ||W||_inf^2).  The check passes if every ratio is <= 1
(arbitrary loose criterion; the point is that the ratio does not grow with D).
Run: python3 -B t2d_sizecheck_hyp.py   (about 20 s)
"""
import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2lib import *
T0 = time.time()
theta = 0.1
def Wb(y): return math.exp(-1/((y - 1)*(2 - y))) if 1 < y < 2 else 0.0
Wsup = max(Wb(1 + i/10000) for i in range(1, 10000))
SQ = primaries(9000, sqfree=True)
XIS = [TRIV, 7]

def energy(Hr, Xv, F, k_xi):
    ns = [n for n in SQ if Xv < N(n) < 2*Xv]
    fs = [f for f in SQ if F <= N(f) < 2*F]
    R = int(2*math.sqrt(Hr)) + 3
    A, B = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing='ij'); A = A.ravel(); B = B.ravel()
    NN = A*A - A*B + B*B; keep = (NN > 0) & (NN <= Hr); A, B = A[keep], B[keep]
    coef = {n: Acm(n)*CHARS[k_xi, cls(n)]*Wb(N(n)/Xv) for n in ns}
    chik = {n: chi_vec(n, A, B) for n in ns}
    tot = 0.0
    for f in fs:
        v = np.zeros(len(A), complex)
        for n in ns:
            c = coef[n]*complex(chi(n, f))**4
            if c != 0: v += c*chik[n]
        tot += float(np.sum(np.abs(v)**2))
    return tot/(Xv*F), len(ns), len(fs), len(A)          # this is E of eq:energy (already divided by XF once)

rows = []
for D in (500.0, 1500.0, 4000.0):
    H = D**(1 + theta)
    for (B, F) in [(1, 1), (1, 7), (3, 1), (2, 13)]:
        Xv = D/(B*F); Hr = D*D/(H*B*B)
        for k in XIS:
            E, nn, nf, nk = energy(Hr, Xv, F, k)
            r = E/(Xv*F*Wsup**2)
            rows.append((D, B, F, k, Hr, Xv, r, nn, nf, nk))
            print(f"  D={D:g} B={B} F={F} xi#{k}: calH={Hr:.1f} ({nk} rows), X={Xv:.1f} ({nn} columns), {nf} f: E = {E:.4g}, E/(XF||W||^2) = {r:.4f}")
worst = max(r[6] for r in rows)
report("eq:b2-hyp", 0.0 if worst <= 1 else worst, 0.0, f"size check only: E/(XF||W||^2) <= {worst:.3f} over {len(rows)} cases in the ranges eq:b2-ranges "
       f"(theta = 1/10, C = 1, D = 500, 1500, 4000)")
print(f"t2d done in {time.time() - T0:.1f} s")
