"""t05: rows, labels and moments (guide L1187-1210, L1291-1322, L1353-1428).

Claims tested
  C45 L1298-1309 the label set X_u (model in a finite group): finite, closed under conjugation, contains the numerator
                 chi(u), the denominator psi conj chi(u) and their conjugates.
  C46 L1313-1321 the two witnesses for a CM member are co-oriented (abar with conj chi, or alpha with chi).
  C47 L1418-1420 int_C alpha(z)^e f(|z|^2) dz = 0 (e != 0), and lattice/coset sums sum_z alpha(z)^e f(N z/X) have no
                 volume main term (they are O(X^{-A}) for smooth f).
  C48 L1202      Prop. cm-completed: OAI78's exponent O/2 + max(M-O, 2M-O-N) (lem:unmarked-completed-moment, L3155-3160)
                 with N = 1 is <= max(M, 2M-1), i.e. the bound Z^eps (Q + Q^2/Z).
  C49 L1202      sanity numbers for sum_{N(m)<=Q} |B_m|^2 against Q + Q^2/Z (CM weight abar^2), no criterion.
  C50 L1389-1390 (MS) of Part I is the special case eps=1, nu=1, no slots, row range <= D^{11/10} of eq:cm-inverse:
                 the conditions r+2z <= m-c_1, 2r+8z <= 3m-c_2 hold with r=1, m=1+theta, z=0, c_1=theta, c_2=1, and
                 Z^{m+eps} normalization reproduces D^{2+theta+eps}.
Run time: about 40 s.
"""
import os, sys, math, cmath, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmcommon import *
from cmreport import claim, cited, info, untested
import numpy as np

# ---------------------------------------------------------------- C45, C46 (finite-group model)
Gm = [(i, j) for i in range(6) for j in range(2)]                 # model of Theta_0 = <eta, T^>: Z/6 x Z/2
neg = lambda th: ((-th[0]) % 6, (-th[1]) % 2)
eta = (1, 0)
def conjc(ch):
    k, th, s = ch
    return (-k, neg(th), -s)
Xu = {(0, th, s) for th in Gm for s in (1, -1)} | {(1, th, -1) for th in Gm} | {(-1, th, 1) for th in Gm}
num = (0, (0, 0), 1); den = (1, eta, -1)                            # chi(u);  psi conj chi(u) = abar eta conj chi(u)
need = [num, den, conjc(num), conjc(den)]
ok45 = all(x in Xu for x in need) and all(conjc(x) in Xu for x in Xu) and len(Xu) == 4*len(Gm)
Xanti = {(0, th, s) for th in Gm for s in (1, -1)} | {(1, th, 1) for th in Gm} | {(-1, th, -1) for th in Gm}
Xhalf = {(0, th, s) for th in Gm for s in (1, -1)} | {(1, th, -1) for th in Gm}
fal45 = (all(x in Xanti for x in need)) or all(conjc(x) in Xhalf for x in Xhalf)
claim("C45", 1299, "label set X_u", "E",
      "X_u = {theta chi(u)^+-1} u {abar theta conj chi(u), alpha theta chi(u)}: finite, conjugation-closed, contains numerator, "
      "denominator and conjugates (model Theta_0 = Z/6 x Z/2)", ok45,
      f"|X_u| = {len(Xu)} = 4|Theta_0|; falsifiers: anti-oriented set misses the denominator, half set not closed",
      falsifier_failed=not fal45)
cm = [x for x in Xu if x[0] != 0]
co = all(x[0]*x[2] == -1 for x in cm)
# witnesses: M_r coefficient mu(n) abar(n) theta(n) conj chi_n(u) -> (k=1, sigma=-1); S_m: abar(l)theta(l)conj chi_l(u)
wit = [(1, eta, -1), (1, eta, -1), conjc((1, eta, -1))]
claim("C46", 1320, "co-orientation", "E", "all CM members of X_u and both witnesses have k*sigma = -1 (co-oriented)",
      co and all(w[0]*w[2] == -1 for w in wit), f"{len(cm)} CM members; falsifier: anti-oriented set has k*sigma = +1",
      falsifier_failed=all(x[0]*x[2] == 1 for x in Xanti if x[0] != 0))

# ---------------------------------------------------------------- C47
def bump(r):
    out = np.zeros_like(r)
    m = (r > 1) & (r < 2)
    out[m] = np.exp(-1/((r[m] - 1)*(2 - r[m])))
    return out
h = 0.004; xs = np.arange(-1.5, 1.5, h); XX, YY = np.meshgrid(xs, xs, indexing="ij")
Z = XX + 1j*YY; R2 = np.abs(Z)**2
with np.errstate(divide="ignore", invalid="ignore"):
    AL = np.where(R2 > 0, Z/np.abs(Z), 0)
F = bump(R2)
ints = {e: abs(np.sum(AL**e*F)*h*h) for e in (1, 2, -1, 3)}
int0 = abs(np.sum(F)*h*h)
# lattice cosets a + mO, alpha^e f(N z / X)
w = cmath.exp(2j*math.pi/3)
def ring(r):
    return np.exp(-(r - 1.5)**2/0.02)
def coset_sum(a, m, e, X, prof):
    mc = m[0] + m[1]*w; ac = a[0] + a[1]*w
    R = int(2*math.sqrt(3*X)/abs(mc)) + 4
    A, B = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    zc = ac + mc*(A + B*w)
    n = np.abs(zc)**2
    f = prof(n/X)
    al = np.zeros_like(zc); nz = n > 0
    al[nz] = zc[nz]/np.sqrt(n[nz])
    return abs(np.sum(al**e*f)), np.sum(f)
cos = [((1, 0), (3, 0)), ((2, 1), (5, 0)), ((1, 0), (2, 3)), ((0, 0), (1, 0))]
worst_ring = 0.0; bump_by_X = {}; vol = []
for (a, m) in cos:
    for e in (1, 2, -1):
        for X in (3e3, 3e4, 3e5):
            s1, s0 = coset_sum(a, m, e, X, bump); bump_by_X[X] = max(bump_by_X.get(X, 0), s1/s0); vol.append(s0)
dec = bump_by_X[3e3] > bump_by_X[3e4] > bump_by_X[3e5]
ok47 = max(ints.values()) < 1e-10*int0 and dec and bump_by_X[3e5] < 1e-9
claim("C47", 1420, "lattice main terms vanish", "E",
      "int alpha^e f(|z|^2) = 0 for e=1,2,-1,3; coset sums sum_{z in a+mO} alpha(z)^e f(N z/X) have no volume term", ok47,
      f"max |int|/int(e=0) {max(ints.values())/int0:.1e}; C_c bump profile, max |sum|/sum(e=0) over 4 cosets (incl. O), "
      f"e=1,2,-1: " +
      ", ".join(f"X={X:g}: {v:.1e}" for X, v in sorted(bump_by_X.items())) +
      f" (decreasing: O(X^-A) is asymptotic); falsifier e=0: int = {int0:.3f}, coset sums ~ volume (up to {max(vol):.3g})",
      falsifier_failed=int0 > 0.01)

# ---------------------------------------------------------------- C48
grid = [Fr(i, 20) for i in range(0, 61)]
ok48 = True; fal48 = False
for M in grid:
    for O in [o for o in grid if o <= 2*M]:
        lhs = O/2 + max(M - O, 2*M - O - 1)
        ok48 &= lhs <= max(M, 2*M - 1)
        fal48 |= (2*O + max(M - O, 2*M - O - 1)) > max(M, 2*M - 1)
claim("C48", 1202, "Prop. cm-completed exponent", "I",
      "O/2 + max(M-O, 2M-O-N) <= max(M, 2M-N) at N=1 for 0<=O<=2M (summing O-dyads gives Z^eps(Q+Q^2/Z))", ok48,
      "exact Fractions on a grid M in [0,3]; falsifier weight 2O instead of O/2 exceeds", falsifier_failed=fal48)

# ---------------------------------------------------------------- C49 sanity: completed rows with abar^2
Zs = 1500.0
def Vb(r):
    return math.exp(4 - 1/((r - 1)*(2 - r))) if 1 < r < 2 else 0.0
pr = [p for p in PR_ODD if N(p) <= 2*Zs]
g2p = {p: gamma_j(2, p) for p in pr}
def gamma2_sf(c):
    f = factor(c); g = 1 + 0j
    for p in f:
        g *= g2p[p]
    for i in range(len(f)):
        for j in range(i + 1, len(f)):
            g *= (chi(f[i], f[j])*chi(f[j], f[i]))**2
    return g
# check twisted multiplicativity (OAI78 L1043) against brute force on a few composites
tw = max(abs(gamma2_sf(c) - gamma_j(2, c)) for c in [mul(pr[1], pr[3]), mul(pr[2], pr[6]), mul((-5, 0), pr[0])])
terms = []
for n in primaries(int((2*Zs)**(1/3)) + 1):
    Nn3 = N(n)**3
    for c in primaries(int(2*Zs/Nn3), sqfree=True):
        wv = Vb(N(c)*Nn3/Zs)
        if wv:
            terms.append((c, n, wv))
R = 70
rowsm = sorted([(x, y) for x in range(-R, R + 1) for y in range(-R, R + 1) if 0 < N((x, y)) <= 3000], key=N)
symtab = {}
import zw
def chivec(p):
    """chi_p(m) for all rows m, via a table on canonical residues mod p (values from Euler's criterion)."""
    if p not in symtab:
        tab = {}
        for r in zw.residues(p):
            j = zw.sym6_prime(r, p)
            tab[r] = 0 if j is None else ZETA[j]
        symtab[p] = np.array([tab[zw.canon_mod(m, p)] for m in rowsm], complex)
    return symtab[p]
B2 = np.zeros(len(rowsm), complex); B1 = np.zeros(len(rowsm), complex)
for c, n, wv in terms:
    v = np.ones(len(rowsm), complex)
    for p in factor(c):
        v *= chivec(p)
    for p in factor(n):
        v *= chivec(p)**3
    A = mul(c, pw(n, 3)); g2 = gamma2_sf(c)
    base = g2/(math.sqrt(N(c))*N(n))*wv
    B2 += base*abar(A)**2*v
    B1 += base*abar(A)*v
norms = np.array([N(m) for m in rowsm])
txt = []
for Q in (30, 100, 300, 1000, 3000):
    sel = norms <= Q
    s2 = np.sum(np.abs(B2[sel])**2); s1 = np.sum(np.abs(B1[sel])**2)
    txt.append(f"Q={Q}: CM {s2/(Q + Q*Q/Zs):.3f}, OAI {s1/(Q + Q*Q/Zs):.3f}")
info("C49", 1202, "Prop. cm-completed (sanity)", "A",
     f"sum_{{N(m)<=Q}}|B_m|^2/(Q+Q^2/Z), Z={Zs:g}, {len(terms)} completed terms, abar^2 (CM) vs abar (OAI78)",
     "; ".join(txt) + f"; twisted multiplicativity of gamma2 vs brute force: {tw:.1e}")

# ---------------------------------------------------------------- C50
ok50 = True; fal50 = False
for th in [Fr(k, 100) for k in range(1, 11)]:
    r, m, z, c1, c2 = Fr(1), 1 + th, Fr(0), th, Fr(1)
    ok50 &= (r + 2*z <= m - c1) and (2*r + 8*z <= 3*m - c2) and c1 > 0 and c2 > 0 and m <= Fr(11, 10)
    # normalization: M(D,H) = sum |A_u|^2 = D * sum |M_u|^2 << D * Z^{m+eps} = D^{2+theta+eps} with Z = D
    ok50 &= (1 + m == 2 + th)
    fal50 |= (r + 2*z <= m - 2*th)          # falsifier: a margin c_1 > theta is NOT available
claim("C50", 1389, "(MS) as a special case", "I",
      "(MS) [theta in (0,1/10]] = eq:cm-inverse with eps=1, nu=1, no slots, r=1, m=1+theta, c_1=theta, c_2=1; m <= 11/10",
      ok50, "exact; note c_1 must be taken <= theta (fixed for each fixed theta); falsifier c_1 = 2 theta fails",
      falsifier_failed=not fal50)

# ---------------------------------------------------------------- C103 sanity: CM inverse moment without slots
# sum_{0<N(u)<=H} |D^{-1/2} sum_n mu(n) alpha(n)^eps chi_n(u)^eps W(N(n)/D)|^2 / H  for D = 400, H = D^{1.1}
Dd = 400.0; Hh = Dd**1.1
cols = [n for n in primaries(int(2*Dd), sqfree=True) if Vb(N(n)/Dd) > 0]
Rr = int(math.sqrt(Hh)*1.2) + 2
urows = [(x, y) for x in range(-Rr, Rr + 1) for y in range(-Rr, Rr + 1) if 0 < N((x, y)) <= Hh]
rowsm = urows; symtab.clear()
tot = {"eps=-1 (co-oriented, CM witness)": 0.0, "eps=+1 (co-oriented, (MS))": 0.0, "anti-oriented alpha*conj chi": 0.0}
acc = {k: np.zeros(len(urows), complex) for k in tot}
for n in cols:
    v = np.ones(len(urows), complex)
    for p in factor(n):
        v *= chivec(p)
    w = mu(n)*Vb(N(n)/Dd)
    acc["eps=-1 (co-oriented, CM witness)"] += w*abar(n)*np.conj(v)
    acc["eps=+1 (co-oriented, (MS))"] += w*alpha(n)*v
    acc["anti-oriented alpha*conj chi"] += w*alpha(n)*np.conj(v)
info("C103", 1363, "eq:cm-inverse (sanity, no slots)", "A",
     f"sum_{{N(u)<=H}}|M_u|^2/H with D={Dd:g}, H=D^1.1={Hh:.0f}, {len(cols)} columns, {len(urows)} rows",
     "; ".join(f"{k}: {np.sum(np.abs(a)**2)/Dd/Hh:.3f}" for k, a in acc.items()) +
     " (all O(1): the anti-oriented family is not large; the guide only says the proof does not extend to it)")
