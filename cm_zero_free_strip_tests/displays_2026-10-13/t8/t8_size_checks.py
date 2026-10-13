"""t8_size_checks.py -- finite size checks for three asymptotic bounds (kind A) in Appendices E and G of
cm-eleven-twelfths-guide.tex.  These bounds have unspecified constants and cannot be proved by computation; the
checks only confirm that the computed sizes are consistent with them (and they would expose a gross error, e.g. a
missing angular factor, which makes the u = 1 row or the sum of A_1 grow like D).

  SZ1  eq:d194 (Prop. prop:s6-ms): M^nu(D, H) = sum_{0<N(u)<=H} |A_u^nu(D)|^2, H = D^{1+theta}, theta = 1/10, for
       D = 100 ... 6400, nu of order 12 modulo pi13 = 4 + 3 omega (S = {2, lambda, pi13}); also nu = 1.
       Reported: M/D^{2+theta} and M/diag, diag = #{u} * sum_n |W(N(n)/D)|^2 (the diagonal term).
       Sanity criterion: M/D^{2+theta} <= 10 and M/diag in [0.5, 2] for every D.
  SZ2  eq:s6-bound: A_1^nu(D) = sum_n mu_K(n) alpha(n) nu(n) W(N(n)/D) for 120 values of D in [10^2, 10^6]
       (mu_K by a vectorized sieve on the norms).  Reported per decade: max |A_1|/D^{1/2}, max |A_1|/D^{11/12}, and
       the trivial size sum |W| (number of terms).  Sanity criterion: |A_1(D)| <= D^{11/12} for D >= 10^3.
       Also the finite inequality of Lemma lem:b1-extract (iii) for the twisted rows (u_0 = 1, N_0 = 1), as used in
       the proof of eq:s6-bound: |A_1(D) - A_{pi^6}(D)| <= (1 + 16 b) ||W||_inf D / Y for pi in P_Y, Y = H^{1/6}
       (this one is an exact inequality: PASS/FAIL).
  SZ3  eq:s4-marked (= eq:cm-inverse): sum_{0<N(u)<=Z^m} |M_u(Z^r; W) Q_u|^2 / Z^m with Z^m = H, for eps = +1 and -1,
       nu = 1 and nu = chi_.(lambda), no slots (r = log_H D with r + 2z <= 1 - c_1) and one slot of primes of norm
       in [13, 26) with unimodular coefficients (r + 2z <= 1 - c_1, 2r + 8z <= 3 - c_2 checked).
       Sanity criterion: ratio <= 10.
Usage: python3 t8_size_checks.py      (about 30 s)
"""
import os, sys, math, cmath, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import *
import numpy as np
random.seed(5)

PI13 = (4, 3); W0 = (-PI13[0]*pow(PI13[1], -1, 13)) % 13
DLOG = {}; x = 1
for k in range(12): DLOG[x] = k; x = (2*x) % 13
DLOG_ARR = np.zeros(13, int)
for r_, k in DLOG.items(): DLOG_ARR[r_] = k
def red13(z): return (z[0] + z[1]*W0) % 13
def nu12(n):
    r = red13(n); return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)
def muk(n): return 1 if n == (1, 0) else mu(n)
def Wb(y):                                             # bump on (1, 2), sup norm 1, complex phase
    y = np.asarray(y, float)
    out = np.zeros(y.shape, complex)
    m = (y > 1) & (y < 2)
    out[m] = np.exp(4 - 1/((y[m] - 1)*(2 - y[m])))*np.exp(0.7j*y[m])
    return out
B_SUPP = 2.0
THETA = 0.1

def lattice_rows(H):
    R = int(2*math.sqrt(H)) + 3
    UA, UB = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing='ij'); UA = UA.ravel(); UB = UB.ravel()
    NU = UA*UA - UA*UB + UB*UB; keep = (NU > 0) & (NU <= H)
    return UA[keep], UB[keep]

ok_all = True
# ------------------------------------------------------------------ SZ1 eq:d194
print("SZ1 eq:d194: twisted mean square, H = D^{1.1}")
worst1 = 0.0; band1 = (1e9, -1e9)
for nuname, nuf in (("order 12 mod pi13", nu12), ("1", lambda n: 1.0 + 0j)):
    for D in (100, 200, 400, 800, 1600, 3200, 6400):
        H = D**(1 + THETA)
        UA, UB = lattice_rows(H)
        cols = [n for n in primaries(2*D, sqfree=True) if D < N(n) < 2*D and nuf(n) != 0]
        Au = np.zeros(len(UA), complex); diag_w = 0.0
        for n in cols:
            c = muk(n)*alpha_c(n)*nuf(n)*Wb(N(n)/D)[()]
            Au += c*chi_vec(n, UA, UB); diag_w += abs(c)**2
        M = float(np.sum(np.abs(Au)**2)); diag = len(UA)*diag_w
        r1 = M/D**(2 + THETA); r2 = M/diag
        worst1 = max(worst1, r1); band1 = (min(band1[0], r2), max(band1[1], r2))
        print(f"      nu={nuname:18s} D={D:5d} H={H:9.1f} rows={len(UA):6d} cols={len(cols):5d}  M/D^(2+theta)={r1:.4f}  M/diag={r2:.4f}")
ok = worst1 <= 10 and 0.5 <= band1[0] and band1[1] <= 2
ok_all &= ok
print(f"[{'SANITY-OK' if ok else 'SANITY-FAIL'}] SZ1 eq:d194: max M/D^(2+theta) = {worst1:.3f}; M/diag in [{band1[0]:.3f}, {band1[1]:.3f}]")

# ------------------------------------------------------------------ SZ2 eq:s6-bound
print("SZ2 eq:s6-bound: A_1^nu(D) up to D = 10^6")
Dmax = 10**6; M2 = int(B_SUPP*Dmax) + 1
R = int(math.sqrt(4*M2/3)) + 3
ia = np.arange((-R - 1)//3, (R - 1)//3 + 2); jb = np.arange(-R//3 - 1, R//3 + 2)
AA, BB = np.meshgrid(3*ia + 1, 3*jb, indexing='ij'); AA = AA.ravel().astype(np.int64); BB = BB.ravel().astype(np.int64)
NN = AA*AA - AA*BB + BB*BB
keep = (NN <= M2) & ~((AA % 2 == 0) & (BB % 2 == 0))         # primary, coprime to 2 (hence to 6), N <= 2 Dmax
AA, BB, NN = AA[keep], BB[keep], NN[keep]
# smallest-prime-factor sieve on the norms
spf = np.zeros(M2 + 1, np.int64)
for p in range(2, int(M2**0.5) + 1):
    if spf[p] == 0:
        blk = spf[p*p::p]; blk[blk == 0] = p
spf[spf == 0] = np.arange(M2 + 1)[spf == 0]
MU = np.ones(len(NN), np.int64); rem = NN.copy()
while True:
    act = rem > 1
    if not act.any(): break
    p = np.where(act, spf[np.maximum(rem, 1)], 1)
    e = np.zeros(len(rem), np.int64)
    while True:
        m = act & (rem % p == 0)
        if not m.any(): break
        rem = np.where(m, rem//p, rem); e += m
    split = act & (p % 3 == 1); inert = act & (p % 3 == 2)
    MU = np.where(split & (e == 1), -MU, MU)
    both = (AA % np.maximum(p, 1) == 0) & (BB % np.maximum(p, 1) == 0)
    MU = np.where(split & (e == 2) & ~both, 0, MU)            # pi^2 | n
    MU = np.where(split & (e >= 3), 0, MU)
    MU = np.where(inert & (e == 2), -MU, MU)
    MU = np.where(inert & (e >= 4), 0, MU)
# cross-check the sieve against mu() of cmlib on small elements
small = np.nonzero(NN <= 3000)[0]
bad = sum(1 for i in small if MU[i] != mu((int(AA[i]), int(BB[i]))) and NN[i] > 1)
r13 = (AA + BB*W0) % 13
NUv = np.where(r13 == 0, 0, np.exp(2j*np.pi*DLOG_ARR[r13]/12))
ALv = (AA + BB*(-0.5 + 0.5j*math.sqrt(3)))/np.sqrt(NN)
COEF = MU*ALv*NUv
sel = COEF != 0
NNs, COEFs = NN[sel].astype(float), COEF[sel]
print(f"      {len(NN)} primary elements with N <= {M2}; sieve mu agrees with cmlib mu on {len(small)} elements: {bad} mismatches")
ok_all &= (bad == 0)
Ds = np.exp(np.linspace(math.log(100), math.log(Dmax), 120))
A1 = np.array([np.sum(COEFs*Wb(NNs/D)) for D in Ds])
triv = np.array([np.sum(np.abs(Wb(NNs/D))) for D in Ds])
viol = 0
for lo in (2, 3, 4, 5):
    m = (Ds >= 10**lo) & (Ds < 10**(lo + 1)) if lo < 5 else (Ds >= 10**5)
    a = np.abs(A1[m]); d = Ds[m]
    print(f"      D in [1e{lo}, 1e{lo+1}]: max|A_1|/D^(1/2) = {np.max(a/d**0.5):.3f}, max|A_1|/D^(11/12) = "
          f"{np.max(a/d**(11/12)):.4f}, trivial size sum|W| / D ~ {np.max(triv[m]/d):.3f}")
    if lo >= 3: viol += int(np.sum(a > d**(11/12)))
ok = viol == 0; ok_all &= ok
print(f"[{'SANITY-OK' if ok else 'SANITY-FAIL'}] SZ2 eq:s6-bound: |A_1^nu(D)| <= D^(11/12) at all {int(np.sum(Ds >= 1e3))} "
      f"sample points D >= 10^3")
# Lemma b1-extract (iii) for the twisted rows: exact inequality
e3 = 0.0; cnt = 0
for D in (1e5, 3e5, 1e6):
    H = D**(1 + THETA); Y = H**(1/6)
    PY = [p for p in PR_ODD if Y/2 < N(p) <= Y and red13(p) != 0]
    A1D = np.sum(COEFs*Wb(NNs/D))
    for p in PY:
        pc = conj(p); np_ = N(p)
        c1 = AA*pc[0] - BB*pc[1]; c2 = AA*pc[1] + BB*pc[0] - BB*pc[1]      # n * conj(p)
        divp = (c1 % np_ == 0) & (c2 % np_ == 0)
        lhs = abs(np.sum((COEF*Wb(NN/D))[divp]))                         # A_1 - A_{p^6} = sum over p | n
        rhs = (1 + 16*B_SUPP)*1.0*D/Y
        cnt += 1; e3 = max(e3, lhs/rhs)
        print(f"      (iii) D={D:.0e} Y={Y:.2f} pi={p} N={np_}: |A_1 - A_pi^6| = {lhs:.2f} <= {rhs:.1f}")
report("SZ2-iii", e3, 1.0, f"Lemma lem:b1-extract (iii) for the twisted rows, {cnt} (D, pi) pairs (err = lhs/rhs, must be <= 1)")

# ------------------------------------------------------------------ SZ3 eq:s4-marked
print("SZ3 eq:s4-marked: sum |M_u Q_u|^2 / Z^m, Z^m = H")
worst3 = 0.0
for H, D, slot in ((20000, 100, False), (20000, 400, False), (60000, 1000, False), (60000, 40, True)):
    UA, UB = lattice_rows(H)
    r = math.log(D)/math.log(H)
    if slot:
        Pl = [p for p in PR_ODD if 13 <= N(p) < 26]; z = math.log(13)/math.log(H)
        assert r + 2*z <= 1 - 0.05 and 2*r + 8*z <= 3 - 0.05, (r, z)
    else:
        Pl = []; z = 0.0
    for eps in (1, -1):
        for nuname, nuf in (("1", lambda n: 1.0 + 0j), ("chi(lambda)", lambda n: complex(chi(n, (1, 2))))):
            cols = [n for n in primaries(2*D, sqfree=True) if D < N(n) < 2*D]
            Mu = np.zeros(len(UA), complex)
            for n in cols:
                c = muk(n)*(alpha_c(n) if eps == 1 else alpha_c(n).conjugate())*nuf(n)*Wb(N(n)/D)[()]
                Mu += c*chi_vec(n, UA, UB, eps % 6)
            Mu /= math.sqrt(D)
            Q = np.ones(len(UA), complex)
            if slot:
                Q = np.zeros(len(UA), complex)
                for p in Pl:
                    a_p = cmath.exp(2j*math.pi*random.random())
                    Q += a_p*nuf(p)*chi_vec(p, UA, UB, eps % 6)*Wb(N(p)/13.0)[()]
                Q /= math.sqrt(13.0)
            ratio = float(np.sum(np.abs(Mu*Q)**2))/H
            worst3 = max(worst3, ratio)
            print(f"      H={H} D={D} r={r:.3f} slot={'yes' if slot else 'no '} z={z:.3f} eps={eps:+d} nu={nuname:11s}: "
                  f"sum|M_u Q_u|^2/H = {ratio:.4f}")
ok = worst3 <= 10; ok_all &= ok
print(f"[{'SANITY-OK' if ok else 'SANITY-FAIL'}] SZ3 eq:s4-marked: max ratio {worst3:.3f}")
report("SZ-sanity", 0.0 if ok_all else 1.0, 0.0, "all sanity criteria of SZ1-SZ3 (and the mu sieve cross-check) hold")
