"""t6 (2026-10-13): Proposition prop:s1-reflection(a), eq:d127 (the Mellin representation of the completed sum), with the
definitions eq:d124 (D_2, L_2, T_2 = D_2 L_2) and eq:d125 (T_V), and the Mellin form of eq:s1-central ((R6) with
Psi = nu chi_.(m)).

Setting: base twist A (Psi_0 = conj((lambda/.)_3), E = (6)), one local prime pi of norm 7 with j = 0..5, as in
tests_cm/phase2_2026-10-12/s1/t_s1_endtoend_active.py.  Two log-Schwartz tests without compact support:
  'G'   V_G(y) = (2 sqrt(pi))^{-1} exp(-(log y)^2/4),  Vhat_G(t) = e^{t^2}   (the guide's Gaussian; Vhat checked too)
  'K'   V(y) = exp(-y-1/y),                          Vhat(t) = 2 K_t(2)
Left side: T_V(X;Psi) from eq:d125, with Psi(n b^3) and conj(alpha(n b^3))^2 evaluated on the whole index n b^3.
Right side: (1/2 pi i) int_{(a)} Vhat(s-1/2) X^{s-1/2} D_2(s,Psi) L_2(s,Psi) ds with D_2 and L_2 computed as two separate
Dirichlet series (trapezoid rule on Re s = a = 1.5 and 2).  Truncating the Dirichlet series only drops terms whose
weight V(N(n)N(b)^3/X) is below 1e-13.
Falsifiers: Vhat(s) X^s (no shift by 1/2); L_2 with N(b)^{-3s} (no +1/2); OAI78's angular powers conj(alpha),
conj(alpha)^3.
eq:s1-central: C^(2)_V(X;m,nu) from its definition (chi_{cn^3}(m) from the prime factorisation of c n^3, with
multiplicity) against the same Mellin integral with Psi = nu chi_.(m), i.e. chi_c(m) in D_2 and chi_n(m)^3 in L_2,
for nu = Psi_0 and three rows m (with a square factor, with a factor at lambda and at 2, and a unit).
eq:d144/eq:d145: the inner integral of the compensated completed row, (1/2 pi i) int_{(4)} X^t e^{t^2} T^(2)_D(t+1/2,nu') dt,
equals the marked completed sum with the Gaussian V_G (D = pi, nu' = Psi_0 chi_.(m)).
Usage: python3 t6_mellin_rep.py      (about 1 minute)
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
from t6_lib import *
from zw import sym6_prime, prim_assoc
import mpmath as mp
t0 = time.time()
S = Setup(25000, 0)
CHECKS = []
def check(name, ok, msg):
    CHECKS.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)
def cplx(x): return complex(cx(x))

mp.mp.dps = 25
# Vhat_G(t) = e^{t^2}: check against the definition int V_G(y) y^{t-1} dy at three points
VG = lambda y: mp.exp(-mp.log(y)**2/4)/(2*mp.sqrt(mp.pi))
err = max(abs(mp.quad(lambda u: VG(mp.e**u)*mp.e**(t*u), [-40, -10, 0, 10, 40])-mp.e**(t*t))/abs(mp.e**(t*t))
          for t in [mp.mpf('0.3'), mp.mpc(1, 2), mp.mpc(-0.7, 0.4)])
errK = max(abs(mp.quad(lambda u: mp.exp(-mp.e**u-mp.e**-u)*mp.e**(t*u), [-7, -2, 0, 2, 7])-2*mp.besselk(t, 2))
           for t in [mp.mpf('0.3'), mp.mpc(1, 2)])
check("Vhat of the two tests", err < 1e-12 and errK < 1e-12, f"Vhat_G(t) = e^(t^2): {mp.nstr(err, 2)}; Vhat(t) = 2K_t(2) for exp(-y-1/y): {mp.nstr(errK, 2)}")

mp.mp.dps = 15
TESTS = {'G': (lambda y: np.exp(-np.log(y)**2/4)/(2*math.sqrt(math.pi)), lambda t: np.exp(t*t), 7.0),
         'K': (lambda y: np.exp(-y-1/y), None, 26.0)}
# ---------------- data: squarefree odd n with factorisations; odd primary b ----------------
NL = []
for cc_, nc, g1_, _, _, _, slam, even, haspi in S.cdata:
    if even: continue
    NL.append((cc_, nc, slam*g1_/abs(cplx(cc_))))                     # gamma_2(n) (fast; brute-force check in t6_theta_translates)
BL = []
for x in range(-60, 61):
    for y in range(-60, 61):
        b = (x, y)
        if primary(b) and not (x % 2 == 0 and y % 2 == 0) and N(b)**3 <= S.B: BL.append(b)
fac_b = {}
for b in BL:
    f = []; r = b
    for p in S.plist:
        if N(r) == 1: break
        while divides(p, r): r = div(r, p); f.append(p)
    fac_b[b] = f
nN = np.array([x[1] for x in NL], float); nC = np.array([cplx(x[0]) for x in NL]); nG = np.array([x[2] for x in NL])
nR = np.array([S.red(x[0]) for x in NL])
nPsi0 = np.array([np.prod([S.psi0(p) for p in S.sqf[x[0]]]) if S.sqf[x[0]] else 1 for x in NL])
bN = np.array([N(b) for b in BL], float); bC = np.array([cplx(b) for b in BL]); bR = np.array([S.red(b) for b in BL])
print(f"{len(NL)} squarefree n, {len(BL)} b; {time.time()-t0:.1f} s")

PAIRS = {}
def pairs(X, ycut):
    """All (n, b) with N(n)N(b)^3/X < ycut: weights and the whole-index data of A = n b^3 (exact)."""
    key = (X, ycut)
    if key not in PAIRS:
        rows = []
        for bi, b in enumerate(BL):
            nb = bN[bi]
            for i in np.nonzero(nN*nb**3/X < ycut)[0]:
                n_ = NL[i][0]; A = mul(n_, pw(b, 3)); Ac = cplx(A)
                rows.append((i, bi, np.conj(S.lam3(A)), S.red(A), np.conj(Ac/abs(Ac))**2, S.sqf[n_]+3*fac_b[b]))
        PAIRS[key] = rows
    return PAIRS[key]
def T_def(X, kind, j, chim=None):
    """eq:d125 with Psi = Psi_0 chi_pi^j (or Psi = Psi_0 chi_.(m) when chim is given), evaluated on the whole index."""
    V = TESTS[kind][0]; ycut = 1e5 if kind == 'G' else 40.0
    rows = pairs(X, ycut)
    i = np.array([r[0] for r in rows]); bi = np.array([r[1] for r in rows])
    lam = np.array([r[2] for r in rows]); rA = np.array([r[3] for r in rows]); aA = np.array([r[4] for r in rows])
    if chim is None: psi = lam*S.chipow(rA, j)                        # Psi_0(nb^3) chi_pi(nb^3)^j, zero retained
    else: psi = lam*np.array([chim(r[5]) for r in rows])               # chi_{nb^3}(m) over the primes of nb^3
    y = nN[i]*bN[bi]**3/X
    return np.sum(nG[i]*aA*psi/(np.sqrt(nN[i])*bN[bi])*V(y))

def mellin(X, kind, j, a, shift=0.5, lexp=0.5, apow=(2, 6), dn=None, db=None):
    V, Vh, U = TESTS[kind]
    h = 0.02; u = np.arange(-U, U+1e-9, h); s = a+1j*u
    if Vh is None:
        key = (round(a-shift, 6), U)
        if key not in VHK: VHK[key] = np.array([complex(2*mp.besselk(mp.mpc(a-shift, uu), 2)) for uu in u])
        vh = VHK[key]
    else:
        vh = Vh(s-shift)
    key = (U, a, j, apow, lexp)
    if dn is not None or key not in DL:
        cn = (dn if dn is not None else nG*np.conj(nC/np.abs(nC))**apow[0]*nPsi0*S.chipow(nR, j))
        cb = (db if db is not None else np.conj(bC/np.abs(bC))**apow[1]*S.chipow(bR, (3*j) % 6))
        D = np.exp(-np.outer(s, np.log(nN))) @ cn
        L = np.exp(-np.outer(3*s-lexp, np.log(bN))) @ cb
        if dn is not None: return np.sum(vh*X**(s-shift)*D*L)*h/(2*np.pi)
        DL[key] = (D, L)
    D, L = DL[key]
    return np.sum(vh*X**(s-shift)*D*L)*h/(2*np.pi)
VHK = {}; DL = {}
worst = 0; fal = {'noshift': 0, 'L2exp': 0, 'zeta-angles': 0}; rows = []
for kind, X in [('G', 0.5), ('G', 3.0), ('K', 40.0), ('K', 300.0)]:
    for j in range(6):
        lhs = T_def(X, kind, j)
        for a in (1.5, 2.0):
            rhs = mellin(X, kind, j, a)
            worst = max(worst, abs(rhs-lhs)/abs(lhs))
        fal['noshift'] = max(fal['noshift'], abs(mellin(X, kind, j, 2.0, shift=0.0)-lhs)/abs(lhs))
        fal['L2exp'] = max(fal['L2exp'], abs(mellin(X, kind, j, 2.0, lexp=0.0)-lhs)/abs(lhs))
        fal['zeta-angles'] = max(fal['zeta-angles'], abs(mellin(X, kind, j, 2.0, apow=(1, 3))-lhs)/abs(lhs))
    rows.append(f"{kind} X={X}: T_V(j=5) = {lhs:.10f}")
    print("   ", rows[-1], f"({time.time()-t0:.1f} s)", flush=True)
check("eq:d127 (with eq:d124, eq:d125)", worst < 1e-9 and min(fal.values()) > 1e-2,
      f"2 tests x 2 X x j=0..5 x a=1.5,2: max rel err {worst:.1e}; falsifiers (largest deviation): "
      + ", ".join(f"{k} {v:.2f}" for k, v in fal.items()))

# ---------------- eq:s1-central: Psi = nu chi_.(m), nu = Psi_0 ----------------
MS = {'pi1^2 pi2': mul(pw((-2, -3), 2), (1, 3)), 'lambda 2 pi3': mul(mul(LAMp, (2, 0)), (4, 3)), 'unit omega': (0, 1)}
worst = 0; fal2 = 0
for name, m in MS.items():
    def chim(primes, m=m):
        r = 1+0j
        for p in primes:
            k = sym6_prime(m, p)
            if k is None: return 0
            r *= cmath.exp(1j*math.pi*k/3)
        return r
    def chim1(n_): return chim(S.sqf[n_])
    cm_n = np.array([chim1(x[0]) for x in NL])
    cm_b = np.array([chim(fac_b[b]) for b in BL])
    dn = nG*np.conj(nC/np.abs(nC))**2*nPsi0*cm_n
    db = np.conj(bC/np.abs(bC))**6*cm_b**3                       # chi_n(m)^3 (Psi_0^3 = 1)
    for kind, X in [('G', 30.0), ('K', 300.0)]:
        lhs = T_def(X, kind, 0, chim=chim)
        rhs = mellin(X, kind, 0, 2.0, dn=dn, db=db)
        worst = max(worst, abs(rhs-lhs)/abs(lhs))
        fal2 = max(fal2, abs(mellin(X, kind, 0, 2.0, dn=dn, db=np.conj(bC/np.abs(bC))**6*cm_b)-lhs)/abs(lhs))
check("eq:s1-central Mellin form (R6), Psi = nu chi_.(m)", worst < 1e-9 and fal2 > 1e-3,
      f"3 rows m (square factor, lambda and 2, unit) x 2 tests: max rel err {worst:.1e}; falsifier chi_n(m) instead of chi_n(m)^3: {fal2:.1e}")
# ---------------- eq:d144 / eq:d145: the inner Mellin integral of the compensated row with a mark D ----------------
# (1/2 pi i) int_{(4)} X^t Phi(t) T^(2)_D(t+1/2, nu') dt, Phi(t) = e^{t^2}, equals the marked completed sum with the
# Gaussian test V_G (Vhat_G = Phi):  sum gamma_2(c) conj(alpha(cn^3))^2 nu'(cn^3) 1_{D | cn^3} V_G(N(c)N(n)^3/X)/(sqrt(N c) N n).
# D = pi (norm 7), nu' = Psi_0 chi_.(m) with m = (1,3) (a prime of norm 7 other than pi), X = 0.7, 2.
m2 = (1, 3) if S.PI != (1, 3) else (-2, -3)
def chim2(primes):
    r = 1+0j
    for p in primes:
        k = sym6_prime(m2, p)
        if k is None: return 0
        r *= cmath.exp(1j*math.pi*k/3)
    return r
worst = 0; falm = 0
for X in (0.7, 2.0):
    rows = pairs(X, 1e5)
    i = np.array([r[0] for r in rows]); bi = np.array([r[1] for r in rows])
    co = nG[i]*np.array([r[4] for r in rows])*np.array([r[2] for r in rows])*np.array([chim2(r[5]) for r in rows])
    mark = np.array([1.0 if (r[3] == 0) else 0.0 for r in rows])          # pi | n b^3  iff  red(n b^3) = 0
    Nc_, Nn_ = nN[i], bN[bi]
    lhs = np.sum(co*mark*np.exp(-np.log(Nc_*Nn_**3/X)**2/4)/(2*math.sqrt(math.pi))/(np.sqrt(Nc_)*Nn_))
    h = 0.01; u = np.arange(-7, 7+1e-9, h); tt = 4+1j*u
    TD = np.exp(-np.outer(tt+0.5, np.log(Nc_))-np.outer(3*(tt+0.5)-0.5, np.log(Nn_))) @ (co*mark)
    rhs = np.sum(X**tt*np.exp(tt*tt)*TD)*h/(2*np.pi)
    worst = max(worst, abs(rhs-lhs)/abs(lhs))
    TDw = np.exp(-np.outer(tt, np.log(Nc_))-np.outer(3*tt-0.5, np.log(Nn_))) @ (co*mark)              # falsifier: T_D(t), no shift
    falm = max(falm, abs(np.sum(X**tt*np.exp(tt*tt)*TDw)*h/(2*np.pi)-lhs)/abs(lhs))
check("eq:d144/eq:d145 inner Mellin integral (mark D = pi)", worst < 1e-8 and falm > 1e-2,
      f"X = 0.7, 2 on the line Re t = 4: max rel err {worst:.1e}; falsifier T_D(t) instead of T_D(t+1/2): {falm:.2f}")
nf = sum(1 for _, ok in CHECKS if not ok)
print(f"{len(CHECKS)} checks, {nf} failed; time {time.time()-t0:.1f} s")
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1': nf += 1; print("CM_TESTS_FORCE_FAIL=1")
print("PASS" if nf == 0 else "FAIL")
raise SystemExit(0 if nf == 0 else 1)
