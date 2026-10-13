"""Block b4 checks (CM guide, Sections 9-10): new identities introduced in drafts_cm_fill/b4.tex.

C1  CM paired identity, forward direction with a twist xi1 (literal sum over the 432 characters of frak G).
C2  CM paired identity, reverse direction (literal sum over the 432 characters).
C3  Second Poisson identity (CM), end to end with the self-dual Gaussian Phi(t) = exp(-2 pi t/sqrt 3).
C4  First Poisson identity of the transfer (CM), end to end with the same Gaussian.
C5  Regrouping: sum of mu(e)mu(w) over the preimages of (r, f', k') equals tau(r) 1_{f' | k'}.
C6  Cube removal: the long part of the inversion equals the b-sum with beta_0(b, f) (exact, conj(alpha)^6).
C7  Hecke-family facts behind the use of [GL]: (x/k)_2 kappa_lambda(x)^{e_k} is trivial on units; reciprocity
    factor depends only on classes mod 4; even columns.
C8  Local-cost bookkeeping in the squarefree-row lemma (exponent enumeration).
Falsifiers (wrong variants that must fail) are printed after each check.
"""
import os, sys, math, cmath, itertools, random
import os; LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'part_B_lines_0550_0997')
sys.path.insert(0, LIB)
from cmlib import *          # numpy as np, zw helpers, chi, mu, factor, primaries, gamma_j_fast, G_fast, chi_vec, ...
random.seed(4)
SQ3 = math.sqrt(3.0)

# Exit status (added 2026-10-12): every [FAIL] line and every falsifier that DID NOT FAIL is counted in NFAILED, and
# the script ends with sys.exit(1) if NFAILED > 0.  CM_TESTS_FORCE_FAIL=1 makes the first check count as failed.
NFAILED = [0]
FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']
def _forced():
    if FORCE[0]:
        FORCE[0] = False
        return True
    return False
def line(name, err, tol, extra=""):
    ok = err <= tol
    if _forced(): ok = False; extra = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + extra
    NFAILED[0] += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: err {err:.3g} (tol {tol:g}) {extra}", flush=True)
def fals(name, err, thr=1e-3):
    ok = err > thr
    if _forced(): ok = False; name = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + name
    NFAILED[0] += not ok
    print(f"      falsifier [{name}]: err {err:.3g} -> {'fails as required' if ok else 'DID NOT FAIL'}", flush=True)

# ------------------------------------------------------------------ the group frak G (primary classes mod 72)
M72 = 72
els = [(a, b) for a in range(M72) for b in range(M72) if a % 3 == 1 and b % 3 == 0 and not (a % 2 == 0 and b % 2 == 0)]
idx = {e: i for i, e in enumerate(els)}
def red(x): return (x[0] % M72, x[1] % M72)
order = len(els)
mt = np.array([[idx[red(mul(x, y))] for y in els] for x in els])
one = idx[(1, 0)]
inv = np.array([int(np.where(mt[i] == one)[0][0]) for i in range(order)])
def centered(v): return v - M72 if v > M72//2 else v
Gtab = np.zeros(order, complex)
for e in els:
    a0, b0 = centered(e[0]), centered(e[1]); best = None
    for i in range(-3, 4):
        for j in range(-3, 4):
            n = (a0 + M72*i, b0 + M72*j)
            if n == (1, 0): best = n; break
            if N(n) >= 2 and mu(n) != 0 and (best is None or N(n) < N(best)): best = n
        if best == (1, 0): break
    Gtab[idx[e]] = 1.0 if best == (1, 0) else G_fast(best)
def cls(n): return idx[red(n)]
def Gquot(a, b):
    """G(a b^{-1}) on frak G."""
    return Gtab[mt[cls(a), inv[cls(b)]]]

# characters of frak G from the regular representation (as in tests_cm/.../t3_gauss.py)
rng = np.random.default_rng(5)
gens = [idx[e] for e in random.sample(els, 6)]
Mreg = np.zeros((order, order), complex)
for g in gens:
    c = rng.normal() + 1j*rng.normal()
    for x in range(order): Mreg[mt[g, x], x] += c
w_, V_ = np.linalg.eig(Mreg)
chars = np.array([np.conj(V_[:, k]/V_[one, k]) for k in range(order)])
chars = np.exp(2j*np.pi*np.round(np.angle(chars)/(2*np.pi)*order)/order)
assert max(np.max(np.abs(chars[:, mt[i, :]] - chars[:, [i]]*chars)) for i in range(order)) < 1e-8

g2c = {}
def g2(n):
    if n == (1, 0): return 1.0 + 0j
    if n not in g2c: g2c[n] = gamma_j_fast(2, n)
    return g2c[n]
def al(n): return alpha_c(n) if n != (1, 0) else 1.0 + 0j
def muu(n): return 1 if n == (1, 0) else mu(n)

sq = primaries(400, sqfree=True)            # contains n = 1

# ------------------------------------------------------------------ C1, C2
k1 = 7; k2 = 11                         # xi1 = chars[k1] (forward twist), xi = chars[k2] (reverse)
Ginv = Gtab[inv]
c_fw = (chars.conj() @ (chars[k1]*Ginv))/order          # xi1(t) G(t^{-1}) = sum c_{xi'} xi'(t)
l_rv = (chars.conj() @ (chars[k2]*Gtab))/order          # (xi G)(t) = sum l_{xi1} xi1(t)
pairs = []
cand = [n for n in sq if N(n) <= 250]
while len(pairs) < 60:
    z1, z2 = random.choice(cand), random.choice(cand)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 4000: continue
    pairs.append((z1, z2))
e1 = e2 = f1 = f1b = f2 = 0
for z1, z2 in pairs:
    X1 = chars[:, cls(z1)]; X2 = chars[:, cls(z2)]
    a1 = al(z1).conjugate()**2*g2(z1); a2 = al(z2).conjugate()**2*g2(z2)
    gam = gauss_psi_general(z1, z2, conj_first=True)
    lhs = muu(z1)*muu(z2)*al(z1).conjugate()*al(z2)*chars[k1, cls(z1)]*np.conj(chars[k1, cls(z2)])*gam
    rhs = np.sum(c_fw*a1*X1*np.conj(a2*X2))
    e1 = max(e1, abs(lhs - rhs))
    za1 = al(z1).conjugate()*g2(z1); za2 = al(z2).conjugate()*g2(z2)
    f1 = max(f1, abs(lhs - np.sum(c_fw*za1*X1*np.conj(za2*X2))))
    lhs_b = muu(z1)*muu(z2)*al(z1).conjugate()*al(z2)*np.conj(chars[k1, cls(z1)])*chars[k1, cls(z2)]*gam
    f1b = max(f1b, abs(lhs_b - rhs))
    gr = gauss_psi_general(z1, z2, conj_first=False)
    L2 = a1*chars[k2, cls(z1)]*np.conj(a2*chars[k2, cls(z2)])*gr
    R2 = muu(z1)*muu(z2)*np.sum(l_rv*al(z1).conjugate()*X1*np.conj(al(z2).conjugate()*X2))
    e2 = max(e2, abs(L2 - R2))
    f2 = max(f2, abs(L2 - muu(z1)*muu(z2)*np.sum(l_rv*X1*np.conj(X2))))
line("C1 forward paired identity with twist xi1 (CM)", e1, 1e-10, f"{len(pairs)} pairs, 432 characters")
fals("zeta coefficients conj(alpha) gamma_2", f1)
fals("xi1 conjugated on the wrong column", f1b)
line("C2 reverse paired identity (CM): Moebius coefficients mu conj(alpha) xi1", e2, 1e-10)
fals("Moebius coefficients without conj(alpha) (zeta form)", f2)

# ------------------------------------------------------------------ smooth weights
def bump(x, a, b):
    x = np.asarray(x, float); out = np.zeros_like(x)
    m = (x > a) & (x < b)
    out[m] = np.exp(-1.0/((x[m] - a)*(b - x[m])))
    return out
def Uw(x, a, b): return bump(x, a, b)*np.exp(1j*np.asarray(x, float))     # complex test weight
def Phi(t): return np.exp(-2*np.pi*np.asarray(t, float)/SQ3)              # self-dual: Phihat = Phi
Phihat = Phi
CUT = 22.0                                                                # Phi(CUT) ~ e^{-80}

def latt(R):
    pts = lattice(R); A = np.array([p[0] for p in pts], np.int64); B = np.array([p[1] for p in pts], np.int64)
    return A, B
def mulvec(c, A, B):
    """c * (A + B omega) for arrays A, B."""
    a, b = c
    return a*A - b*B, a*B + b*A - b*B
def Nvec(A, B): return A*A - A*B + B*B
def coprime(a, b): return N(gcd_primary(a, b)) == 1
def sqfree_list(lo, hi, extra=lambda n: True):
    # primaries(...) already contains n = 1
    return [n for n in primaries(int(hi), sqfree=True) if N(n) >= lo and extra(n)]

# ------------------------------------------------------------------ C3: second Poisson identity
def check_C3(C, t, d, ell, Y, xi1, cm=True, a=0.5, b=2.0, coef_rhs=None):
    """LHS = sum_y Phi(N y/Y) |P(y)|^2 ;  RHS = Z' + sum_{g,e,w,h,xi'} ... (with Sigma_{xi'} c_{xi'} written out
    as xi1(n1) conj(xi1(n2)) G(n2 n1^{-1}))."""
    cols = sqfree_list(a*ell, b*ell, lambda n: coprime(n, t))
    A, B = latt(CUT*Y)
    P = np.zeros(A.shape, complex)
    for n in cols:
        cf = muu(n)*(al(n).conjugate() if cm else 1.0)*xi1(n)*chi(n, C)**4*chi(n, d)*Uw(N(n)/ell, a, b)
        if cf == 0: continue
        P += cf*np.conj(chi_vec(n, A, B))
    lhs = float(np.sum(Phi(Nvec(A, B)/Y)*np.abs(P)**2))
    # zero frequency
    Z = 0.0
    for g in sqfree_list(a*ell, b*ell, lambda n: coprime(n, mul(C, t))):
        pr = 1.0
        for p in factor(g): pr *= 1 - 1/N(p)
        Z += Y*Phihat(0.0)*abs(Uw(N(g)/ell, a, b))**2*pr
    U0 = lambda x: Uw(x, a, b)/np.sqrt(np.asarray(x, float))
    coef = coef_rhs or (lambda n: al(n).conjugate()**2*g2(n))
    main = 0j
    gs = sqfree_list(1, b*ell, lambda n: coprime(n, mul(C, t)))
    for g in gs:
        for e in divisors_sqfree(g):
            for w in sqfree_list(1, b*ell/N(g), lambda n: coprime(n, mul(mul(g, C), t))):
                Ngw = N(g)*N(w)
                r_ = mul(t, div(g, e))
                ncols = sqfree_list(a*ell/Ngw, b*ell/Ngw, lambda n: coprime(n, r_))
                if not ncols: continue
                x = np.array([N(n)*Ngw/ell for n in ncols])
                kappa = Y*N(g)**2/(N(e)*ell**2)
                Hmax = CUT*x.max()**2/kappa
                Ah, Bh = latt(Hmax)
                nz = (Ah != 0) | (Bh != 0)
                Ah, Bh = Ah[nz], Bh[nz]
                if N(w) > 1:   # (h, w) = 1
                    keep = np.ones(Ah.shape, bool)
                    for p in factor(w):
                        keep &= chi_vec(p, Ah, Bh, 6) != 0
                    Ah, Bh = Ah[keep], Bh[keep]
                Nh = Nvec(Ah, Bh)
                de = mul(d, e); Ade, Bde = mulvec(de, Ah, Bh)
                Cew = mul(mul(C, e), w)
                V = np.array([coef(n)*xi1(n)*chi(n, Cew)**4*chi_vec(n, Ade, Bde) for n in ncols])   # rows n, cols h
                u0 = U0(x)
                pre = Y*muu(e)*muu(w)*N(g)/(N(e)*ell)
                for i, n1 in enumerate(ncols):
                    for j, n2 in enumerate(ncols):
                        Fq = Gquot(n2, n1)
                        ker = Phihat(kappa*Nh/(x[i]*x[j]))
                        main += pre*Fq*u0[i]*np.conj(u0[j])*np.sum(V[i]*np.conj(V[j])*ker)
    rhs = Z + main
    return lhs, rhs

xi1a = lambda n: complex(chi(n, (-1, 0)))*complex(cubic3(n, (0, 1)))     # chi_n(-1) (omega/n)_3: a character of frak G
C_, t_ = (1, 3), (4, 3)          # 1+3w (norm 7), 4+3w (norm 13)
C2_ = mul((1, 3), (-5, 0))      # norm 175
cases = [(C_, t_, C_, 30.0, 400.0), (C_, (1, 0), (1, 0), 30.0, 250.0), (C2_, t_, (-5, 0), 40.0, 500.0)]
worst = 0
for (C, t, d, ell, Y) in cases:
    lhs, rhs = check_C3(C, t, d, ell, Y, xi1a)
    worst = max(worst, abs(lhs - rhs)/abs(lhs))
    print(f"      C3 case C={C} t={t} e0={d} ell={ell} Y={Y}: LHS {lhs:.12g}  RHS {rhs.real:.12g}{rhs.imag:+.2g}i")
line("C3 second Poisson identity (CM coefficients mu conj(alpha) xi1 -> a_{xi'})", worst, 1e-9, "relative error")
lhs, rhs = check_C3(C_, t_, C_, 30.0, 400.0, xi1a, cm=False)
fals("zeta Moebius coefficients mu xi1 on the left", abs(lhs - rhs)/abs(lhs), 1e-8)
lhs, rhs = check_C3(C_, t_, C_, 30.0, 400.0, xi1a, coef_rhs=lambda n: al(n).conjugate()*g2(n))
fals("zeta cubic coefficients conj(alpha) gamma_2 on the right", abs(lhs - rhs)/abs(lhs), 1e-8)

# ------------------------------------------------------------------ C4: first Poisson identity of the transfer
def check_C4(L, F, Hs, xi, rev_cm=True, u=1.0, v=2.0):
    fs = sqfree_list(F, 2*F - 1e-9)
    cols = sqfree_list(u*L, v*L)
    A, B = latt(CUT*Hs)
    lhs = 0.0
    for f in fs:
        S = np.zeros(A.shape, complex)
        for n in cols:
            cf = al(n).conjugate()**2*g2(n)*xi(n)*chi(n, f)**4*Uw(N(n)/L, u, v)
            if cf == 0: continue
            S += cf*chi_vec(n, A, B)
        lhs += float(np.sum(Phi(Nvec(A, B)/Hs)*np.abs(S)**2))
    lhs /= (L*F)
    Z = 0.0
    for f in fs:
        for Cc in cols:
            if not coprime(Cc, f): continue
            pr = 1.0
            for p in factor(Cc): pr *= 1 - 1/N(p)
            Z += Hs*Phihat(0.0)*abs(Uw(N(Cc)/L, u, v))**2*pr
    Z /= (L*F)
    main = 0j
    for Cc in sqfree_list(1, v*L):
        for t in sqfree_list(1, v*L/N(Cc), lambda n: coprime(n, Cc)):
            ell = L/(N(Cc)*N(t))
            xs = sqfree_list(u*ell, v*ell, lambda n: coprime(n, t))
            if not xs: continue
            z = np.array([N(x)/ell for x in xs])
            for d in divisors_sqfree(Cc):
                wCd = Hs*N(Cc)/(N(d)*L**2*F)
                kappa = Hs*N(Cc)**2/(L**2*N(d))
                Hmax = CUT*z.max()**2/kappa
                Ah, Bh = latt(Hmax)
                nz = (Ah != 0) | (Bh != 0); Ah, Bh = Ah[nz], Bh[nz]
                if N(t) > 1:
                    keep = np.ones(Ah.shape, bool)
                    for p in factor(t): keep &= chi_vec(p, Ah, Bh, 6) != 0
                    Ah, Bh = Ah[keep], Bh[keep]
                Nh = Nvec(Ah, Bh)
                for f in fs:
                    if not coprime(f, mul(Cc, t)): continue
                    f2 = mul(f, f); Ay, By = mulvec(f2, Ah, Bh)
                    Vm = np.array([muu(x)*(al(x).conjugate() if rev_cm else 1.0)*chi(x, Cc)**4*chi(x, d)
                                   *np.conj(chi_vec(x, Ay, By)) for x in xs])
                    Wz = Uw(z, u, v)/np.sqrt(z)
                    pre = muu(d)*muu(t)*wCd
                    for i, x1 in enumerate(xs):
                        for j, x2 in enumerate(xs):
                            E = xi(x1)*np.conj(xi(x2))*Gquot(x1, x2)          # (xi G)(x1 x2^{-1})
                            ker = Phihat(kappa*Nh/(z[i]*z[j]))
                            main += pre*E*Wz[i]*np.conj(Wz[j])*np.sum(Vm[i]*np.conj(Vm[j])*ker)
    return lhs, Z + main

xia = lambda n: complex(cubic3(n, (0, 1)))
worst = 0
for (L, F, Hs) in [(20.0, 7.0, 25.0), (25.0, 1.0, 30.0)]:
    lhs, rhs = check_C4(L, F, Hs, xia)
    worst = max(worst, abs(lhs - rhs)/abs(lhs))
    print(f"      C4 case L={L} F={F} H={Hs}: LHS {lhs:.12g}  RHS {rhs.real:.12g}{rhs.imag:+.2g}i")
line("C4 first Poisson identity of the transfer (CM): A(W) = Z + sum ... p_{hf^2}", worst, 1e-9, "relative error")
lhs, rhs = check_C4(20.0, 7.0, 25.0, xia, rev_cm=False)
fals("p_y without conj(alpha) (zeta form)", abs(lhs - rhs)/abs(lhs))

# ------------------------------------------------------------------ C5: regrouping over preimages
def sqf_divs(n): return divisors_sqfree(n)
pr_small = [p for p in PR_ODD if N(p) <= 61][:12]
bad = 0; cnt = 0
for trial in range(400):
    ps = random.sample(pr_small, 6)
    r = (1, 0); fp = (1, 0)
    for p in ps[:random.randint(0, 2)]: r = mul(r, p)
    rest = ps[3:]
    for p in rest[:random.randint(0, 3)]: fp = mul(fp, p)
    # k' : random product of primes from f' and others (with multiplicity), times a unit
    kp = random.choice(UNITS)
    for p in random.sample(pr_small, random.randint(0, 4)): kp = mul(kp, p)
    for p in factor(fp) if fp != (1, 0) else []:
        if random.random() < 0.6: kp = mul(kp, p)
    # brute force over all tuples (t, C, e, w, d) with the conditions of the second Poisson identity
    tot = 0
    fpr = list(factor(fp)) if fp != (1, 0) else []
    for t in sqf_divs(r):
        for assign in itertools.product(range(3), repeat=len(fpr)):
            C = e = w = (1, 0)
            for p, s in zip(fpr, assign):
                if s == 0: C = mul(C, p)
                elif s == 1: e = mul(e, p)
                else: w = mul(w, p)
            for d in sqf_divs(C):
                de = mul(d, e)
                if not divides(de, kp): continue
                h = div(kp, de)
                g = mul(e, div(r, t))
                # conditions: (g, C t) = 1, (w, g C t h) = 1, h != 0
                if not coprime(g, mul(C, t)): continue
                ok = all(not divides(p, h) for p in (factor(w) if w != (1, 0) else []))
                if not ok: continue
                tot += muu(e)*muu(w)
    tau_r = len(sqf_divs(r))
    want = tau_r*(1 if divides(fp, kp) else 0)
    bad += (tot != want); cnt += 1
line("C5 regrouping: sum over preimages of mu(e)mu(w) = tau(r) 1_{f'|k'}", bad, 0, f"{cnt} random (r, f', k')")

# ------------------------------------------------------------------ C6: cube removal, long part
Xc = 3000.0; Hc = 2.0
sqX = primaries(int(2*Xc), sqfree=True)
bs = primaries(int((2*Xc)**(1/3)) + 2)
def Wc(y): return bump(y, 1.0, 2.0)*np.exp(0.7j*np.asarray(y, float))
def Tcm(Yv, k, f, xi, bexp=6):
    s = 0j
    for b_ in bs:                      # bs contains b = 1
        Nb = N(b_)
        if Nb**3 > 2*Yv: continue
        Pb = xi(b_)*chi(b_, k)*chi(b_, f)**4 if b_ != (1, 0) else 1.0
        if Pb == 0: continue
        cb = (al(b_).conjugate()**bexp)*Pb**3*math.sqrt(Nb)
        for n in sqX:
            yv = N(n)*Nb**3/Yv
            if not (1 < yv < 2): continue
            s += al(n).conjugate()**2*g2(n)*xi(n)*chi(n, k)*chi(n, f)**4*cb*complex(Wc(yv))
    return s/math.sqrt(Yv)
def Psi(n, k, f, xi): return 1.0 if n == (1, 0) else xi(n)*chi(n, k)*chi(n, f)**4
err6 = 0; fl6 = 0
for (k, f) in [((5, 2), (1, 3)), ((7, -4), (1, 0))]:
    hs = bs
    long_h = 0j; long_hf = 0j
    for h in hs:
        if N(h) <= Hc or N(h)**3 > 2*Xc: continue
        if mu(h) == 0: continue
        Ph = Psi(h, k, f, xia)
        long_h += mu(h)*al(h).conjugate()**6*Ph**3/N(h)*Tcm(Xc/N(h)**3, k, f, xia)
        long_hf += mu(h)*al(h).conjugate()**3*Ph**3/N(h)*Tcm(Xc/N(h)**3, k, f, xia, bexp=3)
    # b-sum
    bsum = 0j
    for b_ in bs:
        Nb = N(b_)
        if Nb <= Hc or Nb**3 > 2*Xc: continue
        s_mu = sum(muu(h) for h in [hh for hh in bs if N(hh) <= Nb] if divides(h, b_) and N(h) > Hc)
        beta0 = al(b_).conjugate()**6*xia(b_)**3*(1 if coprime(b_, f) else 0)*s_mu
        if beta0 == 0: continue
        Lb = Xc/Nb**3
        col = sum(al(n).conjugate()**2*g2(n)*xia(n)*chi(n, k)*chi(n, f)**4*complex(Wc(N(n)/Lb)) for n in sqX if 1 < N(n)/Lb < 2)
        bsum += beta0*chi(b_, k)**3/Nb*col/math.sqrt(Lb)
    err6 = max(err6, abs(long_h - bsum)); fl6 = max(fl6, abs(long_hf - bsum))
    print(f"      C6 case k={k} f={f}: |P_long| = {abs(long_h):.3g}, error {abs(long_h - bsum):.2g}")
line("C6 cube removal: P_long = sum_b beta_0(b,f) chi_b(k)^3 N(b)^{-1} L_b^{-1/2} (column sum at L_b)", err6, 1e-10)
fals("conj(alpha)^3 (zeta) in the inversion and in T", fl6, 1e-8)

# ------------------------------------------------------------------ C7: the Hecke family for [GL]
def quad(n, x): return chi(n, x)**3 if n != (1, 0) else 1.0
def kappa_lam(x):          # nontrivial character mod lambda: O/lambda = F_3, x = a + b omega = a + b mod lambda
    r = (x[0] + x[1]) % 3
    return 0 if r == 0 else (1 if r == 1 else -1)
sq7 = [n for n in primaries(600, sqfree=True)]
bad_u = 0
for k in sq7:
    ek = 0 if N(k) % 4 == 1 else 1
    for u_ in UNITS:
        val = quad(k, u_)*(kappa_lam(u_)**ek)
        bad_u += abs(val - 1) > 1e-12
line("C7a (x/k)_2 kappa_lambda(x)^{e_k} is trivial on the six units", bad_u, 0, f"{len(sq7)} squarefree k")
# reciprocity factor depends only on (k mod 4, n mod 4)
tab = {}; bad_r = 0
for _ in range(3000):
    k, n = random.choice(sq7), random.choice(sq7)
    if not coprime(k, n): continue
    R = quad(k, n)*quad(n, k)        # (n/k)_2 (k/n)_2, real +-1
    key = ((k[0] % 4, k[1] % 4), (n[0] % 4, n[1] % 4))
    if key in tab and abs(tab[key] - R) > 1e-9: bad_r += 1
    tab.setdefault(key, R)
line("C7b quadratic reciprocity factor (n/k)_2 (k/n)_2 depends only on k, n mod 4", bad_r, 0, f"{len(tab)} class pairs")
# product character for k1 = k2 mod 4: trivial on units
bad_p = 0; cntp = 0
for _ in range(2000):
    k1, k2 = random.choice(sq7), random.choice(sq7)
    if k1 == k2 or not coprime(k1, k2) or red(k1) == red(k2) and False: continue
    if (k1[0] - k2[0]) % 4 or (k1[1] - k2[1]) % 4: continue
    cntp += 1
    for u_ in UNITS:
        bad_p += abs(quad(k1, u_)*quad(k2, u_) - 1) > 1e-12
line("C7c k1 = k2 mod 4: (x/k1 k2)_2 trivial on units (so the product is primitive mod k1 k2)", bad_p, 0, f"{cntp} pairs")
# even columns: (-2 n'/k)_2 = (-2/k)_2 (n'/k)_2
bad_e = 0
for _ in range(500):
    k, n = random.choice(sq7), random.choice(sq7)
    bad_e += abs(quad(k, mul((-2, 0), n)) - quad(k, (-2, 0))*quad(k, n)) > 1e-12
line("C7d even columns: (-2n'/k)_2 = (-2/k)_2 (n'/k)_2 with |(-2/k)_2| = 1", bad_e, 0)

# ------------------------------------------------------------------ C8: local costs in the squarefree-row lemma
bad8 = 0
for vs in (0, 1):                       # v_p(s0): p divides the fixed part s0 of the row
    for vg in range(0, 13):
        if vs == 0 and vg == 0: continue
        j = (vs + 4*vg) % 6
        cost = {0: 1, 4: 1}.get(j, 2)    # exponent of q in the cost of a^2 Y (q for j in {0,4}, q^2 otherwise)
        allowed = 2*vs + vg if vs == 1 else vg
        if vs == 1: allowed = 2              # p | s0: cost <= N(p)^2 (and N(s0)^2 pays it)
        bad8 += cost > allowed
line("C8 local cost exponent <= 2 v_p(s0) or v_p(g) (as in a^2 Y << H^2 N(g)/X)", bad8, 0)

# ------------------------------------------------------------------ C10: scale identities of the second transfer (sympy)
import sympy as sp
Hs_, L_, F_, Sg_, cs_, NC, Nt, Ne0, Ng, Ne, Nw, Nh2 = sp.symbols('H L F Sigma c_* NC Nt Ne0 Ng Ne Nw Nh2', positive=True)
ell_ = L_/(NC*Nt); w_ = Hs_*NC/(Ne0*L_**2*F_); Y_ = cs_*Sg_*L_*F_*Ne0/(Hs_*NC**2)
Nr = Nt*Ng/Ne; Nf = NC*Ne*Nw; Nk = Ne0*Ne*Nh2
chk = [sp.simplify(ell_/(Ng*Nw) - L_/(Nr*Nf)),
       sp.simplify(w_*Y_*Ng/(Ne*ell_) - cs_*Sg_*Nr/L_**2),
       sp.simplify(Y_*Nh2*Ng**2/(Ne*ell_**2) - cs_*Sg_*F_*Nk*Nr**2/(Hs_*L_)),
       sp.simplify(w_*Y_*ell_ - cs_*Sg_/(NC**2*Nt)),
       sp.simplify((Sg_*Nr/L_**2)*Nr*(L_/Nr)**2 - Sg_)]
line("C10 scale identities: X0' = L/(N(r)N(f')), weights c_* Sigma N(r)/L^2, Phihat argument, zero-frequency weight, block size",
     sum(1 for c in chk if c != 0), 0, "sympy")
# first-transfer normalization: (1/(LF)) * H/(N(e0) N(t) l) = w_{C,e0}
chk2 = sp.simplify(Hs_/(L_*F_)/(Ne0*Nt*ell_) - w_)
line("C10b first transfer: (1/(LF)) H/(N(e0) N(t) ell) = w_{C,e0}", 0 if chk2 == 0 else 1, 0, "sympy")

# ------------------------------------------------------------------ C11: multiplicity of (f, h) -> y = h f^2
bad11 = 0; mx = 0
fs_ = primaries(200, sqfree=True)
cnt = {}
for f in fs_:
    f2 = mul(f, f)
    for a in range(-12, 13):
        for b in range(-12, 13):
            h = (a, b)
            if h == (0, 0): continue
            y = mul(h, f2)
            cnt[y] = cnt.get(y, 0) + 1
for y, c in cnt.items():
    # tau(y) for the part coprime to 6 bounds the number of squarefree primary f with f^2 | y
    nf = sum(1 for f in fs_ if divides(mul(f, f), y))
    bad11 += c > nf; mx = max(mx, c)
line("C11 the map (f, h) -> h f^2 has multiplicity <= #{f squarefree : f^2 | y} <= tau(y)", bad11, 0, f"max multiplicity {mx}")
sys.exit(1 if NFAILED[0] else 0)
