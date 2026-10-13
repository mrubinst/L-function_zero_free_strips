"""b2 checks: (1) class form of the CM paired identities (forward and reverse);
(2) end-to-end exact form of the CM first Poisson reduction:
    M_D = (1/D) sum_u Phi(N(u)/H) |A_u(D)|^2 = Z + sum_xi c_xi S_xi,
with S_xi in the (r0, f, k) variables (r0 = g/e0, f = e0 v0, k = e0 h), and
sum_xi c_xi xi(m1) conj(xi(m2)) = G(m2 m1^{-1}) evaluated on classes.
Gaussian Phi(t)=exp(-pi t), hatPhi(t)=(2/sqrt3)exp(-4 pi t/3) (the identity holds for any Schwartz Phi)."""
import os, sys, math, cmath, random
sys.path.insert(0, sys.argv[1])
from cmlib import *
random.seed(12)

D = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
H = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0

# ---- G on classes mod 4 (G depends only on n mod 4; checked in B25)
sqall = primaries(3000, sqfree=True)
Gcls = {}
for n in sqall:
    c = (n[0] % 4, n[1] % 4)
    if c not in Gcls:
        Gcls[c] = 1.0 + 0j if n == (1, 0) else G_fast(n)
assert len(Gcls) == 12, len(Gcls)
cls4 = list(Gcls)
def m4(x): return (x[0] % 4, x[1] % 4)
def mul4(x, y): return m4(mul(x, y))
inv4 = {c: [d for d in cls4 if mul4(c, d) == (1, 0)][0] for c in cls4}
def Gq(a, b):
    """G(a b^{-1}) on classes (a, b coprime to 2)."""
    return Gcls[mul4(m4(a), inv4[m4(b)])]

g2c = {}
def g2(n):
    if n not in g2c: g2c[n] = 1.0 + 0j if n == (1, 0) else gamma_j_fast(2, n)
    return g2c[n]
def Acm(n): return alpha_c(n).conjugate()**2*g2(n)        # conj(alpha)^2 gamma_2
def Azeta(n): return alpha_c(n).conjugate()*g2(n)         # zeta coefficient (falsifier)
def muk(n): return 1 if n == (1, 0) else mu(n)

# ---- (1) class form of the paired identities
sq = [n for n in sqall if N(n) <= 400]
e_f = e_r = e_xi = 0; f_f = f_r = 0; cnt = 0
xi1 = lambda n: complex(chi(n, (0, 1)))      # chi_n(omega): a character of the classes mod 72 (depends on N(n) mod 36)
while cnt < 80:
    z1, z2 = random.choice(sq), random.choice(sq)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 5000: continue
    cnt += 1
    gf = gauss_psi_general(z1, z2, conj_first=True)
    gr = gauss_psi_general(z1, z2, conj_first=False)
    lhs = muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*gf
    rhs = Gq(z2, z1)*Acm(z1)*Acm(z2).conjugate()
    e_f = max(e_f, abs(lhs - rhs))
    f_f = max(f_f, abs(lhs - Gq(z2, z1)*Azeta(z1)*Azeta(z2).conjugate()))
    # with a nontrivial class character xi_1 on the Moebius side
    lx = lhs*xi1(z1)*xi1(z2).conjugate()
    rx = xi1(z1)*xi1(z2).conjugate()*Gq(z2, z1)*Acm(z1)*Acm(z2).conjugate()
    e_xi = max(e_xi, abs(lx - rx))
    # reverse: A(z1) conj A(z2) gamma(chi_z1 conj chi_z2) = mu mu conj(alpha(z1)) alpha(z2) G(z1 z2^{-1})
    lr = Acm(z1)*Acm(z2).conjugate()*gr
    rr = muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*Gq(z1, z2)
    e_r = max(e_r, abs(lr - rr))
    f_r = max(f_r, abs(lr - muk(z1)*muk(z2)*alpha_c(z1)*alpha_c(z2).conjugate()*Gq(z1, z2)))
print(f"PAIRED forward class form: max err {e_f:.2e} over {cnt} pairs; zeta-coefficient falsifier {f_f:.2e}")
print(f"PAIRED forward with xi_1 = chi_n(omega): max err {e_xi:.2e}")
print(f"PAIRED reverse class form: max err {e_r:.2e}; falsifier mu alpha (not conj alpha): {f_r:.2e}")

# ---- (2) end-to-end reduction
def Wc(y): return math.exp(-1/((y - 1)*(2 - y)))*(1 + 0.5j*y) if 1 < y < 2 else 0.0
def W0(x): return x**-0.5*Wc(x).conjugate() if x > 0 else 0.0
Phi = lambda t: np.exp(-np.pi*t)
hatPhi = lambda t: 2/SQ3*np.exp(-4*np.pi*t/3)
S = [n for n in primaries(2*D, sqfree=True)]                     # all squarefree primary n coprime to 6, N < 2D
cols = [n for n in S if D < N(n) < 2*D]
# LHS
R = int(2*math.sqrt(14*H)) + 3
UA, UB = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing='ij'); UA = UA.ravel(); UB = UB.ravel()
NU = UA*UA - UA*UB + UB*UB; keep = NU <= 14*H; UA, UB, NU = UA[keep], UB[keep], NU[keep]
Au = np.zeros(len(UA), complex)
for n in cols:
    Au += muk(n)*alpha_c(n)*Wc(N(n)/D)*chi_vec(n, UA, UB)
MD = np.sum(Phi(NU/H)*np.abs(Au)**2)/D
# Z
Z = H/D*hatPhi(0.0)*sum(abs(Wc(N(g)/D))**2*np.prod([1 - 1/N(p) for p in set(factor(g))] or [1.0]) for g in cols)
def coprime(a, b): return N(gcd_primary(a, b)) == 1
def fourth(f): return pw(f, 4)
RHS = {"cm": 0j, "zeta": 0j, "kf2": 0j, "noexcl": 0j, "Gswap": 0j}
npairs = 0
for r0 in S:
    for f in S:
        if N(r0)*N(f) >= 2*D or not coprime(r0, f): continue
        ms_all = [m for m in S if D < N(r0)*N(f)*N(m) < 2*D]
        if not ms_all: continue
        npairs += 1
        Kmax = 12*(2*D/N(r0))**2/H
        Rk = int(2*math.sqrt(Kmax)) + 3
        KA, KB = np.meshgrid(np.arange(-Rk, Rk + 1), np.arange(-Rk, Rk + 1), indexing='ij'); KA = KA.ravel(); KB = KB.ravel()
        NK = KA*KA - KA*KB + KB*KB; kk = (NK > 0) & (NK <= Kmax); KA, KB, NK = KA[kk], KB[kk], NK[kk]
        f4 = fourth(f); f2 = pw(f, 2)
        XA = KA*f4[0] - KB*f4[1]; XB = KA*f4[1] + KB*f4[0] - KB*f4[1]      # k f^4 = (KA + KB w)(c + d w)
        YA = KA*f2[0] - KB*f2[1]; YB = KA*f2[1] + KB*f2[0] - KB*f2[1]      # k f^2 (falsifier)
        for variant in RHS:
            ms = ms_all if variant == "noexcl" else [m for m in ms_all if coprime(m, r0)]
            if not ms: continue
            coef = np.array([(Azeta(m) if variant == "zeta" else Acm(m))*W0(N(r0)*N(f)*N(m)/D) for m in ms])
            if variant == "kf2":                                     # falsifier: chi_m(k f^2)
                CH = np.array([chi_vec(m, YA, YB) for m in ms])
            else:
                CH = np.array([chi_vec(m, XA, XB) for m in ms])
            V = coef[:, None]*CH                                     # [m, k]
            if variant == "Gswap":
                Gm = np.array([[Gq(m1, m2) for m2 in ms] for m1 in ms])
            else:
                Gm = np.array([[Gq(m2, m1) for m2 in ms] for m1 in ms])       # G(m2 m1^{-1})
            Nm = np.array([N(m) for m in ms], float)
            arg = H*NK[None, None, :]/(N(f)**2*Nm[:, None, None]*Nm[None, :, None])
            T = np.einsum('ij,ik,jk,ijk->', Gm, V, np.conj(V), hatPhi(arg))
            RHS[variant] += H/D**2*muk(f)*N(r0)*T
tot = {k: Z + v for k, v in RHS.items()}
err = abs(MD - tot["cm"])/abs(MD)
print(f"REDUCTION D={D:g} H={H:g}: {len(cols)} columns, {npairs} (r0,f) pairs; M_D = {MD:.10f}, Z = {Z:.6f}, Z + sum c S = {tot['cm']:.10f}; rel err {err:.2e}")
for k in ("zeta", "kf2", "noexcl", "Gswap"):
    print(f"   falsifier {k}: rel err {abs(MD - tot[k])/abs(MD):.2e}")

# Exit status (added 2026-10-12): 1 unless every identity error above is <= 1e-10 and every falsifier error is > 1e-6.
# CM_TESTS_FORCE_FAIL=1 counts the first check as failed (self-test of the wiring).
checks = [("paired forward err <= 1e-10", e_f <= 1e-10), ("paired forward with xi_1 err <= 1e-10", e_xi <= 1e-10),
          ("paired reverse err <= 1e-10", e_r <= 1e-10), ("forward falsifier (zeta coefficients) > 1e-6", f_f > 1e-6),
          ("reverse falsifier (mu alpha) > 1e-6", f_r > 1e-6), ("reduction rel err <= 1e-10", err <= 1e-10)]
checks += [(f"reduction falsifier {k} > 1e-6", abs(MD - tot[k])/abs(MD) > 1e-6) for k in ("zeta", "kf2", "noexcl", "Gswap")]
import os
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    checks[0] = (checks[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, c in checks if not c]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(checks), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(checks)))
sys.exit(1 if bad else 0)
