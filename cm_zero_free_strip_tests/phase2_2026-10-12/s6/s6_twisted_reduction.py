"""s6-T1: end-to-end exact form of the first Poisson reduction (Proposition prop:b2-reduction of the guide) for the
TWISTED family of the insert (Appendix app:s6):
    A_u(D) = sum_n mu(n) alpha(n) nu(n) chi_n(u) W(N(n)/D),   n primary squarefree coprime to S = {2, 3, pi13},
with nu of order 12 modulo pi13 = 4 + 3 omega (as in s6_twisted_paired.py).  The claim tested is
    M_D = (1/D) sum_{u in O} Phi(N(u)/H) |A_u(D)|^2 = Z + sum_{xi in hat G_S} c_xi(conj nu) S_xi,
with S_xi in the final (r0, f, k) form (7.17) of the guide, every column-type variable (g, r0, f, m_i) coprime to S,
and the class form  sum_xi c_xi(conj nu) xi(m1) conj(xi(m2)) = conj(nu(m1)) nu(m2) G(m2 m1^{-1})  (checked literally
over the 144 characters of G_S in s6_twisted_paired.py).
Gaussian Phi(t) = exp(-pi t), hatPhi(t) = (2/sqrt3) exp(-4 pi t/3); the identity holds for any radial Schwartz Phi.
Falsifiers (must fail): orientation nu(m1) conj(nu(m2)); no twist on the dual side; r0, f not coprime to pi13;
zeta coefficients conj(alpha) gamma_2.
Usage: python3 s6_twisted_reduction.py [path to tests_cm/part_B_lines_0550_0997] [D] [H]
"""
import os, sys, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')
sys.path.insert(0, os.path.abspath(LIB))
from cmlib import *
D = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
H = float(sys.argv[3]) if len(sys.argv) > 3 else 170.0

PI = (4, 3); P13 = 13
W0 = (-PI[0]*pow(PI[1], -1, P13)) % P13
DLOG = {}; x = 1
for k in range(12): DLOG[x] = k; x = (2*x) % P13
def red13(x): return (x[0] + x[1]*W0) % P13
def nu(n):
    r = red13(n)
    return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)
def copS(n): return red13(n) != 0                 # coprime to pi13 (coprimality to 6 is built into primaries())

sqall = primaries(3000, sqfree=True)
Gcls = {}
for n in sqall:
    c = (n[0] % 4, n[1] % 4)
    if c not in Gcls: Gcls[c] = 1.0 + 0j if n == (1, 0) else G_fast(n)
def m4(x): return (x[0] % 4, x[1] % 4)
def mul4(x, y): return m4(mul(x, y))
cls4 = list(Gcls)
inv4 = {c: [d for d in cls4 if mul4(c, d) == (1, 0)][0] for c in cls4}
def Gq(a, b): return Gcls[mul4(m4(a), inv4[m4(b)])]          # G(a b^{-1})
g2c = {}
def g2(n):
    if n not in g2c: g2c[n] = 1.0 + 0j if n == (1, 0) else gamma_j_fast(2, n)
    return g2c[n]
def Acm(n): return alpha_c(n).conjugate()**2*g2(n)
def Azeta(n): return alpha_c(n).conjugate()*g2(n)
def muk(n): return 1 if n == (1, 0) else mu(n)

def Wc(y): return math.exp(-1/((y - 1)*(2 - y)))*(1 + 0.5j*y) if 1 < y < 2 else 0.0
def W0f(x): return x**-0.5*Wc(x).conjugate() if x > 0 else 0.0
Phi = lambda t: np.exp(-np.pi*t)
hatPhi = lambda t: 2/SQ3*np.exp(-4*np.pi*t/3)
Sall = primaries(2*D, sqfree=True)                         # squarefree primary coprime to 6
S = [n for n in Sall if copS(n)]                           # ... and coprime to pi13
cols = [n for n in S if D < N(n) < 2*D]
# ---------------- left side
R = int(2*math.sqrt(14*H)) + 3
UA, UB = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing='ij'); UA = UA.ravel(); UB = UB.ravel()
NU = UA*UA - UA*UB + UB*UB; keep = NU <= 14*H; UA, UB, NU = UA[keep], UB[keep], NU[keep]
Au = np.zeros(len(UA), complex)
for n in cols:
    Au += muk(n)*alpha_c(n)*nu(n)*Wc(N(n)/D)*chi_vec(n, UA, UB)
MD = np.sum(Phi(NU/H)*np.abs(Au)**2)/D
# ---------------- zero frequency
Z = H/D*hatPhi(0.0)*sum(abs(Wc(N(g)/D))**2*np.prod([1 - 1/N(p) for p in set(factor(g))] or [1.0]) for g in cols)
def coprime(a, b): return N(gcd_primary(a, b)) == 1
RHS = {"twist": 0j, "orient": 0j, "notwist": 0j, "noS": 0j, "zeta": 0j}
npairs = 0
for r0 in Sall:
    for f in Sall:
        if N(r0)*N(f) >= 2*D or not coprime(r0, f): continue
        rf_ok = copS(r0) and copS(f)
        ms_all = [m for m in S if D < N(r0)*N(f)*N(m) < 2*D and coprime(m, r0)]
        if not ms_all: continue
        npairs += rf_ok
        Kmax = 12*(2*D/N(r0))**2/H
        Rk = int(2*math.sqrt(Kmax)) + 3
        KA, KB = np.meshgrid(np.arange(-Rk, Rk + 1), np.arange(-Rk, Rk + 1), indexing='ij'); KA = KA.ravel(); KB = KB.ravel()
        NK = KA*KA - KA*KB + KB*KB; kk = (NK > 0) & (NK <= Kmax); KA, KB, NK = KA[kk], KB[kk], NK[kk]
        f4 = pw(f, 4)
        XA = KA*f4[0] - KB*f4[1]; XB = KA*f4[1] + KB*f4[0] - KB*f4[1]      # k f^4
        CH = np.array([chi_vec(m, XA, XB) for m in ms_all])
        Nm = np.array([N(m) for m in ms_all], float)
        arg = H*NK[None, None, :]/(N(f)**2*Nm[:, None, None]*Nm[None, :, None])
        HP = hatPhi(arg)
        for variant in RHS:
            if variant != "noS" and not rf_ok: continue
            coef = np.array([(Azeta(m) if variant == "zeta" else Acm(m))*W0f(N(r0)*N(f)*N(m)/D) for m in ms_all])
            V = coef[:, None]*CH
            if variant == "orient":
                Gm = np.array([[nu(m1)*nu(m2).conjugate()*Gq(m2, m1) for m2 in ms_all] for m1 in ms_all])
            elif variant == "notwist":
                Gm = np.array([[Gq(m2, m1) for m2 in ms_all] for m1 in ms_all])
            else:
                Gm = np.array([[nu(m1).conjugate()*nu(m2)*Gq(m2, m1) for m2 in ms_all] for m1 in ms_all])
            T = np.einsum('ij,ik,jk,ijk->', Gm, V, np.conj(V), HP)
            RHS[variant] += H/D**2*muk(f)*N(r0)*T
tot = {k: Z + v for k, v in RHS.items()}
err = abs(MD - tot["twist"])/abs(MD)
report("s6-T1", err, 1e-12, f"twisted reduction M_D = Z + sum_xi c_xi(conj nu) S_xi, D={D:g}, H={H:g}: {len(cols)} columns, "
       f"{npairs} (r0,f) pairs; M_D={MD:.10f}, Z={Z:.6f}")
for k, nm in (("orient", "nu(m1) conj nu(m2)"), ("notwist", "no twist on dual side"),
              ("noS", "r0, f not restricted to be coprime to pi13"), ("zeta", "zeta coefficients")):
    falsifier("s6-T1", nm, abs(MD - tot[k])/abs(MD), 1e-6)
