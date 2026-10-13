"""t2b: end-to-end test of the intermediate form eq:b2-Sxi (7.32) of the first Poisson reduction, with Z of eq:d045 (7.29).

    M_D = D^-1 sum_{u in O} Phi(N(u)/H) |A_u(D)|^2,   A_u(D) = sum_n mu_K(n) alpha(n) chi_n(u) W(N(n)/D)   (eq:Au)

is compared with Z + sum_xi c_xi S_xi, where S_xi is computed in the (g, e0, h, z1, z2) variables of eq:b2-Sxi, i.e.
BEFORE the change of variables (r0, f, k) (the final form eq:b2-output is tested in tests_cm/proofs_2026-10-12/b2).
The sum over xi is taken literally: K(z1, z2) = sum_xi c_xi xi(z1) conj(xi(z2)) with c_xi of eq:d037 and the
independent character table of t2lib.  The identity holds for every radial Schwartz Phi; we use the Gaussian
Phi(t) = exp(-pi t), Phihat(t) = (2/sqrt3) exp(-4 pi t/3) (part B, B10), and a complex smooth W on (1, 2).
The h-sums are truncated where Phihat < 1e-20 relative (argument > 12), the u-sum at N(u) <= 14 H.

Also tested:
  eq:d045 (equality): the h = 0 terms of the Poisson formula, (H/D^2) sum_g sum_{e0|g} mu(e0) N(g)/N(e0) Phihat(0)
          |W0(N(g)/D)|^2 * K(1,1), equal Z = (H/D) Phihat(0) sum_g |W(N(g)/D)|^2 prod_{pi|g}(1 - 1/N(pi)).
  eq:d045 (bound, size check): Z/(H ||W||_inf^2) for D = 30 ... 5000 stays bounded.
  eq:d047: the W0 identity numerically on all (g, z1, z2) used (symbolically in t2c).
Falsifiers (must fail): zeta coefficients conj(alpha) gamma_2; twist h e0 instead of h e0^5; c_xi with G(t) for G(t^-1);
  N(g) for N(g)/N(e0); mu(e0) dropped.
Run: python3 -B t2b_sxi_endtoend.py [D H ...]   (default: D = 120 with H = 170 and 40, and D = 300 with H = 360; about 20 s)
"""
import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2lib import *
T0 = time.time()
SQ3 = math.sqrt(3.0)
args = [float(a) for a in sys.argv[1:]]
CASES = [(args[i], args[i + 1]) for i in range(0, len(args), 2)] if args else [(120.0, 170.0), (120.0, 40.0), (300.0, 360.0)]

def Wc(y): return math.exp(-1/((y - 1)*(2 - y)))*(1 + 0.5j*y) if 1 < y < 2 else 0.0
def W0(x): return x**-0.5*Wc(x).conjugate() if x > 0 else 0.0
Phi = lambda t: np.exp(-np.pi*t)
hatPhi = lambda t: 2/SQ3*np.exp(-4*np.pi*t/3)
X = CHARS
CXI_bad = (X.conj() @ GVEC)/432.0                     # falsifier: G(t) in place of G(t^-1)

def lattice_arrays(Rmax, include0=False):
    R = int(2*math.sqrt(Rmax)) + 3
    A, B = np.meshgrid(np.arange(-R, R + 1, dtype=np.int64), np.arange(-R, R + 1, dtype=np.int64), indexing='ij')
    A = A.ravel(); B = B.ravel(); NN = A*A - A*B + B*B
    k = (NN <= Rmax) & ((NN > 0) | include0)
    return A[k], B[k], NN[k]

def run(D, H):
    t0 = time.time()
    S_all = primaries(2*D, sqfree=True)
    cols = [n for n in S_all if D < N(n) < 2*D]
    # ---- left side
    UA, UB, NU = lattice_arrays(14*H, include0=True)
    Au = np.zeros(len(UA), complex)
    for n in cols:
        Au += muk(n)*alpha_c(n)*Wc(N(n)/D)*chi_vec(n, UA, UB)
    MD = float(np.sum(Phi(NU/H)*np.abs(Au)**2))/D
    # ---- Z (eq:d045) and the direct h = 0 terms
    Z = H/D*hatPhi(0.0)*sum(abs(Wc(N(g)/D))**2*math.prod(1 - 1/N(p) for p in set(factor(g))) for g in cols)
    K11 = np.sum(CXI)                                   # K(1,1) = sum_xi c_xi = G(1)
    Zdir = 0j
    for g in cols:
        for e0 in divisors_sqfree(g):
            Zdir += H/D**2*muk(e0)*N(g)/N(e0)*hatPhi(0.0)*abs(W0(N(g)/D))**2*K11
    # ---- S in the (g, e0, h, z1, z2) variables
    var = ["cm", "zeta", "e0", "Gswap", "noNe0", "nomu"]
    S = dict.fromkeys(var, 0j)
    e47 = 0.0; nterms = 0
    for g in S_all:
        zs = [z for z in S_all if coprime(z, g) and D < N(g)*N(z) < 2*D]
        if not zs: continue
        nz = len(zs); Nz = np.array([N(z) for z in zs], float)
        idx = [cls(z) for z in zs]
        Xz = X[:, idx]
        cop = np.array([[1.0 if coprime(a, b) else 0.0 for b in zs] for a in zs])
        Kc = ((CXI[:, None]*Xz).T @ Xz.conj())*cop           # K[i, j] = sum_xi c_xi xi(z_i) conj(xi(z_j)), (z_i, z_j) = 1
        Kb = ((CXI_bad[:, None]*Xz).T @ Xz.conj())*cop
        w0 = np.array([W0(N(g)*N(z)/D) for z in zs])
        Vc = np.array([Acm(z) for z in zs])*w0
        Vz = np.array([Azeta(z) for z in zs])*w0
        # eq:d047 numerically
        for i, a in enumerate(zs):
            for j, b in enumerate(zs):
                l = Wc(N(g)*N(a)/D).conjugate()*Wc(N(g)*N(b)/D)/math.sqrt(N(a)*N(b))
                r = N(g)/D*W0(N(g)*N(a)/D)*W0(N(g)*N(b)/D).conjugate()
                e47 = max(e47, abs(l - r))
        for e0 in divisors_sqfree(g):
            Rmax = 12.0*N(e0)*Nz.max()**2/H
            HA, HB, NH = lattice_arrays(Rmax)
            nterms += nz*nz*len(NH)
            e5 = pw(e0, 5)
            xA = HA*e5[0] - HB*e5[1]; xB = HA*e5[1] + HB*e5[0] - HB*e5[1]           # h e0^5
            yA = HA*e0[0] - HB*e0[1]; yB = HA*e0[1] + HB*e0[0] - HB*e0[1]           # h e0 (falsifier)
            CH = np.array([chi_vec(z, xA, xB) for z in zs])                          # [i, h]
            CHf = np.array([chi_vec(z, yA, yB) for z in zs])
            HP = hatPhi(H*NH[None, None, :]/(N(e0)*Nz[:, None, None]*Nz[None, :, None])) if nz*nz*len(NH) < 4e6 else None
            def bil(K, V, C):
                P = V[:, None]*C
                if HP is not None:
                    return np.einsum('ij,ih,jh,ijh->', K, P, P.conj(), HP)
                tot = 0j
                for i in range(nz):
                    hp = hatPhi(H*NH[None, :]/(N(e0)*Nz[i]*Nz[:, None]))
                    tot += np.sum(K[i, :, None]*P[i][None, :]*P.conj()*hp)
                return tot
            pre = H/D**2*muk(e0)*N(g)/N(e0)
            tc = bil(Kc, Vc, CH)
            S["cm"] += pre*tc
            S["zeta"] += pre*bil(Kc, Vz, CH)
            S["e0"] += pre*bil(Kc, Vc, CHf)
            S["Gswap"] += pre*bil(Kb, Vc, CH)
            S["noNe0"] += H/D**2*muk(e0)*N(g)*tc
            S["nomu"] += H/D**2*N(g)/N(e0)*tc
    tot = {k: Z + v for k, v in S.items()}
    err = abs(MD - tot["cm"])/abs(MD)
    print(f"  D={D:g} H={H:g}: {len(cols)} columns; M_D = {MD:.12f}; Z = {Z:.12f}; Z + sum_xi c_xi S_xi = {tot['cm'].real:.12f}"
          f"{tot['cm'].imag:+.1e}i; {nterms:.3g} (z1,z2,h) terms; {time.time() - t0:.1f} s")
    report(f"eq:b2-Sxi(D={D:g},H={H:g})", err, 1e-12, f"M_D = Z + sum_xi c_xi S_xi in the (g,e0,h,z1,z2) form, literal xi-sum; relative error")
    names = {"zeta": "zeta coefficients conj(alpha) gamma_2", "e0": "twist chi(h e0) for chi(h e0^5)",
             "Gswap": "c_xi with G(t) for G(t^-1)", "noNe0": "N(g) for N(g)/N(e0)", "nomu": "mu(e0) dropped"}
    for k, nm in names.items():
        falsifier(f"eq:b2-Sxi(D={D:g},H={H:g})", nm, abs(MD - tot[k])/abs(MD), 1e-6)
    ez = abs(Z - Zdir)/abs(Z)
    report(f"eq:d045(D={D:g},H={H:g})", ez, 1e-13, f"Z of eq:d045 = the h = 0 Poisson terms summed over e0 | g with K(1,1) = sum c_xi; relative error")
    report(f"eq:d047(D={D:g})", e47, 1e-13, "W0 identity on every (g, z1, z2) of the sum, absolute error")
    return err

for (D, H) in CASES:
    run(D, H)

# ---- eq:d045 bound: size check of Z/(H ||W||^2) = (2/sqrt3) D^-1 sum_g |W(N(g)/D)|^2 prod(1-1/N(pi)) / ||W||_inf^2
Wsup = max(abs(Wc(1 + i/10000)) for i in range(1, 10000))
SQBIG = primaries(10000, sqfree=True)
ratios = []
for D in (30, 120, 500, 2000, 5000):
    s = sum(abs(Wc(N(g)/D))**2*math.prod(1 - 1/N(p) for p in set(factor(g))) for g in SQBIG if D < N(g) < 2*D)
    ratios.append((D, 2/SQ3*s/D/Wsup**2))
print("  eq:d045 size check Z/(H ||W||_inf^2): " + ", ".join(f"D={D}: {r:.4f}" for D, r in ratios))
rmax = max(r for _, r in ratios)
report("eq:d045-bound", 0.0 if rmax < 1.0 else rmax, 0.0, f"size check only: Z/(H||W||^2) <= {rmax:.3f} for D = 30..5000 (bounded, as Z << H||W||^2 says)")
print(f"t2b done in {time.time() - T0:.1f} s")
