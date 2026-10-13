"""t5: reflection of a CM plain sum on an actual row, with deleted and restored Euler factors
(inserts.tex, Lemma s5-rowchar and the CM form of OAI78 L12733-12761), and mask deletion (2.1a) with alpha.

Row character psi_k(n) = tau(n) chi_n(k) on primary n outside S = {lambda, 2}; CM plain coefficient alpha(l)^e psi_k(l).
xi_k = primitive Hecke character inducing alpha^e psi_k^0 (conductor f_k), r_k = primes of S u supp(k) not dividing f_k.
Claim:  T_k(X;W) = eps(xi_k) sum_{d0 | r_k} sum_{rad h0 | r_k} mu(d0) xi_k(d0) conj(xi_k(h0)) (N d0 N h0)^{-1/2}
                   T_{conj xi_k^0}(C_k N d0 / (X N h0); W^sharp),   C_k = 3 N(f_k)/(4 pi^2),
and C_k N(r_k) <= 27 C_k^fin N(R_0k)  (C_k^fin, R_0k: OAI78's data for the finite-order part).
Falsifiers: OAI78's finite-order data (C_k^fin, R_0k) in place of (C_k, r_k); the J_0 dual profile.
"""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from s5lib import (report, falsifier, summary, chi_p, primes_upto, cx, N, mul, OM, divides, primary_assoc, conj)
from hecke import from_primary_function, primitive, gauss, elements

LAM = (1, 2); TWO = (-2, 0)
alpha = lambda a: cx(a) / abs(cx(a))
P = primes_upto(40)
p7 = [p for p in P if N(p) == 7][0]
p13 = [p for p in P if N(p) == 13][0]
p19 = [p for p in P if N(p) == 19][0]
p31 = [p for p in P if N(p) == 31][0]

def cubic_mod2(a):
    return {(1, 0): 1, (0, 1): OM, (1, 1): OM ** 2}[(a[0] % 2, a[1] % 2)]

def W(x): return np.sqrt(x) * np.exp(-x)
def Ws(y): return np.sqrt(y) * np.exp(-y)
def Ws_J0(y): return np.array([float(mp.gamma(1.5) * mp.hyp1f1(1.5, 1, -t)) for t in np.atleast_1d(y)])

# Row characters.  OAI78's psi_k(n) = tau(n) chi_n(k) is, by sextic reciprocity, tau'(n) chi_k(n) with a fixed-ray factor
# absorbed in tau'.  We use the computable form psi(n) = tau(n) prod_p chi_p(n)^{v_p} (residue symbols modulo the moving
# primes, zero extension retained, exponent 6 = redundant zero), which is the same coefficient class.
def make_psi(tau, kfac):
    def psi(a):
        v = tau(a)
        for p, ex in kfac:
            c = chi_p(a, p)
            if c == 0: return 0
            v *= c ** ex
        return v
    return psi

rows = [
    ("psi=chi_p7 chi_p13^6, e=+1", [(p7, 1), (p13, 6)], lambda a: 1, 1),
    ("psi=chi_p7 chi_p13^6, e=-1", [(p7, 1), (p13, 6)], lambda a: 1, -1),
    ("psi=(./2)_3 chi_p19^2 chi_p7^3, e=+1", [(p19, 2), (p7, 3)], cubic_mod2, 1),
    ("psi=chi_p31^5, e=-1", [(p31, 5)], lambda a: 1, -1),
    ("psi=(./2)_3 chi_p13^4, e=-1", [(p13, 4)], cubic_mod2, -1),
]

for name, kfac, tau, e in rows:
    kprimes = [p for p, ex in kfac]
    m = (6, 0)
    for p in kprimes: m = mul(m, p)
    mprimes = [LAM] * 2 + [TWO] + kprimes
    psi = make_psi(tau, kfac)
    def cm(a):
        return alpha(a) ** e * psi(a)
    # build on residues coprime to m, as Hecke characters (values on primary generators)
    ch_cm = primitive(from_primary_function(cm, m, mprimes, e))
    ch_fin = primitive(from_primary_function(psi, m, mprimes, 0))
    f = ch_cm.m; g = ch_fin.m
    fpr = set(ch_cm.primes); gpr = set(ch_fin.primes)
    allp = [LAM, TWO] + kprimes
    r = [p for p in allp if p not in fpr]              # r_k: zero primes of alpha^e psi_k^0 not dividing f_k
    R0 = [p for p in allp if p not in gpr]             # OAI78's R_{0,k}
    C = 3 * N(f) / (4 * math.pi ** 2); Cfin = 3 * N(g) / (4 * math.pi ** 2)
    Nr = math.prod(N(p) for p in r); NR0 = math.prod(N(p) for p in R0)
    eps = e * gauss(ch_cm) * (cx(f) / abs(cx(f))) ** e
    print(f"--- {name}: N(f_k)={N(f)}, N(g_k)={N(g)}, r_k norms={[N(p) for p in r]}, R_0k norms={[N(p) for p in R0]}")
    report(f"5-{name}-conductor-bound", max(0.0, C * Nr / (27 * Cfin * NR0) - 1), 0.0,
           f"C_k N(r_k) = {C*Nr:.3f} <= 27 C_k^fin N(R_0k) = {27*Cfin*NR0:.3f}")

    def T_direct(X, cut=60.0):
        A, B, Z, n = elements(int(cut * X) + 2)
        tot = 0
        sel = ((A - 1) % 3 == 0) & (B % 3 == 0) & ~((A % 2 == 0) & (B % 2 == 0))
        for a, b, nn in zip(A[sel], B[sel], n[sel]):
            v = cm((int(a), int(b)))
            if v != 0:
                tot += v * W(nn / X)
        return tot / math.sqrt(X)

    def xi_val(ch, d):
        z = cx(d); return ch(np.array([d[0]]), np.array([d[1]]))[0] * (z / abs(z)) ** ch.e

    def T_dual(ch, rlist, Y, prof, cut=60.0):
        A, B, Z, n = elements(int(cut * Y) + 2)
        vals = np.conj(ch(A, B) * (Z / np.abs(Z)) ** ch.e)
        for p in rlist:
            c_, d_ = p; x_, y_ = c_ - d_, -d_
            Pp = A * x_ - B * y_; Qp = A * y_ + B * x_ - B * y_
            vals = vals * ~((Pp % N(p) == 0) & (Qp % N(p) == 0))
        ns, idx = np.unique(n, return_inverse=True)
        return np.sum(vals * prof(ns / Y)[idx]) / 6.0 / math.sqrt(Y)

    def rhs(ch, rlist, CC, epsv, X, prof):
        tot = 0
        # h0 = product of powers of the primes in rlist, truncated where the reflected scale is negligible
        for k_ in range(len(rlist) + 1):
            for sub in itertools.combinations(rlist, k_):
                d0 = (1, 0)
                for p in sub: d0 = mul(d0, p)
                c0 = (-1) ** k_ * xi_val(ch, d0) / math.sqrt(N(d0))
                exps = [range(0, 40)] * len(rlist)
                def rec(i, h0):
                    nonlocal tot
                    if i == len(rlist):
                        Y = CC * N(d0) / (X * N(h0))
                        if Y < 1 / 80.0:
                            return
                        tot += c0 * np.conj(xi_val(ch, h0)) / math.sqrt(N(h0)) * T_dual(ch, rlist, Y, prof)
                        return
                    hh = h0
                    for ex in range(0, 60):
                        if CC * N(d0) / (X * N(hh)) < 1 / 80.0:
                            break
                        rec(i + 1, hh)
                        hh = mul(hh, rlist[i])
                rec(0, (1, 0))
        return epsv * tot

    worst = 0; f1 = 0; f2 = 0
    for X in [0.7, 2.0, 5.0]:
        lhs = T_direct(X)
        worst = max(worst, abs(lhs - rhs(ch_cm, r, C, eps, X, Ws)))
        if (C, Nr) != (Cfin, NR0):
            f1 = max(f1, abs(lhs - rhs(ch_cm, R0, Cfin, eps, X, Ws)))
        f2 = max(f2, abs(lhs - rhs(ch_cm, r, C, eps, X, Ws_J0)))
    report(f"5-{name}-reflection", worst, 1e-9, "OAI78 L12736-12744 in CM form: deleted/restored factors, C_k, r_k, eps(xi_k)")
    if (C, Nr) != (Cfin, NR0):
        falsifier(f"5-{name}-reflection", "OAI78's finite-order C_k and R_0k", f1, 1e-4)
    falsifier(f"5-{name}-reflection", "J_0 dual profile", f2, 1e-4)

    # mask deletion (2.1a) with alpha: extra mask R = (p19 or p7) * (inert 5)
    Rx = [q for q in [p19, (-5, 0), p7] if q not in kprimes][:2]
    def S_masked(X, cut=60.0):
        A, B, Z, n = elements(int(cut * X) + 2)
        sel = ((A - 1) % 3 == 0) & (B % 3 == 0) & ~((A % 2 == 0) & (B % 2 == 0))
        tot = 0
        for a, b, nn in zip(A[sel], B[sel], n[sel]):
            l = (int(a), int(b))
            if any(divides(q, l) for q in Rx): continue
            v = cm(l)
            if v != 0: tot += v * W(nn / X)
        return tot / math.sqrt(X)
    err = 0
    for X in [1.0, 4.0]:
        ie = 0
        for k_ in range(len(Rx) + 1):
            for sub in itertools.combinations(Rx, k_):
                d = (1, 0)
                for q in sub: d = mul(d, q)
                ie += (-1) ** k_ * cm(d) / math.sqrt(N(d)) * T_direct(X / N(d))
        err = max(err, abs(S_masked(X) - ie))
    report(f"5-{name}-mask-deletion", err, 1e-10, "(2.1a): S_{k,R}(X) = sum_{d|R} mu(d) alpha(d)^e psi_k^0(d) N(d)^{-1/2} T_k(X/N d)")
summary()
