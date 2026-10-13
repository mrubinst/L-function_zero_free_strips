# b2 tests: new identities and inequalities in b2.tex, with suggested tests

Conventions as in the guide and in `tests_cm/part_B_lines_0550_0997/cmlib.py` (primary = 1 mod 3, lambda = 1 + 2 omega,
chi_n(x) = (x/n)_6 with zero values kept, e(z) = exp(4 pi i Im z/sqrt3), alpha(n) = n/|n|,
a_xi = conj(alpha)^2 gamma_2 xi). "Existing" names the claim ids of `tests_cm/part_B_lines_0550_0997` (B..) and the
earlier independent scripts.

| id | label in b2.tex | statement | class | test | existing coverage | run by b2 |
|---|---|---|---|---|---|---|
| b2-T1 | eq:b2-g1g2 | conj(alpha(n)) gamma_2(n) gamma_1(n) = mu_K(n) G(n) (OAI (4.4)), squarefree primary n coprime to 6 | E | brute-force Gauss sums, 150 squarefree n; falsifier: alpha instead of conj(alpha) | an earlier independent script ("Lemma 2.1(3)", 80 n); QF scripts (eq. (4.12)) | no (covered) |
| b2-T2 | eq:b2-quotient | chi_x(-1) conj(G(x)) G(y) R(x,y) = G(y x^{-1}) on frak G, with R(a,b) = chi_b(a)/chi_a(b) for coprime primes and chi_x(-1) = (-1)^{(N(x)-1)/6} | E | evaluate G on the 12 classes mod 4 from squarefree representatives; check for all pairs of coprime squarefree (x,y) with N <= 600; falsifier: drop chi_x(-1) | implicitly in B27, B73 (paired identities) and t6_identities.py | implicitly (T4) |
| b2-T3 | eq:b2-remove | E_(r) <= tau(r) sum_{j\|r} E(H, X/N(j), F N(j)) | I | the identity eq:exclusionid and the inequality on computed averages | section-10 checks (README of the revised Section 10: exclusion-removal identity to 3e-15, inequality on computed averages) | no (covered) |
| b2-T4 | eq:b2-forward | CM forward paired identity with an arbitrary class character xi_1; class form: mu mu conj(alpha(z1)) alpha(z2) gamma(conj chi_z1 chi_z2) = G(z2 z1^{-1}) A(z1) conj(A(z2)), A = conj(alpha)^2 gamma_2 | E | class form for 80 random coprime squarefree pairs (N(z1 z2) <= 5000), also multiplied by xi_1(z1)conj(xi_1(z2)) with xi_1(n) = chi_n(omega); falsifier: zeta coefficient conj(alpha) gamma_2 | B27 (xi_1 = 1, literal sum over all 432 xi); t6_identities.py (class form) | YES: max err 1.65e-15; falsifier 2.0 |
| b2-T5 | eq:b2-reverse | a_xi(z1) conj(a_xi(z2)) gamma(chi_z1 conj chi_z2) = sum_{xi1} l_{xi1}(xi) beta_{xi1}(z1) conj(beta_{xi1}(z2)), beta = mu_K conj(alpha) xi1; class form: = mu mu conj(alpha(z1)) alpha(z2) G(z1 z2^{-1}) | E | class form, 80 pairs; falsifier: beta = mu_K alpha xi1 | B73 (literal xi expansion) | YES: max err 1.54e-15; falsifier 2.0 |
| b2-T6 | eq:b2-twisted | a_xi(ab) = a_xi(a) a_xi(b) chi_b(a)^4, \|a_xi(a)\| = 1 | E | random coprime pairs, two xi; falsifier chi_b(a)^2 | B30; t6_identities.py | no (covered) |
| b2-T7 | eq:b2-crt | gamma(conj chi_z1 chi_z2) = R(z1,z2) gamma_{-1}(z1) gamma_1(z2), and gamma(chi_z1 conj chi_z2) = R(z1,z2) gamma_1(z1) gamma_{-1}(z2), with R(z1,z2) = chi_z2(z1) conj(chi_z1(z2)) = +-1 | E | brute-force Gauss sums for 80 coprime pairs; check R real and +-1; falsifier: drop R | QF (4.27) checks; implied by B27/B73 | implicitly (T4, T5) |
| b2-T8 | eq:b2-two | mu gamma_{-1} = chi_n(-1) conj(G) a^zeta and mu gamma_1 = G conj(a^zeta), a^zeta = conj(alpha) gamma_2 | E | brute force, 150 n | B23 (first); second is T1 rearranged | no (covered) |
| b2-T9 | Normalizations paragraph | gamma(psi) with modulus generator zeta q equals psi(zeta) gamma(psi); gamma_2(-pi) = gamma_2(pi), alpha(-pi) = -alpha(pi), so gamma_2(-pi)^3 = +alpha(-pi) (no Moebius sign for the Ireland-Rosen normalization) | E | brute force for the 12 smallest primary primes; gamma(psi) for 20 coprime pairs and the 6 units | none | YES: gamma_2(-pi) = gamma_2(pi) to 1.4e-15; gamma_2(-pi)^3 = alpha(-pi) to 3.6e-15; generator change of gamma(psi), 20 pairs x 6 units, to 4.4e-15 |
| b2-T10 | Step 1 | mu(n1)mu(n2) = mu(z1)mu(z2) and conj(alpha(n1))alpha(n2) = conj(alpha(z1))alpha(z2) for n_i = g z_i | E | random squarefree pairs | B28 | no (covered) |
| b2-T11 | Step 1, Z | zero-frequency term Z = (H/D) Phi^(0) sum_g \|W(N(g)/D)\|^2 prod_{pi\|g}(1 - 1/N(pi)) (CM coefficient contributes \|mu alpha\|^2 = 1) | E | part of the end-to-end test (T13) | B17 end-to-end (Poisson step only) | YES (inside T13) |
| b2-T12 | Step 3 | the change of variables (g,e0,v0,h) <-> (r0,f,k): bijection, g v0 = r0 f, h e0^5 v0^4 = k f^4, mu(e0)mu(v0) = mu(f), N(g)/N(e0) = N(r0), N(h)/(N(e0)N(v0)^2) = N(k)/N(f)^2 | E | enumerate all tuples with N(g v0) <= 300, N(h) <= 50 and compare with the triples; falsifier f = e0 v0^2 | B31, B32 (identities on random data; bijection not enumerated) | implicitly (T13) |
| b2-T13 | eq:b2-Sxi, eq:b2-output | END TO END: M_D = (1/D) sum_u Phi(N(u)/H) \|A_u(D)\|^2 equals Z + sum_xi c_xi S_xi with S_xi in the final (r0,f,k) form eq:b2-output (sum_xi c_xi xi(m1) conj(xi(m2)) = G(m2 m1^{-1}) on classes) | E | script below, Gaussian Phi(t) = exp(-pi t), Phi^(t) = (2/sqrt3) exp(-4 pi t/3), complex smooth W on (1,2); falsifiers: zeta coefficients, chi(k f^2), no exclusion (m_i,r0)=1, G(t) instead of G(t^{-1}) | none (B17 checks only the Poisson step; B32/B33 check the algebra symbolically) | YES: D=120, H=170: rel err 1.2e-15; D=120, H=40: rel err 3.0e-15; falsifiers 4.7e-2, 4.1e-4, 3.2e-3, 5.2e-2 (H=170) and 8.7e-2, 9.7e-5, 8.1e-2, 8.3e-2 (H=40) |
| b2-T14 | Step 4 | N(r0)N(f) <= b0 D; H N(k) <= C_Phi b0^2 D^2/N(r0)^2; X = D/(BF), H_rows = C1 D^2/(H B^2) satisfy eq:b2-ranges with C = C1 = max(1, C_Phi b0^2) | I | sympy, as in B32 | B32 | no (covered) |
| b2-T15 | Step 5 | eta = H N(k)/(N(f)^2 X^2) <= C1 for N(k) <= H_rows and N(f) >= F; varrho = N(r0)N(f)/(BF) in [1,4); kernel support [a0/4, b0]^2 for W supported in [a0, b0] | I | sympy: eta <= C1 F^2/N(f)^2 <= C1; sampling of the support | none | no (suggested; one-line algebra) |
| b2-T16 | Step 5 | bookkeeping sum_{N(r0)~B} (H N(r0)/D^2)(D/B)^2 << H, so M(D,H) << D^{2+theta+eps} | I | sympy | B33 | no (covered) |
| b2-T17 | Prop. b2-reduction | the analytic inequality itself (Lemma B.2 separation, hypothesis eq:b2-hyp) | A | not testable on finite data beyond sanity (B33 bookkeeping, B85 sanity of E/Sigma) | B33 | no |

Notes.
- T13 is the main new test. Its Gaussian Phi is not >= 1 on [0,1] and its transform is not compactly supported; the
  identity M_D = Z + sum c_xi S_xi holds for every radial Schwartz Phi, so this tests the exact algebra of Steps 1-3
  (Poisson with the e0-exclusion, the CM paired identity, Moebius inversion over v0, the change of variables, the
  fourth power k f^4 and the exclusion (m_i, r0) = 1). The k-sums are truncated where Phi^ < 1e-16 relative.
- At D = 40 the k f^2 falsifier does not fail (only one column m per (r0, f) with f != 1, so the character cancels in
  |.|^2); use D >= 100.
- Suggested home: a new script `t8_b2_reduction.py` in `tests_cm/part_B_lines_0550_0997/` (or the Section-7 folder of
  the published test suite), registering T4, T5, T9, T13 as claims; run time about 2 s.

## Script used (run as `python3 b2_reduction_check.py <path to tests_cm/part_B_lines_0550_0997> 120 170`)

```python
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
```

## Normalization check (T9)

```python
import sys, math, cmath
sys.path.insert(0, sys.argv[1])
from cmlib import *
# gamma_2 with modulus generator m (ideal (n), n primary): N^{-1/2} sum_{x mod n} chi_n(x)^2 e(x/m)
def g2_gen(n, m):
    mm = cx(m); s = 0
    for x in residues(n):
        c = chi(n, x)
        if c != 0: s += c**2*cmath.exp(4j*math.pi*(cx(x)/mm).imag/math.sqrt(3))
    return s/math.sqrt(N(n))
mx = 0; my = 0
for p in [q for q in PR_ODD if N(q) < 200][:12]:
    a = g2_gen(p, p); b = g2_gen(p, (-p[0], -p[1]))
    mx = max(mx, abs(a - b))
    my = max(my, abs(b**3 - (-alpha_c(p))), abs(a**3 + alpha_c(p)))
print("gamma_2(-pi) = gamma_2(pi):", mx, " gamma_2(pi)^3 = -alpha(pi) = alpha(-pi):", my)
# gamma(psi) with modulus generator zeta*q equals psi(zeta) gamma(psi), psi = conj(chi_z1) chi_z2, q = z1 z2
import random
random.seed(3)
sq = [n for n in primaries(200, sqfree=True) if n != (1, 0)]
mz = 0; cnt = 0
while cnt < 20:
    z1, z2 = random.choice(sq), random.choice(sq)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 2500: continue
    cnt += 1; q = mul(z1, z2)
    def gam_gen(m):
        mm = cx(m); s = 0
        for x in residues(q):
            v = chi(z1, x).conjugate()*chi(z2, x)
            if v != 0: s += v*cmath.exp(4j*math.pi*(cx(x)/mm).imag/math.sqrt(3))
        return s/math.sqrt(N(q))
    g0 = gam_gen(q)
    for u in UNITS:
        psi_u = chi(z1, u).conjugate()*chi(z2, u)
        mz = max(mz, abs(gam_gen(mul(u, q)) - psi_u*g0))
print("gamma(psi) with generator zeta q = psi(zeta) gamma(psi), 20 pairs x 6 units:", mz)
```

Output (Oct 12 2026):

```
PAIRED forward class form: max err 1.65e-15 over 80 pairs; zeta-coefficient falsifier 2.00e+00
PAIRED forward with xi_1 = chi_n(omega): max err 1.61e-15
PAIRED reverse class form: max err 1.54e-15; falsifier mu alpha (not conj alpha): 2.00e+00
REDUCTION D=120 H=170: 34 columns, 109 (r0,f) pairs; M_D = 0.0070864918, Z = 0.006963, Z + sum c S = 0.0070864918+0.0000000000j; rel err 1.22e-15
   falsifier zeta: rel err 4.73e-02
   falsifier kf2: rel err 4.14e-04
   falsifier noexcl: rel err 3.19e-03
   falsifier Gswap: rel err 5.17e-02
REDUCTION D=120 H=40: 34 columns, 109 (r0,f) pairs; M_D = 0.0016153990, Z = 0.001638, Z + sum c S = 0.0016153990-0.0000000000j; rel err 2.95e-15
   falsifier zeta: rel err 8.66e-02
   falsifier kf2: rel err 9.65e-05
   falsifier noexcl: rel err 8.12e-02
   falsifier Gswap: rel err 8.30e-02
gamma_2(-pi) = gamma_2(pi): 1.37e-15  gamma_2(pi)^3 = -alpha(pi) = alpha(-pi): 3.58e-15
gamma(psi) with generator zeta q = psi(zeta) gamma(psi), 20 pairs x 6 units: 4.43e-15
```
