"""t8_signal_init.py -- the CM signal steps of the second canonical Poisson transformation (eq:s4-a0) and of the
initialization (eq:s4-Astar) in Appendix E of cm-eleven-twelfths-guide.tex, in the form in which they are used:
after the ray factor G([z2][z1]^{-1}) is expanded over the characters theta of the ray group.

Claim (eq:s4-Astar; OAI78 L11803-11814 with the CM coefficient): for coprime squarefree primary z1, z2 outside S and the
Moebius class m = mu conj(alpha) nu_0 (the coefficient of R_0, eq:s4-R0),
   m(z1) conj(m(z2)) gamma(conj chi_z1 chi_z2) = sum_theta Ghat(theta) A_*(z1) conj(A_*(z2)),
   A_*(z) = conj(alpha(z))^2 gamma_2(z) nu_*(z),  nu_* = nu_0 conj(theta),  Ghat(theta) = |G|^{-1} sum_t G(t) conj theta(t),
the sum being over the 12 characters theta of (O/4)^x (G is a function of classes mod 4, (G1) of Appendix G).
eq:s4-a0 is the same switch in the second canonical Poisson transformation (OAI78 L10770-10807): the class
mu conj(alpha) nu_i turns into a_0 nu_i conj(theta) with a_0 = conj(alpha)^2 gamma_2 and the same ray factor; the
test with nu_0 = nu_i covers it.  Three twists nu_0: 1, chi_.(lambda) (modulus 9) and an order-12 character modulo
pi13 = 4 + 3 omega (then z1, z2 prime to pi13).

Also (eq:initial-common-extraction of OAI78 with A_*, "the extraction holds by twisted multiplicativity"):
   A_*(sn) conj(chi_sn(d')) chi_sn(h') = A_*(s) conj(chi_s(d')) chi_s(h') * A_*(n) chi_n(d'h') chi_n(d's)^4
for coprime squarefree s, n and arbitrary d', h' (zero extensions included).

Falsifiers (must fail): nu_* = nu_0 theta (theta not conjugated); the zeta class A_* = conj(alpha) gamma_2 nu_*;
Ghat built from G(t^{-1}); the extraction without the CRT factor chi_n(s)^4.
Usage: python3 t8_signal_init.py
"""
import os, sys, math, cmath, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import *
random.seed(13)

def muk(n): return 1 if n == (1, 0) else mu(n)
def al(n): return 1.0 + 0j if n == (1, 0) else alpha_c(n)
def phi(n): return al(n).conjugate()
_g2 = {}
def g2(n):
    if n == (1, 0): return 1.0 + 0j
    if n not in _g2: _g2[n] = gamma_j_fast(2, n)
    return _g2[n]
def chi0(n, x): return 1.0 + 0j if n == (1, 0) else complex(chi(n, x))

# ---- (O/4)^x = C6 x C2 and its 12 characters
def m4(x): return (x[0] % 4, x[1] % 4)
def mul4(x, y): return m4(mul(x, y))
U4 = [(a, b) for a in range(4) for b in range(4) if not (a % 2 == 0 and b % 2 == 0)]
assert len(U4) == 12
def order4(g):
    k, y = 1, g
    while y != (1, 0): y = mul4(y, g); k += 1
    return k
g1 = next(g for g in U4 if order4(g) == 6)
P1 = [(1, 0)]
for _ in range(5): P1.append(mul4(P1[-1], g1))
g2_ = next(g for g in U4 if order4(g) == 2 and g not in P1)
IDX = {}
for i in range(6):
    for j in range(2):
        IDX[P1[i] if j == 0 else mul4(P1[i], g2_)] = (i, j)
assert len(IDX) == 12
CH = [(a, b) for a in range(6) for b in range(2)]
def th(ch, x):
    i, j = IDX[m4(x)]
    return cmath.exp(2j*math.pi*(ch[0]*i/6 + ch[1]*j/2))
inv4 = {c: next(d for d in U4 if mul4(c, d) == (1, 0)) for c in U4}

SQ = primaries(3000, sqfree=True)
Gcls = {}
for n in SQ:
    if n != (1, 0) and m4(n) not in Gcls: Gcls[m4(n)] = G_fast(n)
Gcls[(1, 0)] = 1.0 + 0j          # G(1) = 1 (also the value at squarefree n = 1 mod 4)
assert len(Gcls) == 12
# G is a class function mod 4 on squarefree elements
gerr = max(abs(G_fast(n) - Gcls[m4(n)]) for n in SQ[1:300])
report("T8.6-pre", gerr, 1e-10, "G(n) depends only on n mod 4 (299 squarefree n): the ray group used below")
Ghat = {ch: sum(Gcls[t]*th(ch, t).conjugate() for t in U4)/12 for ch in CH}
Ghat_bad = {ch: sum(Gcls[inv4[t]]*th(ch, t).conjugate() for t in U4)/12 for ch in CH}

# ---- twists nu_0
PI13 = (4, 3); W0 = (-PI13[0]*pow(PI13[1], -1, 13)) % 13
DLOG = {}; x = 1
for k in range(12): DLOG[x] = k; x = (2*x) % 13
def red13(z): return (z[0] + z[1]*W0) % 13
def nu12(n):
    r = red13(n); return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)
NUS = {"1": lambda n: 1.0 + 0j, "chi_.(lambda)": lambda n: chi0(n, (1, 2)), "order 12 mod pi13": nu12}

pool = [n for n in SQ if N(n) <= 300 and red13(n) != 0]
pairs = []
while len(pairs) < 50:
    z1, z2 = random.choice(pool), random.choice(pool)
    if N(gcd_primary(z1, z2)) == 1 and N(z1)*N(z2) <= 8000: pairs.append((z1, z2))
pairs += [((1, 0), pool[3]), (pool[5], (1, 0))]

e = 0.0; f_th = 0.0; f_zeta = 0.0; f_G = 0.0
for name, nu0 in NUS.items():
    for z1, z2 in pairs:
        gs = gauss_psi_general(z1, z2, conj_first=True)            # gamma(conj chi_z1 chi_z2)
        m1 = muk(z1)*phi(z1)*nu0(z1); m2 = muk(z2)*phi(z2)*nu0(z2)
        lhs = m1*m2.conjugate()*gs
        def rhs(cls_power, ntheta, Gh):
            s = 0j
            for ch in CH:
                A1 = phi(z1)**cls_power*g2(z1)*nu0(z1)*(th(ch, z1).conjugate() if ntheta else th(ch, z1))
                A2 = phi(z2)**cls_power*g2(z2)*nu0(z2)*(th(ch, z2).conjugate() if ntheta else th(ch, z2))
                s += Gh[ch]*A1*A2.conjugate()
            return s
        e = max(e, abs(lhs - rhs(2, True, Ghat)))
        f_th = max(f_th, abs(lhs - rhs(2, False, Ghat)))
        f_zeta = max(f_zeta, abs(lhs - rhs(1, True, Ghat)))
        f_G = max(f_G, abs(lhs - rhs(2, True, Ghat_bad)))
report("T8.6", e, 1e-10, f"eq:s4-Astar / eq:s4-a0: m(z1)conj m(z2) gamma(conj chi_z1 chi_z2) = sum_theta Ghat(theta) A_*(z1)conj A_*(z2), "
       f"A_* = conj(alpha)^2 gamma_2 nu_0 conj(theta); 3 twists nu_0 x {len(pairs)} coprime pairs (incl. z = 1)")
falsifier("T8.6", "nu_* = nu_0 theta", f_th, 1e-3)
falsifier("T8.6", "zeta class conj(alpha) gamma_2 nu_*", f_zeta, 1e-3)
falsifier("T8.6", "Ghat from G(t^-1)", f_G, 1e-3)

# ---- extraction identity with A_*
e2 = 0.0; f2 = 0.0; cnt = 0
small = [n for n in pool if N(n) <= 120]
others = [(5, 3), (2, 9), (1, 2), (0, 1), (7, 0), (4, 3), (13, 4), (-1, 6)]
for ch in [(1, 0), (2, 1), (5, 1)]:
    nu0 = NUS["order 12 mod pi13"]
    def Astar(z): return phi(z)**2*g2(z)*nu0(z)*th(ch, z).conjugate()
    for s_, n_ in itertools.combinations(small, 2):
        if N(gcd_primary(s_, n_)) != 1 or N(s_)*N(n_) > 3000: continue
        sn = mul(s_, n_)
        for dp in others[:4]:
            for hp in others[3:]:
                cnt += 1
                lhs = Astar(sn)*chi0(sn, dp).conjugate()*chi0(sn, hp)
                rhs = (Astar(s_)*chi0(s_, dp).conjugate()*chi0(s_, hp)
                       *Astar(n_)*chi0(n_, mul(dp, hp))*chi0(n_, mul(dp, s_))**4)
                bad = (Astar(s_)*chi0(s_, dp).conjugate()*chi0(s_, hp)
                       *Astar(n_)*chi0(n_, mul(dp, hp))*chi0(n_, dp)**4)
                e2 = max(e2, abs(lhs - rhs)); f2 = max(f2, abs(lhs - bad))
report("T8.7", e2, 1e-10, f"eq:s4-Astar extraction (OAI78 eq:initial-common-extraction with A_*): {cnt} cases, 3 ray summands theta")
falsifier("T8.7", "without the CRT factor chi_n(s)^4", f2, 1e-3)
