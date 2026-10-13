"""s6-T2: the CM paired identities (Lemma lem:b2-paired of the guide) over an ENLARGED class group G_S,
with a finite-order twist nu whose conductor is a prime outside 6 and whose order (12) does not divide 6.

Setting (Appendix app:s6 of the insert):
  pi13 = 4 + 3 omega (primary, norm 13), S = {primes over 2, 3} + {pi13}.
  nu(n) = exp(2 pi i dlog_2(n mod pi13)/12) on primary n coprime to S (order 12), 0 otherwise.
  G_S = (O/4)^x x (O/pi13)^x (= the classes of primary elements coprime to S modulo 4 pi13), order 144.
  G(t), R, chi_t(-1) depend on t mod 4; nu on t mod pi13.
  c_xi(xi1) = |G_S|^{-1} sum_t xi1(t) G(t^{-1}) conj(xi(t)),   (xi G)(t) = sum_{xi1} l_{xi1}(xi) xi1(t).

Literal checks (sum over all 144 characters xi):
  forward  (eq:b2-forward, xi1 = conj nu):
     mu(z1)mu(z2) conj(alpha(z1)) alpha(z2) conj(nu(z1)) nu(z2) gamma(conj chi_z1 chi_z2)
        = sum_xi c_xi(conj nu) a_xi(z1) conj(a_xi(z2)),      a_xi = conj(alpha)^2 gamma_2 xi
  reverse  (eq:b2-reverse, for xi = nu and for a random xi):
     a_xi(z1) conj(a_xi(z2)) gamma(chi_z1 conj chi_z2) = sum_{xi1} l_{xi1}(xi) beta_{xi1}(z1) conj(beta_{xi1}(z2)),
     beta_{xi1} = mu conj(alpha) xi1
Falsifiers (must fail): xi1 = nu (wrong orientation); c_xi built with G(t) instead of G(t^{-1});
  untwisted c_xi (xi1 = 1); reverse with beta = mu alpha xi1.
Usage: python3 s6_twisted_paired.py [path to tests_cm/part_B_lines_0550_0997]
"""
import os, sys, math, cmath, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')
sys.path.insert(0, os.path.abspath(LIB))
from cmlib import *          # zw arithmetic, chi, mu, factor, primaries, gamma_j_fast, gauss_psi_general, G_fast, alpha_c
random.seed(20261012)

PI = (4, 3); P13 = 13
assert primary(PI) and N(PI) == P13
W0 = (-PI[0]*pow(PI[1], -1, P13)) % P13          # image of omega in F_13
assert (W0*W0 + W0 + 1) % P13 == 0
DLOG = {}
x = 1
for k in range(12):
    DLOG[x] = k; x = (2*x) % P13
assert len(DLOG) == 12                             # 2 is a primitive root mod 13

def red13(x): return (x[0] + x[1]*W0) % P13
def nu(n):
    r = red13(n)
    return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)

# ---- (O/4)^x: structure C6 x C2
def m4(x): return (x[0] % 4, x[1] % 4)
U4 = [(a, b) for a in range(4) for b in range(4) if not (a % 2 == 0 and b % 2 == 0)]
assert len(U4) == 12
def mul4(x, y): return m4(mul(x, y))
def order4(g):
    k, y = 1, g
    while y != (1, 0): y = mul4(y, g); k += 1
    return k
g1 = next(g for g in U4 if order4(g) == 6)
pows1 = [(1, 0)]
for _ in range(5): pows1.append(mul4(pows1[-1], g1))
g2 = next(g for g in U4 if order4(g) == 2 and g not in pows1)
IDX4 = {}
for i in range(6):
    for j in range(2):
        e = pows1[i] if j == 0 else mul4(pows1[i], g2)
        IDX4[e] = (i, j)
assert len(IDX4) == 12                             # unique representation g1^i g2^j

# ---- G on classes mod 4 (G depends only on n mod 4)
sqall = primaries(3000, sqfree=True)
Gcls = {}
for n in sqall:
    c = m4(n)
    if c not in Gcls: Gcls[c] = 1.0 + 0j if n == (1, 0) else G_fast(n)
assert len(Gcls) == 12
inv4 = {c: next(d for d in U4 if mul4(c, d) == (1, 0)) for c in U4}

# ---- the group G_S = (O/4)^x x (O/pi)^x, elements (c4, r13) with r13 in 1..12
GS = [(c4, r) for c4 in U4 for r in range(1, 13)]
assert len(GS) == 144
CHARS = [(j1, j2, j3) for j1 in range(6) for j2 in range(2) for j3 in range(12)]
def xi_val(ch, t):
    (i1, i2) = IDX4[t[0]]; i3 = DLOG[t[1]]
    j1, j2, j3 = ch
    return cmath.exp(2j*math.pi*(j1*i1/6 + j2*i2/2 + j3*i3/12))
def cls(n): return (m4(n), red13(n))
NU_CH = (0, 0, 1)                                     # nu as a character of G_S
assert all(abs(xi_val(NU_CH, cls(n)) - nu(n)) < 1e-12 for n in sqall[:400] if nu(n) != 0)
def tinv(t): return (inv4[t[0]], pow(t[1], -1, P13))
def Gt(t): return Gcls[t[0]]

def cvec(xi1_fn, Gfun):
    """c_xi(xi1) for all xi: |G_S|^{-1} sum_t xi1(t) G(.) conj(xi(t))."""
    return {ch: sum(xi1_fn(t)*Gfun(t)*xi_val(ch, t).conjugate() for t in GS)/len(GS) for ch in CHARS}
c_twist = cvec(lambda t: xi_val(NU_CH, t).conjugate(), lambda t: Gt(tinv(t)))      # xi1 = conj nu
c_wrong = cvec(lambda t: xi_val(NU_CH, t), lambda t: Gt(tinv(t)))                  # xi1 = nu (falsifier)
c_Gswap = cvec(lambda t: xi_val(NU_CH, t).conjugate(), lambda t: Gt(t))            # G(t) (falsifier)
c_plain = cvec(lambda t: 1.0, lambda t: Gt(tinv(t)))                               # xi1 = 1 (falsifier)
assert max(abs(v) for v in c_twist.values()) <= 1 + 1e-12

g2c = {}
def g2(n):
    if n not in g2c: g2c[n] = 1.0 + 0j if n == (1, 0) else gamma_j_fast(2, n)
    return g2c[n]
def muk(n): return 1 if n == (1, 0) else mu(n)
def a_xi(ch, n): return alpha_c(n).conjugate()**2*g2(n)*xi_val(ch, cls(n))

sq = [n for n in sqall if N(n) <= 400 and nu(n) != 0]
pairs = []
while len(pairs) < 60:
    z1, z2 = random.choice(sq), random.choice(sq)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 6000: continue
    pairs.append((z1, z2))
e_f = f_wrong = f_G = f_plain = 0.0
for z1, z2 in pairs:
    gf = gauss_psi_general(z1, z2, conj_first=True)
    lhs = muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*nu(z1).conjugate()*nu(z2)*gf
    A1 = {ch: a_xi(ch, z1) for ch in CHARS}; A2 = {ch: a_xi(ch, z2) for ch in CHARS}
    def rhs(c): return sum(c[ch]*A1[ch]*A2[ch].conjugate() for ch in CHARS)
    e_f = max(e_f, abs(lhs - rhs(c_twist)))
    f_wrong = max(f_wrong, abs(lhs - rhs(c_wrong)))
    f_G = max(f_G, abs(lhs - rhs(c_Gswap)))
    f_plain = max(f_plain, abs(lhs - rhs(c_plain)))
ok = report("s6-T2a", e_f, 1e-10, f"forward paired identity over G_S (144 chars), xi1 = conj nu, nu of order 12 mod pi13; {len(pairs)} pairs")
falsifier("s6-T2a", "xi1 = nu (orientation)", f_wrong, 1e-3)
falsifier("s6-T2a", "G(t) in place of G(t^-1)", f_G, 1e-3)
falsifier("s6-T2a", "untwisted c_xi (xi1 = 1)", f_plain, 1e-3)

# ---- reverse identity: (xi G)(t) = sum_{xi1} l_{xi1}(xi) xi1(t)
def lvec(ch):
    return {c1: sum(xi_val(ch, t)*Gt(t)*xi_val(c1, t).conjugate() for t in GS)/len(GS) for c1 in CHARS}
e_r = f_r = 0.0
for ch in (NU_CH, random.choice(CHARS)):
    L = lvec(ch)
    for z1, z2 in pairs[:30]:
        gr = gauss_psi_general(z1, z2, conj_first=False)
        lhs = a_xi(ch, z1)*a_xi(ch, z2).conjugate()*gr
        rhs = sum(L[c1]*muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*xi_val(c1, cls(z1))*xi_val(c1, cls(z2)).conjugate()
                  for c1 in CHARS)
        bad = sum(L[c1]*muk(z1)*muk(z2)*alpha_c(z1)*alpha_c(z2).conjugate()*xi_val(c1, cls(z1))*xi_val(c1, cls(z2)).conjugate()
                  for c1 in CHARS)
        e_r = max(e_r, abs(lhs - rhs)); f_r = max(f_r, abs(lhs - bad))
report("s6-T2b", e_r, 1e-10, "reverse paired identity over G_S, xi = nu and a random xi, 30 pairs each")
falsifier("s6-T2b", "beta = mu alpha xi1 (not conj alpha)", f_r, 1e-3)
