"""s6-T3, s6-T4, s6-T5.

s6-T3: exact (Fraction) check of the exponent bookkeeping of OAI78 Part I (paper.tex L5905-6012, L6465-6510,
L6582-6698), which Section s6.3 of the insert claims is unchanged for F_CM.  Includes: C_I(s) = s - 2/3 from
eq:stage-mellin-exponent, the envelope R(delta), the three pieces of E_I(1/2), the comparison
E_I(1/2) - Delta_1 <= -1/24 - Delta_1/2 under the CM bin ceiling delta <= 2 beta_* - 1 = 5/6 + 2 Delta_1, the slope
range [33/50, 62/75], m_c = 1021/25000, the small-row exponent -43/300, the large-row exponent at v = 4000, the loss
budget 46 eps_* < 1/4800, the consistency of E_I with E_{sigma0}(d;R,0) of eq:stage-high-exponent, and the Part II
bootstrap Delta <= 1/24, kappa <= 5/6 obtained from beta_* <= 11/12.
Falsifier: without the bin ceiling (delta = 1) the third piece gives 0, so the comparison fails for small Delta_1.

s6-T4: unit rows.  For each unit u != 1 the numerator character chi_.(u) is nonprincipal on primes outside S
(OAI78 L5502-5514): a primary prime p coprime to 6 with chi_p(u) != 1 is exhibited; and the six units give six
distinct characters.  Falsifier: u = 1 gives the principal character.

s6-T5: the twisted family has the CM Euler product: for primary primes p outside S = {2,3,pi13},
alpha(p) nu(p) = conj(psi((p))) with psi = conj(alpha) eta, eta = conj(nu); psi is not of finite order (its values at
the first primes are not roots of unity of bounded order), and psi((a)) = conj(alpha(a)) for a = 1 mod 3 pi13
(infinity type z-bar/|z|).  Falsifier: psi((a)) = alpha(a).
Usage: python3 s6_margins_and_rows.py [path to tests_cm/part_B_lines_0550_0997]
"""
import os, sys, math, cmath
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')
sys.path.insert(0, os.path.abspath(LIB))
from cmlib import *

# ------------------------------------------------------------------ s6-T3
fails = []
def chk(name, cond):
    if not cond: fails.append(name)
lx = ly = h = F(1, 2); ell = F(0); z0 = F(17, 50)
C = lambda s: s + lx/2 - 1 + h/6
chk("C_I(s) = s - 2/3", C(F(0)) == F(-2, 3) and C(F(11, 12)) == F(1, 4))
lx2, h2, ell2 = F(17, 48), F(13, 16), F(1, 6)
chk("h = 1 - l_x + ell (Part II)", h2 == 1 - lx2 + ell2)
chk("C_II(s) = s - 11/16", F(0) + lx2/2 - 1 + h2/6 == F(-11, 16))
def R(d): return min(F(1), max(1 - d/2, F(4, 3) - d))
def EI(d, delta): return F(-3, 4) + (delta + R(delta))/2 + (d - F(1, 2))*(R(delta) + delta/2 - z0)
grid = [F(k, 600) for k in range(0, 601)]
for dl in grid:
    e = EI(F(1, 2), dl)
    piece = (-F(1, 4) + dl/2) if dl <= F(1, 3) else (F(-1, 12) if dl <= F(2, 3) else -(1 - dl)/4)
    chk(f"E_I(1/2) piece at {dl}", e == piece)
    slope = R(dl) + dl/2 - z0
    chk(f"slope range at {dl}", F(33, 50) <= slope <= F(62, 75))
worst = F(-10)
for k in range(1, 201):
    D1 = F(k, 4800)                               # 0 < Delta_1 <= 1/12 (beta_* <= 1)
    for dl in grid:
        if dl > F(5, 6) + 2*D1: continue          # CM bin ceiling delta <= 2 beta_* - 1
        worst = max(worst, EI(F(1, 2), dl) - D1 + F(1, 24) + D1/2)
chk("E_I(h) - Delta_1 <= -1/24 - Delta_1/2 under the bin ceiling", worst <= 0)
D1 = F(1, 4800)
chk("FALSIFIER: no bin ceiling (delta=1) breaks the comparison", not (EI(F(1, 2), F(1)) - D1 <= -F(1, 24) - D1/2))
mc = F(1, 24) - F(62, 75)*F(1, 1000)
chk("m_c = 1021/25000", mc == F(1021, 25000))
dmin = F(1, 63)
chk("row power 1 + 3/5 - z0 = 63/50", 1 + F(3, 5) - z0 == F(63, 50))
chk("small rows -43/300", h*(z0 - F(1, 6)) - ly/2 + F(63, 50)*dmin == F(-43, 300))
eps = F(1, 10**6); v = 4000
B0 = lx/2 + 1 + ly
chk("B_0 = 7/4", B0 == F(7, 4))
chk("large rows exponent < -1 at v=4000", B0 + F(501, 1000)*(1 + eps) - F(v, 1000) + eps < -1)
chk("margins m_w = 1/40, m_z = 1/1200", ly/20 == F(1, 40) and h/600 == F(1, 1200))
chk("46 eps_* < 1/4800 = m_0/4", 46*eps < F(1, 4800) and F(1, 1200)/4 == F(1, 4800))
chk("principal losses < m_0/4", (F(3, 2) + 1)*eps < F(1, 4800) and 2*eps < F(1, 4800))
# consistency of E_I(d) with the general exponent E_{sigma0}(d;R,g) of eq:stage-high-exponent (sigma0 = 11/12, g = 0)
for dl in grid[::10]:
    a = (1 + dl)/2
    for d in (dmin, F(1, 4), F(1, 2), F(501, 1000)):
        Eg = a - F(11, 12) + h*(z0 - F(1, 6)) - a*ly - ell/2 + 0 + d*(R(dl) + dl/2 - z0)
        chk(f"E_sigma0 = E_I at delta={dl}, d={d}", Eg == EI(d, dl))
# Part II bootstrap from beta_* <= 11/12
for k in range(1, 101):
    beta = F(7, 8) + F(k, 100)*(F(11, 12) - F(7, 8))
    Delta = beta - F(7, 8); kappa = 2*beta - 1
    chk("Delta <= 1/24 and kappa = 3/4 + 2 Delta <= 5/6", 0 < Delta <= F(1, 24) and kappa == F(3, 4) + 2*Delta and kappa <= F(5, 6))
report("s6-T3", float(len(fails)), 0, f"exponent bookkeeping of OAI78 Part I (exact fractions); failed: {fails[:3]}")
falsifier("s6-T3", "no bin ceiling", 1.0 if (EI(F(1, 2), F(1)) - D1 > -F(1, 24) - D1/2) else 0.0, 0.5)

# ------------------------------------------------------------------ s6-T4
UN = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]     # 1, -1, omega, -omega, omega^2, -omega^2
assert all(N(u) == 1 for u in UN)
primes = [p for p in PR_ODD if N(p) < 2000]
sigs = {}
bad = 0
for u in UN:
    sig = tuple(int(round(cmath.phase(chi(p, u))/(math.pi/3))) % 6 for p in primes[:40])
    sigs[u] = sig
    nontriv = any(s != 0 for s in sig)
    if u != (1, 0) and not nontriv: bad += 1
distinct = len(set(sigs.values())) == 6
report("s6-T4", float(bad + (0 if distinct else 1)), 0, "each unit u != 1 has chi_p(u) != 1 at some of the first 40 primary primes; "
       "the six unit characters are distinct")
falsifier("s6-T4", "u = 1 is principal (no p with chi_p(1) != 1)", 0.0 if any(s != 0 for s in sigs[(1, 0)]) else 1.0, 0.5)

# ------------------------------------------------------------------ s6-T5
PI = (4, 3); P13 = 13
W0 = (-PI[0]*pow(PI[1], -1, P13)) % P13
DLOG = {}; x = 1
for k in range(12): DLOG[x] = k; x = (2*x) % P13
def red13(x): return (x[0] + x[1]*W0) % P13
def nu(n):
    r = red13(n)
    return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)
def psi(n):                       # psi((n)) = conj(alpha(n)) eta(n), eta = conj(nu), n primary coprime to S
    return alpha_c(n).conjugate()*nu(n).conjugate()
err = 0.0
for p in primes[:300]:
    if red13(p) == 0: continue
    err = max(err, abs(alpha_c(p)*nu(p) - psi(p).conjugate()))
# psi on a = 1 mod 3 pi13 equals conj(alpha(a)); falsifier alpha(a)
e2 = 0.0; f2 = 0.0; cnt = 0
for a in range(-60, 61):
    for b in range(-60, 61):
        n = (a, b)
        if N(n) <= 1 or not primary(n) or red13(n) != 1: continue
        cnt += 1
        e2 = max(e2, abs(psi(n) - alpha_c(n).conjugate()))
        f2 = max(f2, abs(psi(n) - alpha_c(n)))
# not of finite order: the arguments of psi(p) at primes are not multiples of 2 pi / M for M <= 1000
ang = [cmath.phase(psi(p))/(2*math.pi) for p in primes[:50] if red13(p) != 0]
fin = any(all(abs(t*M - round(t*M)) < 1e-9 for t in ang) for M in range(1, 1001))
report("s6-T5", max(err, e2) + (1.0 if fin else 0.0), 1e-12,
       f"alpha(p)nu(p) = conj psi((p)) at 300 primes; psi((a)) = conj alpha(a) for {cnt} elements a = 1 mod 3 pi13; "
       "psi not of order <= 1000")
falsifier("s6-T5", "infinity type alpha instead of conj(alpha)", f2, 1e-3)
