"""r1 tests of Section 7 (guide lines 945-1180): the finite group frak G of primary invertible classes
mod 72, the functions G, R, chi_x(-1) on it, the quotient identity (7.7) and the paired identities
(7.3)/(7.9)/(7.10) in class form and in literal character-expansion form.

Independent implementation (r1lib.py); run with  python3 -I test_classgroup.py
"""
import os
import sys
import random
import cmath
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *  # noqa

rep = Reporter("test_classgroup")
random.seed(20261012)

# ------------------------------------------------------------------ C1: the group frak G
M = 72
# invertible mod 72 <=> coprime to 2 and to lambda; lambda | a + b w  <=>  a + b = 0 (mod 3) (since w = 1 mod lambda)
units72 = [(a, b) for a in range(M) for b in range(M) if not (a % 2 == 0 and b % 2 == 0) and (a + b) % 3 != 0]
frakG = [x for x in units72 if x[0] % 3 == 1 and x[1] % 3 == 0]


def m72(x):
    return (x[0] % M, x[1] % M)


Gset = set(frakG)
closed = all(m72(mul(x, y)) in Gset for x in frakG[:60] for y in frakG)
fib = {}
for x in frakG:
    fib.setdefault(cls4(x), 0)
    fib[cls4(x)] += 1
rep.check("C1", len(units72) == 2592 and len(frakG) == 432 and closed and len(fib) == 12
          and set(fib.values()) == {36} and set(fib) == set(units_mod4()),
          "|(O/72)^x| = %d, |frak G| = %d (= 2592/6), closed under multiplication: %s; reduction to (O/4)^x "
          "is onto the %d classes, every fibre of size %s" % (len(units72), len(frakG), closed, len(fib),
                                                            sorted(set(fib.values()))))

# every class of frak G contains a squarefree primary n coprime to 6 (so G, xi are defined on frak G)
SQ = primaries(9000, squarefree=True)
hit = {}
for n in SQ:
    c = m72(n)
    if c not in hit:
        hit[c] = n
rep.check("C1b", len(hit) == 432, "squarefree primary n coprime to 6 with N(n) <= 9000 meet %d of the 432 classes"
          % len(hit))

# ------------------------------------------------------------------ C2: G depends only on n mod 4
SQG = [n for n in SQ if N(n) <= 1300]


def G_closed(n):
    a, b = n
    i = 1j
    gq = (1 + i ** ((-b) % 4) + i ** (a % 4) + i ** ((b - a) % 4)) / 2
    return chi(n, (4, 0)).conjugate() * gq


Gval = {}
err_cls = 0.0
err_abs = 0.0
err_closed = 0.0
for n in SQG:
    g = G(n)
    err_abs = max(err_abs, abs(abs(g) - 1))
    err_closed = max(err_closed, abs(g - G_closed(n)))
    c = cls4(n)
    if c in Gval:
        err_cls = max(err_cls, abs(g - Gval[c]))
    else:
        Gval[c] = g
rep.check("C2", err_cls < 1e-9 and err_abs < 1e-9 and err_closed < 1e-9 and len(Gval) == 12,
          "G(n) = conj(chi_n(4)) gamma_3(n) for %d squarefree n (N <= 1300): |G| = 1 (err %.1e), depends only on n mod 4 "
          "(spread %.1e), equals OAI's closed form conj(chi_n(4))(1+i^-b+i^a+i^(b-a))/2 (err %.1e)"
          % (len(SQG), err_abs, err_cls, err_closed))
# falsifier: n mod 2 does not determine G
G2 = {}
spread2 = 0.0
for n in SQG:
    c = (n[0] % 2, n[1] % 2)
    if c in G2:
        spread2 = max(spread2, abs(G(n) - G2[c]))
    else:
        G2[c] = G(n)
rep.falsifier("C2", "G depends only on n mod 2", spread2 > 0.5, "max spread within classes mod 2 = %.2f" % spread2)
# falsifier: gamma_3 alone is NOT G (the factor conj(chi_n(4)) matters)
d3 = max(abs(gamma(3, n) - G(n)) for n in SQG)
rep.falsifier("C2", "drop conj(chi_n(4))", d3 > 0.5, "max |gamma_3 - G| = %.2f" % d3)

Gc = dict(Gval)  # G on the 12 classes mod 4
Gc[(1, 0)] = 1.0 + 0j

# ------------------------------------------------------------------ C3: R
PR = [p for p in primary_primes(3000)]
sample = PR[:120] + random.sample([n for n in SQ if 1 < N(n) <= 3000], 80)
Rcls = {}
err_R = 0.0
bad_pm = 0
pairs = 0
for i in range(len(sample)):
    for j in range(len(sample)):
        a, b = sample[i], sample[j]
        if a == b or N(gcd_primary(a, b)) > 1:
            continue
        if random.random() > 0.15:
            continue
        pairs += 1
        r = chi(b, a) * chi(a, b).conjugate()  # R(a,b) = chi_b(a)/chi_a(b)
        if abs(r - 1) > 1e-9 and abs(r + 1) > 1e-9:
            bad_pm += 1
        key = (cls4(a), cls4(b))
        if key in Rcls:
            err_R = max(err_R, abs(r - Rcls[key]))
        else:
            Rcls[key] = r.real
sym = max(abs(Rcls[(x, y)] - Rcls[(y, x)]) for (x, y) in Rcls if (y, x) in Rcls)
rep.check("C3", bad_pm == 0 and err_R < 1e-9 and len(Rcls) == 144 and sym < 1e-9,
          "R(a,b) = chi_b(a)/chi_a(b) over %d coprime pairs: values +-1 (violations %d), depends only on classes mod 4 "
          "(spread %.1e), all 144 class pairs seen, symmetric (err %.1e)" % (pairs, bad_pm, err_R, len(Rcls) and sym))
cl = units_mod4()
bich = max(abs(Rcls[(mul4(x, xp), y)] - Rcls[(x, y)] * Rcls[(xp, y)]) for x in cl for xp in cl for y in cl)
chim1 = {}
for n in SQG:
    chim1.setdefault(cls4(n), chi(n, (-1, 0)).real)
diag = max(abs(Rcls[(c, c)] - chim1[c]) for c in cl)
Gmult = max(abs(Gc[mul4(x, y)] - Gc[x] * Gc[y] * Rcls[(x, y)]) for x in cl for y in cl)
rep.check("C3b", bich < 1e-9 and diag < 1e-9 and Gmult < 1e-9,
          "on classes mod 4: R is a bicharacter (err %.1e), R(c,c) = chi_c(-1) (err %.1e), G(xy) = G(x)G(y)R(x,y) "
          "(err %.1e) [OAI1112 (4.1) and App. A.1]" % (bich, diag, Gmult))
# falsifier: classes mod 2 do not determine R
R2 = {}
spr = 0.0
for (x, y), v in Rcls.items():
    k2 = ((x[0] % 2, x[1] % 2), (y[0] % 2, y[1] % 2))
    if k2 in R2:
        spr = max(spr, abs(v - R2[k2]))
    else:
        R2[k2] = v
rep.falsifier("C3", "R depends only on classes mod 2", spr > 1, "spread %.0f" % spr)

# ------------------------------------------------------------------ C4: chi_x(-1)
err_m1 = 0
cls_m1 = {}
spread_m1 = 0
for n in primaries(3000):
    v = chi(n, (-1, 0)).real
    if abs(v - (-1) ** ((N(n) - 1) // 6)) > 1e-9:
        err_m1 += 1
    c = cls4(n)
    if c in cls_m1:
        spread_m1 = max(spread_m1, abs(v - cls_m1[c]))
    else:
        cls_m1[c] = v
rep.check("C4", err_m1 == 0 and spread_m1 < 1e-9,
          "chi_x(-1) = (-1)^((N(x)-1)/6) for all primary x coprime to 6, N <= 3000 (violations %d); it depends only on "
          "x mod 4 (spread %.0f)" % (err_m1, spread_m1))
fal = sum(1 for n in primaries(3000) if abs(chi(n, (-1, 0)).real - (-1) ** ((N(n) - 1) // 3)) > 1e-9)
rep.falsifier("C4", "(-1)^((N-1)/3), i.e. chi_x(-1) = 1 (as for the cubic symbol)", fal > 0, "%d violations" % fal)

# ------------------------------------------------------------------ C5: quotient identity (7.7) on frak G
err_q = 0.0
f1 = 0.0
f2 = 0.0
for x in cl:
    for y in cl:
        t = mul4(y, inv4(x))
        lhs = chim1[x] * Gc[x].conjugate() * Gc[y] * Rcls[(x, y)]
        err_q = max(err_q, abs(lhs - Gc[t]))
        f1 = max(f1, abs(Gc[x].conjugate() * Gc[y] * Rcls[(x, y)] - Gc[t]))
        f2 = max(f2, abs(lhs - Gc[mul4(x, inv4(y))]))
rep.check("C5", err_q < 1e-9,
          "(7.7) chi_x(-1) conj(G(x)) G(y) R(x,y) = G(y x^-1) for all 144 pairs of classes mod 4, hence on frak G x frak G "
          "(err %.1e)" % err_q)
rep.falsifier("C5", "drop chi_x(-1)", f1 > 0.5, "err %.2f" % f1)
rep.falsifier("C5", "G(x y^-1) on the right", f2 > 0.5, "err %.2f" % f2)

# ------------------------------------------------------------------ C6: paired identities, class form
xi_omega = lambda n: chi(n, (0, 1))  # chi_n(omega): a class function mod 9, a character of frak G
c0 = (1, 2)
xi_R = lambda n: Rcls[(cls4(n), c0)]  # R(., c0): a character of (O/4)^x (bicharacter), pulled back to frak G
SQP = [n for n in SQ if N(n) <= 400]
cnt = 0
e_f = e_r = e_x = 0.0
fz = fG = fr = 0.0
while cnt < 70:
    z1, z2 = random.choice(SQP), random.choice(SQP)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 2600:
        continue
    cnt += 1
    gf = gauss_char(z1, z2, True)   # gamma(conj chi_z1 chi_z2)
    gr = gauss_char(z1, z2, False)  # gamma(chi_z1 conj chi_z2)
    A = lambda n: alpha(n).conjugate() ** 2 * gamma(2, n)      # CM coefficient conj(alpha)^2 gamma_2
    Az = lambda n: alpha(n).conjugate() * gamma(2, n)          # zeta coefficient
    t21 = mul4(z2, inv4(z1))  # class of z2 z1^-1
    t12 = mul4(z1, inv4(z2))
    lhs = mu(z1) * mu(z2) * alpha(z1).conjugate() * alpha(z2) * gf
    e_f = max(e_f, abs(lhs - Gc[t21] * A(z1) * A(z2).conjugate()))
    fz = max(fz, abs(lhs - Gc[t21] * Az(z1) * Az(z2).conjugate()))
    fG = max(fG, abs(lhs - Gc[t12] * A(z1) * A(z2).conjugate()))
    for xi1 in (xi_omega, xi_R):
        e_x = max(e_x, abs(lhs * xi1(z1) * xi1(z2).conjugate()
                           - xi1(z1) * xi1(z2).conjugate() * Gc[t21] * A(z1) * A(z2).conjugate()))
    lr = A(z1) * A(z2).conjugate() * gr
    rr = mu(z1) * mu(z2) * alpha(z1).conjugate() * alpha(z2) * Gc[t12]
    e_r = max(e_r, abs(lr - rr))
    fr = max(fr, abs(lr - mu(z1) * mu(z2) * Gc[t12]))
rep.check("C6", max(e_f, e_r, e_x) < 1e-9,
          "class form over %d random coprime squarefree pairs (N(z1 z2) <= 2600): forward (7.3)/(7.9) err %.1e, "
          "with xi_1 = chi_n(omega) or a character mod 4 err %.1e, reverse (7.10) err %.1e" % (cnt, e_f, e_x, e_r))
rep.falsifier("C6", "zeta coefficient conj(alpha) gamma_2 in the forward identity", fz > 0.1, "err %.2f" % fz)
rep.falsifier("C6", "G(t) for G(t^-1) in the forward identity", fG > 0.1, "err %.2f" % fG)
rep.falsifier("C6", "reverse identity without conj(alpha(z1)) alpha(z2) (zeta form)", fr > 0.1, "err %.2f" % fr)

# ------------------------------------------------------------------ C7: literal character expansion
# characters of (O/4)^x (order 12) found by brute force; characters of frak G that do not factor through
# (O/4)^x are orthogonal to t -> xi_1(t) G(t^-1) when xi_1 factors through mod 4, so their c_xi vanish.
def order4(x):
    k, y = 1, x
    while y != (1, 0):
        y = mul4(y, x)
        k += 1
    return k


# (O/4)^x = Z/3 x (Z/2)^2: an element g3 of order 3 and two elements h1, h2 of order 2 generate it
g3 = next(g for g in cl if order4(g) == 3)
inv2 = [g for g in cl if order4(g) == 2]
h1 = inv2[0]
h2 = next(g for g in inv2 if g != h1)
words = {}
for i in range(3):
    for j in range(2):
        for k in range(2):
            z = (1, 0)
            for _ in range(i):
                z = mul4(z, g3)
            if j:
                z = mul4(z, h1)
            if k:
                z = mul4(z, h2)
            words[z] = (i, j, k)
assert len(words) == 12
chars = []
for a in range(3):
    for b in range(2):
        for c in range(2):
            chars.append({z: cmath.exp(2j * math.pi * a * i / 3) * (-1) ** (b * j) * (-1) ** (c * k)
                          for z, (i, j, k) in words.items()})
assert all(abs(ch[mul4(x, y)] - ch[x] * ch[y]) < 1e-12 for ch in chars for x in cl for y in cl)
cxi = [sum(Gc[inv4(t)] * ch[t].conjugate() for t in cl) / 12 for ch in chars]
lit = 0.0
cnt = 0
while cnt < 40:
    z1, z2 = random.choice(SQP), random.choice(SQP)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 2600:
        continue
    cnt += 1
    A = lambda n: alpha(n).conjugate() ** 2 * gamma(2, n)
    lhs = mu(z1) * mu(z2) * alpha(z1).conjugate() * alpha(z2) * gauss_char(z1, z2, True)
    rhs = sum(c * (A(z1) * ch[cls4(z1)]) * (A(z2) * ch[cls4(z2)]).conjugate() for c, ch in zip(cxi, chars))
    lit = max(lit, abs(lhs - rhs))
rep.check("C7", len(chars) == 12 and max(abs(c) for c in cxi) <= 1 + 1e-12 and abs(sum(cxi) - 1) < 1e-12 and lit < 1e-9,
          "the 12 characters of (O/4)^x: |c_xi| <= 1 (max %.3f), sum c_xi = G(1) = 1 (err %.1e); literal expansion "
          "sum_xi c_xi a_xi(z1) conj(a_xi(z2)) reproduces the left side of (7.3) for %d pairs (err %.1e)"
          % (max(abs(c) for c in cxi), abs(sum(cxi) - 1), cnt, lit))

sys.exit(1 if rep.summary() else 0)
