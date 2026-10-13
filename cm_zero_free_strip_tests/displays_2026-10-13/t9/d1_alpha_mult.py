"""d1: (F.1) eq:s5-alpha-mult and the allocation scalar used in (F.23) eq:s5-coeff / Lemma s5-coefficients.

F.1  alpha(l l') = alpha(l) alpha(l'), |alpha(l)| = 1 for ideals l, l' outside S (S contains the primes above 6),
     where alpha(l) = l0/|l0| with l0 the primary generator of l.
     Test (exact part): for every pair of primary l0, l0' coprime to 6 with N <= 400, the product l0 l0' (computed
     exactly in Z[omega]) is primary and coprime to 6; the primary generator of the ideal (l l') is found
     independently from a randomly rotated generator u l0 l0'; alpha of that generator equals alpha(l) alpha(l')
     to rounding, and |alpha| = 1.
F.23 consistency: allocating p^{i1} to l1 and p^{i2} to l2 multiplies alpha(l1 l2)^e by alpha(p)^{e(i1+i2)}
     (the frozen unimodular scalar of Lemma s5-coefficients (iii)-(v)), checked with the primary generator of the
     product ideal computed independently.
Falsifiers: Ireland-Rosen's normalization (generator = 2 mod 3) and the first-sextant normalization (argument in
[0, pi/3)) are not multiplicative.
"""
import os, sys, random, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t9lib import (report, falsifier, summary, cx, mul, conj, N, UNITS, primary_assoc, is_primary, divides,
                   LAM, TWO, alpha_el, sextant_generator)

rnd = random.Random(20261013)

def coprime6(a):
    return not divides(LAM, a) and not divides(TWO, a)

prim = []
R = 26
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        x = (a, b)
        if 0 < N(x) <= 400 and is_primary(x) and coprime6(x):
            prim.append(x)
print(f"{len(prim)} primary elements coprime to 6 with N <= 400 (one per ideal outside S = {{lambda, 2}})")

def alpha_ideal(gen):
    """alpha of the ideal (gen): alpha of its primary generator, found from an arbitrary generator."""
    u, g0 = primary_assoc(gen)
    return alpha_el(g0)

def ir_gen(gen):         # Ireland-Rosen: generator = 2 mod 3, i.e. minus the primary one
    u, g0 = primary_assoc(gen); return (-g0[0], -g0[1])

nonprimary = 0; not_coprime = 0; err = 0; errmod = 0; f_ir = 0; f_sext = 0; npairs = 0
for x in prim:
    errmod = max(errmod, abs(abs(alpha_ideal(x)) - 1))
    for y in prim:
        p = mul(x, y)
        npairs += 1
        if not is_primary(p): nonprimary += 1
        if not coprime6(p): not_coprime += 1
        rot = mul(p, UNITS[rnd.randrange(6)])            # arbitrary generator of the product ideal
        lhs = alpha_ideal(rot)
        err = max(err, abs(lhs - alpha_el(x) * alpha_el(y)))
        f_ir = max(f_ir, abs(alpha_el(ir_gen(rot)) - alpha_el(ir_gen(x)) * alpha_el(ir_gen(y))))
        f_sext = max(f_sext, abs(alpha_el(sextant_generator(rot)) - alpha_el(sextant_generator(x)) * alpha_el(sextant_generator(y))))
report("F1-primary-closed", nonprimary + not_coprime, 0,
       f"{npairs} products l0 l0' are primary and coprime to 6 (exact Z[omega] arithmetic)")
report("F1-mult", err, 1e-14, "alpha(l l') = alpha(l) alpha(l') with the primary generator of l l' found independently")
report("F1-modulus", errmod, 1e-15, "|alpha(l)| = 1")
falsifier("F1-mult", "Ireland-Rosen generator (= 2 mod 3)", f_ir, 0.5)
falsifier("F1-mult", "first-sextant generator", f_sext, 0.1)

# F.23 consistency: allocation scalar
err = 0; fals = 0
primes = [x for x in prim if N(x) in (7, 13, 19, 25, 31, 37)]
for _ in range(400):
    p = rnd.choice(primes); l1 = rnd.choice(prim); l2 = rnd.choice(prim)
    i1 = rnd.randrange(0, 4); i2 = rnd.randrange(0, 4); e = rnd.choice((1, -1))
    n1 = l1; n2 = l2
    for _ in range(i1): n1 = mul(n1, p)
    for _ in range(i2): n2 = mul(n2, p)
    rot = mul(mul(n1, n2), UNITS[rnd.randrange(6)])
    lhs = alpha_ideal(rot) ** e
    rhs = alpha_el(p) ** (e * (i1 + i2)) * alpha_ideal(mul(l1, l2)) ** e
    err = max(err, abs(lhs - rhs))
    if i1 + i2 > 0 and abs(alpha_el(p) ** e - 1) > 1e-3:
        fals = max(fals, abs(lhs - alpha_ideal(mul(l1, l2)) ** e))
report("F23-allocation", err, 1e-13, "alpha(p^{i1} l1 p^{i2} l2)^e = alpha(p)^{e(i1+i2)} alpha(l1 l2)^e (frozen unimodular scalar)")
falsifier("F23-allocation", "dropping the extracted scalar alpha(p)^{e(i1+i2)}", fals, 1e-3)
summary()
