"""c1: the presentation eq:s2-pres (C.1), the Gauss sum gamma(phi) of eq:d147 (C.2) and the twisted Gauss sum identity
eq:s2-gauss (C.3) of Appendix C of the CM guide.

Characters: the eight primitive phi of t7lib.test_chars (phi(zeta) = zeta; Q = 7 ... 225; lambda-exponent 0, 1, 2;
primes over 2, inert 5, split 7, 13).

Claims
  T7-C1a eq:s2-pres  psi((a)) := phi(a) conj(alpha(a)) does not depend on the generator a of (a) (all six associates),
         is multiplicative on ideals, and has infinity type zbar/|z|: psi((b))/psi((c)) = conj(alpha(b/c)) whenever
         b = c mod f (b, c coprime to f).  All a with N(a) <= 400; 3000 random pairs b = c mod f.
         Falsifier: phi(a) alpha(a) (infinity type z/|z|, wrong for phi(zeta) = zeta) depends on the generator.
  T7-C1b Lemma s2-pres(iii)  psi = alphabar: phi(a) = the unit congruent to a mod 3, and psi((n)) = conj(alpha(n)) on
         primary n (N(n) <= 400).
  T7-C2  eq:d147  gamma(phi) = Q^{-1/2} sum_{x mod f} phi(x) e(x/f): example gamma(alphabar) = -1, i.e.
         sum_zeta zeta e(zeta/3) = -3 exactly; W(alphabar) = +1; W = -gamma sqrt(Q)/f does not depend on the generator
         of (f) (all six), since gamma changes by phi(zeta) = zeta.
  T7-C3  eq:s2-gauss  sum_{x mod f} phi(x) e(hx/f) = conj(phi(h)) sqrt(Q) gamma(phi) for EVERY h mod f (also h not
         coprime to f, where both sides vanish by primitivity).  Checked EXACTLY in Z[zeta_M], M = lcm(6, Q) (the values
         of phi are sixth roots of unity and e(hx/f) = exp(2 pi i b/Q) with b the omega-coordinate of h x conj(f)),
         and numerically; also |gamma| = 1 exactly (gamma conj(gamma) = 1 in Z[zeta_M]).
         Falsifiers: phi(h) in place of conj(phi(h)) fails; an imprimitive character (u3 viewed mod 3 pi7) fails at
         some h not coprime to f.
"""
import os, sys, math, cmath, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t7lib import *

t0 = time.time()
rng = random.Random(20261013)
chars = test_chars()

# ---------------------------------------------------------------- T7-C1a
pts = lattice(400)
err_gen = err_mult = err_inf = 0.0
fals_gen = 0.0
for ch in chars:
    for a in pts:
        if not coprime(a, ch.f):
            continue
        v = ch.psi_elem(a)
        for u in UNITS:
            err_gen = max(err_gen, abs(ch.psi_elem(mul(u, a)) - v))
            # falsifier: phi(a) alpha(a)
            za, zua = cx(a), cx(mul(u, a))
            fals_gen = max(fals_gen, abs(ch.phi(mul(u, a)) * zua / abs(zua) - ch.phi(a) * za / abs(za)))
    cop = [a for a in pts if coprime(a, ch.f)]
    for _ in range(1500):
        a, b = rng.choice(cop), rng.choice(cop)
        err_mult = max(err_mult, abs(ch.psi_elem(mul(a, b)) - ch.psi_elem(a) * ch.psi_elem(b)))
    for _ in range(3000):
        b = rng.choice(cop)
        y = (rng.randint(-6, 6), rng.randint(-6, 6))
        c = add(b, mul(ch.f, y))
        if N(c) == 0:
            continue
        lhs = ch.psi_elem(b) / ch.psi_elem(c)
        q = cx(b) / cx(c)
        err_inf = max(err_inf, abs(lhs - (q / abs(q)).conjugate()))
report("T7-C1a", "eq:s2-pres", "psi((a)) = phi(a) conj(alpha(a)) is generator-independent, multiplicative, of type zbar/|z|",
       max(err_gen, err_mult, err_inf) < 1e-12,
       f"8 characters; max errors: generators {err_gen:.1e}, multiplicativity {err_mult:.1e}, "
       f"psi((b))/psi((c)) = conj alpha(b/c) for b = c mod f: {err_inf:.1e}; falsifier phi(a)alpha(a): generator "
       f"dependence {fals_gen:.2f}", falsifier_failed=fals_gen > 0.5)

# ---------------------------------------------------------------- T7-C1b
A = chars[0]
e1 = e2 = 0.0
for a in pts:
    if divides(LAM, a):
        continue
    j = A.idx(a)
    e2 = max(e2, 0.0 if divides((3, 0), sub(a, UNITS[j])) else 1.0)
    if primary(a):
        z = cx(a)
        e1 = max(e1, abs(A.psi_elem(a) - (z / abs(z)).conjugate()))
report("T7-C1b", "lem:s2-pres(iii)", "alphabar: phi(a) = unit = a mod 3; psi((n)) = conj(alpha(n)) on primary n",
       e1 < 1e-12 and e2 == 0, f"max error {e1:.1e}; unit mismatches {int(e2)}")

# ---------------------------------------------------------------- T7-C2
gA = sum(cx(u) * cmath.exp(4j * math.pi * (cx(u) / 3).imag / SQ3) for u in UNITS)
# exact: e(zeta/3) = exp(2 pi i b/9) with b the omega-coordinate of zeta*conj(3) = 3 zeta, i.e. e(zeta/3) = exp(2 pi i b_zeta/3)
exact_terms = [(j, (3 * UNITS[j][1]) % 9) for j in range(6)]
M = 18
phiM = cyclo(M)
vec = [0] * M
for j, b in exact_terms:
    vec[(j * 3 + b * 2) % M] += 1           # zeta6^j = zeta18^{3j}, exp(2 pi i b/9) = zeta18^{2b}
vec[0] += 3                                  # + 3
exact_ok = all(c == 0 for c in reduce_cyclo(vec, phiM))
werr = 0.0
for ch in chars:
    W0 = ch.root_number()
    for u in UNITS:
        g2 = Char(ch.name, mul(u, ch.f), ch.k3, ch.comps)
        werr = max(werr, abs(g2.root_number() - W0))
report("T7-C2", "eq:d147", "gamma(alphabar) = -1 (sum zeta e(zeta/3) = -3 exactly), W(alphabar) = 1; W independent of generator f",
       exact_ok and abs(gA + 3) < 1e-12 and abs(A.root_number() - 1) < 1e-12 and werr < 1e-12,
       f"sum = {gA:.12f}; exact in Z[zeta_18]: {exact_ok}; W(alphabar) = {A.root_number():.12f}; "
       f"max |W(f) - W(zeta f)| over 8 characters and 6 units: {werr:.1e}")


# ---------------------------------------------------------------- T7-C3
def gauss_check(ch, exact=True):
    """returns (numeric max error, exact ok, falsifier phi(h) max error)."""
    Q = ch.Q
    M = Q * 6 // math.gcd(Q, 6)
    phiM = cyclo(M) if exact else None
    fc = conj(ch.f)
    idx = {x: ch.idx(x) for x in ch.reps}
    # G1 = sum phi(x) e(x/f) as exponent vector
    G1 = [0] * M
    for x, j in idx.items():
        if j is not None:
            b = mul(x, fc)[1] % Q
            G1[(j * (M // 6) + b * (M // Q)) % M] += 1
    G1num = sum(c * cmath.exp(2j * math.pi * k / M) for k, c in enumerate(G1) if c)
    nerr = 0.0; ex_ok = True; ferr = 0.0
    # |gamma| = 1 exactly: G1 * conj(G1) = Q
    if exact:
        prod = [0] * M
        nz = [(k, c) for k, c in enumerate(G1) if c]
        for k1, c1 in nz:
            for k2, c2 in nz:
                prod[(k1 - k2) % M] += c1 * c2
        prod[0] -= Q
        ex_ok &= all(c == 0 for c in reduce_cyclo(prod, phiM))
    for h in ch.reps:
        vec = [0] * M
        for x, j in idx.items():
            if j is not None:
                b = mul(mul(h, x), fc)[1] % Q
                vec[(j * (M // 6) + b * (M // Q)) % M] += 1
        lhs = sum(c * cmath.exp(2j * math.pi * k / M) for k, c in enumerate(vec) if c)
        jh = ch.idx(h)
        if jh is None:
            rhs = 0
            fal = 0
        else:
            rhs = ZETA6[(-jh) % 6] * G1num
            fal = ZETA6[jh] * G1num
            if exact:
                sh = ((-jh) % 6) * (M // 6)
                for k, c in enumerate(G1):
                    if c:
                        vec[(k + sh) % M] -= c
        if exact:
            ex_ok &= all(c == 0 for c in reduce_cyclo(vec, phiM))
        nerr = max(nerr, abs(lhs - rhs))
        ferr = max(ferr, abs(lhs - fal))
    return nerr, ex_ok, ferr


tot_err = 0.0; tot_ok = True; fal_min = 1e9; lines = []
for ch in chars:
    ne, eo, fe = gauss_check(ch)
    tot_err = max(tot_err, ne); tot_ok &= eo
    if any(ch.idx(h) not in (None, 0, 3) for h in ch.reps):     # phi takes a non-real value
        fal_min = min(fal_min, fe)
    lines.append(f"{ch.name.split(':')[0]}(Q={ch.Q}): num {ne:.1e}, exact {eo}, falsifier {fe:.2f}")
# imprimitive character: u3 viewed modulo 3*pi7 (phi trivial on the pi7 part)
p7 = split_primary(7)[0]
imp = Char("imprimitive u3 mod 3 pi7", mul((3, 0), p7), 1, [])
ne_imp, eo_imp, _ = gauss_check(imp, exact=True)
print("   " + "; ".join(lines))
report("T7-C3", "eq:s2-gauss", "sum_x phi(x) e(hx/f) = conj(phi(h)) sqrt(Q) gamma for all h mod f, exactly in Z[zeta_M]; |gamma| = 1",
       tot_ok and tot_err < 1e-9,
       f"8 characters, all h mod f (incl. h not coprime to f); exact: {tot_ok}; numeric max error {tot_err:.1e}; "
       f"falsifiers: phi(h) for conj(phi(h)) min-over-characters max error {fal_min:.2f}; imprimitive u3 mod 3 pi7: "
       f"exact identity holds? {eo_imp} (numeric error {ne_imp:.2f})",
       falsifier_failed=(fal_min > 0.5 and not eo_imp and ne_imp > 0.5))
print(f"time {time.time() - t0:.1f}s")
finish()
