"""t2a: the finite Gauss-sum and character identities of Section 7 of the CM guide (displays (7.8)-(7.30)).

Ids are the guide labels.  Every check prints RESULT|<label>|PASS/FAIL|...; every falsifier (a deliberately wrong
variant) must fail.  Exit status 1 if any check fails or any falsifier does not fail (cmlib atexit hook).

  eq:b2-g1g2   conj(alpha(n)) gamma_2(n) gamma_1(n) = mu_K(n) G(n); |gamma_j(n)| = 1 for 6 not dividing j
  eq:b2-crt    gamma(conj chi_z1 chi_z2) = R(z1,z2) gamma_-1(z1) gamma_1(z2), and the reversed version
  eq:b2-two    mu gamma_-1 = chi_n(-1) conj(G) a^zeta,  mu gamma_1 = G conj(a^zeta),  a^zeta = conj(alpha) gamma_2
  eq:d040      both equalities (middle term with chi_z1(-1), R; right side G(z2 z1^-1))
  eq:d041      both identities
  eq:d042      reverse class form with a character xi of frak G (6 characters)
  eq:b2-forward, eq:b2-reverse   literal sums over the 432 characters of an independent character table
  eq:d039      the coefficients c_xi(xi_1), l_xi1(xi): expansions, c_xi(1) = c_xi, all |.| <= 1
  eq:d037      c_xi from its definition: sum_xi c_xi xi(t) = G(t^-1), sum_xi c_xi = G(1) = 1, |c_xi| <= 1
  eq:b2-twisted / eq:d038   a_xi(ab) = a_xi(a) a_xi(b) chi_b(a)^4 and |a_xi| = 1
  eq:d044      mu(n1)mu(n2) = mu(z1)mu(z2), conj(alpha(n1))alpha(n2) = conj(alpha(z1))alpha(z2)  (n_i = g z_i)
  eq:d046      mu(n1)mu(n2)conj(alpha(n1))alpha(n2) gamma(psi) = sum_xi c_xi a_xi(z1) conj(a_xi(z2)), g != 1 included;
               z1 = z2 = 1 gives sum_xi c_xi = G(1) = 1
Run: python3 -B t2a_gauss_identities.py   (about 10 s)
"""
import os, sys, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t2lib as L
from t2lib import *
T0 = time.time()
random.seed(20261013)
sq = primaries(3000, sqfree=True)
sq_small = [n for n in sq if N(n) <= 400]
NS = sq[:150]                                   # 150 squarefree n incl. 1 and composites
ncomp = sum(1 for n in NS if len(factor(n)) > 1)

# brute-force control of the vectorized Gauss sums
e_bf = max(abs(gamma_j_fast(j, n) - gamma_j(j, n)) for n in sq[1:30] for j in (-1, 1, 2, 3))
report("helper", e_bf, 1e-12, "vectorized gamma_j = brute-force gamma_j (common.py) for 29 n, j = -1, 1, 2, 3")

# ------------------------------------------------------------------ eq:b2-g1g2
err = 0; e_abs = 0; f1 = f2 = 0
for n in NS:
    lhs = alpha_c(n).conjugate()*gam(2, n)*gam(1, n) if n != (1, 0) else 1.0
    rhs = muk(n)*Gn(n)
    err = max(err, abs(lhs - rhs))
    if n != (1, 0):
        e_abs = max(e_abs, max(abs(abs(gam(j, n)) - 1) for j in (1, 2, 3, 4, 5, -1)))
        f1 = max(f1, abs(alpha_c(n)*gam(2, n)*gam(1, n) - rhs))
        f2 = max(f2, abs(lhs - Gn(n)))
report("eq:b2-g1g2", err, 1e-11, f"{len(NS)} squarefree n ({ncomp} composite), N(n) <= {N(NS[-1])}")
falsifier("eq:b2-g1g2", "alpha(n) in place of conj(alpha(n))", f1, 1e-3)
falsifier("eq:b2-g1g2", "Moebius sign dropped", f2, 1e-3)
report("gamma_j-unimodular", e_abs, 1e-11, f"|gamma_j(n)| = 1 for j = +-1, 2, 3, 4, 5 and the same {len(NS)-1} n (text before eq:b2-g1g2)")

# ------------------------------------------------------------------ pairs
pairs = []
while len(pairs) < 80:
    z1, z2 = random.choice(sq_small), random.choice(sq_small)
    if not coprime(z1, z2) or N(z1)*N(z2) > 5000: continue
    pairs.append((z1, z2))
pairs += [((1, 0), (7, 9)), ((-5, 0), (1, 0)), ((1, 0), (1, 0))]
GF = {p: gauss_psi_general(p[0], p[1], conj_first=True) for p in pairs}
GR = {p: gauss_psi_general(p[0], p[1], conj_first=False) for p in pairs}
nR = sum(1 for z1, z2 in pairs if abs(Rf(z1, z2) + 1) < 1e-9)

# ------------------------------------------------------------------ eq:b2-crt
err = eR = 0; f1 = f2 = 0
for (z1, z2) in pairs:
    R = Rf(z1, z2)
    eR = max(eR, min(abs(R - 1), abs(R + 1)))
    rf = R*gam(-1, z1)*gam(1, z2); rr = R*gam(1, z1)*gam(-1, z2)
    err = max(err, abs(GF[(z1, z2)] - rf), abs(GR[(z1, z2)] - rr))
    f1 = max(f1, abs(GF[(z1, z2)] - gam(-1, z1)*gam(1, z2)))
    f2 = max(f2, abs(GF[(z1, z2)] - R*gam(1, z1)*gam(-1, z2)))
report("eq:b2-crt", max(err, eR), 1e-11, f"{len(pairs)} coprime squarefree pairs, N(z1 z2) <= 5000, brute-force Gauss sums mod z1 z2; "
       f"both directions; R(z1,z2) = chi_z2(z1)/chi_z1(z2) is +-1 ({nR} pairs with R = -1)")
falsifier("eq:b2-crt", "R(z1,z2) dropped", f1, 1e-3)
falsifier("eq:b2-crt", "gamma_1(z1) gamma_-1(z2) (directions swapped)", f2, 1e-3)

# ------------------------------------------------------------------ eq:b2-two
e1 = e2 = 0; f1 = f2 = 0
for n in NS:
    az = alpha_c(n).conjugate()*gam(2, n)
    l1 = muk(n)*gam(-1, n); r1 = complex(chi(n, (-1, 0)))*Gn(n).conjugate()*az
    l2 = muk(n)*gam(1, n); r2 = Gn(n)*az.conjugate()
    e1 = max(e1, abs(l1 - r1)); e2 = max(e2, abs(l2 - r2))
    f1 = max(f1, abs(l1 - complex(chi(n, (-1, 0)))*Gn(n)*az)); f2 = max(f2, abs(l2 - Gn(n)*az))
report("eq:b2-two", max(e1, e2), 1e-11, f"both identities, {len(NS)} squarefree n (errors {e1:.1e}, {e2:.1e})")
falsifier("eq:b2-two", "G(n) for conj(G(n)) in the first", f1, 1e-3)
falsifier("eq:b2-two", "a^zeta unconjugated in the second", f2, 1e-3)

# ------------------------------------------------------------------ eq:d040
ea = eb = 0; f1 = f2 = 0
for (z1, z2) in pairs:
    az1 = alpha_c(z1).conjugate()*gam(2, z1); az2 = alpha_c(z2).conjugate()*gam(2, z2)
    lhs = muk(z1)*muk(z2)*GF[(z1, z2)]
    mid = complex(chi(z1, (-1, 0)))*Gn(z1).conjugate()*Gn(z2)*Rf(z1, z2)*az1*az2.conjugate()
    rhs = Gq(z2, z1)*az1*az2.conjugate()
    ea = max(ea, abs(lhs - mid)); eb = max(eb, abs(mid - rhs))
    f1 = max(f1, abs(lhs - Gq(z1, z2)*az1*az2.conjugate()))
    f2 = max(f2, abs(lhs - mid/complex(chi(z1, (-1, 0)))))
report("eq:d040", max(ea, eb), 1e-11, f"both equalities, {len(pairs)} pairs (first {ea:.1e}, second {eb:.1e}); G(z2 z1^-1) on classes")
falsifier("eq:d040", "G(z1 z2^-1) on the right", f1, 1e-3)
falsifier("eq:d040", "chi_z1(-1) dropped from the middle", f2, 1e-3)

# ------------------------------------------------------------------ eq:d041
ea = eb = 0; f1 = f2 = 0
for n in NS:
    if n == (1, 0): continue
    al = alpha_c(n)
    la = al.conjugate()**2*gam(2, n)*gam(1, n); ra = al.conjugate()*muk(n)*Gn(n)
    lb = (al.conjugate()**2*gam(2, n)).conjugate()*gam(-1, n); rb = al*muk(n)*complex(chi(n, (-1, 0)))*Gn(n).conjugate()
    ea = max(ea, abs(la - ra)); eb = max(eb, abs(lb - rb))
    f1 = max(f1, abs(la - al*muk(n)*Gn(n))); f2 = max(f2, abs(lb - al.conjugate()*muk(n)*complex(chi(n, (-1, 0)))*Gn(n).conjugate()))
report("eq:d041", max(ea, eb), 1e-11, f"both identities, {len(NS)-1} squarefree n != 1 (errors {ea:.1e}, {eb:.1e})")
falsifier("eq:d041", "alpha for conj(alpha) on the right of the first", f1, 1e-3)
falsifier("eq:d041", "conj(alpha) for alpha on the right of the second", f2, 1e-3)

# ------------------------------------------------------------------ eq:d042 (with characters xi of frak G)
XIS = [TRIV] + random.sample(range(432), 5)
def xi(k, n): return CHARS[k, cls(n)]
def axi(k, n): return Acm(n)*xi(k, n)
err = 0; f1 = 0
for k in XIS:
    for (z1, z2) in pairs:
        lhs = axi(k, z1)*axi(k, z2).conjugate()*GR[(z1, z2)]
        t = mulm(red(z1), INV[red(z2)]); ts = mulm(red(z2), INV[red(z1)])
        rhs = muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*CHARS[k, IDX[t]]*Gcls(t)
        err = max(err, abs(lhs - rhs))
        f1 = max(f1, abs(lhs - muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*CHARS[k, IDX[ts]]*Gcls(ts)))
report("eq:d042", err, 1e-11, f"{len(pairs)} pairs x 6 characters xi of frak G (trivial + 5 random)")
falsifier("eq:d042", "(xi G)(z2 z1^-1) on the right", f1, 1e-3)

# ------------------------------------------------------------------ eq:d039 and eq:d037: coefficient tables
X = CHARS                                                   # X[k, t]
# c_xi(xi1) = |G|^-1 sum_t xi1(t) G(t^-1) conj(xi(t));  C[k1, k]
C = (X*GINV[None, :]) @ X.conj().T/432.0
Lm = (X*GVEC[None, :]) @ X.conj().T/432.0                   # Lm[k, k1] = l_{xi_k1}(xi_k) = |G|^-1 sum_t (xi_k G)(t) conj(xi_k1(t))
rec_c = np.max(np.abs(C @ X - X*GINV[None, :]))              # sum_xi c_xi(xi1) xi(t) = xi1(t) G(t^-1)
rec_l = np.max(np.abs(Lm @ X - X*GVEC[None, :]))             # sum_xi1 l_xi1(xi) xi1(t) = (xi G)(t)
e_triv = np.max(np.abs(C[TRIV, :] - CXI))
maxc = max(np.abs(C).max(), np.abs(Lm).max())
fal = np.max(np.abs(((X*GVEC[None, :]) @ X.conj().T/432.0) @ X - X*GINV[None, :]))
report("eq:d039", max(rec_c, rec_l, e_triv, max(0.0, maxc - 1)), 1e-11,
       f"all 432 x 432 coefficients: expansions reproduce xi1(t)G(t^-1) and (xi G)(t) for all t ({rec_c:.1e}, {rec_l:.1e}); "
       f"c_xi(1) = c_xi ({e_triv:.1e}); max |c|, |l| = {maxc:.3f} <= 1")
falsifier("eq:d039", "c_xi(xi1) computed with G(t) for G(t^-1)", fal, 1e-3)
e37 = max(np.max(np.abs(CXI @ X - GINV)), abs(CXI.sum() - 1), max(0.0, np.abs(CXI).max() - 1))
report("eq:d037", e37, 1e-11, f"c_xi from its definition (independent table, basis orders {ORDERS}): sum_xi c_xi xi(t) = G(t^-1) for all 432 t; "
       f"sum c_xi = G(1) = 1; max |c_xi| = {np.abs(CXI).max():.3f}")

# ------------------------------------------------------------------ eq:b2-forward, eq:b2-reverse (literal sums)
XI1S = [TRIV] + random.sample(range(432), 3)
ef = er = 0; f1 = f2 = f3 = 0
for (z1, z2) in pairs[:60]:
    i1, i2 = cls(z1), cls(z2)
    a1, a2 = Acm(z1), Acm(z2); b1, b2 = Azeta(z1), Azeta(z2)
    pre = muk(z1)*muk(z2)*alpha_c(z1).conjugate()*alpha_c(z2)*GF[(z1, z2)]
    for k1 in XI1S:
        lhs = pre*X[k1, i1]*X[k1, i2].conjugate()
        rhs = np.sum(C[k1, :]*X[:, i1]*X[:, i2].conj())*a1*a2.conjugate()
        ef = max(ef, abs(lhs - rhs))
        f1 = max(f1, abs(lhs - np.sum(C[k1, :]*X[:, i1]*X[:, i2].conj())*b1*b2.conjugate()))
        f3 = max(f3, abs(pre*X[k1, i2]*X[k1, i1].conjugate() - rhs))
    for k in XIS[:4]:
        lhs = a1*X[k, i1]*(a2*X[k, i2]).conjugate()*GR[(z1, z2)]
        beta1 = muk(z1)*alpha_c(z1).conjugate()*X[:, i1]; beta2 = muk(z2)*alpha_c(z2).conjugate()*X[:, i2]
        rhs = np.sum(Lm[k, :]*beta1*beta2.conj())
        er = max(er, abs(lhs - rhs))
        g1 = muk(z1)*X[:, i1]; g2 = muk(z2)*X[:, i2]
        f2 = max(f2, abs(lhs - np.sum(Lm[k, :]*g1*g2.conj())))
report("eq:b2-forward", ef, 1e-10, f"60 pairs x 4 characters xi1 (incl. trivial), literal sum over the 432 xi with c_xi(xi1)")
falsifier("eq:b2-forward", "zeta coefficients conj(alpha) gamma_2 in a_xi", f1, 1e-3)
falsifier("eq:b2-forward", "xi1 conjugated on the wrong column", f3, 1e-3)
report("eq:b2-reverse", er, 1e-10, f"60 pairs x 4 characters xi, literal sum over the 432 xi1 with l_xi1(xi), beta = mu conj(alpha) xi1")
falsifier("eq:b2-reverse", "beta = mu xi1 (no conj(alpha), zeta form)", f2, 1e-3)

# ------------------------------------------------------------------ eq:b2-twisted / eq:d038
err = eabs = 0; f1 = 0; cnt = 0
while cnt < 120:
    a, b = random.choice(sq_small), random.choice(sq_small)
    if not coprime(a, b) or N(a)*N(b) > 6000: continue
    cnt += 1; ab = mul(a, b)
    for k in XIS:
        lhs = axi(k, ab); rhs = axi(k, a)*axi(k, b)*complex(chi(b, a))**4
        err = max(err, abs(lhs - rhs)); eabs = max(eabs, abs(abs(axi(k, a)) - 1), abs(abs(lhs) - 1))
        f1 = max(f1, abs(lhs - axi(k, a)*axi(k, b)*complex(chi(b, a))**2))
report("eq:b2-twisted", max(err, eabs), 1e-11, f"120 coprime squarefree pairs (N(ab) <= 6000) x 6 characters xi: twisted multiplicativity and |a_xi| = 1")
falsifier("eq:b2-twisted", "chi_b(a)^2 in place of chi_b(a)^4", f1, 1e-3)
report("eq:d038", max(err, eabs), 1e-11, "same data as eq:b2-twisted with (v0, n) = (a, b)")

# ------------------------------------------------------------------ eq:d044 and eq:d046
gl = [g for g in sq_small if N(g) <= 200]
trip = []
while len(trip) < 70:
    g = random.choice(gl); z1, z2 = random.choice(sq_small), random.choice(sq_small)
    if not (coprime(z1, z2) and coprime(g, mul(z1, z2))) or N(z1)*N(z2) > 4000: continue
    trip.append((g, z1, z2))
trip += [((-5, 0), (1, 0), (1, 0)), (mul((-1, -3), (4, 3)), (1, 0), (7, 9))]
e44 = 0; f44 = 0; e46 = 0; f46a = f46b = 0; ng = sum(1 for g, _, _ in trip if g != (1, 0))
for (g, z1, z2) in trip:
    n1, n2 = mul(g, z1), mul(g, z2)
    e44 = max(e44, abs(muk(n1)*muk(n2) - muk(z1)*muk(z2)),
              abs(alpha_c(n1).conjugate()*alpha_c(n2) - alpha_c(z1).conjugate()*alpha_c(z2)))
    f44 = max(f44, abs(alpha_c(n1)*alpha_c(n2) - alpha_c(z1)*alpha_c(z2)))
    gp = gauss_psi_general(z1, z2, conj_first=True)
    lhs = muk(n1)*muk(n2)*alpha_c(n1).conjugate()*alpha_c(n2)*gp
    i1, i2 = cls(z1), cls(z2)
    rhs = np.sum(CXI*X[:, i1]*X[:, i2].conj())*Acm(z1)*Acm(z2).conjugate()
    e46 = max(e46, abs(lhs - rhs))
    f46a = max(f46a, abs(lhs - np.sum(CXI*X[:, i1]*X[:, i2].conj())*Azeta(z1)*Azeta(z2).conjugate()))
    f46b = max(f46b, abs(muk(n1)*muk(n2)*alpha_c(n1)*alpha_c(n2).conjugate()*gp - rhs))
report("eq:d044", e44, 1e-12, f"{len(trip)} triples (g, z1, z2), g squarefree with 0-2 primes, (g, z1 z2) = (z1, z2) = 1")
falsifier("eq:d044", "alpha(n1) alpha(n2) = alpha(z1) alpha(z2) (no conjugation)", f44, 1e-3)
one = abs(CXI.sum() - 1)
report("eq:d046", max(e46, one), 1e-10, f"{len(trip)} triples ({ng} with g != 1, incl. z1 = 1 and z1 = z2 = 1), literal sum over the 432 xi; "
       f"z1 = z2 = 1: sum c_xi = G(1) = 1 (err {one:.1e})")
falsifier("eq:d046", "zeta coefficients conj(alpha) gamma_2 on the right", f46a, 1e-3)
falsifier("eq:d046", "alpha(n1) conj(alpha(n2)) on the left", f46b, 1e-3)
print(f"t2a done in {time.time() - T0:.1f} s")
