"""t8_coefficients.py -- displays of Appendix E (sec:s4-inverse) of cm-eleven-twelfths-guide.tex that are identities
about the column twist phi = conj(alpha) and the coefficient classes.

  T8.1  eq:d170   phi(n) = conj(alpha(n)) is completely multiplicative on primary elements and |phi| = 1
                  (the sentence after the display).  Falsifier: alpha evaluated at a non-primary associate.
  T8.2  eq:s4-a   a(n) = conj(alpha)^2 gamma_2 nu rho = phi * conj(alpha) gamma_2 nu rho, with nu = chi_.(lambda) and a
                  puncture rho; and |a| = 1 on squarefree n prime to the puncture.
  T8.3  eq:d174   phi(n1) conj(phi(n2)) = phi(A1) conj(phi(A2)) phi(x1) conj(phi(x2)) for n_a = A_a C t' x_a with C, t'
                  common to both sides.  Falsifier: different C on the two sides.
  T8.4  eq:s4-R0  the orientation step: for eps = -1 the residual coefficient mu alpha^eps nu chi_c(u)^eps is the R_0
                  coefficient mu conj(alpha) nu_0 conj(chi_c(u)) with nu_0 = nu; for eps = +1 its conjugate is, with
                  nu_0 = conj(nu); and the R_0 coefficient is phi(c) times OAI78's (eq:initial-residual-polynomial).
                  Includes rows u sharing primes with c (zero extensions).  Falsifier: an anti-oriented sum
                  mu alpha^eps nu chi_c(u)^(-eps) is NOT of this form; it has the twist alpha (checked positively too).
  T8.5  eq:d173   Lemma lem:s4-transfer (column twist) on random finite data: a fibrewise identity stays an identity
                  after multiplying by phi(n1) conj(phi(n2)); the factorization n_a = e_a x_a splits the twist.
                  Falsifier: an identity that holds globally but not fibrewise does not survive the twist.
  T8.12 eq:s4-overlap pointwise (s4 S4.6 checks it summed over the assignments): for squarefree n and one prime from
                  each of 1-3 disjoint slots, with j = (n, P), c = n0 P0, both signs eps, rows sharing primes.
                  Falsifier: surviving slot coefficients without conj(alpha(p)^eps).

Usage: python3 t8_coefficients.py   (helpers of tests_cm/part_B_lines_0550_0997 are imported read-only)
Exit status 1 if a check fails or a falsifier does not fail.
"""
import os, sys, math, cmath, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import *          # zw arithmetic, chi, mu, factor, primaries, gamma_j_fast, report, falsifier, alpha_c
random.seed(20261013)

def al(n): return 1.0 + 0j if n == (1, 0) else alpha_c(n)
def phi(n): return al(n).conjugate()
def muk(n): return 1 if n == (1, 0) else mu(n)
_g2 = {}
def g2(n):
    if n == (1, 0): return 1.0 + 0j
    if n not in _g2: _g2[n] = gamma_j_fast(2, n)
    return _g2[n]
def chi0(n, x):                     # chi_n(x) with chi_1 = 1
    return 1.0 + 0j if n == (1, 0) else complex(chi(n, x))
def chipow(n, x, e):                # chi_n(x)^e with zero extension (e = -1 means conj)
    c = chi0(n, x)
    if c == 0: return 0j
    return c if e == 1 else c.conjugate()
def coprime(a, b): return N(gcd_primary(a, b)) == 1

PRIM = primaries(3000)                         # primary elements coprime to 6 (all, not only squarefree)
SQF = [n for n in PRIM if mu(n) != 0]

# ------------------------------------------------------------------ T8.1 eq:d170
err = 0.0; ferr = 0.0; cnt = 0
small = [n for n in PRIM if N(n) <= 200]
def phi_bad(n):                     # alpha at a non-primary associate depending on n (breaks the normalization)
    u = UNITS[1 + (N(n) % 5)]
    return al(mul(u, n)).conjugate()
for a_ in small:
    for b_ in small:
        ab = mul(a_, b_)
        assert primary(ab)
        cnt += 1
        err = max(err, abs(phi(ab) - phi(a_)*phi(b_)), abs(abs(phi(ab)) - 1))
        ferr = max(ferr, abs(phi_bad(ab) - phi_bad(a_)*phi_bad(b_)))
report("T8.1", err, 1e-12, f"eq:d170: phi = conj(alpha) completely multiplicative on {cnt} pairs of primary elements "
       "(N <= 200 each, not necessarily coprime or squarefree), |phi| = 1")
falsifier("T8.1", "alpha at a non-primary associate", ferr, 1e-3)

# ------------------------------------------------------------------ T8.2 eq:s4-a
LAMB = (1, 2)
PUNCT = [p for p in PR_ODD if N(p) == 37][0]
def nu(n): return chi0(n, LAMB)
def rho(n): return 1 if coprime(n, PUNCT) else 0
err = 0.0; err_mod = 0.0
for n in SQF[:400]:
    a1 = phi(n)**2*g2(n)*nu(n)*rho(n)
    a2 = phi(n)*(phi(n)*g2(n)*nu(n)*rho(n))
    err = max(err, abs(a1 - a2))
    if rho(n): err_mod = max(err_mod, abs(abs(a1) - 1))
report("T8.2", max(err, err_mod), 1e-12, "eq:s4-a: conj(alpha)^2 gamma_2 nu rho = phi conj(alpha) gamma_2 nu rho and |a| = 1 off the "
       "puncture, 400 squarefree n, nu = chi_.(lambda), puncture of norm 37")

# ------------------------------------------------------------------ T8.3 eq:d174
err = 0.0; ferr = 0.0
pool = [n for n in PRIM if N(n) <= 400]
for _ in range(3000):
    A1, A2, C, t, x1, x2, C2 = (random.choice(pool) for _ in range(7))
    n1 = mul(mul(mul(A1, C), t), x1); n2 = mul(mul(mul(A2, C), t), x2)
    lhs = phi(n1)*phi(n2).conjugate()
    rhs = phi(A1)*phi(A2).conjugate()*phi(x1)*phi(x2).conjugate()
    err = max(err, abs(lhs - rhs))
    if C2 != C:
        n2b = mul(mul(mul(A2, C2), t), x2)       # C not common to both sides
        ferr = max(ferr, abs(phi(n1)*phi(n2b).conjugate() - rhs))
report("T8.3", err, 1e-12, "eq:d174: twist on n_a = A_a C t' x_a splits as phi(A1)conj phi(A2) phi(x1)conj phi(x2); 3000 random tuples")
falsifier("T8.3", "C different on the two sides", ferr, 1e-3)

# ------------------------------------------------------------------ T8.4 eq:s4-R0 (orientation)
pr7 = [p for p in PR_ODD if N(p) == 7]
ROWS = [(5, 3), (2, 9), LAMB, (0, 1), (-1, 0), pr7[0], mul(pr7[0], pr7[1]), mul(mul(pr7[0], pr7[0]), (4, 3)), (11, 0)]
csq = [n for n in SQF if N(n) <= 1500]
e_or = 0.0; e_oai = 0.0; f_anti = 0.0; e_anti_alpha = 0.0; nzero = 0
for eps in (1, -1):
    for u in ROWS:
        for c in csq:
            apow = al(c) if eps == 1 else phi(c)
            res = muk(c)*apow*nu(c)*chipow(c, u, eps)                     # residual coefficient mu alpha^eps psi_u
            if res == 0: nzero += 1
            if eps == 1:
                got = res.conjugate(); nu0 = nu(c).conjugate()                # conjugate the whole polynomial
            else:
                got = res; nu0 = nu(c)
            r0 = muk(c)*phi(c)*nu0*chipow(c, u, -1)                         # R_0 coefficient (eq:s4-R0)
            oai = muk(c)*nu0*chipow(c, u, -1)                               # OAI78 L11717 coefficient
            e_or = max(e_or, abs(got - r0))
            e_oai = max(e_oai, abs(r0 - phi(c)*oai))
            # anti-oriented sum: alpha^eps paired with chi^(-eps)
            anti = muk(c)*apow*nu(c)*chipow(c, u, -eps)
            got_a = anti.conjugate() if eps == -1 else anti                 # bring its row character to conj(chi)
            nu0a = nu(c).conjugate() if eps == -1 else nu(c)
            f_anti = max(f_anti, abs(got_a - muk(c)*phi(c)*nu0a*chipow(c, u, -1)))
            e_anti_alpha = max(e_anti_alpha, abs(got_a - muk(c)*al(c)*nu0a*chipow(c, u, -1)))
report("T8.4", max(e_or, e_oai, e_anti_alpha), 1e-12,
       f"eq:s4-R0: eps=+1,-1, {len(ROWS)} rows (units, lambda, rows sharing primes with c), {len(csq)} squarefree c; "
       f"residual coefficient = R_0 coefficient = phi(c) x OAI78's; {nzero} zero-extended values; anti-oriented sum "
       "has twist alpha")
falsifier("T8.4", "anti-oriented sum in the R_0 form with conj(alpha)", f_anti, 1e-3)

# ------------------------------------------------------------------ T8.5 eq:d173 (Lemma lem:s4-transfer)
def column_twist_trial(fibrewise=True):
    nfib = 25
    fib = [(random.choice(pool), random.choice(pool)) for _ in range(nfib)]      # lambda -> (n1, n2)
    taus = []; sigmas = []
    for lam in range(nfib):
        L = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(random.randint(1, 5))]
        R = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(random.randint(1, 5))]
        R[-1] += sum(L) - sum(R)                                          # fibrewise identity
        taus += [(lam, v) for v in L]; sigmas += [(lam, v) for v in R]
    if not fibrewise:                                                     # move mass between fibres 0 and 1
        d = complex(0.7, -0.4)
        sigmas.append((0, d)); sigmas.append((1, -d))
    tw = lambda lam: phi(fib[lam][0])*phi(fib[lam][1]).conjugate()
    glob = abs(sum(v for _, v in taus) - sum(v for _, v in sigmas))
    twisted = abs(sum(tw(l)*v for l, v in taus) - sum(tw(l)*v for l, v in sigmas))
    # factorization n_a = e_a x_a (sigma-dependent): twist splits
    fac = 0.0
    for lam in range(nfib):
        n1, n2 = fib[lam]
        for e1 in [d for d in pool if divides(d, n1)][:4]:
            for e2 in [d for d in pool if divides(d, n2)][:4]:
                x1, x2 = div(n1, e1), div(n2, e2)
                fac = max(fac, abs(tw(lam) - phi(e1)*phi(e2).conjugate()*phi(x1)*phi(x2).conjugate()))
    return glob, twisted, fac
g_err = t_err = f_err = 0.0
for _ in range(20):
    g, t, f = column_twist_trial(True)
    g_err = max(g_err, g); t_err = max(t_err, t); f_err = max(f_err, f)
g_bad = t_bad = 0.0
for _ in range(20):
    g, t, _f = column_twist_trial(False)
    g_bad = max(g_bad, g); t_bad = max(t_bad, t)
report("T8.5", max(g_err, t_err, f_err), 1e-12, "eq:d173: 20 random fibrewise identities (25 fibres, random primary column "
       "pairs) survive the twist phi(n1)conj phi(n2); twist splits over n_a = e_a x_a")
report("T8.5-control", g_bad, 1e-12, "the non-fibrewise variant still satisfies the untwisted identity (control)")
falsifier("T8.5", "identity true globally but not fibrewise, after the twist", t_bad, 1e-3)

# ------------------------------------------------------------------ T8.12 eq:s4-overlap, pointwise
# For squarefree n and one prime p_i from each slot (distinct): j = (n, P), P = prod p_i, n = j n0, P = j P0, c = n0 P0,
# I0 = surviving slots (p_i not dividing n).  Display:
#   mu(n) alpha(n)^eps psi_u(n) prod_I a_i(p_i) psi_u(p_i)
#     = [mu(j) alpha(j)^eps psi_u(j)^2 prod_{i not in I0} a_i(p_i)] mu(c) alpha(c)^eps psi_u(c) prod_{I0} mu(p_i) conj(alpha(p_i)^eps) a_i(p_i)
pr_small = [p for p in PR_ODD if N(p) <= 200]
def apw(x, eps): return al(x) if eps == 1 else phi(x)
e12 = 0.0; f12 = 0.0; cnt = 0; nsurv = [0, 0, 0, 0]
for _ in range(4000):
    eps = random.choice((1, -1)); u = random.choice(ROWS + [(1, 0), (2, -5)])
    k_ = random.randint(1, 3)
    ps = random.sample(pr_small, k_)
    a_ = [complex(random.gauss(0, 1), random.gauss(0, 1)) for _ in ps]
    n = random.choice(csq)
    if random.random() < 0.6:                      # force overlaps with some slot primes
        for p in random.sample(ps, random.randint(1, k_)):
            if not divides(p, n) and N(n)*N(p) <= 40000: n = mul(n, p)
    if mu(n) == 0: continue
    psi_ = lambda x: nu(x)*chipow(x, u, eps) if x != (1, 0) else 1.0 + 0j
    lhs = muk(n)*apw(n, eps)*psi_(n)
    for p, a in zip(ps, a_): lhs *= a*psi_(p)
    j = (1, 0); P0 = (1, 0); assigned = 1.0 + 0j; surv = 1.0 + 0j; surv_bad = 1.0 + 0j; ns = 0
    for p, a in zip(ps, a_):
        if divides(p, n):
            j = mul(j, p); assigned *= a
        else:
            P0 = mul(P0, p); ns += 1
            surv *= muk(p)*apw(p, eps).conjugate()*a
            surv_bad *= muk(p)*a
    n0 = div(n, j); c_ = mul(n0, P0)
    rhs_core = muk(j)*apw(j, eps)*psi_(j)**2*assigned*muk(c_)*apw(c_, eps)*psi_(c_)
    e12 = max(e12, abs(lhs - rhs_core*surv)); cnt += 1; nsurv[ns] += 1
    if ns > 0: f12 = max(f12, abs(lhs - rhs_core*surv_bad))
report("T8.12", e12, 1e-12, f"eq:s4-overlap pointwise: {cnt} cases, eps = +-1, 1-3 slots, rows incl. units and rows sharing primes; "
       f"cases with 0/1/2/3 surviving slots: {nsurv}")
falsifier("T8.12", "surviving slot coefficient without conj(alpha(p)^eps)", f12, 1e-3)
