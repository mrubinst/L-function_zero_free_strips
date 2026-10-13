"""Checks for inserts.tex (agent s4): the CM inverse moment eq:cm-inverse for co-oriented coefficients.

Every identity that inserts.tex writes out is tested here, each with a falsifier (a deliberately wrong variant
that must fail).  Helpers (exact Eisenstein arithmetic, sextic symbols, vectorized Gauss sums) are those of
tests_cm/part_B_lines_0550_0997 (cmlib.py, common.py, zw.py, fastgauss.py).

Conventions (OAI78 Section sec:arithmetic = the guide): primary = 1 mod 3; chi_n(x) = (x/n)_6, zero on nonunits;
e(z) = exp(4 pi i Im z/sqrt 3); gamma_j(n) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n); G(n) = conj(chi_n(4))
gamma_3(n); alpha(n) = n/|n|; phi = conj(alpha) is the column twist.

  S4.1  switching identities (s4-S1)-(s4-S4) pointwise on squarefree u (composites included)
  S4.2  anti-oriented collapse (s4-S5): mu alpha gamma_{-1} = chi(-1) conj(G) gamma_2 (no angular factor left)
  S4.3  twisted multiplicativity of a = conj(alpha)^2 gamma_2 nu, nu a ray character
  S4.4  ray-quotient forms of the two class switches (s4-P1), (s4-P2) for coprime squarefree pairs
  S4.5  CM marked cube inversion (s4-inversion), with nu, a puncture and a product-form mark
  S4.6  CM overlap identity (s4-overlap) of the initialization, both signs epsilon
  S4.7  CM sixth-power column identity (s4-sixth) of the amplification, both signs epsilon
"""
import os, sys, math, cmath, itertools, random
LIB = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   '..', '..', 'part_B_lines_0550_0997'))
sys.path.insert(0, LIB)
from cmlib import *          # numpy as np, zw helpers, chi, mu, factor, primaries, gamma_j_fast, G_fast, ...

random.seed(7)
RESULTS = []
_FORCE = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]   # self-test (added 2026-10-12): first check counts as failed

def line(name, err, tol, extra=""):
    ok = err <= tol
    if _FORCE[0]:
        _FORCE[0] = False; ok = False; extra = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + extra
    RESULTS.append((name, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: err {err:.3g} (tol {tol:g}) {extra}", flush=True)

def fals(name, err, thr=1e-6):
    ok = err > thr
    if _FORCE[0]:
        _FORCE[0] = False; ok = False; name = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + name
    RESULTS.append(("falsifier " + name, ok))
    print(f"      falsifier [{name}]: err {err:.3g} -> {'fails as required' if ok else 'DID NOT FAIL'}", flush=True)

def al(n): return 1.0 + 0j if n == (1, 0) else alpha_c(n)
def phi(n): return al(n).conjugate()                 # the column twist conj(alpha)
def muu(n): return 1 if n == (1, 0) else mu(n)
_gc = {}
def gam(j, n):
    if n == (1, 0): return 1.0 + 0j
    key = (j % 6, n)
    if key not in _gc: _gc[key] = gamma_j_fast(j % 6, n)
    return _gc[key]
def Gf(n): return 1.0 + 0j if n == (1, 0) else chi(n, (4, 0)).conjugate()*gam(3, n)
def chim1(n): return 1.0 + 0j if n == (1, 0) else chi(n, (-1, 0))
def coprime(a, b): return N(gcd_primary(a, b)) == 1

# squarefree primary u coprime to 6, composites included
SQF = [n for n in primaries(1300, sqfree=True) if n != (1, 0)]
COMP = [n for n in SQF if len(factor(n)) >= 2]
print(f"{len(SQF)} squarefree u with N(u) <= 1300, of which {len(COMP)} composite "
      f"(max number of prime factors {max(len(factor(n)) for n in SQF)})")

# ------------------------------------------------------------------ S4.1 switching identities
e1 = e2 = e3 = e4 = 0.0; f1 = f2 = f3 = f4 = 0.0
for u in SQF:
    a = phi(u)**2*gam(2, u)                          # CM cubic coefficient conj(alpha)^2 gamma_2
    m = muu(u)*phi(u)                                # CM Moebius coefficient mu conj(alpha)
    # (S1) first canonical Poisson, side 1: a gamma_1 = mu conj(alpha) G
    e1 = max(e1, abs(a*gam(1, u) - m*Gf(u)))
    f1 = max(f1, abs(phi(u)*gam(2, u)*gam(1, u) - m*Gf(u)))           # zeta class on the left
    # (S2) first canonical Poisson, side 2: conj(a) gamma_{-1} = conj(mu conj(alpha) chi(-1) G)
    e2 = max(e2, abs(a.conjugate()*gam(-1, u) - (m*chim1(u)*Gf(u)).conjugate()))
    f2 = max(f2, abs(a.conjugate()*gam(-1, u) - (m*chim1(u)*Gf(u))))   # forgot to conjugate
    # (S3) second Poisson / initialization, side 1: mu conj(alpha) gamma_{-1} = a chi(-1) conj(G)
    e3 = max(e3, abs(m*gam(-1, u) - a*chim1(u)*Gf(u).conjugate()))
    f3 = max(f3, abs(m*gam(-1, u) - a*chim1(u)*Gf(u)))                # G instead of conj G
    # (S4) second Poisson / initialization, side 2: conj(m) gamma_1 = conj(a) G
    e4 = max(e4, abs(m.conjugate()*gam(1, u) - a.conjugate()*Gf(u)))
    f4 = max(f4, abs(muu(u)*gam(1, u) - a.conjugate()*Gf(u)))          # zeta Moebius class
line("S4.1 (s4-S1) conj(alpha)^2 gamma_2 gamma_1 = mu conj(alpha) G", e1, 1e-10, f"{len(SQF)} u")
fals("S4.1 zeta cubic class conj(alpha) gamma_2 on the left", f1)
line("S4.1 (s4-S2) conj(conj(alpha)^2 gamma_2) gamma_{-1} = conj(mu conj(alpha) chi(-1) G)", e2, 1e-10)
fals("S4.1 side 2 without the conjugation", f2)
line("S4.1 (s4-S3) mu conj(alpha) gamma_{-1} = conj(alpha)^2 gamma_2 chi(-1) conj(G)", e3, 1e-10)
fals("S4.1 G in place of conj(G)", f3)
line("S4.1 (s4-S4) conj(mu conj(alpha)) gamma_1 = conj(conj(alpha)^2 gamma_2) G", e4, 1e-10)
fals("S4.1 zeta Moebius class mu on the left", f4)

# ------------------------------------------------------------------ S4.2 anti-oriented collapse
e5 = 0.0; f5 = 0.0
for u in SQF:
    lhs = muu(u)*al(u)*gam(-1, u)                    # anti-oriented: alpha with conj(chi)
    e5 = max(e5, abs(lhs - gam(2, u)*chim1(u)*Gf(u).conjugate()))
    f5 = max(f5, abs(lhs - phi(u)**2*gam(2, u)*chim1(u)*Gf(u).conjugate()))
line("S4.2 (s4-S5) mu alpha gamma_{-1} = gamma_2 chi(-1) conj(G) (no angular factor)", e5, 1e-10)
fals("S4.2 claim that the anti-oriented class is conj(alpha)^2 gamma_2", f5)

# ------------------------------------------------------------------ S4.3 twisted multiplicativity
LAM = (1, 2)
def nu_ray(n): return 1.0 + 0j if n == (1, 0) else chi(n, LAM)    # n -> chi_n(lambda): a ray character on 3
def acm(n, nu=nu_ray): return phi(n)**2*gam(2, n)*nu(n)
e6 = 0.0; f6 = 0.0; cnt = 0
small = [n for n in SQF if N(n) <= 160]
for a_, b_ in itertools.combinations(small, 2):
    if not coprime(a_, b_) or N(a_)*N(b_) > 1300: continue
    ab = mul(a_, b_); cnt += 1
    e6 = max(e6, abs(acm(ab) - acm(a_)*acm(b_)*chi(b_, a_)**4))
    f6 = max(f6, abs(acm(ab) - acm(a_)*acm(b_)*chi(b_, a_)**2))
line("S4.3 (eq:b2-twisted) a(ab) = a(a) a(b) chi_b(a)^4 for a = conj(alpha)^2 gamma_2 chi_.(lambda)", e6, 1e-10,
     f"{cnt} coprime pairs")
fals("S4.3 chi_b(a)^2 in place of chi_b(a)^4", f6)

# ------------------------------------------------------------------ S4.4 ray-quotient forms (P1), (P2)
M72 = 72
els = [(a, b) for a in range(M72) for b in range(M72) if a % 3 == 1 and b % 3 == 0 and not (a % 2 == 0 and b % 2 == 0)]
idx = {e: i for i, e in enumerate(els)}
def red(x): return (x[0] % M72, x[1] % M72)
order = len(els)
mt = np.array([[idx[red(mul(x, y))] for y in els] for x in els])
one = idx[(1, 0)]
inv = np.array([int(np.where(mt[i] == one)[0][0]) for i in range(order)])
def centered(v): return v - M72 if v > M72//2 else v
Gtab = np.zeros(order, complex)
for e in els:
    a0, b0 = centered(e[0]), centered(e[1]); best = None
    for i in range(-3, 4):
        for j in range(-3, 4):
            n = (a0 + M72*i, b0 + M72*j)
            if n == (1, 0): best = n; break
            if N(n) >= 2 and mu(n) != 0 and (best is None or N(n) < N(best)): best = n
        if best == (1, 0): break
    Gtab[idx[e]] = 1.0 if best == (1, 0) else Gf(best)
def Gq(a, b): return Gtab[mt[idx[red(a)], inv[idx[red(b)]]]]      # G([a][b]^{-1})
# G is a class function: two squarefree representatives of one class give the same value
gerr = 0.0
for n in SQF[:150]:
    gerr = max(gerr, abs(Gf(n) - Gtab[idx[red(n)]]))
line("S4.4 G(n) depends only on the class of n mod 72 (squarefree n)", gerr, 1e-10)
pairs = []
for a_, b_ in itertools.combinations([n for n in SQF if N(n) <= 200], 2):
    if coprime(a_, b_) and N(a_)*N(b_) <= 4000: pairs.append((a_, b_))
random.shuffle(pairs); pairs = pairs[:70] + [((1, 0), n) for n in SQF[:5]] + [(n, (1, 0)) for n in SQF[:5]]
eP1 = eP2 = fP1 = fP2 = 0.0
for z1, z2 in pairs:
    A1, A2 = acm(z1, lambda n: 1.0), acm(z2, lambda n: 1.0)
    M1, M2 = muu(z1)*phi(z1), muu(z2)*phi(z2)
    g_fwd = gauss_psi_general(z1, z2, conj_first=False) if mul(z1, z2) != (1, 0) else 1.0   # gamma(chi_z1 conj chi_z2)
    g_bwd = gauss_psi_general(z1, z2, conj_first=True) if mul(z1, z2) != (1, 0) else 1.0    # gamma(conj chi_z1 chi_z2)
    # (P1) cubic -> Moebius (first canonical Poisson)
    eP1 = max(eP1, abs(A1*A2.conjugate()*g_fwd - M1*M2.conjugate()*Gq(z1, z2)))
    fP1 = max(fP1, abs(phi(z1)*gam(2, z1)*(phi(z2)*gam(2, z2)).conjugate()*g_fwd - M1*M2.conjugate()*Gq(z1, z2)))
    # (P2) Moebius -> cubic (second Poisson, initialization)
    eP2 = max(eP2, abs(M1*M2.conjugate()*g_bwd - A1*A2.conjugate()*Gq(z2, z1)))
    fP2 = max(fP2, abs(M1*M2.conjugate()*g_bwd - A1*A2.conjugate()*Gq(z1, z2)))
line("S4.4 (s4-P1) a(z1)conj(a(z2)) gamma(chi_z1 conj chi_z2) = m(z1)conj(m(z2)) G([z1][z2]^-1)", eP1, 1e-10,
     f"{len(pairs)} pairs")
fals("S4.4 zeta cubic class in (s4-P1)", fP1)
line("S4.4 (s4-P2) m(z1)conj(m(z2)) gamma(conj chi_z1 chi_z2) = a(z1)conj(a(z2)) G([z2][z1]^-1)", eP2, 1e-10)
fals("S4.4 G([z1][z2]^-1) in (s4-P2)", fP2)

# ------------------------------------------------------------------ S4.5 CM marked cube inversion
def bumpW(y):            # smooth annular test on (1/2, 2), complex
    return math.exp(-1/((y - 0.5)*(2.0 - y)))*cmath.exp(0.9j*y) if 0.5 < y < 2.0 else 0.0
def Vstar(y): return math.sqrt(y)*bumpW(y)
PRIMES = [p for p in PR_ODD if N(p) <= 6000]
def check_inversion(X, k, f, r_punct, slots, bexp=6, hexp=6, shift_mark=True):
    """Return |LHS - RHS| and |LHS| for (s4-inversion)."""
    def rho(n): return 1 if n == (1, 0) else (1 if coprime(n, r_punct) else 0)
    def Psi(n): return 1.0 + 0j if n == (1, 0) else nu_ray(n)*rho(n)*chi(n, k)*chi(n, f)**4
    def mark(A):         # product-form mark d(A) = sum over tuples prod a_i(p_i) 1_{p_i | A}
        tot = 1.0 + 0j
        for lst in slots:
            s = 0j
            for p, c in lst:
                if divides(p, A): s += c
            tot *= s
        return tot
    sqf = [n for n in primaries(int(2*X) + 1, sqfree=True)]
    lhs = 0j
    for n in sqf:
        if N(n) >= 2*X: continue
        w = bumpW(N(n)/X)
        if w == 0: continue
        lhs += phi(n)**2*gam(2, n)*Psi(n)*mark(n)*w
    lhs *= X**-0.5
    # right side: sum over h of mu(h) conj(alpha(h))^hexp Psi(h)^3/N(h) T_{d_h}(X/N(h)^3)
    allp = primaries(int((2*X)**(1/3)) + 2)                # candidates for h and b (any primary, N^3 < 2X)
    rhs = 0j
    for h in allp:
        if muu(h) == 0 or N(h)**3 >= 2*X: continue
        Y = X/N(h)**3
        h3 = mul(mul(h, h), h)
        T = 0j
        for b in allp:
            if N(b)**3*N(h)**3 >= 2*X: continue
            b3 = mul(mul(b, b), b)
            for n in sqf:
                if N(n)*N(b)**3 >= 2*Y: continue
                v = Vstar(N(n)*N(b)**3/Y)
                if v == 0: continue
                A = mul(mul(h3, n), b3) if shift_mark else mul(n, b3)
                T += (phi(n)**2*gam(2, n)*Psi(n)*phi(b)**bexp*Psi(b)**3/(math.sqrt(N(n))*N(b))
                      *mark(A)*v)
        rhs += muu(h)*phi(h)**hexp*Psi(h)**3/N(h)*T
    return abs(lhs - rhs), abs(lhs)
P7 = [p for p in PRIMES if N(p) == 7]; P13 = [p for p in PRIMES if N(p) == 13]
p19, p31, p37, p43 = [[p for p in PRIMES if N(p) == q][0] for q in (19, 31, 37, 43)]
slotsA = [[(P7[0], 0.8 + 0.3j), (p19, -0.5 + 0.6j)], [(P13[0], 0.7 - 0.2j), (p31, 0.4 + 0.9j)]]
slotsB = []                                              # empty mark
cases = [(900.0, (5, 3), p43, p37, slotsB), (900.0, (2, 7), (1, 0), (1, 0), slotsB),
         (3000.0, (4, -3), p43, p37, slotsA)]
emax = 0.0; fmax3 = 0.0; fmaxm = 0.0
for X, k, f, rp, sl in cases:
    err, size = check_inversion(X, k, f, rp, sl)
    emax = max(emax, err/max(size, 1e-30))
    ef, _ = check_inversion(X, k, f, rp, sl, hexp=3)     # zeta power conj(alpha)^3 in the inversion
    fmax3 = max(fmax3, ef/max(size, 1e-30))
    if sl:
        em, _ = check_inversion(X, k, f, rp, sl, shift_mark=False)
        fmaxm = max(fmaxm, em/max(size, 1e-30))
    print(f"      X={X:g} k={k} f={f} mark={'yes' if sl else 'no'}: |LHS|={size:.3g}, rel err {err/max(size,1e-30):.2g}")
line("S4.5 (s4-inversion) CM marked cube inversion with conj(alpha)^6, nu, puncture, mark", emax, 1e-11, "relative")
fals("S4.5 conj(alpha(h))^3 (zeta power) in the inversion", fmax3)
fals("S4.5 unshifted mark d in place of d_h(A) = d(h^3 A)", fmaxm)

# ------------------------------------------------------------------ S4.6 overlap identity (initialization)
def check_overlap(eps, u, X, slots, twist_slots=True, eps_alpha=None):
    """Direct M_u Q_u (unnormalized) versus the overlap decomposition (s4-overlap)."""
    if eps_alpha is None: eps_alpha = eps
    nu = nu_ray
    def psi(n): return 1.0 + 0j if n == (1, 0) else nu(n)*chi(n, u)**eps if chi(n, u) != 0 else 0j
    def apow(n, e): return al(n)**e if e >= 0 else phi(n)**(-e)
    sqf = [n for n in primaries(int(2*X) + 1, sqfree=True)]
    W = lambda y: bumpW(y)
    Mu = sum(muu(n)*apow(n, eps_alpha)*psi(n)*W(N(n)/X) for n in sqf if N(n) < 2*X)
    Wi = [lambda y: 1.0 + 0.3*y, lambda y: 2.0 - 0.1*y]       # slot weights (any row-independent functions)
    Q = 1.0 + 0j
    for i, lst in enumerate(slots):
        Q *= sum(c*psi(p)*Wi[i](N(p)) for p, c in lst)
    direct = Mu*Q
    # decomposition: sum over the assigned set A, assigned tuple (j), surviving tuple with P0 | c
    I = list(range(len(slots))); dec = 0j
    for r_ in range(len(I) + 1):
        for A in itertools.combinations(I, r_):
            I0 = [i for i in I if i not in A]
            for assigned in itertools.product(*[slots[i] for i in A]):
                j = (1, 0); outer = 1.0 + 0j
                for i, (p, c) in zip(A, assigned):
                    j = mul(j, p); outer *= c*Wi[i](N(p))
                outer *= muu(j)*apow(j, eps_alpha)*psi(j)**2
                inner = 0j
                # c ranges over squarefree c prime to j with P0 | c; we enumerate it as c = n0 P0 (n0 squarefree),
                # which is a bijection onto that set, and evaluate the coefficient of c itself
                for surv in itertools.product(*[slots[i] for i in I0]):
                    P0 = (1, 0); coef = 1.0 + 0j
                    for i, (p, c) in zip(I0, surv):
                        P0 = mul(P0, p)
                        tw = apow(p, -eps_alpha) if twist_slots else 1.0
                        coef *= muu(p)*tw*c*Wi[i](N(p))
                    for n0 in sqf:
                        if N(j)*N(n0) >= 2*X: continue
                        c_ = mul(n0, P0)
                        if muu(c_) == 0 or not coprime(c_, j): continue
                        base = muu(c_)*apow(c_, eps_alpha)*psi(c_)
                        if base == 0: continue
                        inner += base*coef*W(N(j)*N(n0)/X)
                dec += outer*inner
    return abs(direct - dec), abs(direct)
slotsO = [[(P7[0], 0.8 + 0.3j), (P13[0], -0.5 + 0.6j)], [(P7[1], 0.7 - 0.2j), (p19, 0.4 + 0.9j)]]
eO = 0.0; fO = 0.0
for eps in (1, -1):
    for u in [(5, 3), (7, 1), mul(P7[0], (2, -1))]:
        err, size = check_overlap(eps, u, 300.0, slotsO)
        eO = max(eO, err/size)
        ef, _ = check_overlap(eps, u, 300.0, slotsO, twist_slots=False)
        fO = max(fO, ef/size)
line("S4.6 (s4-overlap) M_u Q_u = sum over assignments, slot coefficients mu(p) conj(alpha(p)^eps) a_i(p)",
     eO, 1e-12, "relative, eps = +1, -1, three rows (one sharing a slot prime)")
fals("S4.6 slot coefficients without conj(alpha(p)^eps)", fO)

# ------------------------------------------------------------------ S4.7 sixth-power column identity
def Mrow(v, D, eps):
    sqf = [n for n in primaries(int(2*D) + 1, sqfree=True)]
    tot = 0j
    for n in sqf:
        if N(n) >= 2*D: continue
        c = chi(n, v)
        if c == 0: continue
        tot += muu(n)*(al(n)**eps if eps > 0 else phi(n))*nu_ray(n)*c**eps*bumpW(N(n)/D)
    return tot*D**-0.5
e7 = 0.0; f7 = 0.0
for eps in (1, -1):
    for u, a_ in [((5, 3), mul(P7[0], P13[0])), ((2, 9), P7[1]), (mul(P7[0], (4, 3)), P7[0])]:
        D = 400.0
        a6 = (1, 0)
        for _ in range(6): a6 = mul(a6, a_)
        ua6 = mul(u, a6)
        lhs = Mrow(u, D, eps)
        rhs = 0j; rhs_f = 0j
        for d in divisors_sqfree(a_):
            psid = 1.0 + 0j if d == (1, 0) else (nu_ray(d)*chi(d, u)**eps if chi(d, u) != 0 else 0j)
            ad = 1.0 + 0j if d == (1, 0) else (al(d)**eps if eps > 0 else phi(d))
            term = muu(d)*psid*N(d)**-0.5*Mrow(ua6, D/N(d), eps)
            rhs += ad*term; rhs_f += term
        e7 = max(e7, abs(lhs - rhs)/abs(lhs)); f7 = max(f7, abs(lhs - rhs_f)/abs(lhs))
line("S4.7 (s4-sixth) M_u(D) = sum_{d | rad a} mu(d) alpha(d)^eps psi_u(d) N(d)^{-1/2} M_{ua^6}(D/N(d))", e7, 1e-12,
     "relative; includes u sharing a prime with a")
fals("S4.7 without alpha(d)^eps", f7)

print()
npass = sum(1 for _, ok in RESULTS if ok); print(f"SUMMARY: {npass}/{len(RESULTS)} as required"
                                                + ("" if npass == len(RESULTS) else "  <-- PROBLEM"))

# A printed failure must also fail the suite.
raise SystemExit(0 if (npass == len(RESULTS)) else 1)
