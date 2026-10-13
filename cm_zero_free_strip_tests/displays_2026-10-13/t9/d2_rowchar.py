"""d2: Lemma s5-rowchar and Steps 3 and 7 of the proof of Proposition s5-plain (CM guide, Appendix F).

F.5  eq:d175   alpha^e psi_k^*((a)) = phi_alpha(a)^e phi_0(a) (a/|a|)^e   for a prime to 3 g_k,
               phi_alpha(a) = u(a)^{-1} (a = u(a) a_0, a_0 primary), phi_0 primitive mod g_k, trivial on units.
               Also: phi_alpha is a character of (O/3)^x with phi_alpha(u) = u^{-1}; phi_k(u) = u^{-e};
               Xi_k = alpha^e psi_k^* on ideals prime to 3 g_k  (Lemma s5-rowchar (a)).
F.4  eq:s5-natural   [alpha^e psi_k^0](l) = Xi_k(l) 1_{(l, r_k) = 1}   for EVERY nonzero ideal l (N(l) <= 2500),
               r_k = primes not dividing f_k that lie in S = {lambda, 2} or at which psi_k^0 vanishes.
F.24 eq:d187   T_k(X;W) = sum_{d0 | r_k} mu(d0) Xi_k(d0) N(d0)^{-1/2} T_{Xi_k}(X/N(d0); W).
F.25 eq:d188   T_{conj Xi_k}(Y;W#) = sum_{rad h0 | r_k} conj Xi_k(h0) N(h0)^{-1/2} T_{conj(alpha^e psi_k^0)}(Y/N(h0); W#).
F.26 eq:d189   Delta_k = Q (S(X1)S(X2) - S(Y1)S(Y2)) = Q (X1X2)^{-1/2} sum_{l1,l2} alpha(l1l2)^e psi_k(l1l2) {W1W2 - W1W2}
               when X1 X2 = Y1 Y2, with OAI78's row character psi_k(n) = tau(n) chi_n(k) evaluated on the product
               l1 l2 (factored independently).
Rows: psi(n) = tau(n) prod_p chi_p(n)^{v_p} (residue symbols modulo the row primes, zero extension kept, v_p = 6 a
redundant zero), as in tests_cm/phase2_2026-10-12/s5/t5; this lies in OAI78's coefficient class (sextic reciprocity).
Characters are built exactly as finite tables on (O/m)^x and reduced to their conductors (hecke.py).
Falsifiers: OAI78's R_{0,k} (computed from g_k) in place of r_k; r_k without the primes of S; no indicator;
phi_alpha(a) = u(a); omitting mu(d0); squarefree h0 only; unequal products X1X2 != Y1Y2; no angular factor.
"""
import os, sys, math, itertools, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t9lib import (report, falsifier, summary, cx, mul, add, conj, N, OM, UNITS, chi_p, primes_upto, primary_assoc,
                   is_primary, divides, exact_div, from_primary_function, primitive, elements, LAM, TWO, alpha_el,
                   cubic_mod2, sextant_generator)

P = primes_upto(4000)
def prime_of_norm(n, idx=0): return [p for p in P if N(p) == n][idx]
p7, p13, p19, p31 = prime_of_norm(7), prime_of_norm(13), prime_of_norm(19), prime_of_norm(31)
one = lambda a: 1

def make_psi(tau, kfac):
    def psi(a):
        v = tau(a)
        for p, ex in kfac:
            c = chi_p(a, p)
            if c == 0: return 0
            v *= c ** ex
        return v
    return psi

rows = [
    ("chi_p7^5, e=+1 (guide example)", [(p7, 5)], one, 1),
    ("chi_p7^5, e=-1", [(p7, 5)], one, -1),
    ("chi_p7 chi_p13^6, e=+1", [(p7, 1), (p13, 6)], one, 1),
    ("chi_p7 chi_p13^6, e=-1", [(p7, 1), (p13, 6)], one, -1),
    ("(./2)_3 chi_p19^2 chi_p7^3, e=+1", [(p19, 2), (p7, 3)], cubic_mod2, 1),
    ("chi_p31^5, e=-1", [(p31, 5)], one, -1),
    ("(./2)_3 chi_p13^4, e=-1", [(p13, 4)], cubic_mod2, -1),
    ("chi_p13^3, e=+1", [(p13, 3)], one, 1),
]

def coprime_to(a, plist): return all(not divides(p, a) for p in plist)

# all nonzero ideals with N <= BI, one generator each
BI = 2500
RR = int(math.ceil(2 * math.sqrt(BI / 3.0))) + 2
ideals = {}
for a in range(-RR, RR + 1):
    for b in range(-RR, RR + 1):
        x = (a, b)
        if 0 < N(x) <= BI:
            g = sextant_generator(x); ideals[g] = 1
ideals = list(ideals)
print(f"{len(ideals)} nonzero ideals with N <= {BI}")

def xi_val(ch, a):
    z = cx(a); return ch(np.array([a[0]]), np.array([a[1]]))[0] * (z / abs(z)) ** ch.e

def T_xi(ch, X, Wf, nmax, conjugate=False):
    """X^{-1/2} sum over nonzero ideals a with N(a) <= nmax of Xi(a) W(N a/X)  (= 1/6 x element sum)."""
    if nmax < 1: return 0.0
    A, B, Z, n = elements(int(nmax) + 1)
    k = n <= nmax
    A, B, Z, n = A[k], B[k], Z[k], n[k]
    v = ch(A, B) * (Z / np.abs(Z)) ** ch.e
    if conjugate: v = np.conj(v)
    return np.sum(v * Wf(n / X)) / 6.0 / math.sqrt(X)

def bump(x):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x); m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out
def W12(x): return np.sqrt(x) * np.exp(-x)
def W12s(y): return np.sqrt(y) * np.exp(-y)          # its W^sharp (closed form, t1 1b)

worst = {k: 0.0 for k in ["F5", "F5-unit", "F5-phik", "F5-agree", "F4", "F24", "F25"]}
fal = {k: 0.0 for k in ["F5-noinv", "F4-R0", "F4-noS", "F4-noind", "F24-R0", "F24-nomu", "F25-sqfree", "F25-R0"]}
fal_used = {k: False for k in fal}

for name, kfac, tau, e in rows:
    kprimes = [p for p, ex in kfac]
    m = (6, 0)
    for p in kprimes: m = mul(m, p)
    mprimes = [LAM, LAM, TWO] + kprimes
    psi = make_psi(tau, kfac)
    def cmf(a, psi=psi, e=e): return alpha_el(a) ** e * psi(a)
    Xi = primitive(from_primary_function(cmf, m, mprimes, e))          # Xi_k, conductor f_k
    fin = primitive(from_primary_function(psi, m, mprimes, 0))          # psi_k^*, conductor g_k (phi_0 = its table)
    f, g = Xi.m, fin.m; fpr, gpr = set(Xi.primes), set(fin.primes)
    allp = [LAM, TWO] + kprimes                                          # S u {zeros of psi_k^0}
    r = [p for p in allp if p not in fpr]
    R0 = [p for p in allp if p not in gpr]
    rnoS = [p for p in kprimes if p not in fpr]
    print(f"--- {name}: N(f_k)={N(f)}, N(g_k)={N(g)}, r_k norms {[N(p) for p in r]}, R_0k norms {[N(p) for p in R0]}")
    # Lemma s5-rowchar (c): r_k squarefree (distinct primes), prime to f_k, and r_k, R_0k differ at most at lambda
    bad = (len(set(r)) != len(r)) or any(p in fpr for p in r) or not (set(r) ^ set(R0)) <= {LAM}
    report(f"F4-rk-{name}", float(bad), 0.0, "r_k squarefree, prime to f_k, differs from R_0k at most at lambda")

    # ---------------- F.5 (eq:d175) and Lemma (a)
    for a in [(x, y) for x in range(-14, 15) for y in range(-14, 15)]:
        if N(a) == 0 or not coprime_to(a, [LAM] + list(gpr)): continue
        u, a0 = primary_assoc(a)
        phi_alpha = 1 / cx(u)
        z = cx(a); ang = z / abs(z)
        lhs = alpha_el(a0) ** e * fin(np.array([a0[0]]), np.array([a0[1]]))[0]       # alpha^e psi^* on the ideal (a)
        rhs = phi_alpha ** e * fin(np.array([a[0]]), np.array([a[1]]))[0] * ang ** e
        worst["F5"] = max(worst["F5"], abs(lhs - rhs))
        fal["F5-noinv"] = max(fal["F5-noinv"], abs(lhs - cx(u) ** e * fin(np.array([a[0]]), np.array([a[1]]))[0] * ang ** e))
        fal_used["F5-noinv"] = True
        # phi_alpha depends only on a mod 3 and is multiplicative
        a3 = add(a, mul((3, 0), (x_ := (a[0] % 5) - 2, 1)))
        if N(a3) and coprime_to(a3, [LAM]):
            worst["F5-unit"] = max(worst["F5-unit"], abs(1 / cx(primary_assoc(a3)[0]) - phi_alpha))
        # Xi_k agrees with alpha^e psi_k^* on ideals prime to 3 g_k
        worst["F5-agree"] = max(worst["F5-agree"], abs(xi_val(Xi, a) - lhs))
    for u in UNITS:
        worst["F5-unit"] = max(worst["F5-unit"], abs(1 / cx(primary_assoc(u)[0]) - 1 / cx(u)))
        worst["F5-phik"] = max(worst["F5-phik"], abs(Xi(np.array([u[0]]), np.array([u[1]]))[0] - cx(u) ** (-e)))
        worst["F5-phik"] = max(worst["F5-phik"], abs(fin(np.array([u[0]]), np.array([u[1]]))[0] - 1))

    # ---------------- F.4 (eq:s5-natural) for every ideal with N <= BI
    e4 = 0; fr0 = 0; fns = 0; fni = 0
    for a in ideals:
        if coprime_to(a, [LAM, TWO]):
            u, a0 = primary_assoc(a); lhs = cmf(a0)
        else:
            lhs = 0
        xv = xi_val(Xi, a)
        e4 = max(e4, abs(lhs - xv * coprime_to(a, r)))
        fr0 = max(fr0, abs(lhs - xv * coprime_to(a, R0)))
        fns = max(fns, abs(lhs - xv * coprime_to(a, rnoS)))
        fni = max(fni, abs(lhs - xv))
    worst["F4"] = max(worst["F4"], e4)
    if set(r) != set(R0):
        fal["F4-R0"] = max(fal["F4-R0"], fr0); fal_used["F4-R0"] = True
    fal["F4-noS"] = max(fal["F4-noS"], fns); fal_used["F4-noS"] = True
    fal["F4-noind"] = max(fal["F4-noind"], fni); fal_used["F4-noind"] = True
    report(f"F4-{name}", e4, 1e-12, f"[alpha^e psi^0](l) = Xi_k(l) 1_(l,r_k)=1 on all {len(ideals)} ideals N <= {BI}")

    # ---------------- F.24 (eq:d187)
    def Tk(X, Wf, nmax):
        tot = 0
        A, B, Z, n = elements(int(nmax) + 1)
        for a_, b_, nn in zip(A, B, n):
            a = (int(a_), int(b_))
            if nn > nmax or not is_primary(a) or not coprime_to(a, [LAM, TWO]): continue
            v = cmf(a)
            if v != 0: tot += v * Wf(np.array([nn / X]))[0]
        return tot / math.sqrt(X)
    def rhs187(X, Wf, nmax, rlist, use_mu=True):
        tot = 0
        for k_ in range(len(rlist) + 1):
            for sub in itertools.combinations(rlist, k_):
                d0 = (1, 0)
                for p in sub: d0 = mul(d0, p)
                nd = N(d0)
                if nd > nmax: continue
                c = ((-1) ** k_ if use_mu else 1) * xi_val(Xi, d0) / math.sqrt(nd)
                tot += c * T_xi(Xi, X / nd, Wf, nmax / nd)
        return tot
    e24 = 0
    for X, Wf, cut in [(3.0, bump, 2.0), (40.0, bump, 2.0), (300.0, bump, 2.0), (0.7, W12, 60.0), (5.0, W12, 60.0)]:
        lhs = Tk(X, Wf, cut * X)
        e24 = max(e24, abs(lhs - rhs187(X, Wf, cut * X, r)))
        if set(r) != set(R0):
            fal["F24-R0"] = max(fal["F24-R0"], abs(lhs - rhs187(X, Wf, cut * X, R0))); fal_used["F24-R0"] = True
        fal["F24-nomu"] = max(fal["F24-nomu"], abs(lhs - rhs187(X, Wf, cut * X, r, use_mu=False))); fal_used["F24-nomu"] = True
    worst["F24"] = max(worst["F24"], e24)
    report(f"F24-{name}", e24, 1e-12, "T_k(X;W) = sum_{d0|r_k} mu Xi_k(d0) N(d0)^{-1/2} T_Xi(X/N d0; W), bump X=3,40,300 and x^{1/2}e^{-x} X=0.7,5")

    # ---------------- F.25 (eq:d188)
    def T_conj_nat(Y, Wf, nmax):
        tot = 0
        A, B, Z, n = elements(int(nmax) + 1)
        for a_, b_, nn in zip(A, B, n):
            a = (int(a_), int(b_))
            if nn > nmax or not is_primary(a) or not coprime_to(a, [LAM, TWO]): continue
            v = cmf(a)
            if v != 0: tot += np.conj(v) * Wf(np.array([nn / Y]))[0]
        return tot / math.sqrt(Y)
    def h0_list(rlist, nmax, sqfree=False):
        out = [(1, 0)]
        for p in rlist:
            new = []
            for h in out:
                hh = h
                while N(hh) <= nmax:
                    new.append(hh)
                    if sqfree and hh != h: break
                    hh = mul(hh, p)
            out = new
        return out
    def rhs188(Y, Wf, nmax, rlist, sqfree=False):
        tot = 0
        for h0 in h0_list(rlist, nmax, sqfree):
            nh = N(h0)
            tot += np.conj(xi_val(Xi, h0)) / math.sqrt(nh) * T_conj_nat(Y / nh, Wf, nmax / nh)
        return tot
    e25 = 0
    for Y in [0.5, 3.0, 12.0]:
        nmax = 60 * Y
        lhs = T_xi(Xi, Y, W12s, nmax, conjugate=True)
        e25 = max(e25, abs(lhs - rhs188(Y, W12s, nmax, r)))
        fal["F25-sqfree"] = max(fal["F25-sqfree"], abs(lhs - rhs188(Y, W12s, nmax, r, sqfree=True))); fal_used["F25-sqfree"] = True
        if set(r) != set(R0):
            fal["F25-R0"] = max(fal["F25-R0"], abs(lhs - rhs188(Y, W12s, nmax, R0))); fal_used["F25-R0"] = True
    worst["F25"] = max(worst["F25"], e25)
    report(f"F25-{name}", e25, 1e-12, "T_conjXi(Y;W#) = sum_{rad h0|r_k} conj Xi(h0) N(h0)^{-1/2} T_conj(alpha^e psi^0)(Y/N h0; W#), Y=0.5,3,12")

report("F5", worst["F5"], 1e-12, "alpha^e psi^*((a)) = phi_alpha(a)^e phi_0(a) (a/|a|)^e for all a prime to 3g_k in a box, all rows")
report("F5-phialpha", worst["F5-unit"], 1e-12, "phi_alpha(a) = u(a)^{-1} depends only on a mod 3; phi_alpha(u) = u^{-1}")
report("F5-units", worst["F5-phik"], 1e-12, "phi_k(u) = u^{-e} and phi_0(u) = 1 on the six units")
report("F5-agree", worst["F5-agree"], 1e-12, "Xi_k = alpha^e psi_k^* on ideals prime to 3 g_k (Lemma s5-rowchar (a))")
falsifier("F5", "phi_alpha(a) = u(a) (no inverse)", fal["F5-noinv"], 1e-3)
for k, lab in [("F4-R0", "OAI78's R_0k (from g_k) in place of r_k"), ("F4-noS", "r_k without the primes of S"),
               ("F4-noind", "no indicator 1_(l,r_k)=1"), ("F24-R0", "R_0k in place of r_k"),
               ("F24-nomu", "omit mu(d0)"), ("F25-sqfree", "squarefree h0 only"), ("F25-R0", "R_0k in place of r_k")]:
    assert fal_used[k], k
    falsifier(k.split("-")[0], lab, fal[k], 1e-3)

# ---------------- F.26 (eq:d189) with OAI78's row character psi_k(n) = tau(n) chi_n(k)
PP = [p for p in P if N(p) <= 40000]
def factor_primary(n):
    """Factor a primary n coprime to 6 into primary primes (trial division)."""
    out = []; rem = n
    for p in PP:
        if N(rem) == 1: break
        if N(p) > N(rem): continue
        while N(rem) % N(p) == 0 and divides(p, rem):
            out.append(p); rem = exact_div(rem, p)
    assert N(rem) == 1, (n, rem)
    return out
def psi_oai(n, k):
    v = cubic_mod2(n)
    for p in factor_primary(n):
        v *= chi_p(k, p)
    return v

prim = []
for a in range(-30, 31):
    for b in range(-30, 31):
        x = (a, b)
        if 0 < N(x) <= 220 and is_primary(x) and coprime_to(x, [LAM, TWO]): prim.append(x)
def bumpab(x, a, b):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x); mm = (x > a) & (x < b); u = (2 * x[mm] - a - b) / (b - a)
    out[mm] = np.exp(-1.0 / (1 - u * u)); return out
W1 = lambda x: bumpab(x, 1.0, 2.0)
W2 = lambda x: bumpab(x, 0.8, 2.5) * (1 + 0.3 * np.asarray(x, dtype=float))
rng = random.Random(5)
err = 0; f_uneq = 0; f_noalpha = 0
for k in [mul((-2, 0), p7), mul(p13, (2, 5)), (11, -4)]:
    for e in (1, -1):
        X1, X2, Y1 = 40.0, 70.0, 56.0; Y2 = X1 * X2 / Y1
        Q = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
        vals = {x: alpha_el(x) ** e * psi_oai(x, k) for x in prim}
        def S(X, Wf):
            return sum(v * Wf(np.array([N(x) / X]))[0] for x, v in vals.items()) / math.sqrt(X)
        lhs = Q * (S(X1, W1) * S(X2, W2) - S(Y1, W1) * S(Y2, W2))
        Y2b = 1.1 * Y2
        lhs_b = Q * (S(X1, W1) * S(X2, W2) - S(Y1, W1) * S(Y2b, W2))
        tot = 0; totb = 0; tot0 = 0
        for l1 in prim:
            for l2 in prim:
                n1, n2 = N(l1), N(l2)
                br = W1(np.array([n1 / X1]))[0] * W2(np.array([n2 / X2]))[0] - W1(np.array([n1 / Y1]))[0] * W2(np.array([n2 / Y2]))[0]
                brb = W1(np.array([n1 / X1]))[0] * W2(np.array([n2 / X2]))[0] - W1(np.array([n1 / Y1]))[0] * W2(np.array([n2 / Y2b]))[0]
                if br == 0 and brb == 0: continue
                l = mul(l1, l2)                                   # the product, factored independently
                c = psi_oai(l, k)
                tot += alpha_el(l) ** e * c * br; totb += alpha_el(l) ** e * c * brb; tot0 += c * br
        rhs = Q / math.sqrt(X1 * X2) * tot
        err = max(err, abs(lhs - rhs))
        f_uneq = max(f_uneq, abs(lhs_b - Q / math.sqrt(X1 * X2) * totb))
        f_noalpha = max(f_noalpha, abs(lhs - Q / math.sqrt(X1 * X2) * tot0))
report("F26", err, 1e-12, "Delta_k = Q (X1X2)^{-1/2} sum alpha(l1l2)^e psi_k(l1l2){...}, psi_k(n) = (n/2)_3 chi_n(k), 3 rows, e = +-1")
falsifier("F26", "unequal products X1X2 != Y1Y2", f_uneq, 1e-4)
falsifier("F26", "no angular factor", f_noalpha, 1e-4)
summary()
