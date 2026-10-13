"""Part D, test 01: the two eta quotients of Section 'Two examples' (lines 1531-1554) and the Part I rows
(lines 1625-1626), plus the gamma factor of the functional equation (line 1605).

Sources of truth used here (independent of the paper):
  * PARI/GP 2.17 (mfinit, mffrometaquo, mfcoefs, mfisCM, mfconductor, lfuncheckfeq);
  * an independent Python computation of the q-expansions from Euler's pentagonal and Jacobi's identities;
  * exact Z[omega] arithmetic (zw.py) for the Hecke-character sides.
Each claim has a falsifier that must fail.  Runtime: about 1 minute.
"""
import os, sys, math, cmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import *

NBIG = 40000      # PARI coefficients (for a_{p^2}, p <= 199, and exact Hecke-character identities)
N0 = 3000         # the paper's range

# ---------------------------------------------------------------- PARI data
out = gp(f"""
F=mffrometaquo([3,8]); mf9=mfinit([9,4],0);
print(mfdim(mf9)); print(mftobasis(mf9,F)); print(#mfeigenbasis(mf9)); print(mfcoefs(mfeigenbasis(mf9)[1],50)==mfcoefs(F,50));
print(mfisCM(F)); print(mfconductor(mf9,F)); print(mfparams(F)[1..3]);
G=mffrometaquo([2,3;6,3]); mf12=mfinit([12,3,-3],0);
print(mfdim(mf12)); print(mftobasis(mf12,G)); print(#mfeigenbasis(mf12)); print(mfcoefs(mfeigenbasis(mf12)[1],50)==mfcoefs(G,50));
print(mfisCM(G)); print(mfconductor(mf12,G)); print(mfparams(G)[1..3]);
print(mftobasis(mfinit([12,3],1),G,1)); print(mftobasis(mfinit([3,4],1),F,1));
an=Vec(mfcoefs(F,4000))[2..4001]; bn=Vec(mfcoefs(G,4000))[2..4001];
print(lfuncheckfeq(lfuncreate([an,0,[0,1],4,9,1])));
print(lfuncheckfeq(lfuncreate([an,0,[0,1],4,27,1])));
print(lfuncheckfeq(lfuncreate([bn,0,[0,1],3,12,1])));
print(lfuncheckfeq(lfuncreate([bn,0,[0,1],3,36,1])));
print(Vec(mfcoefs(F,{NBIG})));
print(Vec(mfcoefs(G,{NBIG})));
""", timeout=900).strip().split("\n")
(d9, tb9, ne9, eig9, cm9, cond9, par9, d12, tb12, ne12, eig12, cm12, cond12, par12,
 wrongG, wrongF, feqF, feqFw, feqG, feqGw, aF, aG) = out
aF = [int(x) for x in aF.strip("[]").split(",")]
aG = [int(x) for x in aG.strip("[]").split(",")]
assert len(aF) == NBIG + 1 and len(aG) == NBIG + 1

# ---------------------------------------------------------------- independent q-expansions in Python
def pent(Nn):           # prod (1-q^n) = sum_k (-1)^k q^{k(3k-1)/2}
    a = [0]*(Nn + 1); k = 0
    while True:
        done = True
        for kk in (k, -k) if k else (0,):
            e = kk*(3*kk - 1)//2
            if e <= Nn:
                a[e] += (-1)**kk; done = False
        if done and k > 0:
            break
        k += 1
    return a
def jac3(Nn):           # prod (1-q^n)^3 = sum_k (-1)^k (2k+1) q^{k(k+1)/2}
    a = [0]*(Nn + 1); k = 0
    while k*(k + 1)//2 <= Nn:
        a[k*(k + 1)//2] += (-1)**k*(2*k + 1); k += 1
    return a
def pmul(a, b, Nn):
    c = [0]*(Nn + 1)
    nz = [i for i in range(Nn + 1) if a[i]]
    for i in nz:
        ai = a[i]
        for j in range(Nn + 1 - i):
            if b[j]:
                c[i + j] += ai*b[j]
    return c
def subs(a, m, Nn):     # a(q) -> a(q^m)
    c = [0]*(Nn + 1)
    for i in range(Nn//m + 1):
        c[i*m] = a[i]
    return c
M = N0
P1, P3 = pent(M), jac3(M)
P8 = pmul(pmul(P3, P3, M), pmul(P1, P1, M), M)          # prod(1-q^n)^8
f4 = [0] + subs(P8, 3, M)[:M]                           # eta(3z)^8 = q prod(1-q^{3n})^8
P7 = pmul(pmul(P3, P3, M), P1, M)
f4_wrong = [0] + subs(P7, 3, M)[:M]                     # eta(3z)^7 (falsifier)
g3 = [0] + pmul(subs(P3, 2, M), subs(P3, 6, M), M)[:M]  # eta(2z)^3 eta(6z)^3 = q prod(1-q^{2n})^3(1-q^{6n})^3

report("D01.1", 1532, "eta products", "q-expansions of eta(3z)^8 and eta(2z)^3eta(6z)^3 to q^3000: PARI = Euler/Jacobi product",
       aF[:M + 1] == f4 and aG[:M + 1] == g3 and aF[:M + 1] != f4_wrong,
       "falsifier eta(3z)^7 differs")
report("D01.2", 1537, "eta(3z)^8", "weight 4, level 9 newform with CM by Q(sqrt-3)",
       d9 == "1" and tb9 == "[1]~" and eig9 == "1" and cm9 == "-3" and cond9 == "9" and par9 == "[9, 4, 1]"
       and wrongF == "[]~",
       f"mfdim S_4^new(9)={d9}, coords={tb9}, CM disc={cm9}, conductor={cond9}; falsifier: not in S_4(Gamma0(3)) ({wrongF})")
report("D01.3", 1543, "eta(2z)^3eta(6z)^3", "weight 3, level 12 newform, nebentypus chi_{-3}, CM by Q(sqrt-3)",
       d12 == "1" and tb12 == "[1]~" and eig12 == "1" and cm12 == "-3" and cond12 == "12" and par12 == "[12, 3, -3]"
       and wrongG == "[]~",
       f"mfparams={par12}, CM disc={cm12}, conductor={cond12}; falsifier: not in S_3(Gamma0(12)) trivial char ({wrongG})")

# ---------------------------------------------------------------- split primes and their generators
P = rational_primes(NBIG)
split = [p for p in P if p % 3 == 1 and p < N0]
report("D01.4", 1537, "count", "number of split rational primes p < 3000 is 207 (414 primary generators)",
       len(split) == 207 and 2*len(split) == 414 and len([p for p in P if p % 3 == 1 and p < 2900]) != 207,
       f"count={len(split)}; falsifier bound 2900 gives {len([p for p in P if p % 3 == 1 and p < 2900])}")

def cxp(x): return complex(cx(x))
bad = bad_f = 0
for p in split:
    for pi in split_primary_pair(p):
        z = cxp(pi)
        if abs(2*(z**3).real - aF[p]) > 1e-6: bad += 1
        if abs(2*((-z)**3).real - aF[p]) > 1e-6: bad_f += 1         # non-primary associate -pi
report("D01.5", 1539, "a_p(eta(3z)^8)", "a_p = 2 Re(pi^3) for both primary generators of all 207 split p < 3000",
       bad == 0 and bad_f == 414, f"mismatches={bad}/414; falsifier -pi (non-primary) mismatches {bad_f}/414")

inert = [p for p in P if p % 3 == 2]
ok1 = all(aF[p] == 0 for p in inert if p <= NBIG)
sq = [p for p in inert if p*p <= NBIG]
ok2 = all(aF[p*p] == -p**3 for p in sq)
fal = all(aF[p*p] == p**3 for p in sq)
report("D01.6", 1539, "inert, eta(3z)^8", "p = 2 mod 3: a_p = 0 (all p <= 40000) and a_{p^2} = -p^3",
       ok1 and ok2 and not fal, f"a_{{p^2}} checked for {len(sq)} primes p <= {sq[-1]}; falsifier +p^3 fails")

# exact identity a_m = sum_{g primary, N(g)=m} g^3  (L(f,s+3/2) = L(s,alpha^3), conductor lambda)
E = eis_elements_upto(NBIG)
def exact_sum(m, fn, cond):
    s = (0, 0)
    for g in E.get(m, []):
        if cond(g):
            s = add(s, fn(g))
    return s
badE = 0; badF = 0
for m in range(1, NBIG + 1):
    if m % 3 == 0:
        tgt = (0, 0)
    else:
        tgt = exact_sum(m, lambda g: pw(g, 3), primary)
    if tgt != (aF[m], 0): badE += 1
    if m <= N0 and m % 3:
        # falsifier: twist alpha^3 by (2/g)_3  (wrong character)
        s = 0
        for g in E.get(m, []):
            if primary(g) and not divides((2, 0), g):
                s += complex(cx(pw(g, 3)))*complex(cx(mod2_class(g)))
            elif primary(g):
                s += 0
        if abs(s - aF[m]) > 1e-6: badF += 1
report("D01.7", 1540, "L(f,s+3/2)=L(s,alpha^3)", "a_m = sum_{g primary, N g = m} g^3 exactly (0 if 3|m), all m <= 40000",
       badE == 0 and badF > 100, f"exact mismatches={badE}; falsifier alpha^3*(2/.)_3 mismatches {badF} of m<=3000")

# ---------------------------------------------------------------- eta(2z)^3 eta(6z)^3
def chi_prime(pi, x):      # chi_pi(x) = (x/pi)_6 for a primary prime pi coprime to 6 (Euler's criterion)
    j = sym6_prime(x, pi)
    return 0 if j is None else cmath.exp(1j*math.pi*j/3)
bad = bad_f = bad_sym = nontriv = 0
for p in split:
    for pi in split_primary_pair(p):
        z = cxp(pi)
        c4 = chi_prime(pi, (4, 0)); c2_3 = cubic_sym_prime((2, 0), pi)
        if abs(c4 - c2_3) > 1e-9: bad_sym += 1
        if abs(c2_3 - 1) > 1e-9: nontriv += 1
        if abs(2*(z*z*c2_3).real - aG[p]) > 1e-6: bad += 1
        if abs(2*(z*z*c2_3.conjugate()).real - aG[p]) > 1e-6: bad_f += 1
report("D01.8", 1545, "a_p(eta(2z)^3eta(6z)^3)", "a_p = 2Re(pi^2 (2/pi)_3) = 2Re(pi^2 chi_pi(4)), both generators, 207 split p",
       bad == 0 and bad_sym == 0 and bad_f == nontriv > 0,
       f"mismatches={bad}/414, (2/pi)_3 != chi_pi(4) at {bad_sym}; falsifier conj((2/pi)_3) mismatches {bad_f}/414"
       f" (= all {nontriv} generators with (2/pi)_3 != 1)")

ok1 = all(aG[p] == 0 for p in inert if p != 2)
sq2 = [p for p in inert if p != 2 and p*p <= NBIG]
ok2 = all(aG[p*p] == p*p for p in sq2)
fal = all(aG[p*p] == -p*p for p in sq2)
report("D01.9", 1545, "inert, eta(2z)^3eta(6z)^3", "p nmid 12, p = 2 mod 3: a_p = 0 and a_{p^2} = p^2 (nebentypus chi_{-3})",
       ok1 and ok2 and not fal, f"{len(sq2)} primes 5 <= p <= {sq2[-1]}; falsifier -p^2 fails")
report("D01.10", 1547, "a_4", "a_4 = 0, not 4 (and a_2 = 0)", aG[4] == 0 and aG[2] == 0 and aG[4] != 4,
       f"a_2={aG[2]}, a_4={aG[4]}; falsifier value 4 rejected")
pw3 = [3**k for k in range(1, 12) if 3**k <= NBIG]
pw2 = [2**k for k in range(1, 20) if 2**k <= NBIG]
ok = aG[3] == -3 and all(aG[q] == (-3)**k for k, q in enumerate(pw3, 1)) and all(aG[q] == 0 for q in pw2)
fal = all(aG[q] == q for q in pw3)
report("D01.11", 1549, "Euler factors at 2,3", "a_3 = -3, a_{3^k} = (-3)^k (factor (1+3^{-s})^{-1}), a_{2^k} = 0 (factor 1)",
       ok and not fal, f"a_3={aG[3]}, 3^k<= {pw3[-1]}, 2^k<= {pw2[-1]}; falsifier (1-3^-s)^-1 fails")

# exact identity for the primitive character: a_m = (1/6) sum_{g: N g = m, (g,2)=1} g^2 eps(g), eps(g) = g mod 2
badE = 0; badF = 0
for m in range(1, NBIG + 1):
    s = (0, 0); sf = 0
    for g in E.get(m, []):
        if divides((2, 0), g):
            continue
        u = mod2_class(g)
        s = add(s, mul(pw(g, 2), u))
        if m <= N0:
            sf += complex(cx(pw(g, 2)))*complex(cx(u)).conjugate()
    if s[1] != 0 or s[0] % 6 or s[0]//6 != aG[m]: badE += 1
    if m <= N0 and abs(sf/6 - aG[m]) > 1e-6: badF += 1
report("D01.12", 1547, "L(f,s+1)=L(s,alpha^2 chi(4))", "primitive character (conductor 2, value -1 at lambda): a_m = (1/6) sum_{(g,2)=1, Ng=m} g^2 eps(g) exactly, m <= 40000",
       badE == 0 and badF > 100, f"exact mismatches={badE}; falsifier conj(eps) mismatches {badF} of m<=3000")

# ---------------------------------------------------------------- Dirichlet series with sextic symbols (L^(6) and Part I rows)
from fastgauss import primary_primes
PR = [pi for pi in primary_primes(N0) if pi != (-2, 0)]   # primary primes coprime to 6, norm <= N0
def euler_series(vals, Nn):
    """Dirichlet coefficients of prod_pi (1 - v_pi N(pi)^{-s})^{-1} up to Nn, vals = list of (N(pi), v_pi)."""
    c = [0j]*(Nn + 1); c[1] = 1
    for q, v in vals:
        if q > Nn: continue
        # multiply by sum_k v^k q^{-ks}
        new = c[:]
        qk, vk = q, v
        while qk <= Nn:
            for n in range(1, Nn//qk + 1):
                if c[n]: new[n*qk] += vk*c[n]
            qk *= q; vk *= v
        c = new
    return c
def inverse_series(vals, Nn):
    """Coefficients of prod (1 - v N^{-s})."""
    c = [0j]*(Nn + 1); c[1] = 1
    for q, v in vals:
        if q > Nn: continue
        new = c[:]
        for n in range(1, Nn//q + 1):
            if c[n]: new[n*q] -= v*c[n]
        c = new
    return c
def al(x): z = complex(cx(x)); return z/abs(z)
vals2 = [(N(pi), al(pi)**2*chi_prime(pi, (4, 0))) for pi in PR]
L6 = euler_series(vals2, N0)
err = 0; errF = 0
for n in range(1, N0 + 1):
    lhs = aG[n]/n
    rhs = 0; k = 0; m = n
    while True:
        rhs += (-1)**k*L6[m]
        if m % 3: break
        m //= 3; k += 1
    err = max(err, abs(lhs - rhs))
    errF = max(errF, abs(lhs - L6[n]))
report("D01.13", 1552, "L(f,s+1)=(1+3^-s)^-1 L^(6)", "Dirichlet coefficients agree for n <= 3000 (chi_pi(4) by Euler's criterion)",
       err < 1e-9 and errF > 0.5, f"max err={err:.1e}; falsifier without (1+3^-s)^-1: err={errF:.2f}")

def part1_rows(k, u0, Nn):
    vals = [(N(pi), al(pi)**k*chi_prime(pi, u0)) for pi in PR]
    return inverse_series(vals, Nn)
def recip_coeffs(a, k, Nn):
    """Dirichlet inverse of b(n) = a_n/n^{k/2}, restricted to the multiplicative part coprime to 6."""
    b = [0.0]*(Nn + 1)
    for n in range(1, Nn + 1):
        if math.gcd(n, 6) == 1: b[n] = a[n]/n**(k/2)
    inv = [0.0]*(Nn + 1); inv[1] = 1.0
    for n in range(2, Nn + 1):
        if math.gcd(n, 6) != 1: continue
        s = 0.0
        for d in range(2, n + 1):
            if n % d == 0 and b[d]: s += b[d]*inv[n//d]
        inv[n] = -s
    return inv
NR = 1500
for name, a, k, u0, u0bad in [("eta(3z)^8", aF, 3, (1, 0), (4, 0)), ("eta(2z)^3eta(6z)^3", aG, 2, (4, 0), (2, 0))]:
    inv = recip_coeffs(a, k, NR)
    rows = part1_rows(k, u0, NR); rowsF = part1_rows(k, u0bad, NR)
    e = max(abs(inv[n] - rows[n]) for n in range(1, NR + 1) if math.gcd(n, 6) == 1)
    eF = max(abs(inv[n] - rowsF[n]) for n in range(1, NR + 1) if math.gcd(n, 6) == 1)
    report(f"D01.14{'a' if k == 3 else 'b'}", 1626, f"Part I row, {name}",
           f"1/L(f,s+k/2) = sum mu_K(n) alpha(n)^{k} chi_n(u0) (coprime to 6), u0 = {u0[0]}, n <= {NR}",
           e < 1e-9 and eF > 0.1, f"max err={e:.1e}; falsifier u0={u0bad[0]}: err={eF:.2f}")

report("D01.15", 1605, "Gamma(s+k/2)", "PARI lfuncheckfeq with gamma factor Gamma_C(s) (= Gamma(s+k/2) after the shift s -> s+k/2), conductor 9 / 12",
       int(feqF) < -100 and int(feqG) < -100 and int(feqFw) > -10 and int(feqGw) > -10,
       f"log2 FE error: eta(3z)^8 {feqF}, eta(2z)^3eta(6z)^3 {feqG}; falsifier conductors 27/36: {feqFw}, {feqGw}")
