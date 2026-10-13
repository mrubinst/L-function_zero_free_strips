"""Independent checks (r3) of the Part II identities that feed Prop. cm-rows and Section 14, and of the
identification L(E_d,s+1/2) = L(s,psi_d^*) up to Euler factors, with falsifiers.
  G1  OAI78 eq:signal: gamma_1 gamma_2 = mu alpha G, gamma_2^3 = mu alpha     (squarefree c)
  G2  guide eq:signal-cm: conj(alpha) gamma_2 gamma_1 = mu G                  (falsifier: alpha)
  G3  companion: mu gamma_{-1} = chi_n(-1) G^{-1} conj(alpha) gamma_2         (falsifier: alpha)
  G4  CM weight: conj(alpha)^2 gamma_2 eta gamma_1 = mu G psi, psi = conj(alpha) eta
  G5  chi_p(4)^3 = 1 at good primes, so b_p^3 = a_p^2                          (falsifier: chi_p(2)^3)
  C1  sqrt(p)(psi*(pi)+psi*(pibar)) = a_p(E_d), split p not dividing 6d       (falsifier: alpha chi(4d))
  C2  psi*(-p) = -1 and a_p = 0 for inert p not dividing 6d                    (falsifier: no angular factor)
  C3  psi*(pibar) = conj(psi*(pi))
  C4  n -> chi_n(4d) is periodic modulo M = 144|d| on primary n prime to 6d (finite-order ray character)
  C5  bad primes: the Euler factor of L(E_d,s+1/2) at p | 6d is 1 or 1/(1+p^{-2s}) (no zeros, poles on Re s=0)
"""
import random, math, cmath, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r3lib import *

random.seed(3)
ok = True
_FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']   # self-test (added 2026-10-12): first check counts as failed
def rep(name, cond, info=''):
    global ok
    if _FORCE[0]:
        _FORCE[0] = False; cond = False; name = '[CM_TESTS_FORCE_FAIL=1] ' + name
    print(('[PASS] ' if cond else '[FAIL] ') + name + ('  :: ' + info if info else ''))
    ok &= bool(cond)

P = primary_primes(400)
split = [p for p in P if p[1] != 0]
inert = [p for p in P if p[1] == 0]
# squarefree c: single primes and products of two or three distinct non-conjugate primes
cs = [[p] for p in P[:45]]
while len(cs) < 75:
    k = random.choice([2, 2, 3])
    pr = random.sample(split[:30] + inert[:2], k)
    if len(set(pr)) < k or any(conj(a) == b for a in pr for b in pr):
        continue
    if norm(prod(pr)) > 15000:
        continue
    cs.append(pr)

e1 = e2 = e3 = e4 = 0; f2 = f3 = 0; e0 = 0
u0 = (2, 1)   # eta = chi_.(2+w), a finite-order character, as a sample finite part
for pr in cs:
    c = prod(pr)
    ga = gauss_all(pr)
    g1, g2, gm1 = ga[1], ga[2], ga[-1]
    m, a = mu(pr), alpha(c)
    Gc = chi(pr, (4, 0)).conjugate() * ga[3]
    e0 = max(e0, abs(g2 ** 3 - m * a))
    e1 = max(e1, abs(g1 * g2 - m * a * Gc))
    e2 = max(e2, abs(a.conjugate() * g2 * g1 - m * Gc))
    f2 = max(f2, abs(a * g2 * g1 - m * Gc))
    chim1 = chi(pr, (-1, 0))
    e3 = max(e3, abs(m * gm1 - chim1 / Gc * a.conjugate() * g2))
    f3 = max(f3, abs(m * gm1 - chim1 / Gc * a * g2))
    eta = chi(pr, u0)
    if eta != 0:
        psi = a.conjugate() * eta
        e4 = max(e4, abs(a.conjugate() ** 2 * g2 * eta * g1 - m * Gc * psi))
rep('G1 eq:signal (OAI78 L965) on %d squarefree c' % len(cs), max(e0, e1) < 1e-11, 'max err %.1e, %.1e' % (e0, e1))
rep('G2 eq:signal-cm', e2 < 1e-11 and f2 > 0.1, 'err %.1e; falsifier alpha: %.2f' % (e2, f2))
rep('G3 companion mu gamma_{-1} = chi(-1) G^{-1} conj(alpha) gamma_2', e3 < 1e-11 and f3 > 0.1,
    'err %.1e; falsifier alpha: %.2f' % (e3, f3))
rep('G4 CM weight conj(alpha)^2 gamma_2 eta gamma_1 = mu G psi', e4 < 1e-11, 'err %.1e' % e4)
bad4 = [p for p in P if chi_prime(p, (4, 0)) not in [0] and cx(chi_prime(p, (4, 0))) ** 3 != 1]
bad4 = [p for p in P if abs(cx(chi_prime(p, (4, 0))) ** 3 - 1) > 1e-12]
bad2 = [p for p in P if abs(cx(chi_prime(p, (2, 0))) ** 3 - 1) > 1e-12]
rep('G5 chi_p(4)^3 = 1 at %d good primes' % len(P), not bad4 and len(bad2) > 0,
    'falsifier chi_p(2)^3 = 1 fails at %d primes' % len(bad2))

# ---- curves -----------------------------------------------------------------------------
ds = [1, 2, 3, -3, 5, 7, -11, 16, -432, 30, -6, 13]
B = 1500
script = "default(parisize,10^8);\n"
for d in ds:
    script += "E=ellinit([0,%d]); forprime(p=2,%d, print(%d,\" \",p,\" \",ellap(E,p)));\n" % (d, B, d)
    script += "print(\"COND \",%d,\" \",ellglobalred(E)[1]);\n" % d
lines = gp(script)
ap = {}; cond = {}
for l in lines:
    t = l.split()
    if t[0] == 'COND':
        cond[int(t[1])] = int(t[2])
    else:
        ap[(int(t[0]), int(t[1]))] = int(t[2])
splitB = [p for p in primary_primes(B) if p[1] != 0]
inertB = [(-p, 0) for p in range(5, B + 1) if is_prime_int(p) and p % 3 == 2]
err1 = 0; fal1 = 0; n1 = 0; err2 = 0; n2 = 0; fal2 = 0; err3 = 0
for d in ds:
    x = (4 * d, 0)
    for pi in splitB:
        p = norm(pi)
        if (6 * d) % p == 0:
            continue
        def psistar(q):
            c = chi_prime(q, x)
            return alpha(q).conjugate() * cx(c)
        def falsi(q):
            return alpha(q) * cx(chi_prime(q, x))
        s = math.sqrt(p) * (psistar(pi) + psistar(conj(pi)))
        err1 = max(err1, abs(s - ap[(d, p)]))
        fal1 = max(fal1, abs(math.sqrt(p) * (falsi(pi) + falsi(conj(pi))) - ap[(d, p)]))
        err3 = max(err3, abs(psistar(conj(pi)) - psistar(pi).conjugate()))
        n1 += 1
    for q in inertB:
        p = -q[0]
        if (6 * d) % p == 0:
            continue
        v = alpha(q).conjugate() * cx(chi_prime(q, x))
        err2 = max(err2, abs(v + 1) + abs(ap[(d, p)]))
        fal2 += abs(cx(chi_prime(q, x)) + 1) < 1e-9
        n2 += 1
rep('C1 split traces, %d pairs (d,pi), both primes above p' % n1, err1 < 1e-9 and fal1 > 1,
    'max err %.1e; falsifier alpha*chi(4d) max err %.1f' % (err1, fal1))
rep('C2 psi*(-p) = -1 and a_p = 0, %d inert pairs' % n2, err2 < 1e-12 and fal2 == 0,
    'max err %.1e; falsifier (no angular factor) equals -1 in %d cases' % (err2, fal2))
rep('C3 psi*(pibar) = conj psi*(pi)', err3 < 1e-12, 'max err %.1e' % err3)

# C4 periodicity of n -> chi_n(4d) modulo M = 144|d| (finite-order ray character), with falsifier M/4
def rand_primary_coprime(dd, lim):
    while True:
        k = random.randint(1, 3)
        pr = random.sample([p for p in P if (6 * dd) % norm(p) != 0 and norm(p) < 200 and
                            (p[1] != 0 or (6 * dd) % (-p[0]) != 0)], k)
        n = prod(pr)
        if norm(n) < lim:
            return pr, n
def factor_primary(n, plist):
    out = []; m = n
    for p in plist:
        while True:
            # test p | m: m * conj(p) divisible by N(p)
            t = mul(m, conj(p)); Np = norm(p)
            if t[0] % Np == 0 and t[1] % Np == 0 and norm(m) % Np == 0:
                m = (t[0] // Np, t[1] // Np); out.append(p)
            else:
                break
        if norm(m) == 1:
            break
    return out, m
plist = primary_primes(40000)
viol = 0; tested = 0; violf = 0
for d in [1, 2, 5, -3, 7]:
    M = 144 * abs(d)
    x = (4 * d, 0)
    for _ in range(60):
        pr, n = rand_primary_coprime(d, 5000)
        v = chi(pr, x)
        for _ in range(3):
            k = (random.randint(-3, 3), random.randint(-3, 3))
            n2 = (n[0] + M * k[0], n[1] + M * k[1])
            f2_, rest = factor_primary(n2, plist)
            if norm(rest) != 1:
                continue
            if any((6 * d) % norm(p) == 0 if p[1] != 0 else (6 * d) % (-p[0]) == 0 for p in f2_):
                continue
            tested += 1
            if abs(chi(f2_, x) - v) > 1e-9:
                viol += 1
            # falsifier: shift by 9 (a multiple of 3, so still primary, but not of the conductor)
            n3 = (n[0] + 9 * k[0], n[1] + 9 * k[1])
            if is_primary(n3) and n3 != n:
                f3_, rest3 = factor_primary(n3, plist)
                if norm(rest3) == 1 and not any((6 * d) % norm(p) == 0 if p[1] != 0 else (6 * d) % (-p[0]) == 0 for p in f3_):
                    violf += abs(chi(f3_, x) - v) > 1e-9
rep('C4 chi_n(4d) periodic mod 144|d| (%d tests)' % tested, viol == 0 and violf > 0,
    'violations %d; falsifier (shift by 9, not a multiple of the modulus) differs in %d cases' % (viol, violf))
print('   conductors N(E_d):', cond)
print('ALL PASS' if ok else 'SOME FAIL')

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
