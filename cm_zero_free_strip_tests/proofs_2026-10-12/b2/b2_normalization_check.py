import sys, math, cmath
sys.path.insert(0, sys.argv[1])
from cmlib import *
# gamma_2 with modulus generator m (ideal (n), n primary): N^{-1/2} sum_{x mod n} chi_n(x)^2 e(x/m)
def g2_gen(n, m):
    mm = cx(m); s = 0
    for x in residues(n):
        c = chi(n, x)
        if c != 0: s += c**2*cmath.exp(4j*math.pi*(cx(x)/mm).imag/math.sqrt(3))
    return s/math.sqrt(N(n))
mx = 0; my = 0
for p in [q for q in PR_ODD if N(q) < 200][:12]:
    a = g2_gen(p, p); b = g2_gen(p, (-p[0], -p[1]))
    mx = max(mx, abs(a - b))
    my = max(my, abs(b**3 - (-alpha_c(p))), abs(a**3 + alpha_c(p)))
print("gamma_2(-pi) = gamma_2(pi):", mx, " gamma_2(pi)^3 = -alpha(pi) = alpha(-pi):", my)
# gamma(psi) with modulus generator zeta*q equals psi(zeta) gamma(psi), psi = conj(chi_z1) chi_z2, q = z1 z2
import random
random.seed(3)
sq = [n for n in primaries(200, sqfree=True) if n != (1, 0)]
mz = 0; cnt = 0
while cnt < 20:
    z1, z2 = random.choice(sq), random.choice(sq)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 2500: continue
    cnt += 1; q = mul(z1, z2)
    def gam_gen(m):
        mm = cx(m); s = 0
        for x in residues(q):
            v = chi(z1, x).conjugate()*chi(z2, x)
            if v != 0: s += v*cmath.exp(4j*math.pi*(cx(x)/mm).imag/math.sqrt(3))
        return s/math.sqrt(N(q))
    g0 = gam_gen(q)
    for u in UNITS:
        psi_u = chi(z1, u).conjugate()*chi(z2, u)
        mz = max(mz, abs(gam_gen(mul(u, q)) - psi_u*g0))
print("gamma(psi) with generator zeta q = psi(zeta) gamma(psi), 20 pairs x 6 units:", mz)

# Exit status (added 2026-10-12): 1 unless every identity error above is <= 1e-10 (no falsifier here).
# CM_TESTS_FORCE_FAIL=1 counts the first check as failed (self-test of the wiring).
checks = [("gamma_2(-pi) = gamma_2(pi) err <= 1e-10", mx <= 1e-10), ("gamma_2(pi)^3 = -alpha(pi) err <= 1e-10", my <= 1e-10),
          ("generator change of gamma(psi) err <= 1e-10", mz <= 1e-10)]
import os
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    checks[0] = (checks[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, c in checks if not c]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(checks), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(checks)))
sys.exit(1 if bad else 0)
