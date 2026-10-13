"""s1 finite checks of the local data of the completed reflection (identical in the CM case, since no derivative
enters), at several good primes p (split and inert):

 (a) OAI78 eq:reflection-local-fourier:  C_{p,j}(h) = q^{-1/2} g^-_{p,j} chi_p(h)^{-j} (j!=0,h!=0), 0 (j!=0,h=0),
     -1/q (j=0,h!=0), 1-1/q (j=0,h=0);  |g^{pm}_{p,m}| = 1 for m != 0 mod 6.
 (b) the local sum (OAI78 L2224-2251), for every unit u = c/p mod p, sigma_p = lambda^2 u, eps_p = -(lambda^3 u sigma_p)^{-1}:
       sum_{h != 0} C_{p,j}(h) chi_p(sigma_p h)^{-2} e(eps_p h^{-1} x / p) = chi_p(sigma_p)^{-2} varpi_{p,j} B_p(x)
     for all x mod p and j = 0..5, with varpi, B_p as in OAI78 L1806-1840.
 (c) cross-prime phase OAI78 eq:reflection-cross-prime-phase:  chi_p(sigma_p)^{-2} varpi_{p,j} / chi_p(u)^{2j+2}
     is independent of the unit u (this is xi_{p,j}).
 (d) pair phase OAI78 eq:reflection-pair-phase:  chi_p(q)^{2j_p+2} chi_q(p)^{2j_q+2} = (q/p)_3^{j_p+j_q+2}.
 (e) marked combination (cor:marked-reflection):  1_{q|A} = 1 - chi_q(A)^0, inactive coefficient 1-(1-1/q) = 1/q,
     active coefficient: -(zero-mask transform) gives scalar -varpi_{q,0} = g^+_{q,2} chi_q(eps_q)^{-2}, local factor
     q^{-1/2} chi_q(x)^{-2}; checked by summing the active local transform of 1_{q|A} directly.
Falsifiers (must fail): (b) with B_p exponent -j instead of -j-2; (c) with exponent 2j instead of 2j+2;
 (d) with exponent j_p+j_q instead of j_p+j_q+2.
Usage: python3 t_s1_local_phases.py   (seconds)
"""
import os, sys, math, cmath, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zw import mul, add, sub, N, cx, divides, residues, sym6_prime, pw, reduce_mod, canon_mod
from fastgauss import primary_primes, cubic_prime
LAM = (1, 2)
def e_add(z): return cmath.exp(4j*math.pi*z.imag/math.sqrt(3))
PR = [p for p in primary_primes(200) if p != (-2, 0) and N(p) <= 61]   # good primes: split of norm <= 61, inert -5
print("primes:", PR)
worst = {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0, 'gauss': 0}
fals = {'b': 1e9, 'c': 1e9, 'd': 1e9}
for p in PR:
    q = N(p); RES = residues(p); P = cx(p)
    def chi(x, m):            # chi_p(x)^m with zero convention (also m = 0 mod 6)
        k = sym6_prime(x, p)
        return 0 if k is None else cmath.exp(1j*math.pi*k*(m % 6)/3)
    def inv(x):
        for y in RES:
            if divides(p, sub(mul(x, y), (1, 0))): return y
        raise ValueError
    def ep(y): return e_add(cx(y)/P)                   # e(y/p)
    def gpm(m, s): return sum(chi(y, m)*e_add(s*cx(y)/P) for y in RES)/math.sqrt(q)
    UNITS = [y for y in RES if not divides(p, y)]
    for m in range(1, 6):
        worst['gauss'] = max(worst['gauss'], abs(abs(gpm(m, 1))-1), abs(abs(gpm(m, -1))-1))
    for j in range(6):
        C = {h: sum(chi(y, j)*ep(mul((-1, 0), mul(h, y))) for y in RES)/q for h in RES}
        for h in RES:
            if j != 0:
                pred = 0 if divides(p, h) else q**-0.5*gpm(j, -1)*chi(h, -j)
            else:
                pred = (1-1/q) if divides(p, h) else -1/q
            worst['a'] = max(worst['a'], abs(C[h]-pred))
        ratios = []
        for u in UNITS[:12]:
            sig = mul(mul(LAM, LAM), u); eps = mul((-1, 0), inv(mul(mul(pw(LAM, 3), u), sig)))
            if j == 4: varpi = gpm(4, -1)
            elif j == 0: varpi = -gpm(2, 1)*chi(eps, -2)
            else: varpi = gpm(j, -1)*gpm(j+2, 1)*chi(eps, -j-2)
            for x in RES:
                lhs = sum(C[h]*chi(mul(sig, h), -2)*ep(mul(mul(eps, inv(h)), x)) for h in UNITS)
                if j == 4: Bx = q**-0.5*(-1+q*(1 if divides(p, x) else 0))
                elif j == 0: Bx = q**-0.5*chi(x, -2)
                else: Bx = chi(x, -j-2)
                rhs = chi(sig, -2)*varpi*Bx
                worst['b'] = max(worst['b'], abs(lhs-rhs))
                if j not in (0, 4):
                    fals['b'] = min(fals['b'], max(abs(lhs-chi(sig, -2)*varpi*chi(x, -j)), 0) if not divides(p, x) and abs(chi(x, 2)-1) > 1e-9 else 1e9)
            ratios.append(chi(sig, -2)*varpi/chi(u, 2*j+2))
            if j == 0:     # (e) active zero-mask transform of 1_{q|A} = 1 - chi^0: minus the active j=0 local sum
                for x in RES:
                    act = -sum(C[h]*chi(mul(sig, h), -2)*ep(mul(mul(eps, inv(h)), x)) for h in UNITS)
                    pred = chi(sig, -2)*(gpm(2, 1)*chi(eps, -2))*q**-0.5*chi(x, -2)
                    worst['e'] = max(worst['e'], abs(act-pred))
        worst['c'] = max(worst['c'], max(abs(r-ratios[0]) for r in ratios))
        wrong = [r*chi(u, 2) for r, u in zip(ratios, UNITS[:12])]          # ratio with exponent 2j instead of 2j+2
        fals['c'] = min(fals['c'], max(abs(r-wrong[0]) for r in wrong))
# (d) pair phase with cubic reciprocity, over pairs of distinct primes (including inert -5, -11, and split primes)
PP = [p for p in primary_primes(400) if p != (-2, 0)][:25]
for p, qq in itertools.permutations(PP, 2):
    if N(p) == N(qq) and p[1] == 0: continue
    for jp in range(6):
        for jq in range(6):
            def chi6(x, pr, m):
                k = sym6_prime(x, pr); return 0 if k is None else cmath.exp(1j*math.pi*k*(m % 6)/3)
            lhs = chi6(qq, p, 2*jp+2)*chi6(p, qq, 2*jq+2)
            rhs = cubic_prime(qq, p)**((jp+jq+2) % 3)
            worst['d'] = max(worst['d'], abs(lhs-rhs))
            if (jp+jq) % 3 != (jp+jq+2) % 3 and abs(cubic_prime(qq, p)-1) > 1e-9:
                fals['d'] = min(fals['d'], abs(lhs-cubic_prime(qq, p)**((jp+jq) % 3)))
print("max errors:", {k: f"{v:.1e}" for k, v in worst.items()})
print("falsifier minimal deviations (must be >> 0):", {k: f"{v:.2e}" for k, v in fals.items()})
ok = max(worst.values()) < 1e-9 and min(fals.values()) > 1e-3
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':   # self-test of the exit status (added 2026-10-12)
    ok = False; print("CM_TESTS_FORCE_FAIL=1: the check below is counted as failed")
print("PASS" if ok else "FAIL")

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
