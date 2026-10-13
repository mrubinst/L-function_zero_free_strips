"""Block b1 tests: every step of the point-count proof (prop:b1-points), exactly in Z[omega]."""
import sys, math
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'part_A_lines_0001_0549'))
from zw import mul, add, sub, neg, conj, N, divides, primary, UNITS, sym6_prime, prim_assoc

def isprime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True

def primary_primes_above(p):
    out = []
    for A in range(-p, p+1):
        for B in range(-p, p+1, 1):
            if A*A - A*B + B*B == p and A % 3 == 1 and B % 3 == 0:
                out.append((A, B))
    return out

Z = (0, 0); ONE = (1, 0)
def zpow(j):  # zeta^j as element of Z[omega], zeta = 1+omega
    return UNITS[j % 6]
def tr(x):  # x + conj x, an integer (as Z[omega] element (t,0))
    s = add(x, conj(x)); assert s[1] == 0; return s[0]

def chi_idx(t, pi):
    return sym6_prime((t, 0), pi)  # index or None

def char_val(t, pi, k):
    """chi_pi(t)^k as Z[omega] element, zero-extended."""
    j = chi_idx(t, pi)
    return Z if j is None else zpow(k*j)

fails = {'aff': 0, 'S': 0, 'Jratio': 0, 'J': 0, 'rk': 0, 'ap': 0, 'hasse': 0}
fals = {'J+pi': 0, 'rk_unconj': 0, 'ap_unconj': 0, 'ap_d': 0}
cases = 0
for p in [q for q in range(7, 400) if isprime(q) and q % 3 == 1]:
    pis = primary_primes_above(p)
    assert len(pis) == 2
    for d in [x for x in range(-30, 31) if x != 0 and x % p != 0]:
        # brute force count
        sq = {}
        for y in range(p): sq[y*y % p] = sq.get(y*y % p, 0) + 1
        aff = sum(sq.get((x**3 + d) % p, 0) for x in range(p))
        ap = p + 1 - (aff + 1)
        for pi in pis:
            cases += 1
            kap = {t: char_val(t, pi, 2) for t in range(p)}
            rho = {t: char_val(t, pi, 3) for t in range(p)}
            # affine count via characters: sum_V (1+rho(V+d))(1+kappa(V)+conj kappa(V))
            tot = Z
            for V in range(p):
                f1 = add(ONE, rho[(V + d) % p]); f2 = add(add(ONE, kap[V]), conj(kap[V]))
                tot = add(tot, mul(f1, f2))
            if tot != (aff, 0): fails['aff'] += 1
            S = Z
            for V in range(p): S = add(S, mul(rho[(V + d) % p], kap[V]))
            if add(add((p, 0), S), conj(S)) != (aff, 0): fails['aff'] += 1
            Jkr = Z; Jkk = Z
            for t in range(p):
                Jkr = add(Jkr, mul(kap[t], rho[(1 - t) % p]))
                Jkk = add(Jkk, mul(kap[t], kap[(1 - t) % p]))
            if S != mul(mul(rho[d % p], kap[d % p]), Jkr): fails['S'] += 1
            if Jkk != mul(conj(kap[4 % p]), Jkr): fails['Jratio'] += 1
            if Jkk != neg(pi): fails['J'] += 1
            if Jkk == pi: fals['J+pi'] += 1
            c4d = char_val(4*d, pi, 1)
            if mul(rho[d % p], kap[(4*d) % p]) != conj(c4d): fails['rk'] += 1
            if mul(rho[d % p], kap[(4*d) % p]) == c4d and c4d not in (ONE, (-1, 0)): fals['rk_unconj'] += 1
            if ap != tr(mul(pi, conj(c4d))): fails['ap'] += 1
            if ap == tr(mul(pi, c4d)) and ap != 0: fals['ap_unconj'] += 1
            cd = char_val(d, pi, 1)
            if ap == tr(mul(pi, conj(cd))): fals['ap_d'] += 1
            if ap*ap > 4*p: fails['hasse'] += 1
print('cases', cases)
print('failures (all must be 0):', fails)
print('falsifier hits (must be < cases; J+pi must be 0):', fals)

# Exit status (added 2026-10-12): 1 unless every failure counter is 0, the falsifier J = +pi never holds, and every
# falsifier misses in at least one case.  CM_TESTS_FORCE_FAIL=1 counts the first check as failed (self-test).
checks = [('failures[%s] == 0' % k, v == 0) for k, v in fails.items()]
checks += [('falsifier J+pi hits == 0', fals['J+pi'] == 0)]
checks += [('falsifier %s hits < cases' % k, v < cases) for k, v in fals.items()]
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    checks[0] = (checks[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, ok in checks if not ok]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(checks), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(checks)))
sys.exit(1 if bad else 0)
