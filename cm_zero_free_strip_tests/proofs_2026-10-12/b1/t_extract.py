"""Block b1 tests for the extraction lemmas (lem:b1-extract, lem:b1-uniform)."""
import math, random, bisect, os, sys
import numpy as np
import sympy as sp

# --- T3: #{m in Z[w]: N(m)<=x} <= 1+8x for all x (check at every attained norm value and just below 1)
R = 400
A, B = np.meshgrid(np.arange(-R-1, R+2), np.arange(-R-1, R+2))
Nv = (A*A - A*B + B*B).ravel()
Nv.sort()
X = 0.75*R*R   # disc fully inside the box (|A|,|B|<=R covers N<=0.75 R^2)
vals = np.unique(Nv[Nv <= X])
counts = np.searchsorted(Nv, vals, side='right')
ok = np.all(counts <= 1 + 8*vals)
fals = np.any(counts > 1 + 3*vals)
print('T3 lattice count <= 1+8x: ', ok, ' max ratio count/(1+8x)=%.4f' % np.max(counts/(1+8*vals)), ' falsifier 1+3x violated:', fals)

# --- T4: multiples count #{n primary, (n,6)=1, pi|n, N(n)<=X} <= 1+8X/N(pi)
def isprime(n):
    return sp.isprime(n)
def primary_primes_upto(Y):
    out = []
    for p in sp.primerange(5, Y+1):
        if p % 3 == 1:
            for a in range(-p, p+1):
                for b in range(-p, p+1, 3):
                    if a % 3 == 1 and a*a - a*b + b*b == p:
                        out.append(((a, b), p))
        elif p*p <= Y:
            out.append(((-p, 0), p*p))
    return out
def mul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)
Xmax = 3000
prim = [(a, b) for a in range(-80, 81) for b in range(-80, 81) if a % 3 == 1 and b % 3 == 0 and not (a % 2 == 0 and b % 2 == 0) and 0 < a*a-a*b+b*b <= Xmax]
# coprime to 3 automatic for primary; coprime to 2 = not both even
okm = True
for (pi, Npi) in primary_primes_upto(200):
    # n = pi*m with N(m) <= X/N(pi): count primary n with pi | n directly
    cnt = 0
    for n in prim:
        a, b = n; c, d = pi
        # pi | n iff n*conj(pi) has both coords divisible by N(pi)
        cpi = (c - d, -d)
        q = mul(n, cpi)
        if q[0] % Npi == 0 and q[1] % Npi == 0:
            cnt += 1
    if cnt > 1 + 8*Xmax/Npi: okm = False
print('T4 multiples count <= 1+8X/N(pi):', okm)

# --- T6: #{pi primary prime | d, N(pi) > Y/2} < 2 log|d| / log(Y/2) for Y>2
def primes_above_count(d, Y):
    c = 0
    for p, e in sp.factorint(abs(d)).items():
        if p in (2, 3): continue
        if p % 3 == 1:
            if p > Y/2: c += 2
        else:
            if p*p > Y/2: c += 1
    return c
random.seed(1)
ok6 = True; fal6 = False
for _ in range(20000):
    d = random.choice([1, -1]) * random.randint(1, 10**7)
    Y = random.uniform(2.5, 10**6)
    c = primes_above_count(d, Y)
    if c > 0 and not (c < 2*math.log(abs(d))/math.log(Y/2)): ok6 = False
for p in [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97]:
    Y = 2*p - 1; c = primes_above_count(p, Y)
    if c >= math.log(p)/math.log(Y/2): fal6 = True
print('T6 exclusion count < 2log|d|/log(Y/2):', ok6, ' falsifier (no factor 2) violated:', fal6)

# --- T7: prime supply P(Y) = #{pi primary prime, pi nmid 6, Y/2 < N(pi) <= Y} >= 2 c1 Y/log Y, sanity with c1 = 1/5
Ymax = 10**6
norms = []
for p in sp.primerange(5, Ymax+1):
    if p % 3 == 1: norms += [p, p]
    elif p*p <= Ymax: norms.append(p*p)
norms.sort()
def P(Y):
    return bisect.bisect_right(norms, Y) - bisect.bisect_right(norms, Y/2)
t7zero = all(P(Y) == 0 for Y in [1, 2, 3, 4, 5, 6, 6.99])
print('T7 P(Y)=0 for Y<7:', t7zero, ' P(7)=', P(7))
c1 = 0.2
bad = [Y for Y in range(7, Ymax+1) if P(Y) < 2*c1*Y/math.log(Y)]
print('T7 largest Y<=1e6 with P(Y) < 0.4 Y/log Y:', max(bad) if bad else None, ' min ratio P/(Y/logY) over [1e3,1e6]: %.3f' % min(P(Y)/(Y/math.log(Y)) for Y in range(1000, Ymax+1, 997)))

# --- T8: case-2 inequality: Y>=Y1>=7 and c1 Y/log Y < log H/(5 log(Y/2)) imply Y < 2 log H/(5 c1)
ok8 = True
for _ in range(200000):
    c1_ = random.uniform(0.01, 1); Y = random.uniform(7, 1e6); logH = random.uniform(0, 1e5)
    if c1_*Y/math.log(Y) < logH/(5*math.log(Y/2)):
        if not (Y < 2*logH/(5*c1_)): ok8 = False
okhalf = all(math.log(Y/2) >= 0.5*math.log(Y) for Y in np.linspace(4, 1e6, 100000))
print('T8 case-2 implication:', ok8, ' log(Y/2)>=log(Y)/2 for Y>=4:', okhalf)

# --- T9: exponent identities with Y^6 = D^{1+th}/N0
D, th, N0 = sp.symbols('D theta N0', positive=True)
Y = (D**(1+th)/N0)**sp.Rational(1, 6)
e1 = sp.simplify(sp.expand_power_base(D**(1+th/2)*Y**sp.Rational(-1, 2) / (N0**sp.Rational(1, 12)*D**(sp.Rational(11, 12)+5*th/12)), force=True))
e2 = sp.simplify(sp.expand_power_base(D/Y / (N0**sp.Rational(1, 6)*D**(sp.Rational(5, 6)-th/6)), force=True))
f1 = sp.simplify(sp.expand_power_base(D**(1+th/2)*Y**sp.Rational(-1, 2) / (N0**sp.Rational(1, 6)*D**(sp.Rational(11, 12)+5*th/12)), force=True))
print('T9 D^{1+th/2}Y^{-1/2} = N0^{1/12}D^{11/12+5th/12}:', e1 == 1, '; D/Y = N0^{1/6}D^{5/6-th/6}:', e2 == 1, '; falsifier N0^{1/6} in first:', f1 != 1)
# Y <= D for D>=1 when N0>=1, theta<=1/10: (1+th)/6 <= 1
print('T9b Y<=H^{1/6}=D^{(1+th)/6}<=D for th<=1/10:', (1+sp.Rational(1, 10))/6 <= 1)

# Exit status (added 2026-10-12): 1 unless every True/False printed above is True (T7's largest-Y line is a sanity
# number with no criterion) and P(7) = 2 (README row b1.20).  CM_TESTS_FORCE_FAIL=1 counts the first check as failed.
checks = [('T3 count <= 1+8x', bool(ok)), ('T3 falsifier 1+3x violated', bool(fals)), ('T4', okm),
          ('T6', ok6), ('T6 falsifier violated', fal6), ('T7 P(Y)=0 for Y<7', t7zero), ('T7 P(7)=2', P(7) == 2),
          ('T8 case-2 implication', ok8), ('T8 log(Y/2)>=log(Y)/2', okhalf), ('T9 first identity', e1 == 1),
          ('T9 second identity', e2 == 1), ('T9 falsifier N0^{1/6}', f1 != 1),
          ('T9b', bool((1+sp.Rational(1, 10))/6 <= 1))]
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    checks[0] = (checks[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, c in checks if not c]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(checks), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(checks)))
sys.exit(1 if bad else 0)
