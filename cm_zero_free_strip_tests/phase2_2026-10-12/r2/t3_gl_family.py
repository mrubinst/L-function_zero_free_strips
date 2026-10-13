"""r2 test T3: the Goldmakher-Louvel Hecke-family axioms for the quadratic family of Section 9.1 (edit E6 of r2).

For squarefree primary k coprime to 6 put chi_(k)(x) = (x/k)_2 kappa_lambda(x)^{e_k}, where kappa_lambda is the
nontrivial character modulo lambda (x -> x mod lambda in {+1,-1}) and N(k) = (-1)^{e_k} mod 4.  We check:
  (0) chi_(k) is trivial on the six units (so it is a Hecke character of trivial infinite type), and
      chi_(k)(n) = (n/k)_2 for primary n;
  (1) values are +-1 or 0;
  (2) reciprocity: for coprime k, k', chi_(k)(k') chi_(k')(k) = (k'/k)_2 (k/k')_2 depends only on (k mod 4, k' mod 4)
      [GL, Definition 1 (2) with G = (O/4O)^x and [(k)] = k mod 4];
  (3) if k1 = k2 mod 4 are coprime, then e_{k1} = e_{k2} and chi_(k1) conj(chi_(k2)) = (./k1k2)_2, which is trivial on
      the units (a Hecke character modulo k1k2) and nontrivial at each prime of k1k2 (primitive);
  (4) even columns: (-2n'/k)_2 = (-2/k)_2 (n'/k)_2.
Falsifiers (must fail): reciprocity factor as a function of the classes modulo 2 (instead of 4); chi_(k) without the
factor kappa_lambda^{e_k} is not trivial on the units.
Usage: python3 t3_gl_family.py     (a few seconds)
"""
import random, math, itertools, os, sys
random.seed(7)
def mul(x, y): a, b = x; c, d = y; return (a*c - b*d, a*d + b*c - b*d)
def sub(x, y): return (x[0]-y[0], x[1]-y[1])
def conj(x): return (x[0]-x[1], -x[1])
def N(x): a, b = x; return a*a - a*b + b*b
def divides(d, x):
    n = N(d); q = mul(x, conj(d)); return q[0] % n == 0 and q[1] % n == 0
def div(x, d):
    n = N(d); q = mul(x, conj(d)); return (q[0]//n, q[1]//n)
def red(y, m):
    n = N(m); q = mul(y, conj(m)); return sub(y, mul(((2*q[0]+n)//(2*n), (2*q[1]+n)//(2*n)), m))
def powmod(x, e, m):
    r = (1, 0); b = red(x, m)
    while e:
        if e & 1: r = red(mul(r, b), m)
        b = red(mul(b, b), m); e >>= 1
    return r
def primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
def prim(x):
    for u in UNITS:
        if primary(mul(u, x)): return mul(u, x)
# primary primes coprime to 6 of norm <= B
B = 3000
P = []
for a in range(-60, 61):
    for b in range(-60, 61):
        n = N((a, b))
        if 1 < n <= B and primary((a, b)) and n % 2 and n % 3:
            p = (a, b)
            # prime iff norm is a rational prime, or norm = q^2 with q = 2 mod 3 prime and p = -q
            def isp(m): return m > 1 and all(m % d for d in range(2, int(m**.5)+1))
            if isp(n) or (b == 0 and isp(-a) and (-a) % 3 == 2):
                P.append(p)
P = sorted(set(P), key=N)
def quad_prime(x, pi):
    if divides(pi, x): return 0
    v = powmod(x, (N(pi)-1)//2, pi)
    if divides(pi, sub(v, (1, 0))): return 1
    if divides(pi, sub(v, (-1, 0))): return -1
    raise ValueError
def factor(k):
    out = []; m = k
    for p in P:
        if N(m) == 1: break
        while divides(p, m): out.append(p); m = div(m, p)
    assert N(m) == 1
    return out
def quad(x, k):           # (x/k)_2, Jacobi-type
    r = 1
    for p in factor(k): r *= quad_prime(x, p)
    return r
def kappa_lam(x):         # x mod lambda in {+1,-1}; lambda | (x - (x0 + x1)) since omega = 1 mod lambda
    s = (x[0] + x[1]) % 3
    return {1: 1, 2: -1, 0: 0}[s]
def e_k(k): return 0 if N(k) % 4 == 1 else 1
def chi(k, x): return quad(x, k)*kappa_lam(x)**e_k(k)
# squarefree primary k coprime to 6, N(k) <= 1500
K = []
for a in range(-45, 46):
    for b in range(-45, 46):
        k = (a, b)
        if primary(k) and 1 < N(k) <= 1500 and N(k) % 2 and N(k) % 3:
            f = factor(k)
            if len(set(f)) == len(f): K.append(k)
print("squarefree primary k coprime to 6:", len(K))
ok = True
# (0), (1)
bad0 = 0; bad0f = 0
for k in K:
    for u in UNITS:
        if chi(k, u) != 1: bad0 += 1
        if quad(u, k) != 1: bad0f += 1
    for _ in range(5):
        n = random.choice(K)
        if chi(k, n) != quad(n, k): bad0 += 1
print("(0) chi_(k) trivial on units and = (n/k)_2 on primary n: violations", bad0, "| falsifier without kappa^e: units with (u/k)_2 != 1:", bad0f)
ok &= bad0 == 0 and bad0f > 0
def mod4(x): return (x[0] % 4, x[1] % 4)
def mod2(x): return (x[0] % 2, x[1] % 2)
tab4 = {}; tab2 = {}; viol4 = 0; viol2 = 0; npairs = 0
for _ in range(6000):
    k, kp = random.sample(K, 2)
    if set(factor(k)) & set(factor(kp)): continue
    v = chi(k, kp)*chi(kp, k); npairs += 1
    key4 = (mod4(k), mod4(kp)); key2 = (mod2(k), mod2(kp))
    if tab4.setdefault(key4, v) != v: viol4 += 1
    if tab2.setdefault(key2, v) != v: viol2 += 1
print(f"(2) reciprocity factor as a function of classes mod 4: {npairs} coprime pairs, {len(tab4)} class pairs, violations {viol4}"
      f" | falsifier mod 2: violations {viol2}")
ok &= viol4 == 0 and viol2 > 0
viol3 = 0; n3 = 0
for _ in range(4000):
    k1, k2 = random.sample(K, 2)
    if mod4(k1) != mod4(k2) or set(factor(k1)) & set(factor(k2)): continue
    n3 += 1
    if e_k(k1) != e_k(k2): viol3 += 1
    k12 = mul(k1, k2)
    for u in UNITS:
        if quad(u, k1)*quad(u, k2) != 1: viol3 += 1
    for _ in range(4):        # compare on x prime to lambda: the product is induced from (./k1k2)_2 (GL's primitive convention)
        x = (random.randint(-500, 500), random.randint(-500, 500))
        if kappa_lam(x) == 0: continue
        if chi(k1, x)*chi(k2, x) != quad(x, k12): viol3 += 1
    for p in set(factor(k12)):     # nontrivial at each prime: some x with (x/p)_2 = -1
        if not any(quad_prime((a, b), p) == -1 for a in range(0, 6) for b in range(0, 6)): viol3 += 1
print(f"(3) {n3} coprime pairs k1 = k2 mod 4: e equal, product = (./k1k2)_2, trivial on units, primitive: violations {viol3}")
ok &= viol3 == 0 and n3 > 20
viol4e = 0
for _ in range(500):
    k = random.choice(K); n1 = random.choice(K)
    if quad(mul((-2, 0), n1), k) != quad((-2, 0), k)*quad(n1, k): viol4e += 1
print("(4) (-2n'/k)_2 = (-2/k)_2 (n'/k)_2: violations", viol4e)
ok &= viol4e == 0
# Exit status (added 2026-10-12): 1 unless ok.  CM_TESTS_FORCE_FAIL=1 forces FAIL (self-test of the wiring).
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    ok = False; print("CM_TESTS_FORCE_FAIL=1: the checks are counted as failed")
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
