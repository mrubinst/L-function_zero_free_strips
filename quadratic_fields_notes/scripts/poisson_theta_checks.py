"""Numerical checks for the sections on Poisson summation and on the cubic theta function
(quadratic-fields.tex, sec:poisson and sec:theta).

Usage: python3 poisson_theta_checks.py
Each check prints one line; the script stops with an AssertionError if a check fails.
Elements of Z[omega] are pairs (a, b) = a + b*omega.  Conventions as in the notes and the
papers: sqrt(-3) = 1 + 2 omega, e(z) = exp(4 pi i Im z / sqrt 3), e~(z) = exp(2 pi i (z + zbar)).
"""
import cmath, math, random, itertools
from math import gcd, sqrt, pi
import mpmath as mp

W = cmath.exp(2j*pi/3)
def mul(x, y): a, b = x; c, d = y; return (a*c - b*d, a*d + b*c - b*d)
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def neg(x): return (-x[0], -x[1])
def conj(x): a, b = x; return (a - b, -b)
def N(x): a, b = x; return a*a - a*b + b*b
def cx(x): return x[0] + x[1]*W
def divides(d, x):
    n = N(d); q = mul(x, conj(d)); return q[0] % n == 0 and q[1] % n == 0
def div(x, d):
    n = N(d); q = mul(x, conj(d)); assert q[0] % n == 0 and q[1] % n == 0; return (q[0]//n, q[1]//n)
def pw(x, e):
    r = (1, 0)
    for _ in range(e): r = mul(r, x)
    return r
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]          # zeta^j, zeta = 1 + omega = e^{i pi/3}
ZETA = [cmath.exp(1j*pi*j/3) for j in range(6)]
SQ3 = (1, 2)                                                        # sqrt(-3) = 1 + 2 omega
assert mul(SQ3, SQ3) == (-3, 0)
def primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
def prim_assoc(x):
    c = [mul(u, x) for u in UNITS if primary(mul(u, x))]; assert len(c) == 1; return c[0]
def isprime_int(n): return n > 1 and all(n % k for k in range(2, int(n**.5) + 1))
def e_add(z): return cmath.exp(4j*pi*z.imag/sqrt(3))
def e_check(z): return cmath.exp(2j*pi*2*z.real)

def hnf(m):
    """Z-basis (A, 0), (B, C) of the ideal mO, with A, C > 0 and 0 <= B < A."""
    v1, v2 = m, mul(m, (0, 1))
    # make the omega-coordinate of the second vector the gcd
    a1, b1 = v1; a2, b2 = v2
    while b1 != 0:
        q = b2 // b1
        a2, b2, a1, b1 = a1, b1, a2 - q*a1, b2 - q*b1
    C = abs(b2); B = a2 if b2 > 0 else -a2
    A = abs(a1)
    A = gcd(A, abs(N(m)) // C) if A else N(m)//C
    assert A*C == N(m), (m, A, C)
    return A, B % A, C
_hnf = {}
def canon(x, m):
    """Canonical representative of x modulo m."""
    if m not in _hnf: _hnf[m] = hnf(m)
    A, B, C = _hnf[m]
    k = x[1] // C
    x0 = x[0] - k*B; y0 = x[1] - k*C
    return (x0 % A, y0)
def residues(m):
    A, B, C = hnf(m)
    return [(x, y) for x in range(A) for y in range(C)]
for m in [(7, 0), (1, 3), (-5, 6), (2, 0), (1, -1), (4, 9)]:
    R = residues(m)
    assert len(R) == N(m) and len({canon(add(r, mul(m, (s, t))), m) for r in R for s in (-2, 3) for t in (1, -4)}) == N(m)
def powmod(x, e, m):
    r = (1, 0); b = canon(x, m)
    while e:
        if e & 1: r = canon(mul(r, b), m)
        b = canon(mul(b, b), m); e >>= 1
    return r

def primes_upto(B):
    out = []
    R = int(2*sqrt(B)) + 3
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            x = (a, b)
            if not primary(x) or N(x) > B: continue
            n = N(x)
            if (isprime_int(n) and n % 3 == 1) or (b == 0 and isprime_int(-a) and (-a) % 3 == 2):
                out.append(x)
    return sorted(out, key=lambda x: (N(x), x))
PR = primes_upto(20000)
PR_ODD = [p for p in PR if p != (-2, 0)]
def factor(n):
    """Primary n (coprime to 3) as a list of primary primes with multiplicity."""
    out = []; m = n; n0 = N(n)
    for p in PR:
        if N(m) == 1: break
        if N(p)*N(p) > n0 and N(m) > 1:
            break
        while divides(p, m):
            out.append(p); m = div(m, p)
    if N(m) > 1:
        out.append(prim_assoc(m))
    return out
def sym6(x, p):
    """(x/p)_6 as an index j (value zeta^j), or None if p | x; p a primary prime, p coprime to 6."""
    x = canon(x, p)
    if x == (0, 0): return None
    v = powmod(x, (N(p) - 1)//6, p)
    for j, u in enumerate(UNITS):
        if canon(u, p) == v: return j
    raise ValueError
def cubic(x, p):
    """(x/p)_3 as an exponent k (value omega^k), or None; p a primary prime coprime to 3 (p = -2 allowed)."""
    x = canon(x, p)
    if x == (0, 0): return None
    v = powmod(x, (N(p) - 1)//3, p)
    for k, u in enumerate([(1, 0), (0, 1), (-1, -1)]):
        if canon(u, p) == v: return k
    raise ValueError
_tab = {}
def chi(n, x, j=1):
    """chi_n(x)^j = (x/n)_6^j for primary n coprime to 6 (0 if not coprime)."""
    r = 1
    for p in factor(n):
        s = sym6(x, p)
        if s is None: return 0
        r *= ZETA[(j*s) % 6]
    return r
def cub(n, x):
    """(x/n)_3 for primary n coprime to 3, as a complex number."""
    r = 1
    for p in factor(n):
        s = cubic(x, p)
        if s is None: return 0
        r *= W**s
    return r
def mu(n):
    f = factor(n); return 0 if len(set(f)) < len(f) else (-1)**len(f)
def alpha(n): z = cx(n); return z/abs(z)

def gauss(psi, m):
    """N(m)^{-1/2} sum_{x mod m} psi(x) e(x/m)."""
    mm = cx(m)
    return sum(psi(x)*e_add(cx(x)/mm) for x in residues(m))/sqrt(N(m))
def gam(j, n): return gauss(lambda x: chi(n, x, j), n)
def ok(name, err, tol):
    assert err <= tol, (name, err)
    print(f"[PASS] {name}: max error {err:.2e}")

# ---------- 1. duality
B = lambda z, w: 2*(z*w).imag/sqrt(3)
gram = [[B(cx(u), cx(v)) for v in [(1, 0), (0, 1)]] for u in [(1, 0), (0, 1)]]
ok("Gram matrix of B on 1, omega is [[0,1],[1,-1]]", max(abs(gram[0][0]), abs(gram[0][1] - 1), abs(gram[1][0] - 1), abs(gram[1][1] + 1)), 1e-12)
err = 0
for z in [0.3 + 1.1j, -2.2 + 0.4j]:
    err = max(err, abs(e_add(z) - e_check(z/cx(SQ3))))
ok("e(z) = e~(z/sqrt(-3))", err, 1e-12)
# trace dual of O is (1/sqrt(-3)) O: Tr(x/sqrt(-3)) in Z for x in O
err = max(abs(2*(cx(x)/cx(SQ3)).real - round(2*(cx(x)/cx(SQ3)).real)) for x in [(1, 0), (0, 1), (3, -5)])
ok("Tr(x/sqrt(-3)) is an integer for x in O", err, 1e-12)
# orthogonality of additive characters
err = 0
for m in [(7, 0), (1, 3), (2, 0), (4, 9), (1, -1)]:
    for h in residues(m):
        s = sum(e_add(cx(mul(h, x))/cx(m)) for x in residues(m))
        err = max(err, abs(s - (N(m) if canon(h, m) == (0, 0) else 0)))
ok("orthogonality sum_{x mod m} e(hx/m) = N(m) 1_{m|h}", err, 1e-9)

# ---------- 2. Gauss sums
sqf = [n for n in [(a, b) for a in range(-40, 41) for b in range(-40, 41)]
       if primary(n) and 1 < N(n) <= 700 and N(n) % 2 and mu(n) != 0]
sqf.sort(key=lambda x: (N(x), x))
err_prim = err_abs = 0
for n in sqf[:60]:
    for j in (1, 2, 3, 4, 5):
        g = gam(j, n)
        err_abs = max(err_abs, abs(abs(g) - 1))
        for h in residues(n)[::7]:
            lhs = sum(chi(n, x, j)*e_add(cx(mul(h, x))/cx(n)) for x in residues(n))
            err_prim = max(err_prim, abs(lhs - chi(n, h, -j)*sqrt(N(n))*g))
ok("primitive Gauss sum: sum psi(x) e(hx/m) = conj psi(h) N(m)^{1/2} gamma(psi), psi = chi_n^j", err_prim, 1e-8)
ok("|gamma_j(n)| = 1 for squarefree n, j = 1..5", err_abs, 1e-10)
# identities of Lemma 4.1 of the paper
e1 = e2 = e3 = e4 = 0
for n in sqf:
    g1, gm1, g2, g3 = gam(1, n), gam(-1, n), gam(2, n), gam(3, n)
    G = chi(n, (4, 0), -1)*g3
    e1 = max(e1, abs(g1*gm1 - chi(n, (-1, 0))))
    e2 = max(e2, abs(g2**3 - mu(n)*alpha(n)))
    e3 = max(e3, abs(g1*g2 - mu(n)*alpha(n)*G))
    e4 = max(e4, abs(mu(n)*gm1 - chi(n, (-1, 0))*alpha(n).conjugate()*g2/G))
ok(f"gamma_1 gamma_-1 = chi_n(-1), {len(sqf)} squarefree n, N(n) <= 700", e1, 1e-9)
ok("gamma_2^3 = mu alpha", e2, 1e-9)
ok("gamma_1 gamma_2 = mu alpha G, G = conj chi_n(4) gamma_3", e3, 1e-9)
ok("Moebius-Gauss identity mu gamma_-1 = chi_n(-1) G^{-1} conj(alpha) gamma_2", e4, 1e-9)
# Hasse-Davenport duplication at primes: G(chi)G(chi^4) = conj chi(4) G(chi^2) G(chi^3)
err = 0
for p in [q for q in PR_ODD if N(q) <= 700]:
    gs = {j: gam(j, p)*sqrt(N(p)) for j in (1, 2, 3, 4)}
    err = max(err, abs(gs[1]*gs[4] - chi(p, (4, 0), -1)*gs[2]*gs[3])/N(p))
ok("Hasse-Davenport duplication at primes of norm <= 700", err, 1e-10)
# CRT factorization gamma_j(ab) = (chi_a(b) chi_b(a))^j gamma_j(a) gamma_j(b); paired characters
err = err2 = 0; pairs = 0
for a, b in itertools.combinations(sqf[:25], 2):
    if gcd(N(a), N(b)) != 1 or N(a)*N(b) > 2500: continue
    ab = mul(a, b); pairs += 1
    for j in (1, 2, 3, -1):
        f = (chi(a, b)*chi(b, a))
        err = max(err, abs(gam(j, ab) - f**j*gam(j, a)*gam(j, b)))
    psi = lambda x: chi(a, x, -1)*chi(b, x)
    R = chi(a, b, -1)*chi(b, a)
    err2 = max(err2, abs(gauss(psi, ab) - R*gam(-1, a)*gam(1, b)), abs(R*R - 1))
ok(f"CRT factorization of gamma_j over {pairs} coprime pairs", err, 1e-9)
ok("gamma(conj chi_z1 chi_z2) = R gamma_-1(z1) gamma_1(z2), R = +-1", err2, 1e-9)

# ---------- 3. Dunn-Radziwill cubic Gauss sums
def gDR(c, r=(1, 0)):
    cc = cx(c)
    return sum(cub(c, x)*e_check(cx(mul(r, x))/cc) for x in residues(c))
err = 0
for p in [q for q in PR_ODD if N(q) <= 400 and q[1] != 0]:
    P = N(p); a, b = p
    gnotes = sum(cub(p, (t, 0))*cmath.exp(2j*pi*t/P) for t in range(1, P))
    err = max(err, abs(gDR(p) - gnotes))
    for t in range(1, 6):
        err = max(err, abs(e_add(t/cx(p)) - cmath.exp(-2j*pi*b*t/P)))
ok("e(t/pi) = exp(-2 pi i b t/p), and DR's g(pi) = g(pi) of Section 3 for split pi", err, 1e-9)
cs = [c for c in [(a, b) for a in range(-25, 26) for b in range(-25, 26)] if primary(c) and 1 < N(c) <= 400]
err = err_r = err_g2 = 0
for c in cs:
    g = gDR(c)
    err = max(err, abs(g**3 - mu(c)*cx(c)**2*cx(conj(c))))
    if mu(c) != 0 and N(c) % 2:
        err_g2 = max(err_g2, abs(gam(2, c) - cub(c, SQ3)*g/abs(cx(c))))
    for r in [(2, 0), (0, 1), (-3, 0), (5, 1)]:
        if gcd(N(r), N(c)) == 1 and mu(c) != 0:
            err_r = max(err_r, abs(gDR(c, r) - cub(c, r).conjugate()*g))
ok(f"g(c)^3 = mu(c) c^2 cbar for {len(cs)} primary c of norm <= 400", err, 1e-6)
ok("g(r,c) = conj (r/c)_3 g(c) for (r,c) = 1", err_r, 1e-8)
ok("gamma_2(n) = (sqrt(-3)/n)_3 g(n)/|n|", err_g2, 1e-9)
ok("g(-2) = 2", abs(gDR((-2, 0)) - 2), 1e-12)
err = 0
sq = [c for c in cs if mu(c) != 0 and N(c) <= 120]
for a, b in itertools.combinations(sq, 2):
    if gcd(N(a), N(b)) != 1: continue
    ab = mul(a, b)
    lhs = gDR(ab)/abs(cx(ab)); rhs = cub(b, a).conjugate()*gDR(a)/abs(cx(a))*gDR(b)/abs(cx(b))
    err = max(err, abs(lhs - rhs))
ok("twisted multiplicativity g~(ab) = conj (a/b)_3 g~(a) g~(b)", err, 1e-9)

# ---------- 4. Poisson summation
def lattice_sum(f, R):
    s = 0
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            s += f((a, b))
    return s
err = 0
for a_ in (0.3, 1.0, 2.5):
    lhs = lattice_sum(lambda x: math.exp(-a_*N(x)), 40)
    rhs = 2*pi/(sqrt(3)*a_)*lattice_sum(lambda y: math.exp(-4*pi*pi*N(y)/(3*a_)), 40)
    err = max(err, abs(lhs - rhs)/lhs)
ok("Poisson for the Gaussian: sum e^{-aN(x)} = (2pi/(sqrt3 a)) sum e^{-4pi^2 N(y)/(3a)}", err, 1e-12)
# twisted Poisson lemma (paper's Lemma 4.2) with Phi(t) = exp(-t), Phihat(t) = (2pi/sqrt3) exp(-4 pi^2 t/3)
Phi = lambda t: math.exp(-t); Phihat = lambda t: 2*pi/sqrt(3)*math.exp(-4*pi*pi*t/3)
z1, z2, g = (1, 3), (4, 3), (-2, 3)                                  # norms 7, 13, 19
m = mul(z1, z2)
psi = lambda x: chi(z1, x, -1)*chi(z2, x)
gp = gauss(psi, m)
for H in (40.0, 300.0):
    R = int(3*sqrt(H)) + 6
    lhs = 0
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            k = (a, b)
            if divides(g, k): continue
            lhs += psi(k)*Phi(N(k)/H)
    rhs = 0
    for d, md in [((1, 0), 1), (g, -1)]:
        inner = 0
        Rh = 12
        for a in range(-Rh, Rh + 1):
            for b in range(-Rh, Rh + 1):
                h = (a, b)
                ph = psi(h).conjugate() if psi(h) != 0 else 0
                if ph: inner += ph*Phihat(H*N(h)/(N(d)*N(m)))
        rhs += md*psi(d)/N(d)*inner
    rhs *= H*gp/sqrt(N(m))
    ok(f"Poisson with a character and an excluded prime, H = {H:g}", abs(lhs - rhs)/max(1, abs(lhs)), 1e-9)

# ---------- 5. the Kubota character is a homomorphism on Gamma(3)
def mmul(A, Bm):
    (a, b), (c, d) = A; (e, f), (g_, h) = Bm
    return ((add(mul(a, e), mul(b, g_)), add(mul(a, f), mul(b, h))),
            (add(mul(c, e), mul(d, g_)), add(mul(c, f), mul(d, h))))
def kappa(Mx):
    (a, b), (c, d) = Mx
    assert primary(a) and primary(d) and c[0] % 3 == 0 and c[1] % 3 == 0
    return 1 if c == (0, 0) else cub(a, c)
def kappa_alt(Mx):
    (a, b), (c, d) = Mx
    return 1 if c == (0, 0) else cub(d, b)
random.seed(1)
# unipotent elements have kappa = 1, so products of them test nothing: build general elements of Gamma(3)
def rand_gamma3_general():
    while True:
        a = (3*random.randint(-4, 4) + 1, 3*random.randint(-4, 4))
        c = (3*random.randint(-4, 4), 3*random.randint(-4, 4))
        if c == (0, 0) or gcd(N(a), N(c)) != 1: continue
        m = mul((3, 0), c)
        inv = [x for x in residues(m) if canon(mul(a, x), m) == (1, 0)]
        if not inv: continue
        d = inv[0]
        bc = sub(mul(a, d), (1, 0))
        b = div(bc, c)
        assert b[0] % 3 == 0 and b[1] % 3 == 0 and primary(d)
        return ((a, b), (c, d))
cnt = nontriv = 0; err = err_alt = 0
for _ in range(400):
    A = rand_gamma3_general(); Bm = rand_gamma3_general()
    AB = mmul(A, Bm)
    if AB[1][0] == (0, 0): continue
    kA, kB, kAB = kappa(A), kappa(Bm), kappa(AB)
    err = max(err, abs(kAB - kA*kB))
    for Mx in (A, Bm, AB):
        if Mx[0][1] != (0, 0): err_alt = max(err_alt, abs(kappa(Mx) - kappa_alt(Mx)))
    cnt += 1; nontriv += abs(kAB - 1) > 1e-9
assert nontriv > cnt//3
ok(f"kappa(AB) = kappa(A) kappa(B) on {cnt} random pairs in Gamma(3) ({nontriv} with kappa(AB) != 1)", err, 1e-9)
maxnorm = 0
random.seed(1)
for _ in range(400):
    A = rand_gamma3_general(); Bm = rand_gamma3_general(); AB = mmul(A, Bm)
    maxnorm = max(maxnorm, max(N(x) for Mx in (A, Bm, AB) for row in Mx for x in row))
print(f"       (largest entry norm in these matrices: {maxnorm})")
# the action on H^3 (quaternion formula) and the cusp coordinates of g^{-1}(a/c + z, v)
def act(Mx, w):
    (a, b), (c, d) = [[cx(x) for x in row] for row in Mx]
    z, v = w
    den = abs(c*z + d)**2 + abs(c)**2*v*v
    return (((a*z + b)*(c*z + d).conjugate() + a*c.conjugate()*v*v)/den, v/den)
def minv(Mx):
    (a, b), (c, d) = Mx
    return ((d, neg(b)), (neg(c), a))
err = err_g = 0
for _ in range(50):
    A = rand_gamma3_general(); Bm = rand_gamma3_general()
    w = (complex(random.uniform(-1, 1), random.uniform(-1, 1)), random.uniform(0.2, 2))
    w1 = act(mmul(A, Bm), w); w2 = act(A, act(Bm, w))
    err = max(err, abs(w1[0] - w2[0]) + abs(w1[1] - w2[1]))
    (a, bg), (c, dl) = A
    z = complex(random.uniform(-.5, .5), random.uniform(-.5, .5)); v = random.uniform(0.1, 1)
    lhs = act(minv(A), (cx(a)/cx(c) + z, v))
    cc = cx(c)
    rhs = (-cx(dl)/cc - z.conjugate()/(cc*cc*(v*v + abs(z)**2)), v/(abs(cc)**2*(v*v + abs(z)**2)))
    err_g = max(err_g, abs(lhs[0] - rhs[0]) + abs(lhs[1] - rhs[1]))
ok("the action is an action: (gh)w = g(hw)", err, 1e-9)
ok("g^{-1}(a/c + z, v) = (-delta/c - zbar/(c^2(v^2+|z|^2)), v/(N(c)(v^2+|z|^2)))", err_g, 1e-9)
ok("kappa = (c/a)_3 = (b/d)_3", err_alt, 1e-9)

# ---------- 6. hyperbolic Laplacian, Bessel integral
mp.mp.dps = 30
def lap(f, x, y, v, h=mp.mpf('1e-6')):
    fxx = (f(x + h, y, v) - 2*f(x, y, v) + f(x - h, y, v))/h**2
    fyy = (f(x, y + h, v) - 2*f(x, y, v) + f(x, y - h, v))/h**2
    fvv = (f(x, y, v + h) - 2*f(x, y, v) + f(x, y, v - h))/h**2
    fv = (f(x, y, v + h) - f(x, y, v - h))/(2*h)
    return v**2*(fxx + fyy + fvv) - v*fv
muv = mp.mpc(0.4, -0.7)
F1 = lambda x, y, v: v*mp.besselk(mp.mpf(1)/3, 4*mp.pi*abs(muv)*v)*mp.exp(2j*mp.pi*2*mp.re(muv*mp.mpc(x, y)))
F0 = lambda x, y, v: v**(mp.mpf(2)/3)
err = 0
for (x, y, v) in [(0.1, 0.2, 0.3), (-0.4, 0.05, 0.8)]:
    x, y, v = mp.mpf(x), mp.mpf(y), mp.mpf(v)
    err = max(err, abs(lap(F1, x, y, v) + mp.mpf(8)/9*F1(x, y, v))/abs(F1(x, y, v)),
              abs(lap(F0, x, y, v) + mp.mpf(8)/9*F0(x, y, v))/abs(F0(x, y, v)))
ok("v^{2/3} and v K_{1/3}(4pi|mu|v) e~(mu z) have Laplace eigenvalue -8/9", float(err), 1e-8)
err = 0
for s in [mp.mpf('1.3'), mp.mpc('0.7', '2.0')]:
    L = mp.mpf('0.37')
    lhs = mp.quad(lambda v: v**(2*s)*mp.besselk(mp.mpf(1)/3, 4*mp.pi*L*v), [0, 1, mp.inf])
    rhs = 2**(2*s - 1)*mp.gamma(s + mp.mpf(1)/3)*mp.gamma(s + mp.mpf(2)/3)*(4*mp.pi*L)**(-2*s - 1)
    err = max(err, abs(lhs - rhs)/abs(rhs))
ok("int_0^oo v^{2s} K_{1/3}(4pi|l|v) dv = 2^{2s-1} Gamma(s+1/3) Gamma(s+2/3) (4pi|l|)^{-2s-1}", float(err), 1e-15)
z = mp.mpc(0, 1)*mp.sqrt(3)
ok("alpha((sqrt -3)^{-3}) = i", float(abs(z**-3/abs(z**-3) - 1j)), 1e-20)

# ---------- 7. the local transform: a sextic twist becomes chi^{-j-2}
err = 0; nchk = 0
for p in [q for q in PR_ODD if N(q) <= 200]:
    res = [x for x in residues(p) if x != (0, 0)]
    inv = {}
    for x in res:
        for y in res:
            if canon(mul(x, y), p) == (1, 0): inv[x] = y; break
    for j in range(6):
        C = {}
        for h in residues(p):
            C[h] = sum((chi(p, y, j) if j else (0 if canon(y, p) == (0, 0) else 1))*e_add(-cx(mul(h, y))/cx(p))
                       for y in residues(p))/N(p)
        for sig, eps in [((2, 1), (1, 0)), ((1, 3), (5, 2))]:
            if divides(p, sig) or divides(p, eps): continue
            ratios = []
            for x in residues(p):
                lhs = sum(C[h]*chi(p, mul(sig, h), -2)*e_add(cx(mul(eps, mul(inv[h], x)))/cx(p)) for h in res)
                if j == 4:
                    Bx = (-1 + (N(p) if canon(x, p) == (0, 0) else 0))/sqrt(N(p))
                elif j == 0:
                    Bx = chi(p, x, -2)/sqrt(N(p))
                else:
                    Bx = chi(p, x, -j - 2)
                if abs(Bx) > 1e-12:
                    ratios.append(lhs/(chi(p, sig, -2)*Bx))
                else:
                    err = max(err, abs(lhs))
            w0 = ratios[0]
            err = max(err, abs(abs(w0) - 1), max(abs(r - w0) for r in ratios))
            nchk += 1
ok(f"local transform sum_h C_pj(h) chi(sigma h)^-2 e(eps h^-1 x/p) = unit * B_pj(x) ({nchk} cases)", err, 1e-9)
print("all checks passed")
