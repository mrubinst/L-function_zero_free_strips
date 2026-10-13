"""Helpers for the r2 tests: exact Z[omega] arithmetic with big integers, residue symbols, and the cusp data of
[OAI1112, App. A.2] with L = 72 (copied verbatim from t2_grouping_cusp_data.py)."""
import cmath, math, random, subprocess, sys, time
random.seed(20261012)
W = cmath.exp(2j*math.pi/3)

# ----------------------------------------------------------------- exact Z[omega] arithmetic, x = a + b*omega
def mul(x, y): a, b = x; c, d = y; return (a*c - b*d, a*d + b*c - b*d)
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def sub(x, y): return (x[0]-y[0], x[1]-y[1])
def neg(x): return (-x[0], -x[1])
def conj(x): return (x[0]-x[1], -x[1])
def N(x): a, b = x; return a*a - a*b + b*b
def rdiv(a, n):            # nearest integer to a/n, exact (n > 0)
    return (2*a + n)//(2*n)
def divides(d, x):
    n = N(d); q = mul(x, conj(d)); return q[0] % n == 0 and q[1] % n == 0
def div(x, d):
    n = N(d); q = mul(x, conj(d)); assert q[0] % n == 0 and q[1] % n == 0, (x, d); return (q[0]//n, q[1]//n)
def red(y, m):
    n = N(m); q = mul(y, conj(m)); return sub(y, mul((rdiv(q[0], n), rdiv(q[1], n)), m))
def canon(y, m):           # canonical representative mod m (for comparing residues)
    n = N(m); q = mul(y, conj(m)); return sub(y, mul((q[0]//n, q[1]//n), m))
def gcd(x, y):
    while y != (0, 0):
        x, y = y, red(x, y)
    return x
def egcd(x, y):            # returns (g, u, v) with u x + v y = g
    u0, v0, u1, v1 = (1, 0), (0, 0), (0, 0), (1, 0)
    while y != (0, 0):
        n = N(y); q = mul(x, conj(y)); qq = (rdiv(q[0], n), rdiv(q[1], n))
        x, y = y, sub(x, mul(qq, y))
        u0, u1 = u1, sub(u0, mul(qq, u1)); v0, v1 = v1, sub(v0, mul(qq, v1))
    return x, u0, v0
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
def inv_mod(x, m):
    g, u, v = egcd(x, m)
    for e in UNITS:
        if g == e:
            return red(mul(u, conj_unit(e)), m)
    raise ValueError(("not invertible", x, m))
def conj_unit(e):          # inverse of a unit
    for f in UNITS:
        if mul(e, f) == (1, 0): return f
def pw(x, k):
    r = (1, 0)
    for _ in range(k): r = mul(r, x)
    return r
def crt(pairs):            # pairs (residue, modulus) with pairwise coprime moduli -> (x, M)
    x, Mo = (0, 0), (1, 0)
    for r_, m_ in pairs:
        # x' = x + Mo * t with x + Mo t = r_ mod m_
        t = mul(sub(r_, x), inv_mod(Mo, m_))
        x = add(x, mul(Mo, red(t, m_))); Mo = mul(Mo, m_); x = red(x, Mo)
    return x, Mo
def primary(x): return x[0] % 3 == 1 and x[1] % 3 == 0
def mod3eq(x, y): d = sub(x, y); return d[0] % 3 == 0 and d[1] % 3 == 0
LAM = (1, 2); TWO = (2, 0)
def val(q, x):
    k = 0
    while x != (0, 0) and divides(q, x): x = div(x, q); k += 1
    return k

# ----------------------------------------------------------------- residue symbols
def cube_root_index(v, pi):          # v = 1, omega, omega^2 mod pi -> 0,1,2
    for k, u in enumerate([(1, 0), (0, 1), (-1, -1)]):
        if divides(pi, sub(v, u)): return k
    raise ValueError
def powmod(x, e, m):
    r = (1, 0); b = red(x, m)
    while e:
        if e & 1: r = red(mul(r, b), m)
        b = red(mul(b, b), m); e >>= 1
    return r
def cubic_sym_prime(x, pi):          # (x/pi)_3 for a prime pi not above 3; 0 if pi | x
    if divides(pi, x): return 0
    return W**cube_root_index(powmod(x, (N(pi)-1)//3, pi), pi)
def sext_sym_prime(x, pi):           # (x/pi)_6
    if divides(pi, x): return 0
    v = powmod(x, (N(pi)-1)//6, pi)
    for k, u in enumerate(UNITS):    # UNITS[k] = zeta^k, zeta = 1 + omega = e^{i pi/3}
        if divides(pi, sub(v, u)): return cmath.exp(1j*math.pi*k/3)
    raise ValueError
def gp_factor(n):
    out = subprocess.run(['gp', '-q', '-f'], input=f'print(factor({n}))\n', capture_output=True, text=True).stdout.strip()
    out = out.replace('Mat(', '').replace(')', '')
    if out.startswith('['): out = out[1:-1]
    res = []
    for row in out.split(';'):
        p, e = row.split(',')
        res.append((int(p), int(e)))
    return res
def sqrt_mod_m3(p):                  # a root w of w^2 + w + 1 = 0 mod p (p = 1 mod 3)
    for g in range(2, p):
        w = pow(g, (p-1)//3, p)
        if w != 1: return w
def prime_factors(d):                # prime elements of O dividing d (d coprime to 3), with multiplicity
    out = []; m = d
    for p, e in gp_factor(N(d)):
        assert p != 3
        if p == 2 or p % 3 == 2:
            q = (p, 0); k = 0
            while divides(q, m): m = div(m, q); k += 1
            out.append((q, k))
        else:
            w = sqrt_mod_m3(p)
            for root in (w, (-1-w) % p):
                pi = gcd((p, 0), (-root, 1))
                k = 0
                while divides(pi, m): m = div(m, pi); k += 1
                if k: out.append((pi, k))
    assert N(m) == 1, (d, m)
    return out
def cubic_sym(x, d):                 # (x/d)_3, d coprime to 3 (Jacobi-type product over prime ideals)
    if N(d) == 1: return 1
    r = 1
    for pi, k in prime_factors(d):
        r *= cubic_sym_prime(x, pi)**k
    return r
# supplements for primary a (used for the (A.10) formulas)
def lam3(a): assert primary(a); return W**((-(a[1]//3)) % 3)          # (lambda/a)_3
def om3(a): return W**(((N(a)-1)//3) % 3)                            # (omega/a)_3
def two3(a): return cubic_sym_prime(a, TWO)                          # (a/2)_3 = (2/a)_3 for primary odd a
def unit3(u, a):                                                    # (u/a)_3 for a unit u, a primary
    k = UNITS.index(u)   # u = zeta^k; zeta = -omega^2 so (zeta/a)_3 = (omega/a)_3^2
    return om3(a)**(2*k)
def sym_c0_over_a(c0, a):            # (c0/a)_3 for c0 = unit * 2^i * lambda^j, a primary odd
    i = val(TWO, c0); j = val(LAM, c0); u = div(c0, mul(pw(TWO, i), pw(LAM, j)))
    assert u in UNITS
    return unit3(u, a)*two3(a)**i*lam3(a)**j

# ----------------------------------------------------------------- additive characters
def e_frac(num, den):                # e(num/den), num in O, den a positive rational integer: exp(2 pi i b/den)
    return cmath.exp(2j*math.pi*((num[1] % den)/den))
def e_K(num, den_el):                # e(num/den_el) for num, den_el in O
    n = N(den_el); return e_frac(mul(num, conj(den_el)), n)
def inO(num, den_el): return divides(den_el, num)

# ----------------------------------------------------------------- fixed data
L = (72, 0)
M = mul(pw(LAM, 12), pw(L, 4)); M2 = mul(M, M)
assert M == (729*72**4, 0)
vM = {TWO: val(TWO, M), LAM: val(LAM, M)}
assert vM == {TWO: 12, LAM: 28}, vM
PI7 = [p for p in [(3, 1), (-2, -3), (1, 3), (-1, -3), (3, 2)] if primary(p) and N(p) == 7][0]
PI13 = [(a, b) for a in range(-6, 7) for b in range(-6, 7) if N((a, b)) == 13 and primary((a, b))][0]
PI_INERT = (-5, 0)
def nonzero_res(pi):
    n = N(pi); seen = {}
    for a in range(n):
        for b in range(n):
            x = canon((a, b), pi)
            if x != (0, 0) and x not in seen: seen[x] = 1
    return list(seen)

def gauss(pi, j):                    # gamma_j(pi) = N^{-1/2} sum_{x != 0 mod pi} chi_pi(x)^j e(x/pi)
    s = 0
    for x in nonzero_res(pi):
        s += sext_sym_prime(x, pi)**j*e_K(x, pi)
    return s/math.sqrt(N(pi))
def Ccoef(pi, j, h):                 # C_{pi,j}(h) = N^{-1} sum_y chi_pi(y)^j e(-h y/pi)  (chi^0 = 1_{pi nmid y})
    s = 0
    for y in nonzero_res(pi):
        s += sext_sym_prime(y, pi)**j * e_K(neg(mul(h, y)), pi)
    return s/N(pi)
def Bloc(pi, j, x):                  # (8.5)
    n = N(pi); ch = sext_sym_prime(x, pi)
    if j not in (0, 4): return ch**(-j-2) if ch != 0 else 0
    if j == 0: return n**-0.5*(ch**-2 if ch != 0 else 0)
    return n**-0.5*(-1 + (n if divides(pi, x) else 0))

# ----------------------------------------------------------------- one translate
def normalize(a, c):
    if divides(LAM, c):
        for u in UNITS:
            if mod3eq(mul(u, a), (1, 0)): return mul(u, a), mul(u, c)
    for u in UNITS:
        if mod3eq(mul(u, c), (1, 0)): return mul(u, a), mul(u, c)
    raise ValueError
def translate_data(h0, A, hs):
    """A: list of primes; hs: lifted h_pi (= 0 mod M^2).  Returns the cusp data of z = lambda^2 (h0/72 + sum h/pi)."""
    r = (1, 0)
    for p in A: r = mul(r, p)
    den = mul(L, r)
    num = mul(h0, r)
    for p, h in zip(A, hs): num = add(num, mul(L, mul(h, div(r, p))))
    num = mul(mul(LAM, LAM), num)
    if num == (0, 0): raise ValueError
    g = gcd(num, den); a, c = div(num, g), div(den, g)
    a, c = normalize(a, c)
    assert divides(r, c); c0 = div(c, r)
    # c0 | 8 lambda^2 up to a unit
    assert divides(c0, mul((8, 0), mul(LAM, LAM))), c0
    # delta' (OAI App. A.2)
    pairs = []
    for q in (TWO, LAM):
        vq0 = val(q, c0)
        if vq0 > 0:
            mod = pw(q, vM[q] + vq0); pairs.append((inv_mod(a, mod), mod))
        else:
            mod = pw(q, vM[q]); pairs.append(((0, 0), mod))
    pairs.append((inv_mod(a, r), r))
    dp, _ = crt(pairs)
    t = sub(mul(a, dp), (1, 0)); assert divides(c, t); bg = div(t, c)
    vl = val(LAM, c)
    if vl >= 2:
        H = ((1, 0), (0, 0), (0, 0), (1, 0)); sigma = '0'; u0 = None
    elif vl == 1:
        u0 = LAM if mod3eq(c, LAM) else neg(LAM); assert mod3eq(c, u0)
        H = ((1, 0), (0, 0), u0, (1, 0)); sigma = '-' if u0 == LAM else '+'
    else:
        u0 = (a[0] % 3, a[1] % 3)
        H = (u0, (-1, 0), (1, 0), (0, 0)); sigma = {0: '0', 1: '-', 2: '+'}[u0[1]]
    # g1 = g H^{-1};  H^{-1} = (d -b; -c a)
    A1, B1, C1, D1 = H; Hi = (D1, neg(B1), neg(C1), A1)
    ga, gb, gc, gd = a, bg, c, dp
    g1 = (add(mul(ga, Hi[0]), mul(gb, Hi[2])), add(mul(ga, Hi[1]), mul(gb, Hi[3])),
          add(mul(gc, Hi[0]), mul(gd, Hi[2])), add(mul(gc, Hi[1]), mul(gd, Hi[3])))
    det = sub(mul(g1[0], g1[3]), mul(g1[1], g1[2])); assert det == (1, 0)
    inG3 = mod3eq(g1[0], (1, 0)) and mod3eq(g1[3], (1, 0)) and mod3eq(g1[1], (0, 0)) and mod3eq(g1[2], (0, 0))
    kap = 1 if g1[1] == (0, 0) else cubic_sym(g1[1], g1[3])     # kappa = (b1/d1)_3 on Gamma(3)
    return dict(a=a, c=c, c0=c0, r=r, dp=dp, bg=bg, sigma=sigma, u0=u0, inG3=inG3, kap=kap, vl=vl)

