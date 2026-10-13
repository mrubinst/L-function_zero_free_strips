"""r2 test T2: Step 2 and the Step 6 grouping of Appendix A (Proposition prop:b3-reflection), with L = 72.

What is tested (exact arithmetic in Z[omega] with Python big integers; residue symbols by factoring norms):
For a group (h0, A) of translates z_h = lambda^2 (h0/72 + sum_{pi in A} h_pi/pi), with the representatives
h_pi = 0 mod M^2, M = lambda^12 72^4 (Step 1), and delta' chosen by the congruences of [OAI1112, App. A.2]
(a delta' = 1 mod q^{v_q(M c0)} for q | c0; delta' = 0 mod q^{v_q(M)} for q | M, q not dividing c0;
a delta' = 1 mod r), we check
  (a)  c = c0 r is the same for the whole group, c0 | 8 lambda^2 up to a unit, and a = sigma_pi h_pi mod pi
       with sigma_pi = lambda^2 c/pi;
  (b)  g = (a b_g; c delta') is in SL_2(O), g1 = g H^{-1} is in Gamma(3) (H as in [QF (5.18)]), and the cusp
       label sigma and delta' mod lambda^3 c0 are the same for the whole group;
  (c)  kappa(g1) = kappa0 prod_{pi in A} chi_pi(a)^2 with kappa0 constant on the group, and kappa0 agrees with
       the three-case formula [OAI1112 (A.10)];
  (d)  (b3-crt): ce(-delta' l/c) = psi(lambda^4 l) prod e(eps_pi h_pi^{-1} lambda^4 l/pi), psi(x) = e(beta x/(lambda^3 c0)),
       beta = -delta'_0 r^{-1} mod lambda^3 c0, eps_pi = -((lambda^3 c/pi) sigma_pi)^{-1} mod pi  (exact, as an
       element of O);
  (e)  the grouped local identity used in Step 6: for every exponent vector (j_pi) and many x in O,
         sum_{h_pi != 0} prod_pi C_{pi,j_pi}(h_pi) * conj(kappa(g1)) * ce(-delta' l/c)        (x = lambda^4 l)
           = conj(kappa0) psi(x) prod_pi chi_pi(sigma_pi)^{-2} varpi_{pi,j} B_{pi,j}(x),
       with varpi_{pi,j} as in [OAI1112 (A.12)] (|varpi| = 1) and B_{pi,j} as in (8.5).
Falsifiers (must fail): the Kubota shift chi_pi(a)^{+2} instead of chi_pi(a)^{-2} in (e) (i.e. kappa instead of
conj(kappa)); B_{pi,j} with exponent -j instead of -j-2; dropping the factor psi.
Kappa: kappa(g1) = (b1/d1)_3 on Gamma(3) [QF, remark after Theorem 22], evaluated by factoring N(d1) (PARI/GP).
Usage: python3 t2_grouping_cusp_data.py          (about 1-2 minutes)
"""
import cmath, math, os, random, subprocess, sys, time
t0 = time.time()
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
print("primes:", PI7, N(PI7), PI13, N(PI13), PI_INERT, N(PI_INERT))
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

def A10_formula(d):
    a, c0, r, bg, u0, vl = d['a'], d['c0'], d['r'], d['bg'], d['u0'], d['vl']
    if vl >= 2:
        return sym_c0_over_a(c0, a)
    if vl == 1:
        m = sub(a, mul(u0, bg)); assert primary(m)
        t1 = unit3((-1, 0), m)*lam3(m)**1 if u0 == LAM else lam3(m)      # (-u0/m)_3 with u0 = +-lambda
        c0u = div(c0, u0)
        return t1*sym_c0_over_a(c0u, a)
    # (a/c0)_3 with c0 = unit * 2^i
    i = val(TWO, c0); return two3(a)**i

# ----------------------------------------------------------------- groups
ok = True; nfail = 0
FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']   # self-test (added 2026-10-12): first check counts as failed
def check(cond, msg):
    global ok, nfail
    if FORCE[0]:
        FORCE[0] = False; cond = False; msg = "[CM_TESTS_FORCE_FAIL=1] " + msg
    if not cond:
        ok = False; nfail += 1; print("   FAIL:", msg)
H0S = [(0, 0), (1, 0), (3, 0), (1, 2), (2, 1), (5, 7), (9, 3), (4, 0), (6, 3), (12, 0), (2, 5), (8, 4)]
H0S += [(random.randrange(72), random.randrange(72)) for _ in range(28)]   # random classes mod 72
GROUPS = [(h0, [PI7]) for h0 in H0S] + [(h0, [PI_INERT]) for h0 in [(1, 0), (1, 2), (3, 0)]] + \
         [(h0, [PI7, PI13]) for h0 in [(1, 0), (1, 2), (3, 0)]]
fals = {'kappa_unconj': 0.0, 'B_exp_-j': 0.0, 'no_psi': 0.0}
worst = 0.0; ntr = 0; types = {}
for h0, A in GROUPS:
    res_lists = [nonzero_res(p) for p in A]
    combos = [[]]
    for rl in res_lists: combos = [cmb + [h] for cmb in combos for h in rl]
    data = []
    for cmb in combos:
        hs = []
        for p, h in zip(A, cmb):
            tt = mul(h, inv_mod(M2, p)); hs.append(mul(M2, red(tt, p)))
            assert divides(p, sub(hs[-1], h)) and divides(M2, hs[-1])
        d = translate_data(h0, A, hs); d['hs'] = hs; d['hres'] = cmb; data.append(d); ntr += 1
    c = data[0]['c']; c0 = data[0]['c0']; r = data[0]['r']; D0 = mul(mul(LAM, mul(LAM, LAM)), c0)
    types[(data[0]['sigma'], data[0]['vl'], N(c0))] = types.get((data[0]['sigma'], data[0]['vl'], N(c0)), 0) + 1
    # (a), (b)
    check(all(d['c'] == c for d in data), f"c not constant, group {h0},{A}")
    check(all(d['inG3'] for d in data), f"g1 not in Gamma(3), group {h0},{A}")
    check(len({d['sigma'] for d in data}) == 1, f"sigma not constant, group {h0},{A}")
    check(len({canon(d['dp'], D0) for d in data}) == 1, f"delta' mod lambda^3 c0 not constant, group {h0},{A}")
    sig = {p: div(mul(mul(LAM, LAM), c), p) for p in A}
    for d in data:
        for p, h in zip(A, d['hs']):
            check(divides(p, sub(d['a'], mul(sig[p], h))), "a != sigma_pi h_pi mod pi")
    # (c) kappa0
    k0s = []
    for d in data:
        f = 1
        for p in A: f *= cubic_sym_prime(d['a'], p)      # chi_pi(a)^2 = (a/pi)_3
        k0s.append(d['kap']/f)
    check(max(abs(k - k0s[0]) for k in k0s) < 1e-9, f"kappa0 not constant, group {h0},{A}: {set(round(k.real,3)+1j*round(k.imag,3) for k in k0s)}")
    kap0 = k0s[0]
    fA10 = A10_formula(data[0])
    check(abs(fA10 - kap0) < 1e-9, f"(A.10) formula {fA10} != kappa0 {kap0}, group {h0},{A}, vl={data[0]['vl']}")
    # (d) CRT for the additive character
    dp0 = data[0]['dp']; rinv = inv_mod(r, D0); beta = red(neg(mul(dp0, rinv)), D0)
    eps = {p: red(neg(inv_mod(mul(div(mul(pw(LAM, 3), c), p), sig[p]), p)), p) for p in A}
    l3c = mul(pw(LAM, 3), c)
    for d in data:
        # element  -delta'/(lambda^3 c) - beta/(lambda^3 c0) - sum eps h^{-1}/pi  must lie in O
        den = mul(l3c, (1, 0)); num = neg(d['dp'])
        num = sub(num, mul(beta, r))                                   # beta/(lambda^3 c0) = beta r/(lambda^3 c)
        for p, h in zip(A, d['hs']):
            num = sub(num, mul(mul(eps[p], inv_mod(h, p)), div(l3c, p)))
        check(inO(num, den), f"(b3-crt) fails, group {h0},{A}")
    # (e) grouped local identity, all exponent vectors
    xs = [(random.randint(-10**6, 10**6), random.randint(-10**6, 10**6)) for _ in range(40)]
    xs += [mul(A[0], (random.randint(-99, 99), random.randint(-99, 99))) for _ in range(8)]
    jvecs = [[j] for j in range(6)] if len(A) == 1 else [[j1, j2] for j1 in range(6) for j2 in (0, 1, 4)]
    for jv in jvecs:
        Cvals = {}
        for p, j in zip(A, jv):
            for h in nonzero_res(p): Cvals[(p, h)] = Ccoef(p, j, h)
        varpi = {}
        for p, j in zip(A, jv):
            if j == 4: varpi[p] = gauss(p, 4)
            elif j == 0: varpi[p] = -gauss(p, 2)*sext_sym_prime(eps[p], p)**-2
            else: varpi[p] = sext_sym_prime((-1, 0), p)**j*gauss(p, j)*gauss(p, (j+2) % 6)*sext_sym_prime(eps[p], p)**(-j-2)
        for x in xs:
            lhs = 0; lhs_k = 0
            for d in data:
                coef = 1
                for p, h in zip(A, d['hres']): coef *= Cvals[(p, h)]
                ph = e_K(neg(mul(d['dp'], x)), l3c)
                kp = complex(d['kap'])
                lhs += coef*kp.conjugate()*ph
                lhs_k += coef*kp*ph
            psi = e_K(mul(beta, x), D0)
            rhs = kap0.conjugate()*psi; rhs_B = kap0.conjugate()*psi; rhs_np = kap0.conjugate()
            for p, j in zip(A, jv):
                fac = sext_sym_prime(sig[p], p)**-2*varpi[p]
                rhs *= fac*Bloc(p, j, x); rhs_np *= fac*Bloc(p, j, x)
                ch = sext_sym_prime(x, p); bwrong = (ch**(-j) if ch != 0 else 0) if j not in (0, 4) else Bloc(p, j, x)
                rhs_B *= fac*bwrong
            err = abs(lhs - rhs); worst = max(worst, err)
            fals['kappa_unconj'] = max(fals['kappa_unconj'], abs(lhs_k - rhs))
            fals['B_exp_-j'] = max(fals['B_exp_-j'], abs(lhs - rhs_B))
            fals['no_psi'] = max(fals['no_psi'], abs(lhs - rhs_np))
    print(f"group h0={h0} A={[N(p) for p in A]}: c={c} N(c)={N(c)} c0={c0} sigma={data[0]['sigma']} v_lam(c)={data[0]['vl']} "
          f"kappa0={complex(kap0):.3f} (A.10)={complex(fA10):.3f} translates={len(data)}  t={time.time()-t0:.0f}s")
print(f"translates {ntr}; cusp types (sigma, v_lam(c), N(c0)): {dict(sorted(types.items()))}")
print(f"(e) worst abs error over all groups, exponent vectors and x: {worst:.2e}")
print("falsifiers (max abs error, must be >= 0.1):", {k: f"{v:.2f}" for k, v in fals.items()})
passed = ok and worst < 1e-9 and min(fals.values()) > 0.1
print("PASS" if passed else f"FAIL ({nfail} failed checks)", f"time {time.time()-t0:.0f}s")
sys.exit(0 if passed else 1)   # exit status added 2026-10-12
