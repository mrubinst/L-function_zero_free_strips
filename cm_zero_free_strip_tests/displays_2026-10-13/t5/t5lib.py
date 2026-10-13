"""Helpers for the t5 display tests (Appendix A of cm-eleven-twelfths-guide.tex, displays (A.1)-(A.41)).

Conventions are the guide's: omega = e^{2 pi i/3}; a + b omega is the integer pair (a, b); primary = 1 mod 3;
lambda = 1 + 2 omega = i sqrt 3; chi_pi(x) = (x/pi)_6; ce(z) = exp(2 pi i (z + conj z)); e(x) = ce(x/lambda)
= exp(4 pi i Im(x)/sqrt 3); alpha(n) = n/|n|.

Read-only imports: zw.py and fastgauss.py from tests_cm/part_B_lines_0550_0997 (exact Z[omega] arithmetic, primary
primes, cubic symbols, prime cubic Gauss sums).  Nothing in other suites is modified.

Theta coefficients: d_0, d_+, d_- are built from [QF, Theorem 23] ((5.7), (5.9), (5.10)) exactly as in the guide's
(A.4): t_0 = tau, t_+(l) = omega tau_2(omega l) ce(l), t_-(l) = omega^2 tau_1(omega^2 l) ce(l), d_s(l) = conj(t_s(-l)).
The Gauss sums g(r, c), c squarefree primary (possibly even), come from the prime Gauss sums and twisted
multiplicativity (fastgauss.py); for every table entry the index is kept EXACTLY as m = lambda^4 l in O together with
its decomposition l = u lambda^k c d^3.
Cusp data: own implementation of [QF, Section 5.5, item 3] ((5.18)): normalization of (a, c), delta, H, sigma,
g1 = g H^{-1} in Gamma(3), kappa(g1) = (c1/a1)_3 [QF (5.4)], cross-checked against (b1/d1)_3.
"""
import os, sys, math, cmath, time
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, '..', '..'))
PARTB = os.path.join(TESTS, 'part_B_lines_0550_0997')
if PARTB not in sys.path:
    sys.path.insert(0, PARTB)
import numpy as np
import fastgauss as fg
from zw import (mul, add, sub, neg, conj, N, cx, divides, div, pw, primary, UNITS, LAM, reduce_mod, canon_mod,
                powmod, isprime_int)

W = cmath.exp(2j*math.pi/3)
E9 = cmath.exp(2j*math.pi/9)
LAMc = complex(cx(LAM))
SQ3 = math.sqrt(3.0)
ZETA6 = [cmath.exp(1j*math.pi*k/3) for k in range(6)]     # UNITS[k] = zeta^k, zeta = 1 + omega

# ------------------------------------------------------------------------------------------------ reporting
_FAILS = []
_COUNT = [0]
def report(tid, label, ok, detail):
    """One line per check: [PASS]/[FAIL] tid (label) detail.  Failures are collected for the exit status."""
    _COUNT[0] += 1
    if os.environ.get('CM_TESTS_FORCE_FAIL') == '1' and _COUNT[0] == 1:
        ok = False; detail = 'CM_TESTS_FORCE_FAIL=1 (forced failure); ' + detail
    if not ok:
        _FAILS.append(tid)
    print(f"[{'PASS' if ok else 'FAIL'}] {tid} ({label}) {detail}", flush=True)
    return ok

def falsifier(tid, name, err, must_exceed):
    """A deliberately wrong variant: its error must exceed must_exceed."""
    worked = err > must_exceed
    if not worked:
        _FAILS.append(tid + ' falsifier ' + name)
    print(f"      falsifier for {tid} [{name}]: err={err:.3g} -> {'fails as required' if worked else 'DID NOT FAIL (bad falsifier)'}",
          flush=True)
    return worked

def finish(t0):
    print(f"time {time.time() - t0:.1f}s")
    if _FAILS:
        print(f"EXIT STATUS 1: {len(_FAILS)} failed: {_FAILS}")
        sys.exit(1)
    print(f"EXIT STATUS 0: all {_COUNT[0]} checks passed")
    sys.exit(0)

# ------------------------------------------------------------------------------------------------ arithmetic
def ce(z): return np.exp(4j*np.pi*np.real(z))
def e_add(z): return np.exp(4j*np.pi*np.imag(z)/SQ3)

def val(q, x):
    k = 0
    while x != (0, 0) and divides(q, x):
        x = div(x, q); k += 1
    return k

def gcd_O(x, y):
    while y != (0, 0):
        x, y = y, reduce_mod(x, y)
    return x

def egcd(x, y):
    """(g, u, v) with u x + v y = g."""
    u0, v0, u1, v1 = (1, 0), (0, 0), (0, 0), (1, 0)
    while y != (0, 0):
        n = N(y); q = mul(x, conj(y)); qq = (round(q[0]/n), round(q[1]/n))
        x, y = y, sub(x, mul(qq, y))
        u0, u1 = u1, sub(u0, mul(qq, u1)); v0, v1 = v1, sub(v0, mul(qq, v1))
    return x, u0, v0

def unit_inv(e):
    for f in UNITS:
        if mul(e, f) == (1, 0): return f

def inv_mod(x, m):
    g, u, _ = egcd(x, m)
    assert g in UNITS, ("not invertible", x, m)
    return reduce_mod(mul(u, unit_inv(g)), m)

def mod3eq(x, y):
    d = sub(x, y); return d[0] % 3 == 0 and d[1] % 3 == 0

def coords(z):
    """Exact coordinates (a, b) of a complex number z = a + b omega known to lie in O (rounded)."""
    b = z.imag/(SQ3/2); a = z.real + b/2
    return (int(round(a)), int(round(b)))

_PRIMES = fg.primary_primes(200000)            # primary primes of norm <= 2e5, including -2
def factor_prim(m):
    """Primary prime factors (with multiplicity) of m coprime to lambda; m up to a unit.  Returns (unit, [primes])."""
    out = []; r = m
    for p in _PRIMES:
        if N(r) == 1: break
        if N(p) > N(r):
            break
        while divides(p, r):
            out.append(p); r = div(r, p)
    if N(r) > 1:                                # leftover prime of norm > 2e5 (only if N(r) is prime)
        n = N(r); assert isprime_int(n) or (math.isqrt(n)**2 == n and isprime_int(math.isqrt(n))), (m, r)
        p = fg.prim_assoc(r); out.append(p); r = div(r, p)
    return r, out

def cubic_sym(x, a):
    """(x/a)_3 for a coprime to 3 (product over the prime factors of a; the unit part of a is ignored)."""
    if N(a) == 1: return 1
    _, ps = factor_prim(a)
    r = 1
    for p in ps:
        r *= fg.cubic_prime(x, p)
    return r

def sext_prime(x, pi):
    """chi_pi(x) = (x/pi)_6 for a primary prime pi coprime to 6; 0 if pi | x."""
    if divides(pi, x): return 0
    v = powmod(x, (N(pi) - 1)//6, pi)
    for k, u in enumerate(UNITS):
        if divides(pi, sub(v, u)): return ZETA6[k]
    raise ValueError((x, pi))

def lam3_euler(x):
    """(lambda/x)_3 for primary x coprime to 3 (possibly even), by Euler's criterion at each prime factor of x."""
    return cubic_sym(LAM, x)

def lam3_supp(x):
    """(lambda/x)_3 for primary x = a + b omega by the supplement: omega^{-b/3} (used in t_b3_prop_endtoend.py)."""
    return W**((-(x[1]//3)) % 3)

def om3(x):
    """(omega/x)_3 = omega^{(N(x)-1)/3} for x coprime to 3."""
    return W**(((N(x) - 1)//3) % 3)

def psi4(x):
    """A character of the primary classes mod 4 (coprime to 2): x = omega^k (1 + 2 y) mod 4, psi4 = (-1)^{y_0}."""
    for k in range(3):
        y = mul(x, pw((-1, -1), k))                    # x * omega^{-k} = x * omega^{2k}
        if (y[0] - 1) % 2 == 0 and y[1] % 2 == 0:
            y0 = ((y[0] - 1)//2) % 2
            return -1 if y0 else 1
    return 0

def prime_coprime6(x):
    return primary(x) and not divides((2, 0), x)

# ------------------------------------------------------------------------------------------------ theta tables
class ThetaTables:
    """Coefficients d_0, d_+, d_- on lambda^{-4} O with N(lambda^4 l) <= B, from [QF, Theorem 23]."""
    def __init__(self, B, verbose=True):
        t0 = time.time()
        self.B = B
        plist = fg.primary_primes(B)
        gp = {p: fg.g_prime(p) for p in plist}
        sq = fg.squarefree_list(B, plist)
        L2 = (-3, 0); wL2 = mul((0, 1), L2); w2L2 = mul((-1, -1), L2)
        self.cdata = []
        for c, f in sq:
            g = 1 + 0j; s1 = s2 = s3 = 1
            for p in f:
                g *= gp[p]; s1 *= fg.cubic_prime(L2, p); s2 *= fg.cubic_prime(wL2, p); s3 *= fg.cubic_prime(w2L2, p)
            for i in range(len(f)):
                for j in range(i + 1, len(f)):
                    g *= fg.pairsym(f[i], f[j])*fg.pairsym(f[j], f[i])
            # g(1,c), g(lambda^2, c), g(omega lambda^2, c), g(omega^2 lambda^2, c)   [QF (5.5), (5.24)]
            self.cdata.append((c, N(c), g, np.conj(s1)*g, np.conj(s2)*g, np.conj(s3)*g))
        R = int(round(B**(1/3))) + 2
        self.dl = [(x, y) for x in range(-2*R, 2*R + 1) for y in range(-2*R, 2*R + 1)
                   if primary((x, y)) and N((x, y))**3 <= B]
        self.gtab = {c: (g1, gl2, gwl2, gw2l2) for c, nc, g1, gl2, gwl2, gw2l2 in self.cdata}
        if verbose:
            print(f"   theta tables: {len(plist)} primes, {len(self.cdata)} squarefree c, B = {B}, {time.time()-t0:.1f}s", flush=True)
        self._cache = {}

    def table(self, sigma, Bm=None):
        """Arrays (m, l, d, k, c, b, u): m = lambda^4 l exactly (list of pairs), l complex, d = d_sigma(l), and the
        decomposition l = u lambda^k c b^3 (k = lambda-exponent of l), for all l with N(m) <= Bm."""
        Bm = self.B if Bm is None else Bm
        assert Bm <= self.B
        key = (sigma, Bm)
        if key in self._cache: return self._cache[key]
        M = []; Ls = []; D = []; K = []; Cs = []; Bs = []; Us = []
        for c, nc, g1, gl2, gwl2, gw2l2 in self.cdata:
            if nc > Bm: continue
            for d in self.dl:
                nd = N(d); base = nc*nd**3
                if base > Bm: continue
                absdc = abs(cx(d))/abs(cx(c))
                cd3 = mul(c, pw(d, 3))
                if sigma == '0':
                    for k in range(-3, 60):
                        if base*3**(k + 4) > Bm: break
                        r = k % 3
                        if r == 1: continue
                        if r == 2:                       # k = 3 nu - 4, nu >= 1, all six units
                            nu = (k + 4)//3
                            rows = [((1, 0), np.conj(gl2)), ((-1, 0), np.conj(gl2)), ((0, 1), E9**-1*np.conj(gwl2)),
                                    ((0, -1), E9**-1*np.conj(gwl2)), ((-1, -1), E9*np.conj(gw2l2)), ((1, 1), E9*np.conj(gw2l2))]
                            fac = absdc*3**(nu/2 + 2)
                        else:                            # k = 3 nu - 3, nu >= 0, units +-1
                            nu = (k + 3)//3
                            rows = [((1, 0), np.conj(g1)), ((-1, 0), np.conj(g1))]
                            fac = absdc*3**(nu/2 + 2.5)
                        for u, tv in rows:
                            m = mul(u, mul(pw(LAM, k + 4), cd3))
                            # d_0(l) = conj(tau(-l)) = conj(tau(l)) (tau is even: both signs are in the table)
                            M.append(m); D.append(np.conj(tv*fac)); K.append(k); Cs.append(c); Bs.append(d); Us.append(u)
                else:
                    if sigma == '+':
                        # d_+(l) = omega^2 conj(tau_2(-omega l)) ce(l); tau_2 at -omega^j lambda^-4 c d^3 [QF (5.10)]
                        rows = [((1, 0), 9*W*W*np.conj(gl2)), ((0, 1), 9*E9**-1*np.conj(gwl2)), ((-1, -1), 9*W*E9*np.conj(gw2l2))]
                        for up, t2 in rows:
                            u = mul(up, (-1, -1))                       # l = omega^{j+2} lambda^-4 c d^3
                            m = mul(u, cd3); lc = complex(cx(m))/9
                            M.append(m); D.append(W*W*np.conj(t2*absdc)*complex(ce(lc))); K.append(-4); Cs.append(c); Bs.append(d); Us.append(u)
                    else:
                        # d_-(l) = omega conj(tau_1(-omega^2 l)) ce(l); tau_1 at omega^j lambda^-4 c d^3 [QF (5.9)]
                        rows = [((1, 0), 9*W*np.conj(gl2)), ((0, 1), 9*E9**-1*W*W*np.conj(gwl2)), ((-1, -1), 9*E9*np.conj(gw2l2))]
                        for up, t1 in rows:
                            u = neg(mul(up, (0, 1)))                    # l = -omega^{j+1} lambda^-4 c d^3
                            m = mul(u, cd3); lc = complex(cx(m))/9
                            M.append(m); D.append(W*np.conj(t1*absdc)*complex(ce(lc))); K.append(-4); Cs.append(c); Bs.append(d); Us.append(u)
        Ls = np.array([complex(cx(m))/9 for m in M])
        out = dict(m=M, l=Ls, d=np.array(D), k=np.array(K), c=Cs, b=Bs, u=Us)
        self._cache[key] = out
        return out

# ------------------------------------------------------------------------------------------------ cusp data
def normalize(a, c):
    """Common unit so that a = 1 mod 3 if lambda | c, and c = 1 mod 3 otherwise (Step 2(a); [QF Sec. 5.5])."""
    for u in UNITS:
        aa, cc = mul(u, a), mul(u, c)
        if divides(LAM, c):
            if mod3eq(aa, (1, 0)): return aa, cc
        elif mod3eq(cc, (1, 0)): return aa, cc
    raise ValueError((a, c))

def lowest_terms(num, den):
    """num/den in lowest terms (a, c), normalized; 0 -> (0, 1)."""
    if num == (0, 0): return (0, 0), (1, 0)
    g = gcd_O(num, den)
    return normalize(div(num, g), div(den, g))

def matmul(A, B):
    (a, b), (c, d) = A; (e, f), (g, h) = B
    return ((add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h))), (add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h))))

def matinv(A):
    (a, b), (c, d) = A
    return ((d, neg(b)), (neg(c), a))

def kubota(g):
    """kappa(g) for g in Gamma(3): (c/a)_3 if c != 0, else 1 [QF (5.4)]."""
    (a, b), (c, d) = g
    return 1 if c == (0, 0) else cubic_sym(c, a)

def kubota_bd(g):
    """The alternative form (b/d)_3 on Gamma(3) [QF, remark after Theorem 22]."""
    (a, b), (c, d) = g
    return 1 if b == (0, 0) else cubic_sym(b, d)

def in_gamma3(g):
    (a, b), (c, d) = g
    return mod3eq(a, (1, 0)) and mod3eq(d, (1, 0)) and mod3eq(b, (0, 0)) and mod3eq(c, (0, 0))

def cusp_data(a, c):
    """Cusp data of a/c ((a, c) coprime and normalized) following [QF, Sec. 5.5, item 3]:
    delta with a delta = 1 mod 3c (lambda | c) or a delta = 1 mod c and delta = 0 mod 3 (lambda not dividing c);
    g = (a, (a delta - 1)/c; c, delta); H of (5.18); sigma; g1 = g H^{-1} in Gamma(3); kappa(g1)."""
    vl = val(LAM, c)
    if vl >= 1:
        mod = mul((3, 0), c); delta = inv_mod(a, mod)
    else:
        # delta = a^{-1} mod c and delta = 0 mod 3 (CRT; c is coprime to 3)
        if N(c) == 1:
            delta = (0, 0)
        else:
            ai = inv_mod(a, c)
            t = mul(ai, inv_mod((3, 0), c))              # delta = 3 t, 3 t = a^{-1} mod c
            delta = mul((3, 0), reduce_mod(t, c))
    t = sub(mul(a, delta), (1, 0)); assert divides(c, t); bg = div(t, c)
    g = ((a, bg), (c, delta))
    one, zero = (1, 0), (0, 0)
    if vl >= 2:
        H = ((one, zero), (zero, one)); sigma = '0'
    elif vl == 1:
        u = LAM if mod3eq(c, LAM) else neg(LAM); assert mod3eq(c, u)
        H = ((one, zero), (u, one)); sigma = '-' if u == LAM else '+'
    else:
        u = (a[0] % 3, a[1] % 3)
        H = ((u, (-1, 0)), (one, zero)); sigma = {0: '0', 1: '-', 2: '+'}[u[1]]
    g1 = matmul(g, matinv(H))
    (A1, B1), (C1, D1) = g1
    assert sub(mul(A1, D1), mul(B1, C1)) == (1, 0)
    return dict(a=a, c=c, delta=delta, b=bg, g=g, H=H, sigma=sigma, g1=g1, inG3=in_gamma3(g1),
                kappa=kubota(g1) if in_gamma3(g1) else None, kappa_bd=kubota_bd(g1) if in_gamma3(g1) else None, vl=vl)

GAMMA = {'0': ((1, 0), (0, 0), (0, 0), (1, 0)), '+': ((1, 0), (0, 0), (0, 1), (1, 0)), '-': ((1, 0), (0, 0), (-1, -1), (1, 0))}

def act(g, z, v):
    """Action of g = (a, b, c, d) (complex entries) on H^3 [QF (5.1)]:
    (z, v) -> (((az+b) conj(cz+d) + a conj(c) v^2)/D, v/D), D = |cz+d|^2 + |c|^2 v^2."""
    a, b, c, d = g
    den = abs(c*z + d)**2 + abs(c)**2*v*v
    return (((a*z + b)*np.conj(c*z + d) + a*np.conj(c)*v*v)/den, v/den)

def mat_c(g):
    (a, b), (c, d) = g
    return tuple(complex(cx(x)) for x in (a, b, c, d))
