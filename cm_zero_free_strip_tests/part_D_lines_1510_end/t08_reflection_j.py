"""Part D, test 08: end-to-end single-cusp theta reflection with j zbar-derivatives, j = 1, 2, 3.

This tests, on actual theta coefficients, the general-j claims of lines 1585-1600 (reflected Gamma factors,
cusp factor conj(alpha(c))^{2j}, dual angular factor alpha(mu)^j, kernel R_j, argument N(mu)X/N(c)^2 independent of
j), and for j = 2 the structure of eq. (reflection), line 1681 (alpha(l)^2, V# with kernel R_2, argument
N(l)X/N(c)^2, coefficients d_sigma from the appendix tables).

For a cusp a/C (a, C coprime) the identity tested is
  LHS_j(X) = sum_{0 != l in lambda^-3 O} d_0(l) ce(l a/C) conj(alpha(l))^j N(l)^{-1/2} V(N(l)/X)
  RHS_j(X) = (-1)^j conj(kappa(g1)) conj(alpha(C))^{2j}
             * sum_{0 != l in lambda^-4 O} d_sigma(l) alpha(l)^j ce(-delta' l/C) N(l)^{-1/2} V#_j((2pi)^4 N(l) X/N(C)^2),
  V#_j(y) = (1/2 pi i) int_{Re t=0} V^(-t) R_j(t) y^{-t} dt,   R_j as in line 1589,
which follows from the chain rule (L1571) and the two Mellin transforms (L1579, L1585) applied to one translate
of conj(theta) (sign (-1)^j from (-(cv)^{-2})^j; the i^j and 1/4 cancel).  The cusp normalizer H, the matrix g,
g1 = g H^{-1} = I mod 3, and kappa(g1) = (C1/A1)_3 are as in tests_guide/test_sec7_reflection_j1.py (OpenAI's
11/12 paper, Appendix A); the coefficient tables d_0, d_+, d_- are built here directly from the appendix formulas
of the CM paper (t_+(l) = omega tau_2(omega l) ce(l), t_-(l) = omega^2 tau_1(omega^2 l) ce(l), d_s(l) = conj(t_s(-l))).
Weight V(y) = exp(-(log y)^2/(2 s^2)), s = 0.5 (V^ Gaussian).

Controls (each must fail, by an amount of size ~1, at some cusp; at real C, e.g., conj(alpha(C))^2 = 1 cannot be detected): kernel R_{j+1} in place of R_j; conj(alpha(C))^{2j+2};
sign (-1)^{j+1}.  For j = 1 this reproduces the guide's test 7.6.
PASS: relative error < 1e-6 at every nondegenerate cusp for each j, |RHS| < 1e-8 where LHS vanishes, controls > 1e-2.
Runtime: a few minutes (tables with B = 80000).
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import report, info
from fastgauss import *
import appendix_theta as AT
from scipy.special import loggamma
t00 = time.time()
LAMp = (1, 2); LAMc = cx(LAMp)
s_w = 0.5; K = 6.5; YCUT = 4e5
def V(y): return np.exp(-np.log(y)**2/(2*s_w*s_w))
def Rj(j, t):
    a = (j + 1)/2
    return np.exp(loggamma(a + 1/6 + t) + loggamma(a - 1/6 + t) - loggamma(a + 1/6 - t) - loggamma(a - 1/6 - t))
TAU = np.arange(-14/s_w, 14/s_w + 1e-9, 0.05); TT = 1j*TAU
def FVj(j): return np.sqrt(2*np.pi)*s_w*np.exp(s_w*s_w*TT*TT/2)*Rj(j, TT)*0.05/(2*np.pi)
FV = {j: FVj(j) for j in range(1, 5)}
def Vnat(y, j):
    y = np.asarray(y, float); out = np.empty(len(y), complex)
    for i in range(0, len(y), 4000):
        ly = np.log(y[i:i + 4000]); out[i:i + 4000] = np.exp(-np.outer(ly, TT))@FV[j]
    return out
def ce(z): return np.exp(4j*np.pi*np.real(z))

# ---------------------------------------------------------------- cusp data (as in the guide's j=1 script)
def mod3eq(x, y): d = sub(x, y); return d[0] % 3 == 0 and d[1] % 3 == 0
def vlam(x):
    k = 0
    while divides(LAMp, x): x = div(x, LAMp); k += 1
    return k
def cubic_general(x, a, plist):
    assert primary(a)
    r = 1; m = a
    for p in plist:
        if N(m) == 1: break
        while divides(p, m):
            r *= cubic_prime(x, p); m = div(m, p)
    assert N(m) == 1, (a, m)
    return r
def matmul(A, B):
    (a, b), (c, d) = A; (e, f), (g, h) = B
    return ((add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h))), (add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h))))
def inv(A):
    (a, b), (c, d) = A; return ((d, (-b[0], -b[1])), ((-c[0], -c[1]), a))
def cusp_data(a, C, plist):
    v = vlam(C); one = (1, 0); zero = (0, 0)
    if v >= 2:
        assert mod3eq(a, one); H = ((one, zero), (zero, one)); sigma = '0'
    elif v == 1:
        assert mod3eq(a, one)
        u0 = LAMp if mod3eq(C, LAMp) else (-1, -2); assert mod3eq(C, u0)
        H = ((one, zero), (u0, one)); sigma = '-' if u0 == LAMp else '+'
    else:
        assert primary(C)
        u0 = (a[0] % 3, a[1] % 3)
        eps = {0: 0, 1: 1, 2: -1}[u0[1]]
        H = ((u0, (-1, 0)), (one, zero)); sigma = {0: '0', 1: '-', -1: '+'}[eps]
    found = {}
    R = 14
    for x in range(-R, R + 1):
        for y in range(-R, R + 1):
            dp = (x, y); t = sub(mul(a, dp), one)
            if not divides(C, t): continue
            bg = div(t, C); g = ((a, bg), (C, dp)); g1 = matmul(g, inv(H))
            (A1, B1), (C1, D1) = g1
            if mod3eq(A1, one) and mod3eq(D1, one) and mod3eq(B1, zero) and mod3eq(C1, zero):
                if N(A1) > 4000: continue
                kap = 1 if C1 == (0, 0) else cubic_general(C1, A1, plist)
                key = (round(complex(kap).real, 6), round(complex(kap).imag, 6))
                if key not in found or N(dp) < N(found[key][0]): found[key] = (dp, kap)
    return [dict(sigma=sigma, delta=dp, kappa=kap, v=v) for dp, kap in found.values()]

# ---------------------------------------------------------------- tables from the appendix formulas
B = int(sys.argv[1]) if len(sys.argv) > 1 else 80000
plist = primary_primes(B)
gp_ = {p: g_prime(p) for p in plist}
sq = squarefree_list(B, plist)
L2 = (-3, 0); wL2 = mul((0, 1), L2); w2L2 = mul((-1, -1), L2)
# check of the fast Gauss sums against the brute-force g(r,c) of the appendix definition (small c)
fast_chk = 0
cdata = []
for c, f in sq:
    g = 1 + 0j; s1 = s2 = s3 = 1
    for p in f:
        g *= gp_[p]; s1 *= cubic_prime(L2, p); s2 *= cubic_prime(wL2, p); s3 *= cubic_prime(w2L2, p)
    for i in range(len(f)):
        for j in range(i + 1, len(f)): g *= pairsym(f[i], f[j])*pairsym(f[j], f[i])
    gr = {'1': g, 'L2': np.conj(s1)*g, 'wL2': np.conj(s2)*g, 'w2L2': np.conj(s3)*g}   # g(r,c) = conj((r/c)_3) g(1,c)
    if N(c) <= 120:
        for key, r in [('1', (1, 0)), ('L2', L2), ('wL2', wL2), ('w2L2', w2L2)]:
            fast_chk = max(fast_chk, abs(gr[key] - AT.g(r, c)))
    cdata.append((cx(c), N(c), gr))
info(f"tables: {len(cdata)} squarefree c, fast vs brute-force g(r,c) max diff {fast_chk:.1e}, {time.time() - t00:.0f}s")
Dmax = int(round(B**(1/3))) + 2
bl = [(cx((x, y)), N((x, y))) for x in range(-2*Dmax, 2*Dmax + 1) for y in range(-2*Dmax, 2*Dmax + 1)
      if primary((x, y)) and N((x, y))**3 <= B]
Z9 = cmath.exp(2j*math.pi/9)

def tau_entries(Nmax):
    """All (mu, tau(mu)) with mu in lambda^-3 O (tau's support), N(mu) <= Nmax, from the table of line 1711."""
    L = []; Cf = []
    for cc, nc, gr in cdata:
        if nc > 27*Nmax: continue
        for bb, nb in bl:
            base = nc*nb**3
            if base > 27*Nmax: continue
            f = abs(bb)/abs(cc)
            for k in range(-3, 40):
                if base*3.0**k > Nmax*(1 + 1e-12): break
                r = k % 3
                if r == 1: continue
                lk = LAMc**k*cc*bb**3
                if r == 2:
                    nu = (k + 4)//3; A = 3**(nu/2 + 2)*f
                    vals = [(1, A*np.conj(gr['L2'])), (-1, A*np.conj(gr['L2'])),
                            (W, A*Z9**-1*np.conj(gr['wL2'])), (-W, A*Z9**-1*np.conj(gr['wL2'])),
                            (W*W, A*Z9*np.conj(gr['w2L2'])), (-W*W, A*Z9*np.conj(gr['w2L2']))]
                else:
                    nu = (k + 3)//3; A = 3**(nu/2 + 5/2)*f
                    vals = [(1, A*np.conj(gr['1'])), (-1, A*np.conj(gr['1']))]
                for u, t in vals:
                    L.append(u*lk); Cf.append(t)
    return np.array(L), np.array(Cf)

def table_d0(Nmax):
    mu, t = tau_entries(Nmax)
    return -mu, np.conj(t)                  # d_0(l) = conj(tau(-l)):  l = -mu

def table_dpm(sign, Nmax):
    L = []; Cf = []
    for cc, nc, gr in cdata:
        if nc > 81*Nmax: continue
        for bb, nb in bl:
            base = nc*nb**3
            if base > 81*Nmax: continue
            f = 9*abs(bb)/abs(cc); m0 = LAMc**-4*cc*bb**3
            if sign == '+':
                # tau_2(-omega^j lam^-4 c b^3); t_+(l) = omega tau_2(omega l) ce(l) with omega l = -omega^j m0
                for wj, t2 in [(1, f*W*W*np.conj(gr['L2'])), (W, f*Z9**-1*np.conj(gr['wL2'])), (W*W, f*W*Z9*np.conj(gr['w2L2']))]:
                    l = -wj*m0/W; tp = W*t2*ce(l)
                    L.append(-l); Cf.append(np.conj(tp))          # d_+(-l) = conj(t_+(l))
            else:
                # tau_1(omega^j lam^-4 c b^3); t_-(l) = omega^2 tau_1(omega^2 l) ce(l) with omega^2 l = omega^j m0
                for wj, t1 in [(1, f*W*np.conj(gr['L2'])), (W, f*Z9**-1*W*W*np.conj(gr['wL2'])), (W*W, f*Z9*np.conj(gr['w2L2']))]:
                    l = wj*m0/(W*W); tm = W*W*t1*ce(l)
                    L.append(-l); Cf.append(np.conj(tm))
    return np.array(L), np.array(Cf)
info(f"tables ready {time.time() - t00:.0f}s")

cusps = [((1, 0), (-3, 0)), ((4, 0), (-3, 0)), ((1, 0), pw(LAMp, 3)), ((7, 0), (-6, 0)), ((-2, -3), (-6, 0)),
         ((1, 0), LAMp), ((1, 0), (-1, -2)), ((1, 0), mul((-2, 0), LAMp)), ((7, 0), mul((2, 0), LAMp)), ((1, 0), mul(LAMp, (-2, -3))),
         ((0, 0), (1, 0)), ((0, 1), (1, 0)), ((0, -1), (1, 0)), ((1, 0), (-2, 0)), ((0, 1), (-2, 0)), ((0, -1), (-2, 0)),
         ((1, 1), (-2, -3)), ((2, 0), (-2, -3)), ((0, 2), (-2, -3))]
JS = (1, 2, 3)
worst = {j: 0 for j in JS}; ctrl = {j: [0.0, 0.0, 0.0] for j in JS}; zero_cases = []; ncfg = 0; sig_seen = set()
for a, C in cusps:
    for cd in cusp_data(a, C, plist):
        NC = N(C); Cc = cx(C); ac = cx(a)/Cc
        X = 3.16*NC if cd['sigma'] == '0' else 5.46*NC
        NL = X*math.exp(K*s_w); NR = YCUT*NC**2/((2*math.pi)**4*X)
        if 27*NL > B or (27 if cd['sigma'] == '0' else 81)*NR > B:
            info(f"skip (B too small) {a} {C}"); continue
        Ll, Cl = table_d0(NL)
        Lr, Cr = table_d0(NR) if cd['sigma'] == '0' else table_dpm(cd['sigma'], NR)
        dl_ = cx(cd['delta']); nl = np.abs(Lr)**2; kap = np.conj(complex(cd['kappa']))
        ncfg += 1; sig_seen.add(cd['sigma'])
        for j in JS:
            lhs = np.sum(Cl*ce(Ll*ac)*np.conj(Ll/np.abs(Ll))**j*np.abs(Ll)**-1*V(np.abs(Ll)**2/X))
            base = Cr*(Lr/np.abs(Lr))**j*ce(-dl_*Lr/Cc)*nl**-0.5
            y = (2*np.pi)**4*nl*X/NC**2
            core = kap*np.conj(Cc/abs(Cc))**(2*j)*np.sum(base*Vnat(y, j))
            rhs = (-1)**j*core
            if abs(lhs) < 1e-8:
                zero_cases.append(abs(rhs)); continue
            err = abs(lhs - rhs)/abs(lhs); worst[j] = max(worst[j], err)
            c1 = abs(lhs - (-1)**j*kap*np.conj(Cc/abs(Cc))**(2*j)*np.sum(base*Vnat(y, j + 1)))/abs(lhs)
            c2 = abs(lhs - (-1)**j*kap*np.conj(Cc/abs(Cc))**(2*j + 2)*np.sum(base*Vnat(y, j)))/abs(lhs)
            c3 = abs(lhs + rhs)/abs(lhs)
            ctrl[j] = [max(ctrl[j][0], c1), max(ctrl[j][1], c2), max(ctrl[j][2], c3)]
        print(f"   a={a} C={C} sigma={cd['sigma']} kappa={complex(cd['kappa']):.3f}: rel.err j=1,2,3 so far "
              + ", ".join(f"{worst[j]:.1e}" for j in JS), flush=True)
info(f"{ncfg} cusp configurations (sigma types {sorted(sig_seen)}), degenerate (LHS=0) cases {len(zero_cases)}, max |RHS| there "
     f"{max(zero_cases) if zero_cases else 0:.1e}, {time.time() - t00:.0f}s")
zok = (not zero_cases) or max(zero_cases) < 1e-8
for j in JS:
    report(f"D08.{j}", 1585 if j != 2 else 1681, f"reflection j={j}",
           f"single-cusp theta reflection with {j} derivative(s): sign (-1)^j, conj(alpha(C))^(2j), alpha(l)^j, kernel R_{j}, argument (2pi)^4 N(l)X/N(C)^2"
           + (" [= structure of eq. (reflection)]" if j == 2 else ""),
           worst[j] < 1e-6 and zok and fast_chk < 1e-9 and min(ctrl[j]) > 1e-2,
           f"worst rel err {worst[j]:.1e} over {ncfg} cusps; controls (largest error over cusps; a wrong variant must fail somewhere): R_{j+1} {ctrl[j][0]:.2f}, "
           f"conj(alpha(C))^{2*j+2} {ctrl[j][1]:.2f}, sign flipped {ctrl[j][2]:.2f}")
