"""b3 end-to-end test of Proposition prop:b3-reflection for P = {} and three characters xi_0 (before grouping).

Twists (Psi(n) = xi_0(n) on primary n coprime to 6; phi(x) = (lambda/x)_3 xi_0(x)):
  '1'   : xi_0 = 1,                                   phi = (lambda/x)_3   (period 18; cusps with lambda^2 | c, sigma = 0)
  'A'   : xi_0 = conj((lambda/.)_3),                   phi = 1              (period 6;  c | 2, sigma = 0, +, -)
  'B'   : xi_0 = conj((lambda/.)_3) (omega/.)_3,        phi = (omega/x)_3    (period 2 lambda^3; c | 2 lambda, all sigma)
All three xi_0 are characters of the group of primary classes modulo 72.

Identity tested (Steps 1-6 of the proof, with the finite Fourier expansion of phi modulo L = 18):
  T(X;Psi) = -(1/81) sum_{h0 mod 18} phihat(h0) conj(kappa(g1)) conj(alpha(c))^4
             * sum_{0 != l in lambda^-4 O} d_sigma(l) alpha(l)^2 N(l)^{-1/2} ce(-delta l/c) V#(N(l) X/N(c)^2),
where Psi(n) = 1 on primary n coprime to 6, T is eq:T (with V(y) = sqrt(y) W(y) replaced by the log-Gaussian
V(y) = exp(-(log y)^2/(2 s^2)), whose Mellin transform is entire and Gaussian on vertical lines), phi(x) =
(lambda/x)_3 on primary x coprime to 6 (period 18), z_h = lambda^2 h0/18 = a/c in lowest terms, and
(sigma, delta, kappa) are the cusp data of a/c (any g = (a *; c delta) with g H^{-1} = I mod 3, H as in
OAI1112 App. A.2).  The left side is computed directly from eq:T with brute-force-checked gamma_2.
Controls (must fail): conj(alpha(c))^2 instead of ^4; zeta prefactor -i/81 instead of -1/81; one-derivative kernel.
Helpers fastgauss.py, zw.py are copies of tests_cm/part_B_lines_0550_0997 (unchanged).
Usage: python3 t_b3_prop_endtoend.py      (about 3-5 minutes)
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastgauss import *
from zw import canon_mod, residues, reduce_mod
from scipy.special import loggamma
import numpy as np
t0 = time.time()
LAMp = (1, 2); LAMc = cx(LAMp)
s_w = 0.5; KCUT = 7.5; YCUT = 4e5
X = 2400.0; B = 110000
def V(y): return np.exp(-np.log(y)**2/(2*s_w*s_w))
def Q2(t): return np.exp(loggamma(4/3+t)+loggamma(5/3+t)-loggamma(4/3-t)-loggamma(5/3-t))
def Q1(t): return np.exp(loggamma(7/6+t)+loggamma(5/6+t)-loggamma(7/6-t)-loggamma(5/6-t))
TAU = np.arange(-14/s_w, 14/s_w+1e-9, 0.05); TT = 1j*TAU
FV2 = np.sqrt(2*np.pi)*s_w*np.exp(s_w*s_w*TT*TT/2)*Q2(TT)*0.05/(2*np.pi)
FV1 = np.sqrt(2*np.pi)*s_w*np.exp(s_w*s_w*TT*TT/2)*Q1(TT)*0.05/(2*np.pi)
def Vnat(y, FVv):
    y = np.asarray(y, float); out = np.empty(len(y), complex)
    for i in range(0, len(y), 4000):
        out[i:i+4000] = np.exp(-np.outer(np.log(y[i:i+4000]), TT)) @ FVv
    return out
def ce(z): return np.exp(4j*np.pi*np.real(z))
def e_add(z): return np.exp(4j*np.pi*np.imag(z)/math.sqrt(3))
# ---------------- tables (as in t5_reflection_j2.py) ----------------
plist = primary_primes(B)
gp = {p: g_prime(p) for p in plist}
print("prime Gauss sums", len(plist), round(time.time()-t0, 1))
sq = squarefree_list(B, plist)
L2 = (-3, 0); wL2 = mul((0, 1), L2); w2L2 = mul((-1, -1), L2)
cdata = []
for c, f in sq:
    g = 1+0j; s1 = s2 = s3 = 1; slam = 1
    for p in f:
        g *= gp[p]; s1 *= cubic_prime(L2, p); s2 *= cubic_prime(wL2, p); s3 *= cubic_prime(w2L2, p)
        slam *= cubic_prime(LAMp, p)
    for i in range(len(f)):
        for j in range(i+1, len(f)): g *= pairsym(f[i], f[j])*pairsym(f[j], f[i])
    cdata.append((cx(c), N(c), g, np.conj(s1)*g, np.conj(s2)*g, np.conj(s3)*g, slam, (-2, 0) in f))
print("squarefree c:", len(cdata), round(time.time()-t0, 1))
Dmax = int(round(B**(1/3)))+2
dl = []
for x in range(-2*Dmax, 2*Dmax+1):
    for y in range(-2*Dmax, 2*Dmax+1):
        if primary((x, y)) and N((x, y))**3 <= B: dl.append((cx((x, y)), N((x, y)), (x % 2 == 0 and y % 2 == 0)))
E9 = cmath.exp(2j*math.pi/9)
def table_d0(Nmax):
    L = []; Cf = []
    for cc, nc, g1_, gl2, gwl2, gw2l2, _, _ in cdata:
        if nc > 27*Nmax: continue
        for dd, nd, _ in dl:
            base = nc*nd**3
            if base > 27*Nmax: continue
            absdc = abs(dd)/abs(cc)
            for k in range(-3, 40):
                if base*3.0**k > Nmax*(1+1e-12): break
                r = k % 3
                if r == 1: continue
                if r == 2:
                    nu = (k+4)//3
                    vals = [(1, np.conj(gl2)), (-1, np.conj(gl2)), (W, E9**-1*np.conj(gwl2)), (-W, E9**-1*np.conj(gwl2)), (W*W, E9*np.conj(gw2l2)), (-W*W, E9*np.conj(gw2l2))]
                    fac = absdc*3**(nu/2+2)
                else:
                    nu = (k+3)//3
                    vals = [(1, np.conj(g1_)), (-1, np.conj(g1_))]
                    fac = absdc*3**(nu/2+5/2)
                lk = LAMc**k*cc*dd**3
                for u, tv in vals:
                    L.append(u*lk); Cf.append(np.conj(tv*fac))
    return np.array(L), np.array(Cf)
def table_dpm(sign, Nmax):
    L = []; Cf = []
    for cc, nc, g1_, gl2, gwl2, gw2l2, _, _ in cdata:
        if nc > 81*Nmax: continue
        for dd, nd, _ in dl:
            base = nc*nd**3
            if base > 81*Nmax: continue
            absdc = abs(dd)/abs(cc); m0 = LAMc**-4*cc*dd**3
            if sign == '+':
                for up, t2 in [(1, 9*W*W*np.conj(gl2)), (W, 9*E9**-1*np.conj(gwl2)), (W*W, 9*W*E9*np.conj(gw2l2))]:
                    l = up*W*W*m0; L.append(l); Cf.append(W*W*np.conj(t2*absdc)*ce(l))
            else:
                for up, t1 in [(1, 9*W*np.conj(gl2)), (W, 9*E9**-1*W*W*np.conj(gwl2)), (W*W, 9*E9*np.conj(gw2l2))]:
                    l = -up*W*m0; L.append(l); Cf.append(W*np.conj(t1*absdc)*ce(l))
    return np.array(L), np.array(Cf)
# ---------------- cusp data (as in t5_reflection_j2.py) ----------------
def mod3eq(x, y): d = sub(x, y); return d[0] % 3 == 0 and d[1] % 3 == 0
def vlam(x):
    k = 0
    while divides(LAMp, x): x = div(x, LAMp); k += 1
    return k
def cubic_general(x, a):
    r = 1; m = a
    for p in plist:
        if N(m) == 1: break
        while divides(p, m):
            r *= cubic_prime(x, p); m = div(m, p)
    assert N(m) == 1
    return r
def matmul(A, Bm):
    (a, b), (c, d) = A; (e, f), (g, h) = Bm
    return ((add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h))), (add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h))))
def inv(A):
    (a, b), (c, d) = A; return ((d, (-b[0], -b[1])), ((-c[0], -c[1]), a))
def cusp_data(a, C):
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
    R = 14
    for x in range(-R, R+1):
        for y in range(-R, R+1):
            dp = (x, y); t = sub(mul(a, dp), one)
            if not divides(C, t): continue
            bg = div(t, C); g = ((a, bg), (C, dp)); (A1, B1), (C1, D1) = matmul(g, inv(H))
            if mod3eq(A1, one) and mod3eq(D1, one) and mod3eq(B1, zero) and mod3eq(C1, zero) and N(A1) <= 4000:
                kap = 1 if C1 == (0, 0) else cubic_general(C1, A1)
                return sigma, dp, kap
    raise ValueError((a, C))
def gcd_O(x, y):
    while y != (0, 0):
        x, y = y, reduce_mod(x, y)
    return x
UNITS6 = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]
def normalize(a, c):
    if divides(LAMp, c):
        for u in UNITS6:
            if mod3eq(mul(u, a), (1, 0)): return mul(u, a), mul(u, c)
    else:
        for u in UNITS6:
            if mod3eq(mul(u, c), (1, 0)): return mul(u, a), mul(u, c)
    raise ValueError
# ---------------- checks of the supplement and of gamma_2 ----------------
bad = 0
for p in plist[:300]:
    if p == (-2, 0): continue
    bad += abs(cubic_prime(LAMp, p)-W**((-(p[1]//3)) % 3)) > 1e-9
print("supplement (lambda/p)_3 = omega^{-b/3}: failures", bad)
# ---------------- left side and right side for each twist ----------------
def lam3(x):            # (lambda/x)_3 for primary x (supplement): omega^{-b/3}
    return W**((-(x[1]//3)) % 3)
def om3(x):             # (omega/x)_3 = omega^{(N(x)-1)/3} for x coprime to 3
    return W**(((N(x)-1)//3) % 3)
XI = {'1': lambda x: 1, 'A': lambda x: np.conj(lam3(x)), 'B': lambda x: np.conj(lam3(x))*om3(x)}
# for the direct sum we need xi_0(n) for squarefree n from cdata: store the factorizations
sqf = {cx(c): f for c, f in sq}
def xi_n(tw, cc):
    r = 1
    for p in sqf[cc]:
        r *= XI[tw](p)
    return r
L18 = (18, 0); RES = residues(L18)
tabcache = {}
results = {}
for tw in ('1', 'A', 'B'):
    Tdir = 0j; nterms = 0
    for cc, nc, g1_, _, _, _, slam, even in cdata:
        if even or nc*1.0 > X*math.exp(KCUT*s_w): continue
        gam2 = slam*g1_/abs(cc); an = cc/abs(cc); xin = xi_n(tw, cc)
        for dd, nd, deven in dl:
            if deven: continue
            y = nc*nd**3/X
            if y > math.exp(KCUT*s_w): continue
            ab = dd/abs(dd)
            Tdir += np.conj(an)**2*gam2*xin*np.conj(ab)**6/(math.sqrt(nc)*nd)*V(y); nterms += 1   # Psi(b)^3 = 1
    def phi(x):
        if not (x[0] % 3 == 1 and x[1] % 3 == 0): return 0
        if x[0] % 2 == 0 and x[1] % 2 == 0: return 0
        return lam3(x)*XI[tw](x)
    phivals = [(cx(x), phi(x)) for x in RES]
    tot = {'main': 0j, 'a2': 0j, 'zeta': 0j, 'k1': 0j}; ncusp = 0; types = {}
    for h0 in RES:
        ph = sum(v*e_add(-cx(h0)*xc/18) for xc, v in phivals if v != 0)/324
        if abs(ph) < 1e-13: continue
        num = mul(mul(LAMp, LAMp), h0); den = L18
        if num == (0, 0): a, c = (0, 0), (1, 0)
        else:
            g = gcd_O(num, den); a, c = div(num, g), div(den, g)
        a, c = normalize(a, c)
        sigma, dp, kap = cusp_data(a, c)
        NC = N(c); key = (sigma, NC)
        if key not in tabcache:
            NR = 27*YCUT*NC**2/((2*math.pi)**4*X)
            assert (27 if sigma == '0' else 81)*NR <= B, key
            Lr, Cr = table_d0(NR) if sigma == '0' else table_dpm(sigma, NR)
            nl = np.abs(Lr)**2; arg = (2*np.pi)**4*nl*X/(27*NC**2)
            common = Cr*nl**-0.5*(Lr/np.abs(Lr))**2
            tabcache[key] = (Lr, common*Vnat(arg, FV2), common*Vnat(arg, FV1))
        Lr, base2, base1 = tabcache[key]
        Cc = cx(c); aC = np.conj(Cc/abs(Cc)); phase = ce(-cx(dp)*Lr/Cc)
        S2 = np.sum(base2*phase); S1 = np.sum(base1*phase)
        tot['main'] += -(1/81)*ph*np.conj(kap)*aC**4*S2
        tot['a2'] += -(1/81)*ph*np.conj(kap)*aC**2*S2
        tot['zeta'] += -(1j/81)*ph*np.conj(kap)*aC**4*S2
        tot['k1'] += -(1/81)*ph*np.conj(kap)*aC**4*S1
        ncusp += 1; types[(sigma, NC)] = types.get((sigma, NC), 0)+1
    err = abs(tot['main']-Tdir)/abs(Tdir)
    ctrl = {k: abs(tot[k]-Tdir)/abs(Tdir) for k in ('a2', 'zeta', 'k1')}
    results[tw] = (err, ctrl)
    print(f"twist {tw}: direct terms {nterms}; cusps {ncusp}, types (sigma,N(c)) {dict(sorted(types.items()))}")
    print(f"   T direct = {Tdir:.12f}   T via Prop = {tot['main']:.12f}   rel.err {err:.2e}")
    print("   controls (rel.err):", {k: f"{v:.2e}" for k, v in ctrl.items()}, round(time.time()-t0, 1))
worst = max(r[0] for r in results.values())
cmax = {k: max(r[1][k] for r in results.values()) for k in ('a2', 'zeta', 'k1')}
print("worst rel.err", f"{worst:.2e}", " control max over twists", {k: f"{v:.2e}" for k, v in cmax.items()})
print("PASS" if worst < 1e-6 and min(cmax.values()) > 1e-2 else "FAIL", "time", round(time.time()-t0, 1))

# Exit status (added 2026-10-12): 1 unless the supplement check had 0 failures, worst rel.err < 1e-6 and every control
# fails (rel.err > 1e-2), the criteria of the PASS/FAIL line above.  CM_TESTS_FORCE_FAIL=1 counts the first check as
# failed (self-test of the wiring).
checks = [("supplement failures == 0", bad == 0), ("worst rel.err < 1e-6", worst < 1e-6),
          ("every control rel.err > 1e-2", min(cmax.values()) > 1e-2)]
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    checks[0] = (checks[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, c in checks if not c]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(checks), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(checks)))
sys.exit(1 if bad else 0)
