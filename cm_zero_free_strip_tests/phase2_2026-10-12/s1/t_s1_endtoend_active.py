"""s1 end-to-end test of the CM completed reflection (Proposition prop:s1-reflection) with an ACTIVE local prime,
its grouped form (local factors B_p, scalar zeta^CM, cross-prime bookkeeping kappa_F, vartheta), and the marked
corollary (cor:s1-marked), for two non-compact log-Schwartz tests V.

Data.  Base twist 'A' of t_b3_prop_endtoend.py:  Psi_0(n) = conj((lambda/n)_3) on primary n prime to 6, so that
phi(x) = (lambda/x)_3 Psi_0(x) = 1 on primary x prime to 6 (period 6, used with L = 18).  One local prime
pi = PI (primary, norm 7) with exponent j in {0,...,5}:  Psi(n) = Psi_0(n) chi_pi(n)^j  (chi^0 = 1_{pi not | n}).

Left side (definition, OAI78 normalization, CM weight):
   T_V(X;Psi) = sum_{n sf, b; (nb,6)=1} gamma_2(n) conj(alpha(n))^2 Psi(n) conj(alpha(b))^6 Psi(b)^3 / (sqrt(N n) N b)
                * V(N(n) N(b)^3 / X).
Right side, UNGROUPED (Steps 1-6 of the proof, every translate h = (h0 mod 18, h_pi mod pi) separately):
   -(1/81) sum_h phihat(h0) C_{pi,j}(h_pi) conj(kappa(G_h)) conj(alpha(c_h))^4
           * sum_{l in lambda^-4 O} d_sigma(l) alpha(l)^2 N(l)^{-1/2} ce(-delta_h l / c_h) V#(N(l) X / N(c_h)^2),
   z_h = lambda^2 (h0/18 + h_pi/pi) = a/c reduced; cusp data (sigma, delta, kappa) by brute force (any valid choice).
Right side, GROUPED (statement of the proposition): for each h0 and active pi (c = c_F pi),
   zeta*vartheta(x) * B_pi(x),   zeta = -(1/81) conj(alpha(c))^4 phihat(h0) conj(kappa_F) chi_pi(sigma_pi)^{-2} varpi_{pi,j},
   where conj(kappa_F) vartheta(x) is extracted from ONE active translate, and the test checks that it is the same for all
   h_pi != 0 (this is OAI78 eq:reflection-multiplier + eq:reflection-additive-phase, the cross-prime bookkeeping);
   plus, for j = 0, the inactive branch with factor (1 - 1/N(pi)).
Marked: T_Q with Q = {pi} (insert 1_{pi | n b^3}, base twist A) against cor:s1-marked: inactive branch N(pi)^{-1} times
   the P = {} branch; active branch zeta_Q = -(1/81)conj(alpha(c))^4 phihat conj(kappa_F) chi(sigma)^{-2} (-varpi_{pi,0}),
   local factor N(pi)^{-1/2} chi_pi(x)^{-2}.
Tests V:  'gauss'  V(y) = exp(-(log y)^2/(2 s^2)), s = 1/2,  Vhat(t) = sqrt(2 pi) s exp(s^2 t^2/2);
          'bessel' V(y) = exp(-y - 1/y),                     Vhat(t) = 2 K_t(2)   (exponential, not Gaussian, decay).
Controls (must fail): conj(alpha(c))^2;  zeta prefactor -i/81;  one-derivative kernel R_1;  B_pi with chi^{-j}
   instead of chi^{-j-2} (multiplier chi_p(a)^2 forgotten);  varpi with eps_p replaced by 1.
Usage: python3 t_s1_endtoend_active.py [gauss|bessel] [pi_index]     (a few minutes)
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastgauss import *
from zw import canon_mod, residues, reduce_mod, sym6_prime
from scipy.special import loggamma
import numpy as np, mpmath as mp
t0 = time.time()
VKIND = sys.argv[1] if len(sys.argv) > 1 else 'gauss'
PIDX = int(sys.argv[2]) if len(sys.argv) > 2 else 0
LAMp = (1, 2); LAMc = cx(LAMp)
X = 3000.0; B = 150000; YCUT = 4e5
s_w = 0.5
if VKIND == 'gauss':
    YMAX = math.exp(7.5*s_w)
    def V(y): return np.exp(-np.log(y)**2/(2*s_w*s_w))
    TAU = np.arange(-14/s_w, 14/s_w+1e-9, 0.05)
    VHAT = np.sqrt(2*np.pi)*s_w*np.exp(-s_w*s_w*TAU*TAU/2)          # Vhat(-i tau)
else:
    YMAX = 45.0
    def V(y): return np.exp(-y-1/y)
    TAU = np.arange(-40, 40+1e-9, 0.05)
    VHAT = np.array([float(2*mp.re(mp.besselk(-1j*t, 2))) for t in TAU])   # Vhat(-i tau) = 2 K_{-i tau}(2) (real)
TT = 1j*TAU
def R2(t): return np.exp(loggamma(4/3+t)+loggamma(5/3+t)-loggamma(4/3-t)-loggamma(5/3-t))
def R1(t): return np.exp(loggamma(7/6+t)+loggamma(5/6+t)-loggamma(7/6-t)-loggamma(5/6-t))
FV2 = VHAT*R2(TT)*0.05/(2*np.pi); FV1 = VHAT*R1(TT)*0.05/(2*np.pi)
def Vnat(y, FVv):
    y = np.asarray(y, float); out = np.empty(len(y), complex)
    for i in range(0, len(y), 4000):
        out[i:i+4000] = np.exp(-np.outer(np.log(y[i:i+4000]), TT)) @ FVv
    return out
def ce(z): return np.exp(4j*np.pi*np.real(z))
def e_add(z): return np.exp(4j*np.pi*np.imag(z)/math.sqrt(3))
# ---------------- tables ----------------
plist = primary_primes(B)
PI = [p for p in plist if N(p) == 7][PIDX]; NPI = 7; PIc = cx(PI)
print("V =", VKIND, " pi =", PI)
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
    cdata.append((c, N(c), g, np.conj(s1)*g, np.conj(s2)*g, np.conj(s3)*g, slam, (-2, 0) in f, PI in f))
print("squarefree c:", len(cdata), round(time.time()-t0, 1))
Dmax = int(round(B**(1/3)))+2
dl = []
for x in range(-2*Dmax, 2*Dmax+1):
    for y in range(-2*Dmax, 2*Dmax+1):
        if primary((x, y)) and N((x, y))**3 <= B: dl.append(((x, y), N((x, y)), (x % 2 == 0 and y % 2 == 0)))
E9 = cmath.exp(2j*math.pi/9)
UE = {1: (1, 0), -1: (-1, 0), W: (0, 1), -W: (0, -1), W*W: (-1, -1), -W*W: (1, 1)}   # complex unit -> exact
def ukey(u):
    for k, v in UE.items():
        if abs(k-u) < 1e-9: return v
    raise ValueError(u)
def table_d0(Nmax):
    L = []; Cf = []; XE = []
    for cc_, nc, g1_, gl2, gwl2, gw2l2, _, _, _ in cdata:
        if nc > 27*Nmax: continue
        cc = cx(cc_)
        for dd_, nd, _ in dl:
            base = nc*nd**3
            if base > 27*Nmax: continue
            dd = cx(dd_); absdc = abs(dd)/abs(cc); cd3 = mul(cc_, pw(dd_, 3))
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
                lk = LAMc**k*cc*dd**3; xk = mul(pw(LAMp, k+4), cd3)        # lambda^4 * lambda^k c d^3 (exact)
                for u, tv in vals:
                    L.append(u*lk); Cf.append(np.conj(tv*fac)); XE.append(mul(ukey(u), xk))
    return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64)
def table_dpm(sign, Nmax):
    L = []; Cf = []; XE = []
    for cc_, nc, g1_, gl2, gwl2, gw2l2, _, _, _ in cdata:
        if nc > 81*Nmax: continue
        cc = cx(cc_)
        for dd_, nd, _ in dl:
            base = nc*nd**3
            if base > 81*Nmax: continue
            dd = cx(dd_); absdc = abs(dd)/abs(cc); m0 = LAMc**-4*cc*dd**3; cd3 = mul(cc_, pw(dd_, 3))
            if sign == '+':
                for up, t2 in [(1, 9*W*W*np.conj(gl2)), (W, 9*E9**-1*np.conj(gwl2)), (W*W, 9*W*E9*np.conj(gw2l2))]:
                    l = up*W*W*m0; L.append(l); Cf.append(W*W*np.conj(t2*absdc)*ce(l)); XE.append(mul(ukey(up*W*W), cd3))
            else:
                for up, t1 in [(1, 9*W*np.conj(gl2)), (W, 9*E9**-1*W*W*np.conj(gwl2)), (W*W, 9*E9*np.conj(gw2l2))]:
                    l = -up*W*m0; L.append(l); Cf.append(W*np.conj(t1*absdc)*ce(l)); XE.append(mul(ukey(-up*W), cd3))
    return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64)
# ---------------- residue arithmetic modulo pi (Z[omega]/pi = F_7) ----------------
w0 = (-PI[0]*pow(PI[1], -1, NPI)) % NPI          # image of omega in F_7
def red(x): return (x[0]+x[1]*w0) % NPI          # exact element -> F_7
SIX = {}                                          # F_7^* -> sextic symbol (x/pi)_6 as complex
for r in range(1, NPI):
    SIX[r] = cmath.exp(1j*math.pi*sym6_prime((r, 0), PI)/3)
CHI = np.zeros(NPI, complex)
for r in range(1, NPI): CHI[r] = SIX[r]
def chipow(r, j):                                 # chi_pi(r)^j with the zero convention (also for j = 0)
    r = np.asarray(r) % NPI
    out = np.where(r == 0, 0, CHI[r]**j)
    return out
EPI = np.array([e_add(t/PIc) for t in range(NPI)])   # e(t/pi), t in F_7 (e(y/pi) depends on y mod pi)
def inv7(r): return pow(int(r) % NPI, -1, NPI)
# ---------------- cusp data ----------------
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
    best = None
    for Rb in (14, 30):
        for x in range(-Rb, Rb+1):
            for y in range(-Rb, Rb+1):
                dp = (x, y); t = sub(mul(a, dp), one)
                if not divides(C, t): continue
                bg = div(t, C); g = ((a, bg), (C, dp)); (A1, B1), (C1, D1) = matmul(g, inv(H))
                if mod3eq(A1, one) and mod3eq(D1, one) and mod3eq(B1, zero) and mod3eq(C1, zero):
                    if best is None or N(A1) < best[0]: best = (N(A1), dp, A1, C1)
        if best is not None and best[0] <= B: break
    _, dp, A1, C1 = best
    kap = 1 if C1 == (0, 0) else cubic_general(C1, A1)
    return sigma, dp, kap
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
# ---------------- twist A ----------------
def lam3(x): return W**((-(x[1]//3)) % 3)               # (lambda/x)_3 for primary x
def psi0(x): return np.conj(lam3(x))                     # Psi_0 = conj((lambda/.)_3)
sqf = {c: f for c, f in sq}
L18 = (18, 0); RES = residues(L18); RESP = residues(PI)
def phi(x):
    if not (x[0] % 3 == 1 and x[1] % 3 == 0): return 0
    if x[0] % 2 == 0 and x[1] % 2 == 0: return 0
    return 1                                             # lam3(x)*psi0(x) = 1
phivals = [(cx(x), phi(x)) for x in RES]
PHIHAT = {h0: sum(v*e_add(-cx(h0)*xc/18) for xc, v in phivals if v != 0)/324 for h0 in RES}
def Cpj(h, j):                                           # C_{pi,j}(h) = N(pi)^{-1} sum_y chi(y)^j e(-h y/pi)
    return sum(chipow(red(y), j)*e_add(-cx(h)*cx(y)/PIc) for y in RESP)/NPI
def gauss_pm(m, sgn):                                    # N(pi)^{-1/2} sum_y chi(y)^m e(sgn y/pi)
    return sum(chipow(red(y), m % 6)*e_add(sgn*cx(y)/PIc) for y in RESP)/math.sqrt(NPI)
# ---------------- direct sums ----------------
def direct(j, mark=False):
    T = 0j
    for cc_, nc, g1_, _, _, _, slam, even, haspi in cdata:
        if even or nc > X*YMAX: continue
        cc = cx(cc_); gam2 = slam*g1_/abs(cc); an = cc/abs(cc)
        xin = 1
        for p in sqf[cc_]: xin *= psi0(p)
        if mark: chin = 1                                    # Psi = Psi_0, mark 1_{pi | n b^3}
        else: chin = chipow(red(cc_), j)                     # chi_pi(n)^j, zero if pi | n
        for dd_, nd, deven in dl:
            if deven: continue
            y = nc*nd**3/X
            if y > YMAX: continue
            ab = cx(dd_)/math.sqrt(nd)
            if mark:
                if not (haspi or red(dd_) == 0): continue
                chib = 1
            else:
                chib = chipow(red(dd_), (3*j) % 6)          # Psi(b)^3 = chi_pi(b)^{3j}, zero at pi | b
            T += np.conj(an)**2*gam2*xin*chin*np.conj(ab)**6*chib/(math.sqrt(nc)*nd)*V(y)
    return T
# Psi(b)^3 = Psi_0(b)^3 chi_pi(b)^{3j} = chi_pi(b)^{3j} (Psi_0 is cubic); zero at pi | b for every j (incl. 0).
# ---------------- translates ----------------
tabcache = {}
def table(sigma, NC):
    key = (sigma, NC)
    if key not in tabcache:
        NR = 27*YCUT*NC**2/((2*math.pi)**4*X)
        Lr, Cr, XE = table_d0(NR) if sigma == '0' else table_dpm(sigma, NR)
        nl = np.abs(Lr)**2; arg = (2*np.pi)**4*nl*X/(27*NC**2)
        common = Cr*nl**-0.5*(Lr/np.abs(Lr))**2
        xr = (XE[:, 0]+XE[:, 1]*w0) % NPI
        tabcache[key] = (Lr, common*Vnat(arg, FV2), common*Vnat(arg, FV1), xr, Cr)
    return tabcache[key]
def translate(h0, hp):
    num = mul(mul(LAMp, LAMp), add(mul(h0, PI), mul((18, 0), hp))); den = mul(L18, PI)
    if num == (0, 0): a, c = (0, 0), (1, 0)
    else:
        g = gcd_O(num, den); a, c = div(num, g), div(den, g)
    a, c = normalize(a, c)
    a = add(reduce_mod(sub(a, (1, 0)), mul((3, 0), c)), (1, 0)) if divides(LAMp, c) else reduce_mod(a, mul((3, 0), c))
    sigma, dp, kap = cusp_data(a, c)
    return a, c, sigma, dp, kap
TR = {}
H0S = [h0 for h0 in RES if abs(PHIHAT[h0]) > 1e-13]
for h0 in H0S:
    for hp in RESP:
        TR[(h0, hp)] = translate(h0, hp)
print("translates:", len(TR), " h0 classes with phihat != 0:", len(H0S), round(time.time()-t0, 1))
types = {}
for (h0, hp), (a, c, sigma, dp, kap) in TR.items(): types[(sigma, N(c))] = types.get((sigma, N(c)), 0)+1
print("cusp types (sigma, N(c)):", dict(sorted(types.items())))
# per-translate dual sums (ungrouped) and per-l data for the grouped check
def dual(h0, hp, which=2):
    a, c, sigma, dp, kap = TR[(h0, hp)]
    Lr, b2, b1, xr, Cr = table(sigma, N(c))
    Cc = cx(c); phase = ce(-cx(dp)*Lr/Cc)
    return np.conj(kap)*phase, (b2 if which == 2 else b1), Cc, xr, Cr
results = {}
for j in range(6):
    Tdir = direct(j)
    U = {'main': 0j, 'a2': 0j, 'zeta': 0j, 'k1': 0j}
    G = {'main': 0j, 'Bj': 0j, 'eps1': 0j}
    worst_lam = 0.0
    for h0 in H0S:
        ph = PHIHAT[h0]
        # ungrouped
        for hp in RESP:
            Cf_ = Cpj(hp, j)
            if abs(Cf_) < 1e-14: continue
            kp, b2, Cc, xr, Cr = dual(h0, hp, 2); _, b1, _, _, _ = dual(h0, hp, 1)
            aC = np.conj(Cc/abs(Cc))
            S2 = np.sum(b2*kp); S1 = np.sum(b1*kp)
            U['main'] += -(1/81)*ph*Cf_*aC**4*S2
            U['a2'] += -(1/81)*ph*Cf_*aC**2*S2
            U['zeta'] += -(1j/81)*ph*Cf_*aC**4*S2
            U['k1'] += -(1/81)*ph*Cf_*aC**4*S1
        # grouped: inactive branch (j = 0 only)
        if j == 0:
            kp, b2, Cc, xr, Cr = dual(h0, (0, 0), 2); aC = np.conj(Cc/abs(Cc))
            val = -(1/81)*ph*(1-1/NPI)*aC**4*np.sum(b2*kp)
            for k in G: G[k] += val
        # grouped: active branch, kappa_F*vartheta from each h_pi != 0, must agree
        lams = []
        for hp in RESP:
            if hp == (0, 0) or red(hp) == 0: continue
            a, c, sigma, dp, kap = TR[(h0, hp)]
            kp, b2, Cc, xr, Cr = dual(h0, hp, 2)
            cp = div(c, PI); sig_p = red(mul(mul(LAMp, LAMp), cp)); l3 = red(mul(pw(LAMp, 3), cp))
            eps_p = (-inv7(l3*sig_p)) % NPI; hinv = inv7(red(hp))
            Lam = kp*CHI[(sig_p*red(hp)) % NPI]**2*EPI[(-eps_p*hinv*xr) % NPI]
            lams.append((Lam, c, sigma, sig_p, eps_p, b2, xr, Cr))
        Lam0 = lams[0][0]; supp = np.abs(lams[0][7]) > 1e-12
        for Lam, *_ in lams[1:]:
            worst_lam = max(worst_lam, float(np.max(np.abs(Lam-Lam0)[supp])))
        Lam, c, sigma, sig_p, eps_p, b2, xr, Cr = lams[0]
        aC = np.conj(cx(c)/abs(cx(c)))
        # varpi_{p,j} (OAI78 L1832-1840), with g^{pm}_{p,m}; chi_p(sigma_p)^{-2}
        if j == 4: varpi = gauss_pm(4, -1)
        elif j == 0: varpi = -gauss_pm(2, +1)*CHI[eps_p]**-2
        else: varpi = gauss_pm(j, -1)*gauss_pm(j+2, +1)*CHI[eps_p]**(-j-2)
        e3 = (3*eps_p) % NPI                               # control: eps_p replaced by 3 eps_p (3 generates F_7^*)
        if j == 4: varpi_e1 = varpi
        elif j == 0: varpi_e1 = -gauss_pm(2, +1)*CHI[e3]**-2
        else: varpi_e1 = gauss_pm(j, -1)*gauss_pm(j+2, +1)*CHI[e3]**(-j-2)
        pref = -(1/81)*ph*aC**4*CHI[sig_p]**-2
        if j == 4: Bx = NPI**-0.5*(-1+NPI*(xr == 0))
        elif j == 0: Bx = NPI**-0.5*chipow(xr, -2 % 6)
        else: Bx = chipow(xr, (-j-2) % 6)
        Bwrong = chipow(xr, (-j) % 6) if j not in (0, 4) else Bx*chipow(xr, 2)
        G['main'] += pref*varpi*np.sum(b2*Lam*Bx)
        G['Bj'] += pref*varpi*np.sum(b2*Lam*Bwrong)
        G['eps1'] += pref*varpi_e1*np.sum(b2*Lam*Bx)
    eu = abs(U['main']-Tdir)/abs(Tdir); eg = abs(G['main']-Tdir)/abs(Tdir)
    ctrl = {k: abs(U[k]-Tdir)/abs(Tdir) for k in ('a2', 'zeta', 'k1')}
    ctrl.update({k: abs(G[k]-Tdir)/abs(Tdir) for k in ('Bj', 'eps1') if not (k == 'eps1' and j == 4)})
    results[j] = (eu, eg, worst_lam, ctrl)
    print(f"j={j}: T direct {Tdir:.12f}\n     ungrouped {U['main']:.12f} rel.err {eu:.2e}\n     grouped   {G['main']:.12f} rel.err {eg:.2e}"
          f"   max|kappaF*vartheta(h_pi) - (h_pi=1)| = {worst_lam:.1e}")
    print("     controls:", {k: f"{v:.2e}" for k, v in ctrl.items()}, round(time.time()-t0, 1))
# ---------------- marked corollary, Q = {pi}, base twist A ----------------
TQ = direct(0, mark=True)
M = {'main': 0j, 'noact': 0j}
for h0 in H0S:
    ph = PHIHAT[h0]
    kp, b2, Cc, xr, Cr = dual(h0, (0, 0), 2); aC = np.conj(Cc/abs(Cc))
    val = -(1/81)*ph*aC**4*np.sum(b2*kp)/NPI                # inactive: N(pi)^{-1} x (P = {} branch)
    M['main'] += val; M['noact'] += val
    hp1 = [hp for hp in RESP if red(hp) != 0][0]
    a, c, sigma, dp, kap = TR[(h0, hp1)]
    kp, b2, Cc, xr, Cr = dual(h0, hp1, 2)
    cp = div(c, PI); sig_p = red(mul(mul(LAMp, LAMp), cp)); l3 = red(mul(pw(LAMp, 3), cp))
    eps_p = (-inv7(l3*sig_p)) % NPI; hinv = inv7(red(hp1))
    Lam = kp*CHI[(sig_p*red(hp1)) % NPI]**2*EPI[(-eps_p*hinv*xr) % NPI]
    aC = np.conj(cx(c)/abs(cx(c)))
    zetaQ = -(1/81)*ph*aC**4*CHI[sig_p]**-2*(gauss_pm(2, +1)*CHI[eps_p]**-2)    # -varpi_{q,0} = g^+_{q,2} chi(eps)^{-2}
    M['main'] += zetaQ*np.sum(b2*Lam*NPI**-0.5*chipow(xr, 4))
em = abs(M['main']-TQ)/abs(TQ); emc = abs(M['noact']-TQ)/abs(TQ)
print(f"marked Q={{pi}}: T_Q direct {TQ:.12f}  via cor {M['main']:.12f}  rel.err {em:.2e}  (control: active branch dropped {emc:.2e})")
worst = max(max(r[0], r[1]) for r in results.values()); wl = max(r[2] for r in results.values())
cmin = min(min(r[3].values()) for r in results.values())
ok = worst < 1e-7 and wl < 1e-9 and em < 1e-7 and cmin > 1e-3 and emc > 1e-3
print(f"worst rel.err {worst:.2e}; worst kappa_F*vartheta spread {wl:.1e}; marked {em:.2e}; smallest control error {cmin:.2e}")
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':   # self-test of the exit status (added 2026-10-12)
    ok = False; print("CM_TESTS_FORCE_FAIL=1: the check below is counted as failed")
print("PASS" if ok else "FAIL", "time", round(time.time()-t0, 1))

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
