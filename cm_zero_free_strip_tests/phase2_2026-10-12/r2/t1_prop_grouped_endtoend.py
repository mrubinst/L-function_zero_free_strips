"""r2 test T1: Proposition prop:b3-reflection end to end, in its GROUPED form, with one prime in P.

Twist: Psi(n) = xi_0(n) chi_pi(n)^j with xi_0 = conj((lambda/.)_3) (a character of the primary classes mod 72),
N(pi) = 7, and every j = 0,...,5.  Then phi(x) = (lambda/x)_3 xi_0(x) = 1 on primary x prime to 6.
Left side: T(X;Psi) computed from its definition (8.1), with the log-normal V(y) = exp(-(log y)^2/(2 s^2)).
Right side: the sum over the groups (h0, A), h0 in O/72O and {pi : j_pi != 0} c A c {pi}, of
   C * sum_{0 != l in lambda^-4 O} d_sigma(l) alpha(l)^2 N(l)^{-1/2} psi(lambda^4 l) prod_{A} B_{pi,j}(lambda^4 l)
       * V#(N(l) X/N(c)^2),
   C = -(conj(alpha(c))^4/81) phihat(h0) conj(kappa0) prod_{inactive}(1 - 1/N(pi)) prod_{A} chi_pi(sigma_pi)^{-2} varpi_{pi,j},
exactly as in (A.1) and (A.16), with the cusp data (c, sigma, kappa0, beta, sigma_pi, eps_pi) computed from
[OAI1112, App. A.2] with L = 72 and M = lambda^12 72^4 (helpers r2lib.py, as in t2_grouping_cusp_data.py), the
coefficients d_sigma from [QF, Theorem 23] and the CM kernel of (8.4).  This tests the proposition as stated (with
local factors and the constant (A.16)), and in particular conj(alpha(c))^4 at cusps c = c0*pi that are not real.
Falsifiers (must fail): alpha(c)^4 (unconjugated); conj(alpha(c))^2; zeta prefactor -i/81; kappa0 unconjugated;
B_{pi,j} with exponent -j instead of -j-2; one-derivative kernel.
Theta tables and V# follow tests_cm/proofs_2026-10-12/b3/t_b3_prop_endtoend.py (helpers fastgauss.py, zw.py imported
read-only from there).  Usage: python3 -B t1_prop_grouped_endtoend.py     (a few minutes)
"""
import os, sys, time, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
B3 = os.path.join(HERE, '..', '..', 'proofs_2026-10-12', 'b3')
sys.path.insert(0, B3)
import fastgauss as fg
import zw
sys.path.insert(0, HERE)
import r2lib as R
import numpy as np
from scipy.special import loggamma
t0 = time.time()
W = fg.W
LAMp = (1, 2); LAMc = fg.cx(LAMp)
s_w = 0.5; KCUT = 7.5; YCUT = 4e5
X = 3200.0; B = 150000
PI = R.PI7; NPI = R.N(PI)
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
# ---------------- theta tables (as in t_b3_prop_endtoend.py) ----------------
plist = fg.primary_primes(B)
gp = {p: fg.g_prime(p) for p in plist}
print("prime Gauss sums", len(plist), round(time.time()-t0, 1), flush=True)
sq = fg.squarefree_list(B, plist)
L2 = (-3, 0); wL2 = fg.mul((0, 1), L2); w2L2 = fg.mul((-1, -1), L2)
cdata = []
for c, f in sq:
    g = 1+0j; s1 = s2 = s3 = 1; slam = 1
    for p in f:
        g *= gp[p]; s1 *= fg.cubic_prime(L2, p); s2 *= fg.cubic_prime(wL2, p); s3 *= fg.cubic_prime(w2L2, p)
        slam *= fg.cubic_prime(LAMp, p)
    for i in range(len(f)):
        for j in range(i+1, len(f)): g *= fg.pairsym(f[i], f[j])*fg.pairsym(f[j], f[i])
    cdata.append((c, fg.cx(c), fg.N(c), g, np.conj(s1)*g, np.conj(s2)*g, np.conj(s3)*g, slam, (-2, 0) in f))
print("squarefree c:", len(cdata), round(time.time()-t0, 1), flush=True)
Dmax = int(round(B**(1/3)))+2
dl = []
for x in range(-2*Dmax, 2*Dmax+1):
    for y in range(-2*Dmax, 2*Dmax+1):
        if fg.primary((x, y)) and fg.N((x, y))**3 <= B: dl.append(((x, y), fg.cx((x, y)), fg.N((x, y)), (x % 2 == 0 and y % 2 == 0)))
E9 = cmath.exp(2j*math.pi/9)
def table_d0(Nmax):
    L = []; Cf = []
    for _, cc, nc, g1_, gl2, gwl2, gw2l2, _, _ in cdata:
        if nc > 27*Nmax: continue
        for _, dd, nd, _ in dl:
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
    for _, cc, nc, g1_, gl2, gwl2, gw2l2, _, _ in cdata:
        if nc > 81*Nmax: continue
        for _, dd, nd, _ in dl:
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
# ---------------- chi_pi on residues mod pi (pi split, N = 7) ----------------
p7 = NPI; w0 = [w for w in range(p7) if (w*w + w + 1) % p7 == 0 and R.divides(PI, (-w, 1))][0]
chi_tab = {}
for r_ in range(p7):
    chi_tab[r_] = complex(R.sext_sym_prime((r_, 0), PI)) if r_ else 0j
def res_mod_pi(a, b): return (a + b*w0) % p7          # image of a + b*omega in F_7 = O/pi
def chi_pi_int(x): return chi_tab[res_mod_pi(x[0], x[1])]
def coords(z):                                         # complex -> (a, b) with z = a + b*omega (rounded)
    b = np.rint(np.imag(z)/(math.sqrt(3)/2)).astype(np.int64); a = np.rint(np.real(z) + b/2).astype(np.int64)
    return a, b
chi_arr = np.array([chi_tab[r_] for r_ in range(p7)])
# ---------------- left side ----------------
JS = range(6)
Tdir = {j: 0j for j in JS}; Tdir_b3 = {j: 0j for j in JS}
for c, cc, nc, g1_, _, _, _, slam, even in cdata:
    if even or nc*1.0 > X*math.exp(KCUT*s_w): continue
    gt = g1_/abs(cc); an = cc/abs(cc); chn = chi_pi_int(c)       # gamma_2(n) xi_0(n) = gtilde(n)
    for b_, dd, nd, deven in dl:
        if deven: continue
        y = nc*nd**3/X
        if y > math.exp(KCUT*s_w): continue
        ab = dd/abs(dd); chb = chi_pi_int(b_)
        w = np.conj(an)**2*gt*V(y)/(math.sqrt(nc)*nd)
        for j in JS:
            pn = chn**j if j else (1.0 if chn != 0 else 0.0)
            pb = chb**(3*j) if j else (1.0 if chb != 0 else 0.0)     # Psi(b)^3; xi_0(b)^3 = 1
            Tdir[j] += w*pn*np.conj(ab)**6*pb
            Tdir_b3[j] += w*pn*np.conj(ab)**3*pb                    # falsifier on the left: conj(alpha(b))^3
print("left side done", round(time.time()-t0, 1), flush=True)
# ---------------- groups and cusp data ----------------
L72 = (72, 0)
RES72 = zw.residues(L72)
def phi(x):
    if not (x[0] % 3 == 1 and x[1] % 3 == 0): return 0
    if x[0] % 2 == 0 and x[1] % 2 == 0: return 0
    return 1
phivals = [(fg.cx(x), phi(x)) for x in RES72 if phi(x)]
def phihat(h0):
    hc = fg.cx(h0)
    return sum(e_add(-hc*xc/72) for xc, _ in phivals)/5184
def group_data(h0, A):
    if A:
        hres = R.nonzero_res(PI)[0]
        hs = [R.mul(R.M2, R.red(R.mul(hres, R.inv_mod(R.M2, PI)), PI))]
        d = R.translate_data(h0, A, hs)
        kap0 = complex(d['kap'])/complex(R.cubic_sym_prime(d['a'], PI))
    else:
        if h0 == (0, 0):
            d = dict(a=(0, 0), c=(1, 0), c0=(1, 0), r=(1, 0), dp=(0, 0), sigma='0', vl=0); kap0 = 1.0
        else:
            d = R.translate_data(h0, [], []); kap0 = complex(d['kap'])
    c = d['c']; c0 = d['c0']; r = d['r']
    D0 = R.mul(R.pw(R.LAM, 3), c0)
    beta = R.red(R.neg(R.mul(d['dp'], R.inv_mod(r, D0))), D0) if R.N(D0) > 1 else (0, 0)
    out = dict(c=c, c0=c0, sigma=d['sigma'], kap0=kap0, beta=beta)
    if A:
        sig = R.div(R.mul(R.mul(R.LAM, R.LAM), c), PI)
        eps = R.red(R.neg(R.inv_mod(R.mul(R.div(R.mul(R.pw(R.LAM, 3), c), PI), sig), PI)), PI)
        out.update(sigpi=sig, eps=eps)
    return out
tabcache = {}
def dual_base(sigma, NC):
    key = (sigma, NC)
    if key not in tabcache:
        NR = 27*YCUT*NC**2/((2*math.pi)**4*X)
        assert (27 if sigma == '0' else 81)*NR <= B, key
        Lr, Cr = table_d0(NR) if sigma == '0' else table_dpm(sigma, NR)
        nl = np.abs(Lr)**2; arg = (2*np.pi)**4*nl*X/(27*NC**2)
        common = Cr*nl**-0.5*(Lr/np.abs(Lr))**2
        a9, b9 = coords(9*Lr)                                   # x = lambda^4 l = 9 l in O
        assert np.max(np.abs(a9 + b9*W - 9*Lr)) < 1e-6
        resx = (a9 + b9*w0) % p7
        tabcache[key] = (Lr, common*Vnat(arg, FV2), common*Vnat(arg, FV1), resx)
        print(f"   table {key}: {len(Lr)} coefficients, t={time.time()-t0:.0f}s", flush=True)
    return tabcache[key]
def Bvec(j, resx, wrong=False):
    ch = chi_arr[resx]; nz = resx != 0
    out = np.zeros(len(resx), complex)
    if j not in (0, 4):
        e = -j if wrong else -j-2
        out[nz] = ch[nz]**e
    elif j == 0:
        out[nz] = NPI**-0.5*ch[nz]**-2
    else:
        out[:] = NPI**-0.5*(-1 + NPI*(~nz))
    return out
gam = {j: R.gauss(PI, j) for j in range(6)}
def varpi(j, eps):
    ce_ = complex(R.sext_sym_prime(eps, PI))
    if j == 4: return gam[4]
    if j == 0: return -gam[2]*ce_**-2
    return complex(R.sext_sym_prime((-1, 0), PI))**j*gam[j]*gam[(j+2) % 6]*ce_**(-j-2)
variants = ['main', 'alpha_c_unconj', 'abar_c^2', 'zeta_-i', 'kappa_unconj', 'B_-j', 'k1', 'left_abar3']
tot = {j: {v: 0j for v in variants} for j in JS}
ngroups = 0; types = {}
for h0 in RES72:
    ph = phihat(h0)
    if abs(ph) < 1e-13: continue
    for Akind in ('inactive', 'active'):
        A = [PI] if Akind == 'active' else []
        gd = group_data(h0, A); ngroups += 1
        c = gd['c']; cc = fg.cx(c); NC = fg.N(c)
        types[(gd['sigma'], NC)] = types.get((gd['sigma'], NC), 0) + 1
        Lr, base2, base1, resx = dual_base(gd['sigma'], NC)
        psi = e_add(fg.cx(gd['beta'])*LAMc*Lr/fg.cx(gd['c0']))          # psi(lambda^4 l) = e(beta lambda l / c0)
        ac = cc/abs(cc)
        for j in JS:
            if Akind == 'inactive' and j != 0: continue                    # inactive branch only for j = 0
            if Akind == 'inactive':
                Bf = np.ones(len(Lr)); Bw = Bf; scal = (1 - 1/NPI)
            else:
                Bf = Bvec(j, resx); Bw = Bvec(j, resx, wrong=True)
                scal = complex(R.sext_sym_prime(gd['sigpi'], PI))**-2*varpi(j, gd['eps'])
            S2 = np.sum(base2*psi*Bf); S1 = np.sum(base1*psi*Bf); S2w = np.sum(base2*psi*Bw)
            k0 = gd['kap0']
            pre = ph*scal/81
            tot[j]['main'] += -np.conj(ac)**4*pre*np.conj(k0)*S2
            tot[j]['alpha_c_unconj'] += -ac**4*pre*np.conj(k0)*S2
            tot[j]['abar_c^2'] += -np.conj(ac)**2*pre*np.conj(k0)*S2
            tot[j]['zeta_-i'] += -1j*np.conj(ac)**4*pre*np.conj(k0)*S2
            tot[j]['kappa_unconj'] += -np.conj(ac)**4*pre*k0*S2
            tot[j]['B_-j'] += -np.conj(ac)**4*pre*np.conj(k0)*S2w
            tot[j]['k1'] += -np.conj(ac)**4*pre*np.conj(k0)*S1
            tot[j]['left_abar3'] += -np.conj(ac)**4*pre*np.conj(k0)*S2
print(f"groups {ngroups}; cusp types (sigma, N(c)): {dict(sorted(types.items()))}")
worst = 0; ctrl = {v: 0 for v in variants if v != 'main'}
for j in JS:
    err = abs(tot[j]['main']-Tdir[j])/abs(Tdir[j]); worst = max(worst, err)
    for v in ctrl:
        ref = Tdir_b3[j] if v == 'left_abar3' else Tdir[j]
        ctrl[v] = max(ctrl[v], abs(tot[j][v]-ref)/abs(ref))
    print(f"j={j}: T direct = {Tdir[j]:.12f}   T via Proposition (grouped) = {tot[j]['main']:.12f}   rel.err {err:.2e}")
print("controls (max rel.err over j, must be >= 1e-2):", {k: f"{v:.2e}" for k, v in ctrl.items()})
# Exit status (added 2026-10-12): 1 unless the PASS criterion below holds.  CM_TESTS_FORCE_FAIL=1 forces FAIL (self-test).
passed = worst < 1e-8 and min(ctrl.values()) > 1e-2
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    passed = False; print("CM_TESTS_FORCE_FAIL=1: the check below is counted as failed")
print("PASS" if passed else "FAIL", "time", round(time.time()-t0, 1))
sys.exit(0 if passed else 1)
