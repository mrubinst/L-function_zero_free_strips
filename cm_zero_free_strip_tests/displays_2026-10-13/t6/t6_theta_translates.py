"""t6 (2026-10-13): direct pointwise checks of the theta-side displays of the proof of Proposition prop:s1-reflection
(Appendix app:s1-general, Section sec:s1-proof), at the 252 translates of the s1 end-to-end setting (L = 18, base
twist A, one local prime pi of norm 7) and at extra cusps a/c with lambda | c.

 T1  eq:d129   finite Fourier inversion: sum_h C(h) ce(l z_h) = phi(lambda^3 l) chi_pi(lambda^3 l)^j for all
               l in lambda^-3 O (all x = lambda^3 l mod 18 pi), j = 0..5, with z_h = lambda^2 (h0/18 + h_pi/pi),
               C(h) = phihat(h0) C_{pi,j}(h_pi); and the constant term sum_h C(h) = 0.
               Falsifiers: z_h without lambda^2; C_{pi,j} with e(+h y).
 T2  eq:d130   cusp expansion: conj(theta(H(z',v'))) computed from the identity-cusp expansion at the point H(z',v')
               equals C_H(v') + sum_l d_sigma(l) v' K_{1/3}(4 pi |l| v') ce(l z') (d_sigma from the tables; C_H =
               (3^{5/2}/2) v'^{2/3} for sigma = 0, else 0), for every matrix H occurring.  Falsifier: d_+ <-> d_-.
 T3  eq:s1-automorphy (with eq:d130): conj(theta(a/c+z, v)) [identity expansion] = conj(kappa(G)) *
               [C_H(v') + sum_l d_sigma(l) v'K ce(l z')] at (z',v') = g^{-1}(a/c+z, v), g = (a b'; c delta').
               Falsifiers: kappa(G) instead of its conjugate (where kappa != 1); d_+ <-> d_-.
 T4  eq:s1-cusp and eq:d131: d_zbar^2 conj(theta(a/c+z,v))|_{z=0} = conj(kappa(G)) (cv)^{-4} G_h(1/(N(c)v)),
               G_h(v') = -4 pi^2 sum_l d_sigma(l) N(l) alpha(l)^2 v'K_{1/3}(4 pi|l|v') ce(-delta' l/c);  the left side is
               the termwise derivative sum_mu d_0(mu)(2 pi i conj mu)^2 vK ce(mu a/c).
               Falsifiers: (cv)^{-2} (one-derivative power); conj(alpha(l))^2 in G_h.
 T5  eq:d132   the majorant: every table entry l = u lambda^k n b^3 has |d_sigma(l)| <= 27 3^{k/6} |b| (so the
               termwise majorant holds), at most 6 units per (k,n,b), and the partial sums of the left side stay
               below the (finite) right side at Re s' = 2 and 1.6.  Falsifier: constant 9 instead of 27.
 T6  eq:s1-fe (with eq:s1-dual, eq:d133): for Re s < 1/2, J_2(s) = int_0^oo f(v) v^{2s-1} dv computed by quadrature
               (f from the identity expansion of F for v >= v0 and from the 252 cusp expansions eq:s1-cusp for v < v0,
               two values of v0) equals sum_h C(h) conj(kappa) conj(alpha(c))^4 N(c)^{2-2s} J^vee_h(2-s), J^vee from
               eq:s1-dual (dual Dirichlet series).  Falsifiers: conj(alpha(c))^2; N(c)^{1-2s}.
               The same quadrature at Re s > 3/2 is compared with eq:s1-direct (D_2 L_2 from Dirichlet series).
Falsifier deviations are reported as the largest deviation over all cases (a falsifier is detected when the
check with the wrong formula would fail).
Usage: python3 t6_theta_translates.py [pi_index]      (about 1-2 minutes)
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
from t6_lib import *
from zw import prim_assoc
from scipy.special import kv, gamma as sgamma, loggamma
import mpmath as mp
t0 = time.time()
PIDX = int(sys.argv[1]) if len(sys.argv) > 1 else 0
S = Setup(25000, PIDX)
CHECKS = []
def check(name, ok, msg):
    CHECKS.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)
SIG0 = 3**2.5/2
NMAX = 300.0                     # norms of indices l kept in the expansions (heights >= 0.18)
TAB = {s: S.table(s, NMAX) for s in ('0', '+', '-')}
L0, D0, X0 = TAB['0']
print(f"tables: d_0 {len(L0)}, d_+ {len(TAB['+'][0])}, d_- {len(TAB['-'][0])} entries; {time.time()-t0:.1f} s")

def theta_bar(z, v):             # identity-cusp expansion of conj(theta)
    return SIG0*v**(2/3)+np.sum(D0*v*kv(1/3, 4*np.pi*np.abs(L0)*v)*ce(L0*z))
def cusp_series(sigma, z, v, swap=False):
    s = sigma if not swap else {'0': '0', '+': '-', '-': '+'}[sigma]
    Lr, Dr, _ = TAB[s]
    return (SIG0*v**(2/3) if sigma == '0' else 0)+np.sum(Dr*v*kv(1/3, 4*np.pi*np.abs(Lr)*v)*ce(Lr*z))
def cplx(x): return complex(cx(x))

# ---------------- T1: finite Fourier inversion eq:d129 ----------------
RX = residues(mul(S.L18, S.PI))
xs = np.array([cplx(x) for x in RX]); ls = xs/LAMc**3
h0s = np.array([cplx(h) for h in S.RES]); phs = np.array([S.PHIHAT[h] for h in S.RES])
hps = np.array([cplx(h) for h in S.RESP])
worst = 0; f1 = f2 = 0; csum = 0
for j in range(6):
    Cp = np.array([S.Cpj(h, j) for h in S.RESP])
    Cp_wrong = np.array([sum(S.chipow(S.red(y), j)*e_add(+cplx(h)*cplx(y)/S.PIc) for y in S.RESP)/7 for h in S.RESP])
    zh = LAMc**2*(h0s[:, None]/18+hps[None, :]/S.PIc)                    # 324 x 7 translates
    Ch = phs[:, None]*Cp[None, :]
    lhs = np.array([np.sum(Ch*ce(l*zh)) for l in ls])
    rhs = np.array([S.phi(x)*S.chipow(S.red(x), j) for x in RX])
    worst = max(worst, np.max(np.abs(lhs-rhs)))
    lw1 = np.array([np.sum(Ch*ce(l*zh/LAMc**2)) for l in ls]); f1 = max(f1, np.max(np.abs(lw1-rhs)))
    Chw = phs[:, None]*Cp_wrong[None, :]
    if j not in (0, 3):
        lw2 = np.array([np.sum(Chw*ce(l*zh)) for l in ls]); f2 = max(f2, np.max(np.abs(lw2-rhs)))
    csum = max(csum, abs(np.sum(Ch)))
check("T1 eq:d129 Fourier inversion", worst < 1e-12 and f1 > 0.1 and f2 > 0.1 and csum < 1e-13,
      f"{len(RX)} classes x mod 18pi, j=0..5: max err {worst:.1e}; |sum_h C(h)| = {csum:.1e}; "
      f"falsifiers: z_h without lambda^2 {f1:.2f}, e(+hy) in C_pi,j {f2:.2f}")

# ---------------- T2: cusp expansion eq:d130 for every H that occurs ----------------
def matc(M): return tuple(cplx(e) for row in M for e in row)
Hs = {}
for t in S.TR.values(): Hs[(t['H'], t['sigma'])] = 1
EXTRA = []                                   # cusps with lambda | c (sigma from v_lambda(c) = 1, 2, 3)
for (a, c) in [((1, 0), LAMp), ((-2, 0), LAMp), ((1, 0), mul(LAMp, LAMp)), ((4, 0), mul(LAMp, LAMp)), ((1, 3), pw(LAMp, 3)),
               ((1, 0), mul((-2, 0), LAMp)), ((7, 3), mul((-2, -3), LAMp)), ((4, 3), mul((-2, 0), mul(LAMp, LAMp)))]:
    g = S.gcd_O(a, c); assert N(g) == 1, (a, c)
    a2, c2 = S.normalize(a, c)
    a2 = add(reduce_mod(sub(a2, (1, 0)), mul((3, 0), c2)), (1, 0))
    sigma, dp, kap, bp, H = S.cusp_data(a2, c2)
    EXTRA.append(dict(a=a2, c=c2, sigma=sigma, dp=dp, kap=kap, bp=bp, H=H)); Hs[(H, sigma)] = 1
worst = 0; fsw = 0; npts = 0
for (H, sigma) in Hs:
    Hc = matc(H)
    for zp, vp in [(0.05+0.11j, 0.62), (-0.13+0.02j, 0.71), (0.21-0.07j, 0.55)]:
        z1, v1 = act(Hc, zp, vp)
        lhs = theta_bar(z1, v1); rhs = cusp_series(sigma, zp, vp)
        worst = max(worst, abs(lhs-rhs)/abs(lhs)); npts += 1
        if sigma != '0': fsw = max(fsw, abs(lhs-cusp_series(sigma, zp, vp, swap=True))/abs(lhs))
check("T2 eq:d130 cusp expansions", worst < 1e-9 and fsw > 1e-2,
      f"{len(Hs)} matrices H (sigma in {sorted(set(s for _, s in Hs))}), {npts} points: max rel err {worst:.1e}; "
      f"falsifier d_+ <-> d_-: {fsw:.2e}")

# ---------------- T3, T4: automorphy and second derivative at the translates ----------------
def ginv_point(t, z, v):
    a, c, bp, dp = cplx(t['a']), cplx(t['c']), cplx(t['bp']), cplx(t['dp'])
    return act((dp, -bp, -c, a), a/c+z, v)
w3 = w3e = 0; fk = fs = 0; w4 = 0; f4a = f4b = 0; nk = 0
items = list(S.TR.values())+EXTRA
for t in items:
    a, c = cplx(t['a']), cplx(t['c']); Nc = abs(c)**2; kap = t['kap']
    v = Nc**-0.5
    for z in [0, 0.03+0.02j]:
        zp, vp = ginv_point(t, z, v)
        assert min(v, vp) > 0.15, (v, vp)
        lhs = theta_bar(a/c+z, v)
        rhs = np.conj(kap)*cusp_series(t['sigma'], zp, vp)
        err = abs(lhs-rhs)/abs(lhs)
        if t in EXTRA: w3e = max(w3e, err)
        else: w3 = max(w3, err)
        if abs(kap-1) > 1e-6:
            fk = max(fk, abs(lhs-kap*cusp_series(t['sigma'], zp, vp))/abs(lhs)); nk += 1
        if t['sigma'] != '0':
            fs = max(fs, abs(lhs-np.conj(kap)*cusp_series(t['sigma'], zp, vp, swap=True))/abs(lhs))
    # second derivative at z = 0 (eq:s1-cusp, eq:d131)
    lhs = np.sum(D0*(2j*np.pi*np.conj(L0))**2*v*kv(1/3, 4*np.pi*np.abs(L0)*v)*ce(L0*a/c))
    Lr, Dr, _ = TAB[t['sigma']]; vp = 1/(Nc*v); dpc = cplx(t['dp'])
    base = -4*np.pi**2*Dr*np.abs(Lr)**2*vp*kv(1/3, 4*np.pi*np.abs(Lr)*vp)*ce(-dpc*Lr/c)
    al2 = (Lr/np.abs(Lr))**2
    G = np.sum(base*al2); Gw = np.sum(base*np.conj(al2))
    rhs = np.conj(kap)*(c*v)**-4*G
    w4 = max(w4, abs(lhs-rhs)/abs(lhs))
    f4a = max(f4a, abs(lhs-np.conj(kap)*(c*v)**-2*G)/abs(lhs))
    f4b = max(f4b, abs(lhs-np.conj(kap)*(c*v)**-4*Gw)/abs(lhs))
check("T3 eq:s1-automorphy (+eq:d130)", max(w3, w3e) < 1e-9 and fk > 1e-2 and fs > 1e-2,
      f"{len(S.TR)} translates (sigma=0,+,-; N(c)=1,4,7,28) and {len(EXTRA)} cusps with lambda|c, 2 points each: "
      f"max rel err {w3:.1e} / {w3e:.1e}; falsifiers kappa (not conj, {nk} cases) {fk:.2f}, d_+<->d_- {fs:.2f}")
check("T4 eq:s1-cusp, eq:d131", w4 < 1e-9 and f4a > 1e-2 and f4b > 1e-2,
      f"{len(items)} cusps: max rel err {w4:.1e}; falsifiers (cv)^-2: {f4a:.2f}, conj(alpha(l))^2: {f4b:.2f}")

# ---------------- T5: the majorant eq:d132 ----------------
viol = violf = 0; maxratio = 0; cnt = 0; units = {}; supp_bad = 0
# metadata of each entry l = u lambda^k c d^3 (c squarefree primary, d primary): then n = c, b = d (n b^3 = c d^3)
for sigma in ('0', '+', '-'):
    Lr, Dr, XE, MK = S.table(sigma, NMAX, meta=True)
    for l, d, xe, (k, nc, nd) in zip(Lr, Dr, XE, MK):
        if abs(d) < 1e-12: continue
        if k < -4 or abs(abs(l)**2-3.0**k*nc*nd**3) > 1e-9*abs(l)**2: supp_bad += 1
        bound = 27*3**(k/6)*math.sqrt(nd)
        cnt += 1; maxratio = max(maxratio, abs(d)/bound)
        if abs(d) > bound*(1+1e-9): viol += 1
        if abs(d) > 9*3**(k/6)*math.sqrt(nd)*(1+1e-9): violf += 1
        x = (int(xe[0]), int(xe[1]))
        for _ in range(int(k)+4): x = div(x, LAMp)
        key = (sigma, int(k), prim_assoc(x))
        units[key] = units.get(key, 0)+1
maxu = max(units.values())
sums = {}
for sp in (2.0, 1.6):
    lhs = 0
    for sigma in ('0', '+', '-'):
        Lr, Dr, _ = TAB[sigma]; lhs = max(lhs, np.sum(np.abs(Dr)*np.abs(Lr)**(2*(0.5-sp))))
    # right side, each sum over all primary elements (n squarefree) up to norm 2e5 plus an integral tail bound
    sk = sum(3**(k*(2/3-sp)) for k in range(-4, 400))
    nn = [N(c) for c, _ in squarefree_list(200000, primary_primes(200000))]
    sn = sum(x**(0.5-sp) for x in nn)
    bb = [N((x, y)) for x in range(-600, 601) for y in range(-600, 601) if primary((x, y)) and N((x, y)) <= 200000]
    sbv = sum(x**(2-3*sp) for x in bb)
    rhs = 6*27*sk*sn*sbv
    sums[sp] = (lhs, rhs)
okmaj = all(l <= r for l, r in sums.values())
check("T5 eq:d132 majorant", viol == 0 and violf > 0 and maxu <= 6 and okmaj and supp_bad == 0,
      f"{cnt} nonzero entries (all sigma), support u lambda^k n b^3 (k>=-4, n sf) violated {supp_bad} times; |d| <= 27 3^(k/6)|b| violated {viol} times, max ratio {maxratio:.4f}; units per (k,n,b) <= {maxu}; "
      + "; ".join(f"Re s'={sp}: partial LHS (one sigma) {l:.4g} <= RHS (truncated, finite) {r:.4g}" for sp, (l, r) in sums.items())
      + f"; falsifier constant 9: {violf} violations")
print(f"[time {time.time()-t0:.1f} s]")

# ---------------- T6: eq:s1-fe and eq:s1-direct by quadrature of f(v) ----------------
# F-coefficients at the identity cusp: d_0(mu) phi(lambda^3 mu) chi_pi(lambda^3 mu)^j (by T1); f(v) = d_zbar^2 F at z=0.
Lf, Df, Xf = S.table_d0(900.0)                                   # complete: 27*900 <= B
lam3mu = [(int(x[0]), int(x[1])) for x in Xf]                   # lambda^4 mu (exact); lambda^3 mu = that / lambda
lam3mu = [div(x, LAMp) if divides(LAMp, x) else None for x in lam3mu]
phimask = np.array([0 if x is None else S.phi(x) for x in lam3mu], float)
redx = np.array([0 if x is None else S.red(x) for x in lam3mu])
sel = phimask != 0
Lf, Df, redx = Lf[sel], Df[sel], redx[sel]
order = np.argsort(np.abs(Lf)); Lf, Df, redx = Lf[order], Df[order], redx[order]; aLf = np.abs(Lf)
def f_id(vv, j):
    co = Df*S.chipow(redx, j)*(2j*np.pi*np.conj(Lf))**2
    out = []
    for v in vv:
        k = np.searchsorted(aLf, 60/(4*np.pi*v))
        assert k < len(aLf) or 4*np.pi*aLf[-1]*v > 45
        out.append(np.sum(co[:k]*v*kv(1/3, 4*np.pi*aLf[:k]*v)))
    return np.array(out)
# cusp side: group the translates by (sigma, N(c)); A_class(l) = sum_h C(h) conj(kappa) c^-4 (-4pi^2) d N alpha^2 ce(-delta' l/c)
def cusp_classes(j, alpha_pow=4):
    cls = {}
    for (h0, hp), t in S.TR.items():
        Ch = S.PHIHAT[h0]*S.Cpj(hp, j)
        if abs(Ch) < 1e-15: continue
        c = cplx(t['c']); Nc = round(abs(c)**2); Lr, Dr, _ = TAB[t['sigma']]
        key = (t['sigma'], Nc)
        vec = Ch*np.conj(t['kap'])*c**-4*(-4*np.pi**2)*Dr*np.abs(Lr)**2*(Lr/np.abs(Lr))**2*ce(-cplx(t['dp'])*Lr/c)
        cls[key] = cls.get(key, 0)+vec
    out = {}
    for key, A in cls.items():
        aL = np.abs(TAB[key[0]][0]); o = np.argsort(aL); out[key] = (A[o], aL[o])
    return out
def f_cusp(vv, cls):
    out = np.zeros(len(vv), complex)
    for (sigma, Nc), (A, aL) in cls.items():
        for i, v in enumerate(vv):
            vp = 1/(Nc*v); k = np.searchsorted(aL, 60/(4*np.pi*vp))
            out[i] += v**-4*np.sum(A[:k]*vp*kv(1/3, 4*np.pi*aL[:k]*vp))
    return out
xg, wg = np.polynomial.legendre.leggauss(400)
def quad_log(fun, lo, hi):                 # int_lo^hi g(v) dv/v... returns nodes and weights in v for int dv
    u = 0.5*(np.log(hi)-np.log(lo))*xg+0.5*(np.log(hi)+np.log(lo)); w = 0.5*(np.log(hi)-np.log(lo))*wg
    v = np.exp(u); return v, w*v
def J2_quad(s, j, v0, cls):
    v1, w1 = quad_log(None, v0, 25.0); v2, w2 = quad_log(None, 4e-4, v0)
    return np.sum(f_id(v1, j)*v1**(2*s-1)*w1)+np.sum(f_cusp(v2, cls)*v2**(2*s-1)*w2)
def J2_fe(s, j, apow=4, npow=2):
    sp = 2-s; tot = 0
    for (h0, hp), t in S.TR.items():
        Ch = S.PHIHAT[h0]*S.Cpj(hp, j)
        if abs(Ch) < 1e-15: continue
        c = cplx(t['c']); Nc = abs(c)**2; Lr, Dr, _ = TAB[t['sigma']]
        dual = -(np.pi/2)*(2*np.pi)**(-2*sp)*complex(mp.gamma(sp+1/3)*mp.gamma(sp+2/3))*np.sum(
            Dr*(Lr/np.abs(Lr))**2*np.abs(Lr)**(2*(0.5-sp))*ce(-cplx(t['dp'])*Lr/c))
        tot += Ch*np.conj(t['kap'])*np.conj(c/abs(c))**apow*Nc**(npow-2*s)*dual
    return tot
w6 = 0; f6a = f6b = 0; v0dep = 0
for j in (0, 1, 4):
    cls = cusp_classes(j)
    for s in (-2.0, -2.5+1.5j):
        rhs = J2_fe(s, j)
        lA = J2_quad(s, j, 0.22, cls); lB = J2_quad(s, j, 0.15, cls)
        w6 = max(w6, abs(lA-rhs)/abs(rhs), abs(lB-rhs)/abs(rhs)); v0dep = max(v0dep, abs(lA-lB)/abs(lA))
        f6a = max(f6a, abs(lA-J2_fe(s, j, apow=2))/abs(lA)); f6b = max(f6b, abs(lA-J2_fe(s, j, npow=1.5))/abs(lA))
    print(f"   j={j}: eq:s1-fe at s=-2, -2.5+1.5i done ({time.time()-t0:.1f} s)", flush=True)
check("T6a eq:s1-fe (+eq:s1-dual)", w6 < 1e-8 and f6a > 1e-2 and f6b > 1e-2,
      f"j=0,1,4, s=-2,-2.5+1.5i, v0=0.22,0.15: max rel err {w6:.1e} (v0-dependence {v0dep:.1e}); "
      f"falsifiers conj(alpha(c))^2: {f6a:.2f}, N(c)^(1-2s): {f6b:.2f}")
# eq:s1-direct at Re s > 3/2: D_2 and L_2 as separate Dirichlet series in w = s - 1/2
def gam2_fast(cc_, g1_, slam): return slam*g1_/abs(cplx(cc_))
DN = []; DC = []; DR = []; DC1 = []
for cc_, nc, g1_, _, _, _, slam, even, haspi in S.cdata:
    if even: continue
    cc = cplx(cc_); xin = 1
    for p in S.sqf[cc_]: xin *= S.psi0(p)
    DN.append(nc); DC.append(gam2_fast(cc_, g1_, slam)*np.conj(cc/abs(cc))**2*xin); DR.append(S.red(cc_))
    DC1.append(gam2_fast(cc_, g1_, slam)*np.conj(cc/abs(cc))*xin)
DN, DC, DR, DC1 = np.array(DN, float), np.array(DC), np.array(DR), np.array(DC1)
BN = []; BA = []; BR = []
for x in range(-400, 401):
    for y in range(-400, 401):
        b = (x, y)
        if not primary(b) or (x % 2 == 0 and y % 2 == 0) or N(b) > 160000: continue
        BN.append(N(b)); BA.append(cplx(b)/abs(cplx(b))); BR.append(S.red(b))
BN, BA, BR = np.array(BN, float), np.array(BA), np.array(BR)
def D2(w, j, one=False): return np.sum((DC1 if one else DC)*S.chipow(DR, j)*DN**(-w))                    # Psi(n) = Psi_0(n) chi_pi(n)^j
def L2(w, j, apow=6): return np.sum(np.conj(BA)**apow*S.chipow(BR, (3*j) % 6)*BN**(-3*w+0.5))   # Psi_0(b)^3 = 1
w7 = 0; f7 = 0; f7b = 0
for j in (0, 1, 4):
    cls = cusp_classes(j)
    for s in (4.0, 3.5+2j):
        lhs = J2_quad(s, j, 0.22, cls)
        T2 = D2(s-0.5, j)*L2(s-0.5, j)
        pre = (3*np.pi/2)*(27/(2*np.pi)**2)**s*complex(mp.gamma(s+1/3)*mp.gamma(s+2/3))
        w7 = max(w7, abs(lhs-pre*T2)/abs(lhs))
        f7b = max(f7b, abs(lhs-pre*D2(s-0.5, j, one=True)*L2(s-0.5, j, apow=3))/abs(lhs))
        f7 = max(f7, abs(lhs-(3**2.5/4)*(27/(2*np.pi)**2)**s*complex(mp.gamma(s+1/3)*mp.gamma(s+2/3))*T2)/abs(lhs))
check("T6b eq:s1-direct by quadrature", w7 < 1e-7 and f7 > 1e-2 and f7b > 1e3*max(w7, 1e-12),
      f"j=0,1,4, s=4, 3.5+2i: J_2 by quadrature of f vs (3pi/2)(27/(2pi)^2)^s GG D_2 L_2: max rel err {w7:.1e}; "
      f"falsifiers one-derivative constant 3^(5/2)/4: {f7:.2f}, zeta angular powers conj(alpha), conj(alpha)^3: {f7b:.1e}")
# brute-force gamma_2 (definition eq:d036) against the fast values used in D_2
from zw import chi as chi6
bad = 0; nchk = 0
for cc_, nc, g1_, _, _, _, slam, even, haspi in S.cdata:
    if even or nc == 1 or nc > 400: continue
    cc = cplx(cc_); s = 0
    for x in residues(cc_):
        cv = chi6(cc_, x)
        if cv != 0: s += cv**2*e_add(cplx(x)/cc)
    nchk += 1
    if abs(s/math.sqrt(nc)-gam2_fast(cc_, g1_, slam)) > 1e-9: bad += 1
check("T6c gamma_2 brute force", bad == 0 and nchk > 50, f"{nchk} squarefree n (N(n) <= 400, odd): fast gamma_2 = definition eq:d036 ({bad} mismatches)")

nf = sum(1 for _, ok in CHECKS if not ok)
print(f"{len(CHECKS)} checks, {nf} failed; time {time.time()-t0:.1f} s")
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1': nf += 1; print("CM_TESTS_FORCE_FAIL=1")
print("PASS" if nf == 0 else "FAIL")
raise SystemExit(0 if nf == 0 else 1)
