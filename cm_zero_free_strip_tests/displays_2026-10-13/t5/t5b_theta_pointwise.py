"""t5b: pointwise checks of Appendix A with the actual cubic theta function of [QF, Theorem 23].

Displays tested (labels of cm-eleven-twelfths-guide.tex, numbers of the 2026-10-13 build):
  (A.4)  eq:d098       t_sigma, d_sigma: the expansions F_sigma built from them equal conj(theta(gamma_sigma P))
  (A.14) eq:b3-translates  Theta_Psi(z,v) = sum_h c_F(h) conj(theta(z + z_h, v)), pointwise, L = 72; sum c_F = 0;
                       the coefficient identity sum_h c_F(h) ce(l z_h) = phi(x) prod chi_pi(x)^j (x = lambda^3 l)
  (A.15) eq:d105       g in SL_2(O), g1 = g H^{-1} in Gamma(3), theta(HP) = theta(gamma_sigma P) (exact decomposition
                       H gamma_sigma^{-1} = S U, S in SL_2(Z), U in Gamma(3) with kappa(U) = 1, and numerically)
  (A.16) eq:b3-automorphy  conj(theta(a/c+z, v)) = conj(kappa(g1)) F_sigma(g^{-1}(a/c+z, v)), pointwise
  (A.25) eq:d109       g^{-1}(a/c+z, v) by the formula vs the action (5.1)
  (A.28) eq:b3-cusp    d_zbar^2 conj(theta(a/c+z,v))|_{z=0} = conj(kappa(g1)) (cv)^{-4} G(1/(N(c)v))
  (A.29) eq:d112       the expansion of F_sigma (constant term only for sigma = 0) [used in A.16, A.28]
  (A.30) eq:b3-G       G(v') = -4 pi^2 sum d_sigma(l) N(l) alpha(l)^2 v'K(4 pi|l|v') ce(-delta l/c)  [used in A.28]
  (A.31) eq:d113       f(v) = sum_h c_F(h) conj(kappa(g1,h)) (c_h v)^{-4} G_h(1/(N(c_h)v)) (f from the twisted series)
  (A.35) eq:b3-fe      J(s) = sum_h c_F conj(kappa) conj(alpha(c))^4 N(c)^{2-2s} J^vee_h(2-s) at Re s < 1/2, with
                       J(s) = int_0^oo f(v) v^{2s-1} dv computed by quadrature in v
Twists Psi = xi_0 prod chi_pi^{j_pi} (phi(x) = (lambda/x)_3 xi_0(x)):
  T1: phi = 1 (xi_0 = conj((lambda/.)_3)), P empty;  T2: phi = (lambda/x)_3 (xi_0 = 1);  T3: phi = (omega/x)_3;
  T4: phi = psi4 (a character mod 4);  T5_j: phi = 1, P = {pi}, N(pi) = 7, j = 0, 1, 2, 4.
Cusp data: own implementation of [QF Sec. 5.5 (5.18)] in t5lib.py (kappa(g1) = (c1/a1)_3, checked = (b1/d1)_3).
Falsifiers (must fail): theta(HP) vs theta(gamma_{-sigma}P); kappa instead of conj(kappa); (conj(c) v)^{-4};
conj(alpha(l))^2 in G; z - z_h in (A.14); conj(alpha(c))^2 in (A.35).
Usage: python3 -B t5b_theta_pointwise.py      (about 2-3 minutes)
"""
import os, sys, time, math, cmath
t0 = time.time()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t5lib import *
from scipy.special import kv, gamma as cgamma

B = 50000
TT = ThetaTables(B)
RTAB = math.sqrt(B/81.0)*0.999            # every l with |l| <= RTAB is in the tables
TAB = {}
for s in '0+-':
    tb = TT.table(s)
    order = np.argsort(np.abs(tb['l']))
    TAB[s] = dict(l=tb['l'][order], d=tb['d'][order], absl=np.abs(tb['l'][order]),
                  m=[tb['m'][i] for i in order])
CUT = 3.6                                  # keep 4 pi |l| v <= 45.2: K_{1/3} below e^-45
SIG0 = 3**2.5/2

def terms(sig, v):
    t = TAB[sig]; R = CUT/v
    assert R <= RTAB, ("table too short", sig, v, R, RTAB)
    n = np.searchsorted(t['absl'], R)
    return t['l'][:n], t['d'][:n]

def F(sig, z, v, deriv=None):
    """F_sigma(z, v) = conj(theta(gamma_sigma(z, v))) from (A.29) (eq:d112); deriv 'zbar2' or 'z2' applies the operator
    d_zbar^2 or d_z^2 termwise."""
    l, d = terms(sig, v)
    val = d*v*kv(1/3, 4*np.pi*np.abs(l)*v)*ce(l*z)
    if deriv == 'zbar2': val = val*(2j*np.pi*np.conj(l))**2
    elif deriv == 'z2': val = val*(2j*np.pi*l)**2
    out = np.sum(val)
    if deriv is None and sig == '0': out += SIG0*v**(2/3)
    return out

def theta(z, v): return np.conj(F('0', z, v))

def relerr(a, b): return abs(a - b)/max(abs(b), 1e-300)

# ---------------------------------------------------------------------------------------- (A.15): theta(HP) = theta(gamma_sigma P)
def find_S(Mx):
    """S in SL_2(Z) with S = Mx mod 3 (entries in -2..2), so that U = S^{-1} Mx is in Gamma(3)."""
    rng = range(-2, 3)
    for a in rng:
        for b in rng:
            for c in rng:
                for d in rng:
                    if a*d - b*c != 1: continue
                    S = (((a, 0), (b, 0)), ((c, 0), (d, 0)))
                    U = matmul(matinv(S), Mx)
                    if in_gamma3(U): return S, U
    return None, None
one, zero = (1, 0), (0, 0)
Hlist = [(((one, zero), (zero, one)), '0', 'I')]
Hlist += [(((one, zero), (LAM, one)), '-', '(1 0; lambda 1)'), (((one, zero), (neg(LAM), one)), '+', '(1 0; -lambda 1)')]
for x in range(3):
    for y in range(3):
        Hlist.append((((((x, y)), (-1, 0)), (one, zero)), {0: '0', 1: '-', 2: '+'}[y], f'(u -1; 1 0), u={x}+{y}w'))
ok_exact = True; worst = 0; fals = 0
pts = [(0.13 + 0.21j, 0.83), (-0.27 + 0.05j, 1.1)]
for H, sig, name in Hlist:
    gs = ((one, zero), (GAMMA[sig][2], one))
    Mx = matmul(H, matinv(gs))
    S, U = find_S(Mx)
    good = S is not None and kubota(U) == 1 and kubota_bd(U) == 1
    ok_exact &= good
    for z, v in pts:
        zH, vH = act(mat_c(H), z, v)
        lhs = theta(zH, vH)
        rhs = np.conj(F(sig, z, v))                                  # theta(gamma_sigma P) = conj(F_sigma(P))
        worst = max(worst, relerr(lhs, rhs))
        other = {'0': '+', '+': '-', '-': '+'}[sig]
        fals = max(fals, relerr(lhs, np.conj(F(other, z, v))))
report("t5b.A15", "eq:d105 theta(HP)=theta(gamma_sigma P)", ok_exact and worst < 1e-10,
       f"{len(Hlist)} matrices H of [QF (5.18)] (I, (1 0; +-lambda 1), (u -1; 1 0) for the 9 u mod 3): H gamma_sigma^-1 = S U "
       f"with S in SL_2(Z), U in Gamma(3), kappa(U) = 1 (exact): {ok_exact}; theta(HP) vs theta(gamma_sigma P) at 2 points, "
       f"max rel err {worst:.1e}")
falsifier("t5b.A15", "label sigma swapped (0->+, +->-, -->+)", fals, 1e-3)

# ---------------------------------------------------------------------------------------- cusps for (A.16), (A.25), (A.28)
cusp_list = []
for c in [(1, 0), (-2, 0), (1, 2), (-1, -2), (-3, 0), (2, 4), (3, 1), (-6, 0), (-2, -3), (4, 3), (5, 3), (-5, 0)]:
    for a in [(0, 0), (1, 0), (0, 1), (2, 1), (1, 3), (-1, 1), (4, 1)]:
        if a == (0, 0) and c != (1, 0): continue
        if a != (0, 0) and gcd_O(a, c) not in UNITS: continue
        aa, cc = lowest_terms(a, c) if a != (0, 0) else ((0, 0), (1, 0))
        if (aa, cc) in [(x['a'], x['c']) for x in cusp_list]: continue
        cusp_list.append(cusp_data(aa, cc))
assert all(cd['inG3'] for cd in cusp_list)
kap_agree = all(abs(cd['kappa'] - cd['kappa_bd']) < 1e-12 for cd in cusp_list)
types = sorted(set((cd['sigma'], cd['vl'], N(cd['c'])) for cd in cusp_list))
print(f"   {len(cusp_list)} cusps a/c; (sigma, v_lambda(c), N(c)) types: {types}; kappa values "
      f"{sorted(set(round(cmath.phase(cd['kappa']), 3) for cd in cusp_list))}; (c1/a1)_3 = (b1/d1)_3: {kap_agree}", flush=True)

# (A.25) eq:d109 numerically
e25 = 0
for cd in cusp_list:
    ac = complex(cx(cd['a'])); cc = complex(cx(cd['c'])); dl = complex(cx(cd['delta']))
    ginv = mat_c(matinv(cd['g']))
    for z in [0.031 + 0.017j, -0.05 + 0.08j]:
        for v in [0.3, 0.9]:
            zp, vp = act(ginv, ac/cc + z, v)
            zf = -dl/cc - np.conj(z)/(cc**2*(v*v + abs(z)**2)); vf = v/(abs(cc)**2*(v*v + abs(z)**2))
            e25 = max(e25, abs(zp - zf) + abs(vp - vf))
report("t5b.A25", "eq:d109 g^{-1}(a/c+z,v)", e25 < 1e-12 and kap_agree,
       f"{len(cusp_list)} cusps x 4 points: formula vs the action [QF (5.1)] of g^{{-1}}, max abs err {e25:.1e}")

# (A.16) + (A.29): conj(theta(a/c+z, v)) = conj(kappa) F_sigma(g^{-1}(a/c+z, v))
e16 = 0; f16k = 0; f16s = 0
for cd in cusp_list:
    ac = complex(cx(cd['a'])); cc = complex(cx(cd['c'])); NC = abs(cc)**2
    ginv = mat_c(matinv(cd['g'])); kap = cd['kappa']
    for z in [0.04 + 0.02j, -0.03 + 0.05j]:
        v = 1/math.sqrt(NC) if NC > 1 else 0.7
        lhs = F('0', ac/cc + z, v)
        zp, vp = act(ginv, ac/cc + z, v)
        Fs = F(cd['sigma'], zp, vp)
        e16 = max(e16, relerr(np.conj(kap)*Fs, lhs))
        f16k = max(f16k, relerr(kap*Fs, lhs))
        other = {'0': '+', '+': '-', '-': '+'}[cd['sigma']]
        f16s = max(f16s, relerr(np.conj(kap)*F(other, zp, vp), lhs))
report("t5b.A16", "eq:b3-automorphy (with eq:d112)", e16 < 1e-9,
       f"{len(cusp_list)} cusps, 2 points each (v = N(c)^-1/2): max rel err {e16:.1e}; F_sigma from (A.29) with d_sigma of (A.4)")
falsifier("t5b.A16", "kappa instead of conj(kappa)", f16k, 1e-3)
falsifier("t5b.A16", "wrong cusp label", f16s, 1e-3)

# (A.28) + (A.30): d_zbar^2 conj(theta(a/c+z,v))|_{z=0} = conj(kappa) (cv)^{-4} G(1/(N(c)v))
def G(sig, vp, delta, c, variant=None):
    """(A.30) eq:b3-G."""
    l, d = terms(sig, vp)
    ang = (l/np.abs(l))**2 if variant != 'conj' else np.conj(l/np.abs(l))**2
    return -4*np.pi**2*np.sum(d*np.abs(l)**2*ang*vp*kv(1/3, 4*np.pi*np.abs(l)*vp)*ce(-delta*l/c))
e28 = 0; f28c = 0; f28a = 0; f28k = 0
for cd in cusp_list:
    ac = complex(cx(cd['a'])); cc = complex(cx(cd['c'])); NC = abs(cc)**2; dl = complex(cx(cd['delta'])); kap = cd['kappa']
    for v in ([1/math.sqrt(NC)] if NC > 1 else [0.7, 0.45]):
        lhs = F('0', ac/cc, v, deriv='zbar2')
        Gv = G(cd['sigma'], 1/(NC*v), dl, cc)
        rhs = np.conj(kap)*(cc*v)**-4*Gv
        e28 = max(e28, relerr(rhs, lhs))
        f28c = max(f28c, relerr(np.conj(kap)*(np.conj(cc)*v)**-4*Gv, lhs))
        f28a = max(f28a, relerr(np.conj(kap)*(cc*v)**-4*G(cd['sigma'], 1/(NC*v), dl, cc, 'conj'), lhs))
        f28k = max(f28k, relerr(kap*(cc*v)**-4*Gv, lhs))
report("t5b.A28", "eq:b3-cusp, eq:b3-G", e28 < 1e-9,
       f"{len(cusp_list)} cusps: second zbar-derivative of conj(theta) at a/c from the expansion at infinity vs "
       f"conj(kappa)(cv)^-4 G(1/(N(c)v)), G by (A.30): max rel err {e28:.1e}")
falsifier("t5b.A28", "(conj(c) v)^-4", f28c, 1e-3)
falsifier("t5b.A28", "conj(alpha(l))^2 in G", f28a, 1e-3)
falsifier("t5b.A28", "kappa unconjugated", f28k, 1e-3)

# ---------------------------------------------------------------------------------------- twists, (A.14), (A.31), (A.35)
PI = [p for p in fg.primary_primes(7) if N(p) == 7][0]
def phi_base(x):
    """1 on primary x coprime to 6, else 0."""
    return 1 if prime_coprime6(x) else 0
TWISTS = [
    ('T1 phi=1', lambda x: np.conj(lam3_euler(x)), None),
    ('T2 phi=(lambda/x)_3', lambda x: 1, None),
    ('T3 phi=(omega/x)_3', lambda x: np.conj(lam3_euler(x))*om3(x), None),
    ('T4 phi=psi4', lambda x: np.conj(lam3_euler(x))*psi4(x), None),
] + [(f'T5 phi=1, P={{pi7}}, j={j}', lambda x: np.conj(lam3_euler(x)), j) for j in (0, 1, 2, 4)]

L72 = 72
RES72 = [(a, b) for a in range(72) for b in range(72)]
XC72 = np.array([complex(cx(x)) for x in RES72])
XA = np.array([x[0] for x in RES72]); XB = np.array([x[1] for x in RES72])
def phihat_all(phi_vals):
    """phihat(h0) = N(72)^-1 sum_x phi(x) e(-h0 x/72) for all h0 (A.13), exact phases: e(y/72) = exp(2 pi i Q/72),
    Q = omega-coordinate of y (y = h0 x)."""
    out = np.zeros(len(RES72), complex)
    for i0 in range(0, len(RES72), 400):
        ha = XA[i0:i0 + 400, None]; hb = XB[i0:i0 + 400, None]
        Q = (XA[None, :]*hb + XB[None, :]*(ha - hb)) % 72
        out[i0:i0 + 400] = (np.exp(-2j*np.pi*Q/72) @ phi_vals)/5184
    return out
# convention check: e(-h0 x/72) by the complex formula agrees with the exact phase exp(-2 pi i Q/72)
econv = 0
for i, k in [(5, 17), (777, 2500), (3001, 4000), (5183, 1)]:
    Q = (XA[k]*XB[i] + XB[k]*(XA[i] - XB[i])) % 72
    econv = max(econv, abs(e_add(-XC72[i]*XC72[k]/72) - np.exp(-2j*np.pi*Q/72)))
RESPI = sorted(set(canon_mod((a, b), PI) for a in range(7) for b in range(7)))
def Ccoef(j, h):
    """C_{pi,j}(h) = N(pi)^-1 sum_{y mod pi} chi_pi(y)^j e(-h y/pi), chi^0 = 1_{pi nmid y} (A.13)."""
    s = 0
    for y in RESPI:
        ch = sext_prime(y, PI)
        if ch == 0: continue
        s += ch**j*e_add(-complex(cx(mul(h, y)))/complex(cx(PI)))
    return s/N(PI)

def twisted_coeff(xi0, j, x):
    """phi(x) prod chi_pi(x)^j for x = lambda^3 l (the twist of (A.11), eq:d103)."""
    if not prime_coprime6(x): return 0
    val = lam3_euler(x)*xi0(x)
    if j is not None:
        ch = sext_prime(x, PI)
        val *= ch**j if ch != 0 else 0
    return val

tw_summary = []
e14c = 0; e14p = 0; f14 = 0; s14 = 0; e31 = 0; f31k = 0; f31c = 0; e35 = 0; f35a = 0; f35n = 0; q35 = 0
for name, xi0, j in TWISTS:
    t1 = time.time()
    phi_vals = np.array([(lam3_euler(x)*xi0(x) if prime_coprime6(x) else 0) for x in RES72])
    ph = phihat_all(phi_vals)
    # translates: h = (h0, h_pi); c_F(h) = phihat(h0) C_{pi,j}(h_pi) (A.14)
    hpis = [(0, 0)] if j is None else RESPI
    Cv = {h: (1.0 if j is None else Ccoef(j, h)) for h in hpis}
    trans = []
    for i, h0 in enumerate(RES72):
        if abs(ph[i]) < 1e-12: continue
        for hp in hpis:
            cF = ph[i]*Cv[hp]
            if abs(cF) < 1e-13: continue
            if j is None:
                num = mul(mul(LAM, LAM), h0); den = (72, 0)
            else:
                num = mul(mul(LAM, LAM), add(mul(h0, PI), mul((72, 0), hp))); den = mul((72, 0), PI)
            a, c = lowest_terms(num, den)
            cd = cusp_data(a, c)
            assert cd['inG3']
            zh = complex(cx(LAM))**2*(complex(cx(h0))/72 + (complex(cx(hp))/complex(cx(PI)) if j is not None else 0))
            assert abs(zh - complex(cx(a))/complex(cx(c))) < 1e-9
            trans.append(dict(cF=cF, zh=zh, cd=cd, h0=h0, hp=hp))
    sumcF = abs(sum(t['cF'] for t in trans))
    s14 = max(s14, sumcF)
    # coefficient identity: sum_h c_F(h) e(x (h0/72 + h_pi/pi)) = phi(x) prod chi^j(x), all x with N(x) <= 1500
    xs = [(a_, b_) for a_ in range(-45, 46) for b_ in range(-45, 46) if N((a_, b_)) <= 1500]
    xc = np.array([complex(cx(x)) for x in xs])
    lhs = np.zeros(len(xs), complex)
    for t in trans:
        lhs += t['cF']*ce(xc*t['zh']/LAMc**3)                      # ce(l z_h), l = lambda^-3 x
    rhs = np.array([twisted_coeff(xi0, j, x) for x in xs])
    e14c = max(e14c, np.max(np.abs(lhs - rhs)))
    # pointwise (A.14): Theta_Psi(z, v) = sum_h c_F(h) conj(theta(z + z_h, v))
    tw = np.array([twisted_coeff(xi0, j, div(m, LAM)) if divides(LAM, m) else 0 for m in TAB['0']['m']])
    for z, v in [(0.137 + 0.052j, 0.41), (-0.21 + 0.33j, 0.55)]:
        l, d = terms('0', v); n = len(l)
        Th = np.sum(d*tw[:n]*v*kv(1/3, 4*np.pi*np.abs(l)*v)*ce(l*z))
        R = sum(t['cF']*F('0', z + t['zh'], v) for t in trans)
        Rf = sum(t['cF']*F('0', z - t['zh'], v) for t in trans)
        e14p = max(e14p, relerr(R, Th)); f14 = max(f14, relerr(Rf, Th))
    # (A.31): f(v) = d_zbar^2 Theta_Psi(z, v)|_{z=0} vs the cusp sum
    NCmax = max(abs(complex(cx(t['cd']['c'])))**2 for t in trans)
    v0 = 1/math.sqrt(NCmax)
    groups = {}
    for t in trans:
        cd = t['cd']; key = (cd['sigma'], cd['c'])
        groups.setdefault(key, []).append(t)
    GA = {}
    for (sig, c), ts in groups.items():            # per (sigma, c): A(l) = sum_h c_F(h) conj(kappa_h) ce(-delta_h l/c)
        cc = complex(cx(c)); l = TAB[sig]['l']
        A = np.zeros(len(l), complex); Ak = np.zeros(len(l), complex)
        for t in ts:
            ph_ = t['cF']*ce(-complex(cx(t['cd']['delta']))*l/cc)
            A += np.conj(t['cd']['kappa'])*ph_; Ak += t['cd']['kappa']*ph_
        GA[(sig, c)] = (cc, A, Ak)
    def f_direct(v):
        l, d = terms('0', v); n = len(l)
        return np.sum(d*tw[:n]*(2j*np.pi*np.conj(l))**2*v*kv(1/3, 4*np.pi*np.abs(l)*v))
    def f_cusp(v, variant=None):
        tot = 0
        for (sig, c), (cc, A, Ak) in GA.items():
            NC = abs(cc)**2; vp = 1/(NC*v)
            l, d = terms(sig, vp); n = len(l)
            AA = Ak[:n] if variant == 'kappa' else A[:n]
            Gv = -4*np.pi**2*np.sum(AA*d*np.abs(l)**2*(l/np.abs(l))**2*vp*kv(1/3, 4*np.pi*np.abs(l)*vp))
            pre = (np.conj(cc)*v)**-4 if variant == 'cbar' else (cc*v)**-4
            tot += pre*Gv
        return tot
    for v in [v0, 1.15*v0, 0.87*v0]:
        if CUT/v > RTAB or CUT*NCmax*v > RTAB: continue
        fd = f_direct(v); fc = f_cusp(v)
        e31 = max(e31, relerr(fc, fd))
        f31k = max(f31k, relerr(f_cusp(v, 'kappa'), fd)); f31c = max(f31c, relerr(f_cusp(v, 'cbar'), fd))
    # (A.35): J(s) at s = -3 + 0.5 i, by quadrature in log v (cusp sum for v < v0, direct series for v >= v0)
    s = -3 + 0.5j
    def J_quad(step):
        ts = np.arange(-9.0, 2.6, step); tot = 0
        for tt in ts:
            v = math.exp(tt)
            fv = f_cusp(v) if v < v0 else f_direct(v)
            tot += fv*v**(2*s)                       # f(v) v^{2s-1} dv = f(v) v^{2s} dt
        return tot*step
    Jq = J_quad(0.02); Jq2 = J_quad(0.025)
    q35 = max(q35, relerr(Jq2, Jq))
    sp_ = 2 - s
    def Jvee(sig, delta, c, sp):
        tb = TAB[sig]; l = tb['l']; d = tb['d']
        return (-(np.pi/2)*(2*np.pi)**(-2*sp)*cgamma(sp + 1/3)*cgamma(sp + 2/3)
                * np.sum(d*(l/np.abs(l))**2*np.abs(l)**(-2*(sp - 0.5))*ce(-delta*l/c)))
    rhs = 0; rhs_a2 = 0; rhs_n = 0
    for t in trans:
        cd = t['cd']; cc = complex(cx(cd['c'])); NC = abs(cc)**2; ac = cc/abs(cc)
        jv = Jvee(cd['sigma'], complex(cx(cd['delta'])), cc, sp_)
        base = t['cF']*np.conj(cd['kappa'])*jv
        rhs += base*np.conj(ac)**4*NC**(2 - 2*s)
        rhs_a2 += base*np.conj(ac)**2*NC**(2 - 2*s)
        rhs_n += base*np.conj(ac)**4*NC**(1 - 2*s)
    e35 = max(e35, relerr(rhs, Jq)); f35a = max(f35a, relerr(rhs_a2, Jq)); f35n = max(f35n, relerr(rhs_n, Jq))
    tps = sorted(set((t['cd']['sigma'], N(t['cd']['c'])) for t in trans))
    print(f"   {name}: {len(trans)} translates, cusp types (sigma, N(c)) {tps}; |sum c_F| = {sumcF:.1e}; "
          f"(A.31) at v0={v0:.3f}: rel err {relerr(f_cusp(v0), f_direct(v0)):.1e}; (A.35) J(-3+0.5i) quad {Jq:.6e} vs "
          f"fe {rhs:.6e} rel err {relerr(rhs, Jq):.1e}  [{time.time()-t1:.1f}s]", flush=True)
    tw_summary.append((name, len(trans), tps))
report("t5b.A14", "eq:b3-translates", e14c < 1e-11 and e14p < 1e-10 and s14 < 1e-12 and econv < 1e-9,
       f"{len(TWISTS)} twists (L = 72, phihat over all 5184 h0; P = {{pi7}} with j = 0,1,2,4): coefficient identity at all x, "
       f"N(x) <= 1500: max abs err {e14c:.1e}; Theta_Psi(z,v) = sum_h c_F conj(theta(z+z_h,v)) at 2 points: max rel err "
       f"{e14p:.1e}; |sum_h c_F(h)| <= {s14:.1e} (constant terms cancel)")
falsifier("t5b.A14", "z - z_h in place of z + z_h", f14, 1e-3)
report("t5b.A31", "eq:d113", e31 < 1e-9,
       f"f(v) = d_zbar^2 Theta_Psi at z = 0 from the twisted series vs the sum over translates of "
       f"c_F conj(kappa)(c v)^-4 G(1/(N(c)v)), {len(TWISTS)} twists, 2 heights each: max rel err {e31:.1e}")
falsifier("t5b.A31", "kappa unconjugated", f31k, 1e-3)
falsifier("t5b.A31", "(conj(c) v)^-4", f31c, 1e-3)
report("t5b.A35", "eq:b3-fe", e35 < 1e-7 and q35 < 1e-9,
       f"J(s) at s = -3+0.5i: quadrature of int f(v) v^(2s-1) dv (trapezoid in log v; steps 0.02/0.025 agree to {q35:.1e}) "
       f"vs the right side of (A.35) with J^vee from (A.32) over the tables: max rel err {e35:.1e} over {len(TWISTS)} twists")
falsifier("t5b.A35", "conj(alpha(c))^2 in place of ^4", f35a, 1e-3)
falsifier("t5b.A35", "N(c)^(1-2s) in place of N(c)^(2-2s)", f35n, 1e-3)
finish(t0)
