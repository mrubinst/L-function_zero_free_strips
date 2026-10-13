"""Section 'The completed sums and the second theta derivative' (lines 705-808): claims B34-B52, B54, B56-B61.
(B53 and B55, the reflected weight and the cusp prefactor, are tested end to end in t5_reflection_j2.py.)

The theta coefficients are implemented here directly from the tables of the paper's appendix (lines 1700-1760):
tau(+-lambda^{3nu-4} c b^3) etc., g(r,c) = sum_{x mod c} (x/c)_3 ce(rx/c) with (x/(-2))_3 = x mod 2, and
d_sigma(l) = conj(t_sigma(-l)).  B37 compares this transcription with tests_guide/theta.py (a transcription of
Dunn-Radziwill whose automorphy is checked there, and spot-checked again here).

B34 b = 1 term of T = X^{-1/2} * column sum of eq:energy.      B35 the (n, b) sums are finite.
B36 lambda = i sqrt 3.                                          B37 appendix table = theta.py; automorphy spot check.
B38 eq:theta-coeff (falsifier: (lambda/n)_3 unconjugated).      B39 selection of indices (falsifier: 'x squarefree').
B40 support of d_0, d_+, d_- in u lambda^m n b^3, m >= -4, n squarefree.
B41 d_0(lambda^{-4}) = 0 (falsifiers: d_+(lambda^{-4}), d_-(lambda^{-4}) are nonzero).
B42 (lambda/n)_3 depends only on n mod 9 (falsifier: mod 3).   B43 (lambda/n b^3)_3 = (lambda/n)_3.
B44 d_zbar^2 ce(lz) = (2 pi i conj l)^2 ce(lz) (falsifier: (2 pi i l)^2).
B45 Bessel-Mellin with two derivatives leaves conj(alpha(l))^2 (falsifier: conj(alpha)^1).
B46 conj(alpha(lambda^{-3}))^2 = -1 and the factorization of conj(alpha(l))^2 (falsifier: exponent 1).
B47 the T summand equals -(1/81) d_0(l)(lambda/lambda^3 l)_3 conj(alpha(l))^2 N(l)^{-1/2} Psi(lambda^3 l) V(27 N(l)/X).
B48-B52 chain rule (sympy).  B54 gamma quotient shifts.  B56, B59 Psi = xi_0 prod chi_pi^{j_pi}.
B57 local factors (finite sums).  B58 B_{pi,1} = (x/pi)_2.  B60 Q^2/X < X iff Q < X.  B61 decay of V#.
Run time: about 1 minute.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import *
import theta as TH
TH.factor_primary = lambda n, plist=None: list(factor(n))       # same factorization, faster
from fastgauss import cubic_prime
import sympy as sp, mpmath as mp
from fractions import Fraction as Fr
random.seed(808)
LAMc = cx(LAM)
Wc = cmath.exp(2j*math.pi/3)
Z9 = cmath.exp(2j*math.pi/9)
sq = primaries(3000, sqfree=True)

# ------------------------------------------------------------------ appendix tables (own transcription)
def cub_appendix(x, c):
    """(x/c)_3 for squarefree primary c, possibly even: at -2 the cube root of unity = x mod 2; else exponent (N-1)/3."""
    r = 1 + 0j; m = c
    if divides((2, 0), m):
        m = div(m, (-2, 0))
        if divides((2, 0), x): return 0
        for val, rep in [(1, (1, 0)), (Wc, (0, 1)), (Wc*Wc, (-1, -1))]:
            if divides((2, 0), sub(x, rep)): r *= val; break
    for p in (factor(m) if N(m) > 1 else []):
        r *= cubic_prime(x, p)
    return r
_gc = {}
def g_app(r, c):
    """g(r,c) = sum_{x mod c} (x/c)_3 ce(rx/c).  Even c (or N(c)<=300): brute force with cub_appendix (Euler criterion,
    exponent (N-1)/3, and x mod 2 at -2).  Odd c with N(c)>300: vectorized, (x/c)_3 = chi_c(x)^2 (checked in t3, B22)."""
    key = (r, c)
    if key not in _gc:
        if divides((2, 0), c) or N(c) <= 300:
            cc = cx(c); s = 0
            for x in residues(c):
                k = cub_appendix(x, c)
                if k: s += k*cmath.exp(4j*math.pi*(cx(mul(r, x))/cc).real)
        else:
            A, B = hnf_reps(c)
            s = complex(np.sum(chi_vec(c, A, B, 2)*ce_phase_vec(A, B, c, r)))
        _gc[key] = s
    return _gc[key]
def lam_split(m):
    k = 0
    while divides(LAM, m): m = div(m, LAM); k += 1
    return k, m
PRIM_ALL = [(-2, 0)] + list(PR_ODD)
def ucb(m):
    """m coprime to lambda: (u, c, b) with m = u c b^3, c squarefree primary, b primary (may be even); None if some
    prime exponent is 2 mod 3."""
    rest = m; c = (1, 0); b = (1, 0)
    for p in PRIM_ALL:
        if N(rest) == 1: break
        if N(p) > N(rest): break
        e = 0
        while divides(p, rest): rest = div(rest, p); e += 1
        if e:
            q, r = divmod(e, 3)
            if r == 2: return None
            c = mul(c, pw(p, r)); b = mul(b, pw(p, q))
    assert N(rest) == 1, (m, rest)
    return rest, c, b                                         # rest is the unit u
def tau_app(m):
    """tau(l) at l = lambda^{-3} m (m in O), from the first table of the appendix."""
    if m == (0, 0): return 0
    k, mp_ = lam_split(m); kk = k - 3
    dec = ucb(mp_)
    if dec is None: return 0
    u, c, b = dec; bc = abs(cx(b))/abs(cx(c)); L2 = mul(LAM, LAM)
    if kk % 3 == 2:                                           # kk = 3 nu - 4, nu >= 1
        nu = (kk + 4)//3
        if nu < 1: return 0
        if u in [(1, 0), (-1, 0)]: return 3**(nu/2 + 2)*g_app(L2, c).conjugate()*bc
        if u in [(0, 1), (0, -1)]: return 3**(nu/2 + 2)*Z9**-1*g_app(mul((0, 1), L2), c).conjugate()*bc
        return 3**(nu/2 + 2)*Z9*g_app(mul((-1, -1), L2), c).conjugate()*bc
    if kk % 3 == 0:                                           # kk = 3 nu - 3, nu >= 0, units +-1 only
        nu = (kk + 3)//3
        if nu < 0: return 0
        if u in [(1, 0), (-1, 0)]: return 3**(nu/2 + 2.5)*g_app((1, 0), c).conjugate()*bc
    return 0
def d0_lam3(m):
    """d_0(l) = conj(tau(-l)) at l = lambda^{-3} m."""
    return complex(tau_app(neg(m))).conjugate()
def tau12_app(which, m4):
    """tau_1(omega^j lambda^{-4} c b^3) or tau_2(-omega^j lambda^{-4} c b^3), argument l = lambda^{-4} m4, m4 in O."""
    if m4 == (0, 0): return 0
    k, mp_ = lam_split(m4)
    if k != 0: return 0
    dec = ucb(mp_)
    if dec is None: return 0
    u, c, b = dec; bc = abs(cx(b))/abs(cx(c)); L2 = mul(LAM, LAM)
    gs = {0: g_app(L2, c), 1: g_app(mul((0, 1), L2), c), 2: g_app(mul((-1, -1), L2), c)}
    if which == 1:
        js = {(1, 0): 0, (0, 1): 1, (-1, -1): 2}
        if u not in js: return 0
        j = js[u]; pref = {0: 9*Wc, 1: 9*Z9**-1*Wc*Wc, 2: 9*Z9}[j]
    else:
        js = {(-1, 0): 0, (0, -1): 1, (1, 1): 2}
        if u not in js: return 0
        j = js[u]; pref = {0: 9*Wc*Wc, 1: 9*Z9**-1, 2: 9*Wc*Z9}[j]
    return pref*gs[j].conjugate()*bc
L4c = LAMc**-4
def ce_c(z): return cmath.exp(4j*math.pi*z.real)
def d_pm(sign, m4):
    """d_+(l), d_-(l) at l = lambda^{-4} m4: t_+(l) = omega tau_2(omega l) ce(l), t_-(l) = omega^2 tau_1(omega^2 l) ce(l),
    d_sigma(l) = conj(t_sigma(-l))."""
    l = L4c*cx(m4); ml = neg(m4)
    if sign == '+':
        t = Wc*tau12_app(2, mul((0, 1), ml))*ce_c(-l)
    else:
        t = Wc*Wc*tau12_app(1, mul((-1, -1), ml))*ce_c(-l)
    return complex(t).conjugate()

# ---------------- B36
report("B36", abs(LAMc - 1j*SQ3) + abs(LAMc - (1 + 2*Wc)), 1e-14, "lambda = 1 + 2 omega = i sqrt(3)")

# ---------------- B37: own transcription vs theta.py; automorphy spot check
err = 0; cnt = 0; nz = 0
for a in range(-11, 12):
    for b in range(-11, 12):
        m = (a, b)
        if m == (0, 0): continue
        t1 = tau_app(m); t2 = TH.tau(m)
        err = max(err, abs(t1 - t2)); cnt += 1; nz += abs(t1) > 1e-9
th = TH.Theta(R=4.5); errA = 0
for (z, v) in [(0.13 + 0.21j, 0.93), (-0.21 + 0.05j, 0.97)]:
    a0 = th(z, v); z2, v2 = TH.act(((0, 0), (-1, 0), (1, 0), (0, 0)), z, v)
    errA = max(errA, abs(th(z + 1, v) - a0)/abs(a0), abs(th(z2, v2) - a0)/abs(a0))
report("B37", max(err, errA), 1e-8, f"appendix tau table = theta.py tau on {cnt} indices lambda^-3 m, |m|<=11 ({nz} nonzero, even c included): "
       f"max diff {err:.1e}; automorphy of theta (z->z+1, inversion): {errA:.1e}", note="cited input (DR); spot check only")

# ---------------- B38: d_0(lambda^{-3} n b^3) = 3^{5/2}|b| conj((lambda/n)_3) gamma_2(n)
bs = [(1, 0), (-2, -3), (1, 3), (4, 3), (-5, 0), mul((-2, -3), (-2, -3)), mul((-2, -3), (1, 3))]
err = 0; f38 = 0; cnt = 0
for n in sq[:30]:
    for b in bs:
        x = mul(n, pw(b, 3))
        d0 = d0_lam3(x)
        g2 = gamma_j_fast(2, n)
        rhs = 3**2.5*abs(cx(b))*complex(cubic3(n, LAM)).conjugate()*g2
        err = max(err, abs(d0 - rhs)/abs(rhs)); cnt += 1
        f38 = max(f38, abs(d0 - 3**2.5*abs(cx(b))*complex(cubic3(n, LAM))*g2)/abs(rhs))
report("B38", err, 1e-10, f"{cnt} indices (30 squarefree n, b in 7 values incl. b non-squarefree and inert)")
falsifier("B38", "(lambda/n)_3 not conjugated", f38, 1e-3)

# ---------------- B39: selection
from collections import Counter
bad = 0; cnt = 0; fsel = 0
for x in primaries(1200):
    cn = Counter(factor(x)); form = all(e % 3 != 2 for e in cn.values())
    nonzero = abs(d0_lam3(x)) > 1e-9
    bad += (form != nonzero); cnt += 1
    fsel += (nonzero != (mu(x) != 0))
report("B39", bad, 0, f"{cnt} primary x coprime to 6 (N<=1200): d_0(lambda^-3 x) != 0 iff x = n b^3 with n squarefree")
falsifier("B39", "'iff x squarefree'", fsel, 0)

# ---------------- B40, B41: support of the three coefficient sequences on lambda^{-4} O
bad = 0; cnt = 0; nzs = Counter()
for a in range(-9, 10):
    for b in range(-9, 10):
        m4 = (a, b)
        if m4 == (0, 0): continue
        vals = {'0': (d0_lam3(div(m4, LAM)) if divides(LAM, m4) else 0), '+': d_pm('+', m4), '-': d_pm('-', m4)}
        k, mp_ = lam_split(m4)                               # l = lambda^{k-4} * (unit c b^3 part)
        for s_, v in vals.items():
            if abs(v) > 1e-9:
                nzs[s_] += 1; cnt += 1
                bad += (k - 4 < -4) or (ucb(mp_) is None)
report("B40", bad, 0, f"{cnt} nonzero coefficients d_sigma(l), l in lambda^-4 O, |m|<=9 ({dict(nzs)}): all of the form u lambda^m n b^3, m>=-4, n squarefree")
d0_l4 = 0.0 if not divides(LAM, (1, 0)) else d0_lam3(div((1, 0), LAM))   # l = lambda^{-4} * 1: m4 = 1 is not divisible by lambda
# direct: tau at -lambda^{-4}: exponent -4 = 3nu-4 with nu = 0 < 1, and -4 != 3nu-3: zero in the table
k_exp = -4; row1 = (k_exp % 3 == 2 and (k_exp + 4)//3 >= 1); row2 = (k_exp % 3 == 0 and (k_exp + 3)//3 >= 0)
report("B41", abs(d0_l4) + (1 if (row1 or row2) else 0), 0, "d_0(lambda^-4) = conj(tau(-lambda^-4)) = 0: lambda^-4 is not in the support of the tau table "
       f"(rows need nu>=1 or exponent 3nu-3); d_+(lambda^-4) = {abs(d_pm('+', (1, 0))):.3g}, d_-(lambda^-4) = {abs(d_pm('-', (1, 0))):.3g}")
falsifier("B41", "the same claim for d_+ (d_+(lambda^-4) = 0)", abs(d_pm('+', (1, 0))), 1e-6)

# ---------------- B42, B43
tab9 = {}; tab3 = {}; bad = 0; bad3 = 0
for n in primaries(3000):
    v = complex(cubic3(n, LAM))
    k9 = (n[0] % 9, n[1] % 9); k3 = (n[0] % 3, n[1] % 3)
    if k9 in tab9 and abs(tab9[k9] - v) > 1e-9: bad += 1
    if k3 in tab3 and abs(tab3[k3] - v) > 1e-9: bad3 += 1
    tab9.setdefault(k9, v); tab3.setdefault(k3, v)
report("B42", bad, 0, f"all primary n coprime to 6 with N(n)<=3000: {len(tab9)} classes mod 9")
falsifier("B42", "depends only on n mod 3", bad3, 0)
err = max(abs(cubic3(mul(n, pw(b, 3)), LAM) - cubic3(n, LAM)) for n in sq[:60] for b in bs[1:5])
report("B43", err, 1e-12, "(lambda/n b^3)_3 = (lambda/n)_3 for 60 n and 4 b; so (lambda/lambda^3 l)_3 d_0(l) = 3^{5/2}|b|gamma_2(n)")

# ---------------- B44, B48-B52: symbolic
z, zb, v, c, cb, dl, dlb, lv, lb = sp.symbols('z zb v c cb delta deltab l lb')
ceexpr = sp.exp(2*sp.pi*sp.I*(lv*z + lb*zb))
ok44 = sp.simplify(sp.diff(ceexpr, zb, 2) - (2*sp.pi*sp.I*lb)**2*ceexpr) == 0
f44 = sp.simplify(sp.diff(ceexpr, zb, 2) - (2*sp.pi*sp.I*lv)**2*ceexpr) != 0
report("B44", 0 if ok44 else 1, 0, "sympy: d^2/dzbar^2 exp(2 pi i (l z + conj(l) zbar)) = (2 pi i conj(l))^2 exp(...)")
falsifier("B44", "(2 pi i l)^2", 1 if f44 else 0, 0)
# B48: (z', v') = g^{-1}(a/c + z, v), g = (a b; c delta) in SL_2, with independent conjugate symbols
a_, b_, ab_, bb_ = sp.symbols('a b ab bb')
def act(A, B, C, Dm, Ab, Bb, Cb, Db, Z, Zb, V):
    den = (C*Z + Dm)*(Cb*Zb + Db) + C*Cb*V**2
    return ((A*Z + B)*(Cb*Zb + Db) + A*Cb*V**2)/den, ((Ab*Zb + Bb)*(C*Z + Dm) + Ab*C*V**2)/den, V/den
Zp, Zpb, Vp = act(dl, -b_, -c, a_, dlb, -bb_, -cb, ab_, a_/c + z, ab_/cb + zb, v)
sub_det = {b_: (a_*dl - 1)/c, bb_: (ab_*dlb - 1)/cb}
claimZ = -dl/c - zb/(c**2*(v**2 + z*zb)); claimV = v/(c*cb*(v**2 + z*zb)); claimZb = -dlb/cb - z/(cb**2*(v**2 + z*zb))
ok48 = all(sp.simplify(sp.together(e).subs(sub_det)) == 0 for e in [Zp - claimZ, Vp - claimV, Zpb - claimZb])
Zg, _, Vg = act(a_, b_, c, dl, ab_, bb_, cb, dlb, a_/c + z, ab_/cb + zb, v)     # falsifier: g instead of g^{-1}
f48 = sp.simplify(sp.together(Zg - claimZ).subs(sub_det)) != 0
# numeric spot check with an SL_2(Z[omega]) matrix and theta.act
A0, C0 = (3, 1), (2, 0)          # coprime (norms 7 and 4); find B0, Dl with A0 Dl - B0 C0 = 1
Dl = None
for B0 in [(x, y) for x in range(-8, 9) for y in range(-8, 9)]:
    t = add((1, 0), mul(B0, C0))
    if divides(A0, t):
        Dl = div(t, A0); break
assert Dl is not None and sub(mul(A0, Dl), mul(B0, C0)) == (1, 0)
num_err = 0
for zz, vv in [(0.17 - 0.08j, 0.61), (-0.3 + 0.2j, 1.3)]:
    zp, vp = TH.act((cx(Dl), -cx(B0), -cx(C0), cx(A0)), cx(A0)/cx(C0) + zz, vv)       # g^{-1} = (delta -b; -c a)
    num_err = max(num_err, abs(zp - (-cx(Dl)/cx(C0) - zz.conjugate()/(cx(C0)**2*(vv**2 + abs(zz)**2)))),
                  abs(vp - vv/(N(C0)*(vv**2 + abs(zz)**2))))
report("B48", (0 if ok48 else 1) + num_err, 1e-12, f"sympy (ad-bc=1) and a numeric SL_2(Z[omega]) example (error {num_err:.1e})")
falsifier("B48", "g in place of g^{-1}", 1 if f48 else 0, 0)
# B49
d_list = []
bad49 = 0
for k in range(1, 6):
    dz = sp.simplify(sp.diff(claimZ, zb, k).subs(z, 0).subs(zb, 0))
    dzb = sp.simplify(sp.diff(claimZb, zb, k).subs(z, 0).subs(zb, 0))
    dv = sp.simplify(sp.diff(claimV, zb, k).subs(z, 0).subs(zb, 0))
    bad49 += (dzb != 0) + (dv != 0)
    bad49 += (k == 1 and sp.simplify(dz + 1/(c*v)**2) != 0) + (k >= 2 and dz != 0)
report("B49", bad49, 0, "sympy, orders 1..5: d^k zbar'=d^k v'=0, d z'=-(cv)^-2, d^k z'=0 (k>=2) at z=zbar=0")
# B50, B51, B52 with a generic polynomial F(X, Y, V) (X = z', Y = zbar', V = v')
X_, Y_, V_ = sp.symbols('X Y V'); co = sp.symbols('a0:45')
Fp = sum(co[i]*X_**(i % 5)*Y_**((i//5) % 3)*V_**(i//15) for i in range(45))
expr = Fp.subs({X_: claimZ, Y_: claimZb, V_: claimV})
lhs = sp.diff(expr, zb, 2).subs(z, 0).subs(zb, 0)
pt = {X_: -dl/c, Y_: -dlb/cb, V_: 1/(c*cb*v)}
rhs = (c*v)**-4*sp.diff(Fp, X_, 2).subs(pt)
ok50 = sp.simplify(lhs - rhs) == 0
f50a = sp.simplify(lhs - (c*v)**-4*sp.diff(Fp, Y_, 2).subs(pt)) != 0
f50b = sp.simplify(lhs - (cb*v)**-4*sp.diff(Fp, X_, 2).subs(pt)) != 0
report("B50", 0 if ok50 else 1, 0, "sympy, generic polynomial F of degree (4,2,2) with symbolic coefficients")
falsifier("B50", "d_{zbar'}^2 F in place of d_{z'}^2 F", 1 if f50a else 0, 0)
falsifier("B50", "(conj(c) v)^-4 in place of (cv)^-4", 1 if f50b else 0, 0)
# B51: the Faa di Bruno coefficients of all other derivatives of F vanish
d1 = {s_: sp.diff(e, zb).subs(z, 0).subs(zb, 0) for s_, e in [('X', claimZ), ('Y', claimZb), ('V', claimV)]}
d2 = {s_: sp.diff(e, zb, 2).subs(z, 0).subs(zb, 0) for s_, e in [('X', claimZ), ('Y', claimZb), ('V', claimV)]}
coeffs = {'F_X': d2['X'], 'F_Y': d2['Y'], 'F_V': d2['V'], 'F_XY': 2*d1['X']*d1['Y'], 'F_XV': 2*d1['X']*d1['V'],
          'F_YY': d1['Y']**2, 'F_YV': 2*d1['Y']*d1['V'], 'F_VV': d1['V']**2}
nonzero = [k for k, e in coeffs.items() if sp.simplify(e) != 0]
report("B51", len(nonzero), 0, "coefficients of F_X, F_Y, F_V, F_XY, F_XV, F_YY, F_YV, F_VV in d_zbar^2 F(z',zbar',v') at 0 all vanish; "
       f"F_XX coefficient = {sp.simplify(d1['X']**2)}")
phi = sum(co[i]*V_**i for i in range(8))
ok52 = sp.simplify(sp.diff(phi.subs(V_, claimV), zb, 2).subs(z, 0).subs(zb, 0)) == 0
ok52b = sp.simplify(sp.diff(phi.subs(V_, v), zb, 2)) == 0
report("B52", 0 if ok52 and ok52b else 1, 0, "d_zbar^2 of a function of v' (resp. v) alone vanishes at z=0 (sympy, generic degree-7 phi)",
       note="'only the identity cusp has a constant term' is cited (DR (5.9),(5.15)); not tested")

# ---------------- B45: Bessel-Mellin with two derivatives
mp.mp.dps = 25; third = mp.mpf(1)/3
worst = 0; f45 = 0
for l in [mp.mpc(0.7, -1.3), mp.mpc(-0.4, 0.9), mp.mpc(2.1, 0.3)]:
    for s in [mp.mpc(1.4, 0.8), mp.mpc(0.9, -2.2), mp.mpf(2.3)]:
        al = l/abs(l)
        num = mp.quad(lambda vv: (2j*mp.pi*mp.conj(l))**2*vv*mp.besselk(third, 4*mp.pi*abs(l)*vv)*vv**(2*s), [0, 0.5, 2, mp.inf])
        cl = -mp.mpf(1)/4*(2*mp.pi)**(-2*s)*mp.gamma(s + mp.mpf(5)/6)*mp.gamma(s + mp.mpf(7)/6)*mp.conj(al)**2*abs(l)**(-2*s)
        worst = max(worst, abs(num - cl)/abs(cl))
        f45 = max(f45, abs(num - cl/mp.conj(al))/abs(cl))
report("B45", float(worst), 1e-15, "int_0^oo (2 pi i conj l)^2 v K_1/3(4 pi|l|v) v^{2s} dv = -(1/4)(2pi)^{-2s} Gamma(s+5/6)Gamma(s+7/6) conj(alpha(l))^2 |l|^{-2s}, 9 cases")
falsifier("B45", "conj(alpha(l))^1", float(f45), 1e-3)

# ---------------- B46
aL3 = (LAMc**-3)/abs(LAMc**-3)
e46 = abs(aL3.conjugate()**2 + 1); fac = 0
for n in sq[:40]:
    for b in bs:
        l = LAMc**-3*cx(n)*cx(b)**3
        fac = max(fac, abs((l/abs(l)).conjugate()**2 - aL3.conjugate()**2*alpha_c(n).conjugate()**2*alpha_c(b).conjugate()**6))
report("B46", e46 + fac, 1e-12, f"conj(alpha(lambda^-3))^2 = {complex(aL3.conjugate()**2):.3g}; factorization over 280 (n, b)")
falsifier("B46", "zeta exponent: conj(alpha(lambda^-3)) = -1", abs(aL3.conjugate() + 1), 1e-3)

# ---------------- B34, B35, B47: the completed sum T and its theta form
def Wb(y): return math.exp(-1/((y - 1)*(2 - y))) if 1 < y < 2 else 0.0
Xs = 1400.0; kk, ff = (5, 3), (7, 9)                         # cols complete: 2X < 3000
def xi(n): return complex(cubic3(n, (0, 1)))                   # (omega/n)_3: a character of n mod 9
def Psi(n): return xi(n)*chi(n, kk)*chi(n, ff)**4
cols = [n for n in sq if N(n) < 2*Xs]
bsall = primaries(int((2*Xs)**(1/3)) + 2)                      # all primary b coprime to 6 with N(b)^3 < 2X
terms = []; bad35 = 0
for b in bsall:
    for n in cols:
        y = N(n)*N(b)**3/Xs
        if Wb(y) == 0: continue
        t = alpha_c(n).conjugate()**2*gamma_j_fast(2, n)*Psi(n)*alpha_c(b).conjugate()**6*Psi(b)**3*math.sqrt(N(b))*Wb(y)/math.sqrt(Xs)
        terms.append((n, b, t))
bad35 = sum(1 for b in primaries(60) if N(b)**3 >= 2*Xs and any(Wb(N(n)*N(b)**3/Xs) for n in cols[:5]))
colsum = sum(alpha_c(n).conjugate()**2*gamma_j_fast(2, n)*xi(n)*chi(n, kk)*chi(n, ff)**4*Wb(N(n)/Xs) for n in cols)
b1 = sum(t for n, b, t in terms if b == (1, 0))
report("B34", abs(b1 - colsum/math.sqrt(Xs)), 1e-12, f"X={Xs:g}, k={kk}, f={ff}, xi=(omega/n)_3: b=1 part {complex(b1):.6g}")
report("B35", bad35, 0, f"{len(terms)} nonzero (n,b) terms; none with N(b)^3 >= 2X")
# termwise identity on a grid of (n, b) (it does not need the support of W): use W(y) = y e^{-y}, V(y) = sqrt(y) W(y)
def Wg(y): return y*math.exp(-y)
X47 = 60.0; err47 = 0; f47 = 0; cnt47 = 0
for n in [(1, 0)] + sq[1:40]:
    for b in bs:
        x = mul(n, pw(b, 3)); y = N(n)*N(b)**3/X47
        t = alpha_c(n).conjugate()**2*gamma_j_fast(2, n)*Psi(n)*alpha_c(b).conjugate()**6*Psi(b)**3*math.sqrt(N(b))*Wg(y)/math.sqrt(X47)
        l = LAMc**-3*cx(n)*cx(b)**3; Nl = abs(l)**2
        V = math.sqrt(27*Nl/X47)*Wg(27*Nl/X47)
        pred = -(1/81)*d0_lam3(x)*complex(cubic3(x, LAM))*(l/abs(l)).conjugate()**2*Nl**-0.5*Psi(x)*V
        sc = max(abs(t), 1e-300)
        err47 = max(err47, abs(t - pred)/sc if abs(t) > 1e-14 else abs(pred)); f47 = max(f47, abs(t + pred)/sc if abs(t) > 1e-14 else 0)
        cnt47 += 1
nzero = sum(1 for n in [(1, 0)] + sq[1:40] for b in bs if abs(Psi(mul(n, pw(b, 3)))) < 1e-12)
report("B47", err47, 1e-10, f"{cnt47} pairs (n,b) (40 n incl. n=1, 7 b incl. non-squarefree and inert; {nzero} with Psi=0), relative error; "
       "the constant is -1/81 = -1/(3^(5/2) 27^(1/2)), the sign from conj(alpha(lambda^-3))^2 = -1")
falsifier("B47", "+1/81", f47, 1e-3)

# ---------------- B54
ok54 = Fr(7, 6) + Fr(1, 2) == Fr(5, 3) and Fr(5, 6) + Fr(1, 2) == Fr(4, 3)
s_, t_ = sp.symbols('s t')
kern = sp.gamma(1 - s_ + sp.Rational(5, 6))*sp.gamma(1 - s_ + sp.Rational(7, 6))/(sp.gamma(s_ + sp.Rational(5, 6))*sp.gamma(s_ + sp.Rational(7, 6)))
claim = sp.gamma(sp.Rational(4, 3) + t_)*sp.gamma(sp.Rational(5, 3) + t_)/(sp.gamma(sp.Rational(4, 3) - t_)*sp.gamma(sp.Rational(5, 3) - t_))
ok54b = sp.simplify(kern.subs(s_, sp.Rational(1, 2) - t_) - claim) == 0
report("B54", 0 if ok54 and ok54b else 1, 0, "7/6+1/2=5/3, 5/6+1/2=4/3; the j=2 Bessel-Mellin quotient Gamma(11/6-s)Gamma(13/6-s)/(Gamma(s+5/6)Gamma(s+7/6)) at s=1/2-t is the kernel of eq:weight")

# ---------------- B56, B59: Psi = xi_0(n) prod chi_pi(n)^{j_pi}, xi_0 a class function mod 72
def vals_primes(k):
    out = {}
    for p in PR_ODD:
        if N(p) > N(k): break
        e = 0; m = k
        while divides(p, m): m = div(m, p); e += 1
        if e: out[p] = e
    return out
P7 = [p for p in PR_ODD if N(p) == 7][0]; P13 = [p for p in PR_ODD if N(p) == 13][0]; P19 = [p for p in PR_ODD if N(p) == 19][0]
P31 = [p for p in PR_ODD if N(p) == 31][0]
tests = [(mul(mul((2, 0), P7), mul(P13, P13)), (1, 0)), (mul(mul(LAM, (0, 1)), P19), P7), (mul((1, 1), (-5, 0)), mul(P31, P13)),
         (mul(P7, P7), P7), (mul(pw(P13, 3), (-1, 0)), P19), (P31, (1, 0))]
bad = 0; nclass = []; bad59 = 0
for k, f in tests:
    vk = vals_primes(k); vf = vals_primes(f)
    Pset = sorted(set(vk) | set(vf))
    jp = {p: (vk.get(p, 0) + 4*vf.get(p, 0)) % 6 for p in Pset}
    for p in Pset:
        if vk.get(p, 0) == 1 and vf.get(p, 0) == 0: bad59 += jp[p] != 1
    tab = {}
    for n in primaries(2500):
        Ps = xi(n)*chi(n, k)*chi(n, f)**4
        prod = 1 + 0j
        for p in Pset:
            cp = chi(p, n)
            prod *= (cp**jp[p] if cp != 0 else 0)
        if (abs(Ps) > 1e-9) != (abs(prod) > 1e-9): bad += 1; continue
        if abs(Ps) < 1e-9: continue
        r = Ps/prod; key = (n[0] % 72, n[1] % 72)
        if key in tab and abs(tab[key] - r) > 1e-9: bad += 1
        tab.setdefault(key, r)
    nclass.append(len(tab))
report("B56", bad, 0, f"6 pairs (k,f) incl. units, lambda, 2, squares and cubes in k: Psi(n)/prod chi_pi(n)^j_pi is a function of n mod 72 "
       f"(classes hit: {nclass}), and the zero sets agree (all primary n, N<=2500)")
report("B59", bad59, 0, "every prime pi || k with pi not dividing 6f has j_pi = 1")

# ---------------- B57, B58: local factors as finite sums
err57 = 0; f57 = 0; ncase = 0; ratios57 = {}
for p in [q for q in PR_ODD if 7 <= N(q) <= 100][:7]:
    res = residues(p); Np = N(p); pc = cx(p)
    inv = {h: next(y for y in res if canon_mod(sub(mul(y, h), (1, 0)), p) == (0, 0)) for h in res if chi(p, h) != 0}
    def gs(k): return sum(chi(p, y)**k*e_add(cx(y)/pc) for y in inv)
    for eps in [y for y in list(inv)[:3]]:
        for j in range(6):
            for x in res:
                cxv = chi(p, x)
                if j == 0:     # active branch of 1_{p nmid x}: coefficient -1/N, Kubota factor chi(h)^-2
                    S = -sum(chi(p, h)**-2*e_add(cx(mul(mul(eps, inv[h]), x))/pc) for h in inv)/Np
                    B = Np**-0.5*(cxv**-2 if cxv != 0 else 0)
                else:          # chi^j = g(chi^-j)^{-1} sum_h chi(h)^-j e(hx/p); reflected with chi(h)^-2 and e(eps h^-1 x/p)
                    S = sum(chi(p, h)**(-j - 2)*e_add(cx(mul(mul(eps, inv[h]), x))/pc) for h in inv)/gs(-j)
                    if j == 4: B = Np**-0.5*(-1 + Np*(1 if cxv == 0 else 0))
                    else: B = (cxv**(-j - 2) if cxv != 0 else 0)
                ncase += 1
                if abs(B) < 1e-12:
                    err57 = max(err57, abs(S))
                else:
                    ratio = S/B
                    key = (p, eps, j)
                    if key in ratios57: err57 = max(err57, abs(ratio - ratios57[key]))
                    ratios57.setdefault(key, ratio)
                    err57 = max(err57, abs(abs(ratio) - 1))
                    if j not in (0, 4) and cxv != 0:
                        f57 = max(f57, abs(abs(S/(cxv**(-j - 1))) - 1) + abs(S/(cxv**(-j - 1)) - S/(chi(p, (1, 0))**(-j - 1))))
report("B57", err57, 1e-10, f"{ncase} (prime, eps, j, x) cases, 7 primes: the reflected local sums equal a phase (|.|=1, independent of x) times B_pi,j(x) "
       "for every j in 0..5, including the N^-1/2 normalizations for j=0 (active) and j=4",
       note="the extra chi(h)^-2 (Kubota character) is the cited theta input; the finite-sum algebra is what is tested")
# falsifier: exponent -j-1 instead of -j-2 (ratio is not constant in x)
fal = 0
p = [q for q in PR_ODD if N(q) == 13][0]; res = residues(p); pc = cx(p)
inv = {h: next(y for y in res if canon_mod(sub(mul(y, h), (1, 0)), p) == (0, 0)) for h in res if chi(p, h) != 0}
gsm = sum(chi(p, y)**-1*e_add(cx(y)/pc) for y in inv)
rat = [sum(chi(p, h)**-3*e_add(cx(mul(inv[h], x))/pc) for h in inv)/gsm/(chi(p, x)**-2) for x in list(inv)[:6]]
falsifier("B57", "B_{pi,1} = chi^(-j-1) = chi^-2 (ratio not constant in x)", max(abs(r_ - rat[0]) for r_ in rat), 1e-3)
err58 = 0
for p in [q for q in PR_ODD if N(q) < 200][:10]:
    for x in residues(p)[:40]:
        cv = chi(p, x)
        if cv == 0: continue
        q2 = 1 if powmod(canon_mod(x, p), (N(p) - 1)//2, p) == (1, 0) else -1
        err58 = max(err58, abs(cv**-3 - cv**3), abs(cv**3 - q2))
report("B58", err58, 1e-12, "10 primes x 40 residues: chi^-3 = chi^3 = quadratic symbol by Euler's criterion")

# ---------------- B60
viol = sum(1 for Q in range(1, 120) for X in range(1, 200) if (Fr(Q*Q, X) < X) != (Q < X))
report("B60", viol, 0, "Q^2/X < X iff Q < X for all integers 1<=Q<120, 1<=X<200")
falsifier("B60", "'iff Q^2 < X'", sum(1 for Q in range(1, 60) for X in range(1, 100) if (Fr(Q*Q, X) < X) != (Q*Q < X)), 0)

# ---------------- B61: decay of V#(x), V(y) = exp(-(log y)^2/(2 s^2)) (V^ Gaussian), by contour shifts
sw = 0.5
def R2(t): return mp.gamma(mp.mpf(4)/3 + t)*mp.gamma(mp.mpf(5)/3 + t)/(mp.gamma(mp.mpf(4)/3 - t)*mp.gamma(mp.mpf(5)/3 - t))
def Vhat(sv): return mp.sqrt(2*mp.pi)*sw*mp.e**(sw*sw*sv*sv/2)
def Vsharp(x, c):
    K = (2*mp.pi)**4*x/27
    f = lambda u: Vhat(-(c + 1j*u))*R2(c + 1j*u)*K**(-(c + 1j*u))
    return mp.quad(f, [-40, -10, 0, 10, 40])/(2*mp.pi)
mp.mp.dps = 20
shift_err = max(abs(Vsharp(x, 0) - Vsharp(x, 1.5))/max(abs(Vsharp(x, 0)), 1e-30) for x in [0.05, 0.5, 2.0])
def Cbound(c):
    f = lambda u: abs(Vhat(-(c + 1j*u))*R2(c + 1j*u))
    return mp.quad(f, [-40, -10, 0, 10, 40])/(2*mp.pi)
viol = 0; rows = []
for c in [1.0, 3.0, 5.0]:
    Cc = Cbound(c)
    for x in [1, 10, 100, 1000, 10000]:
        K = (2*mp.pi)**4*x/27; v0 = abs(Vsharp(x, 0)) if x <= 100 else abs(Vsharp(x, c))
        viol += v0 > Cc*K**(-c)*(1 + 1e-12)
        if c == 3.0: rows.append(f"x={x}: |V#|={float(v0):.2e} <= {float(Cc*K**(-c)):.2e}")
report("B61", float(shift_err) + viol, 1e-8,
       f"V(y)=exp(-(log y)^2/0.5): contour shift Re t=0 -> 1.5 agrees ({float(shift_err):.1e}); |V#(x)| <= C_c (Kx)^-c (K=(2pi)^4/27) for c=1,3,5 "
       f"at x=1..1e4; c=3: " + "; ".join(rows),
       note="class A: rapid decay is qualitative; tested via the contour-shift bound, finite for every c")
