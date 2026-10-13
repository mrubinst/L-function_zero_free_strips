"""t03: the reflection table for the CM target (guide L1213-1255), against OAI78 Prop. prop:completed-reflection
(L1677-1880) and its proof (L2281-2420).

Conventions as in OAI78: check-e(z) = exp(2 pi i (z + zbar)), d/dzbar check-e(mu z) = 2 pi i mubar check-e(mu z),
theta expansion terms d(mu) v K_{1/3}(4 pi |mu| v) check-e(mu z), K = (2 pi)^4/27, lambda = 1 + 2 omega = i sqrt 3.

Claims tested
  C28 L1222-1223 direct Mellin factor: int_0^oo dbar^j[v K_{1/3}(4pi|mu|v) e(mu z)] v^{2s-1} dv equals
                 (i/4)(2pi)^{-2s} G(s+1/3)G(s+2/3) abar(mu) q^{-s} (j=1, OAI78 L2299) and
                 -(pi/2)(2pi)^{-2s} G(s+1/3)G(s+2/3) abar(mu)^2 q^{1/2-s} (j=2), i.e. Gamma(w+5/6)Gamma(w+7/6), w = s-1/2.
  C29 L1251-1252 "relative error 3e-26" of the Bessel-Mellin identities (reproduced at 30 digits).
  C30 L1215, L2319-2331 chain rule: dbar^j of theta-bar(g^{-1}(z+a/c,v)) at z=0 is (-(cv)^{-2})^j d_w^j theta-bar(...).
  C31 L1224-1226 reflected kernel and cusp scalar: the full reflected/direct ratio equals
                 -(i/81) abar(c)^2 q_mu^{-1/2} (K q_mu/q_c^2)^{-t} R_1(t) alpha(mu)       (OAI78, j=1, calibration)
                 -(1/81) abar(c)^4 q_mu^{-1/2} (K q_mu/q_c^2)^{-t} R_CM(t) alpha(mu)^2   (CM, j=2),
                 R_CM(t) = G(4/3+t)G(5/3+t)/(G(4/3-t)G(5/3-t)), with the direct constants 3^{5/2}/4 and 3 pi/2.
  C32 L1228      dual-side angular factor alpha(mu)^2 (part of C31).
  C33h L1233     holomorphy: CM kernel holomorphic for Re t > -4/3, OAI78's for Re t > -5/6 (exact pole sets).
  C34 L1232      |R(sigma+iu)| << (1+|u|)^{4 sigma} on -1/4 <= sigma <= A, both kernels (numerics + Stirling exponent).
  C35 L1233-1234 the reflected weight V#(x) << min(x^{1/4}, x^{-A}) for both kernels (Gaussian test, numerics).
  C36 L1254      the Stirling exponents of the kernels (4 sigma) and of the direct/reflected gamma products.
Run time: about 60 s.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmreport import claim, cited, info, untested
import mpmath as mp
import sympy as sp
from fractions import Fraction as Fr

mp.mp.dps = 30
third = mp.mpf(1)/3
w3 = mp.exp(2j*mp.pi/3)

def bm_mode(j, mu, s):
    f = lambda v: (2j*mp.pi*mp.conj(mu))**j*v*mp.besselk(third, 4*mp.pi*abs(mu)*v)*v**(2*s - 1)
    return mp.quad(f, [0, 0.05, 0.3, 1, mp.inf])

# ---------------------------------------------------------------- C28, C29
errs = []; ferr = []
for mu, s in [(mp.mpc(0.7, -1.3), mp.mpc(1.4, 0.8)), (mp.mpc(-2.1, 0.4), mp.mpc(0.9, -2.5)), (mp.mpc(0.3, 0.2), mp.mpc(2.2, 0.0))]:
    al = mu/abs(mu); q = abs(mu)**2
    G = mp.gamma(s + third)*mp.gamma(s + 2*third)
    oai = (1j/4)*(2*mp.pi)**(-2*s)*G*mp.conj(al)*q**(-s)
    cm = -(mp.pi/2)*(2*mp.pi)**(-2*s)*G*mp.conj(al)**2*q**(mp.mpf(1)/2 - s)
    w = s - mp.mpf(1)/2
    cm_w = -(mp.pi/2)*(2*mp.pi)**(-2*s)*mp.gamma(w + mp.mpf(5)/6)*mp.gamma(w + mp.mpf(7)/6)*mp.conj(al)**2*q**(-w)
    m1 = bm_mode(1, mu, s); m2 = bm_mode(2, mu, s)
    errs += [abs(m1 - oai)/abs(oai), abs(m2 - cm)/abs(cm), abs(m2 - cm_w)/abs(cm_w)]
    # falsifier: Gamma(s+5/6)Gamma(s+7/6) without the half shift
    fals = -(mp.pi/2)*(2*mp.pi)**(-2*s)*mp.gamma(s + mp.mpf(5)/6)*mp.gamma(s + mp.mpf(7)/6)*mp.conj(al)**2*q**(-w)
    ferr.append(abs(m2 - fals)/abs(fals))
E = max(errs)
claim("C28", 1223, "table: direct Mellin factor", "E",
      "Bessel-Mellin for one and two zbar-derivatives; CM factor Gamma(w+5/6)Gamma(w+7/6) in OAI78's variable w=s-1/2",
      E < 1e-20, f"max rel err {mp.nstr(E, 3)} (3 points); falsifier unshifted Gamma(s+5/6)Gamma(s+7/6): rel err {mp.nstr(min(ferr), 3)}",
      falsifier_failed=min(ferr) > 1e-3)
claim("C29", 1251, "rel. err 3e-26", "E", "Bessel-Mellin identities agree to ~1e-26 or better at 30 digits", E < 1e-25,
      f"max rel err {mp.nstr(E, 3)}; falsifier tolerance 1e-40 not met", falsifier_failed=not (E < 1e-40))

# ---------------------------------------------------------------- C30 (sympy, Wirtinger derivatives)
x, y = sp.symbols('x y', real=True)
c = sp.Rational(2) + 3*sp.Rational(-1, 2) + 3*sp.sqrt(3)*sp.I/2          # c = 2 + 3 omega
dp = 1 - (sp.Rational(-1, 2) + sp.sqrt(3)*sp.I/2)                      # delta' = 1 - omega
v = sp.Rational(7, 10)
qc = sp.expand(c*sp.conjugate(c))
p1, p2, p3, p4 = sp.Rational(3, 10) + sp.I/5, sp.Rational(-1, 4) + sp.I/3, sp.Rational(1, 7), sp.Rational(1, 9) - sp.I/11
z = x + sp.I*y; zb = x - sp.I*y
den = v**2 + z*zb
f = -dp/c - zb/(c**2*den)
fb = -sp.conjugate(dp)/sp.conjugate(c) - z/(sp.conjugate(c)**2*den)
h = v/(qc*den)
Fcomp = sp.exp(p1*f + p2*fb + p3*h + p4*f*fb)
dzb = lambda e: (sp.diff(e, x) + sp.I*sp.diff(e, y))/2
D1 = dzb(Fcomp); D2 = dzb(D1)
at0 = {x: 0, y: 0}
w0 = -dp/c; wb0 = sp.conjugate(w0); h0 = 1/(qc*v)
F0 = sp.exp(p1*w0 + p2*wb0 + p3*h0 + p4*w0*wb0)
Fw = (p1 + p4*wb0)*F0; Fww = (p1 + p4*wb0)**2*F0
lhs1 = complex(sp.N(D1.subs(at0), 30)); rhs1 = complex(sp.N(-(c*v)**-2*Fw, 30))
lhs2 = complex(sp.N(D2.subs(at0), 30)); rhs2 = complex(sp.N((c*v)**-4*Fww, 30))
fal2 = complex(sp.N(-(c*v)**-4*Fww, 30))
e30 = max(abs(lhs1 - rhs1)/abs(rhs1), abs(lhs2 - rhs2)/abs(rhs2))
f30 = abs(lhs2 - fal2)/abs(fal2)
claim("C30", 1215, "two derivatives through the cusp map", "E",
      "dbar of the composite at z=0 is -(cv)^-2 F_w (OAI78 L2324) and dbar^2 is +(cv)^-4 F_ww (generic non-holomorphic F)",
      e30 < 1e-20, f"rel err {e30:.1e}; falsifier -(cv)^-4: rel err {f30:.2f}", falsifier_failed=f30 > 0.1)

# ---------------------------------------------------------------- C31, C32
K = (2*mp.pi)**4/27
lam = 1 + 2*w3
ab_lam3 = mp.conj((lam**-3)/abs(lam**-3))
D1c = (1j/4)*ab_lam3*mp.mpf(3)**2.5          # coefficient of 27^s (2pi)^{-2s}
D2c = -(mp.pi/2)*ab_lam3**2*mp.mpf(27)**-0.5*mp.mpf(3)**2.5
calib1 = abs(D1c - mp.mpf(3)**2.5/4)
calib2 = abs(D2c - 3*mp.pi/2)

def Rker(t, a1, a2):
    return mp.gamma(a1 + t)*mp.gamma(a2 + t)/(mp.gamma(a1 - t)*mp.gamma(a2 - t))

def reflected(j, s, cc, mu):
    """Numerical reflected Mellin integral for one cusp term (unit coefficient d_H(mu)=1, kappa = 1, phase 1):
       int_0^oo (-(c v)^-2)^j (2 pi i mu)^j v' K_{1/3}(4 pi |mu| v') v^{2s-1} dv,  v' = 1/(q_c v)."""
    qc = abs(cc)**2
    f = lambda vv: ((-(cc*vv)**-2)**j*(2j*mp.pi*mu)**j*(1/(qc*vv))*mp.besselk(third, 4*mp.pi*abs(mu)/(qc*vv))
                    *vv**(2*s - 1))
    return mp.quad(f, [0, 0.01, 0.1, 1, 10, mp.inf])

e31 = []; f31 = []
cc = mp.mpc(2, 0) + 3*w3          # c = 2 + 3 omega
for mu in [mp.mpc(0.6, 0.9), mp.mpc(-1.7, 0.3)]:
    qmu = abs(mu)**2; amu = mu/abs(mu); qc = abs(cc)**2; ac = cc/abs(cc)
    for s in [mp.mpc(-0.6, 0.7), mp.mpc(-1.1, -1.9)]:
        # OAI78, j = 1:  T(s) = J(s)/(D1 (27/(2pi)^2)^s Gamma Gamma), t = 1/2 - s
        r1 = reflected(1, s, cc, mu)/(D1c*(27/(2*mp.pi)**2)**s*mp.gamma(s + third)*mp.gamma(s + 2*third))
        t = mp.mpf(1)/2 - s
        pred1 = -(1j/81)*mp.conj(ac)**2*qmu**-0.5*(K*qmu/qc**2)**(-t)*Rker(t, mp.mpf(7)/6, mp.mpf(5)/6)*amu
        # CM, j = 2:  J_2(s) = D2 (27/(2pi)^2)^s Gamma(s+1/3)Gamma(s+2/3) T_2(s - 1/2);  t = 1/2 - (s - 1/2) = 1 - s
        r2 = reflected(2, s, cc, mu)/(D2c*(27/(2*mp.pi)**2)**s*mp.gamma(s + third)*mp.gamma(s + 2*third))
        t2 = 1 - s
        pred2 = -(mp.mpf(1)/81)*mp.conj(ac)**4*qmu**-0.5*(K*qmu/qc**2)**(-t2)*Rker(t2, mp.mpf(4)/3, mp.mpf(5)/3)*amu**2
        e31 += [abs(r1/pred1 - 1), abs(r2/pred2 - 1)]
        # falsifiers: OAI78's scalar -i/81 abar(c)^2 and OAI78's kernel for the CM term; and alpha(mu) instead of alpha(mu)^2
        f31 += [abs(r2/(pred2*1j*ac**-2) - 1),
                abs(r2/(pred2*Rker(t2, mp.mpf(7)/6, mp.mpf(5)/6)/Rker(t2, mp.mpf(4)/3, mp.mpf(5)/3)) - 1),
                abs(r2/(pred2/amu) - 1)]
E31 = max(e31)
claim("C31", 1224, "table: reflected kernel and cusp scalar", "E",
      "reflected/direct = -(1/81)abar(c)^4 q^-1/2 (Kq/q_c^2)^-t R_CM(t) alpha(mu)^2; OAI78 j=1 analogue -(i/81)abar(c)^2 R_1(t)"
      " reproduced as calibration", E31 < 1e-15 and calib1 < 1e-25 and calib2 < 1e-25,
      f"max |ratio-1| {mp.nstr(E31, 3)} (2 mu x 2 s); direct constants 3^(5/2)/4 and 3pi/2 err {mp.nstr(max(calib1, calib2), 2)};"
      f" falsifiers (OAI scalar, OAI kernel, alpha^1): min |ratio-1| {mp.nstr(min(f31), 3)}",
      falsifier_failed=min(f31) > 1e-3)
claim("C32", 1228, "table: dual angular factor", "E", "dual side carries alpha(mu)^2 (CM) vs alpha(mu) (OAI78)", E31 < 1e-15,
      "from C31 (falsifier alpha(mu)^1 fails there)", falsifier_failed=min(f31) > 1e-3)

# ---------------------------------------------------------------- C33h
def rightmost_pole(shifts):
    # poles of prod Gamma(a + t) at t = -a - k
    return max(-a for a in shifts)
cm_pole = rightmost_pole([Fr(4, 3), Fr(5, 3)]); oai_pole = rightmost_pole([Fr(7, 6), Fr(5, 6)])
ok = cm_pole == Fr(-4, 3) and oai_pole == Fr(-5, 6) and Fr(-1, 4) > max(cm_pole, oai_pole)
claim("C33h", 1233, "holomorphy of kernels", "I", "CM kernel holomorphic for Re t > -4/3, OAI78's for Re t > -5/6; strip "
      "-1/4 <= Re t contains no pole", ok, f"rightmost poles {cm_pole}, {oai_pole} (denominators 1/Gamma entire); "
      f"falsifier: CM pole at -5/6? {cm_pole == Fr(-5, 6)}", falsifier_failed=not (cm_pole == Fr(-5, 6)))

# ---------------------------------------------------------------- C34, C36
mp.mp.dps = 20
R1 = lambda t: Rker(t, mp.mpf(7)/6, mp.mpf(5)/6)
Rc = lambda t: Rker(t, mp.mpf(4)/3, mp.mpf(5)/3)
mx = 0; mn = 1e9; slopes = []; fmx = 0
for sig in [-0.25, 0, 0.3, 1, 2, 3]:
    for RR in (R1, Rc):
        vals = []
        for u in [0, 0.5, 1, 3, 10, 30, 100, 1e3, 1e4, 1e5, 1e6]:
            r = abs(RR(mp.mpc(sig, u)))/(1 + u)**(4*sig)
            vals.append(r)
            fmx = max(fmx, abs(RR(mp.mpc(sig, u)))/(1 + u)**(2*sig))
        mx = max(mx, max(vals)); mn = min(mn, min(vals))
        if sig != 0:
            u1, u2 = 1e5, 1e6
            slopes.append((mp.log(abs(RR(mp.mpc(sig, u2)))) - mp.log(abs(RR(mp.mpc(sig, u1)))))/mp.log(u2/u1)/(4*sig))
ok34 = mx < 50 and all(abs(sl - 1) < 1e-3 for sl in slopes)
claim("C34", 1232, "kernel growth", "A",
      "|R(sigma+iu)|/(1+|u|)^{4 sigma} bounded for sigma in {-1/4,0,0.3,1,2,3}, u in [0,1e6], both kernels", ok34,
      f"ratio in [{mp.nstr(mn, 3)}, {mp.nstr(mx, 3)}]; fitted exponent/(4 sigma) in [{mp.nstr(min(slopes), 6)}, {mp.nstr(max(slopes), 6)}]; "
      f"falsifier (1+|u|)^{{2 sigma}}: ratio reaches {mp.nstr(fmx, 3)}; the constant is not tested beyond u=1e6",
      falsifier_failed=fmx > 1e6)
# direct/reflected gamma products on Re s < 0: |Gamma(5/3-s)Gamma(4/3-s)/(Gamma(s+1/3)Gamma(s+2/3))| ~ |t|^{2-4 sigma} (OAI)
# and |Gamma(7/3-s)Gamma(8/3-s)/(Gamma(s+1/3)Gamma(s+2/3))| ~ |t|^{4-4 sigma} (CM, in the J-variable)
sl = []
for sig in [-0.5, -1.0]:
    a = [abs(mp.gamma(mp.mpf(5)/3 - mp.mpc(sig, u))*mp.gamma(mp.mpf(4)/3 - mp.mpc(sig, u))/(mp.gamma(mp.mpc(sig, u) + third)*mp.gamma(mp.mpc(sig, u) + 2*third))) for u in (1e5, 1e6)]
    b = [abs(mp.gamma(mp.mpf(7)/3 - mp.mpc(sig, u))*mp.gamma(mp.mpf(8)/3 - mp.mpc(sig, u))/(mp.gamma(mp.mpc(sig, u) + third)*mp.gamma(mp.mpc(sig, u) + 2*third))) for u in (1e5, 1e6)]
    sl.append((mp.log(a[1]/a[0])/mp.log(10), 2 - 4*sig)); sl.append((mp.log(b[1]/b[0])/mp.log(10), 4 - 4*sig))
e36 = max(abs(p - q) for p, q in sl)
claim("C36", 1254, "Stirling exponents", "E",
      "gamma-quotient exponents: kernels 4 sigma; reflected/direct products 2-4 sigma (OAI78) and 4-4 sigma (CM, J-variable)",
      e36 < 1e-3 and ok34, f"max exponent deviation {mp.nstr(e36, 3)} at u=1e5..1e6; falsifier exponent 2 sigma rejected in C34",
      falsifier_failed=fmx > 1e6)

# ---------------------------------------------------------------- C35
mp.mp.dps = 25
def Vsharp(xx, RR):
    f = lambda u: mp.e**((mp.mpc(0, u))**2)*RR(mp.mpc(0, u))*(K*xx)**(-mp.mpc(0, u))
    return mp.re(mp.quad(f, [-mp.inf, -3, 0, 3, mp.inf])/(2*mp.pi))
rows = []; ok35 = True; fal = []
for RR, name, first_pole in [(R1, "OAI", mp.mpf(5)/6), (Rc, "CM", mp.mpf(4)/3)]:
    ratios = []; small = []
    for xx in [mp.mpf(10)**k for k in (-8, -6, -4, -2, 0, 2, 4)]:
        V = abs(Vsharp(xx, RR))
        ratios.append(V/min(xx**0.25, xx**-2))
        if xx < 1e-3:
            small.append(V/xx**(first_pole + mp.mpf(1)/6))     # falsifier: decay faster than the first pole allows
    rows.append(f"{name}: max |V#|/min(x^1/4,x^-2) = {mp.nstr(max(ratios), 3)}; |V#|/x^(pole+1/6) at x=1e-8: {mp.nstr(small[0], 3)}")
    ok35 &= max(ratios) < 10
    fal.append(small[0] > 10)
claim("C35", 1234, "reflected weight bound", "A",
      "V#(x) << min(x^{1/4}, x^{-A}) for both kernels (Gaussian test, A=2, x in 1e-8..1e4)", ok35,
      "; ".join(rows) + "; falsifier: a bound x^(first pole + 1/6) as x->0 fails (ratio large)",
      falsifier_failed=all(fal))
