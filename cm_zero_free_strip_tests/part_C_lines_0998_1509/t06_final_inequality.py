"""t06: the exponent bookkeeping of Part II (guide L1074-1092, L1246-1248, L1306-1311, L1336-1342, L1393-1407,
L1430-1460), verified exactly with Fractions / sympy, against OAI78 L505-516, L8641/L8794, L12531-12570, L15393-15405,
L15726-15733, L15920-16113.

Claims tested
  C51 L1092, L1247  sigma_0 = 7/8 in (1/2,1), C(s) = s - 11/16, C(7/8) = 3/16 = exponent of the compensated low estimate;
                    the design values (b, ell, h, l_x, l_y) = (1/8, 1/6, 13/16, 17/48, 23/48) reproduce C and 3/16.
  C52 L1246, L1455  balanced probe: C_I(s) = s - 2/3, C_I(11/12) = 1/4 = exponent of the balanced low estimate.
  C53 L1311, L1454  bin a <= beta_* <= 11/12 gives delta = 2a-1 <= 5/6.
  C54 L1338         |Q_i| << P_i^{a-1/2} = P_i^{delta/2}; q = x delta with 0 <= x <= 1/2.
  C55 L1342, L1393, L1407, L1423  kappa = 3/4 = 2(7/8)-1 lies in [3/4,1]; the plain-moment hypothesis
                    beta <= (1+kappa)/2 = 7/8 is exactly OpenAI's zero-free region; kappa=3/4 < 2 beta_*-1 under beta_*>7/8.
  C56 L1431-1432    OAI78 eq:common-high-exponent: its two forms agree (sympy), C_0 = -1/48.
  C57 L1437         numerator size on the reflected line: U^{(1-2(1-a))/2} = U^{delta/2}.
  C58 L1442-1446    endpoint certificate -E(h) >= 49/440640 on 0<=delta<=5/6, 0<=x<=1/2, q = x delta, R = R_*:
                    (i) derivation of -E_* (sympy); (ii) OAI78's completed-square identity; (iii) the nonnegativity argument;
                    (iv) exact rational grid; (v) true minimum (numerics).
  C59 L1447-1451    adaptive endpoint: 1 - h/4 = 51/64 and E_actual - Delta <= -49/440640 - (1-h/4)Delta + eps < 0.
  C60 L1452-1453    the term h Delta/4 is the capacity loss from kappa = 2 beta_*-1 (OAI78 eq:capacity-comparison);
                    with kappa = 3/4 it is absent, which only lowers the bound.
  C61 L1076-1077, L1458-1460  logic of the contradiction: margins uniform, criterion hypotheses met for every Delta>0.
Run time: about 10 s.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmreport import claim, cited, info, untested
import sympy as sp
import numpy as np

# ---------------------------------------------------------------- C51
b, ell, h, lx, ly = Fr(1, 8), Fr(1, 6), Fr(13, 16), Fr(17, 48), Fr(23, 48)
s0 = Fr(7, 8)
C = lambda s: s - Fr(11, 16)
ok = (Fr(1, 2) < s0 < 1) and C(s0) == Fr(3, 16)
ok &= (-Fr(5, 6) + lx/3 + ell/6 == -Fr(11, 16))            # G78: C(s) = s - 5/6 + l_x/3 + ell/6
ok &= (lx/2 - 1 + h/6 == -Fr(11, 16))                      # OAI78 L8794
ok &= (lx/2 + b/12 == Fr(3, 16))                           # low estimate exponent (G78; OAI78 prop:probe-low)
ok &= (h == 1 - lx + ell) and (lx + ly + ell == 1) and (ly == lx + b)
ok &= (1 - h/6 + b/12 == s0) and (Fr(11, 12) - ell/4 == s0)
fals = (Fr(7, 8) - Fr(2, 3) == Fr(3, 16))                   # C_I instead of C_II
claim("C51", 1092, "sigma_0, C(s), 3/16", "I",
      "sigma_0=7/8, C(s)=s-11/16, C(7/8)=3/16 = l_x/2+b/12; design b=1/8, ell=1/6, h=13/16, l_x=17/48, l_y=23/48 consistent",
      ok, "exact Fractions (OAI78 L514, L8641, L8794; G78 eq:geometry)", falsifier_failed=not fals)

# ---------------------------------------------------------------- C52
CI = lambda s: s - Fr(2, 3)
ok = CI(Fr(11, 12)) == Fr(1, 4) and Fr(1, 2)*Fr(1, 2) == Fr(1, 4)       # X = Z^{1/2}, low estimate X^{1/2}
claim("C52", 1246, "balanced low estimate Z^{1/4}", "I", "C_I(11/12) = 11/12 - 2/3 = 1/4 = (1/2)(1/2) (X^{1/2}, X=Z^{1/2})", ok,
      "exact (OAI78 L509-510)", falsifier_failed=not (CI(Fr(7, 8)) == Fr(1, 4)))

# ---------------------------------------------------------------- C53, C54
ok53 = all(2*a - 1 <= Fr(5, 6) for a in [Fr(k, 1200) for k in range(612, 1101)]) and 2*Fr(11, 12) - 1 == Fr(5, 6)
claim("C53", 1454, "bin ceiling delta <= 5/6", "I", "a <= beta_* <= 11/12 implies delta = 2a-1 <= 5/6", ok53,
      "exact; falsifier a = 11/12 + 1/1200 gives delta > 5/6", falsifier_failed=(2*(Fr(11, 12) + Fr(1, 1200)) - 1 > Fr(5, 6)))
a = sp.symbols('a')
ok54 = sp.simplify((a - sp.Rational(1, 2)) - (2*a - 1)/2) == 0
claim("C54", 1338, "P_i^{delta/2}", "I", "OAI78 lem:bin-prime-bound exponent a-1/2 equals delta/2; q<=delta/2 <=> x in [0,1/2]",
      ok54, "sympy; falsifier exponent delta (=2a-1) differs", falsifier_failed=sp.simplify((a - sp.Rational(1, 2)) - (2*a - 1)) != 0)

# ---------------------------------------------------------------- C55
kap = 2*Fr(7, 8) - 1
ok55 = kap == Fr(3, 4) and Fr(3, 4) <= kap <= 1 and (1 + kap)/2 == Fr(7, 8)
D = sp.symbols('Delta', positive=True)
k_oai = 2*(sp.Rational(7, 8) + D) - 1
ok55 &= sp.simplify(k_oai - (sp.Rational(3, 4) + 2*D)) == 0
# smaller kappa enlarges the plain capacity (1-2m)/(6 kappa) and weakens n1+n2+6 kappa z <= M
m = sp.symbols('m')
ok55 &= sp.simplify((1 - 2*m)/(6*sp.Rational(3, 4)) - (1 - 2*m)/(6*k_oai) - (1 - 2*m)*12*D/(sp.Rational(9, 2)*(sp.Rational(9, 2) + 12*D))) == 0
claim("C55", 1342, "kappa = 3/4", "I",
      "kappa=2(7/8)-1=3/4 in [3/4,1]; (1+kappa)/2 = 7/8 (OpenAI's theorem); OAI's kappa = 3/4+2Delta > 3/4 under beta_*>7/8",
      ok55, "exact/sympy (OAI78 L12565, L15940-15960)", falsifier_failed=not ((1 + Fr(5, 6))/2 == Fr(7, 8)))

# ---------------------------------------------------------------- C56
dl, q, R, d = sp.symbols('delta q R d')
aa = (1 + dl)/2
hS, ellS, lyS = sp.Rational(13, 16), sp.Rational(1, 6), sp.Rational(23, 48)
E1 = aa - sp.Rational(7, 8) + hS*(sp.Rational(17, 50) - sp.Rational(1, 6)) - aa*lyS - (1 - aa)*ellS - (dl/2 - q)*ellS \
     + d*(R + dl/2 - sp.Rational(17, 50))
E2 = -sp.Rational(1, 48) + sp.Rational(2, 3)*dl + q/6 - hS*(1 - R) + (d - hS)*(R + dl/2 - sp.Rational(17, 50))
ok56 = sp.simplify(E1 - E2) == 0
E2f = E2 + sp.Rational(1, 96)
claim("C56", 1431, "E(d) two forms", "E", "OAI78 eq:common-high-exponent: first form = second form with C_0 = -1/48", ok56,
      "sympy, symbolic in delta, q, R, d", falsifier_failed=sp.simplify(E1 - E2f) != 0)

# ---------------------------------------------------------------- C57
ok57 = sp.simplify((1 - 2*(1 - a))/2 - (2*a - 1)/2) == 0
claim("C57", 1437, "numerator U^{delta/2}", "I", "FE on Re w = 1-a: conductor power U^{(1-2(1-a))/2} = U^{delta/2}", ok57,
      "sympy; falsifier exponent (1-a)/2", falsifier_failed=sp.simplify((1 - a)/2 - (2*a - 1)/2) != 0)

# ---------------------------------------------------------------- C58 certificate
x, y = sp.symbols('x y')
al = sp.Rational(5, 6)
Dx = 3 - sp.Rational(17, 9)*x
Px = (2 - sp.Rational(8, 9)*x)*(1 - x)
J = (al - dl)*Dx + dl*Px
t = 1 + dl*Px/(2*J)
Rstar = 1 - dl + (al - dl)*(t - 1)
Rshort = 1 - dl + dl*Px/Dx*(sp.Rational(3, 2) - t)
Eh = E2.subs({d: hS, q: x*dl, R: Rstar})
mEstar_oai = (1 + (3 + 8*y)*dl)/48 - sp.Rational(13, 32)*(al - dl)*dl*Px/J
c1 = sp.simplify((Rstar - Rshort)) == 0
c2 = sp.simplify((-Eh) - mEstar_oai.subs(y, sp.Rational(1, 2) - x)) == 0
yy = sp.Rational(1, 2) - x
py = 7 + 18*y + 8*y**2; jy = 185 + 170*y + (-138 + 12*y + 96*y**2)*dl; v = 51 + 41*y
c3 = sp.simplify(Dx - ((37 + 34*y)/18).subs(y, yy)) == 0 and sp.simplify(Px - (py/9).subs(y, yy)) == 0 \
     and sp.simplify(J - (jy/108).subs(y, yy)) == 0
quad = 10*(37 + 34*y) - 8*(237 + 377*y + 26*y**2)*dl + 48*(51 + 131*y + 94*y**2 + 32*y**3)*dl**2
c4 = sp.expand(2*jy*(1 + (3 + 8*y)*dl) - 468*(al - dl)*dl*py - quad) == 0
rhs = (3 + 5*y)*((4*v*dl - 79)**2 + 49) + 4*y*(4*v*dl*((1 + 3*y)*(15 + 32*y)*dl + 9 - 13*y) + 265 + 3485*y)
c5 = sp.expand(v*quad - rhs) == 0
c6 = sp.simplify(10368*J*(-Eh) - quad.subs(y, yy)) == 0
# nonnegativity argument: y in [0,1/2]: 9-13y >= 5/2; v <= 17(3+5y); J <= 5/2 on the box
c7 = all(9 - 13*Fr(k, 200) >= Fr(5, 2) for k in range(101)) and all(51 + 41*Fr(k, 200) <= 17*(3 + 5*Fr(k, 200)) for k in range(101))
Jl = sp.lambdify((x, dl), J)
c8 = True; Jmax = Fr(0)
for i in range(0, 101):
    for k in range(0, 101):
        xv = Fr(i, 200); dv = Fr(5, 6)*Fr(k, 100)
        Jv = (Fr(5, 6) - dv)*(3 - Fr(17, 9)*xv) + dv*(2 - Fr(8, 9)*xv)*(1 - xv)
        Jmax = max(Jmax, Jv); c8 &= Jv <= Fr(5, 2) and Jv > 0
c9 = Fr(49, 176256)*Fr(2, 5) == Fr(49, 440640) and Fr(49, 440640) > Fr(1, 10000)
# (iv) exact rational grid of -E_*
def mE(xv, dv):
    yv = Fr(1, 2) - xv
    Pxv = (2 - Fr(8, 9)*xv)*(1 - xv)
    Jv = (Fr(5, 6) - dv)*(3 - Fr(17, 9)*xv) + dv*Pxv
    return (1 + (3 + 8*yv)*dv)/48 - Fr(13, 32)*(Fr(5, 6) - dv)*dv*Pxv/Jv
gmin = min(mE(Fr(i, 100), Fr(5, 6)*Fr(k, 300)) for i in range(51) for k in range(301))
# (v) true minimum (float, fine grid then local refine)
f = lambda X, Dl: float(mE(Fr(X).limit_denominator(10**9), Fr(Dl).limit_denominator(10**9)))
xs = np.linspace(0, 0.5, 201); ds = np.linspace(0, 5/6, 801)
Xg, Dg = np.meshgrid(xs, ds, indexing="ij")
Ygr = 0.5 - Xg; Pg = (2 - 8*Xg/9)*(1 - Xg); Jg = (5/6 - Dg)*(3 - 17*Xg/9) + Dg*Pg
Mg = (1 + (3 + 8*Ygr)*Dg)/48 - 13/32*(5/6 - Dg)*Dg*Pg/Jg
i0 = np.unravel_index(np.argmin(Mg), Mg.shape)
tmin = Mg[i0]; xw, dw = Xg[i0], Dg[i0]
ok58 = all([c1, c2, c3, c4, c5, c6, c7, c8, c9]) and gmin >= Fr(49, 440640) and tmin >= 49/440640
fal58 = tmin >= 3e-4          # falsifier: a certificate margin 3e-4 would be false
claim("C58", 1444, "endpoint certificate", "I",
      "-E(h) >= 49/440640 for 0<=delta<=5/6, 0<=x<=1/2 (q=x delta, R=R_*): derivation, OAI78 identity, positivity, grid",
      ok58, f"R_*=R_short(t) {c1}; -E_* formula {c2}; D_x,P_x,J substitutions {c3}; quadratic {c4}; completed square {c5}; "
      f"10368 J(-E_*) {c6}; 9-13y>=5/2 & v<=17(3+5y) {c7}; 0<J<=5/2 (max {float(Jmax):.4f}) {c8}; 49/176256*2/5=49/440640>1e-4 {c9}; "
      f"exact grid min {float(gmin):.6e}; true min {tmin:.4e} at x={xw:.3f}, delta={dw:.4f} vs 49/440640={49/440640:.4e}; "
      f"falsifier margin 3e-4 false", falsifier_failed=not fal58)

# ---------------------------------------------------------------- C59
eps_tot = sp.symbols('epsilon_total', nonnegative=True)
cert = sp.Rational(49, 440640)
c_a = (1 - hS/4) == sp.Rational(51, 64)
bound = -cert - (1 - hS/4)*D + eps_tot
# E_actual(h) <= E_* + h*Delta/4 + eps  (row count R <= R_* + Delta/4 + eps_row; coefficient of R in E is d = h)
c_b = sp.simplify((-cert + hS*D/4 + eps_tot - D) - bound) == 0
grid_ok = all(bound.subs({D: sp.Rational(k, 24*100), eps_tot: sp.Rational(48, 440640)}) < 0 for k in range(1, 101))
grid_f = any(bound.subs({D: sp.Rational(k, 10**9), eps_tot: sp.Rational(50, 440640)}) >= 0 for k in range(1, 3))
claim("C59", 1450, "eq:adaptive-endpoint", "I",
      "1-h/4 = 51/64 (h=13/16); E_actual-Delta <= -49/440640-(1-h/4)Delta+eps_total < 0 for 0<Delta<=1/24, eps_total<49/440640",
      c_a and c_b and grid_ok, f"exact (sympy); 1-h/4 = {1 - hS/4}; falsifier eps_total = 50/440640 with Delta -> 0 gives >= 0",
      falsifier_failed=grid_f)

# ---------------------------------------------------------------- C60
mm, qq = sp.symbols('m q', nonnegative=True)
lhs = 2*qq*(1 - 2*mm)*(sp.Rational(2, 9) - 1/(sp.Rational(9, 2) + 12*D))
ident = sp.simplify(lhs - 24*qq*(1 - 2*mm)*D/(sp.Rational(9, 2)*(sp.Rational(9, 2) + 12*D))) == 0
worst = (24*sp.Rational(1, 2)*sp.Rational(1, 3)*D/(sp.Rational(9, 2)*sp.Rational(9, 2)))   # 2q<=1, 1-2m<=1/3, D>=0
le = sp.simplify(D/4 - worst) > 0
baseline = sp.Rational(2, 9) == 1/(6*sp.Rational(3, 4))
only_helps = sp.simplify((-cert - D + eps_tot) - bound) <= 0      # = -h D/4 <= 0
claim("C60", 1452, "kappa=3/4 removes h Delta/4", "I",
      "capacity loss 2q(1-2m)(2/9 - 1/(9/2+12 Delta)) <= Delta/4 is the only Delta-dependent row loss; the baseline 2/9 = 1/(6*3/4)"
      " is kappa = 3/4; dropping h Delta/4 lowers the bound", ident and bool(le) and baseline and bool(only_helps),
      f"identity {ident}; worst case {sp.nsimplify(worst/D)}*Delta <= Delta/4; baseline 2/9 = 1/(6 kappa) at kappa=3/4: {baseline}; "
      "falsifier: kappa = 5/6 baseline 1/(6 kappa) = 1/5 != 2/9", falsifier_failed=(1/(6*sp.Rational(5, 6)) != sp.Rational(2, 9)))

# ---------------------------------------------------------------- C61
# continuation criterion: need 0 < omega < Delta_0 = beta_* - 7/8 and sigma > 0, independent of the target.
# low: |J| << Z^{C(7/8)+eps}  -> omega = eps;  high: |J - f| << Z^{C(7/8)+E_actual} = Z^{C(beta_*) + (E_actual - Delta)}
ok61 = True
for Dv in [Fr(1, 10**6), Fr(1, 1000), Fr(1, 100), Fr(1, 24)]:
    eps = min(Dv/2, Fr(48, 440640))
    sigma = Fr(49, 440640) + Fr(51, 64)*Dv - eps          # = -(upper bound for E_actual - Delta), with kappa=2beta_*-1
    sigma_k34 = Fr(49, 440640) + Dv - eps                  # with kappa = 3/4
    ok61 &= (0 < eps < Dv) and sigma > 0 and sigma_k34 >= sigma
# falsifier: without the subtraction of Delta (comparing to C(7/8) instead of C(beta_*)), the row loss h Delta/4 alone
# can exceed the certificate margin, so "E_actual(h) < 0 for every Delta" is not what is proved
ff = Fr(49, 440640) - Fr(13, 16)*Fr(1, 24)/4 < 0
claim("C61", 1460, "contradiction logic", "I",
      "for every Delta in (0,1/24]: omega=eps<Delta and sigma = 49/440640+(51/64)Delta-eps > 0; kappa=3/4 gives larger sigma",
      ok61, "exact; falsifier: measured against C(7/8) (no -Delta), the row loss h Delta/4 alone can exceed the certificate "
      f"(49/440640 - h(1/24)/4 < 0: {ff})", falsifier_failed=ff)
