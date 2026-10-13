"""Part D, test 02: the j-th derivative chain rule (lines 1568-1574) and its mirror image (lines 1641-1642).

With z' = -delta/c - zbar/(c^2 (v^2 + z zbar)),  zbar' = -deltabar/cbar - z/(cbar^2 (v^2 + z zbar)),
v' = v/(c cbar (v^2 + z zbar)) (the coordinates of eq. (chain), line 754), the paper claims, exactly,

  d^j/dzbar^j F(z', zbar', v') |_{z=zbar=0} = (-(c v)^{-2})^j (d^j_{z'} F)(-delta/c, 1/(N(c) v))      (L1571)
  d^j/dz^j    F(z', zbar', v') |_{z=zbar=0} = (-(cbar v)^{-2})^j (d^j_{zbar'} F)(same point)          (L1641)

and that every derivative of order j >= 1 kills a constant term (a function of v' only) (L1574).

Method: sympy, with z, zbar, c, cbar, delta, deltabar, v independent symbols.  Both sides are linear in F
and depend only on the j-jet of F at the image point, so it suffices (and is a proof for that j) to check
every monomial (X-X0)^a (Y-Y0)^b (V-V0)^d with a+b+d <= j.  The paper checked j <= 4 (zbar) and j <= 3 (z);
here j = 1..6 for both.

Falsifiers: (i) the coefficient (c v)^{-2j} without the sign (must fail for odd j); (ii) evaluating at a
point with z != 0 (z = 1/5) instead of z = 0 (must fail); (iii) for the constant term, j = 0 (must not vanish).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import report
import sympy as sp

z, zb, v, c, cb, dl, dlb = sp.symbols('z zb v c cb delta deltab')
X, Y, V = sp.symbols('X Y V')
zp = -dl/c - zb/(c**2*(v**2 + z*zb))
zbp = -dlb/cb - z/(cb**2*(v**2 + z*zb))
vp = v/(c*cb*(v**2 + z*zb))
X0, Y0, V0 = -dl/c, -dlb/cb, 1/(c*cb*v)
JMAX = 6

def monomials(j):
    out = []
    for a in range(j + 1):
        for b in range(j + 1 - a):
            for d in range(j + 1 - a - b):
                out.append((X - X0)**a*(Y - Y0)**b*(V - V0)**d)
    return out

def lhs(m, var, j, zval=0):
    e = m.subs({X: zp, Y: zbp, V: vp}, simultaneous=True)
    e = sp.diff(e, var, j)
    return e.subs({z: zval, zb: 0})

ok_zb = True; ok_z = True; fal_sign = []; fal_point = []
for j in range(1, JMAX + 1):
    good_zb = good_z = True; sign_fails = False; point_fails = False
    for m in monomials(j):
        L = lhs(m, zb, j)
        R = (-(c*v)**-2)**j*sp.diff(m, X, j).subs({X: X0, Y: Y0, V: V0}, simultaneous=True)
        if sp.simplify(L - R) != 0:
            good_zb = False
        if j % 2 == 1:
            Rf = ((c*v)**-2)**j*sp.diff(m, X, j).subs({X: X0, Y: Y0, V: V0}, simultaneous=True)
            if sp.simplify(L - Rf) != 0:
                sign_fails = True
        if m == (X - X0)**(j - 1)*(V - V0):      # point falsifier: z = 1/5 instead of 0
            Lp = lhs(m, zb, j, zval=sp.Rational(1, 5))
            if sp.simplify(Lp - R) != 0:
                point_fails = True
        Lz = lhs(m, z, j)
        Rz = (-(cb*v)**-2)**j*sp.diff(m, Y, j).subs({X: X0, Y: Y0, V: V0}, simultaneous=True)
        if sp.simplify(Lz - Rz) != 0:
            good_z = False
    ok_zb &= good_zb; ok_z &= good_z
    if j % 2 == 1: fal_sign.append(sign_fails)
    fal_point.append(point_fails)
    print(f"   j={j}: zbar-rule {good_zb}, z-rule {good_z}, sign falsifier fails: {sign_fails if j % 2 else 'n/a'},"
          f" z!=0 falsifier fails: {point_fails}", flush=True)

report("D02.1", 1571, "chain rule d_zbar^j", f"exact for every monomial of degree <= j, j = 1..{JMAX} (paper: j <= 4)",
       ok_zb and all(fal_sign) and all(fal_point),
       "falsifiers: unsigned coefficient fails for odd j; z=1/5 fails for all j")
report("D02.2", 1641, "chain rule d_z^j", f"(-(cbar v)^-2)^j d_zbar'^j, exact for j = 1..{JMAX} (paper: j <= 3)", ok_z,
       "same monomial basis")

# constant term killed: F = h(V) arbitrary function of v' only
h = sp.Function('h')
kill = all(sp.simplify(sp.diff(h(vp), zb, j).subs({z: 0, zb: 0})) == 0 for j in range(1, JMAX + 1))
nonzero0 = sp.simplify(h(vp).subs({z: 0, zb: 0})) != 0
# but at a point z != 0 the first derivative of h(v') does NOT vanish (so the evaluation at z = 0 matters)
nz = sp.simplify(sp.diff(h(vp), zb, 1).subs({z: sp.Rational(1, 5), zb: 0})) != 0
report("D02.3", 1574, "constant term", f"d_zbar^j h(v') = 0 at z = zbar = 0 for j = 1..{JMAX}", kill and nonzero0 and nz,
       "falsifiers: j=0 gives h(1/(N(c)v)) != 0; at z=1/5 the first derivative is nonzero")

# the structural facts quoted at L1568: z' affine in zbar at z=0; zbar', v' independent of zbar at z=0
aff = all(sp.simplify(sp.diff(zp.subs(z, 0), zb, k)) == 0 for k in range(2, 6)) and \
      sp.simplify(sp.diff(zp, zb).subs({z: 0, zb: 0}) + (c*v)**-2) == 0 and \
      sp.simplify(sp.diff(zbp.subs(z, 0), zb)) == 0 and sp.simplify(sp.diff(vp.subs(z, 0), zb)) == 0 and \
      sp.simplify(vp.subs({z: 0}) - 1/(c*cb*v)) == 0
fal = sp.simplify(sp.diff(zp, zb, 2).subs({z: sp.Rational(1, 5), zb: 0})) != 0
report("D02.4", 1568, "coordinates", "at z=0: z' affine in zbar with slope -(cv)^-2; zbar', v' independent of zbar; v' = 1/(N(c) v)",
       aff and fal, "falsifier: at z=1/5, d^2 z'/dzbar^2 != 0")

# L1640-1641: d_z^n ce(l z) = (2 pi i l)^n ce(l z), so d_z^{k-1} produces alpha(l)^{k-1} (times |l|^{k-1}); d_zbar gives conj(l)
l, lb = sp.symbols('l lb')
ce_ = sp.exp(2*sp.pi*sp.I*(l*z + lb*zb))
ok = all(sp.simplify(sp.diff(ce_, z, n) - (2*sp.pi*sp.I*l)**n*ce_) == 0 for n in range(1, 7)) and \
     all(sp.simplify(sp.diff(ce_, zb, n) - (2*sp.pi*sp.I*lb)**n*ce_) == 0 for n in range(1, 7))
fal = sp.simplify(sp.diff(ce_, z, 2) - (2*sp.pi*sp.I*lb)**2*ce_) != 0
report("D02.5", 1640, "d_z on ce(lz)", "d_z^n ce(lz) = (2 pi i l)^n ce(lz) and d_zbar^n ce(lz) = (2 pi i lbar)^n ce(lz), n=1..6",
       ok and fal, "falsifier: d_z giving lbar fails")
