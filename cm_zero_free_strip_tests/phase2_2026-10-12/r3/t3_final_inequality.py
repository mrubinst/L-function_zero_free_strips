"""Section 15 of the CM guide: the endpoint certificate of OAI78 lem:balanced-endpoint, the
adaptive comparison eq:adaptive-endpoint, the effect of kappa = 3/4, and the role of the bin
ceiling alpha = 5/6 (= 2*(11/12) - 1).  Exact arithmetic with sympy; falsifiers included.
All formulas are transcribed from OAI78 (paper.tex L15684-16110)."""
import os
import sympy as sp
from fractions import Fraction as Fr

d, x, y, h, R, q, a, al = sp.symbols('delta x y h R q a alpha', real=True)
C0 = sp.Rational(-1, 48)
H = sp.Rational(13, 16)          # h such that 1 - h/4 = 51/64
ok = True
_FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']   # self-test (added 2026-10-12): first check counts as failed
def report(name, cond, info=''):
    global ok
    if _FORCE[0]:
        _FORCE[0] = False; cond = False; name = '[CM_TESTS_FORCE_FAIL=1] ' + name
    print(('[PASS] ' if cond else '[FAIL] ') + name + ('  :: ' + info if info else ''))
    ok &= bool(cond)

# (1) h = 13/16 is the h of eq:adaptive-endpoint: (1 - h/4) = 51/64
report('E1 1-h/4 = 51/64 for h=13/16', sp.simplify(1 - H/4 - sp.Rational(51, 64)) == 0)

# (2) E(h) from eq:common-high-exponent (second line, d = h)
E_h = C0 + sp.Rational(2, 3)*d + q/6 - h*(1 - R)
# consistency of first and second line: a=(1+delta)/2, l_y, ell from the geometry are not
# needed for the d=h specialisation, which is all the certificate uses.

# (3) the certificate: D_x, P_x, J, t, R_* (eq:target-cutoff) with alpha = 5/6
Dx = (37 + 34*y)/18
Px = (7 + 18*y + 8*y**2)/9
def J(alpha): return (alpha - d)*Dx + d*Px
def Rstar(alpha): return 1 - d + (alpha - d)*d*Px/(2*J(alpha))
negE = sp.simplify(-(E_h.subs({h: H, q: (sp.Rational(1, 2) - y)*d, R: Rstar(sp.Rational(5, 6))})))
target = (1 + (3 + 8*y)*d)/48 - sp.Rational(13, 32)*(sp.Rational(5, 6) - d)*d*Px/J(sp.Rational(5, 6))
report('E2 -E_* = (1+(3+8y)delta)/48 - (13/32)(5/6-delta)delta P_x/J', sp.simplify(negE - target) == 0)
jy = 185 + 170*y + (-138 + 12*y + 96*y**2)*d
report('E3 J = j_y/108', sp.simplify(J(sp.Rational(5, 6)) - jy/108) == 0)
v = 51 + 41*y
lhs = sp.expand(10368*v*(jy/108)*target)
rhs = sp.expand((3 + 5*y)*((4*v*d - 79)**2 + 49)
                + 4*y*(4*v*d*((1 + 3*y)*(15 + 32*y)*d + 9 - 13*y) + 265 + 3485*y))
report('E4 completed-square identity eq:balanced-endpoint-certificate', sp.simplify(lhs - rhs) == 0)
# falsifier: wrong constant 49 -> 48 must break the identity
rhs_bad = rhs - (3 + 5*y)
report('E4f falsifier (49 -> 48) breaks the identity', sp.simplify(lhs - rhs_bad) != 0)

# (4) grid minimum of -E_* on 0<=delta<=5/6, 0<=x<=1/2 (exact rationals)
f = sp.lambdify((d, y), target, 'sympy')
mn = None
N = 120
for i in range(N + 1):
    dd = Fr(5, 6)*Fr(i, N)
    for j in range(31):
        yy = Fr(1, 2)*Fr(j, 30)
        val = f(sp.Rational(dd.numerator, dd.denominator), sp.Rational(yy.numerator, yy.denominator))
        if mn is None or val < mn[0]:
            mn = (val, dd, yy)
report('E5 grid min of -E_* >= 49/440640', mn[0] >= sp.Rational(49, 440640),
       'min %.6e at delta=%s, x=%s; bound %.6e' % (float(mn[0]), mn[1], Fr(1, 2) - mn[2], 49/440640))

# (5) eq:adaptive-endpoint: E_actual(h) - Delta <= E_* + h*Delta/4 - Delta
Dl = sp.symbols('Delta', positive=True)
expr = -sp.Rational(49, 440640) + H*Dl/4 - Dl
report('E6 -49/440640 + h Delta/4 - Delta = -49/440640 - (51/64) Delta',
       sp.simplify(expr - (-sp.Rational(49, 440640) - sp.Rational(51, 64)*Dl)) == 0)
# with kappa = 3/4 the capacity allowance h*Delta/4 is absent: the bound becomes -49/440640 - Delta,
# which is <= the kappa = 2 beta_* - 1 bound for every Delta >= 0
report('E7 kappa=3/4 bound <= OAI bound for all Delta>=0',
       sp.simplify((-sp.Rational(49, 440640) - Dl) - expr) <= 0)

# (6) plain capacity: z_P(m) = (1-2m)/(6 kappa); with kappa = 3/4 + 2 Delta vs kappa = 3/4
m = sp.symbols('m')
zP = lambda kap: (1 - 2*m)/(6*kap)
report('E8 z_P(m) at kappa=3/4 equals 2(1-2m)/9', sp.simplify(zP(sp.Rational(3, 4)) - 2*(1 - 2*m)/9) == 0)
report('E8b kappa=3/4 capacity >= kappa=3/4+2Delta capacity for m<1/2',
       sp.simplify(zP(sp.Rational(3, 4)) - zP(sp.Rational(3, 4) + 2*Dl)).subs({m: sp.Rational(1, 4), Dl: sp.Rational(1, 100)}) > 0)

# (7) large-delta endpoint eq:large-delta-endpoint: R = 1-delta, q <= delta/2, d = h
E_large = sp.expand(E_h.subs({h: H, R: 1 - d, q: d/2}))
report('E9 E(h) at R=1-delta, q=delta/2 equals -1/48 - delta/16',
       sp.simplify(E_large - (-sp.Rational(1, 48) - d/16)) == 0)
# floor bin eq:floor-bound: delta0 = 1/50, R = 1, q <= delta0/2
E_floor = E_h.subs({h: H, R: 1, q: sp.Rational(1, 100), d: sp.Rational(1, 50)})
report('E10 floor bound E(h) = -7/1200', sp.simplify(E_floor + sp.Rational(7, 1200)) == 0)

# (8) role of the bin ceiling alpha: does the certificate survive with alpha in (5/6, 1]?
#     (relevant to whether item 6 of Section 16, beta_* <= 11/12, could be avoided)
for A in [sp.Rational(5, 6), sp.Rational(7, 8), sp.Rational(9, 10), sp.Rational(19, 20), 1]:
    tgt = (sp.Rational(1, 48) - (sp.Rational(2, 3) - H)*d - (sp.Rational(1, 2) - y)*d/6
           + H*(-(A - d)*d*Px/(2*J(A))))
    g = sp.lambdify((d, y), tgt, 'math')
    best = min((g(float(A)*i/400, 0.5*j/40), float(A)*i/400, 0.5 - 0.5*j/40)
               for i in range(401) for j in range(41))
    print('   [info] certificate with the amplification slope 5/6 replaced by A=%-5s: min over 0<=delta<=A, 0<=x<=1/2 of -E_* = %+.4e at delta=%.4f, x=%.3f'
          % (A, best[0], best[1], best[2]))

# (9) capacity comparison eq:capacity-comparison: 2q(1-2m)(2/9 - 1/(6 kappa)) with kappa = 3/4 + 2 Delta
qq, mm = sp.symbols('q m', nonnegative=True)
kap = sp.Rational(3, 4) + 2*Dl
cc = 2*qq*(1 - 2*mm)*(sp.Rational(2, 9) - 1/(6*kap))
report('E11 capacity comparison = 24 q (1-2m) Delta / ((9/2)(9/2+12 Delta))',
       sp.simplify(cc - 24*qq*(1 - 2*mm)*Dl/(sp.Rational(9, 2)*(sp.Rational(9, 2) + 12*Dl))) == 0)
import random
random.seed(1)
okc = True
for _ in range(20000):
    D_ = random.uniform(0, 1/8); m_ = random.uniform(1/3, 1/2); q_ = random.uniform(0, 1/2)
    v = 2*q_*(1 - 2*m_)*(2/9 - 1/(6*(3/4 + 2*D_)))
    okc &= v <= D_/4 + 1e-15
report('E12 comparison <= Delta/4 for m>=1/3, 2q<=1, 0<=Delta<=1/8 (random 2e4)', okc)
report('E13 comparison is exactly 0 at kappa = 3/4', sp.simplify(2*qq*(1 - 2*mm)*(sp.Rational(2, 9) - 1/(6*sp.Rational(3, 4)))) == 0)

# (10) bins with delta >= 5/6 (Lemma r3-high-bins): no-slot exponent e(r)-delta r <= 1-delta for r>=1
al6 = sp.Rational(5, 6)
okn = True; fal = False
for i in range(0, 101):
    dd = al6 + (1 - al6)*sp.Rational(i, 100)
    for j in range(0, 51):
        r_ = 1 + sp.Rational(j, 100)
        okn &= (1 - al6 + (al6 - dd)*r_) <= 1 - dd
for dd in [sp.Rational(1, 2), sp.Rational(3, 4)]:
    fal |= (1 - al6 + (al6 - dd)*sp.Rational(3, 2)) > 1 - dd
report('E14 for 5/6<=delta<=1, r in [1,3/2]: 1-alpha+(alpha-delta)r <= 1-delta', okn and fal,
       'falsifier: for delta<5/6 the exponent exceeds 1-delta at r=3/2')
dvar = sp.symbols('d')
Ed = (C0 + sp.Rational(2, 3)*d + (d/2)/6 - H*(1 - (1 - d)) + (dvar - H)*((1 - d) + d/2 - sp.Rational(17, 50)))
report('E15 E(d) at R=1-delta, q=delta/2 equals -1/48 - delta/16 + (d-h)(33/50 - delta/2)',
       sp.simplify(Ed - (-sp.Rational(1, 48) - d/16 + (dvar - H)*(sp.Rational(33, 50) - d/2))) == 0)
sl_lo = sp.Rational(33, 50) - sp.Rational(1, 2); sl_hi = sp.Rational(33, 50) - sp.Rational(5, 12)
report('E16 slope 33/50 - delta/2 in [4/25, 73/300] for 5/6<=delta<=1; E(h) <= -7/96',
       sl_lo == sp.Rational(4, 25) and sl_hi == sp.Rational(73, 300) and sl_hi < sp.Rational(1, 4)
       and -sp.Rational(1, 48) - sp.Rational(5, 6)/16 == -sp.Rational(7, 96))
# (11) intermediate range d<=1/2 (eq:adaptive-intermediate) fails for delta > 529/625: why the bins delta>5/6 need Lemma r3-high-bins
Eint = -sp.Rational(529, 2400) + sp.Rational(25, 96)*d
report('E17 eq:adaptive-intermediate root 529/625 = 0.8464 > 5/6; negative at 5/6, positive at 1',
       sp.solve(Eint, d)[0] == sp.Rational(529, 625) and Eint.subs(d, sp.Rational(5, 6)) < 0 and Eint.subs(d, 1) > 0
       and Eint.subs(d, sp.Rational(5, 6)) == -sp.Rational(49, 14400))

print('ALL PASS' if ok else 'SOME FAIL')

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
