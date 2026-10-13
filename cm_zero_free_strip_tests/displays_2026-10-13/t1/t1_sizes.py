"""Agent t1, new test T1-S: size checks for the uniform bound eq:d023 (5.8) = eq:d025 (5.10) (class A: the implied
constant is not explicit, so only finite sizes and the explicit sub-claims are checked).

T1-S1 (exact) N_0^(1/12) = |4d|^(5/6) for N_0 = |4d|^10 (sympy, exact), falsifier |4d|^(5/12).
T1-S2 (explicit inequality behind "when Y < 7 the right side exceeds the trivial bound"): with Y^6 = D^(1+theta)/N_0,
      Y < 7 implies N_0^(1/12) D^(11/12+5theta/12) > 7^(-1/2) D^(1+theta/2) >= D/sqrt 7.  Checked exactly (sympy) as the
      identity N_0^(1/12) D^(11/12+5theta/12) = D^(1+theta/2) Y^(-1/2), plus a grid; falsifier: "> D^(1+theta/2)".
T1-S3 (sizes): S_d(D) for 0 < |d| <= 30, D in {1e2, 1e3, 1e4}, the complex weight W1 = bump(1,2) e^{iy}, theta = 1/10,
      epsilon = 0: the ratio |S_d(D)| / (N_0^(1/12) D^(11/12+5theta/12) + N_0^(1/6) D^(5/6-theta/6)) is printed (sanity,
      no criterion), and the explicit trivial bound |S_d(D)| <= (1 + 8bD) ||W||_inf of Lemma lem:b1-extract(i) is checked.
"""
import math
import numpy as np
import sympy as sp
from t1lib import check, sanity, gp_ap, qd_array

# T1-S1
d, D, th, N0, Y = sp.symbols('d D theta N0 Y', positive=True)
ok1 = sp.simplify((((4*d)**10)**sp.Rational(1, 12)) - (4*d)**sp.Rational(5, 6)) == 0
bad1 = sp.simplify((((4*d)**10)**sp.Rational(1, 12)) - (4*d)**sp.Rational(5, 12)) != 0
check('T1-S1', 'eq:d023', '5.8', ok1, 'sympy: (|4d|^10)^(1/12) = |4d|^(5/6) exactly; falsifier |4d|^(5/12)', falsifier=bad1)

# T1-S2
Yexpr = (D**(1 + th)/N0)**sp.Rational(1, 6)
first = N0**sp.Rational(1, 12)*D**(sp.Rational(11, 12) + 5*th/12)
ok2 = sp.simplify(sp.powsimp(first/(D**(1 + th/2)*Yexpr**sp.Rational(-1, 2)), force=True)) == 1
viol = 0; wrong_ok = 0
for Dv in [1.0, 2.0, 10.0, 1e3, 1e6, 1e12]:
    for tv in [0.01, 0.05, 0.1]:
        for Yv in [0.1, 1.0, 3.0, 6.99]:
            N0v = Dv**(1 + tv)/Yv**6
            f = N0v**(1/12)*Dv**(11/12 + 5*tv/12)
            viol += not (f > 7**-0.5*Dv**(1 + tv/2)*(1 - 1e-12) and 7**-0.5*Dv**(1 + tv/2) >= Dv/math.sqrt(7)*(1 - 1e-12))
            wrong_ok += f > Dv**(1 + tv/2)
check('T1-S2', 'eq:d023', '5.8', ok2 and viol == 0,
      f'sympy: first term = D^(1+theta/2) Y^(-1/2) exactly; Y<7 gives first term > D^(1+theta/2)/sqrt7 >= D/sqrt7 '
      f'on a grid of 72 (D, theta, Y) (violations {viol}); falsifier "> D^(1+theta/2)" holds only in {wrong_ok}/72',
      falsifier=wrong_ok < 72)

# T1-S3
theta = 0.1
W = lambda y: np.where((y > 1) & (y < 2), np.exp(-1.0/np.clip((y - 1)*(2 - y), 1e-300, None) + 4), 0)*np.exp(1j*y)
Wmax = float(np.max(np.abs(W(np.linspace(1, 2, 200001)))))
b = 2.0
Ds = [100, 1000, 10000]
Mq = 2*max(Ds) + 2
worst = 0.0; worst_d = None; trivial_viol = 0; ratios_by_D = {Dv: 0.0 for Dv in Ds}
for dv in [x for x in range(-30, 31) if x != 0]:
    ap = gp_ap(dv, Mq)
    q = qd_array(dv, Mq, ap)
    N0v = float(abs(4*dv))**10
    for Dv in Ds:
        m = np.arange(1, int(b*Dv) + 1)
        S = complex(np.sum(q[m]*W(m/Dv)))
        bound = N0v**(1/12)*Dv**(11/12 + 5*theta/12) + N0v**(1/6)*Dv**(5/6 - theta/6)
        r = abs(S)/bound
        ratios_by_D[Dv] = max(ratios_by_D[Dv], r)
        if r > worst:
            worst, worst_d = r, (dv, Dv)
        trivial_viol += abs(S) > (1 + 8*b*Dv)*Wmax
check('T1-S3a', 'eq:d023', '5.8', trivial_viol == 0,
      f'explicit trivial bound |S_d(D)| <= (1+8bD)||W||_inf (Lemma b1-extract(i)) for 60 values of d, D in {Ds} '
      f'(violations {trivial_viol})')
sanity('T1-S3b', 'eq:d023', '5.8',
       'max |S_d(D)|/(N_0^(1/12)D^(11/12+5theta/12)+N_0^(1/6)D^(5/6-theta/6)), theta=0.1, eps=0, 0<|d|<=30: '
       + ', '.join(f'D={Dv}: {ratios_by_D[Dv]:.2e}' for Dv in Ds) + f' (max at d,D={worst_d})')
