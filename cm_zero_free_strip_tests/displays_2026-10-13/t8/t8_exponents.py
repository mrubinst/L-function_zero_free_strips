"""t8_exponents.py -- exact (sympy) checks of the exponent relations attached to four displays of
cm-eleven-twelfths-guide.tex (Appendices E and G).

  EX1  eq:d172 (hypotheses r + 2z <= m - c_1, 2r + 8z <= 3m - c_2 of Prop. prop:s4-marked, identical to
       OAI78 Lemma lem:marked L9252-9254): the application stated in the guide (Section sec:78moments, item 3) --
       eq:MS is eq:cm-inverse with Z = H = D^{1+theta}, m = 1, r = 1/(1+theta), no slots, c_1 = theta/(1+theta) --
       satisfies both hypotheses for 0 < theta <= 1/10 (with c_2 = 1); and OAI78 L12345-12348: H >= D^{1+c} gives
       1 - r >= c/(1+c).  Falsifier: theta = 0 (no margin c_1 > 0).
  EX2  eq:d190 (OAI78 Props. prop:balanced-low/high, L3673-3681, L6583-6589, L6690-6694; C_I(s) = s - 2/3 at L5548,
       omega = Delta_1/2 at L6606): Z^{1/4+omega} = Z^{C_I(11/12)+omega}, 0 < omega < Delta_1, and the low bound is
       below the signal exponent C_I(beta_*) by exactly Delta_1/2.
  EX3  eq:d198 = eq:s1-cm-completed at M = 1: Z^{max(M, 2M-1)} = Z^1, and O/2 + max(M-O, 2M-O-1) <= max(M, 2M-1) for
       O >= 0 (the dyad sum after eq:s1-cm-completed), checked symbolically on the two branches.
  EX4  eq:s6-bound: the extraction exponents (as in eq:exponents with N_0 = 1): sqrt(D^{2+theta}/Y) with Y = H^{1/6},
       H = D^{1+theta} gives 11/12 + 5 theta/12, sqrt(D^2/Y^2) gives 5/6 - theta/6 <= 11/12; so D^{11/12+eps} after
       theta -> 0.  Falsifier: Y = H^{1/3} would give 5/6 + ... (a different, wrong exponent for sixth powers).
Usage: python3 t8_exponents.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import report, falsifier
import sympy as sp

th, c, M, O, D1, beta = sp.symbols('theta c M O Delta_1 beta', real=True)

# ---- EX1
m = 1; r = 1/(1 + th); z = 0; c1 = th/(1 + th); c2 = 1
h1 = sp.simplify((m - c1) - (r + 2*z))            # must be >= 0
h2 = sp.simplify((3*m - c2) - (2*r + 8*z))        # must be >= 0
ok1 = (h1 == 0) and all(h2.subs(th, t) > 0 for t in [sp.Rational(k, 1000) for k in range(1, 101)])
ok1 &= sp.simplify(h2 - 2*th/(1 + th)) == 0       # = 2 theta/(1+theta) > 0
ok1 &= all(c1.subs(th, t) > 0 for t in [sp.Rational(k, 1000) for k in range(1, 101)])
# OAI78 L12345-12348: H >= D^{1+c}  =>  r = log_H D <= 1/(1+c)  =>  1 - r >= c/(1+c)
ok1 &= sp.simplify((1 - 1/(1 + c)) - c/(1 + c)) == 0
report("EX1", 0.0 if ok1 else 1.0, 0.0, "eq:d172: the MS application m=1, r=1/(1+theta), z=0, c_1=theta/(1+theta), c_2=1 satisfies "
       "r+2z <= m-c_1 (equality) and 2r+8z <= 3m-c_2 (margin 2theta/(1+theta)) for 0<theta<=1/10")
falsifier("EX1", "theta = 0 (no positive margin c_1)", 0.0 if c1.subs(th, 0) > 0 else 1.0, 0.5)

# ---- EX2
CI = lambda s: s - sp.Rational(2, 3)
omega = D1/2
ok2 = CI(sp.Rational(11, 12)) == sp.Rational(1, 4)
ok2 &= sp.simplify((CI(sp.Rational(11, 12)) + omega) - (sp.Rational(1, 4) + omega)) == 0
ok2 &= sp.simplify(CI(sp.Rational(11, 12) + D1) - (sp.Rational(1, 4) + omega) - D1/2) == 0   # gap Delta_1/2
ok2 &= all((0 < omega.subs(D1, t) < t) for t in [sp.Rational(k, 1200) for k in range(1, 101)])
report("EX2", 0.0 if ok2 else 1.0, 0.0, "eq:d190: C_I(11/12) = 1/4, omega = Delta_1/2 in (0, Delta_1), and C_I(beta_*) - (1/4 + omega) "
       "= Delta_1/2")

# ---- EX3
E = sp.Max(M, 2*M - 1)
ok3 = E.subs(M, 1) == 1
# O/2 + max(M-O, 2M-O-1) <= max(M, 2M-1), O >= 0: branchwise
br1 = sp.simplify(M - (O/2 + M - O))               # = O/2 >= 0
br2 = sp.simplify((2*M - 1) - (O/2 + 2*M - O - 1))  # = O/2 >= 0
ok3 &= sp.simplify(br1 - O/2) == 0 and sp.simplify(br2 - O/2) == 0
# a grid check of the full max-inequality
for Mv in [sp.Rational(k, 10) for k in range(0, 31)]:
    for Ov in [sp.Rational(k, 10) for k in range(0, 31)]:
        lhs = Ov/2 + max(Mv - Ov, 2*Mv - Ov - 1); rhs = max(Mv, 2*Mv - 1)
        ok3 &= lhs <= rhs
report("EX3", 0.0 if ok3 else 1.0, 0.0, "eq:d198: eq:s1-cm-completed at M = 1 gives Z^{1+eps}; O/2+max(M-O,2M-O-1) <= max(M,2M-1) "
       "(symbolic on both branches, grid 31x31)")

# ---- EX4
H = 1 + th                                         # exponents of D
Y = H/6
e_main = sp.simplify(((2 + th) - Y)/2)
e_err = sp.simplify((2 - 2*Y)/2)
ok4 = sp.simplify(e_main - (sp.Rational(11, 12) + 5*th/12)) == 0
ok4 &= sp.simplify(e_err - (sp.Rational(5, 6) - th/6)) == 0
ok4 &= all(e_err.subs(th, t) <= sp.Rational(11, 12) for t in [sp.Rational(k, 1000) for k in range(1, 101)])
ok4 &= sp.limit(e_main, th, 0) == sp.Rational(11, 12)
report("EX4", 0.0 if ok4 else 1.0, 0.0, "eq:s6-bound: extraction exponents 11/12 + 5theta/12 and 5/6 - theta/6 <= 11/12 "
       "(H = D^{1+theta}, Y = H^{1/6}); limit theta -> 0 is 11/12")
Ybad = H/3
falsifier("EX4", "Y = H^{1/3} in place of H^{1/6}",
          float(abs(sp.simplify(((2 + th) - Ybad)/2 - (sp.Rational(11, 12) + 5*th/12)).subs(th, sp.Rational(1, 20)))), 1e-6)
