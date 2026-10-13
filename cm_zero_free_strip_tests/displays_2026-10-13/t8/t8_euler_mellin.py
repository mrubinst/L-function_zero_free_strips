"""t8_euler_mellin.py -- eq:d196 of cm-eleven-twelfths-guide.tex (Appendix G, proof of Theorem thm:s6-1112):
   F_W(s) = What(s) sum_n mu_K(n) alpha(n) nu(n) N(n)^{-s}
          = What(s) prod_{pi not in S} (1 - conj(psi((pi))) N(pi)^{-s}) = What(s) / L^S(s, conj psi),
with What(s) = int_0^infty W(y) y^s dy/y, F_W(s) = int_0^infty A_1^nu(D) D^{-s} dD/D.

Setting (as in tests_cm/phase2_2026-10-12/s6): pi13 = 4 + 3 omega (norm 13), S = {2, lambda, pi13},
nu(n) = exp(2 pi i dlog_2(n mod pi13)/12) on primary n prime to S (order 12), 0 otherwise; eta = conj(nu);
psi((n)) = conj(alpha(n)) eta((n)) for primary n prime to S (eq:d191, eq:d192), so conj(psi((n))) = alpha(n) nu(n).

  T8.9   exact: the Euler product prod_{N(pi) <= B} (1 - alpha(pi) nu(pi) N(pi)^{-s}), expanded as a Dirichlet polynomial
         (ideals = primary generators), has coefficient mu_K(n) alpha(n) nu(n) at every primary n with N(n) <= B; and
         the Dirichlet convolution of mu alpha nu with conj(psi) = alpha nu is the identity (so the series is 1/L^S).
  T8.10  numerical, s = 3 + 2i and s = 2.2 - 5i: sum_{N(n)<=B} mu alpha nu N^{-s} against prod_{N(pi)<=B}(1 - conj psi(pi)N^{-s})
         and against 1 / sum_{N(n)<=B} conj(psi(n)) N^{-s}, with the tail bound sum_{N(n) > B} N(n)^{-Re s} as tolerance.
  T8.11  Mellin convention and the first equality: F_W(s) computed by quadrature in D of the exact finite sums A_1^nu(D)
         (W a complex bump on (1, 2)) equals What(s) sum_n mu alpha nu N^{-s} at s = 4 + 3i (relative error); and with
         W(y) = y^{-s0} phi(y), What(s0) = int_1^2 phi(y) dy/y > 0 (s0 = 0.95 + 14i).
Falsifiers (must fail): the Euler product with psi in place of conj(psi); What(-s) in place of What(s).
Usage: python3 t8_euler_mellin.py
"""
import os, sys, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import *
import numpy as np
TRAP = getattr(np, 'trapezoid', None) or np.trapz

PI13 = (4, 3); W0 = (-PI13[0]*pow(PI13[1], -1, 13)) % 13
DLOG = {}; x = 1
for k in range(12): DLOG[x] = k; x = (2*x) % 13
def red13(z): return (z[0] + z[1]*W0) % 13
def nu(n):
    if n == (1, 0): return 1.0 + 0j
    r = red13(n); return 0j if r == 0 else cmath.exp(2j*math.pi*DLOG[r]/12)
def al(n): return 1.0 + 0j if n == (1, 0) else alpha_c(n)
def muk(n): return 1 if n == (1, 0) else mu(n)
def cpsi(n): return al(n)*nu(n)                    # conj(psi((n))) = alpha(n) nu(n), n primary prime to S
def psi(n): return cpsi(n).conjugate()

# ------------------------------------------------------------------ T8.9 exact Dirichlet coefficients
B = 4000
ELS = [n for n in primaries(B) if red13(n) != 0]           # primary, prime to S (2, lambda, pi13)
IDX = {n: i for i, n in enumerate(ELS)}
PRS = [p for p in PR_ODD if N(p) <= B and red13(p) != 0]
coef = {(1, 0): 1.0 + 0j}                                  # expand prod (1 - cpsi(p) T^{N p}) on primary elements
coef_bad = {(1, 0): 1.0 + 0j}
for p in PRS:
    new = dict(coef); newb = dict(coef_bad)
    for n, c in coef.items():
        if N(n)*N(p) <= B:
            m = mul(n, p); new[m] = new.get(m, 0) - c*cpsi(p)
    for n, c in coef_bad.items():
        if N(n)*N(p) <= B:
            m = mul(n, p); newb[m] = newb.get(m, 0) - c*psi(p)
    coef, coef_bad = new, newb
err = 0.0; ferr = 0.0
for n in ELS:
    target = muk(n)*al(n)*nu(n)
    err = max(err, abs(coef.get(n, 0) - target)); ferr = max(ferr, abs(coef_bad.get(n, 0) - target))
# convolution with conj(psi): sum_{d | n} mu(d) cpsi(d) cpsi(n/d) = [n = 1]
cerr = 0.0
for n in ELS[:600]:
    s = sum(muk(d)*cpsi(d)*cpsi(div(n, d)) for d in ELS if N(d) <= N(n) and divides(d, n))
    cerr = max(cerr, abs(s - (1 if n == (1, 0) else 0)))
report("T8.9", max(err, cerr), 1e-12, f"eq:d196: Euler product prod (1 - conj psi(pi) N(pi)^-s) over {len(PRS)} primes has Dirichlet "
       f"coefficient mu alpha nu at all {len(ELS)} primary n prime to S with N(n) <= {B}; convolution with conj psi = delta (600 n)")
falsifier("T8.9", "Euler product with psi in place of conj(psi)", ferr, 1e-3)

# ------------------------------------------------------------------ T8.10 numerical values
NS = np.array([N(n) for n in ELS], float)
MUAN = np.array([muk(n)*al(n)*nu(n) for n in ELS])
CPS = np.array([cpsi(n) for n in ELS])
e10 = 0.0
for s in (3 + 2j, 2.2 - 5j):
    dser = np.sum(MUAN*NS**(-s))
    eul = np.prod([1 - cpsi(p)*N(p)**(-s) for p in PRS])
    Ls = np.sum(CPS*NS**(-s))
    tail = 0.5*sum(1.0/(m**s.real) for m in range(B + 1, 400*B)) + 0.5*(400*B)**(1 - s.real)/(s.real - 1)
    rel = max(abs(dser - eul), abs(dser - 1/Ls))/max(tail, 1e-300)
    e10 = max(e10, rel)
    print(f"      s={s}: sum mu alpha nu N^-s = {dser:.12f}, Euler product = {eul:.12f}, 1/L^S = {1/Ls:.12f}, tail bound {tail:.2e}")
report("T8.10", e10, 3.0, "eq:d196 numerically at s = 3+2i, 2.2-5i: Dirichlet series, Euler product and 1/L^S(s, conj psi) "
       "agree within a few times the tail bound (err = discrepancy/tail bound)")

# ------------------------------------------------------------------ T8.11 Mellin convention
def Wc(y): return np.where((y > 1) & (y < 2), np.exp(-1/np.clip((y - 1)*(2 - y), 1e-300, None))*(1 + 0.5j*y), 0)
def mellin(Wf, s, a=1.0, b=2.0, M=200001):
    y = np.linspace(a, b, M); f = Wf(y)*y**(s - 1)
    return TRAP(f, y)
s = 4 + 3j
What = mellin(Wc, s); What_bad = mellin(Wc, -s)
sel = NS <= 1000
# F_W(s) = int A_1(D) D^{-s} dD/D, A_1(D) = sum mu alpha nu W(N/D); W(N/D) != 0 for N/2 < D < N; D in (1/2, 1000)
t = np.linspace(math.log(0.5), math.log(1000.0), 60001); Dg = np.exp(t)
A1 = np.zeros(len(Dg), complex)
for c, Nn in zip(MUAN[sel], NS[sel]):
    A1 += c*Wc(Nn/Dg)
FW = TRAP(A1*Dg**(-s), t)
rhs = What*np.sum(MUAN[sel]*NS[sel]**(-s)); rhs_bad = What_bad*np.sum(MUAN[sel]*NS[sel]**(-s))
e11 = abs(FW - rhs)/abs(rhs); f11 = abs(FW - rhs_bad)/abs(rhs)
# positivity at s0 with W(y) = y^{-s0} phi(y)
s0 = 0.95 + 14j
phi_ = lambda y: np.where((y > 1) & (y < 2), np.exp(-1/np.clip((y - 1)*(2 - y), 1e-300, None)), 0.0)
W_s0 = lambda y: y**(-s0)*phi_(y)
val = mellin(W_s0, s0); target = mellin(lambda y: phi_(y) + 0j, 0.0)
e11b = abs(val - target)
report("T8.11", max(e11, e11b), 1e-8, f"eq:d196 first equality: F_W(s) by quadrature of A_1^nu(D) (columns N(n) <= 1000) = "
       f"What(s) sum mu alpha nu N^-s at s = 4+3i (rel err {e11:.2e}); What(s0) = int phi dy/y = {target.real:.6f} > 0 for "
       f"W = y^-s0 phi (err {e11b:.1e})")
falsifier("T8.11", "What(-s) in place of What(s)", f11, 1e-3)
