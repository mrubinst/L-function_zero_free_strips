"""c2: the planar Fourier transforms eq:d151 (C.7) and eq:d152 (C.8) of Appendix C (proof of Lemma s2-fe(i)).

Conventions: e(z) = exp(4 pi i Im z / sqrt 3); for w = x + iy and eta = a + ib, Im(w eta) = x b + y a.

Claims
  T7-C7  eq:d151  int_C exp(-2 pi v |w|^2/sqrt3) e(-w eta) dx dy = sqrt3/(2v) exp(-2 pi |eta|^2/(sqrt3 v)):
         (i) symbolic (sympy: the integral factors into two Gaussian integrals, evaluated in closed form);
         (ii) numerically by 1-D quadratures (mpmath) at 5 values (v, eta).
         Falsifier: the convention e(z) = exp(2 pi i Im z) gives a different transform (numerically).
  T7-C8a the Wirtinger step: d/d(conj eta) e(-w eta) = (2 pi/sqrt3) conj(w) e(-w eta)  (symbolic).
  T7-C8b eq:d152  Ghat(eta) = int conj(w) exp(-2 pi v|w|^2/sqrt3) e(-w eta) dx dy
             = (sqrt3/(2 pi)) d/d(conj eta)[sqrt3/(2v) exp(-2 pi eta conj(eta)/(sqrt3 v))]
             = -sqrt3/(2 v^2) eta exp(-2 pi |eta|^2/(sqrt3 v)):
         the second equality symbolically; the first numerically by quadrature at 5 values (v, eta).
         Falsifiers: d/d(eta) in place of d/d(conj eta) (gives conj(eta)); the sign +.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t7lib import report, finish
import sympy as sp
import mpmath as mp

t0 = time.time()
x, y, a, b = sp.symbols("x y a b", real=True)
v = sp.symbols("v", positive=True)
c = 2 * sp.pi * v / sp.sqrt(3)
k = 4 * sp.pi / sp.sqrt(3)
# e(-w eta) = exp(-i k Im(w eta)) = exp(-i k (x b + y a))
Ix = sp.integrate(sp.exp(-c * x**2 - sp.I * k * b * x), (x, -sp.oo, sp.oo))
Iy = sp.integrate(sp.exp(-c * y**2 - sp.I * k * a * y), (y, -sp.oo, sp.oo))
lhs = sp.simplify(Ix * Iy)
rhs = sp.sqrt(3) / (2 * v) * sp.exp(-2 * sp.pi * (a**2 + b**2) / (sp.sqrt(3) * v))
sym1 = sp.simplify(lhs / rhs) == 1


def Iquad(vv, eta, mom_x=0, mom_y=0, kk=None):
    """int x^mom_x y^mom_y exp(-2 pi v (x^2+y^2)/sqrt3) exp(-i kk Im(w eta)) dx dy via two 1-D quadratures."""
    kk = 4 * mp.pi / mp.sqrt(3) if kk is None else kk
    cc = 2 * mp.pi * vv / mp.sqrt(3)
    A, B = mp.re(eta), mp.im(eta)
    fx = mp.quad(lambda t: t**mom_x * mp.exp(-cc * t * t - 1j * kk * B * t), [-mp.inf, 0, mp.inf])
    fy = mp.quad(lambda t: t**mom_y * mp.exp(-cc * t * t - 1j * kk * A * t), [-mp.inf, 0, mp.inf])
    return fx * fy


mp.mp.dps = 25
samples = [(mp.mpf("0.37"), mp.mpc("0.3", "-0.8")), (mp.mpf(1), mp.mpc("1.1", "0.4")), (mp.mpf("2.5"), mp.mpc("-0.6", "1.7")),
           (mp.mpf("0.9"), mp.mpc(0, 0)), (mp.mpf("1.6"), mp.mpc("2", "-1"))]
err1 = 0; fal1 = 1e9
for vv, eta in samples:
    num = Iquad(vv, eta)
    ex = mp.sqrt(3) / (2 * vv) * mp.exp(-2 * mp.pi * abs(eta)**2 / (mp.sqrt(3) * vv))
    err1 = max(err1, abs(num - ex) / abs(ex))
    if eta != 0:
        fal1 = min(fal1, abs(Iquad(vv, eta, kk=2 * mp.pi) - ex) / abs(ex))
report("T7-C7", "eq:d151", "Fourier transform of the Gaussian with e(z) = exp(4 pi i Im z/sqrt3)", sym1 and err1 < 1e-15,
       f"symbolic (sympy) equality: {sym1}; quadrature at 5 (v, eta): max rel. error {mp.nstr(err1, 3)}; falsifier "
       f"e(z)=exp(2 pi i Im z): min rel. error {mp.nstr(fal1, 3)}", falsifier_failed=fal1 > 1e-2)

# Wirtinger derivative d/d(conj eta) = (d/da + i d/db)/2
wx, wy = sp.symbols("wx wy", real=True)
E = sp.exp(-sp.I * k * (wx * b + wy * a))
dbar = lambda F: (sp.diff(F, a) + sp.I * sp.diff(F, b)) / 2
dhol = lambda F: (sp.diff(F, a) - sp.I * sp.diff(F, b)) / 2
sym2 = sp.simplify(dbar(E) - 2 * sp.pi / sp.sqrt(3) * (wx - sp.I * wy) * E) == 0
report("T7-C8a", "proof of lem:s2-fe(i)", "d/d(conj eta) e(-w eta) = (2 pi/sqrt3) conj(w) e(-w eta)", sym2,
       f"symbolic: {sym2}")

etaexpr = a + sp.I * b
G = sp.sqrt(3) / (2 * v) * sp.exp(-2 * sp.pi * (a**2 + b**2) / (sp.sqrt(3) * v))
mid = sp.sqrt(3) / (2 * sp.pi) * dbar(G)
right = -sp.sqrt(3) / (2 * v**2) * etaexpr * sp.exp(-2 * sp.pi * (a**2 + b**2) / (sp.sqrt(3) * v))
sym3 = sp.simplify(mid - right) == 0
sym3f = sp.simplify(sp.sqrt(3) / (2 * sp.pi) * dhol(G) - right) == 0        # falsifier d/d eta
err2 = 0; fal2 = 1e9; fal3 = 1e9
for vv, eta in samples:
    # conj(w) = x - i y
    num = Iquad(vv, eta, mom_x=1) - 1j * Iquad(vv, eta, mom_y=1)
    ex = -mp.sqrt(3) / (2 * vv**2) * eta * mp.exp(-2 * mp.pi * abs(eta)**2 / (mp.sqrt(3) * vv))
    if eta == 0:
        err2 = max(err2, abs(num))
        continue
    err2 = max(err2, abs(num - ex) / abs(ex))
    fal2 = min(fal2, abs(num - (-mp.sqrt(3) / (2 * vv**2) * mp.conj(eta) * mp.exp(-2 * mp.pi * abs(eta)**2 / (mp.sqrt(3) * vv)))) / abs(ex))
    fal3 = min(fal3, abs(num + ex) / abs(ex))
report("T7-C8b", "eq:d152", "Ghat(eta) = (sqrt3/2pi) dbar(...) = -sqrt3/(2v^2) eta exp(-2 pi|eta|^2/(sqrt3 v))",
       sym3 and err2 < 1e-15,
       f"second equality symbolic: {sym3}; quadrature of int conj(w) e^(-2 pi v|w|^2/sqrt3) e(-w eta) at 5 (v, eta): max error "
       f"{mp.nstr(err2, 3)} (relative; absolute at eta = 0); falsifiers: d/d eta symbolic equality {sym3f}, conj(eta) min rel. error {mp.nstr(fal2, 3)}, "
       f"sign + min rel. error {mp.nstr(fal3, 3)}", falsifier_failed=(not sym3f and fal2 > 1e-2 and fal3 > 1e-2))
print(f"time {time.time() - t0:.1f}s")
finish()
