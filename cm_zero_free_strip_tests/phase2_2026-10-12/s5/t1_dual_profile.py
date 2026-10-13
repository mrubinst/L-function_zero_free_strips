"""t1: the CM dual profile W^sharp (inserts.tex, Lemma s5-sharp and eq. s5-sharp-def).

Claims tested
 1a  int_0^oo J_1(2 sqrt y) y^{s-1} dy = Gamma(s+1/2)/Gamma(3/2-s)      (-1/2 < Re s < 3/4)
     control: int_0^oo J_0(2 sqrt y) y^{s-1} dy = Gamma(s)/Gamma(1-s)   (OAI78's finite-order kernel)
 1b  closed forms: W(x)=x^{1/2}e^{-x} -> W^sharp(y)=y^{1/2}e^{-y};  W(x)=x^{3/2}e^{-x} -> y^{1/2}e^{-y}(2-y)
 1c  Mellin: M W^sharp(s) = M W(1-s) Gamma(s+1/2)/Gamma(3/2-s) for an annular bump W, at several s
 1d  small y: W^sharp(y) = M W(3/2) y^{1/2} + O(y^{3/2}); central-normalized annular coefficients O(2^{-j})
     (finite-order kernel gives only O(2^{-j/2}))
 1e  Stirling: |Gamma(s+1/2)/Gamma(3/2-s)| ~ |t|^{2 sigma - 1}, same exponent as |Gamma(s)/Gamma(1-s)|
Falsifiers: wrong gamma quotient in 1a/1c, wrong small-y exponent in 1d, exponent 2 sigma in 1e.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import numpy as np
from scipy import special, integrate
from s5lib import report, falsifier, summary

mp.mp.dps = 30

def G_cm(s): return mp.gamma(s + 0.5) / mp.gamma(1.5 - s)
def G_fin(s): return mp.gamma(s) / mp.gamma(1 - s)

# ---------------------------------------------------------------- 1a Bessel-Mellin
def mellin_bessel(nu, s):
    # int_0^oo J_nu(2 sqrt y) y^{s-1} dy = 2^{1-2s} int_0^oo J_nu(x) x^{2s-1} dx  (y = x^2/4).
    # On [0,1] subtract the leading Taylor term (x/2)^nu/nu! and integrate it exactly; on [1,oo) use quadosc.
    lead = lambda x: (x / 2) ** nu / mp.factorial(nu)
    a = mp.quad(lambda x: (mp.besselj(nu, x) - lead(x)) * x ** (2 * s - 1), [0, 1])
    a += 2 ** (-nu) / mp.factorial(nu) / (2 * s + nu)
    b = mp.quadosc(lambda x: mp.besselj(nu, x) * x ** (2 * s - 1), [1, mp.inf], period=2 * mp.pi)
    return 2 ** (1 - 2 * s) * (a + b)

worst = 0; worst_f = 0; fals = 1e9
for s in [mp.mpf('0.1'), mp.mpc('0.2', '0.7'), mp.mpc('-0.3', '1.5'), mp.mpc('0.0', '-2.2'), mp.mpc('-0.4', '0.3')]:
    v1 = mellin_bessel(1, s)
    worst = max(worst, abs(v1 - G_cm(s)) / abs(G_cm(s)))
    fals = min(fals, abs(v1 - G_fin(s)) / abs(v1))
    if mp.re(s) > 0:
        v0 = mellin_bessel(0, s)
        worst_f = max(worst_f, abs(v0 - G_fin(s)) / abs(G_fin(s)))
report("1a", float(worst), 1e-12, "J_1 kernel Mellin = Gamma(s+1/2)/Gamma(3/2-s) at 5 points")
report("1a-control", float(worst_f), 1e-12, "J_0 kernel Mellin = Gamma(s)/Gamma(1-s) (OAI78 finite-order case)")
falsifier("1a", "finite-order quotient for the J_1 kernel", float(fals), 1e-3)

# ---------------------------------------------------------------- 1b closed forms
def sharp_num(W, y):
    return mp.quad(lambda x: W(x) * mp.besselj(1, 2 * mp.sqrt(x * y)), [0, 1, 4, 10, 30, mp.inf])

err = 0; errw = 0
for y in [mp.mpf('0.03'), mp.mpf('0.5'), mp.mpf('1.7'), mp.mpf('4.2')]:
    a = sharp_num(lambda x: mp.sqrt(x) * mp.e ** (-x), y); b = mp.sqrt(y) * mp.e ** (-y)
    c = sharp_num(lambda x: x ** 1.5 * mp.e ** (-x), y); d = mp.sqrt(y) * mp.e ** (-y) * (2 - y)
    err = max(err, abs(a - b), abs(c - d))
    # wrong: J_0 kernel for the same W
    a0 = mp.quad(lambda x: mp.sqrt(x) * mp.e ** (-x) * mp.besselj(0, 2 * mp.sqrt(x * y)), [0, 1, 4, 10, 30, mp.inf])
    errw = max(errw, abs(a0 - b))
report("1b", float(err), 1e-20, "x^{1/2}e^{-x} -> y^{1/2}e^{-y}; x^{3/2}e^{-x} -> y^{1/2}e^{-y}(2-y)")
falsifier("1b", "J_0 kernel instead of J_1", float(errw), 1e-3)

# ---------------------------------------------------------------- 1c Mellin of W^sharp for an annular bump
def bump(x):
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    m = (x > 1) & (x < 2)
    u = (2 * x[m] - 3)          # in (-1,1)
    out[m] = np.exp(-1.0 / (1 - u * u))
    return out

xg, xw = np.polynomial.legendre.leggauss(400)
xs = 1.5 + 0.5 * xg; ws = 0.5 * xw
Wx = bump(xs)

def sharp_bump(y, kernel=special.j1):
    y = np.atleast_1d(y)
    return (Wx[None, :] * kernel(2 * np.sqrt(np.outer(y, xs)))) @ ws

def MW_bump(s):
    return complex(np.sum(Wx * xs ** (s - 1) * ws))

# Mellin inversion: W^sharp(y) = (1/2 pi) int M W(1-s) G(s) y^{-s} dt on Re s = sigma, for several sigma
# (no poles of G(s) = Gamma(s+1/2)/Gamma(3/2-s) in Re s > -1/2, so the result is independent of sigma).
xg2, xw2 = np.polynomial.legendre.leggauss(2000)
xs2 = 1.5 + 0.5 * xg2; ws2 = 0.5 * xw2; Wx2 = bump(xs2)
tt = np.arange(-1600, 1600.0001, 0.05)
def MWvec(svec):
    out = np.empty(len(svec), dtype=complex)
    for i in range(0, len(svec), 4000):
        sv = svec[i:i + 4000]
        out[i:i + 4000] = (xs2[None, :] ** (sv[:, None] - 1) * (Wx2 * ws2)[None, :]).sum(axis=1)
    return out
def gamma_ratio(svec, fin=False):
    sv = svec
    if fin:
        return np.exp(special.loggamma(sv) - special.loggamma(1 - sv))
    return np.exp(special.loggamma(sv + 0.5) - special.loggamma(1.5 - sv))
ytest = np.array([0.05, 0.4, 1.3, 3.7, 9.0])
direct = sharp_bump(ytest)
err = 0; fals = 0
for sigma in [-0.3, 0.4, 1.0]:
    s = sigma + 1j * tt
    MW1 = MWvec(1 - s)
    for fin in (False, True):
        integrand = MW1 * gamma_ratio(s, fin)
        vals = np.array([np.real(np.sum(integrand * y ** (-s)) * 0.05 / (2 * np.pi)) for y in ytest])
        if fin:
            if sigma > 0:
                fals = max(fals, np.abs(vals - direct).max())
        else:
            err = max(err, np.abs(vals - direct).max())
report("1c", err, 1e-8, "Mellin inversion of M W(1-s)Gamma(s+1/2)/Gamma(3/2-s) on Re s = -0.3, 0.4, 1.0 reproduces the Hankel transform of a C_c^oo bump")
falsifier("1c", "finite-order gamma quotient", fals, 1e-3)

# ---------------------------------------------------------------- 1d small-y behaviour and annular coefficients
yy = np.array([1e-2, 1e-3, 1e-4, 1e-5])
r = sharp_bump(yy) / np.sqrt(yy)
target = MW_bump(1.5).real
rel = np.abs(r - target) / target
# O(y^{3/2}) remainder: (W^sharp - MW(3/2) y^{1/2}) / y^{3/2} bounded
rem = (sharp_bump(yy) - target * np.sqrt(yy)) / yy ** 1.5
report("1d-leading", float(rel[-1]), 1e-4, f"W^sharp(y)/sqrt(y) -> MW(3/2)={target:.6g}; rel err {rel}")
report("1d-remainder", float(np.ptp(rem) / abs(rem).max()), 0.05, f"(W^sharp - MW(3/2)y^1/2)/y^3/2 = {rem}")
# central normalization on annulus [2^{-j-1}, 2^{-j}] of the main scale: coefficient (2^{-j})^{1/2} sup |W^sharp|
coef_cm = []; coef_fin = []
for j in range(4, 16):
    yj = np.linspace(2.0 ** (-j - 1), 2.0 ** (-j), 50)
    coef_cm.append(2.0 ** (-j / 2) * np.abs(sharp_bump(yj)).max() * 2.0 ** j)
    coef_fin.append(2.0 ** (-j / 2) * np.abs(sharp_bump(yj, special.j0)).max() * 2.0 ** j)
coef_cm = np.array(coef_cm); coef_fin = np.array(coef_fin)
report("1d-annuli", float(coef_cm.max() / coef_cm.min() - 1), 0.5,
       "2^j x (annular coefficient) is bounded for the CM kernel: coefficients are O(2^-j)")
falsifier("1d-annuli", "the same O(2^-j) claim for OAI78's J_0 kernel", float(coef_fin.max() / coef_fin.min()), 10.0)

# ---------------------------------------------------------------- 1e Stirling exponents
err = 0; fals = 0
for sig in [-0.25, 0.0, 0.5, 1.0, 3.0]:
    for t in [1e3, 1e4, 1e5]:
        s = mp.mpc(sig, t)
        rc = abs(G_cm(s)) / t ** (2 * sig - 1)
        rf = abs(G_fin(s)) / t ** (2 * sig - 1)
        err = max(err, abs(rc - 1), abs(rf - 1))
        fals = max(fals, abs(abs(G_cm(s)) / t ** (2 * sig) - 1))
report("1e", float(err), 1e-2, "|G(sigma+it)| / t^{2 sigma -1} -> 1 for both quotients (t = 1e3..1e5)")
falsifier("1e", "exponent 2 sigma", float(fals), 0.5)

summary()
