"""t2: planar Fourier transform of an angular profile (inserts.tex, eq. s5-hankel in the proof of Lemma s5-reflection).

With e(z) = exp(4 pi i Im z/sqrt3) and fhat(w) = int_C f(z) e(-zw) d^2z (the convention of [QF, Thm. 19]):
  g(z) = (z/|z|)^e W(|z|^2),  e = +-1   ==>
  ghat(w) = pi * e * (conj(w)/|w|)^e * int_0^oo W(x) J_1(4 pi |w| sqrt(x/3)) dx ,   ghat(0) = 0.
Tested by direct polar quadrature for an annular C_c^oo bump and for W(x) = x^{1/2} e^{-x}.
Falsifiers: prefactor -i instead of e; (w/|w|)^e instead of its conjugate; J_0 instead of J_1.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import special
from s5lib import report, falsifier, summary, SQ3

def bump(x):
    x = np.asarray(x, dtype=float); out = np.zeros_like(x)
    m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out

profiles = {
    "bump[1,2]": (bump, (1.0, np.sqrt(2.0)), (1.0, 2.0)),
    "x^1/2 e^-x": (lambda x: np.sqrt(x) * np.exp(-x), (0.0, 7.0), (0.0, 49.0)),
}

th = np.linspace(0, 2 * np.pi, 512, endpoint=False); dth = 2 * np.pi / 512
rg, rw = np.polynomial.legendre.leggauss(600)

def ghat_direct(W, rr, e, w):
    r0, r1 = rr
    r = 0.5 * (r1 - r0) * rg + 0.5 * (r1 + r0); wr = 0.5 * (r1 - r0) * rw
    Z = r[:, None] * np.exp(1j * th)[None, :]
    g = (Z / np.abs(Z)) ** e * W(r ** 2)[:, None]
    ph = np.exp(-4j * np.pi * (Z * w).imag / SQ3)
    return np.sum(g * ph * (r * wr)[:, None]) * dth

def hankel(W, xx, w, kernel=special.j1):
    x0, x1 = xx
    x = 0.5 * (x1 - x0) * rg + 0.5 * (x1 + x0); wx = 0.5 * (x1 - x0) * rw
    return np.sum(W(x) * kernel(4 * np.pi * abs(w) * np.sqrt(x / 3)) * wx)

ws = [0.0, 0.3 + 0.1j, -0.7 + 0.45j, 1.1 - 0.9j, 0.05j, 2.2 + 0.0j]
for name, (W, rr, xx) in profiles.items():
    for e in (1, -1):
        err = 0; f1 = f2 = f3 = 0; zero = None
        for w in ws:
            d = ghat_direct(W, rr, e, w)
            if w == 0:
                zero = abs(d); continue
            ang = (np.conj(w) / abs(w)) ** e
            H1 = hankel(W, xx, w)
            pred = np.pi * e * ang * H1
            err = max(err, abs(d - pred))
            f1 = max(f1, abs(d - np.pi * (-1j) * ang * H1))
            f2 = max(f2, abs(d - np.pi * e * (w / abs(w)) ** e * H1))
            f3 = max(f3, abs(d - np.pi * e * ang * hankel(W, xx, w, special.j0)))
        report(f"2-{name}-e{e:+d}", err, 1e-10, "ghat(w) = pi e (conj w/|w|)^e int W(x) J_1(4 pi |w| sqrt(x/3)) dx")
        report(f"2-{name}-e{e:+d}-zero", zero, 1e-12, "ghat(0) = int (z/|z|)^e W(|z|^2) d^2z = 0 (vanishing lattice main term)")
        falsifier(f"2-{name}-e{e:+d}", "prefactor -i", f1, 1e-4)
        falsifier(f"2-{name}-e{e:+d}", "unconjugated angle", f2, 1e-4)
        falsifier(f"2-{name}-e{e:+d}", "J_0 kernel", f3, 1e-4)
# the finite-order case for comparison: radial W has ghat(0) = pi * int W != 0
W = profiles["bump[1,2]"][0]
d0 = ghat_direct(lambda x: W(x), (1.0, np.sqrt(2.0)), 0, 0.0)
falsifier("2-zero", "radial (finite-order) profile has a nonzero zero frequency", abs(d0), 1e-3)
summary()
