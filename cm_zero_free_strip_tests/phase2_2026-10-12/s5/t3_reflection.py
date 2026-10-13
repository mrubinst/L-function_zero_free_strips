"""t3: the CM reflection of plain sums (inserts.tex, Lemma s5-reflection).

For a primitive Hecke character xi of conductor f and infinity type (a/|a|)^e (finite part phi primitive mod f,
phi(u) = u^{-e}), and W with (z/|z|)^e W(|z|^2) Schwartz,
   X^{-1/2} sum_a xi(a) W(Na/X) = eps(xi) Y^{-1/2} sum_b conj(xi(b)) W^sharp(Nb/Y),
   Y = C/X,  C = 3 N(f)/(4 pi^2),  W^sharp(y) = int_0^oo W(x) J_1(2 sqrt(xy)) dx,  eps(xi) = e gamma(phi) (f/|f|)^e.
Sums over nonzero ideals; computed as (1/6) x sums over nonzero elements.
Characters: alpha, conj(alpha), alpha chi_p^5 (N p = 7, conductor p), alpha^e chi_p^j (N p = 13, 19), a product of two
residue characters, and alpha times the cubic character mod 2.  Profiles: x^{1/2}e^{-x}, x^{3/2}e^{-x} (closed-form
W^sharp) and a log-normal profile (W^sharp by Mellin-Barnes).
Falsifiers: (i) the conductor scale of the finite-order part (no change at lambda), (ii) OAI78's finite-order dual
profile (J_0 kernel), (iii) root number without the factor e, (iv) root number without (f/|f|)^e.
"""
import sys, os, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from scipy import special
from s5lib import (report, falsifier, summary, chi_p, primes_upto, cx, N, mul, conj, OM, UNITS)
from hecke import from_primary_function, primitive, gauss, elements, Char

LAM = (1, 2)
TWO = (-2, 0)
alpha = lambda a: cx(a) / abs(cx(a))

P = primes_upto(40)
p7 = [p for p in P if N(p) == 7][0]
p13 = [p for p in P if N(p) == 13][0]
p19 = [p for p in P if N(p) == 19][0]

def cubic_mod2(a):
    # (a/2)_3 = a mod 2 in F_4^x = {1, w, w^2}
    r = ((a[0] % 2), (a[1] % 2))
    return {(1, 0): 1, (0, 1): OM, (1, 1): OM ** 2}[r]

cases = []
cases.append(("alpha", lambda a: alpha(a), (3, 0), [LAM, LAM], 1))
cases.append(("conj alpha", lambda a: 1 / alpha(a), (3, 0), [LAM, LAM], -1))
cases.append(("alpha chi_p7^5", lambda a: alpha(a) * chi_p(a, p7) ** 5, mul((3, 0), p7), [LAM, LAM, p7], 1))
for e in (1, -1):
    for j in (1, 2, 3):
        cases.append((f"alpha^{e} chi_p13^{j}", (lambda e, j: lambda a: alpha(a) ** e * chi_p(a, p13) ** j)(e, j),
                      mul((3, 0), p13), [LAM, LAM, p13], e))
cases.append(("alpha chi_p7 chi_p19^3", lambda a: alpha(a) * chi_p(a, p7) * chi_p(a, p19) ** 3,
              mul(mul((3, 0), p7), p19), [LAM, LAM, p7, p19], 1))
cases.append(("conj alpha (a/2)_3", lambda a: cubic_mod2(a) / alpha(a), (6, 0), [LAM, LAM, TWO], -1))
cases.append(("alpha^-1 chi_p7^2 (a/2)_3", lambda a: chi_p(a, p7) ** 2 * cubic_mod2(a) / alpha(a),
              mul((6, 0), p7), [LAM, LAM, TWO, p7], -1))

# ------------------------------------------------------------------ profiles
def W1(x): return np.sqrt(x) * np.exp(-x)
def W1s(y): return np.sqrt(y) * np.exp(-y)
def W2(x): return x ** 1.5 * np.exp(-x)
def W2s(y): return np.sqrt(y) * np.exp(-y) * (2 - y)
def W1s_J0(y):   # OAI78's finite-order dual of W1: int x^{1/2} e^{-x} J_0(2 sqrt(xy)) dx = Gamma(3/2) 1F1(3/2;1;-y)
    return np.array([float(mp.gamma(1.5) * mp.hyp1f1(1.5, 1, -t)) for t in np.atleast_1d(y)])

SIG = 0.6
def W3(x): return np.exp(-np.log(x) ** 2 / (2 * SIG ** 2))
_tt = np.linspace(-30, 30, 4001)
def W3s(y):
    # Mellin-Barnes on Re s = 0.3: W^sharp(y) = (1/2pi) int MW(1-s) G(s) y^{-s} dt, MW(s) = sqrt(2pi) SIG e^{SIG^2 s^2/2}
    s = 0.3 + 1j * _tt
    MW1 = math.sqrt(2 * math.pi) * SIG * np.exp(SIG ** 2 * (1 - s) ** 2 / 2)
    G = np.exp(special.loggamma(s + 0.5) - special.loggamma(1.5 - s))
    base = MW1 * G
    y = np.atleast_1d(y).astype(float)
    out = np.empty(len(y))
    for i in range(0, len(y), 200):
        yy = y[i:i + 200]
        out[i:i + 200] = np.real((base[None, :] * np.exp(-np.outer(np.log(yy), s))).sum(axis=1)) * (_tt[1] - _tt[0]) / (2 * np.pi)
    return out

profiles = [("x^1/2 e^-x", W1, W1s, 60.0, 60.0), ("x^3/2 e^-x", W2, W2s, 70.0, 70.0), ("lognormal", W3, W3s, 200.0, 1500.0)]

def plain_sum(ch, Wfun, X, cut, conjugate=False):
    A, B, Z, n = elements(int(cut * X) + 2)
    vals = ch(A, B) * (Z / np.abs(Z)) ** ch.e
    if conjugate:
        vals = np.conj(vals)
    return np.sum(vals * Wfun(n / X)) / 6.0 / math.sqrt(X)

def plain_sum_dual(ch, Ws, Y, cut):
    A, B, Z, n = elements(int(cut * Y) + 2)
    ns = np.unique(n)
    vals = np.conj(ch(A, B) * (Z / np.abs(Z)) ** ch.e)
    # group by norm to evaluate W^sharp once per norm
    idx = np.searchsorted(ns, n)
    wsv = Ws(ns / Y)
    return np.sum(vals * wsv[idx]) / 6.0 / math.sqrt(Y)

for name, xi, m, primes, e in cases:
    ch0 = from_primary_function(xi, m, primes, e)
    ch = primitive(ch0)
    f = ch.m; Nf = N(f)
    g = gauss(ch)
    eps = e * g * (cx(f) / abs(cx(f))) ** e
    C = 3 * Nf / (4 * math.pi ** 2)
    # conductor of the finite-order part: phi * (u-character)^{-e}... computed as the primitive version of the
    # finite-order character xi(a)/alpha(a)^e (trivial infinity type)
    chfin = primitive(from_primary_function(lambda a: xi(a) / alpha(a) ** e, m, primes, 0))
    Cfin = 3 * N(chfin.m) / (4 * math.pi ** 2)
    print(f"--- {name}: e={e}, conductor N(f)={Nf} (finite part {N(chfin.m)}), |gamma|={abs(g):.15f}, eps={eps:.6f}")
    for pname, W, Ws, cut, dcut in profiles:
        if pname == "lognormal" and name not in ("alpha", "alpha chi_p7^5", "alpha^-1 chi_p13^2"):
            continue
        worst = 0; f1 = f2 = f3 = f4 = 0
        for X in [0.6, 1.5, 4.0]:
            lhs = plain_sum(ch, W, X, cut)
            Y = C / X
            rhs0 = plain_sum_dual(ch, Ws, Y, dcut)
            worst = max(worst, abs(lhs - eps * rhs0))
            if Cfin != C:
                Yf = Cfin / X
                f1 = max(f1, abs(lhs - eps * plain_sum_dual(ch, Ws, Yf, dcut)))
            else:
                f1 = 1.0   # not applicable
            if pname == "x^1/2 e^-x":
                f2 = max(f2, abs(lhs - eps * plain_sum_dual(ch, W1s_J0, Y, 60.0)))
            f3 = max(f3, abs(lhs - g * (cx(f) / abs(cx(f))) ** e * rhs0)) if e == -1 else 1.0
            f4 = max(f4, abs(lhs - e * g * rhs0)) if abs((cx(f) / abs(cx(f))) ** e - 1) > 1e-6 else 1.0
        report(f"3-{name}-{pname}", worst, 1e-9, f"reflection identity, X in (0.6,1.5,4), C = 3N(f)/(4pi^2) = {C:.5f}")
        if Cfin != C:
            falsifier(f"3-{name}-{pname}", "finite-part conductor scale", f1, 1e-4)
        if pname == "x^1/2 e^-x":
            falsifier(f"3-{name}-{pname}", "J_0 (finite-order) dual profile", f2, 1e-4)
        if e == -1:
            falsifier(f"3-{name}-{pname}", "root number without the factor e", f3, 1e-4)
        if abs((cx(f) / abs(cx(f))) ** e - 1) > 1e-6:
            falsifier(f"3-{name}-{pname}", "root number without (f/|f|)^e", f4, 1e-4)
    report(f"3-{name}-|gamma|", abs(abs(g) - 1), 1e-12, "|eps| = |gamma(phi)| = 1")
    def v3(n):
        v = 0
        while n % 3 == 0: n //= 3; v += 1
        return v, n
    vf, rf = v3(Nf); vg, rg = v3(N(chfin.m))
    bad = (rf != rg) or (vf > max(2, vg)) or (vg > max(2, vf))
    report(f"3-{name}-conductor", 1.0 if bad else 0.0, 0.5,
           f"conductor changes only at lambda: v_lambda(f)={vf}, v_lambda(g)={vg}, prime-to-3 parts {rf}={rg}")
summary()
