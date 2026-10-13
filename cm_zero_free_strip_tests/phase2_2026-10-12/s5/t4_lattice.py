"""t4: masked lattice sums with an angular factor (inserts.tex, Lemma s5-lattice; CM form of OAI78
Lemma lem:centered-lattice-cancellation, L14545-14680).

   L(X,t) = sum_{l outside S} alpha(l)^e theta(l) 1_{(l,R)=1} N(l)^{it} W(N(l)/X),   S = {lambda, 2},
l over primary generators, theta a fixed ray character with conductor in S, W a C_c^oo bump on [1,2].
Claims: (a) for e = +-1 the zero-frequency (volume) term is absent: |L(X,t)| stays bounded (indeed tends to 0 for the
unmasked sum) as X grows, uniformly in the mask; (b) the masked sum equals the inclusion-exclusion
sum_{d | R} mu(d) alpha(d)^e theta(d) N(d)^{it} L^0(X/N(d), t) exactly; (c) the centered and uncentered products obey
(2.18f) of OAI78 a fortiori.
Falsifier: e = 0 (OAI78's finite-order case) has the volume term c_theta X I(t): |L|/X tends to
c = pi/(6 sqrt3) * |int W(y) y^{it} dy| for theta = 1.
"""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from s5lib import report, falsifier, summary, OM, SQ3, primes_upto, N as Nrm, mul, cx, primary_assoc

def bump(x):
    out = np.zeros_like(x); m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out

def I_t(t):
    y = np.linspace(1, 2, 200001)
    return np.trapezoid(bump(y) * y ** (1j * t), y)

def primary_elements(Nmax):
    R = int(math.ceil(2 * math.sqrt(Nmax / 3.0))) + 2
    a = np.arange(-R, R + 1)
    A, B = np.meshgrid(a, a, indexing='ij'); A = A.ravel(); B = B.ravel()
    keep = ((A - 1) % 3 == 0) & (B % 3 == 0) & ~((A % 2 == 0) & (B % 2 == 0))
    A = A[keep]; B = B[keep]; n = A * A - A * B + B * B
    k2 = (n > 0) & (n <= Nmax)
    return A[k2], B[k2], n[k2]

def theta_val(name, A, B, n):
    if name == "1":
        return np.ones(len(A))
    if name == "chi(-1)":                      # chi_l(-1) = +1 iff N(l) = 1 mod 4 (l coprime to 6)
        return np.where(n % 4 == 1, 1.0, -1.0)
    if name == "(l/2)_3":                       # cubic character mod 2: l mod 2 in F_4^x = {1, w, w^2}
        a = A % 2; b = B % 2
        return np.where((a == 1) & (b == 0), 1.0 + 0j, np.where((a == 0) & (b == 1), OM, OM ** 2))
    if name == "chi(-1)(l/2)_3":
        return theta_val("chi(-1)", A, B, n) * theta_val("(l/2)_3", A, B, n)

def divisible(A, B, p):
    c, d = p; q = Nrm(p)
    # p | l  iff  l * conj(p) = 0 mod N(p) coordinatewise; conj(p) = (c - d) - d w;
    # exact product (a + b w)(x + y w) = (a x - b y) + (a y + b x - b y) w with x = c - d, y = -d
    x, y = c - d, -d
    P = A * x - B * y; Q = A * y + B * x - B * y
    return (P % q == 0) & (Q % q == 0)

def Lsum(e, th, t, X, A, B, n, mask_primes=()):
    sel = (n <= 2 * X) & (n >= X)
    a, b, nn = A[sel], B[sel], n[sel]
    Z = a + b * OM
    v = (Z / np.abs(Z)) ** e * theta_val(th, a, b, nn) * nn.astype(float) ** (1j * t) * bump(nn / X)
    for p in mask_primes:
        v = v * ~divisible(a, b, p)
    return v.sum()

Xs = [30.0, 300.0, 3000.0, 3e4, 3e5]
A, B, n = primary_elements(int(2 * Xs[-1]) + 2)
dens = math.pi / (6 * SQ3)

for th in ["1", "chi(-1)", "(l/2)_3", "chi(-1)(l/2)_3"]:
    for t in [0.0, 3.7]:
        It = I_t(t)
        for e in (1, -1):
            vals = [abs(Lsum(e, th, t, X, A, B, n)) for X in Xs]
            # bounded uniformly in X (the constant depends on the fixed modulus of theta), and decaying for large X
            report(f"4a-theta={th}-t={t}-e={e:+d}", max(max(vals) / 10.0, vals[-1] / 0.05), 1.0,
                   "unmasked |L(X,t)| for X=" + ",".join(f"{X:g}" for X in Xs) + ": " + ", ".join(f"{v:.2e}" for v in vals))
        if th == "1":
            v0 = [abs(Lsum(0, th, t, X, A, B, n)) / X for X in Xs]
            rel = abs(v0[-1] - dens * abs(It)) / (dens * abs(It))
            report(f"4a-falsifier-calibration-t={t}", rel, 1e-2,
                   f"e=0: |L|/X -> pi/(6sqrt3)|I(t)| = {dens*abs(It):.5f}; got " + ", ".join(f"{v:.5f}" for v in v0))
            falsifier(f"4a-theta=1-t={t}", "no angular factor (e=0): volume main term", abs(Lsum(0, th, t, Xs[-1], A, B, n)), 1e3)

# ---------------------------------------------------------------- masked sums and inclusion-exclusion
P = primes_upto(40)
R = [p for p in P if Nrm(p) in (7, 13, 19, 31)][:4] + [(-5, 0)]
print("mask primes:", R, [Nrm(p) for p in R])
worst_id = 0; worst_mask = 0
for th in ["1", "(l/2)_3"]:
    for e in (1, -1):
        for t in (0.0, 2.5):
            for X in (300.0, 3000.0, 3e4):
                direct = Lsum(e, th, t, X, A, B, n, R)
                ie = 0
                for k in range(len(R) + 1):
                    for sub in itertools.combinations(R, k):
                        d = (1, 0)
                        for p in sub: d = mul(d, p)
                        u, d0 = primary_assoc(d)            # products of primary are primary; d0 == d
                        dn = Nrm(d0); z = cx(d0)
                        thd = theta_val(th, np.array([d0[0]]), np.array([d0[1]]), np.array([dn]))[0]
                        coef = (-1) ** k * (z / abs(z)) ** e * thd * dn ** (1j * t)
                        Xd = X / dn
                        if 2 * Xd >= 1:
                            ie += coef * Lsum(e, th, t, Xd, A, B, n)
                worst_id = max(worst_id, abs(direct - ie))
                worst_mask = max(worst_mask, abs(direct))
report("4b-inclusion-exclusion", worst_id, 1e-9, "masked sum = sum_{d|R} mu(d) alpha(d)^e theta(d) N(d)^{it} L^0(X/N(d))")
report("4b-masked-bounded", worst_mask, 2.0 ** len(R), f"masked |L| <= (number of divisors of R) x O(1); max {worst_mask:.3g}")

# ---------------------------------------------------------------- (2.18f): centered and uncentered products
T = 9e6; worst = 0; worst0 = 0
for (U1, V1) in [(3000.0, 300.0), (30000.0, 1000.0), (1000.0, 9000.0)]:
    U2 = T / U1; V2 = T / V1
    for e in (1, -1):
        diff = Lsum(e, "1", 1.3, U1, A, B, n) * Lsum(e, "1", 1.3, U2, A, B, n) \
             - Lsum(e, "1", 1.3, V1, A, B, n) * Lsum(e, "1", 1.3, V2, A, B, n)
        worst = max(worst, abs(diff) / math.sqrt(T))
    d0 = Lsum(0, "1", 1.3, U1, A, B, n) * Lsum(0, "1", 1.3, U2, A, B, n) \
       - Lsum(0, "1", 1.3, V1, A, B, n) * Lsum(0, "1", 1.3, V2, A, B, n)
    p0 = abs(Lsum(0, "1", 1.3, U1, A, B, n) * Lsum(0, "1", 1.3, U2, A, B, n)) / math.sqrt(T)
    worst0 = max(worst0, p0)
    print(f"   finite-order: product/sqrt(T) = {p0:.3g}, centered difference/sqrt(T) = {abs(d0)/math.sqrt(T):.3g}")
report("4c-centered", worst, 1e-3, "CM: T^{-1/2}|L(U1)L(U2) - L(V1)L(V2)| is tiny (no main term to cancel), T = 9e6")
falsifier("4c-centered", "finite-order uncentered product has size T^{1/2}", worst0, 5.0)
summary()
