"""Part D, test 05: the remaining identities and numbers of Section 'CM newforms of higher weight' and the
counting statements of the appendix.

  L1521-1524, 1542, 1548  strip arithmetic (exact fractions)
  L1559, 1630, 1638       the weight gamma_2(c) conj(alpha(cn^3)) psi(cn^3) for psi = conj(alpha)^k eta
  L1608                   Poisson transform of zbar^k e^{-pi v |z|^2} is a multiple of xibar^k e^{-pi|xi|^2/v} v^{-1-k}
  L1612-1613              co-oriented classes mu conj(alpha)^k <-> conj(alpha)^{k+1} gamma_2 (from eq. (gauss)),
                          and the cube factor conj(alpha(b))^{3(k+1)}
  L1616                   int alpha(z)^{+-k} f(|z|^2) dz = 0  (k >= 1)
  L1622                   kappa = 3/4
  L1635                   L(s, psi N^{it}) = L(s - it, psi)
  L1692                   at most 5184 * 2^{|P|} terms; c_0 | 8 lambda^2 (consistency with OAI78's count N(L) 2^{|P|}, L = 72)
  L1750-1751              (b^3/k)_2 = (b/k)_2
"""
import os, sys, math, cmath, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import *
from fastgauss import primary_primes
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

random.seed(11)
PR = [p for p in primary_primes(5000) if p != (-2, 0)]

def factor(n):
    out = []; m = n
    for p in PR:
        if N(m) == 1: break
        while divides(p, m):
            out.append(p); m = div(m, p)
    assert N(m) == 1, (n, m)
    return out

def chi6(n, x):
    r = 1 + 0j
    for p in factor(n):
        j = sym6_prime(x, p)
        if j is None: return 0
        r *= cmath.exp(1j*math.pi*j/3)
    return r

def mu(n):
    f = factor(n); return 0 if len(set(f)) < len(f) else (-1)**len(f)

def al(x):
    z = complex(cx(x)); return z/abs(z)

def gamma_j(j, n):
    nn = complex(cx(n)); s = 0
    for x in residues(n):
        c = chi6(n, x)
        if c != 0:
            s += c**j*e_add(complex(cx(x))/nn)
    return s/math.sqrt(N(n))

# ---------------------------------------------------------------- strips (exact)
ok = True
for k in range(1, 30):
    w = k + 1
    ok &= (Fr(w - 1, 2), Fr(w + 1, 2)) == (Fr(k, 2), Fr(k, 2) + 1)              # critical strip of weight k+1
    ok &= Fr(k, 2) + Fr(7, 8) == Fr(7, 8) + Fr(k, 2)
ok &= (Fr(3, 2) + Fr(7, 8), Fr(3, 2) + 1, Fr(3, 2) + Fr(1, 2)) == (Fr(19, 8), Fr(5, 2), Fr(2))
ok &= (Fr(2, 2) + Fr(7, 8), Fr(2, 2) + 1, Fr(2, 2) + Fr(1, 2)) == (Fr(15, 8), Fr(2), Fr(3, 2))
fal = Fr(4, 2) + Fr(7, 8) == Fr(19, 8)        # using the weight instead of k
report("D05.1", 1542, "strips", "critical strip k/2<Re s<k/2+1; boundaries k/2+7/8 (19/8 for k=3, 15/8 for k=2), centres 2 and 3/2; 11/12 analogue k/2+11/12",
       ok and not fal, "falsifier: weight/2+7/8 = 23/8 != 19/8")

# ---------------------------------------------------------------- weight identities L1559 / L1630 / L1638
sq = [n for n in (prim_assoc((a, b)) for a in range(-25, 26) for b in range(-25, 26)
                  if (a, b) != (0, 0) and N((a, b)) % 3 and N((a, b)) % 2 and N((a, b)) < 400) if mu(n) != 0]
sq = sorted(set(sq), key=N)
u_eta = (7, 3)                     # eta = chi_.(u_eta), a finite-order character
err = 0; errf = 0; err0 = 0; errc = 0
for _ in range(60):
    c, n = random.choice(sq), random.choice(sq)
    m = mul(c, pw(n, 3)); k = random.randint(0, 5)
    eta = chi6(m, u_eta)
    psi = np.conj(al(m))**k*eta
    lhs = np.conj(al(m))*psi
    err = max(err, abs(lhs - np.conj(al(m))**(k + 1)*eta))
    errf = max(errf, abs(lhs - np.conj(al(m))**k*eta))
    if k == 0: err0 = max(err0, abs(lhs - np.conj(al(m))*eta))
    # conjugate target alpha^k conj(eta): weight conj(alpha) alpha^k conj(eta) = alpha^{k-1} conj(eta)
    errc = max(errc, abs(np.conj(al(m))*al(m)**k*np.conj(eta) - al(m)**(k - 1)*np.conj(eta)))
report("D05.2", 1559, "weight", "conj(alpha(cn^3)) psi(cn^3) = conj(alpha(cn^3))^{k+1} eta(cn^3); k=0 gives conj(alpha) eta (L1630); conjugate target leaves alpha^{k-1} (L1638)",
       max(err, err0, errc) < 1e-12 and errf > 0.1, f"max err={max(err, err0, errc):.1e}; falsifier exponent k: {errf:.2f}")

# ---------------------------------------------------------------- Poisson transform L1608
mp.mp.dps = 20
def I1(a, xi, v):
    return mp.quad(lambda x: x**a*mp.exp(-mp.pi*v*x*x - 2j*mp.pi*x*xi), [-mp.inf, 0, mp.inf])
def ft(k, xi1, xi2, v):
    return sum(mp.binomial(k, a)*(-1j)**(k - a)*I1(a, xi1, v)*I1(k - a, xi2, v) for a in range(k + 1))
ratios = {}; ratios_f = {}; zero_ok = True
for k in range(0, 6):
    for (xi1, xi2, v) in [(0.3, -0.2, 0.7), (0.5, 0.4, 1.3), (-0.25, 0.6, 2.1), (0.1, 0.15, 0.5)]:
        F = ft(k, xi1, xi2, v)
        base = (xi1 - 1j*xi2)**k*mp.exp(-mp.pi*(xi1**2 + xi2**2)/v)
        ratios.setdefault(k, []).append(F/(base*v**(-1 - k)))
        ratios_f.setdefault(k, []).append(F/(base*v**(-k)))
    if k >= 1:
        zero_ok &= abs(ft(k, 0, 0, 0.9)) < 1e-15
spread = max(max(abs(r - (1j)**(-k)) for r in rs) for k, rs in ratios.items())
spread_f = min(max(abs(r - rs[0]) for r in rs) for k, rs in ratios_f.items())
report("D05.3", 1608, "Poisson transform", "FT of zbar^k e^{-pi v|z|^2} = i^{-k} v^{-1-k} xibar^k e^{-pi|xi|^2/v} (k=0..5), = 0 at xi=0 for k>=1",
       spread < 1e-12 and zero_ok and spread_f > 0.1, f"max |ratio - i^-k|={mp.nstr(spread, 3)}; falsifier v^-k: ratio varies by {mp.nstr(spread_f, 3)}")

# ---------------------------------------------------------------- co-oriented classes L1612-1613
def G(n): return np.conj(chi6(n, (4, 0)))*gamma_j(3, n)
err = 0; errf = 0; cnt = 0
for n in sq[:40]:
    g2 = gamma_j(2, n); gm1 = gamma_j(-1, n)
    for k in range(0, 5):
        lhs = mu(n)*np.conj(al(n))**k*gm1
        rhs = chi6(n, (-1, 0))/G(n)*np.conj(al(n))**(k + 1)*g2
        err = max(err, abs(lhs - rhs)); errf = max(errf, abs(lhs - chi6(n, (-1, 0))/G(n)*np.conj(al(n))**k*g2)); cnt += 1
# cube factor: conj(alpha(l))^{k+1} for l = lambda^-3 n b^3 contains conj(alpha(b))^{3(k+1)}
errb = 0; errbf = 0
lam3 = complex(cx(LAM))**-3
for _ in range(50):
    n, b = random.choice(sq), random.choice(sq); k = random.randint(1, 5)
    l = lam3*complex(cx(n))*complex(cx(b))**3; a = l/abs(l)
    rest = np.conj(lam3/abs(lam3))**(k + 1)*np.conj(al(n))**(k + 1)
    errb = max(errb, abs(np.conj(a)**(k + 1) - rest*np.conj(al(b))**(3*(k + 1))))
    errbf = max(errbf, abs(np.conj(a)**(k + 1) - rest*np.conj(al(b))**(3*k)))
report("D05.4", 1613, "co-oriented classes", "mu conj(alpha)^k gamma_{-1} = chi(-1) G^{-1} conj(alpha)^{k+1} gamma_2 (k=0..4, 40 squarefree n); cube factor conj(alpha(b))^{3(k+1)}",
       err < 1e-10 and errb < 1e-12 and errf > 0.1 and errbf > 0.1,
       f"{cnt} checks, max err={err:.1e}, cube err={errb:.1e}; falsifiers conj(alpha)^k: {errf:.2f}, conj(alpha(b))^(3k): {errbf:.2f}")

# ---------------------------------------------------------------- int alpha^{+-k} f(|z|^2) = 0
def grid_int(h):
    xs = np.arange(-7 + h/2, 7, h)              # midpoint grid (avoids z = 0)
    X, Y = np.meshgrid(xs, xs)
    r2 = X*X + Y*Y; f = np.exp(-r2)*(1 + r2)
    A = (X + 1j*Y)/np.sqrt(r2)
    return {k: abs(np.sum(A**k*f)*h*h) for k in range(-6, 7)}
# the midpoint grid has 90-degree symmetry, so k = 4, 8 cancel only in the limit h -> 0 (error O(h^2)):
v1, v2 = grid_int(0.01), grid_int(0.005)
e1 = max(v for k, v in v1.items() if k); e2 = max(v for k, v in v2.items() if k)
ok = e2 < 1e-4 and e1/e2 > 3 and abs(v2[0] - 2*math.pi) < 1e-6
report("D05.5", 1616, "angular main term", "int_C alpha(z)^{+-k} f(|z|^2) dxdy = 0 for k=1..6 (f = e^{-r^2}(1+r^2)), midpoint grids h=0.01, 0.005",
       ok, f"max |int| (k!=0): {e1:.1e} (h=0.01) -> {e2:.1e} (h=0.005), O(h^2) -> 0; falsifier k=0 gives {v2[0]:.6f} = 2 pi")

# ---------------------------------------------------------------- kappa
report("D05.6", 1622, "kappa", "kappa = 2*(7/8) - 1 = 3/4 (OpenAI's finite-order 7/8 theorem)",
       2*Fr(7, 8) - 1 == Fr(3, 4) and 2*Fr(11, 12) - 1 != Fr(3, 4), "falsifier 2*(11/12)-1 = 5/6")

# ---------------------------------------------------------------- L(s, psi N^{it}) = L(s - it, psi)
mp.mp.dps = 20
vals = [(N(p), al(p)**2*chi6(p, (4, 0))) for p in PR if N(p) < 3000]
def L(s, t=0):
    r = mp.mpf(1)
    for q, v in vals:
        r *= 1/(1 - v*mp.mpf(q)**(1j*t)*mp.mpf(q)**(-s))
    return r
s = mp.mpc(2.5, 0.3); t = 1.7
e = abs(L(s, t) - L(s - 1j*t)); ef = abs(L(s, t) - L(s + 1j*t))
report("D05.7", 1635, "vertical shift", "L(s, psi N^{it}) = L(s-it, psi) (truncated Euler products, psi = alpha^2 chi(4), s=2.5+0.3i, t=1.7)",
       e < 1e-15 and ef > 1e-4, f"err={mp.nstr(e, 3)}; falsifier L(s+it): {mp.nstr(ef, 3)}")

# ---------------------------------------------------------------- term count, c_0 | 8 lambda^2 (L = 72)
L72 = (72, 0); NL = N(L72)
lam2 = pw(LAM, 2)                               # lambda^2 = -3
eight_lam2 = mul((8, 0), lam2)                  # 8 lambda^2 = -24
# lambda^2 h0 / 72 = -h0/24, so for every h0 mod 72 the reduced denominator c_0 divides 24 = -8 lambda^2
den_ok = lam2 == (-3, 0) and eight_lam2 == (-24, 0)
for h0 in residues(L72):
    num = mul(lam2, h0)
    den_ok &= divides((72, 0), mul(num, (24, 0)))  # 24 * (lambda^2 h0/72) is integral, i.e. c_0 | 24
report("D05.8", 1692, "term count", "5184 = N(72) = number of h_0 mod L for L = 72 (OAI78 counts terms by h_0 mod L and active subsets: N(L) 2^{|P|}); "
       "lambda^2 h_0/72 = -h_0/24 so c_0 | 24 = 8|lambda|^2 (= -8 lambda^2)",
       NL == 5184 and den_ok and N((36, 0)) != 5184, f"N(72)={NL}; falsifier L=36 gives {N((36, 0))}. |C| <= 1/81 is OAI78's |zeta| <= 1/81 (cited)")

# ---------------------------------------------------------------- (b^3/k)_2 = (b/k)_2
bad = 0; badf = 0; cnt = 0
for _ in range(400):
    k = random.choice(sq); b = (random.randint(-40, 40), random.randint(-40, 40))
    if b == (0, 0): continue
    q1 = chi6(k, b)**3; q3 = chi6(k, pw(b, 3))**3; q2 = chi6(k, pw(b, 2))**3
    if abs(q1 - q3) > 1e-9: bad += 1
    if abs(q1 - q2) > 1e-9: badf += 1
    cnt += 1
report("D05.9", 1751, "(b^3/k)_2", "(b^3/k)_2 = (b/k)_2 (quadratic symbol = sextic symbol cubed), 400 random pairs",
       bad == 0 and badf > 50, f"mismatches={bad}/{cnt}; falsifier (b^2/k)_2: {badf} mismatches")

# ---------------------------------------------------------------- cusp scalars of the table (line 1226) from the general-j rule
# The single-cusp reflection with j derivatives carries (-1)^j conj(alpha(c))^{2j} (L1585 and test D08).  Writing the
# completed sum in the variables n, b uses conj(alpha(l))^j = conj(alpha(lambda^-3))^j conj(alpha(n b^3))^j, so the cusp
# scalar of the completed reflection has j-dependence (-1)^j conj(alpha(lambda^-3))^{-j} conj(alpha(c))^{2j}.
lam = complex(cx(LAM)); a3 = np.conj((lam**-3)/abs(lam**-3))
ph = {j: (-1)**j*a3**(-j) for j in (1, 2)}
ok = abs(a3 - (-1j)) < 1e-12 and abs(ph[1] - (-1j)) < 1e-12 and abs(ph[2] - (-1)) < 1e-12 and abs(a3**2 - (-1)) < 1e-12
fal = abs((-1)**1 - (-1j)) > 0.5
report("D05.10", 1585, "cusp scalar phases", "(-1)^j conj(alpha(lambda^-3))^{-j} = -i (j=1), -1 (j=2): matches the table's -(i/81)conj(alpha(c))^2 and -(1/81)conj(alpha(c))^4 (line 1226); conj(alpha(lambda^-3))^2 = -1 (line 748)",
       ok and fal, f"conj(alpha(lambda^-3))={a3:.3f}; falsifier (-1)^j alone gives -1 for j=1, not -i")
