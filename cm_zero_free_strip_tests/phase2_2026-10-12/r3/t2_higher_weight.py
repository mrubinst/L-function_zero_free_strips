"""Section 17 (higher weight) of the CM guide: independent checks (r3) of every stated identity.
  H1 chain rule: d_zbar^j F(z',zbar',v') at z=zbar=0 equals (-(cv)^-2)^j (d_z'^j F)(-delta/c, 1/(N(c)v)), j<=5;
     and the symmetric statement for d_z^j with (-(cbar v)^-2)^j d_zbar'^j; falsifier: exponent (cv)^-2j without sign
  H2 direct Mellin transform with j derivatives, j=1..5 (mpmath quadrature); falsifier: Gamma shifted by 1/2
  H3 R_j(1) = (3j-4)(3j-2)(3j+2)(3j+4)/1296, j=1..12 (exact); R_1, R_2 = the two kernels of Section 12
  H4 rightmost pole of R_j at t = 1/6-(j+1)/2 (j=1: -5/6, j=2: -4/3); Stirling |R_j(s+iu)| ~ u^{4s}
  H5 Hecke gamma quotient |Gamma(1-s+k/2)/Gamma(s+k/2)| ~ |t|^{1-2 sigma}, independent of k (k=0..5);
     at sigma=-1/10, k=1: exponent 6/5 (Lemma cm-growth); plain-moment profile Gamma(s+1/2)/Gamma(3/2-s) ~ |t|^{2s-1}
  H6 symmetric powers: local Euler polynomials of Sym^m (m=1..8) at split and inert p equal the product of the
     factors of L(s,psi^{m-2i}) (i<m/2) and of L(s,chi_{-3}^{m/2}) (m even); falsifier: trivial middle character
  H7 eta(3z)^8: a_p = 2 Re(pi^3) (split), a_p = 0, a_{p^2} = -p^3 (p = 2 mod 3); eta(2z)^3 eta(6z)^3:
     a_p = 2 Re(pi^2 chi_pi(4)) (split), a_p = 0, a_{p^2} = p^2 (p = 2 mod 3, p not 2), a_3 = -3, a_2 = a_4 = 0
  H8 expected strips: k/2+7/8 and k/2+1 for k=3 (19/8, 5/2) and k=2 (15/8, 2)
"""
import sys, os, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import mpmath as mp
from r3lib import primary_primes, chi_prime, cx, norm, conj, is_prime_int

ok = True
_FORCE = [os.environ.get('CM_TESTS_FORCE_FAIL') == '1']   # self-test (added 2026-10-12): first check counts as failed
def rep(name, cond, info=''):
    global ok
    if _FORCE[0]:
        _FORCE[0] = False; cond = False; name = '[CM_TESTS_FORCE_FAIL=1] ' + name
    print(('[PASS] ' if cond else '[FAIL] ') + name + ('  :: ' + info if info else ''))
    ok &= bool(cond)

# H1 -------------------------------------------------------------------------------------
z, zb, v, c, cb, dl, dlb = sp.symbols('z zb v c cb delta deltab')
Z1 = -dl / c - zb / (c**2 * (v**2 + z * zb))
Z2 = -dlb / cb - z / (cb**2 * (v**2 + z * zb))
V1 = v / (c * cb * (v**2 + z * zb))
A, B, C = sp.symbols('A B C')
F = sp.Function('F')
res = True; fal = True
for j in range(1, 6):
    expr = F(Z1, Z2, V1)
    d = sp.diff(expr, zb, j).subs({z: 0, zb: 0})
    target = (-(c * v)**-2)**j * sp.diff(F(A, B, C), A, j).subs({A: -dl / c, B: -dlb / cb, C: 1 / (c * cb * v)})
    res &= sp.simplify(d - target) == 0
    wrong = ((c * v)**-2)**j * sp.diff(F(A, B, C), A, j).subs({A: -dl / c, B: -dlb / cb, C: 1 / (c * cb * v)})
    if j % 2 == 1:
        fal &= sp.simplify(d - wrong) != 0
    d2 = sp.diff(expr, z, j).subs({z: 0, zb: 0})
    target2 = (-(cb * v)**-2)**j * sp.diff(F(A, B, C), B, j).subs({A: -dl / c, B: -dlb / cb, C: 1 / (c * cb * v)})
    res &= sp.simplify(d2 - target2) == 0
rep('H1 chain rules d_zbar^j and d_z^j, j=1..5', res and fal, 'falsifier (no sign) fails for odd j')

# H2 -------------------------------------------------------------------------------------
mp.mp.dps = 30
worst = 0; fworst = 1e9
for j in range(1, 6):
    for (l, s) in [(complex(0.7, 0.4), mp.mpc(0.9, 0.3)), (complex(-1.1, 0.5), mp.mpc(1.3, -0.7))]:
        L = mp.mpc(l)
        absl = abs(L)
        f = lambda t: (2j * mp.pi * mp.conj(L))**j * t * mp.besselk(mp.mpf(1) / 3, 4 * mp.pi * absl * t) * t**(2 * s + j - 2)
        lhs = mp.quad(f, [0, 0.05, 0.5, 2, mp.inf])
        al = L / absl
        rhs = (1j**j) / 4 * (2 * mp.pi)**(-2 * s) * mp.gamma(s + mp.mpf(j) / 2 - mp.mpf(1) / 6) * \
            mp.gamma(s + mp.mpf(j) / 2 + mp.mpf(1) / 6) * mp.conj(al)**j * (absl**2)**(-s)
        rhs_bad = (1j**j) / 4 * (2 * mp.pi)**(-2 * s) * mp.gamma(s + mp.mpf(j - 1) / 2 - mp.mpf(1) / 6) * \
            mp.gamma(s + mp.mpf(j - 1) / 2 + mp.mpf(1) / 6) * mp.conj(al)**j * (absl**2)**(-s)
        worst = max(worst, float(abs(lhs - rhs) / abs(rhs)))
        fworst = min(fworst, float(abs(lhs - rhs_bad) / abs(rhs)))
rep('H2 direct Mellin with j derivatives, j=1..5', worst < 1e-20 and fworst > 1e-2,
    'max rel err %.1e; falsifier min rel err %.2f' % (worst, fworst))

# H3 -------------------------------------------------------------------------------------
t = sp.symbols('t')
def Rj(j, tt):
    return sp.prod([sp.gamma(sp.Rational(j + 1, 2) + tt + e) / sp.gamma(sp.Rational(j + 1, 2) - tt + e)
                    for e in (sp.Rational(1, 6), -sp.Rational(1, 6))])
good = True; badf = True
for j in range(1, 13):
    val = sp.gammasimp(Rj(j, 1))
    good &= sp.simplify(val - sp.Rational((3*j - 4) * (3*j - 2) * (3*j + 2) * (3*j + 4), 1296)) == 0
    badf &= sp.simplify(val - sp.Rational((3*j - 4) * (3*j - 1) * (3*j + 2) * (3*j + 4), 1296)) != 0
k1 = sp.gamma(1 + t + sp.Rational(1, 6)) * sp.gamma(1 + t - sp.Rational(1, 6)) / (sp.gamma(1 - t + sp.Rational(1, 6)) * sp.gamma(1 - t - sp.Rational(1, 6)))
k2 = sp.gamma(sp.Rational(4, 3) + t) * sp.gamma(sp.Rational(5, 3) + t) / (sp.gamma(sp.Rational(4, 3) - t) * sp.gamma(sp.Rational(5, 3) - t))
same = sp.simplify(Rj(1, t) - k1) == 0 and sp.simplify(Rj(2, t) - k2) == 0
rep('H3 R_j(1) formula j=1..12; R_1, R_2 = table kernels', good and badf and same)

# H4 -------------------------------------------------------------------------------------
poles_ok = True
for j in range(1, 6):
    t0 = sp.Rational(1, 6) - sp.Rational(j + 1, 2)
    # numerator Gamma((j+1)/2 + t - 1/6) has its rightmost pole at t0
    poles_ok &= (sp.Rational(j + 1, 2) + t0 - sp.Rational(1, 6)) == 0
    poles_ok &= t0 < -sp.Rational(1, 4)
rep('H4a rightmost pole of R_j at 1/6-(j+1)/2 < -1/4 (j=1: -5/6, j=2: -4/3)',
    poles_ok and sp.Rational(1, 6) - 1 == -sp.Rational(5, 6) and sp.Rational(1, 6) - sp.Rational(3, 2) == -sp.Rational(4, 3))
mp.mp.dps = 20
st_ok = True; info = []
for j in range(1, 6):
    for sg in [-0.25, 0.0, 0.5, 2.0]:
        def R(u):
            tt = mp.mpc(sg, u)
            return abs(mp.gamma((j + 1) / 2 + tt + 1 / 6) * mp.gamma((j + 1) / 2 + tt - 1 / 6) /
                       (mp.gamma((j + 1) / 2 - tt + 1 / 6) * mp.gamma((j + 1) / 2 - tt - 1 / 6)))
        e = float(mp.log(R(2e5) / R(1e5)) / mp.log(2))
        st_ok &= abs(e - 4 * sg) < 1e-3
info = 'exponent log2(|R(2e5)|/|R(1e5)|) = 4 sigma within 1e-3'
rep('H4b Stirling |R_j(sigma+iu)| ~ u^{4 sigma}, j=1..5, sigma in {-1/4,0,1/2,2}', st_ok, info)

# H5 -------------------------------------------------------------------------------------
h5 = True
for k in range(0, 6):
    for sg in [-0.1, 0.25, 0.5, 1.1]:
        q = lambda T: abs(mp.gamma(1 - mp.mpc(sg, T) + k / 2) / mp.gamma(mp.mpc(sg, T) + k / 2))
        e = float(mp.log(q(2e5) / q(1e5)) / mp.log(2))
        h5 &= abs(e - (1 - 2 * sg)) < 1e-3
q = lambda T: abs(mp.gamma(mp.mpf(3) / 2 - mp.mpc(-0.1, T)) / mp.gamma(mp.mpf(1) / 2 + mp.mpc(-0.1, T)))
e65 = float(mp.log(q(2e5) / q(1e5)) / mp.log(2))
q2 = lambda T: abs(mp.gamma(mp.mpc(0.3, T) + 0.5) / mp.gamma(1.5 - mp.mpc(0.3, T)))
eprof = float(mp.log(q2(2e5) / q2(1e5)) / mp.log(2))
rep('H5 gamma quotient exponents', h5 and abs(e65 - 1.2) < 1e-3 and abs(eprof - (2 * 0.3 - 1)) < 1e-3,
    'k=0..5 exponent 1-2sigma; sigma=-1/10,k=1: %.4f; profile at sigma=0.3: %.4f' % (e65, eprof))

# H6 -------------------------------------------------------------------------------------
X = sp.symbols('X')
h6 = True; h6f = False
for m in range(1, 9):
    # split p: normalized Frobenius eigenvalues a, 1/a (|a|=1); Sym^m eigenvalues a^{m-2j}
    a = sp.symbols('a')
    lhs = sp.prod([1 - a**(m - 2 * j) * X for j in range(m + 1)])
    rhs = sp.prod([(1 - a**(m - 2 * i) * X) * (1 - a**(-(m - 2 * i)) * X) for i in range(m) if 2 * i < m])
    if m % 2 == 0:
        rhs *= (1 - 1 * X)          # chi_{-3}(p)^{m/2} = 1 at split p
    h6 &= sp.simplify(sp.expand(lhs - rhs)) == 0
    # inert p: eigenvalues +-i of the normalized Frobenius (a_p = 0); psi^{r}(-p) = (-1)^r, N(-p) = p^2
    lhs = sp.expand(sp.prod([1 - (sp.I)**(m - j) * (-sp.I)**j * X for j in range(m + 1)]))
    rhs = sp.prod([(1 - (-1)**(m - 2 * i) * X**2) for i in range(m) if 2 * i < m])
    if m % 2 == 0:
        rhs_good = rhs * (1 - (-1)**(m // 2) * X)
        rhs_bad = rhs * (1 - X)
        h6f |= sp.expand(lhs - rhs_bad) != 0
        rhs = rhs_good
    h6 &= sp.expand(lhs - rhs) == 0
rep('H6 Sym^m local factorization m=1..8 (split and inert)', h6 and h6f,
    'falsifier: middle factor trivial instead of chi_{-3}^{m/2} fails for m = 2 mod 4')

# H7 -------------------------------------------------------------------------------------
Nq = 3000
def series_prod(factors, N):
    """prod over (step, power) of prod_n (1-q^{step n})^power, as a list of ints up to q^N."""
    c = [0] * (N + 1); c[0] = 1
    for step, pw in factors:
        for _ in range(pw):
            for n in range(step, N + 1, step):
                # multiply by (1 - q^n)
                for i in range(N, n - 1, -1):
                    c[i] -= c[i - n]
    return c
f1 = series_prod([(3, 8)], Nq)           # eta(3z)^8 = q * f1(q)
a1 = [0] * (Nq + 2)
for i in range(Nq + 1):
    if i + 1 <= Nq + 1:
        a1[i + 1] = f1[i]
f2 = series_prod([(2, 3), (6, 3)], Nq)   # eta(2z)^3 eta(6z)^3 = q * f2(q)
a2 = [0] * (Nq + 2)
for i in range(Nq + 1):
    a2[i + 1] = f2[i]
P = primary_primes(Nq)
mis1 = mis2 = 0; nsp = 0
for pi in P:
    if pi[1] == 0:
        continue
    p = norm(pi)
    if p > Nq:
        continue
    z3 = cx(pi)**3
    mis1 += abs(2 * z3.real - a1[p]) > 1e-6
    c4 = cx(chi_prime(pi, (4, 0)))
    mis2 += abs(2 * (cx(pi)**2 * c4).real - a2[p]) > 1e-6
    nsp += 1
inert_ok = all(a1[p] == 0 for p in range(2, Nq + 1) if is_prime_int(p) and p % 3 == 2)
inert_sq = all(a1[p * p] == -p**3 for p in range(2, 55) if is_prime_int(p) and p % 3 == 2)
inert2 = all(a2[p] == 0 for p in range(5, Nq + 1) if is_prime_int(p) and p % 3 == 2) and \
    all(a2[p * p] == p**2 for p in range(5, 55) if is_prime_int(p) and p % 3 == 2)
rep('H7 eta-product coefficients (%d split primary generators, q^%d)' % (nsp, Nq),
    mis1 == 0 and mis2 == 0 and inert_ok and inert_sq and inert2 and a2[3] == -3 and a2[2] == 0 and a2[4] == 0,
    'mismatches %d, %d; a_3=%d a_2=%d a_4=%d (eta(2z)^3eta(6z)^3); a_4(eta(3z)^8)=%d' % (mis1, mis2, a2[3], a2[2], a2[4], a1[4]))

# H8 -------------------------------------------------------------------------------------
h8 = all([sp.Rational(3, 2) + sp.Rational(7, 8) == sp.Rational(19, 8), sp.Rational(3, 2) + 1 == sp.Rational(5, 2),
          1 + sp.Rational(7, 8) == sp.Rational(15, 8), 1 + 1 == 2, sp.Rational(1, 2) + sp.Rational(7, 8) == sp.Rational(11, 8)])
rep('H8 expected strips 19/8<Re s<5/2, 15/8<Re s<2; weight 2: 11/8', h8)
print('ALL PASS' if ok else 'SOME FAIL')

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
