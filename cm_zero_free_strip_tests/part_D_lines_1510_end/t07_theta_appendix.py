"""Part D, test 07: the explicit theta coefficients of the appendix (lines 1700-1754).

Claims tested
  L1711-1722  the table of tau(l) (Patterson's coefficients at infinity, as in [DR, (5.7)])
  L1723-1735  the tables of tau_1, tau_2 (other cusps, [DR, (5.13)-(5.14)])
  L1738-1744  t_0 = tau, t_+(l) = omega tau_2(omega l) ce(l), t_-(l) = omega^2 tau_1(omega^2 l) ce(l), d_s(l) = conj(t_s(-l))
              "are precisely the coefficients of the conjugated theta functions at the three cusps"
  L1744-1748  support in {u lambda^m n b^3: m >= -4, n squarefree, n, b primary} and |d_s(u lambda^m n b^3)| <= 27 3^{m/6}|b|
  L1704-1707  the cubic symbol at -2 and |g(r,c)| = |c| (used for the bound)
  (cross-check) eq. (theta-coeff), line 730: d_0(lambda^-3 n b^3) = 3^{5/2}|b| conj((lambda/n)_3) gamma_2(n)

Independent checks ("direct computation"):
  (a) theta(w) = 3^{5/2}/2 v^{2/3} + sum tau(mu) v K_{1/3}(4 pi|mu| v) ce(mu z) is invariant under z -> z+1 and
      under the inversion E = (0 -1; 1 0) in SL_2(Z) (Kubota multiplier 1 on SL_2(Z) [DR, Sec. 5.1]);
  (b) theta(gamma_10 w) and theta(gamma_19 w), gamma_10 = (1 0; omega 1), gamma_19 = (1 0; omega^2 1) (the cusps
      F_10, F_19 of [DR]; Theta_2, Theta_1 of OAI78), computed from the expansion at infinity, equal
      sum t_+(mu) v K ce(mu z) and sum t_-(mu) v K ce(mu z) respectively, with no constant term;
  (c) the conjugated expansions have coefficients d_s.
Falsifiers: zeta_9 <-> zeta_9^{-1} in the tau_1 table; t_+ and t_- swapped; the bound with 9 in place of 27;
the symbol at -2 replaced by its conjugate.
"""
import os, sys, math, cmath, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pd_common import report, info
import appendix_theta as AT
from zw import mul, neg, N, cx, divides, div, pw, primary, LAM, residues, W, sym6_prime, UNITS
import numpy as np
from scipy.special import kv

random.seed(5)
SIG0 = 3**2.5/2

def lattice(maxnorm):
    R = int(2*math.sqrt(maxnorm/3)) + 3
    return [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if 0 < a*a - a*b + b*b <= maxnorm]

# ---------------------------------------------------------------- coefficient lists
# at infinity: mu in lambda^-3 O, mu = n / lambda^3, m = lambda^4 mu = lambda n
INF = []
for n in lattice(260):
    m = mul(LAM, n)
    t = AT.tau(m)
    if t != 0:
        INF.append((complex(cx(n))/AT.LAMc**3, t))
# other cusps: mu in lambda^-4 O, mu = m/9
def cusp_list(sig, maxnorm=760, variant=None):
    out = []
    for m in lattice(maxnorm):
        t = AT.t_coeff(sig, m, variant)
        if t != 0:
            out.append((complex(cx(m))/9, t))
    return out
CP, CM = cusp_list('+'), cusp_list('-')
CMf = cusp_list('-', variant='z9swap')
info(f"coefficients: {len(INF)} at infinity, {len(CP)} (sigma=+), {len(CM)} (sigma=-)")

def series(coefs, z, v, const=0.0):
    mu = np.array([c[0] for c in coefs]); co = np.array([c[1] for c in coefs])
    return const*v**(2/3) + np.sum(co*v*kv(1/3, 4*np.pi*np.abs(mu)*v)*np.exp(4j*np.pi*(mu*z).real))

def theta(z, v):
    return series(INF, z, v, SIG0)

def act(g, z, v):
    a, b, c, d = g
    den = abs(c*z + d)**2 + abs(c)**2*v*v
    return ((a*z + b)*np.conj(c*z + d) + a*np.conj(c)*v*v)/den, v/den

w = complex(cx((0, 1)))
# ---- (a) automorphy at infinity
pts = [(0.13 + 0.21j, 0.97), (-0.21 + 0.05j, 1.02), (0.31 - 0.17j, 0.93)]
err = 0
for z, v in pts:
    t0 = theta(z, v)
    err = max(err, abs(theta(z + 1, v) - t0)/abs(t0))
    z2, v2 = act((0, -1, 1, 0), z, v)
    err = max(err, abs(theta(z2, v2) - t0)/abs(t0))
# falsifier: swap zeta_9 and zeta_9^{-1} in the tau table rows 2,3 (rebuild the list)
INFf = []
for n in lattice(260):
    m = mul(LAM, n)
    t = AT.tau(m)
    if t != 0:
        dec = AT.decompose(m); u = dec[1]
        if u in [(0, 1), (0, -1)]: t *= AT.Z9**2
        if u in [(-1, -1), (1, 1)]: t *= AT.Z9**-2
        INFf.append((complex(cx(n))/AT.LAMc**3, t))
z, v = pts[0]
z2, v2 = act((0, -1, 1, 0), z, v)
errf = abs(series(INFf, z2, v2, SIG0) - series(INFf, z, v, SIG0))/abs(series(INFf, z, v, SIG0))
report("D07.1", 1711, "tau table", "theta from the tau table is invariant under z->z+1 and (z,v)->(-zbar,v)/(|z|^2+v^2) (E in SL_2(Z))",
       err < 1e-9 and errf > 1e-3, f"max rel err={err:.1e} over {len(pts)} points; falsifier zeta_9<->zeta_9^-1 in rows 2,3: {errf:.2e}")

# ---- (b) the other two cusps
g10 = (1, 0, w, 1); g19 = (1, 0, w*w, 1)
errp = errm = 0; fsw = fz9 = 1e9
for (z0, v) in [(-w*w + 0.07 - 0.04j, 1.0), (-w*w - 0.11 + 0.06j, 0.95), (-w*w + 0.02 + 0.13j, 1.05)]:
    z1, v1 = act(g10, z0, v)
    lhs = theta(z1, v1); rhs = series(CP, z0, v)
    errp = max(errp, abs(lhs - rhs)/abs(lhs))
    fsw = min(fsw, abs(lhs - series(CM, z0, v))/abs(lhs))
for (z0, v) in [(-w + 0.07 - 0.04j, 1.0), (-w - 0.11 + 0.06j, 0.95), (-w + 0.02 + 0.13j, 1.05)]:
    z1, v1 = act(g19, z0, v)
    lhs = theta(z1, v1); rhs = series(CM, z0, v)
    errm = max(errm, abs(lhs - rhs)/abs(lhs))
    fz9 = min(fz9, abs(lhs - series(CMf, z0, v))/abs(lhs))
report("D07.2", 1739, "t_+, t_- (cusps)", "theta((1 0; omega 1)w) = sum t_+(mu) vK ce(mu z), theta((1 0; omega^2 1)w) = sum t_-(mu) vK ce(mu z), no constant term",
       max(errp, errm) < 1e-9 and fsw > 1e-3 and fz9 > 1e-3,
       f"rel err +:{errp:.1e}, -:{errm:.1e}; falsifiers: t_+ <-> t_- {fsw:.2e}, zeta_9 <-> zeta_9^-1 in tau_1 {fz9:.2e}")

# ---- (c) conjugated expansions: conj(theta(gamma w)) = C(v) + sum d_s(mu) vK ce(mu z)
def dlist(sig, maxnorm):
    out = []
    for m in lattice(maxnorm):
        d = AT.d_coeff(sig, m)
        if d != 0:
            out.append((complex(cx(m))/9, d))
    return out
D0, DP, DM = dlist('0', 760), dlist('+', 760), dlist('-', 760)
err = 0; errf = 0
for (gm, Dl, z0, v) in [(None, D0, 0.13 + 0.21j, 0.97), (g10, DP, -w*w + 0.07 - 0.04j, 1.0), (g19, DM, -w + 0.02 + 0.13j, 1.05)]:
    if gm is None:
        lhs = np.conj(theta(z0, v)); rhs = series(Dl, z0, v, SIG0)
        rhs_f = series([(mu, np.conj(dd)) for mu, dd in Dl], z0, v, SIG0)
    else:
        z1, v1 = act(gm, z0, v); lhs = np.conj(theta(z1, v1)); rhs = series(Dl, z0, v)
        rhs_f = series([(mu, np.conj(dd)) for mu, dd in Dl], z0, v)
    err = max(err, abs(lhs - rhs)/abs(lhs)); errf = max(errf, abs(lhs - rhs_f)/abs(lhs))
report("D07.3", 1741, "d_sigma", "d_s(l) = conj(t_s(-l)) are the coefficients of the conjugated theta functions at the three cusps",
       err < 1e-9 and errf > 1e-3, f"max rel err={err:.1e}; falsifier d_s(l) = t_s(-l) (no conjugation): {errf:.2e}")

# ---- (d) support and the bound |d| <= 27 3^{m/6} |b|
viol = 0; viol_f = 0; supp_bad = 0; cnt = 0; ratio_max = 0
for sig in ('0', '+', '-'):
    for m in lattice(3000):
        d = AT.d_coeff(sig, m)
        if d == 0: continue
        cnt += 1
        dec = AT.decompose(m)
        if dec is None:
            supp_bad += 1; continue
        k, u, c, b = dec
        mm = k - 4
        if mm < -4: supp_bad += 1
        bound = 27*3**(mm/6)*abs(complex(cx(b)))
        ratio_max = max(ratio_max, abs(d)/bound)
        if abs(d) > bound*(1 + 1e-9): viol += 1
        if abs(d) > 9*3**(mm/6)*abs(complex(cx(b)))*(1 + 1e-9): viol_f += 1
report("D07.4", 1748, "support/bound", "every nonzero d_s(l), N(9l) <= 3000, has l = u lambda^m n b^3 with m >= -4, n sf primary, and |d_s| <= 27 3^{m/6}|b|",
       viol == 0 and supp_bad == 0 and viol_f > 0 and abs(ratio_max - 1) < 1e-9,
       f"{cnt} nonzero coefficients, violations={viol}, bad support={supp_bad}, max |d|/bound={ratio_max:.6f} (attained: row 4); falsifier 9*3^(m/6)|b|: {viol_f} violations")

# ---- (e) |g(r,c)| = |c| for squarefree primary c (incl. even), r in {1, lambda^2, omega lambda^2, omega^2 lambda^2}
from fastgauss import primary_primes, g_prime
PR = primary_primes(200)
sqf = [(1, 0)]
for p in PR:
    for c in list(sqf):
        cc = mul(c, p)
        if N(cc) <= 300: sqf.append(cc)
bad = 0; fast_err = 0
for c in sqf:
    for r in [(1, 0), AT.L2, mul((0, 1), AT.L2), mul((-1, -1), AT.L2)]:
        if abs(abs(AT.g(r, c)) - math.sqrt(N(c))) > 1e-8: bad += 1
    if c in PR:
        fast_err = max(fast_err, abs(AT.g((1, 0), c) - g_prime(c)))
# falsifier for the symbol at -2: with (x/(-2))_3 replaced by its conjugate, theta is no longer automorphic
_old = AT.cubic_sym_prime
def _conj_at_2(x, pi):
    v = _old(x, pi)
    return np.conj(v) if pi == (-2, 0) else v
AT.cubic_sym_prime = _conj_at_2; AT._gcache.clear(); AT._chartab.clear()
INF2 = []
for n in lattice(260):
    m = mul(LAM, n); t = AT.tau(m)
    if t != 0: INF2.append((complex(cx(n))/AT.LAMc**3, t))
AT.cubic_sym_prime = _old; AT._gcache.clear(); AT._chartab.clear()
z, v = pts[0]; z2, v2 = act((0, -1, 1, 0), z, v)
f2 = abs(series(INF2, z2, v2, SIG0) - series(INF2, z, v, SIG0))/abs(series(INF2, z, v, SIG0))
report("D07.5", 1706, "symbol at -2 / |g|", "|g(r,c)| = |c| for all squarefree primary c (incl. even), N(c) <= 300, four twists r; (x/(-2))_3 = x mod 2",
       bad == 0 and fast_err < 1e-9 and f2 > 1e-7,
       f"{len(sqf)} moduli, |g| != |c| in {bad} cases; brute force = fast prime Gauss sums (err {fast_err:.1e}); falsifier conj symbol at -2 breaks automorphy: {f2:.2e} (numerical noise 4e-16)")

# ---- (f) cross-check with eq. (theta-coeff), line 730
def chi6(n, x):
    r = 1 + 0j
    for p in AT.factor_any(n):
        j = sym6_prime(x, p)
        if j is None: return 0
        r *= cmath.exp(1j*math.pi*j/3)
    return r
def gamma2(n):
    nn = complex(cx(n)); s = 0
    for x in residues(n):
        c = chi6(n, x)
        if c != 0: s += c**2*cmath.exp(4j*math.pi*(complex(cx(x))/nn).imag/math.sqrt(3))
    return s/math.sqrt(N(n))
odd_sqf = [c for c in sqf if not divides((2, 0), c) and N(c) > 1][:25]
err = 0; errf = 0
for n in [(1, 0)] + odd_sqf:
    for b in [(1, 0), (-2, -3), (1, 3), (4, 3)]:
        m = mul(pw(LAM, 1), mul(n, pw(b, 3)))           # 9 * lambda^-3 n b^3 = lambda n b^3
        d0 = AT.d_coeff('0', m)
        cub = chi6(n, LAM)**2 if N(n) > 1 else 1
        g2 = gamma2(n) if N(n) > 1 else 1
        rhs = 3**2.5*math.sqrt(N(b))*np.conj(cub)*g2
        err = max(err, abs(d0 - rhs)/abs(rhs)); errf = max(errf, abs(d0 - 3**2.5*math.sqrt(N(b))*cub*g2)/abs(rhs))
report("D07.6", 730, "cross-check eq.(theta-coeff)", "row 4 of the tau table gives d_0(lambda^-3 n b^3) = 3^{5/2}|b| conj((lambda/n)_3) gamma_2(n)",
       err < 1e-9 and errf > 1e-3, f"max rel err={err:.1e} (26 n x 4 b); falsifier (lambda/n)_3 unconjugated: {errf:.2f}")

# ---- (g) transcription against Dunn-Radziwill (references/DR_2109.07463/src/output.tex, eqs. cubiccoeff, cubiccoeff1,
#      cubiccoeff2 and Appendix rows 1, 10, 19): implemented here independently in DR's notation and compared.
def dr_tau1(m):
    dec = AT.decompose(m)
    if dec is None or dec[0] != 0: return 0
    k, u, c, b = dec; f = 9*AT._absbc(b, c)
    if u == (1, 0): return 9*W*np.conj(AT.g(AT.L2, c))*AT._absbc(b, c)
    if u == (0, 1): return 9*cmath.exp(-2j*math.pi/9)*W*W*np.conj(AT.g(mul((0, 1), AT.L2), c))*AT._absbc(b, c)
    if u == (-1, -1): return 9*cmath.exp(2j*math.pi/9)*np.conj(AT.g(mul((-1, -1), AT.L2), c))*AT._absbc(b, c)
    return 0
def dr_tau2(m):
    dec = AT.decompose(m)
    if dec is None or dec[0] != 0: return 0
    k, u, c, b = dec
    if u == (-1, 0): return 9*W*W*np.conj(AT.g(AT.L2, c))*AT._absbc(b, c)
    if u == (0, -1): return 9*cmath.exp(-2j*math.pi/9)*np.conj(AT.g(mul((0, 1), AT.L2), c))*AT._absbc(b, c)
    if u == (1, 1): return 9*W*cmath.exp(2j*math.pi/9)*np.conj(AT.g(mul((-1, -1), AT.L2), c))*AT._absbc(b, c)
    return 0
def dr_d(j, m):           # DR appendix: d_1 = tau(mu), d_10 = omega tau_2(omega mu) ce(mu), d_19 = omega^2 tau_1(omega^2 mu) ce(mu)
    l = complex(cx(m))/9
    if j == 1: return AT.tau(m)
    if j == 10: return W*dr_tau2(mul((0, 1), m))*AT.ce(l)
    if j == 19: return W*W*dr_tau1(mul((-1, -1), m))*AT.ce(l)
mx = 0
for m in lattice(1500):
    mx = max(mx, abs(dr_d(10, m) - AT.t_coeff('+', m)), abs(dr_d(19, m) - AT.t_coeff('-', m)),
             abs(dr_tau1(m) - AT.tau1(m)), abs(dr_tau2(m) - AT.tau2(m)))
report("D07.7", 1724, "transcription vs DR", "tau_1, tau_2 and t_+ = DR d_10, t_- = DR d_19 agree with DR's (5.13)-(5.14) and Appendix rows 10, 19 (all 9l with N <= 1500)",
       mx < 1e-12, f"max |difference| = {mx:.1e}; tau table identical to DR (5.7) by inspection (and D07.1)")

# ---- (h) absolute convergence of the reflected series (line 1696), partial: with the bound V#(y) << min(y^{1/4}, y^{-2})
#      (D04.4) the majorant sum_l |d_s(l)| N(l)^{-1/2} min(y^{1/4}, y^{-2}), y = N(l) X/N(c)^2 (here X/N(c)^2 = 1),
#      has a tail that becomes negligible, while without the decaying weight it grows like sqrt(Y) (falsifier).
terms = []
for sig in ('0', '+', '-'):
    for m in lattice(3000):
        d = AT.d_coeff(sig, m)
        if d != 0:
            nl = N(m)/81
            terms.append((nl, abs(d)/math.sqrt(nl)))
def S(Y, weighted=True):
    return sum(t*(min(nl**0.25, nl**-2) if weighted else 1) for nl, t in terms if nl <= Y)
Ymax = 3000/81
St = S(Ymax)
t1 = (S(Ymax/2) - S(Ymax/4))/St; t2 = (S(Ymax) - S(Ymax/2))/St
u1 = S(Ymax/2, False) - S(Ymax/4, False); u2 = S(Ymax, False) - S(Ymax/2, False)
report("D07.8", 1696, "absolute convergence", "majorant sum |d_s(l)| N(l)^-1/2 min(y^1/4, y^-2): dyadic tail pieces decay geometrically (Y = 3000/81)",
       t2 < t1/2 and t2 < 1e-2 and u2 > u1,
       f"dyadic pieces {t1:.1e} -> {t2:.1e} (ratio {t2/t1:.2f}, expected ~2^-1.5); falsifier (no decaying weight): pieces grow {u1:.0f} -> {u2:.0f}. Partial: finite range")
