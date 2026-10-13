"""t04: Lemma lem:cm-growth and the analytic inputs around it (guide L1082-1092, L1259-1289, L1410-1414).

Claims tested
  C37 L1261-1265 conductor in the functional equation: for psi*_d (t01), 3 N(f_psi) = conductor of E_d, where f_psi is the
                 conductor of psi*_d computed by brute force (smallest modulus g with psi((a)) = abar(a) for a = 1 mod g).
  C38 L1263-1265 Lambda(s,psi) = (3Q)^{s/2}(2pi)^{-s} Gamma(s+1/2) L(s,psi) = eps Lambda(1-s, conj psi), |eps| = 1
                 (L(s,psi) = L(E_d, s+1/2) from PARI; eps = root number).
  C39 L1279-1282 Stirling: |Gamma(3/2-s)/Gamma(1/2+s)| and |Gamma(1-s)/Gamma(s)| are both ~ (3+|t|)^{6/5} on Re s = -1/10;
                 (3Q)^{(1-2 sigma)/2} = (3Q)^{3/5} there.
  C40 L1266, L1282-1284 strip bound |L(sigma+it,psi)| << Q^{3/5}(3+|t|)^2 on [-1/10, 11/10]: PL bookkeeping exact,
                 growth numbers for E_d (sanity only).
  C41 L1273-1277 theta Poisson step: the planar Fourier transform of zbar exp(-pi a|z|^2) is -i a^{-2} xibar exp(-pi|xi|^2/a),
                 which vanishes at xi = 0 (so no constant term, no pole).
  C42 L1087      |1/L^S(s,psi)| <= zeta_K(2) on Re s >= 2.
  C43 L1411-1412 dual profile: MW(1-s) Gamma(s+1/2)/Gamma(3/2-s) (CM) vs MW(1-s) Gamma(s)/Gamma(1-s); both ~ |t|^{2 sigma-1}.
  C44 L1413      the conductor of alpha^e * eta differs from that of eta only at the prime above 3 (eta = chi_.(4d)).
Run time: about 60 s.
"""
import os, sys, math, cmath, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmcommon import *
from cmreport import claim, cited, info, untested
import numpy as np
import mpmath as mp
import zw

# ---------------------------------------------------------------- conductor by brute force
U2 = None          # optional extra finite-order twist chi_.(U2) (U2 supported on lambda), used for C99
def psi_ideal(a, d, e):
    """value at the ideal (a), a coprime to 6d: abar(a')^e chi_{a'}(4d) [chi_{a'}(U2)], a' the primary generator."""
    ap = zw.prim_assoc(a)
    v = abar(ap)**e*chi(ap, (4*d, 0))
    return v*chi(ap, U2) if U2 is not None else v

def is_modulus(g, d, e, B):
    Ng = N(g); R = int(2*math.sqrt(B/Ng)) + 3; tested = 0
    for kx in range(-R, R + 1):
        for ky in range(-R, R + 1):
            a = add((1, 0), mul(g, (kx, ky)))
            na = N(a)
            if na == 0 or na > B or math.gcd(na, 6*abs(d)) != 1:
                continue
            target = abar(a)**e
            if abs(psi_ideal(a, d, e) - target) > 1e-9:
                return False, tested
            tested += 1
    return True, tested

def odd_prime_elems(d):
    out = []
    for p in range(5, abs(d) + 1):
        if abs(d) % p == 0 and all(p % r for r in range(2, int(p**.5) + 1)):
            if p % 3 == 1:
                pis = [q for q in split_primes(p) if N(q) == p]
                out += [(q, 2) for q in pis]
            else:
                out.append(((-p, 0), 1))
    return out

def conductor(d, e, maxnorm=3000):
    gens = [((2, 0), 3), (zw.LAM, 6)] + odd_prime_elems(d)
    cands = []
    for exps in itertools.product(*[range(m + 1) for _, m in gens]):
        g = (1, 0)
        for (q, _), k in zip(gens, exps):
            g = mul(g, pw(q, k))
        if N(g) <= maxnorm:
            cands.append((N(g), g, exps))
    cands.sort()
    valid = []
    for ng, g, exps in cands:
        ok, tested = is_modulus(g, d, e, max(4000, 60*ng))
        if ok and tested >= 30:
            valid.append((ng, g, exps))
    if not valid:
        return None
    f = valid[0]
    # every valid modulus must be a multiple of the minimal one (consistency of the brute force)
    cons = all(all(x >= y for x, y in zip(v[2], f[2])) for v in valid)
    return f, cons, len(valid)

DS = [1, -1, 2, 4, -4, 16, -432, 5, -3, 9, 27, 8, -8]
out = gp("\n".join(f"E=ellinit([0,0,0,0,{d}]); print(ellglobalred(E)[1]); print(ellrootno(E));" for d in DS))
NE = {d: int(out[2*i]) for i, d in enumerate(DS)}; EPS = {d: int(out[2*i + 1]) for i, d in enumerate(DS)}
rows = []; ok37 = True; f37 = True; Qd = {}; FC = {}
for d in DS:
    r = conductor(d, 1)
    if r is None:
        ok37 = False; rows.append(f"d={d}: none"); continue
    (ng, g, exps), cons, nv = r
    Qd[d] = ng; FC[d] = exps
    ok37 &= cons and (3*ng == NE[d])
    f37 &= (4*ng == NE[d])                         # falsifier: |d_K| replaced by 4
    rows.append(f"d={d}:N(f)={ng},N_E={NE[d]}")
claim("C37", 1263, "3Q in Lambda", "E",
      f"3 N(f_psi*) = N(E_d) for {len(DS)} curves (conductor of psi* by brute force over moduli 2^i lambda^j prod p^k)",
      ok37, "; ".join(rows) + f"; falsifier 4N(f) = N_E holds for all? {f37}", falsifier_failed=not f37)

# ---------------------------------------------------------------- C38 FE
mp.mp.dps = 30
pts = ["0.3+2.1*I", "-0.4+5.3*I", "0.15+0.7*I", "1.6-3*I"]
lines = []
for d in [1, 2, 5, 16, -3]:
    lines.append(f"E=ellinit([0,0,0,0,{d}]);" + "".join(f"print(lfun(E,{p}+1/2));print(lfun(E,1-({p})+1/2));" for p in pts))
out = gp("\n".join(lines))
e38 = 0.0; f38a = 1e9; f38b = 1e9; k = 0
for d in [1, 2, 5, 16, -3]:
    Q3 = 3*Qd[d]
    for p in pts:
        s = complex(eval(p.replace("*I", "j")))
        L1 = gpc(out[k]); L2 = gpc(out[k + 1]); k += 2
        lam = lambda s, L, gs, Qn: mp.mpf(Qn)**(s/2)*(2*mp.pi)**(-s)*gs(s)*L
        g_cm = lambda s: mp.gamma(s + mp.mpf(1)/2)
        g_fo = lambda s: mp.gamma(s)
        A = lam(mp.mpc(s), L1, g_cm, Q3); Bv = lam(1 - mp.mpc(s), L2, g_cm, Q3)
        e38 = max(e38, abs(A - EPS[d]*Bv)/abs(A))
        A2 = lam(mp.mpc(s), L1, g_fo, Q3); B2 = lam(1 - mp.mpc(s), L2, g_fo, Q3)
        f38a = min(f38a, abs(A2 - EPS[d]*B2)/abs(A2))
        A3 = lam(mp.mpc(s), L1, g_cm, 4*Qd[d]); B3 = lam(1 - mp.mpc(s), L2, g_cm, 4*Qd[d])
        f38b = min(f38b, abs(A3 - EPS[d]*B3)/abs(A3))
claim("C38", 1265, "functional equation", "E",
      "Lambda(s)=(3Q)^{s/2}(2pi)^{-s}Gamma(s+1/2)L(s,psi) = eps Lambda(1-s), eps=+-1, 5 curves x 4 points (PARI values parsed in double precision)", e38 < 1e-12,
      f"max rel err {float(e38):.1e}; falsifiers Gamma(s) (finite-order factor): min rel err {float(f38a):.2f}; 4Q: {float(f38b):.2f}",
      falsifier_failed=f38a > 1e-3 and f38b > 1e-3)

# ---------------------------------------------------------------- C39
mp.mp.dps = 20
sig = mp.mpf(-1)/10
def gq(t, kind):
    s = mp.mpc(sig, t)
    return abs(mp.gamma(mp.mpf(3)/2 - s)/mp.gamma(mp.mpf(1)/2 + s)) if kind == "cm" else abs(mp.gamma(1 - s)/mp.gamma(s))
ratios = [gq(t, kd)/(3 + t)**(mp.mpf(6)/5) for kd in ("cm", "fo") for t in [0, 1, 5, 20, 100, 1e3, 1e4, 1e5, 1e6]]
slope = [mp.log(gq(1e6, kd)/gq(1e5, kd))/mp.log(10) for kd in ("cm", "fo")]
q35 = (1 - 2*Fr(-1, 10))/2
ok39 = all(abs(sl - mp.mpf(6)/5) < 1e-6 for sl in slope) and max(ratios) < 10 and q35 == Fr(3, 5)
claim("C39", 1280, "exponent 6/5 and Q^{3/5}", "E",
      "|Gamma(3/2-s)/Gamma(1/2+s)| and |Gamma(1-s)/Gamma(s)| ~ (3+|t|)^{6/5} on Re s=-1/10; (1-2sigma)/2 = 3/5", ok39,
      f"fitted exponents {[mp.nstr(x, 8) for x in slope]}; ratio to (3+|t|)^1.2 in [{mp.nstr(min(ratios), 3)},{mp.nstr(max(ratios), 3)}]; "
      f"falsifier exponent 1: deviation {mp.nstr(abs(slope[0] - 1), 3)}", falsifier_failed=abs(slope[0] - 1) > 0.1)

# ---------------------------------------------------------------- C40
# exact PL bookkeeping: 6/5 <= 2, and |s+2| <= (31/30)(3+|t|) for -1/10 <= sigma <= 11/10
ok40 = Fr(6, 5) <= 2 and all(abs(complex(float(sg) + 2, t)) <= 31/30*(3 + t) + 1e-12
                             for sg in (Fr(-1, 10), Fr(11, 10)) for t in np.linspace(0, 1000, 2001))
# sanity numbers with PARI: |L(sigma+it)| / (Q^{3/5}(3+|t|)^2)
gscript = []
for d in [1, 2, 5]:
    gscript.append(f"E=ellinit([0,0,0,0,{d}]);L=lfuninit(E,[1,0.7,105]);m=0;m2=0;"
                   f"for(i=0,24,sg=-0.1+i*0.05;for(t=0,100,v=abs(lfun(L,sg+1/2+I*t))/({Qd[d]}^(3/5)*(3+t)^2);if(v>m,m=v)));"
                   f"for(t=0,100,v=abs(lfun(L,-0.1+1/2+I*t))/({Qd[d]}^(3/5)*(3+t)^(6/5));if(v>m2,m2=v));print(m);print(m2)")
out = gp("\n".join(gscript), prec=19)
nums = [f"d={d}: max|L|/(Q^.6(3+|t|)^2)={float(out[2*i]):.3g}, on Re s=-1/10 /(..)^1.2={float(out[2*i+1]):.3g}"
        for i, d in enumerate([1, 2, 5])]
claim("C40", 1266, "strip growth (PL bookkeeping)", "I",
      "6/5 <= 2 and |s+2| <= (31/30)(3+|t|) on the strip, so PL on Q^{-3/5}L/(s+2)^2 gives Q^{3/5}(3+|t|)^2", ok40,
      "; ".join(nums) + " (sanity: bounded, constants not claimed); falsifier exponent 2 replaced by 1: 6/5<=1 false",
      falsifier_failed=not (Fr(6, 5) <= 1))

# ---------------------------------------------------------------- C41 planar Fourier transform
h = 0.01; X = np.arange(-7, 7, h); XX, YY = np.meshgrid(X, X, indexing="ij")
e41 = 0.0; f41 = 0.0; zero = 0.0
for a in (0.5, 1.0, 2.0):
    g = (XX - 1j*YY)*np.exp(-np.pi*a*(XX**2 + YY**2))
    for xi in [(0.0, 0.0), (0.3, -0.2), (0.7, 0.5), (-0.4, 0.9)]:
        val = np.sum(g*np.exp(-2j*np.pi*(XX*xi[0] + YY*xi[1])))*h*h
        xib = xi[0] - 1j*xi[1]; r2 = xi[0]**2 + xi[1]**2
        pred = -1j*a**-2*xib*np.exp(-np.pi*r2/a)
        if r2 == 0:
            zero = max(zero, abs(val))
        else:
            e41 = max(e41, abs(val - pred)/abs(pred))
            f41 = max(f41, abs(val - (-1j*a**-1*xib*np.exp(-np.pi*r2/a)))/abs(pred))
claim("C41", 1276, "theta: no constant term", "E",
      "FT of zbar e^{-pi a|z|^2} = -i a^{-2} xibar e^{-pi|xi|^2/a} (a=1/2,1,2), zero at xi=0", e41 < 1e-8 and zero < 1e-10,
      f"max rel err {e41:.1e}, |FT(0)| {zero:.1e}; falsifier a^{{-1}} scaling: rel err {f41:.2f}", falsifier_failed=f41 > 0.1)

# ---------------------------------------------------------------- C42
zK2 = float(mp.zeta(2)*(mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9)
from fastgauss import sieve
plist = []
for p in map(int, sieve(200000)):
    if True:
        if p == 3:
            plist.append(3)
        elif p % 3 == 1:
            plist += [p, p]
        elif p*p < 4e10:
            plist.append(p*p)
P = np.array(plist, float)
sup = float(np.prod(1 + P**-2.0))
# actual values for psi*_1 at s = 2 + it (primes of norm <= 2e5 from the helper list, plus lambda, 2 deleted (S))
pr = [p for p in PR_ODD if N(p) <= 200000]
vals = np.array([complex(abar(p)*chi(p, (4, 0))) for p in pr]); norms = np.array([N(p) for p in pr], float)
mx = max(abs(np.prod(1 - vals*norms**(-complex(2, t)))) for t in np.linspace(0, 60, 121))
ok42 = sup < zK2 and mx <= zK2
claim("C42", 1087, "|1/L^S| <= zeta_K(2)", "I",
      "sup over unimodular psi of |prod(1-psi(p)N(p)^-s)| on Re s>=2 is prod(1+N(p)^-2) = zeta_K(2)/zeta_K(4) < zeta_K(2)",
      ok42, f"prod(1+Np^-2) = {sup:.6f} < zeta_K(2) = {zK2:.6f}; max |1/L^S(2+it, psi*_1)| = {mx:.6f}; falsifier bound 1: violated",
      falsifier_failed=sup > 1 and mx > 1)

# ---------------------------------------------------------------- C43
mp.mp.dps = 20
sl = []; fs = []
for sg in [-0.5, 0.0, 0.3, 1.0, 2.0]:
    for kind in ("cm", "fo"):
        def q(t):
            s = mp.mpc(sg, t)
            return abs(mp.gamma(s + mp.mpf(1)/2)/mp.gamma(mp.mpf(3)/2 - s)) if kind == "cm" else abs(mp.gamma(s)/mp.gamma(1 - s))
        e = mp.log(q(1e6)/q(1e5))/mp.log(10)
        sl.append(abs(e - (2*sg - 1))); fs.append(abs(e - 2*sg))
# FE ratio form L(s)/L(1-s) = eps (3Q)^{1/2-s}(2pi)^{2s-1}Gamma(3/2-s)/Gamma(s+1/2), tested at d=1 via PARI (same data as C38)
out = gp("E=ellinit([0,0,0,0,1]);print(lfun(E,0.3+2.1*I+1/2));print(lfun(E,1-(0.3+2.1*I)+1/2))")
s = mp.mpc(0.3, 2.1); L1 = gpc(out[0]); L2 = gpc(out[1])
rhs = EPS[1]*(3*Qd[1])**(mp.mpf(1)/2 - s)*(2*mp.pi)**(2*s - 1)*mp.gamma(mp.mpf(3)/2 - s)/mp.gamma(s + mp.mpf(1)/2)
erat = abs(L1/L2 - rhs)/abs(rhs)
claim("C43", 1411, "dual profile", "E",
      "L(s)/L(1-s) = eps(3Q)^{1/2-s}(2pi)^{2s-1}Gamma(3/2-s)/Gamma(s+1/2), so MW#(s)=MW(1-s)Gamma(s+1/2)/Gamma(3/2-s); "
      "both quotients ~|t|^{2 sigma-1}", max(sl) < 1e-6 and erat < 1e-12,
      f"FE-ratio rel err {float(erat):.1e}; exponent deviations max {float(max(sl)):.1e}; falsifier |t|^(2 sigma): min dev {float(min(fs)):.2f}",
      falsifier_failed=min(fs) > 0.5)

# ---------------------------------------------------------------- C44
rows = []; ok44 = True; f44 = True
for d in [1, -1, 2, 5, 16, -3, 4]:
    r0 = conductor(d, 0); r1 = conductor(d, 1); r2 = conductor(d, -1)
    e0 = r0[0][2]; e1 = r1[0][2]; e2 = r2[0][2]
    same_away = all(e0[i] == e1[i] == e2[i] for i in range(len(e0)) if i != 1)     # index 1 = lambda
    ok44 &= same_away
    f44 &= (e0 == e1)                                                            # falsifier: conductors equal everywhere
    rows.append(f"d={d}: eta {e0}, abar*eta {e1}, alpha*eta {e2} (exponents of 2, lambda, p|d)")
claim("C44", 1413, "conductor change only at lambda", "E",
      "conductors of eta=chi_.(4d), abar*eta and alpha*eta agree away from lambda", ok44,
      "; ".join(rows) + f"; falsifier 'identical conductors' holds for all? {f44}", falsifier_failed=not f44)

# ---------------------------------------------------------------- C99 closure under finite-order twists (L1065-1067)
# psi*_d times theta = chi_.(lambda) (finite order, ramified only at lambda and 2): find a modulus g with
# (psi theta)((a)) = abar(a)^k for all a = 1 mod g.  k = 1 must succeed (infinity type zbar/|z|); the falsifier k = 0
# (i.e. "psi theta has finite order") must find no modulus.
def find_modulus(val, k, d, maxnorm=3000):
    gens = [((2, 0), 3), (zw.LAM, 6)] + odd_prime_elems(d)
    cands = []
    for exps in itertools.product(*[range(m + 1) for _, m in gens]):
        g = (1, 0)
        for (q, _), kk in zip(gens, exps):
            g = mul(g, pw(q, kk))
        if N(g) <= maxnorm:
            cands.append((N(g), g, exps))
    for ng, g, exps in sorted(cands):
        B = max(4000, 60*ng); Rr = int(2*math.sqrt(B/ng)) + 3; ok = True; tested = 0
        for kx in range(-Rr, Rr + 1):
            for ky in range(-Rr, Rr + 1):
                a = add((1, 0), mul(g, (kx, ky))); na = N(a)
                if na == 0 or na > B or math.gcd(na, 6*abs(d)) != 1:
                    continue
                if abs(val(a) - abar(a)**k) > 1e-9:
                    ok = False; break
                tested += 1
            if not ok:
                break
        if ok and tested >= 30:
            return ng, exps
    return None
rows = []; ok99 = True; f99 = True
for d in [1, -1, 4]:
    val = lambda a, d=d: (lambda ap: abar(ap)*chi(ap, (4*d, 0))*chi(ap, zw.LAM))(zw.prim_assoc(a))
    r1 = find_modulus(val, 1, d); r0 = find_modulus(val, 0, d)
    ok99 &= r1 is not None
    f99 &= r0 is None
    rows.append(f"d={d}: modulus for psi*theta with abar: {r1}; with trivial infinity type: {r0}")
claim("C99", 1065, "closure under finite-order twists", "E",
      "psi*_d * chi_.(lambda) has infinity type zbar/|z| (a modulus g with value abar(a) on a = 1 mod g exists), "
      "so its primitive inducer lies in F_CM", ok99,
      "; ".join(rows) + " (exponents of 2, lambda); falsifier 'finite order' finds no modulus", falsifier_failed=f99)

# ---------------------------------------------------------------- C101 theta series (L1274-1275)
# phi(a) = psi((a))/abar(a) on elements; phi(zeta a) conj(zeta) = phi(a); sum_{a in O} phi(a) abar(a)|a| N(a)^{-s}
# = 6 sum_{n primary} psi(n) N(n)^{1/2-s}:  the theta series is that of L(s,psi), shifted by 1/2.
d = 1; s_ = complex(2.5, 1.3)
R = 45; lhs = 0; rhs = 0; e101 = 0.0; cnt = 0; lhs_f = 0
for x in range(-R, R + 1):
    for y in range(-R, R + 1):
        a = (x, y); na = N(a)
        if na == 0 or na > 1500 or math.gcd(na, 6) != 1:
            continue
        ap = zw.prim_assoc(a)
        psi = abar(ap)*chi(ap, (4*d, 0))
        phi = psi/abar(a)
        for z in UNITS:
            za = mul(z, a)
            e101 = max(e101, abs(psi/abar(za)*cx(z).conjugate() - phi)); cnt += 1
        lhs += phi*cx(a).conjugate()*na**(-s_)
        lhs_f += psi*cx(a).conjugate()*na**(-s_)          # falsifier: phi(a) := psi((a)) (no unit compensation)
        if primary(a):
            rhs += 6*psi*na**(0.5 - s_)
ok101 = e101 < 1e-12 and abs(lhs - rhs) < 1e-12*abs(rhs)
claim("C101", 1274, "theta series of psi", "E",
      "phi(zeta a)zetabar = phi(a) for phi = psi/abar; sum_{a in O} phi(a) abar(a) N(a)^-s = 6 sum_n psi(n) N(n)^{1/2-s} (N<=1500)",
      ok101, f"unit relation max err {e101:.1e} ({cnt} cases); series rel err {abs(lhs - rhs)/abs(rhs):.1e}; "
      f"falsifier phi = psi((a)) (no compensation): series sums to {abs(lhs_f):.1e} vs {abs(rhs):.2f}",
      falsifier_failed=abs(lhs_f) < 1e-6*abs(rhs))

# ---------------------------------------------------------------- C102 conductor bound Q_psi << N(u) (L1307-1308), INFO
info("C102", 1308, "Q_psi << N(u)", "A", "conductor norm of psi*_d (row u = (4d)^5 ~ 4d) relative to N(4d)",
     "; ".join(f"d={d}: N(f)/N(4d) = {Qd[d]/(16*d*d):.3g}" for d in DS) + " (bounded: the extra factors are supported on 2, 3)")
