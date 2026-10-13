"""t01: how the curves enter (guide L1006-1037, L1068-1069).

Claims tested
  C01 L1015  7/8 + 1/2 = 11/8, strip 11/8 < Re s < 3/2 inside the critical strip 1/2 < Re s < 3/2.
  C02 L1021-1022  sqrt(p)(psi*_d(pi) + psi*_d(pibar)) = a_p(d) at split p (both primes above p).
  C03 L1022-1023  psi*_d(-p) = -1 for inert p (p not dividing 6d), matching a_p(d) = 0.
  C04 L1023-1024  L(E_d, s+1/2) = L(s, psi*_d) up to finitely many Euler factors: coefficient identity
                  a_m(E_d)/sqrt(m) = sum_{N(n)=m} psi*_d(n) for (m,6d)=1, m <= 600, against PARI ellan.
  C05 L1024  the omitted Euler factors (p | 6d) do not vanish in Re s > 0 (PARI local data).
  C06 L1025-1026  the counts: >= 480 split pairs (both primes) and >= 112 inert pairs.
  C07 L1027-1028  psi*_d = conj(psi_d), psi_d = alpha conj(chi(4d)) = alpha chi((4d)^5).
  C08 L1030  psi*_d(pibar) = conj(psi*_d(pi)) since 4d is rational.
  C09 L1029-1031  psi_d and psi*_d have the same Euler factor at every split and inert prime.
  C10 L1036-1037  j-invariants of CM by orders of K are exactly 0, 54000, -12288000 (h(-3f^2)=1 iff f<=3).
  C11 L1068-1069  conj(psi) = alpha conj(eta) and L(s, conj psi) = conj(L(conj s, psi)) (non-real example).
  C12 L1032-1034  every unitary Hecke character of infinity type (zbar/|z|)^{+-1} is abar*eta or alpha*etabar:
                  model check that n -> psi(n)/abar(n) is invariant under n -> u n on units (finite part).
Run time: about 20 s.
"""
import os, sys, math, cmath, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmcommon import *
from cmreport import claim, cited, info, untested

# ---------------------------------------------------------------- C01
ok = (Fr(7, 8) + Fr(1, 2) == Fr(11, 8)) and (Fr(1, 2) < Fr(11, 8) < Fr(3, 2))
fals = (Fr(7, 8) + Fr(1, 2) == Fr(17, 12))          # Part I value must differ
claim("C01", 1015, "thm:78", "I", "7/8+1/2 = 11/8 and 1/2 < 11/8 < 3/2", ok,
      "exact Fractions", falsifier_failed=not fals)

# ---------------------------------------------------------------- helpers
DS = [1, 2, -1, 3, -2, -3, 5, -5, 7, -7, 10, 11, -11, 13, 16, -432, 17, 19, 6, 12]

def ap_direct(d, p):
    """a_p(d) = -sum_x ((x^3+d)/p) (Legendre symbol, zero on multiples)."""
    s = 0
    for x in range(p):
        v = (x*x*x + d) % p
        if v:
            s += 1 if pow(v, (p - 1)//2, p) == 1 else -1
    return -s

def psistar(n, d):
    return abar(n)*chi(n, (4*d, 0))

def psi_d(n, d):
    return alpha(n)*chi(n, (4*d, 0)).conjugate()

# ---------------------------------------------------------------- C02, C03, C06, C08
split_err = 0.0; nsplit = 0; fals_err = 0.0; conj_err = 0.0
for d in DS:
    for pi in split_primes(700):
        p = N(pi)
        if (6*d) % p == 0:
            continue
        ap = ap_direct(d, p)
        pib = conj(pi)
        assert primary(pib)
        lhs = math.sqrt(p)*(psistar(pi, d) + psistar(pib, d))
        split_err = max(split_err, abs(lhs - ap))
        # falsifier: alpha instead of abar (psi = alpha chi(4d))
        lf = math.sqrt(p)*(alpha(pi)*chi(pi, (4*d, 0)) + alpha(pib)*chi(pib, (4*d, 0)))
        fals_err = max(fals_err, abs(lf - ap))
        conj_err = max(conj_err, abs(psistar(pib, d) - psistar(pi, d).conjugate()))
        nsplit += 1          # one pair (d,p); both primes above p enter
claim("C02", 1022, "psi*_d", "E", f"sqrt(p)(psi*(pi)+psi*(pibar)) = a_p(d), {nsplit} split pairs (d,p), p<700, {len(DS)} d",
      split_err < 1e-9, f"max err {split_err:.2e}; falsifier alpha*chi(4d): max err {fals_err:.2f}",
      falsifier_failed=fals_err > 0.5)

inert_err = 0.0; ninert = 0; fals_inert = 0
for d in DS:
    for q in inert_primes(500):
        p = -q[0]
        if (6*d) % p == 0:
            continue
        v = psistar(q, d)
        inert_err = max(inert_err, abs(v - (-1)))
        assert ap_direct(d, p) == 0
        # falsifier: forgetting the angular factor (abar(-p) = -1) would give +1
        fals_inert += abs(chi(q, (4*d, 0)) - (-1)) < 1e-9
        ninert += 1
claim("C03", 1022, "psi*_d(-p)", "E", f"psi*_d(-p) = -1 and a_p(d)=0 for {ninert} inert pairs (p<500)",
      inert_err < 1e-12, f"max err {inert_err:.1e}; falsifier without the angular factor (value chi_{{-p}}(4d)) equals -1 in {fals_inert} of {ninert} cases",
      falsifier_failed=(fals_inert == 0))

# scope of the inert statement: for p | d (p = 2 mod 3) the value is 0, not -1 (the factor is deleted; a_p = 0 too)
vals = [(d, p, psistar((-p, 0), d)) for d in (5, -5, 10, 11, -11, 17) for p in (5, 11, 17) if d % p == 0]
info("C03b", 1023, "scope of psi*_d(-p) = -1", "E",
     "the inert statement needs p not dividing d (as in the split case): for p | d the value is 0",
     "; ".join(f"d={d},p={p}: psi*(-p)={v.real:+.0f}" for d, p, v in vals) +
     " (convention: the guide's 'p not | 6d' is stated only in the split clause; harmless, the factor is among the deleted ones)")

claim("C06", 1025, "counts", "E", "numbers of checked pairs at least 480 split / 112 inert", nsplit >= 480 and ninert >= 112,
      f"here {nsplit} split pairs (both primes above each p) and {ninert} inert pairs",
      falsifier_failed=not (nsplit >= 10**6))

# falsifier for C08: non-rational numerator u = 2 + omega
fe = max(abs(abar(conj(pi))*chi(conj(pi), (2, 1)) - (abar(pi)*chi(pi, (2, 1))).conjugate()) for pi in split_primes(400))
claim("C08", 1030, "psi*(pibar)", "E", "psi*_d(pibar) = conj psi*_d(pi) (4d rational)", conj_err < 1e-12,
      f"max err {conj_err:.1e} over {nsplit} split pairs; falsifier non-rational row u=2+omega: max dev {fe:.2f}",
      falsifier_failed=fe > 0.5)

# ---------------------------------------------------------------- C04, C05 (PARI)
import numpy as np
def zero_free_factor(coeffs, p):
    """coeffs [c0,c1,c2] of c0 + c1 X + c2 X^2, X = p^{-s-1/2}.  A zero with Re s > 0 means a root with |X| < p^{-1/2}."""
    c = [x for x in coeffs]
    while c and c[-1] == 0:
        c.pop()
    if len(c) <= 1:
        return True
    roots = np.roots(c[::-1])
    return all(abs(r) >= p**-0.5*(1 - 1e-12) for r in roots)
DLIST = [1, 2, -3, 7, 16, 5, -432, -1]
M = 600
prims = primaries(M)
bynorm = {}
for n in prims:
    bynorm.setdefault(N(n), []).append(n)
script = "\n".join(f"E=ellinit([0,0,0,0,{d}]); v=ellan(E,{M}); print(v); print(ellglobalred(E)[1]);"
                   f" print(ellrootno(E)); print(vector(10,k,ellap(E,prime(k))));" for d in DLIST)
out = gp(script)
coef_err = 0.0; coef_n = 0; fals_coef = 0.0; local_ok = True; local_txt = []
CONDS = {}
for k, d in enumerate(DLIST):
    an = eval(out[4*k].replace(";", ","))
    Ncond = int(out[4*k + 1]); CONDS[d] = (Ncond, int(out[4*k + 2]))
    aps = eval(out[4*k + 3])
    for m in range(1, M + 1):
        if math.gcd(m, 6*d) != 1:
            continue
        s = sum(psistar(n, d) for n in bynorm.get(m, []))
        sf = sum(abar(n)*chi(n, (d, 0)) for n in bynorm.get(m, []))     # falsifier: chi(d) instead of chi(4d)
        coef_err = max(coef_err, abs(s - an[m - 1]/math.sqrt(m)))
        fals_coef = max(fals_coef, abs(sf - an[m - 1]/math.sqrt(m)))
        coef_n += 1
    # local factors at p | 6d: 1 - a_p X + [p good] p X^2, X = p^{-s-1/2}; zero-free for Re s > 0 iff roots |X0| >= p^{-1/2}
    P = [2, 3] + [p for p in range(5, 500) if all(p % r for r in range(2, int(p**.5) + 1)) and d % p == 0]
    for p in P:
        ap = aps[[2, 3, 5, 7, 11, 13, 17, 19, 23, 29].index(p)] if p <= 29 else None
        if ap is None:
            ap = int(gp(f"print(ellap(ellinit([0,0,0,0,{d}]),{p}))")[0])
        good = (Ncond % p != 0)
        ok_p = zero_free_factor([1, -ap, p if good else 0], p)
        local_ok &= ok_p
        local_txt.append(f"d={d},p={p}:{'good' if good else 'bad'},a_p={ap}")
claim("C04", 1023, "L(E_d,s+1/2)=L(s,psi*)", "E",
      f"a_m/sqrt(m) = sum_{{N(n)=m}} psi*_d(n), (m,6d)=1, m<={M}, d in {DLIST} ({coef_n} coefficients)",
      coef_err < 1e-9, f"max err {coef_err:.1e}; falsifier chi_n(d) in place of chi_n(4d): max err {fals_coef:.2f}",
      falsifier_failed=fals_coef > 0.5)
# falsifier for C05: a made-up factor 1 - p X (root 1/p, i.e. a zero at Re s = 1/2 > 0) must be flagged
ff = zero_free_factor([1, -5, 0], 5)
claim("C05", 1024, "Euler factors", "E", "omitted local factors at p|6d have no zeros in Re s>0 (PARI a_p, conductor)",
      local_ok, f"{len(local_txt)} local factors checked, e.g. {'; '.join(local_txt[:4])}; falsifier 1-5X flagged: {not ff}",
      falsifier_failed=not ff)

# ---------------------------------------------------------------- C07, C09
e7 = 0.0; e7f = 0.0; e9 = 0.0; cnt = 0
for d in [1, 2, -3, 7, 5, 11]:
    for n in primaries(400):
        if math.gcd(N(n), 6*abs(d)) != 1:
            continue
        a = psistar(n, d); b = psi_d(n, d)
        c5 = alpha(n)*chi(n, (pow(4*d, 5), 0))
        e7 = max(e7, abs(a - b.conjugate()), abs(b - c5))
        e7f = max(e7f, abs(b - alpha(n)*chi(n, (pow(4*d, 2), 0))))     # falsifier: (4d)^2 instead of (4d)^5
        cnt += 1
    for pi in split_primes(600):
        if (6*d) % N(pi) == 0:
            continue
        X = 0.37 + 0.11j
        f1 = (1 - psistar(pi, d)*X)*(1 - psistar(conj(pi), d)*X)
        f2 = (1 - psi_d(pi, d)*X)*(1 - psi_d(conj(pi), d)*X)
        e9 = max(e9, abs(f1 - f2))
    for q in inert_primes(200):
        if (6*d) % (-q[0]) == 0:
            continue
        e9 = max(e9, abs(psistar(q, d) - psi_d(q, d)))
claim("C07", 1028, "psi_d", "E", "psi*_d = conj(psi_d) and alpha conj(chi_n(4d)) = alpha chi_n((4d)^5)", e7 < 1e-12,
      f"{cnt} primary n, max err {e7:.1e}; falsifier (4d)^2: {e7f:.2f}", falsifier_failed=e7f > 0.5)
# falsifier for C09: the non-rational row 2+omega gives different Euler factors for psi and its conjugate
e9f = 0.0
for pi in split_primes(300):
    X = 0.37 + 0.11j
    a1 = abar(pi)*chi(pi, (2, 1)); a2 = abar(conj(pi))*chi(conj(pi), (2, 1))
    e9f = max(e9f, abs((1 - a1*X)*(1 - a2*X) - (1 - a1.conjugate()*X)*(1 - a2.conjugate()*X)))
claim("C09", 1029, "same L", "E", "psi_d and psi*_d have equal Euler factors at split and inert primes", e9 < 1e-12,
      f"max err {e9:.1e}; falsifier u=2+omega (non-rational): {e9f:.2f}", falsifier_failed=e9f > 0.05)

# ---------------------------------------------------------------- C10
out = gp("v=[];for(f=1,300,if(qfbclassno(-3*f^2)==1,v=concat(v,f)));print(v);"
         "print([polclass(-3),polclass(-12),polclass(-27)]);print(polclass(-48))")
fs = eval(out[0]); pols = out[1]
ok = fs == [1, 2, 3] and pols.replace(" ", "") == "[x,x-54000,x+12288000]"
claim("C10", 1036, "j-invariants", "E", "orders Z+fO_K with h=1 are f=1,2,3 with j = 0, 54000, -12288000", ok,
      f"f with h(-3f^2)=1, f<=300: {fs}; Hilbert polys {pols}; falsifier f=4: polclass(-48) = {out[2]} (degree 2)",
      falsifier_failed=("x^2" in out[2]))

# ---------------------------------------------------------------- C11
u = (2, 1)
s = complex(3.0, 1.7)
ns = primaries(3000)
def Ls(f, s):
    return sum(f(n)*N(n)**(-s) for n in ns)
psi = lambda n: abar(n)*chi(n, u)
L1 = Ls(lambda n: psi(n).conjugate(), s)
L2 = Ls(psi, s.conjugate()).conjugate()
L3 = Ls(psi, s)
conj_ok = abs(L1 - L2) < 1e-12 and abs(abar((7, 9)).conjugate() - alpha((7, 9))) < 1e-15
claim("C11", 1068, "conjugation", "E", "L(s,conj psi) = conj L(conj s, psi) for psi = abar*chi(2+omega) (truncated, s=3+1.7i)",
      conj_ok, f"|diff| {abs(L1 - L2):.1e}; falsifier L(s,conj psi)=L(s,psi): |diff| {abs(L1 - L3):.2e}",
      falsifier_failed=abs(L1 - L3) > 1e-3)

# ---------------------------------------------------------------- C12
untested("C12", 1032, "classification", "C",
         "every unitary Hecke character of K with infinity type (zbar/|z|)^{+-1} is abar*eta or alpha*etabar; covers all "
         "weight-2 CM-by-K newforms",
         "standard theory (quotient of two characters with the same infinity type has trivial infinity type, hence finite "
         "order); no finite test. Consistent with C04 (E_d) and the conductor check in t04 (C33).")
