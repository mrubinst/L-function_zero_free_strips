"""T1. Theta transformation, functional equation and root number for CM characters.

Claims tested (insert, Lemma s2-fe):
 (a) Theta_phi(v) = - gamma(phi)/(f sqrt Q) * v^{-2} * conj(Theta_phi(1/(Q v)))       (v>0)
     where Theta_phi(v) = sum_{x in O} phi(x) conj(x) exp(-2 pi v N x / sqrt3).
 (b) Lambda(s,psi) = (3Q)^{s/2} (2pi)^{-s} Gamma(s+1/2) L(s,psi) satisfies
     Lambda(s,psi) = W(psi) Lambda(1-s, conj psi),  W(psi) = -gamma(phi) sqrt(Q)/f, |W|=1.
     Tested via the split-integral formula, which is independent of the split point v0
     iff (a) holds, and which must reproduce the Dirichlet series at s=3.
 (c) W(alphabar) = +1 (root number of 27a = y^2=x^3+16) and the coefficients of L(s,alphabar)
     are a_p(E_16)/sqrt p.
Falsifiers: drop the conjugation, use conj(gamma), flip the sign, use conj(f).
"""
import math
import cmath
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from cmhecke import *

mp.mp.dps = 30
C2 = 2 * math.pi / math.sqrt(3)
results = []
FALS_SOME, FALS_ALL = {}, {}   # falsifier errors (added 2026-10-12): must exceed 1e-6 for some / for every character
_FORCE = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]   # self-test (added 2026-10-12): first check counts as failed


def rep(name, err, tol):
    ok = err < tol
    if _FORCE[0]:
        _FORCE[0] = False; ok = False; name = "[CM_TESTS_FORCE_FAIL=1] " + name
    results.append((name, ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: err={err:.3e} (tol {tol:g})")


def theta(ch, v, pts):
    s = 0
    for x in pts:
        ph = ch.phi(x)
        if ph != 0:
            z = cval(x)
            s += ph * z.conjugate() * math.exp(-C2 * v * norm(x))
    return s


def main():
    chars = make_test_chars()
    pts = lattice(1200)
    for ch in chars:
        Q = ch.Q
        g = ch.gauss()
        fc = cval(ch.f)
        eps0 = g / (fc * math.sqrt(Q))
        rep(f"[{ch.name}] |gamma(phi)| = 1", abs(abs(g) - 1), 1e-12)
        rep(f"[{ch.name}] |W| = 1", abs(abs(ch.root_number()) - 1), 1e-12)
        worst = 0
        worst_f = {"no conj": 1e9, "conj gamma": 1e9, "+ sign": 1e9, "conj f": 1e9}
        for c in [0.35, 0.6, 1.0, 1.7]:
            v = c / math.sqrt(Q)
            lhs = theta(ch, v, pts)
            th2 = theta(ch, 1 / (Q * v), pts)
            rhs = -eps0 * v ** -2 * th2.conjugate()
            scale = max(abs(lhs), 1e-300)
            worst = max(worst, abs(lhs - rhs) / scale)
            fals = {
                "no conj": -eps0 * v ** -2 * th2,
                "conj gamma": -(g.conjugate() / (fc * math.sqrt(Q))) * v ** -2 * th2.conjugate(),
                "+ sign": +eps0 * v ** -2 * th2.conjugate(),
                "conj f": -(g / (fc.conjugate() * math.sqrt(Q))) * v ** -2 * th2.conjugate(),
            }
            for k, val in fals.items():
                worst_f[k] = min(worst_f[k], abs(lhs - val) / scale)
        rep(f"[{ch.name}] theta transformation (rel)", worst, 1e-11)
        for k, val in worst_f.items():
            # falsifier must FAIL: its best error must be large.  (For a real character some
            # falsifiers coincide with the truth; then we only report.)
            print(f"      falsifier '{k}': min rel err = {val:.3e}")
            FALS_SOME.setdefault(f"theta '{k}'", []).append(val)

    # (b) functional equation via the split formula
    pts2 = lattice(400)
    for ch in chars:
        Q = ch.Q
        # everything in mpmath precision (phi takes sixth roots of unity; x = a + b*omega),
        # since at height |t| ~ 11 dividing by Gamma(s+1/2) amplifies double-precision errors.
        OM = mp.exp(2j * mp.pi / 3)
        W6 = mp.exp(2j * mp.pi / 6)

        def root6(c):
            k = round(cmath.phase(c) / (math.pi / 3)) % 6
            assert abs(c - cmath.exp(1j * math.pi * k / 3)) < 1e-9
            return W6 ** k
        fm = ch.f[0] + ch.f[1] * OM
        gm = mp.mpc(0)
        for x in ch.reps:
            val = ch.phi(x)
            if val != 0:
                zz = (x[0] + x[1] * OM) / fm
                gm += root6(val) * mp.exp(4j * mp.pi * mp.im(zz) / mp.sqrt(3))
        gm /= mp.sqrt(Q)
        eps0 = gm / (fm * mp.sqrt(Q))
        data = [(x, root6(ch.phi(x))) for x in pts2 if ch.phi(x) != 0]
        c2 = 2 * mp.pi / mp.sqrt(3)

        def L_split(s, v0, sign=-1, conj_dual=True):
            s = mp.mpc(s)
            v0 = mp.mpf(v0)
            A = mp.mpc(0)
            B = mp.mpc(0)
            for x, ph in data:
                n = norm(x)
                z = x[0] + x[1] * OM
                A += ph * mp.conj(z) * (c2 * n) ** (-(s + 0.5)) * mp.gammainc(s + 0.5, c2 * n * v0)
                dual = (mp.conj(ph) * z) if conj_dual else (ph * mp.conj(z))
                B += dual * (c2 * n) ** (-(1.5 - s)) * mp.gammainc(1.5 - s, c2 * n / (Q * v0))
            I = A + sign * eps0 * mp.mpf(Q) ** (1.5 - s) * B
            return I / (6 * c2 ** (-(s + 0.5)) * mp.gamma(s + 0.5))

        # Dirichlet series at s=3 directly
        big = lattice(6000)
        Ldir = 0
        for x in big:
            v = ch.psi_elem(x)
            if v != 0:
                Ldir += v * norm(x) ** -3.0
        Ldir /= 6
        L3 = L_split(3, 1 / math.sqrt(Q))
        rep(f"[{ch.name}] split formula = Dirichlet series at s=3", abs(complex(L3) - Ldir), 2e-7)
        for s in [0.5 + 3j, 0.2 - 7j, -0.1 + 11j, 0.8 + 0.3j]:
            l1 = L_split(s, mp.mpf(0.7) / mp.sqrt(Q))
            l2 = L_split(s, mp.mpf(1.45) / mp.sqrt(Q))
            rel = abs(l1 - l2) / abs(l1)
            rep(f"[{ch.name}] split-point independence at s={s}", float(rel), 1e-12)
        # falsifier: wrong sign of W -> split point dependence
        s = 0.5 + 3j
        l1 = L_split(s, 0.7 / math.sqrt(Q), sign=+1)
        l2 = L_split(s, 1.45 / math.sqrt(Q), sign=+1)
        fw = float(abs(l1-l2)/abs(l1))
        print(f"      falsifier 'W -> -W': split dependence rel = {fw:.3e} (must be large)")
        FALS_ALL.setdefault("split 'W -> -W'", []).append(fw)
        l1 = L_split(s, 0.7 / math.sqrt(Q), conj_dual=False)
        l2 = L_split(s, 1.45 / math.sqrt(Q), conj_dual=False)
        fd = float(abs(l1-l2)/abs(l1))
        print(f"      falsifier 'dual = psi instead of conj psi': rel = {fd:.3e}")
        FALS_SOME.setdefault("split 'dual = psi'", []).append(fd)

    # (c) alphabar and 27a
    chA = chars[0]
    rep("W(alphabar) = +1 (root number of 27a)", abs(chA.root_number() - 1), 1e-12)
    bad = 0
    count = 0
    for p in range(5, 400):
        if any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
            continue
        # a_p of y^2 = x^3 + 16
        npts = 1 + sum(1 for x in range(p) for y in range(p) if (y * y - x ** 3 - 16) % p == 0)
        ap = p + 1 - npts
        if p % 3 == 1:
            pi = primary_primes_of_norm(p)[0]
            coef = (chA.psi_elem(pi) + chA.psi_elem(conj(pi))) * math.sqrt(p)
            bad = max(bad, abs(coef - ap))
        else:
            # inert: primary generator of (p) is -p; Euler factor (1 - psi(-p) p^{-2s})^{-1}
            v = chA.psi_elem((-p, 0))
            bad = max(bad, abs(v - (-1)) + abs(ap))
        count += 1
    rep(f"L(s,alphabar) = L(E_16, s+1/2) coefficientwise at {count} primes 5<=p<400", bad, 1e-9)
    # Falsifier verdict (added 2026-10-12).  For the real character alphabar several wrong variants coincide with the
    # truth (see above), so each variant must fail (rel err > 1e-6) for at least one character; 'W -> -W' for all.
    fals_ok = all(max(v) > 1e-6 for v in FALS_SOME.values()) and all(min(v) > 1e-6 for v in FALS_ALL.values())
    results.append(("falsifiers", fals_ok))
    print("falsifiers: each wrong variant fails (rel err > 1e-6) for some character, 'W -> -W' for every character:",
          "yes" if fals_ok else "NO")
    print("ALL PASS" if all(ok for _, ok in results) else "SOME FAILED")
    raise SystemExit(0 if all(ok for _, ok in results) else 1)


if __name__ == "__main__":
    main()
