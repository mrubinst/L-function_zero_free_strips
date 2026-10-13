"""T2. Gamma-quotient (Stirling) exponents and the explicit constants of the strip bound.

Claims (insert, Lemma s2-gamma and Lemma cm-growth):
 (a) G_sigma(t) := |Gamma(3/2-sigma-it)/Gamma(1/2+sigma+it)| satisfies
        c1 (3+|t|)^{1-2 sigma} <= G_sigma(t) <= c2 (3+|t|)^{1-2 sigma}
     uniformly for sigma in a compact subset of (-1/2, 3/2); in particular exponent 6/5 on sigma=-1/10,
     the same as |Gamma(1-s)/Gamma(s)| (finite order, OAI78) on that line, and exponent
     2a-1+12e on the reflected line sigma = 1-a-6e of OAI78 lem:buffered-bins.
 (b) explicit constants: with C_* = max(C_left, C_right),
        |L(s,psi)| <= C_* Q^{3/5} |s+2|^2   (-1/10 <= sigma <= 11/10), all psi in F_CM,
     C_left = 3^{3/5}(2pi)^{-6/5} zeta_K(11/10) sup_t G_{-1/10}(t)/(1.9^2+t^2),
     C_right = zeta_K(11/10)/3.1^2.
 (c) actual values of L(s,psi) for the four test characters respect (b) on a grid.
Falsifiers: exponent (1-2 sigma) +- 0.1 must drift; the bound with Q^{3/5} replaced by Q^{1/5}
is not a falsifier of a true statement, so instead we check that the left-line identity
|L(-1/10+it)| = (3Q)^{3/5}(2pi)^{-6/5} G(t) |L(11/10-it, conj psi)| holds (FE) and fails with
exponent 1/2 in place of 3/5.
"""
import math
import cmath
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from cmhecke import *

results = []
FALS = []   # (name, error) of wrong variants that must fail with error > 1e-6 (added 2026-10-12)
_FORCE = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]   # self-test (added 2026-10-12): first check counts as failed


def rep(name, ok, detail=""):
    if _FORCE[0]:
        _FORCE[0] = False; ok = False; name = "[CM_TESTS_FORCE_FAIL=1] " + name
    results.append((name, ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name} {detail}")


def G(sig, t):
    return abs(mp.gamma(mp.mpc(1.5 - sig, -t)) / mp.gamma(mp.mpc(0.5 + sig, t)))


def Gfin(sig, t):
    return abs(mp.gamma(mp.mpc(1 - sig, -t)) / mp.gamma(mp.mpc(sig, t)))


def main():
    mp.mp.dps = 30
    ts = [0] + [10 ** (k / 8) for k in range(-8, 33)]  # 0.1 .. 1e4
    # (a) exponents
    for sig in [-0.25, -0.1, -0.006, 0.0, 0.25, 0.49, 0.5, 0.75, 1.1, 1.25]:
        r = [float(G(sig, t) / (3 + t) ** (1 - 2 * sig)) for t in ts]
        lo, hi = min(r), max(r)
        tail = float(G(sig, 1e4) / mp.mpf(1e4) ** (1 - 2 * sig))
        fin = float(G(sig, 1e4) / Gfin(sig, 1e4)) if sig > 0 else float('nan')
        rep(f"sigma={sig:+.3f}: G/(3+|t|)^(1-2s) in [{lo:.4f},{hi:.4f}]; G/|t|^(1-2s) at t=1e4: {tail:.6f}",
            0 < lo and hi < float('inf') and abs(tail - 1) < 1e-3,
            f"; CM/finite-order quotient at t=1e4: {fin:.6f}" if sig > 0 else "")
        for d in (+0.1, -0.1):
            r1 = float(G(sig, 1e2) / (3 + 1e2) ** (1 - 2 * sig + d))
            r2 = float(G(sig, 1e4) / (3 + 1e4) ** (1 - 2 * sig + d))
            if not (r2 / r1 > 1.5 or r2 / r1 < 1 / 1.5):
                rep(f"   falsifier exponent shift {d:+} at sigma={sig} did not drift", False)
    print("   (falsifiers: exponent (1-2 sigma)+-0.1 drift by a factor >1.5 between t=1e2 and 1e4: checked)")
    # finite-order comparison on sigma=-1/10 (OAI78 exponent 6/5)
    r = float(G(-0.1, 1e4) / Gfin(-0.1, 1e4))
    rep(f"on sigma=-1/10, CM quotient / finite-order quotient at t=1e4 = {r:.6f} (-> 1)", abs(r - 1) < 1e-3)

    # (b) explicit constants
    zK = mp.zeta(1.1) * (mp.zeta(1.1, 1 / 3.) - mp.zeta(1.1, 2 / 3.)) / mp.mpf(3) ** 1.1
    # L(s,chi_{-3}) = 3^{-s}(zeta(s,1/3)-zeta(s,2/3))
    tt = [k / 20 for k in range(0, 4001)] + [200 + k for k in range(0, 9801, 7)]
    supG = max(float(G(-0.1, t) / (3 + t) ** 1.2) for t in tt)
    supL = max(float(G(-0.1, t) / (1.9 ** 2 + t ** 2)) for t in tt)
    Cl = float(mp.mpf(3) ** 0.6 * (2 * mp.pi) ** -1.2 * zK) * supL
    Cr = float(zK) / 3.1 ** 2
    Cstar = max(Cl, Cr)
    print(f"   zeta_K(11/10) = {float(zK):.6f}; sup G_(-1/10)/(3+|t|)^(6/5) = {supG:.6f}; "
          f"sup G/(1.9^2+t^2) = {supL:.6f}")
    print(f"   C_left = {Cl:.6f}, C_right = {Cr:.6f}, C_* = {Cstar:.6f}; "
          f"=> |L| <= {Cstar*(3.1/3)**2:.4f} Q^(3/5) (3+|t|)^2")
    rep("explicit constant C_* < 1", Cstar < 1)

    # (c) actual L values for the test characters
    mp.mp.dps = 60
    chars = make_test_chars()
    OM = mp.exp(2j * mp.pi / 3)
    W6 = mp.exp(2j * mp.pi / 6)
    c2 = 2 * mp.pi / mp.sqrt(3)

    def root6(c):
        k = round(cmath.phase(c) / (math.pi / 3)) % 6
        return W6 ** k
    pts = lattice(300)
    worst = 0
    for ch in chars:
        Q = ch.Q
        fm = ch.f[0] + ch.f[1] * OM
        gm = mp.mpc(0)
        for x in ch.reps:
            val = ch.phi(x)
            if val != 0:
                zz = (x[0] + x[1] * OM) / fm
                gm += root6(val) * mp.exp(4j * mp.pi * mp.im(zz) / mp.sqrt(3))
        gm /= mp.sqrt(Q)
        eps0 = gm / (fm * mp.sqrt(Q))
        data = [(x, root6(ch.phi(x))) for x in pts if ch.phi(x) != 0]
        v0 = 1 / mp.sqrt(Q)

        def L(s, conjpsi=False):
            s = mp.mpc(s)
            A = mp.mpc(0)
            B = mp.mpc(0)
            for x, ph in data:
                n = norm(x)
                z = x[0] + x[1] * OM
                if conjpsi:
                    cA, cB = mp.conj(ph) * z, ph * mp.conj(z)
                    e = mp.conj(eps0)
                else:
                    cA, cB = ph * mp.conj(z), mp.conj(ph) * z
                    e = eps0
                A += cA * (c2 * n) ** (-(s + 0.5)) * mp.gammainc(s + 0.5, c2 * n * v0)
                B += cB * (c2 * n) ** (-(1.5 - s)) * mp.gammainc(1.5 - s, c2 * n / (Q * v0))
            I = A - e * mp.mpf(Q) ** (1.5 - s) * B
            return I / (6 * c2 ** (-(s + 0.5)) * mp.gamma(s + 0.5))
        wloc = 0
        for sig in [-0.1, 0.25, 0.5, 0.75, 1.1]:
            for t in [0, 2, -5, 10, -20, 30]:
                s = mp.mpc(sig, t)
                val = abs(L(s))
                ratio = float(val / (mp.mpf(Q) ** 0.6 * abs(s + 2) ** 2))
                wloc = max(wloc, ratio)
        worst = max(worst, wloc)
        print(f"   [{ch.name}] max |L(s)|/(Q^(3/5)|s+2|^2) on grid = {wloc:.4f}")
        # left-line identity (FE) and falsifier
        errs, ferrs = [], []
        m10 = -mp.mpf(1) / 10
        for t in [0.5, 7, -13]:
            s = mp.mpc(m10, t)
            lhs = abs(L(s))
            Gq = abs(mp.gamma(1.5 - s) / mp.gamma(0.5 + s))
            rhs = (3 * Q) ** (mp.mpf(3) / 5) * (2 * mp.pi) ** (-mp.mpf(6) / 5) * Gq * abs(L(1 - s, conjpsi=True))
            frhs = (3 * Q) ** (mp.mpf(1) / 2) * (2 * mp.pi) ** (-mp.mpf(6) / 5) * Gq * abs(L(1 - s, conjpsi=True))
            errs.append(float(abs(lhs - rhs) / rhs))
            ferrs.append(float(abs(lhs - frhs) / rhs))
        rep(f"   [{ch.name}] left-line FE identity", max(errs) < 1e-20,
            f"(err {max(errs):.1e}; falsifier (3Q)^(1/2): err {min(ferrs):.2e})")
        FALS.append((f"[{ch.name}] (3Q)^(1/2)", min(ferrs)))
    rep(f"actual values respect C_* (worst ratio {worst:.4f} <= C_* = {Cstar:.4f})", worst <= Cstar)
    fals_bad = [n for n, e in FALS if not e > 1e-6]          # falsifier verdict (added 2026-10-12)
    results.append(("falsifiers", not fals_bad))
    print(f"falsifiers: the (3Q)^(1/2) variant fails (rel err > 1e-6) for all {len(FALS)} characters:",
          "yes" if not fals_bad else f"NO {fals_bad}")
    print("ALL PASS" if all(ok for _, ok in results) else "SOME FAILED")
    raise SystemExit(0 if all(ok for _, ok in results) else 1)


if __name__ == "__main__":
    main()
