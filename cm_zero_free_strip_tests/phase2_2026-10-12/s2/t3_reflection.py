"""T3. Reflection of a CM plain sum (insert, Lemma s2-reflection).

Claim: for psi in F_CM with conductor norm Q, root number W(psi), C = 3Q/(2pi)^2, X > 0, Y = C/X,
   X^{-1/2} sum_a psi(a) W(Na/X) = W(psi) Y^{-1/2} sum_a conj(psi)(a) W#(Na/Y),
   M W#(u) = M W(1-u) Gamma(u+1/2)/Gamma(3/2-u),   M W(s) = int_0^oo W(y) y^s dy/y.
Test weight: W(y) = exp(-(log y)^2), M W(s) = sqrt(pi) exp(s^2/4) (entire, Gaussian decay).
Also: W#(y) = M W(3/2) y^{1/2} + O(y^{3/2}) as y -> 0 (no constant term, unlike OAI78's Gamma(s)/Gamma(1-s)).
Falsifiers: OAI78's finite-order kernel Gamma(u)/Gamma(1-u); W(psi) replaced by conj W(psi).
"""
import math
import cmath
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from cmhecke import *

mp.mp.dps = 20
results = []
FALS = []   # (name, error) of wrong variants that must fail with error > 1e-6 (added 2026-10-12)
_FORCE = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]   # self-test (added 2026-10-12): first check counts as failed


def rep(name, ok, detail=""):
    if _FORCE[0]:
        _FORCE[0] = False; ok = False; name = "[CM_TESTS_FORCE_FAIL=1] " + name
    results.append((name, ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name} {detail}")


def MW(s):
    return mp.sqrt(mp.pi) * mp.exp(s * s / 4)


def Wsharp(y, kernel="cm", c=0.5):
    y = mp.mpf(y)

    def f(tau):
        u = mp.mpc(c, tau)
        if kernel == "cm":
            k = mp.gamma(u + 0.5) / mp.gamma(1.5 - u)
        else:
            k = mp.gamma(u) / mp.gamma(1 - u)
        return MW(1 - u) * k * mp.exp(-u * mp.log(y))
    v = mp.quad(f, [-60, -20, -5, 0, 5, 20, 60]) / (2 * mp.pi)
    return mp.re(v) if False else v


def main():
    chars = make_test_chars()
    for ch in [chars[1], chars[2], chars[3]]:
        Q = ch.Q
        Wpsi = ch.root_number()
        C = 3 * Q / (2 * math.pi) ** 2
        Y = 0.45   # choose X so that the dual length is small: W# decays only like exp(-(log y)^2 + O(log y loglog y))
        X = C / Y
        # LHS
        pts = lattice(int(X * math.exp(6.5)) + 1)
        lhs = 0
        for x in pts:
            v = ch.psi_elem(x)
            if v != 0:
                lhs += v * math.exp(-math.log(norm(x) / X) ** 2)
        lhs = lhs / 6 / math.sqrt(X)
        # RHS: group by norm
        nmax = int(Y * math.exp(9.5)) + 1
        pts2 = lattice(nmax)
        bynorm = {}
        for x in pts2:
            v = ch.psi_elem(x)
            if v != 0:
                bynorm[norm(x)] = bynorm.get(norm(x), 0) + v.conjugate()
        rhs = 0
        rhs_f = 0
        rhs_cut8 = 0
        for n, coef in sorted(bynorm.items()):
            if abs(coef) < 1e-12:
                continue
            ws = complex(Wsharp(n / Y))
            rhs += coef * ws
            if n <= Y * math.exp(8.0):
                rhs_cut8 += coef * ws
                rhs_f += coef * complex(Wsharp(n / Y, kernel="fin"))
        print(f"      dual tail beyond Y e^8 (relative): {abs(rhs-rhs_cut8)/abs(rhs):.2e}")
        rhs = Wpsi * rhs / 6 / math.sqrt(Y)
        rhs_f = Wpsi * rhs_f / 6 / math.sqrt(Y)
        rhs_cw = Wpsi.conjugate() * (rhs / Wpsi)
        err = abs(lhs - rhs) / abs(lhs)
        rep(f"[{ch.name}] reflection identity, X={X:.4f}, Y={Y:.4f}", err < 1e-11,
            f"lhs={lhs:.12f} rhs={rhs:.12f} rel err {err:.2e}; falsifiers: finite-order kernel "
            f"{abs(lhs-rhs_f)/abs(lhs):.2e}, conj W {abs(lhs-rhs_cw)/abs(lhs):.2e}")
        FALS.append((f"[{ch.name}] finite-order kernel", abs(lhs-rhs_f)/abs(lhs)))
        FALS.append((f"[{ch.name}] conj W", abs(lhs-rhs_cw)/abs(lhs)))
    # behaviour at zero
    lim = MW(1.5)
    for y in [1e-2, 1e-3, 1e-4]:
        w = Wsharp(y, c=0.0)
        r = abs(w / mp.sqrt(y) - lim) / abs(lim)
        wf = Wsharp(y, kernel="fin", c=0.5)
        print(f"   y={y:g}: W#(y)/y^(1/2) = {complex(w/mp.sqrt(y)):.8f}, M W(3/2) = {float(lim):.8f}, rel diff {float(r):.2e};"
              f" finite-order W#(y) = {complex(wf):.6f} (-> M W(1) = {float(MW(1)):.6f})")
    r = abs(Wsharp(1e-4, c=0.0) / mp.sqrt(1e-4) - lim) / abs(lim)
    rep("W#(y) = M W(3/2) y^(1/2) + O(y^(3/2)) as y->0", r < 1e-3, f"(rel diff at 1e-4: {float(r):.1e})")
    fals_bad = [n for n, e in FALS if not e > 1e-6]          # falsifier verdict (added 2026-10-12)
    results.append(("falsifiers", not fals_bad))
    print(f"falsifiers: finite-order kernel and conj W fail (rel err > 1e-6) for all {len(FALS)//2} characters:",
          "yes" if not fals_bad else f"NO {fals_bad}")
    print("ALL PASS" if all(ok for _, ok in results) else "SOME FAILED")
    raise SystemExit(0 if all(ok for _, ok in results) else 1)


if __name__ == "__main__":
    main()
