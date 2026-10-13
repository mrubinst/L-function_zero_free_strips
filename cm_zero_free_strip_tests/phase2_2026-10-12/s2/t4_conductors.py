"""T4. Lemma s2-pres and Lemma s2-rowcond, by brute force on (O/m)^x.

 (a) alphabar has phi(x) = (unit = x mod 3) and conductor exactly (3) (not (lambda), not (1)).
 (b) For finite-order psi_0 with finite part phi_0 (a character of (O/m0)^x trivial on units),
     psi* = alphabar psi_0 has finite part phi = phi_alphabar * phi_0; its conductor divides lcm(3, f0),
     and v_p(f(psi*)) = v_p(f0) at every prime p != lambda.  (Exponent at lambda may drop or rise.)
 (c) phi(zeta) = zeta for the CM character, phi_0(zeta) = 1 for the finite-order one.
Characters used: phi_0 = zeta(x)^k chi_pi(x)^j with k + j (N pi - 1)/6 = 0 mod 6, pi primary of norm 7, 13, 19,
and k, j ranging over all residues mod 6.
Falsifier: the analogous statement with lcm((lambda), f0) in place of lcm((3), f0) must fail for some phi_0.
"""
import os
import sys
import itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmhecke import *

LAM = (1, 2)


def val(m, p):
    """valuation of element m at prime element p"""
    v = 0
    while divisible(m, p):
        # m = m / p exactly
        mc = cval(m) / cval(p)
        a = round(mc.real + mc.imag / math.sqrt(3))
        b = round(2 * mc.imag / math.sqrt(3))
        assert mul((a, b), p) == m
        m = (a, b)
        v += 1
    return v


def divide(m, p):
    mc = cval(m) / cval(p)
    a = round(mc.real + mc.imag / math.sqrt(3))
    b = round(2 * mc.imag / math.sqrt(3))
    assert mul((a, b), p) == m
    return (a, b)


def conductor(phi, m, primes):
    """smallest divisor g of m such that phi (a char of (O/m)^x) is trivial on {x = 1 mod g, (x,m)=1}."""
    reps = [x for x in residue_reps(m) if phi(x) != 0]
    g = m
    changed = True
    while changed:
        changed = False
        for p in primes:
            if divisible(g, p):
                h = divide(g, p)
                if all(abs(phi(x) - 1) < 1e-9 for x in reps if divisible(sub(x, (1, 0)), h)):
                    g = h
                    changed = True
                    break
    return g


def main():
    ok = True
    phiA = lambda x: (0 if unit_mod3(x) is None else cval(UNITS[unit_mod3(x)]))
    gA = conductor(phiA, (3, 0), [LAM])
    print("conductor of alphabar: norm", norm(gA), "(expect 9)")
    ok &= norm(gA) == 9
    count = 0
    lam_change = {}
    fals_fail = 0
    for p in [7, 13, 19]:
        pi = primary_primes_of_norm(p)[0]
        chi = SexticSymbol(pi)
        e = (p - 1) // 6
        m = mul((3, 0), pi)
        for k, j in itertools.product(range(6), range(6)):
            if (k + j * e) % 6 != 0:
                continue

            def phi0(x, k=k, j=j, chi=chi):
                i = unit_mod3(x)
                c = chi(x)
                if i is None or c == 0:
                    return 0
                return cval(UNITS[i]) ** k * c ** j

            def phi(x, k=k, j=j, chi=chi):
                i = unit_mod3(x)
                c = chi(x)
                if i is None or c == 0:
                    return 0
                return cval(UNITS[i]) ** (k + 1) * c ** j
            for u in UNITS:
                ok &= abs(phi0(u) - 1) < 1e-9 and abs(phi(u) - cval(u)) < 1e-9
            f0 = conductor(phi0, m, [LAM, pi])
            f1 = conductor(phi, m, [LAM, pi])
            # f1 | lcm(3, f0): since 3 = -lambda^2, lcm exponent at lambda is max(2, v_lam(f0))
            cond_div = val(f1, LAM) <= max(2, val(f0, LAM)) and val(f1, pi) <= val(f0, pi)
            same_away = val(f1, pi) == val(f0, pi)
            ok &= cond_div and same_away
            # falsifier: lcm((lambda), f0)
            if not (val(f1, LAM) <= max(1, val(f0, LAM))):
                fals_fail += 1
            lam_change[(val(f0, LAM), val(f1, LAM))] = lam_change.get((val(f0, LAM), val(f1, LAM)), 0) + 1
            count += 1
    print(f"{count} pairs (psi_0, alphabar psi_0) checked; (v_lambda(f0), v_lambda(f*)) counts: {lam_change}")
    print(f"falsifier lcm((lambda),f0) fails for {fals_fail} of them (must be > 0)")
    ok &= fals_fail > 0
    if os.environ.get("CM_TESTS_FORCE_FAIL") == "1":   # self-test of the exit status (added 2026-10-12)
        ok = False; print("CM_TESTS_FORCE_FAIL=1: the checks are counted as failed")
    print("ALL PASS" if ok else "SOME FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
