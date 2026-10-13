#!/usr/bin/env python3
"""Numerical check of the deferred identity in CM/Hecke.lean (`invCoeff_eq_heckeSum`).

Claim (cm-eleven-twelfths-guide.tex, eqs. (ud) and (normgroup)), unnormalized form: for n >= 1,

    invCoeff d n = sum_{c primary, N(c) = n, (c, 6d) = 1}  mu_K(c) * sqrt(N c) * psi_d(c),

    psi_d(c) = alpha(c) * chi_c(u_d),   u_d = (4d)^5,   alpha(c) = c / |c|,
    chi_c(u) = (u / c)_6   (multiplicative in c; at a prime pi it is the sixth root of
                            unity congruent to u^((N pi - 1)/6) mod pi, or 0 if pi | u),

with omega = e^{2 pi i / 3}, c = A + B omega primary iff A = 1, B = 0 (mod 3), and
`invCoeff d` the coefficients of prod_{p not | 6d} (1 - a_p p^{-s} + p^{1-2s}) (CM/Basic.lean):
multiplicative, invCoeff(p) = -a_p, invCoeff(p^2) = p, invCoeff(p^k) = 0 (k >= 3), and 0 on n
with gcd(n, 6d) > 1. Here a_p = p - #{(x, y) in F_p^2 : y^2 = x^3 + d}.

Note sqrt(N c) * psi_d(c) = c * chi_c((4d)^5) = c * conj(chi_c(4d)).

a_p comes from PARI/GP `ellap` and is cross-checked against a direct point count for p < 300.

Falsifiers (each must FAIL for some d):
  F1  alpha(c) replaced by conj(alpha(c));
  F2  u_d = (4d)^5 replaced by 4d (no conjugation of the sextic symbol);
  F3  u_d = (4d)^5 replaced by d^5 (the factor 4 dropped);
  F4  "primary" taken as c = 2 (mod 3) (Ireland-Rosen's convention), i.e. c -> -c.

Exit status 0 iff the claim holds for every tested (d, n) and every falsifier fails somewhere.
"""

import cmath
import math
import subprocess
import sys

DS = [1, 2, 3, -1, -2, -3, 5, 7, -7, 11, 12, 16, -16, 30, -432, 4, 9, 25]
NMAX = 600
OMEGA = cmath.exp(2j * math.pi / 3)


# ---------- Z[omega] arithmetic: (a, b) <-> a + b*omega, omega^2 = -1 - omega ----------

def emul(x, y):
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def econj(x):
    a, b = x
    return (a - b, -b)


def enorm(x):
    a, b = x
    return a * a - a * b + b * b


def ediv_exact(x, y):
    """x / y if y divides x in Z[omega], else None."""
    n = enorm(y)
    a, b = emul(x, econj(y))
    if a % n or b % n:
        return None
    return (a // n, b // n)


def to_c(x):
    return x[0] + x[1] * OMEGA


def is_primary(x):
    return x[0] % 3 == 1 and x[1] % 3 == 0


# ---------- rational helpers ----------

def factor(n):
    f = {}
    m = n
    p = 2
    while p * p <= m:
        while m % p == 0:
            f[p] = f.get(p, 0) + 1
            m //= p
        p += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def is_prime(n):
    return n >= 2 and factor(n) == {n: 1}


def primes_upto(n):
    return [p for p in range(2, n + 1) if is_prime(p)]


# ---------- primary primes and the sextic residue symbol ----------

_split_cache = {}


def primary_primes_over(p):
    """Primary primes of Z[omega] above the rational prime p >= 5."""
    if p % 3 == 2:
        return [(-p, 0)]  # inert; -p = 1 (mod 3)
    if p not in _split_cache:
        r = int(math.isqrt(4 * p // 3)) + 2
        sols = [(A, B) for A in range(-r, r + 1) for B in range(-r, r + 1)
                if enorm((A, B)) == p and is_primary((A, B))]
        assert len(sols) == 2 and sols[1] == econj(sols[0]), (p, sols)
        _split_cache[p] = sols
    return _split_cache[p]


def sextic_symbol_prime(u, pi):
    """(u / pi)_6 in C for a rational integer u and a primary prime pi not dividing 6."""
    N = enorm(pi)
    if pi[1] == 0:  # inert pi = -p, residue field F_{p^2} = F_p[omega]
        p = -pi[0]
        if u % p == 0:
            return 0
        # compute u^((p^2-1)/6) in F_p[omega], compare with the six roots of unity (1+omega)^k
        e = (N - 1) // 6
        res, base = (1, 0), (u % p, 0)
        while e:
            if e & 1:
                res = tuple(t % p for t in emul(res, base))
            base = tuple(t % p for t in emul(base, base))
            e >>= 1
        z, zk = (1, 1), (1, 0)
        for k in range(6):
            if (zk[0] - res[0]) % p == 0 and (zk[1] - res[1]) % p == 0:
                return cmath.exp(1j * math.pi * k / 3)
            zk = emul(zk, z)
        raise AssertionError("not a sixth root of unity")
    p = N
    if u % p == 0:
        return 0
    A, B = pi
    r = (-A * pow(B, -1, p)) % p  # omega = r (mod pi)
    assert (r * r + r + 1) % p == 0
    t = pow(u % p, (p - 1) // 6, p)
    zeta = (1 + r) % p  # 1 + omega = -omega^2 = e^{i pi/3}
    for k in range(6):
        if pow(zeta, k, p) == t:
            return cmath.exp(1j * math.pi * k / 3)
    raise AssertionError("not a sixth root of unity")


def primary_factorization(c):
    """[(pi, e)] with c = prod pi^e, pi primary primes (c primary, coprime to 3)."""
    out = []
    rest = c
    for p in factor(enorm(c)):
        assert p != 3
        for pi in primary_primes_over(p):
            e = 0
            while True:
                q = ediv_exact(rest, pi)
                if q is None:
                    break
                rest, e = q, e + 1
            if e:
                out.append((pi, e))
    assert rest == (1, 0), (c, rest)  # primary => unit part is exactly 1
    return out


def mu_K(fact):
    if any(e >= 2 for _, e in fact):
        return 0
    return (-1) ** len(fact)


def sextic_symbol(u, fact):
    v = 1
    for pi, e in fact:
        v *= sextic_symbol_prime(u, pi) ** e
    return v


# ---------- a_p: PARI/GP and a direct count ----------

def ellap_gp(d, plist):
    prog = "E=ellinit([0,0,0,0,%d]); print(vector(%d,i,ellap(E,%s[i])))\n" % (
        d, len(plist), str(plist).replace(" ", ""))
    out = subprocess.run(["gp", "-q", "-f"], input=prog, capture_output=True, text=True,
                         check=True).stdout.strip()
    vals = [int(t) for t in out.strip("[]").split(",")]
    return dict(zip(plist, vals))


def ap_direct(d, p):
    sq = [0] * p
    for y in range(p):
        sq[y * y % p] += 1
    return p - sum(sq[(x * x * x + d) % p] for x in range(p))


def inv_coeff(d, n, ap):
    if math.gcd(n, 6 * abs(d)) != 1:
        return 0
    v = 1
    for p, k in factor(n).items():
        v *= {1: -ap[p], 2: p}.get(k, 0)
    return v


# ---------- the Z[omega] side ----------

def primary_of_norm(n):
    r = int(math.isqrt(4 * n // 3)) + 1
    return [(A, B) for A in range(-r, r + 1) for B in range(-r, r + 1)
            if enorm((A, B)) == n and is_primary((A, B))]


def hecke_sum(d, n, variant=None):
    """sum over primary c, N c = n, (c, 6d) = 1, of mu_K(c) * sqrt(N c) * psi_d(c)."""
    if math.gcd(n, 6 * abs(d)) != 1:  # (c, 6d) = 1 in Z[omega]  <=>  gcd(N c, 6d) = 1
        return 0
    u = (4 * d) ** 5
    if variant == "F2":
        u = 4 * d
    if variant == "F3":
        u = d ** 5
    s = 0
    for c in primary_of_norm(n):
        fact = primary_factorization(c)
        m = mu_K(fact)
        if m == 0:
            continue
        alpha = to_c(c) / math.sqrt(n)
        if variant == "F1":
            alpha = alpha.conjugate()
        sign = -1 if variant == "F4" else 1  # -c is Ireland-Rosen primary; alpha(-c) = -alpha(c)
        s += m * sign * math.sqrt(n) * alpha * sextic_symbol(u, fact)
    return s


def main():
    ok = True
    plist = primes_upto(NMAX)
    # sanity: the paper's examples mu_K(-5) = -1, mu_K(7) = +1
    assert mu_K(primary_factorization((-5, 0))) == -1
    assert mu_K(primary_factorization((7, 0))) == 1
    print("Claim: invCoeff d n = sum_{c primary, N c = n, (c,6d)=1} mu_K(c) sqrt(N c) psi_d(c),")
    print("       psi_d(c) = alpha(c) (u_d / c)_6, u_d = (4d)^5;  n = 1..%d" % NMAX)
    fal = {k: 0 for k in ["F1", "F2", "F3", "F4"]}
    fal_d = {k: [] for k in fal}
    total = 0
    for d in DS:
        ap = ellap_gp(d, plist)
        for p in plist:
            if p < 300 and (6 * d) % p != 0:
                assert ap[p] == ap_direct(d, p), ("PARI vs direct count", d, p)
        maxerr, bad = 0.0, []
        nonzero = 0
        for n in range(1, NMAX + 1):
            lhs = inv_coeff(d, n, ap)
            rhs = hecke_sum(d, n)
            err = abs(lhs - rhs)
            maxerr = max(maxerr, err)
            nonzero += lhs != 0
            if err > 1e-6:
                bad.append((n, lhs, rhs))
            total += 1
        for k in fal:
            nb = sum(1 for n in range(1, NMAX + 1)
                     if abs(inv_coeff(d, n, ap) - hecke_sum(d, n, k)) > 1e-6)
            fal[k] += nb
            if nb:
                fal_d[k].append(d)
        status = "OK " if not bad else "FAIL"
        ok &= not bad
        print("  d = %5d: %s  max|lhs-rhs| = %.2e  (%d nonzero coefficients)%s" % (
            d, status, maxerr, nonzero, "" if not bad else "  first bad: %s" % (bad[:3],)))
    print("Identity: %s on %d pairs (d, n)." % ("TRUE" if ok else "FALSE", total))
    print("Falsifiers (must fail):")
    names = {"F1": "conj(alpha) instead of alpha",
             "F2": "u = 4d instead of (4d)^5",
             "F3": "u = d^5 instead of (4d)^5",
             "F4": "primary = 2 mod 3 (c -> -c)"}
    for k in fal:
        fails = fal[k] > 0
        ok &= fails
        print("  %s %-32s: %5d mismatches, fails for d in %s  -> %s" % (
            k, names[k], fal[k], fal_d[k], "fails as required" if fails else "DID NOT FAIL"))
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
