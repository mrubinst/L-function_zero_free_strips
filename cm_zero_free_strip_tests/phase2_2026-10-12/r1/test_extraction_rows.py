"""r1 tests of Section 5 (guide lines 621-778): the rows u_0 pi^6 of the extraction.

E1  For a prime pi | d (pi coprime to 6) the row u_0 pi^6 equals the row u_0 exactly:
    chi_n(u_0 pi^6) = chi_n(u_0) for every column n, so A_{u_0 pi^6}(D) = A_{u_0}(D).
    (Supports the optional simplification R1-6 of the report: the condition pi \\nmid d in P_Y is not needed.)
    Falsifier: for pi \\nmid d the rows differ.
E2  Lemma 2(iii) with its explicit constant: |A_{u_0}(D) - A_{u_0 pi^6}(D)| <= (1+16b)||W||_oo D/Y for every
    pi in P_Y, and the count #{n : pi | n, N(n) <= bD} <= 1 + 16 bD/Y used in its proof.
    Sensitivity check: the count is NOT <= bD/(10 Y) (so the test is not vacuous).
E3  The rows u_0 pi^6, pi in P_Y, are pairwise distinct and have norm N_0 N(pi)^6 <= N_0 Y^6 (Lemma 2(ii)).

Run with  python3 -I test_extraction_rows.py
"""
import os
import sys
import math
import cmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *  # noqa

rep = Reporter("test_extraction_rows")


def W(y):
    # smooth, complex, supported in [1, 2]  (a = 1, b = 2)
    if not (1.0 < y < 2.0):
        return 0j
    return math.exp(-1.0 / ((y - 1.0) * (2.0 - y))) * (1.0 + 0.3j * y)


b = 2.0
Winf = max(abs(W(1 + k / 10000.0)) for k in range(1, 10000))

Dmax = 6000
COLS = primaries(int(b * Dmax), squarefree=True)


def row_coeffs(u, D):
    """dict n -> mu(n) alpha(n) chi_n(u) W(N n/D) over the columns."""
    out = {}
    for n in COLS:
        w = W(N(n) / D)
        if w == 0:
            continue
        out[n] = mu(n) * alpha(n) * chi(n, u) * w
    return out


def A(u, D):
    return sum(row_coeffs(u, D).values())


# ---------------------------------------------------------------- E1
d = 7 * 13 * 19
u0 = ((4 * d) ** 5, 0)
primes_d = [p for p in primary_primes(400) if divides(p, (d, 0))]
D = 3000
worst = 0.0
nz_cols = 0
for p in primes_d:
    up = mul(u0, pw(p, 6))
    for n in COLS:
        if N(n) > 2 * D:
            break
        c1 = chi(n, u0)
        c2 = chi(n, up)
        worst = max(worst, abs(c1 - c2))
        if abs(c1) > 0:
            nz_cols += 1
rep.check("E1", worst == 0 and len(primes_d) == 6,
          "d = 7*13*19: for each of the %d primes pi | d, chi_n(u_d pi^6) = chi_n(u_d) for all %d columns with N(n) <= %d "
          "(max difference %.1e), so A_{u_d pi^6}(D) = A_{u_d}(D) exactly" % (len(primes_d), sum(1 for n in COLS if N(n) <= 2 * D), 2 * D, worst))
# falsifier: a prime not dividing d of comparable norm
p_out = next(p for p in primary_primes(400) if N(p) > 13 and not divides(p, (d, 0)))
diff = abs(A(u0, D) - A(mul(u0, pw(p_out, 6)), D))
rep.falsifier("E1", "pi = %s does not divide d" % (p_out,), diff > 1e-6, "|A_u0 - A_u0pi^6| = %.3f" % diff)

# ---------------------------------------------------------------- E2 and E3
worst_ratio = 0.0
worst_count_ratio = 0.0
min_count_ratio_fals = 1e9
ncases = 0
distinct_ok = True
for d in (1, 2, -3, 7):
    u0 = ((4 * d) ** 5, 0)
    N0 = N(u0)
    for D in (1500, Dmax):
        base = row_coeffs(u0, D)
        Abase = sum(base.values())
        for Y in (30, 120, 600):
            if Y > D:
                continue
            PY = [p for p in primary_primes(Y) if Y / 2 < N(p) <= Y and not divides(p, (d, 0))][:12]
            rows = [mul(u0, pw(p, 6)) for p in PY]
            if len(set(rows)) != len(rows) or any(N(r) != N0 * N(p) ** 6 or N(p) ** 6 > Y ** 6 for r, p in zip(rows, PY)):
                distinct_ok = False
            for p in PY:
                # the difference is the sum over columns divisible by p (chi_n(pi)^6 = 1_{pi not | n})
                diffsum = sum(v for n, v in base.items() if divides(p, n))
                direct = Abase - A(mul(u0, pw(p, 6)), D) if ncases < 6 else diffsum
                assert abs(direct - diffsum) < 1e-9, (d, D, Y, p, direct, diffsum)
                bound = (1 + 16 * b) * Winf * D / Y
                worst_ratio = max(worst_ratio, abs(diffsum) / bound)
                cnt = sum(1 for n in primaries(int(b * D)) if divides(p, n))  # all primary n coprime to 6, pi | n
                worst_count_ratio = max(worst_count_ratio, cnt / (1 + 16 * b * D / Y))
                min_count_ratio_fals = min(min_count_ratio_fals, cnt / (b * D / (10 * Y)))
                ncases += 1
rep.check("E2", worst_ratio <= 1 and worst_count_ratio <= 1,
          "%d cases (d in {1,2,-3,7}, D in {1500,6000}, Y in {30,120,600}, up to 12 primes each): "
          "max |A_u0 - A_u0pi^6| / ((1+16b)||W|| D/Y) = %.3f; max count / (1+16bD/Y) = %.3f; "
          "the difference equals the sum over pi | n (checked directly in 6 cases)" % (ncases, worst_ratio, worst_count_ratio))
rep.falsifier("E2", "count <= bD/(10Y) (sensitivity)", min_count_ratio_fals > 1,
              "min count / (bD/(10Y)) = %.2f" % min_count_ratio_fals)
rep.check("E3", distinct_ok, "rows u_0 pi^6 pairwise distinct with N = N_0 N(pi)^6 <= N_0 Y^6 in every case")

sys.exit(1 if rep.summary() else 0)
