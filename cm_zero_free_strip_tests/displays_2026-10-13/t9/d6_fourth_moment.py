"""d6: size check of (F.3) eq:s5-plain-bound, case (1) (z = 0, Q = 1, no restriction on the lengths), on actual rows.

   sum_{k} |S^(e)_{psi_k}(n1;W) S^(e)_{psi_k}(n2;W)|^2 << Z^{M+eps},   S^(e)_psi(n;W) = Z^{-n/2} sum_l alpha(l)^e psi(l) W(N l/Z^n).
Rows: all ideals k coprime to 6 with N(k) <= K = 1500 (one primary generator each, k = 1 included), row character
psi_k(n) = (n/k)_6 = prod_{pi^v || k} (n/pi)_6^v on primary n coprime to 6 (residue symbols modulo the row, the
computable form used in tests_cm/phase2_2026-10-12/s5/t5; tau = 1, S = {lambda, 2}).  W = C_c^oo bump on [1,2].
Here Z^M ~ number of rows; the claim predicts that the mean (1/#rows) sum_k |S S|^2 stays O(Z^eps) for ALL lengths,
including X1 X2 up to 1e9 >> K (far outside the large-sieve range X1 X2 <= K), and for every row, including the
principal row k = 1 (Remark 1 of Appendix F: in the CM case case (1) holds for all rows).
This is a finite sanity check, not a proof: the constant is not computable.
Falsifier: e = 0 (finite-order case) on the same rows, where the principal row k = 1 has |S(X)|^2 ~ c X, so the mean
fourth moment grows like X1 X2 / #rows.
The fast residue symbols are cross-checked against s5lib.chi_p.
"""
import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t9lib import (report, falsifier, summary, cx, mul, N, OM, chi_p, primes_upto, is_primary, divides, exact_div,
                   residue_reps, canon_index, LAM, TWO)

K = 1500
P = [p for p in primes_upto(K) if N(p) <= K]
SIX = np.exp(1j * np.pi * np.arange(6) / 3)

def symbol_table(p):
    """Function (A, B) -> (A + B omega / pi)_6 for numpy arrays, pi = p a primary prime of norm <= K."""
    q = N(p)
    if p[1] != 0:                                    # split prime: O/pi = F_q, omega -> w0
        c, d = p
        w0 = (-c * pow(d, -1, q)) % q
        assert (w0 * w0 + w0 + 1) % q == 0
        zeta = (1 + w0) % q                          # 1 + omega = e^{i pi/3}
        val = {}
        z = 1
        for j in range(6):
            val[z] = SIX[j]; z = z * zeta % q
        e6 = (q - 1) // 6
        tab = np.zeros(q, dtype=complex)
        for r in range(1, q):
            tab[r] = val[pow(r, e6, q)]
        return lambda A, B: tab[(A + B * w0) % q]
    else:                                            # inert prime -q: table on O/q via s5lib
        reps = residue_reps(p)
        tab = np.array([chi_p(r, p) for r in reps], dtype=complex)
        return lambda A, B: tab[canon_index(A, B, p)]
SYM = {p: symbol_table(p) for p in P}

# cross-check against s5lib.chi_p
rng = np.random.default_rng(11)
bad = 0
for p in P[:40] + P[-10:]:
    for _ in range(30):
        a = (int(rng.integers(-200, 200)), int(rng.integers(-200, 200)))
        if abs(SYM[p](np.array([a[0]]), np.array([a[1]]))[0] - complex(chi_p(a, p))) > 1e-12: bad += 1
report("F3-symbols", bad, 0, "fast sextic symbols (F_p discrete powers) agree with s5lib.chi_p on 1500 random elements")

# rows: primary k coprime to 6 with N(k) <= K, factored into the primary primes of P
def factor(k):
    out = {}; rem = k
    for p in P:
        if N(rem) == 1: break
        while N(rem) % N(p) == 0 and divides(p, rem):
            out[p] = out.get(p, 0) + 1; rem = exact_div(rem, p)
    assert N(rem) == 1
    return out
rows = []
R = int(math.ceil(2 * math.sqrt(K / 3))) + 2
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        k = (a, b)
        if 0 < N(k) <= K and is_primary(k) and not divides(TWO, k):
            rows.append((k, factor(k)))
print(f"{len(rows)} rows (ideals coprime to 6, N(k) <= {K}, k = 1 included)")

def bump(x):
    out = np.zeros_like(x); m = (x > 1) & (x < 2); u = 2 * x[m] - 3
    out[m] = np.exp(-1.0 / (1 - u * u)); return out
Xs = [10.0, 30.0, 100.0, 300.0, 1e3, 3e3, 1e4, 3e4]
Nmax = int(2 * Xs[-1]) + 1
RR = int(math.ceil(2 * math.sqrt(Nmax / 3))) + 2
a = np.arange(-RR, RR + 1); A, B = np.meshgrid(a, a, indexing='ij'); A = A.ravel(); B = B.ravel()
n = A * A - A * B + B * B
keep = (n > 0) & (n <= Nmax) & ((A - 1) % 3 == 0) & (B % 3 == 0) & ~((A % 2 == 0) & (B % 2 == 0))
A, B, n = A[keep], B[keep], n[keep]
Z = A + B * OM
ang = Z / np.abs(Z)
Wx = {X: bump(n / X) for X in Xs}
sel = {X: Wx[X] != 0 for X in Xs}

t0 = time.time()
S = {e: np.zeros((len(rows), len(Xs)), dtype=complex) for e in (1, -1, 0)}
for i, (k, fac) in enumerate(rows):
    psi = np.ones(len(A), dtype=complex)
    for p, v in fac.items():
        psi *= SYM[p](A, B) ** v
    for e in (1, -1, 0):
        c = psi * ang ** e
        for j, X in enumerate(Xs):
            m = sel[X]
            S[e][i, j] = np.sum(c[m] * Wx[X][m]) / math.sqrt(X)
print(f"   plain sums computed in {time.time() - t0:.1f} s")

def mean4(e):
    M = np.zeros((len(Xs), len(Xs)))
    for j1 in range(len(Xs)):
        for j2 in range(len(Xs)):
            M[j1, j2] = np.mean(np.abs(S[e][:, j1] * S[e][:, j2]) ** 2)
    return M
def growth(diag):
    """Least-squares exponent of X1 X2 = X^2 in the diagonal means, and the ratio of the largest mean at the upper
    half of the scales to the largest at the lower half.  (Codex review, 2026-10-13: the earlier criterion compared
    the maximum of M with the maxima of its diagonal blocks, which Cauchy-Schwarz makes automatic.)"""
    slope = np.polyfit(np.log(np.array(Xs) ** 2), np.log(diag), 1)[0]
    h = len(Xs) // 2
    return slope, diag[h:].max() / diag[:h].max()
ok = True; worst = (-np.inf, 0.0)
for e in (1, -1):
    M = mean4(e)
    diag = np.diag(M)
    print(f"   e={e:+d}: mean_k |S(X)S(X)|^2 for X = " + ", ".join(f"{X:g}" for X in Xs) + ": " + ", ".join(f"{v:.3g}" for v in diag))
    slope, ratio = growth(diag)
    print(f"          fitted exponent of X1 X2: {slope:.3f}; largest mean at X >= {Xs[len(Xs)//2]:g} / largest at smaller X: {ratio:.3g}")
    ok &= slope <= 0.05 and ratio <= 3
    worst = (max(worst[0], slope), max(worst[1], ratio))
report("F3-sanity", 0.0 if ok else 1.0, 0.0,
       f"mean fourth moment over {len(rows)} rows, X1 X2 from 1e2 to 9e8, e = +-1: fitted exponent of X1 X2 at most "
       f"{worst[0]:.3f} (required <= 0.05), upper-half/lower-half maximum at most {worst[1]:.2f} (required <= 3)")
M0 = mean4(0)
print("   e= 0: mean_k |S(X)S(X)|^2: " + ", ".join(f"{v:.3g}" for v in np.diag(M0)))
i1 = [i for i, (k, fac) in enumerate(rows) if k == (1, 0)][0]
print(f"          principal row k = 1, e = 0: |S(X)|^2 = " + ", ".join(f"{abs(S[0][i1, j]) ** 2:.3g}" for j in range(len(Xs))))
print(f"          principal row k = 1, e = +1: |S(X)|^2 = " + ", ".join(f"{abs(S[1][i1, j]) ** 2:.3g}" for j in range(len(Xs))))
rest = [i for i in range(len(rows)) if i != i1]
print("   e= 0 without the principal row: " + ", ".join(f"{np.mean(np.abs(S[0][rest, j]) ** 4):.3g}" for j in range(len(Xs))))
slope0, ratio0 = growth(np.diag(M0))
print(f"   e= 0: fitted exponent of X1 X2: {slope0:.3f}; upper/lower ratio {ratio0:.3g}")
falsifier("F3-sanity", "e = 0 (no angular factor): principal-row main term, growth ~ X1 X2 (fitted exponent > 0.05)",
          max(slope0 - 0.05, 0.0) + (1.0 if ratio0 > 3 else 0.0), 1e-9)
summary()
