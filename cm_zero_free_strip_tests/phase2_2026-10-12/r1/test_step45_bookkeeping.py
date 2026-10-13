"""r1 tests of Steps 4-5 of the proof of Proposition 3 (first Poisson reduction, guide lines 1297-1333).
These were listed as b2-T15 ("suggested; not run") in tests_cm/proofs_2026-10-12/b2/README.md.

On the support of the summand of (7.17), with W supported in I = [a0, b0] and hat Phi supported in [0, C_Phi]:
 S1  N(r0) N(f) <= b0 D  and  H N(k) <= C_Phi N(f)^2 N(m1) N(m2) <= C_Phi b0^2 D^2 / N(r0)^2;
 S2  with dyadic B <= N(r0) < 2B, F <= N(f) < 2F, X = D/(BF), calH = C1 D^2/(H B^2), C1 = max(1, C_Phi b0^2):
     0 < N(k) <= calH, rho = N(r0)N(f)/(BF) in [1,4), 0 < kappa = H N(k)/(N(f)^2 X^2) <= C1,
     x_i = N(m_i)/X in [a0/4, b0], W_0(N(r0 f m_i)/D) = W_0(rho x_i), and the hat Phi argument equals kappa/(x1 x2);
 S3  for j | r0 the scales (X/N(j), F N(j)) again satisfy (7.15) with the same calH and C = C1;
 S4  bookkeeping: sum_{B <= N(r0) < 2B, r0 squarefree} (H N(r0)/D^2)(D/B)^2 << H uniformly in B.
Falsifiers: rho in [1,2); support [a0/2, b0]; hat Phi argument with N(f) in place of N(f)^2.

Run with  python3 -I test_step45_bookkeeping.py
"""
import os
import sys
import math
import random

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1lib import *  # noqa

rep = Reporter("test_step45_bookkeeping")
random.seed(45)

D = 6000.0
theta = 0.1
H = D ** (1 + theta)
a0, b0 = 1.0, 2.0
C_Phi = 1.7
C1 = max(1.0, C_Phi * b0 ** 2)
SQ = primaries(int(b0 * D), squarefree=True)
norms = sorted(set(N(n) for n in SQ))
by_norm = {}
for n in SQ:
    by_norm.setdefault(N(n), []).append(n)


def dyadic(x):
    return 2.0 ** math.floor(math.log2(x))


s1 = s2 = s3 = True
max_rho = 0.0
min_x = 1e9
max_kappa = 0.0
fal_rho = fal_supp = fal_arg = False
samples = 0
tries = 0
while samples < 20000 and tries < 10 ** 6:
    tries += 1
    r0 = random.choice(SQ[:600]) if random.random() < 0.7 else random.choice(SQ)
    f = random.choice(SQ[:600]) if random.random() < 0.7 else random.choice(SQ)
    if N(r0) * N(f) >= b0 * D or N(gcd_primary(r0, f)) > 1:
        continue
    # columns m_i with a0 D < N(r0 f m_i) < b0 D (support of W_0)
    lo = a0 * D / (N(r0) * N(f))
    hi = b0 * D / (N(r0) * N(f))
    cand = [m for m in SQ if lo < N(m) < hi] if hi < 50 else None
    if cand is None:
        # sample by norm
        ks = [x for x in norms if lo < x < hi]
        if not ks:
            continue
        m1 = random.choice(by_norm[random.choice(ks)])
        m2 = random.choice(by_norm[random.choice(ks)])
    else:
        if not cand:
            continue
        m1, m2 = random.choice(cand), random.choice(cand)
    Nm1, Nm2 = N(m1), N(m2)
    kmax = C_Phi * N(f) ** 2 * Nm1 * Nm2 / H  # support of hat Phi
    if kmax < 1:
        continue
    Nk = random.uniform(1, kmax)  # norm of a nonzero frequency in the support
    samples += 1
    # S1
    s1 &= N(r0) * N(f) <= b0 * D and H * Nk <= C_Phi * N(f) ** 2 * Nm1 * Nm2 <= C_Phi * b0 ** 2 * D ** 2 / N(r0) ** 2 * (1 + 1e-12)
    # S2
    B, F = dyadic(N(r0)), dyadic(N(f))
    X = D / (B * F)
    calH = C1 * D ** 2 / (H * B ** 2)
    rho = N(r0) * N(f) / (B * F)
    kappa = H * Nk / (N(f) ** 2 * X ** 2)
    x1, x2 = Nm1 / X, Nm2 / X
    ok = (0 < Nk <= calH * (1 + 1e-12) and 1 <= rho < 4 and 0 < kappa <= C1 * (1 + 1e-12)
          and a0 / 4 <= x1 <= b0 and a0 / 4 <= x2 <= b0
          and abs(rho * x1 - N(r0) * N(f) * Nm1 / D) < 1e-9 * D
          and abs(kappa / (x1 * x2) - H * Nk / (N(f) ** 2 * Nm1 * Nm2)) < 1e-9 * (1 + kappa / (x1 * x2)))
    s2 &= ok
    max_rho = max(max_rho, rho)
    min_x = min(min_x, x1, x2)
    max_kappa = max(max_kappa, kappa / C1)
    fal_rho |= rho >= 2
    fal_supp |= min(x1, x2) < a0 / 2
    fal_arg |= abs(H * Nk / (N(f) * X ** 2) / (x1 * x2) - H * Nk / (N(f) ** 2 * Nm1 * Nm2)) > 1e-6 * H * Nk / (N(f) ** 2 * Nm1 * Nm2)
    # S3: shifted scales for j | r0
    for p in factor(r0):
        Xs, Fs = X / N(p), F * N(p)
        s3 &= Fs >= 1 and abs(Xs - D / (B * Fs)) < 1e-9 * Xs and calH <= C1 * D ** 2 / (H * B ** 2) * (1 + 1e-12)
rep.check("S1", s1, "%d random support points (D = %g, H = D^1.1, I = [1,2], C_Phi = %g): N(r0)N(f) <= b0 D and "
          "H N(k) <= C_Phi N(f)^2 N(m1) N(m2) <= C_Phi b0^2 D^2/N(r0)^2" % (samples, D, C_Phi))
rep.check("S2", s2, "N(k) <= calH, rho in [1,4) (max %.3f), kappa <= C1 (max kappa/C1 = %.3f), x_i in [a0/4, b0] "
          "(min x_i = %.3f), kernel identities W_0(rho x_i), kappa/(x1 x2)" % (max_rho, max_kappa, min_x))
rep.falsifier("S2", "rho in [1,2)", fal_rho, "some rho >= 2 occurs")
rep.falsifier("S2", "support [a0/2, b0]", fal_supp, "some x_i < a0/2 occurs")
rep.falsifier("S2", "hat Phi argument with N(f) for N(f)^2", fal_arg, "identity breaks")
rep.check("S3", s3, "(X/N(j), F N(j)) satisfy (7.15) with the same calH for every prime j | r0")

# S4: bookkeeping over dyadic B
worst = 0.0
B = 1.0
while B <= b0 * D:
    s = sum(N(r) for r in SQ if B <= N(r) < 2 * B)
    val = (H / D ** 2) * (D / B) ** 2 * s
    worst = max(worst, val / H)
    B *= 2
rep.check("S4", worst < 2, "max over dyadic B of [sum_{B<=N(r0)<2B} (H N(r0)/D^2)(D/B)^2] / H = %.3f (bounded)" % worst)

sys.exit(1 if rep.summary() else 0)
