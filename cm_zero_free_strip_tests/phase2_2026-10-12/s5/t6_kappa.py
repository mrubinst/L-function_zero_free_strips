"""t6: the kappa-dependent inequalities of OAI78's proof of lem:plain at the endpoint kappa = 3/4
(inserts.tex, Step 5 of the proof of Proposition s5-plain).

OAI78 states uniformity for kappa in [3/4, 1]; the CM application uses kappa = 3/4, where 6 kappa - 1 = 7/2 and
s_kappa = (1 + kappa)/2 = 7/8.  Random sampling of the admissible region checks:
  (3.6)   A >= z, A + (6k-1) z <= M  ==>  z <= M/(6k) <= 2M/9,  k z <= M/6
  L12947  A > 5M/6, (6k-1) z <= M - A  ==>  z < M/21,  A - z > 2L  (L = M/4)
  (3.11)/(comparison-uncentered-margin), with b <= L:  A_comp <= 3M/2 - A + 2z <= 23M/30,
          A_comp + (6k-1) z <= 14M/15
Falsifier: kappa = 0.70 (below the range) violates z < M/21 on the sampled region.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5lib import report, falsifier, summary

def check(kappa, n=200000, seed=7):
    rnd = random.Random(seed)
    viol = {"3.6": 0, "z<M/21": 0, "A-z>2L": 0, "Acomp": 0, "Acomp+slots": 0}
    for _ in range(n):
        M = rnd.uniform(0.05, 3.0)
        z = rnd.uniform(0, M)
        A = rnd.uniform(z, 2 * M)
        if A + (6 * kappa - 1) * z > M or z <= 0:
            continue
        if not (z <= M / (6 * kappa) + 1e-12 and M / (6 * kappa) <= 2 * M / 9 + 1e-12 and kappa * z <= M / 6 + 1e-12):
            viol["3.6"] += 1
        if A > 5 * M / 6:
            L = M / 4
            if not z < M / 21: viol["z<M/21"] += 1
            if not A - z > 2 * L: viol["A-z>2L"] += 1
            Acomp = 3 * M / 2 - A + 2 * z      # upper bound (3.11) without xi
            if not Acomp <= 23 * M / 30 + 1e-12: viol["Acomp"] += 1
            if not Acomp + (6 * kappa - 1) * z <= 14 * M / 15 + 1e-12: viol["Acomp+slots"] += 1
    return viol

v = check(0.75)
report("6-kappa=3/4", sum(v.values()), 0, f"violations {v}")
v1 = check(1.0)
report("6-kappa=1", sum(v1.values()), 0, f"violations {v1}")
report("6-s_kappa", abs((1 + 0.75) / 2 - 7 / 8), 0, "s_kappa = (1+kappa)/2 = 7/8 at kappa = 3/4 (OAI78 thm:main: no finite-order zeros in Re s > 7/8)")
vf = check(0.70)
falsifier("6-kappa=3/4", "kappa = 0.70", vf["z<M/21"] + vf["3.6"], 0)
summary()
