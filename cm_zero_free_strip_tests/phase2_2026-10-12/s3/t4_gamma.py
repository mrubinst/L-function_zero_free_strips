"""t4: the gamma quotient on the reflected line and the exponent bookkeeping of the detector (inserts.tex,
Lemma lem:s3-labels (h) and the proof of Proposition prop:s3-detector).

Claims
  S3-T4a  For 51/100 <= a <= 1, 0 < e <= 1/1000 and real t,
              c <= |Gamma(1/2+a+6e-it) / Gamma(3/2-a-6e+it)| / (3+|t|)^{2a-1+12e} <= C
          with C = 1.1 and c = 0.15 (grid: 50 values of a, 4 of e, 400 log-spaced |t| up to 1e12, both signs).
          Falsifier: the same ratio with exponent 2a-1+12e-0.1 is unbounded (exceeds 10 on the grid).
  S3-T4b  OAI78's finite-order quotient |Gamma(a+6e-it)/Gamma(1-a-6e+it)| obeys the same normalized bound <= 1.1
          (so the constant C of [OAI78, L4319] can be kept), and for |t| >= 10 the CM/finite-order ratio lies in
          [1/2, 2] (same Stirling exponent).  Near t = 0, a = 1 the finite-order quotient is small (about 6e, from the
          pole of Gamma at 0), the CM one is not; only upper bounds are used.
          Falsifier: the CM quotient against |t|^{1/2} times the finite-order one is not bounded.
  S3-T4c  Exponent bookkeeping: (i) reflected estimate (a-1/2+6e) + 6e = a-1/2+12e;
          (ii) witness contour [OAI78, L4568-4589]: max over 1 <= t <= 3/2, 51/100 <= sigma <= 1 of -3 + t(5/4 - sigma)
          is -189/100 (attained at t = 3/2, sigma = 51/100), and Q^{3/5} <= (5184 U)^{3/5} <= U^2 for U >= 5184^{3/7};
          falsifier: the claimed maximum -2 is wrong.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3lib import report, summary
import mpmath as mp
from fractions import Fraction as Fr

mp.mp.dps = 30
A = [mp.mpf(51)/100 + k*(mp.mpf(49)/100)/49 for k in range(50)]
E = [mp.mpf(1)/100000, mp.mpf(1)/10000, mp.mpf(5)/10000, mp.mpf(1)/1000]
T = [mp.mpf(0)] + [mp.power(10, mp.mpf(k)/33) for k in range(-66, 400)]   # 1e-2 .. 1e12


def cm_q(a, e, t):
    return abs(mp.gamma(mp.mpf(1)/2 + a + 6*e - 1j*t)/mp.gamma(mp.mpf(3)/2 - a - 6*e + 1j*t))


def fo_q(a, e, t):
    return abs(mp.gamma(a + 6*e - 1j*t)/mp.gamma(1 - a - 6*e + 1j*t))


hi = mp.mpf(0); lo = mp.mpf(10); fhi = mp.mpf(0); r_hi = mp.mpf(0); r_lo = mp.mpf(10); f2 = mp.mpf(0); fo_hi = mp.mpf(0)
for a in A:
    for e in E:
        ex = 2*a - 1 + 12*e
        for t in T:
            for tt in (t, -t):
                q = cm_q(a, e, tt); w = (3 + abs(tt))**ex
                r = q/w
                hi = max(hi, r); lo = min(lo, r)
                fhi = max(fhi, q/(3 + abs(tt))**(ex - mp.mpf(1)/10))
                fq = fo_q(a, e, tt)
                fo_hi = max(fo_hi, fq/w)
                if abs(tt) >= 10:
                    rr = q/fq
                    r_hi = max(r_hi, rr); r_lo = min(r_lo, rr)
                    g = q/(fq*(1 + abs(tt))**mp.mpf(0.5))
                    f2 = max(f2, g, 1/g)
okA = hi <= 1.1 and lo >= 0.15
report("S3-T4a", "lem:s3-labels(h)", "|Gamma(1/2+a+6e-it)/Gamma(3/2-a-6e+it)| / (3+|t|)^{2a-1+12e} in [0.15, 1.1]", okA,
       f"grid {len(A)}x{len(E)}x{2*len(T)}: max {mp.nstr(hi, 6)}, min {mp.nstr(lo, 6)}; falsifier (exponent - 0.1): "
       f"max {mp.nstr(fhi, 6)}", falsifier_failed=fhi > 10)
okB = r_hi <= 2 and r_lo >= 0.5 and fo_hi <= 1.1
report("S3-T4b", "lem:s3-labels(h)", "finite-order quotient obeys the same bound; CM/finite ratio in [1/2,2] for |t|>=10", okB,
       f"finite-order normalized max {mp.nstr(fo_hi, 6)}; ratio max {mp.nstr(r_hi, 6)}, min {mp.nstr(r_lo, 6)}; falsifier (extra factor |t|^(1/2)): max deviation factor "
       f"{mp.nstr(f2, 6)}", falsifier_failed=f2 > 100)

# exponent bookkeeping (exact)
a, e = Fr(3, 4), Fr(1, 2000)
ok_i = (a - Fr(1, 2) + 6*e) + 6*e == a - Fr(1, 2) + 12*e
grid_t = [Fr(1) + Fr(k, 100) for k in range(51)]
grid_s = [Fr(51, 100) + Fr(k, 100) for k in range(50)]
mx = max(-3 + t*(Fr(5, 4) - s) for t in grid_t for s in grid_s)
ok_ii = mx == Fr(-189, 100)
import math
U0 = 5184**(3/7)
ok_iii = all((5184*U)**0.6 <= U**2 for U in [U0*1.0001, 1e3, 1e6, 1e12])
report("S3-T4c", "prop:s3-detector proof", "reflected exponent a-1/2+12e; witness contour exponent max -189/100; "
       "Q^{3/5} <= U^2", ok_i and ok_ii and ok_iii,
       f"max = {mx}; falsifier 'max = -2' true? {mx == Fr(-2)}", falsifier_failed=(mx != Fr(-2)))
summary()
