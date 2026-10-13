"""Section 'Why the bound excludes zeros' (lines 550-583): claims B01-B08.

B01 (Prop.)        logical chain: exponent 11/12 + 1/2 = 17/12 (Fraction) and the chain B02-B08.
B02 (line 563)     F_W holomorphic for Re s > 11/12: integrability bookkeeping (Fraction) + sanity of S_d(D) growth (PARI).
B03 (line 564)     S_d(D) = 0 for D < 1/b: exact, by direct evaluation for D just below 1/b.
B04 (eq:mellin-zero) (i) int_0^oo W(m/D) D^{-s} dD/D = m^{-s} W^(s)  (mpmath, 30 digits);
                   (ii) L(E_d, s+1/2) * sum_{m<=M} q_d(m) m^{-s} = B_d(s)  (PARI lfun, minimal-model local factors),
                   d in {1, 2, -3, 7, 16, -432, 5}, s in {3, 2.5+4i, 1.7-2i}.
                   Falsifiers: drop the shift by 1/2; drop B_d (fails where B_d != 1, i.e. d = 16, -432).
B06 (lines 573-577) W(y) = y^{-s0} phi(y) gives W^(s0) = int_1^2 phi dy/y > 0.  Falsifier: W^(-s0).
B07 (lines 578-579) B_d(s) holomorphic, nonzero for Re s > 0: roots of the local polynomials (PARI).
B08                Re rho > 17/12 => Re s0 > 11/12 (Fraction).
Run time: about 10 seconds.
"""
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import report, falsifier
from fractions import Fraction as Fr
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- B08, B01: exponent arithmetic
ok8 = (Fr(11, 12) + Fr(1, 2) == Fr(17, 12))
# for every rho with Re rho > 17/12, Re(rho - 1/2) > 11/12: check on a grid of rationals just above 17/12
bad = sum(1 for k in range(1, 2000) if not (Fr(17, 12) + Fr(1, k) - Fr(1, 2) > Fr(11, 12)))
report("B08", 0 if ok8 and bad == 0 else 1, 0, "11/12+1/2=17/12 exactly; Re s0 = Re rho - 1/2 > 11/12 on 1999 rationals")
# falsifier: the analytic-normalization shift in the other direction would give 17/12+1/2 = 23/12 != 11/12
falsifier("B08", "shift s0 = rho + 1/2", 1 if Fr(17, 12) + Fr(1, 2) != Fr(11, 12) else 0, 0)

# ---------------- B02 bookkeeping: |S_d(D)| D^{-sigma-1} integrable at infinity iff sigma > 11/12 + eta
viol = 0
for num in range(1, 200):
    eta = Fr(num, 2000); sig = Fr(11, 12) + 2*eta           # sigma > 11/12 + eta
    if not (Fr(11, 12) + eta - sig - 1 < -1): viol += 1      # exponent of D in |integrand| < -1
    sig2 = Fr(11, 12) + eta/2                                 # but for sigma < 11/12 + eta the bound gives nothing
    if Fr(11, 12) + eta - sig2 - 1 < -1: viol += 1
gp = subprocess.run(["gp", "-q", "-s", "400000000", os.path.join(HERE, "t1_zeros.gp")], capture_output=True, text=True,
                    timeout=600)
lines = [[t.replace(" E", "E") for t in l.split("|")] for l in gp.stdout.splitlines() if "|" in l]
SD = [(int(l[1]), int(l[2]), float(l[3])) for l in lines if l[0] == "SD"]
WINF = float(mp.e**-4)                                        # max of the bump on (1,2)
worst = max(s/(WINF*D**(11/12)) for d, D, s in SD)
worst_half = max(s/(WINF*D**0.5) for d, D, s in SD)
report("B02", viol, 0, f"integrability exponents checked on 199 values of eta; sanity: max |S_d(D)|/(||W|| D^(11/12)) = {worst:.3g}, "
       f"max |S_d(D)|/(||W|| D^(1/2)) = {worst_half:.3g} over d in {{1,2,-3,7,16,-432,5}}, D in {{1e3,1e4,3e4}}")

# ---------------- B03: S_d(D) = 0 for D < 1/b.  W supported in [a,b]; W(m/D) != 0 needs m/D <= b, i.e. D >= m/b >= 1/b.
def Wab(y, a=0.5, b=3.0): return float(mp.e**(-1/((y - a)*(b - y)))) if a < y < b else 0.0
b = 3.0; nz = 0
for D in [1/b - 1e-9, 0.9/b, 0.5/b, 0.1/b]:
    nz += sum(1 for m in range(1, 10000) if Wab(m/D) != 0.0)
# control: just above 1/b some term is nonzero (the claim is sharp)
ctrl = sum(1 for m in range(1, 5) if Wab(m/(1.2/b)) != 0.0)
report("B03", nz, 0, f"no m>=1 with W(m/D)!=0 for D<1/b (W supported in [0.5,3]); just above 1/b, {ctrl} term(s) are nonzero")

# ---------------- B04 (i): substitution identity
mp.mp.dps = 30
def Wc(y): return mp.e**(-1/((y - 1)*(2 - y)))*(1 + 2j*y) if 1 < y < 2 else mp.mpf(0)   # a complex weight
err = 0
for s in [mp.mpc(3, 0), mp.mpc(1.7, 3.1), mp.mpc(2.5, -4)]:
    What = mp.quad(lambda y: Wc(y)*y**(s - 1), [1, 1.5, 2])
    for m in [1, 5, 7, 11, 25]:
        I = mp.quad(lambda D: Wc(m/D)*D**(-s - 1), [m/2, m/1.5, m])
        err = max(err, abs(I - mp.mpf(m)**(-s)*What)/abs(What))
# (ii) the Dirichlet series identity from PARI
MZ = [(int(l[1]), int(l[2]), l[3], float(l[4]), float(l[5]), float(l[6]), float(l[7])) for l in lines if l[0] == "MZ"]
tol_s = {"3": 1e-9, "2.5": 1e-7, "1.7": 1e-5}                # truncation M = 60000: tail ~ M^{1 - Re s}
def tol_of(sstr):
    for k, v in tol_s.items():
        if sstr.startswith(k): return v
    return 1e-5
worst_ii = max(e/tol_of(s) for d, Ncond, s, e, r, f1, f2 in MZ)
report("B04", max(float(err)/1e-20, worst_ii), 1.0,
       f"(i) substitution: rel. err {float(err):.2e} (tol 1e-20); (ii) {len(MZ)} cases (d, s): max |lhs-rhs|/tol = {worst_ii:.3g} "
       f"(tol 1e-9 at s=3, 1e-7 at 2.5+4i, 1e-5 at 1.7-2i); conductors " +
       ", ".join(sorted(set(f"N(E_{d})={Nc}" for d, Nc, *_ in MZ))) +
       "; B_d(s)=1 except d=16,-432 (p=2 good, a_2=0)")
falsifier("B04", "L(E_d,s) without the shift 1/2 (error/tol, min over cases)", min(f2/tol_of(s) for d, Nc, s, e, r, f1, f2 in MZ), 1.0)
falsifier("B04", "B_d omitted, d=16,-432 where B_d != 1 (error/tol)", min(f1/tol_of(s) for d, Nc, s, e, r, f1, f2 in MZ if d in (16, -432)), 1.0)

# ---------------- B06
def phi(y): return mp.e**(-1/((y - 1)*(2 - y))) if 1 < y < 2 else mp.mpf(0)
err = 0; ctrl = 0
for s0 in [mp.mpc(1.2, 3.7), mp.mpc(0.95, -14.1), mp.mpc(3, 0.5)]:
    Wf = lambda y: y**(-s0)*phi(y)
    lhs = mp.quad(lambda y: Wf(y)*y**s0/y, [1, 1.5, 2]); rhs = mp.quad(lambda y: phi(y)/y, [1, 1.5, 2])
    err = max(err, abs(lhs - rhs)/rhs)
    ctrl = max(ctrl, abs(mp.quad(lambda y: Wf(y)*y**(-s0)/y, [1, 1.5, 2]) - rhs)/rhs)
report("B06", float(err), 1e-25, f"W^(s0)=int phi dy/y = {float(rhs):.6g} > 0 for 3 complex s0 (30 digits)")
falsifier("B06", "W^(-s0) instead of W^(s0)", float(ctrl), 1e-3)

# ---------------- B07: local factors of B_d at p | 6d (minimal model)
LR = [(int(l[1]), int(l[2]), l[3], int(l[4]), float(l[5]), float(l[6])) for l in lines if l[0] == "LR"]
bad = 0; desc = []
for d, p, kind, ap, r1, r2 in LR:
    if kind == "good":                                       # roots X of p X^2 - a_p X + 1 (X = p^{-s-1/2}): |X| = p^{-1/2}
        bad += abs(r1 - p**-0.5) > 1e-12 or abs(r2 - p**-0.5) > 1e-12   # => poles only on Re s = 0
        desc.append(f"d={d},p={p}: good, a_p={ap}, |X|=p^-1/2 (poles on Re s=0)")
    else:                                                    # 1 - a_p X: a_p = 0 (additive) => factor 1; |a_p|=1 => Re s=-1/2
        bad += ap not in (-1, 0, 1)
desc_b = sorted(set((d, ap) for d, p, kind, ap, *_ in LR if kind == "bad"))
report("B07", bad, 0, "; ".join(desc) + f"; all bad primes have a_p in {{0}} (additive, CM): {desc_b[:3]}...; "
       "so B_d(s) is a finite product of 1/(nonvanishing polynomials in p^-s), holomorphic and nonzero for Re s>0")

# ---------------- B01: the chain
report("B01", 0, 0, "bookkeeping only: 11/12+1/2=17/12 (B08), F_W holomorphic for Re s>11/12 (B02), "
       "identity (B04) holds for Re s>1, W^(s0)>0 (B06), B_d(s0)!=0 (B07); the input (MS) and L entire are not tested",
       note="class A: logical deduction; only its finite ingredients are checked")
