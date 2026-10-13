"""t8_beta_block.py -- eq:s4-beta of cm-eleven-twelfths-guide.tex (Appendix E, Step 2 of the proof of Lemma
lem:s4-canonical): when the completion on the right of (s4-inversion) is reopened, with the cube index b = h h' and
the dyads of h and h' isolated, the coefficient of the block (OAI78 eq.(A), L9984) is
      beta^CM(b, f) = conj(alpha(b))^3 beta(b, f),
where beta is OAI78's coefficient (cube coefficient conj(alpha(b))^3 Psi_k(b)^3 there, conj(alpha(b))^6 Psi_k(b)^3 in
the CM completed sum eq:s4-T).  beta^CM is independent of the row k and of n, and contains the mask 1_{(b,f)=1}.

Test.  Exact finite sums.  For each pair (I1, I2) of norm ranges for (h, h'), the corresponding part of the right side
of (s4-inversion),
   sum_{h in I1} mu(h) conj(alpha(h))^6 Psi_k(h)^3 / N(h) * T_{d_h}(X/N(h)^3; k) restricted to h' in I2,
with T of eq:s4-T (computed from its definition, n squarefree, h' any primary), is compared with the block form
   sum_b beta^CM(b,f) chi_b(k)^3 * sum_n a(n) chi_n(k) chi_n(f)^4 d(n b^3) W(N(n) N(b)^3 / X),
   beta^CM(b,f) = conj(alpha(b))^3 * beta(b,f),
   beta(b,f) = X^{-1/2} N(b)^{1/2} c_{I1,I2}(b) conj(alpha(b))^3 nu(b)^3 rho(b) 1_{(b,f)=1},   (OAI78 L9967-9990)
   c_{I1,I2}(b) = sum over squarefree h | b with h in I1 and b/h in I2 of mu(h),
with ONE beta^CM for all rows k (so it is independent of k and of n).  W is a complex bump on (1/2, 2),
V_*(y) = y^{1/2} W(y), nu = chi_.(lambda), puncture at a prime of norm 37, f of norm 43.  Two scales:
  X = 3000:   dyads N <= 4 / N > 4 for h and for h' (a partition), product mark with two slots;
  X = 200000: h = 1 / N(h) >= 7, and 7 <= N(h') <= 20 / N(h') > 20 (h' != 1, so the n-sums are short), one-slot mark;
              here b runs over 22 elements != 1, including b = p^2 with N(p) = 7.
Errors are relative to each nonzero piece.  Falsifiers (must fail on every piece containing a b != 1):
beta^CM = beta (zeta cube factor); beta^CM = conj(alpha(b))^6 beta.
Also checked: at X = 3000 the four pieces add up to the left side of (s4-inversion) (so the pieces are a partition).
Usage: python3 t8_beta_block.py
"""
import os, sys, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'part_B_lines_0550_0997')))
from cmlib import *

def al(n): return 1.0 + 0j if n == (1, 0) else alpha_c(n)
def phi(n): return al(n).conjugate()
def muk(n): return 1 if n == (1, 0) else mu(n)
_g2 = {}
def g2(n):
    if n == (1, 0): return 1.0 + 0j
    if n not in _g2: _g2[n] = gamma_j_fast(2, n)
    return _g2[n]
def chi0(n, x): return 1.0 + 0j if n == (1, 0) else complex(chi(n, x))
def coprime(a, b): return N(gcd_primary(a, b)) == 1
def bumpW(y):
    return math.exp(-1/((y - 0.5)*(2.0 - y)))*cmath.exp(0.9j*y) if 0.5 < y < 2.0 else 0.0
def Vstar(y): return math.sqrt(y)*bumpW(y)

PRIMES = [p for p in PR_ODD if N(p) <= 200]
P7 = [p for p in PRIMES if N(p) == 7]; P13 = [p for p in PRIMES if N(p) == 13]
p19, p31, p37, p43 = [[p for p in PRIMES if N(p) == q][0] for q in (19, 31, 37, 43)]
LAMB = (1, 2)
def nu(n): return chi0(n, LAMB)
SLOTS_A = [[(P7[0], 0.8 + 0.3j), (p19, -0.5 + 0.6j)], [(P13[0], 0.7 - 0.2j), (p31, 0.4 + 0.9j)]]   # two slots
SLOTS_B = [[(P7[0], 0.8 + 0.3j), (P13[0], -0.5 + 0.6j)]]                                           # one slot
SLOTS = list(SLOTS_A)
def mark(A):
    tot = 1.0 + 0j
    for lst in SLOTS:
        tot *= sum(c for p, c in lst if divides(p, A))
    return tot

f = p43; rp = p37
def rho(n): return 1 if coprime(n, rp) else 0
def Psi(n, k): return 1.0 + 0j if n == (1, 0) else nu(n)*rho(n)*chi0(n, k)*chi0(n, f)**4
ROWS = [(4, -3), (5, 3), (2, 7)]

def run(X, I1, I2, full):
    """All pieces (I1 x I2) at scale X.  full=True: the dyads cover everything, so the pieces must add up to the
    left side of (s4-inversion).  Returns (max rel err of the claim, falsifier errs, partition err, #b, #nonsqf b)."""
    SQF = primaries(int(2*X/(1 if full else 343)) + 1, sqfree=True)
    ALLP = primaries(int((2*X)**(1/3)) + 2)
    def lhs_inversion(k):
        return X**-0.5*sum(phi(n)**2*g2(n)*Psi(n, k)*mark(n)*bumpW(N(n)/X) for n in SQF if N(n) < 2*X)
    def piece_direct(k, i1, i2):
        tot = 0j
        for h in ALLP:
            if muk(h) == 0 or not i1(h) or N(h)**3 >= 2*X: continue
            Y = X/N(h)**3; h3 = pw(h, 3); T = 0j
            for hp in ALLP:
                if not i2(hp) or N(hp)**3 >= 2*Y: continue
                hp3 = pw(hp, 3)
                for n in SQF:
                    if N(n)*N(hp)**3 >= 2*Y: break
                    v = Vstar(N(n)*N(hp)**3/Y)
                    if v == 0: continue
                    T += (phi(n)**2*g2(n)*Psi(n, k)*phi(hp)**6*Psi(hp, k)**3/(math.sqrt(N(n))*N(hp))
                          *mark(mul(mul(h3, n), hp3))*v)
            tot += muk(h)*phi(h)**6*Psi(h, k)**3/N(h)*T
        return tot
    def c_coef(b, i1, i2):
        return sum(muk(h) for h in ALLP if muk(h) != 0 and i1(h) and divides(h, b) and i2(div(b, h)))
    def beta_oai(b, i1, i2):      # OAI78's beta(b,f): cube coefficient conj(alpha(b))^3 nu^3 rho 1_{(b,f)=1}, chi_b(k)^3 removed
        mask = 1 if coprime(b, f) else 0
        return X**-0.5*math.sqrt(N(b))*c_coef(b, i1, i2)*phi(b)**3*nu(b)**3*rho(b)*mask
    def piece_block(k, i1, i2, power, skip_b1=False):
        tot = 0j
        for b in ALLP:
            if N(b)**3 >= 2*X or (skip_b1 and b == (1, 0)): continue
            bet = phi(b)**power*beta_oai(b, i1, i2)      # power 3: the claim eq:s4-beta
            if bet == 0: continue
            b3 = pw(b, 3); inner = 0j
            for n in SQF:
                if N(n)*N(b)**3 >= 2*X: break
                w = bumpW(N(n)*N(b)**3/X)
                if w == 0: continue
                inner += phi(n)**2*g2(n)*nu(n)*rho(n)*chi0(n, k)*chi0(n, f)**4*mark(mul(n, b3))*w
            tot += bet*chi0(b, k)**3*inner
        return tot
    err = perr = 0.0; f0 = f6 = float("inf")
    for k in ROWS:
        total = 0j
        sizes = []
        ds = {}
        for a1, i1 in I1.items():
            for a2, i2 in I2.items():
                ds[(a1, a2)] = piece_direct(k, i1, i2); total += ds[(a1, a2)]
        for (a1, a2), d in ds.items():
            i1, i2 = I1[a1], I2[a2]
            pb = piece_block(k, i1, i2, 3)
            if abs(d) < 1e-14:                       # empty piece: both sides must vanish
                err = max(err, abs(pb)); continue
            err = max(err, abs(d - pb)/abs(d))      # relative to this piece
            if abs(piece_block(k, i1, i2, 3, skip_b1=True)) < 1e-14: continue   # only b = 1: no angular factor
            f0 = min(f0, abs(d - piece_block(k, i1, i2, 0))/abs(d))   # every nonzero piece must fail
            f6 = min(f6, abs(d - piece_block(k, i1, i2, 6))/abs(d))
            print(f"      X={X:g} k={k} (h,h') in ({a1},{a2}): |piece| = {abs(d):.3e}")
        if full:
            perr = max(perr, abs(total - lhs_inversion(k))/abs(lhs_inversion(k)))
    bs = [b for b in ALLP if N(b)**3 < 2*X and b != (1, 0)]
    nsq = sum(1 for b in bs if len(set(factor(b))) < len(factor(b)))
    return err, f0, f6, perr, len(bs), nsq

# run 1: X = 3000, two dyads for h and for h' covering everything (pieces add up to the left side)
r1 = run(3000.0, {"N<=4": lambda n: N(n) <= 4, "N>4": lambda n: N(n) > 4},
         {"N<=4": lambda n: N(n) <= 4, "N>4": lambda n: N(n) > 4}, True)
SLOTS[:] = SLOTS_B          # one-slot mark, so that b = p^2 with p in the slot is not killed by the mark
# run 2: X = 200000, h' != 1 (so b != 1 and the n-sums are short), dyads chosen so that b = p^2 (norm 49) and
# products of two primes occur
r2 = run(200000.0, {"h=1": lambda n: N(n) == 1, "N(h)>=7": lambda n: N(n) >= 7},
         {"7<=N<=20": lambda n: 7 <= N(n) <= 20, "N>20": lambda n: N(n) > 20}, False)
err = max(r1[0], r2[0]); f0 = min(r1[1], r2[1]); f6 = min(r1[2], r2[2])
report("T8.8", err, 1e-11, f"eq:s4-beta: block coefficient beta^CM = conj(alpha(b))^3 beta(b,f), one beta for 3 rows k; "
       f"X = 3000 (4 dyad pairs, {r1[4]} b != 1) and X = 200000 ({r2[4]} b != 1, {r2[5]} not squarefree); relative to each "
       "nonzero piece (falsifiers must fail on every piece with a b != 1 term); two-slot mark (X = 3000) and one-slot mark (X = 200000), puncture, nu = chi_.(lambda)")
report("T8.8-partition", r1[3], 1e-11, "X = 3000: the four pieces add up to the left side of (s4-inversion) (relative)")
falsifier("T8.8", "beta^CM = beta (zeta cube factor)", f0, 1e-3)
falsifier("T8.8", "beta^CM = conj(alpha(b))^6 beta", f6, 1e-3)
