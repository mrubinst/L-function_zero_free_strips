"""t1: the label set Lambda and the members of X_u as functions (inserts.tex, Lemma lem:s3-labels (a), (b), (c), and the
expansion of the truncated inverse in the proof of Proposition prop:s3-detector (iv), [OAI78, L4592-4600]).

Model data: Theta_0 = <chi_.(lambda), chi_.(2), chi_.(omega), chi_.(-1)> (Kummer characters with conductor over 2, 3;
orders 6, 6, 3, 2), target psi = abar * eta with eta = chi_.(lambda).  Members are compared as functions on a sample of
good primary n, including n sharing primes with u (zero convention).

Claims
  S3-T1a (a)  Lambda = {(nu, s): nu in Theta_0 u alpha^s Theta_0} has 4|Theta_0| labels giving 4|Theta_0| distinct
              functions, and conj(phi_{u,nu,s}) = phi_{u,conj nu,-s}: X_u is closed under conjugation.
              Falsifier: the set without the labels (alpha theta, +1) is not closed under conjugation.
  S3-T1b (b)  X_u contains chi(u), conj chi(u), psi conj chi(u) = phi_{u, abar eta, -1}, its conjugate, and the slot
              characters theta conj chi(u); the anti-oriented twist abar theta chi(u) is NOT a member (for any theta),
              i.e. X_u differs from the audit's Theta_0 u abar Theta_0 u alpha Theta_0 (both orientations).
  S3-T1c (c)  complete multiplicativity phi(mn) = phi(m) phi(n) on good primary m, n, zero convention included
              (2000 random pairs, finite and CM members).  Falsifier: replacing the zero convention by
              chi_n(u)^{-1} := 1 at nonunits breaks multiplicativity.
  S3-T1d (iv) for a CM member and D* = 300: sum_{lk = n} mu(l) phi(l) phi(k) V(q_l/D*) equals phi(n) times
              sum_{l | n} mu(l) V(q_l/D*), which is 1 at n = 1 and 0 for 1 < q_n <= D* (all good primary n with
              q_n <= 3 D*).  Falsifier: with the broken zero convention the vanishing fails.
"""
import os, sys, random, itertools, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3lib import *

rng = random.Random(2026)
OM = (0, 1)
p7, p13, p19 = split_prime(7), split_prime(13), split_prime(19)
U_F = [(p7, 1), (p13, 2)]                     # row u = pi7 pi13^2
GENS = [(LAM, 6), ((2, 0), 6), (OM, 3), ((-1, 0), 2)]
THETA0 = [tuple((b, e) for (b, _), e in zip(GENS, ex) if e) for ex in itertools.product(*[range(o) for _, o in GENS])]

# sample of good primary n, including multiples of the primes of u
SAMPLE = []
for x in range(-40, 41):
    for y in range(-40, 41):
        n = (x, y)
        if primary(n) and 1 < N(n) <= 1500 and coprime_to(n, [LAM, (-2, 0)]):
            SAMPLE.append(n)
SAMPLE = SAMPLE[::3] + [mul(p7, p19), mul(p13, (-5, 0)), pw(p7, 2), mul(p13, p13)]
SAMPLE = [prim_assoc(n) for n in SAMPLE]


def vec(m):
    return tuple(complex(round(z.real, 9), round(z.imag, 9)) for z in (m(n) for n in SAMPLE))


def cvec(v):
    return tuple(complex(round(z.real, 9), round(-z.imag, 9)) for z in v)


members = {}
for th in THETA0:
    for s in (1, -1):
        members[(th, s, 0)] = vec(Member(U_F, theta=th, sigma=s, k=0))
        members[(th, s, 1)] = vec(Member(U_F, theta=th, sigma=s, k=1))
vals = set(members.values())
closed = all(cvec(v) in vals for v in vals)
labels_conj = all(members[(tuple((b, (-e) % o) for (b, e), o in zip(th, [dict(GENS)[b] for b, _ in th])), -s, k)]
                  == cvec(members[(th, s, k)]) if all(e for _, e in th) else True
                  for (th, s, k) in members)
part = set(v for (th, s, k), v in members.items() if not (k == 1 and s == 1))
fals_closed = all(cvec(v) in part for v in part)
report("S3-T1a", "lem:s3-labels(a)", "|Lambda| = 4|Theta_0| distinct members; X_u closed under conjugation", closed and
       len(vals) == 4*len(THETA0), f"|Theta_0| = {len(THETA0)}, distinct member functions {len(vals)}, sample of "
       f"{len(SAMPLE)} good primary n; falsifier (drop alpha*theta*chi(u)) closed? {fals_closed}",
       falsifier_failed=not fals_closed)

# membership
eta = ((LAM, 1),)
num = vec(Member(U_F, theta=(), sigma=1, k=0))
den = tuple(complex(round(z.real, 9), round(z.imag, 9)) for z in
            (alpha(n).conjugate()*chi(n, LAM)*chipow(chi(n, Member(U_F).u), -1) for n in SAMPLE))
slots = [vec(Member(U_F, theta=th, sigma=-1, k=0)) for th in THETA0 if all(b in ((2, 0), OM, (-1, 0)) for b, _ in th)]
anti = [tuple(complex(round(z.real, 9), round(z.imag, 9)) for z in
              (alpha(n).conjugate()*Member(U_F, theta=th, sigma=1, k=0).theta_val(n)*chi(n, Member(U_F).u)
               for n in SAMPLE)) for th in THETA0]
okB = num in vals and cvec(num) in vals and den in vals and cvec(den) in vals and all(v in vals for v in slots)
anti_in = sum(1 for v in anti if v in vals)
report("S3-T1b", "lem:s3-labels(b)", "X_u contains numerator, denominator psi conj chi(u), conjugates, slot characters; "
       "anti-oriented abar theta chi(u) not in X_u", okB and anti_in == 0,
       f"{len(slots)} slot-type characters checked; anti-oriented twists found in X_u: {anti_in} of {len(anti)}",
       falsifier_failed=None)

# multiplicativity
class Broken(Member):
    def __call__(self, n):
        c = chi(n, self.u)
        cc = 1 if c == 0 else c**self.sigma            # broken zero convention
        return alpha(n)**self.aexp*self.theta_val(n)*cc

pool = SAMPLE + [mul(p7, n) for n in SAMPLE[:40]] + [mul(p13, n) for n in SAMPLE[:40]]
pool = [prim_assoc(n) for n in pool]
worst = 0.0; broken_worst = 0.0
for _ in range(2000):
    m_, n_ = rng.choice(pool), rng.choice(pool)
    th = rng.choice(THETA0); s = rng.choice((1, -1)); k = rng.choice((0, 1))
    f = Member(U_F, theta=th, sigma=s, k=k); g = Broken(U_F, theta=th, sigma=s, k=k)
    mn = mul(m_, n_)
    worst = max(worst, abs(f(mn) - f(m_)*f(n_)))
    broken_worst = max(broken_worst, abs(g(mn) - g(m_)*g(n_)))
report("S3-T1c", "lem:s3-labels(c)", "phi(mn) = phi(m)phi(n) on good primary m, n (zero convention), finite and CM",
       worst < 1e-9, f"max error {worst:.1e} over 2000 pairs; falsifier (chi^-1 := 1 at nonunits) max error "
       f"{broken_worst:.2f}", falsifier_failed=broken_worst > 0.1)


# truncated inverse coefficient identity
def V(x):
    if x <= 1:
        return 1.0
    if x >= 2:
        return 0.0
    t = x - 1
    a = math.exp(-1/(1 - t)) if t < 1 else 0.0
    b = math.exp(-1/t) if t > 0 else 0.0
    return a/(a + b)


def divisors_primary(n):
    fac = factor_elt(n)
    out = []
    for ex in itertools.product(*[range(v + 1) for _, v in fac]):
        d = (1, 0); mu_ = 1
        for (P, _), e in zip(fac, ex):
            d = mul(d, pw(P, e))
            if e >= 2:
                mu_ = 0
            elif e == 1:
                mu_ = -mu_
        out.append((d, mu_))
    return out


Dst = 300
f = Member(U_F, theta=((LAM, 1), ((2, 0), 2)), sigma=-1, k=1)
g = Broken(U_F, theta=((LAM, 1), ((2, 0), 2)), sigma=-1, k=1)
errs = 0.0; ferr = 0.0; cnt = 0
for x in range(-70, 71):
    for y in range(-70, 71):
        n = (x, y)
        if not (primary(n) and 0 < N(n) <= 3*Dst and coprime_to(n, [LAM, (-2, 0)])):
            continue
        cnt += 1
        tot = 0j; totb = 0j; sv = 0.0
        for d, mu_ in divisors_primary(n):
            if mu_ == 0:
                continue
            k_ = div(n, d)
            tot += mu_*f(d)*f(k_)*V(N(d)/Dst)
            totb += mu_*g(d)*g(k_)*V(N(d)/Dst)
            sv += mu_*V(N(d)/Dst)
        target = f(n)*sv
        errs = max(errs, abs(tot - target))
        if N(n) == 1:
            errs = max(errs, abs(tot - 1))
        elif N(n) <= Dst:
            errs = max(errs, abs(tot))
            ferr = max(ferr, abs(totb))
report("S3-T1d", "prop:s3-detector(iv)", "truncated inverse: coefficient 1 at n = 1, 0 for 1 < q_n <= D*, CM member",
       errs < 1e-9, f"{cnt} good primary n with q_n <= {3*Dst}; max error {errs:.1e}; falsifier (broken zero convention) "
       f"max |coefficient| on 1 < q_n <= D*: {ferr:.2f}", falsifier_failed=ferr > 0.1)
summary()
