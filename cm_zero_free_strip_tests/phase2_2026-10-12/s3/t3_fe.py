"""t3: functional equation and deletion identity for CM members of X_u (inserts.tex, Lemma lem:s3-labels (d), (e), (h)).

For several CM members psi = alpha^sigma theta chi_.(u)^sigma (rows u with one to three good primes, valuations 1..3,
theta in the model group of t2) we
  * compute the conductor f of the inducing character psi* (local exponent search, as in t2), Q = N(f);
  * compute psi* on every prime ideal: on good primes prime to u from the presentation; 0 on primes dividing f;
    at the primes 2, lambda when they do not divide f by psi*((p)) = alpha(x)^sigma / psi(y), where y is good primary,
    prime to u, and x = eps*p*y = 1 mod f for a unit eps;
  * build the Dirichlet coefficients a(m) = sum_{N(I) = m} psi*(I) and ask PARI (lfuncreate/lfunrootres/lfuncheckfeq)
    whether Lambda(s) = (3Q)^{s/2} Gamma_R(s+1/2) Gamma_R(s+3/2) L(s, psi*)  (a constant times (2pi)^{-s} Gamma(s+1/2))
    satisfies Lambda(s) = w Lambda(1-s, conj psi*) with |w| = 1 and no poles.

Claims
  S3-T3a (h)  the functional equation of Lemma lem:cm-growth holds for these psi* with N = 3Q, gamma factor
              Gamma(s+1/2), |w| = 1, and no poles (lfuncheckfeq <= -25 bits).  For sigma = +1 the inducing character is in
              conj(F_CM); the same equation holds (PARI uses the conjugate dual).
              Falsifiers: gamma factor of a finite-order character (Vga = [0, 1]), conductor Q instead of 3Q, and 12Q:
              each gives lfuncheckfeq >= -12.
  S3-T3b      control: a finite-order member theta chi_.(u) has the finite-order equation (Vga = [0,1], N = 3Q).
  S3-T3c (e)  deletion identity: L_orig(s, psi) = sum over good primary n of psi(n) N(n)^{-s} equals
              L(s, psi*) * prod_{p in S} (1 - psi*(p) N(p)^{-s}) at s = 2.5 + 1.3i (truncated at N(n) <= 4e4; tail < 1e-4),
              with no deleted factor at the primes of u (psi*(p) = 0 there).  Falsifier: omitting the factor at the
              unramified S-prime(s) changes the value by more than 1e-3 (only for cases with an unramified S-prime).
"""
import os, sys, random, subprocess, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3lib import *
from sympy import primerange

SCR = __import__('tempfile').mkdtemp(prefix='s3fe_')
os.makedirs(SCR, exist_ok=True)
rng = random.Random(7)
t0 = time.time()


def prime_ideals(B):
    """(generator, norm) for every prime ideal of norm <= B: lambda, -2 (norm 4), inert -q (norm q^2), split primary."""
    out = [(LAM, 3)]
    if 4 <= B:
        out.append(((-2, 0), 4))
    for p in primerange(5, B + 1):
        p = int(p)
        if p % 3 == 2:
            if p*p <= B:
                out.append(((-p, 0), p*p))
        else:
            pi = split_prime(p)
            out.append((pi, p)); out.append((prim_assoc(conj(pi)), p))
    return out


def s_prime_value(member, f, P, sigma, k):
    """psi*((P)) for P in {lambda, -2} not dividing f."""
    for a in range(-60, 61):
        for b in range(-60, 61):
            y = (a, b)
            if not primary(y) or not coprime_to(y, member.bad):
                continue
            for eps in UNITS:
                x = mul(mul(eps, P), y)
                if divides(f, sub(x, (1, 0))):
                    return alpha(x)**(k*sigma)/member(y)
    raise RuntimeError("no auxiliary y found")


def psistar_table(member, f, sigma, k, B):
    tab = []
    for P, q in prime_ideals(B):
        if divides(P, f):
            v = 0
        elif P in (LAM, (-2, 0)):
            v = s_prime_value(member, f, P, sigma, k)
        else:
            v = member(P)          # good primary prime, prime to u (primes of u divide f)
        tab.append((P, q, v))
    return tab


def coefficients(tab, B):
    a = [0j]*(B + 1); a[1] = 1
    for P, q, v in tab:
        if v == 0:
            continue
        for m in range(q, B + 1):
            if m % q == 0:
                a[m] += v*a[m//q]
    return a


def gp_run(script):
    out = subprocess.run(["gp", "-q", "-D", "realprecision=38", "-D", "parisizemax=4000000000"], input=script, capture_output=True, text=True,
                         timeout=1200)
    errs = [l for l in (out.stdout + out.stderr).split("\n") if "***" in l and "Warning" not in l]
    if out.returncode != 0 or errs:
        raise RuntimeError(out.stdout + out.stderr)
    return [l for l in out.stdout.split("\n") if l.strip() and "***" not in l]


def fe_check(a, Ncond, vga, B):
    path = os.path.join(SCR, "coef.gp")
    with open(path, "w") as fh:
        fh.write("a=[" + ",".join(f"{z.real:.17g}{z.imag:+.17g}*I" for z in a[1:B + 1]) + "];\n")
    script = (f"read(\"{path}\");\n"
              f"L0=lfuncreate([a,1,{vga},1,{Ncond},0]); w=lfunrootres(L0)[3];\n"
              f"L=lfuncreate([a,1,{vga},1,{Ncond},w]); print(lfuncheckfeq(L)); print(abs(w));\n")
    out = gp_run(script)
    return int(out[0]), float(out[1])


CASES = [  # (name, u_factors, unit, theta, sigma, k)
    ("abar chibar(pi7)", [(split_prime(7), 1)], (1, 0), (), -1, 1),
    ("abar chi(lam) chibar(pi7 pi13^2)", [(split_prime(7), 1), (split_prime(13), 2)], (1, 0), ((LAM, 1),), -1, 1),
    ("alpha chi(2) chi(-5)", [((-5, 0), 1)], (1, 0), (((2, 0), 1),), 1, 1),
    ("abar chibar(pi7^3)", [(split_prime(7), 3)], (1, 0), (), -1, 1),
    ("alpha chi(omega) chi(omega pi19)", [(split_prime(19), 1)], (0, 1), (((0, 1), 1),), 1, 1),
    ("abar chi(2)^3 chibar(pi13)", [(split_prime(13), 1)], (1, 0), (((2, 0), 3),), -1, 1),
    ("abar chi(lam)^3 chi(-1) chibar(pi7)  [unramified at 2]", [(split_prime(7), 1)], (1, 0),
     ((LAM, 3), ((-1, 0), 1)), -1, 1),
    ("alpha chi(lam)^3 chi(2)^2 chi(-1) chi(pi7)  [unramified at lambda]", [(split_prime(7), 1)], (1, 0),
     ((LAM, 3), ((2, 0), 2), ((-1, 0), 1)), 1, 1),
    ("control: chi(pi7) (finite order)", [(split_prime(7), 1)], (1, 0), (), 1, 0),
]
okA = True; falsA = True; okB = True; okC = True; falsC = True; lines = []; ctl = ""
for name, uf, unit, th, sigma, k in CASES:
    m = Member(uf, unit=unit, theta=th, sigma=sigma, k=k)
    tgt = (lambda x, s=sigma: alpha(x)**s) if k == 1 else (lambda x: 1)
    ex, Q, f = conductor(m, tgt, rng=rng)
    Ncond = 3*Q
    B = int(min(60000, max(6000, 40*math.sqrt(Ncond)*math.log(Ncond + 10))))
    tab = psistar_table(m, f, sigma, k, B)
    a = coefficients(tab, B)
    if k == 1:
        b, w = fe_check(a, Ncond, "[1/2,3/2]", B)
        b1, _ = fe_check(a, Ncond, "[0,1]", B)
        b2, _ = fe_check(a, Q, "[1/2,3/2]", B)
        b3, _ = fe_check(a, 12*Q, "[1/2,3/2]", B)
        okA &= (b <= -25 and abs(w - 1) < 1e-6)
        falsA &= (b1 >= -12 and b2 >= -12 and b3 >= -12)
        lines.append(f"{name}: f exps {ex}, Q={Q}, B={B}: feq {b} bits, |w|={w:.6f}; falsifiers Vga[0,1] {b1}, N=Q {b2}, "
                     f"N=12Q {b3}")
    else:
        b, w = fe_check(a, Ncond, "[0,1]", B)
        b1, _ = fe_check(a, Ncond, "[1/2,3/2]", B)
        okB &= (b <= -25 and abs(w - 1) < 1e-6) and b1 >= -12
        ctl = f"{name}: Q={Q}: feq {b} bits, |w|={w:.6f}; with CM gamma {b1}"
    # deletion identity at s0
    s0 = complex(2.5, 1.3); X = 40000
    lorig = 0j
    for n in [(x, y) for x in range(-240, 241) for y in range(-240, 241)]:
        if primary(n) and 0 < N(n) <= X and coprime_to(n, [LAM, (-2, 0)]):
            lorig += m(n)*N(n)**(-s0)
    Lstar = 1 + 0j; Sfac = 1 + 0j; unram_S = False
    for P, q, v in tab:
        if q <= X:
            Lstar /= (1 - v*q**(-s0))
        if P in (LAM, (-2, 0)):
            Sfac *= (1 - v*q**(-s0))
            unram_S |= (v != 0)
    rhs = Lstar*Sfac
    err = abs(lorig - rhs)
    okC &= err < 1e-4
    if unram_S:
        falsC &= abs(lorig - Lstar) > 1e-3
    lines.append(f"   deletion identity at s=2.5+1.3i: |L_orig - L(psi*) prod_S| = {err:.1e}; unramified S-prime: {unram_S}")

for l in lines:
    print("   " + l)
report("S3-T3a", "lem:s3-labels(h), lem:cm-growth", "FE with (3Q)^{s/2} Gamma(s+1/2), |w|=1, no poles, for CM members of X_u",
       okA, f"{sum(1 for c in CASES if c[5] == 1)} CM members (both sigma)", falsifier_failed=falsA)
report("S3-T3b", "control", "finite-order member satisfies the finite-order FE (Gamma(s)) and not the CM one", okB, ctl,
       falsifier_failed=None)
report("S3-T3c", "lem:s3-labels(e)", "L_orig = L(s,psi*) prod_{p in S}(1 - psi*(p)q^-s); psi*(p) = 0 for p | u", okC,
       "see lines above", falsifier_failed=falsC)
print(f"time {time.time() - t0:.1f}s")
summary()
