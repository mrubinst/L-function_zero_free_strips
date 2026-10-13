"""t2c: the change of variables of Step 3 of Proposition prop:b2-reduction (displays (7.33)-(7.35)) and the W0 identity (7.31).

eq:d048  r0 = g/e0, f = e0 v0, k = e0 h;  inverse e0 = (f, k), v0 = f/e0, g = r0 e0, h = k/e0.
         Enumerated exactly in Z[omega]: the map is a bijection from
             {(g, e0, v0, h): g, v0 squarefree primary coprime to 6, e0 | g, (v0, g) = 1, h != 0, (v0, h) = 1}
         onto {(r0, f, k): r0, f squarefree primary coprime to 6, (r0, f) = 1, k != 0},
         restricted to N(g v0) <= P and N(e0 h) <= K on the left and N(r0 f) <= P, N(k) <= K on the right
         (these bounds correspond under the map, by eq:d049-d050), and the stated inverse is its inverse.
eq:d049  g v0 = r0 f,  h e0^5 v0^4 = k f^4 (exact in Z[omega]),  mu(e0) mu(v0) = mu(f).
eq:d050  N(g)/N(e0) = N(r0),  N(h)/(N(e0) N(v0)^2) = N(k)/N(f)^2 (exact rationals).
eq:d047  conj(W(N(g z1)/D)) W(N(g z2)/D)/sqrt(N(z1 z2)) = (N(g)/D) W0(N(g z1)/D) conj(W0(N(g z2)/D)), W0(x) = x^-1/2 conj(W(x)),
         symbolically (sympy, positive norms and complex values of W).
Falsifiers: f = e0 v0^2; k = h; the condition (v0, h) = 1 dropped; W0(x) = x^{1/2} conj(W(x)).
Run: python3 -B t2c_change_of_variables.py   (about 10 s)
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2lib import *
import sympy as sp
T0 = time.time()
P, K = 300, 60

SQ = primaries(P, sqfree=True)                                    # includes 1
HS = [(a, b) for a in range(-20, 21) for b in range(-20, 21) if 0 < N((a, b)) <= K]
def primes_of(n): return set(factor(n))
def coprime_any(v, h):
    """(v, h) = 1 for squarefree primary v coprime to 6 and any nonzero h."""
    return all(not divides(p, h) for p in primes_of(v))
def gcd_fk(f, k):
    e = (1, 0)
    for p in primes_of(f):
        if divides(p, k): e = mul(e, p)
    return e

def domain(require_vh=True):
    out = []
    for g in SQ:
        for v0 in SQ:
            if N(g)*N(v0) > P or not coprime(g, v0): continue
            for e0 in divisors_sqfree(g):
                for h in HS:
                    if N(e0)*N(h) > K: continue
                    if require_vh and not coprime_any(v0, h): continue
                    out.append((g, e0, v0, h))
    return out
def codomain():
    out = set()
    for r0 in SQ:
        for f in SQ:
            if N(r0)*N(f) > P or not coprime(r0, f): continue
            for k in [(a, b) for a in range(-30, 31) for b in range(-30, 31) if 0 < N((a, b)) <= K]:
                out.add((r0, f, k))
    return out
def fwd(t, variant="cm"):
    g, e0, v0, h = t
    r0 = div(g, e0)
    f = mul(e0, v0) if variant != "fv2" else mul(e0, mul(v0, v0))
    k = mul(e0, h) if variant != "k=h" else h
    return (r0, f, k)
def inv(s):
    r0, f, k = s
    e0 = gcd_fk(f, k); v0 = div(f, e0)
    return (mul(r0, e0), e0, v0, div(k, e0))

Dom = domain(); Cod = codomain()
img = [fwd(t) for t in Dom]
bad_in = sum(1 for s in img if s not in Cod)
inj = len(set(img)) == len(img)
bij = inj and bad_in == 0 and len(set(img)) == len(Cod)
round1 = sum(1 for t, s in zip(Dom, img) if inv(s) != t)
round2 = sum(1 for s in Cod if fwd(inv(s)) != s)
report("eq:d048", (0 if bij else 1) + round1 + round2, 0,
       f"P = {P}, K = {K}: {len(Dom)} tuples (g,e0,v0,h) -> {len(Cod)} triples (r0,f,k); image in codomain: {bad_in == 0}; "
       f"injective: {inj}; inverse e0 = (f,k) etc. fails on {round1} + {round2} elements")
def fal_count(variant=None, require_vh=True):
    dom = Dom if require_vh else domain(require_vh=False)
    im = [fwd(t, variant or "cm") for t in dom]
    return sum(1 for s in im if s not in Cod) + (len(im) - len(set(im))) + abs(len(set(im)) - len(Cod))
falsifier("eq:d048", "f = e0 v0^2", fal_count("fv2"), 0)
falsifier("eq:d048", "k = h (e0 not moved into k)", fal_count("k=h"), 0)
falsifier("eq:d048", "condition (v0, h) = 1 dropped from the domain", fal_count(None, require_vh=False), 0)

# ------------------------------------------------------------------ eq:d049, eq:d050 on every tuple
e49 = e50 = 0; f49 = 0
for t, s in zip(Dom, img):
    g, e0, v0, h = t; r0, f, k = s
    ok1 = mul(g, v0) == mul(r0, f)
    ok2 = mul(h, mul(pw(e0, 5), pw(v0, 4))) == mul(k, pw(f, 4))
    ok3 = muk(e0)*muk(v0) == muk(f)
    e49 += (not ok1) + (not ok2) + (not ok3)
    f49 += mul(h, mul(e0, pw(v0, 4))) != mul(k, pw(f, 4))
    ok4 = Fr(N(g), N(e0)) == N(r0)
    ok5 = Fr(N(h), N(e0)*N(v0)**2) == Fr(N(k), N(f)**2)
    e50 += (not ok4) + (not ok5)
report("eq:d049", e49, 0, f"g v0 = r0 f, h e0^5 v0^4 = k f^4 exactly in Z[omega], mu(e0)mu(v0) = mu(f), on all {len(Dom)} tuples")
falsifier("eq:d049", "h e0 v0^4 = k f^4 (e0 to the first power)", f49, 0)
report("eq:d050", e50, 0, f"N(g)/N(e0) = N(r0) and N(h)/(N(e0)N(v0)^2) = N(k)/N(f)^2 in exact rationals, all {len(Dom)} tuples")

# ------------------------------------------------------------------ eq:d047 symbolically
a, b1, b2, D = sp.symbols('a b1 b2 D', positive=True)
w1, w2 = sp.symbols('w1 w2')                                      # complex values W(a b1/D), W(a b2/D)
def W0sym(x, w, p=sp.Rational(-1, 2)): return x**p*sp.conjugate(w)
lhs = sp.conjugate(w1)*w2/sp.sqrt(b1*b2)
rhs = a/D*W0sym(a*b1/D, w1)*sp.conjugate(W0sym(a*b2/D, w2))
d = sp.simplify(lhs - rhs)
rhs_bad = a/D*W0sym(a*b1/D, w1, sp.Rational(1, 2))*sp.conjugate(W0sym(a*b2/D, w2, sp.Rational(1, 2)))
dbad = sp.simplify(lhs - rhs_bad)
report("eq:d047", 0 if d == 0 else 1, 0, f"sympy: difference simplifies to {d}")
falsifier("eq:d047", "W0(x) = x^{1/2} conj(W(x))", 0 if dbad == 0 else 1, 0)
print(f"t2c done in {time.time() - T0:.1f} s")
