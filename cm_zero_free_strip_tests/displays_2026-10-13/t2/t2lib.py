"""Helpers for the t2 display tests (Sections 7-9 of cm-eleven-twelfths-guide.tex, 2026-10-13).

Imports the part-B helpers read-only (cmlib, common, zw, fastgauss in tests_cm/part_B_lines_0550_0997), and adds an
independent character table of the group frak G = {x in (O/72O)^x : x = 1 mod 3} of the guide (Section 7.1).

The character table is built from an explicit basis of frak G: a greedy search for elements g_1, ..., g_r whose
cyclic subgroups intersect trivially and together generate frak G, so that every element is uniquely
prod g_i^{e_i}, and xi_k(x) = exp(2 pi i sum_i k_i e_i/ord(g_i)).  This is a different construction from the
regular-representation eigenvectors of part B (t3_gauss.py, B26), so the two tables are independent.

Run scripts with `python3 -B` so that nothing is written into the part-B folder.
"""
import os, sys, math, cmath, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LIBB = os.path.normpath(os.path.join(HERE, "..", "..", "part_B_lines_0550_0997"))
sys.path.insert(0, LIBB)
from cmlib import *          # numpy as np, zw helpers, chi, mu, factor, primaries, gamma_j, gamma_j_fast, G_fast, chi_vec, report, falsifier
import numpy as np

M = 72

def red(x, m=M):
    return (x[0] % m, x[1] % m)

def mulm(x, y, m=M):
    return red(mul(x, y), m)

# ------------------------------------------------------------------ the group frak G
ELS = sorted(red((a, b)) for a in range(M) for b in range(M)
             if a % 3 == 1 and b % 3 == 0 and (a*a - a*b + b*b) % 2 == 1)
ONE = (1, 0)
assert len(ELS) == 432, len(ELS)
IDX = {e: i for i, e in enumerate(ELS)}

def order(g):
    k, x = 1, g
    while x != ONE:
        x = mulm(x, g); k += 1
    return k

def cyclic(g):
    out, x = [ONE], g
    while x != ONE:
        out.append(x); x = mulm(x, g)
    out = [ONE] + out[1:]
    return out

def _basis():
    S = {ONE}; basis = []
    cand = sorted(ELS, key=lambda g: -order(g))
    while len(S) < len(ELS):
        for g in cand:
            if g in S: continue
            C = cyclic(g)
            if any(c in S for c in C[1:]): continue      # need <g> cap S = {1}
            newS = {mulm(s, c) for s in S for c in C}
            if len(newS) == len(S)*len(C):
                basis.append((g, len(C))); S = newS; break
        else:
            raise RuntimeError("greedy basis failed")
    return basis

BASIS = _basis()
ORDERS = [o for _, o in BASIS]
assert math.prod(ORDERS) == 432
# exponent vectors of every element
EXPS = {}
for es in itertools.product(*[range(o) for o in ORDERS]):
    x = ONE
    for (g, _), e in zip(BASIS, es):
        for _ in range(e): x = mulm(x, g)
    assert x not in EXPS
    EXPS[x] = es
assert len(EXPS) == 432
EXPARR = np.array([EXPS[e] for e in ELS])                       # [element, i]
KS = list(itertools.product(*[range(o) for o in ORDERS]))      # character labels
KARR = np.array(KS)
PH = (KARR[:, None, :]*EXPARR[None, :, :]/np.array(ORDERS)[None, None, :]).sum(axis=2)
CHARS = np.exp(2j*np.pi*PH)                                     # CHARS[k, element index] = xi_k(element)
TRIV = KS.index(tuple(0 for _ in ORDERS))
INV = {e: [f for f in ELS if mulm(e, f) == ONE][0] for e in ELS}

def cls(n):
    """Index in frak G of the class of a primary n coprime to 2 (n = 1 mod 3 is automatic)."""
    return IDX[red(n)]

# ------------------------------------------------------------------ G on frak G (via classes mod 4)
def _Gmod4():
    reps = {}
    for n in primaries(800, sqfree=True):
        c = (n[0] % 4, n[1] % 4)
        reps.setdefault(c, []).append(n)
    tab = {}; spread = 0.0
    for c, ns in reps.items():
        vals = [(1.0 + 0j) if n == (1, 0) else G_fast(n) for n in ns[:3]]
        spread = max(spread, max(abs(v - vals[0]) for v in vals))
        tab[c] = vals[0]
    return tab, spread

GMOD4, GSPREAD = _Gmod4()
assert len(GMOD4) == 12 and GSPREAD < 1e-10, (len(GMOD4), GSPREAD)
GVEC = np.array([GMOD4[(e[0] % 4, e[1] % 4)] for e in ELS])            # G(t), t in frak G
GINV = np.array([GMOD4[(INV[e][0] % 4, INV[e][1] % 4)] for e in ELS])  # G(t^{-1})

def Gcls(t):
    """G of an element t of frak G given as a residue pair mod 72."""
    return GMOD4[(t[0] % 4, t[1] % 4)]

def Gq(a, b):
    """G(a b^{-1}) for primary a, b coprime to 6 (on classes)."""
    return Gcls(mulm(red(a), INV[red(b)]))

# c_xi = |G|^-1 sum_t G(t^{-1}) conj(xi(t))  (guide (7.8), literally)
CXI = (CHARS.conj() @ GINV)/432.0

# ------------------------------------------------------------------ coefficients
_g = {}
def gam(j, n):
    """gamma_j(n) (vectorized Gauss sum of part B; checked against brute force in t3 and in t2a)."""
    if n == (1, 0): return 1.0 + 0j
    key = (j % 6, n)
    if key not in _g: _g[key] = gamma_j_fast(j, n)
    return _g[key]

def muk(n): return 1 if n == (1, 0) else mu(n)
def Acm(n): return alpha_c(n).conjugate()**2*gam(2, n)      # conj(alpha)^2 gamma_2
def Azeta(n): return alpha_c(n).conjugate()*gam(2, n)       # zeta coefficient conj(alpha) gamma_2
def Gn(n): return 1.0 + 0j if n == (1, 0) else complex(chi(n, (4, 0))).conjugate()*gam(3, n)
def Rf(a, b):
    """R(a, b) = chi_b(a)/chi_a(b) for coprime primary a, b coprime to 6 (chi_b(a) = R(a,b) chi_a(b))."""
    return complex(chi(b, a))/complex(chi(a, b))
def coprime(a, b): return N(gcd_primary(a, b)) == 1
