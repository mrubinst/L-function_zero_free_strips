"""Figures and checks for the subsection on ideals of Z[omega] (quadratic-fields.tex, sec:ideals).

Usage: python3 ideal_figures.py <output directory>
Writes ideal_p7.tex, ideal_ram_inert.tex and ideal_products.tex (TikZ pictures) and asserts every
numerical claim made in that subsection.  Elements are pairs (a, b) = a + b*omega, drawn at
(a - b/2, (sqrt3/2) b) as in lattice.tex.
"""
import math, os, sys, itertools
from math import gcd

def mul(x, y): a, b = x; c, d = y; return (a*c - b*d, a*d + b*c - b*d)
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def conj(x): a, b = x; return (a - b, -b)
def N(x): a, b = x; return a*a - a*b + b*b
def xy(x): a, b = x; return (a - b/2, math.sqrt(3)/2*b)
def divides(d, x):
    n = N(d); q = mul(x, conj(d)); return q[0] % n == 0 and q[1] % n == 0
W = (0, 1)
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
def chi(n): return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)

def window(R):
    """All elements of Z[omega] at distance <= R from 0."""
    K = int(2*R) + 2
    return [(a, b) for a in range(-K, K + 1) for b in range(-K, K + 1) if math.hypot(*xy((a, b))) <= R + 1e-9]

def in_I(m, r, z, c=1):
    """z in c * I_{m,r}, where I_{m,r} = {x + y omega : x + r y = 0 mod m}."""
    x, y = z
    if x % c or y % c: return False
    x //= c; y //= c
    return (x + r*y) % m == 0

# ---------------- checks of the proposition on ideals in coordinates
def rho(m): return sum(1 for r in range(m) if (r*r + r + 1) % m == 0)
for m in range(1, 200):
    for r in range(m):
        if (r*r + r + 1) % m: continue
        # I_{m,r} is closed under multiplication by omega, contains m and omega - r, and has index m
        for z in window(6):
            if in_I(m, r, z): assert in_I(m, r, mul(W, z))
        assert in_I(m, r, (m, 0)) and in_I(m, r, (-r, 1))
        assert len({(x + r*y) % m for x in range(m) for y in range(m)}) == m
        # the map x + y omega -> x + r y mod m is a ring homomorphism
        f = lambda z: (z[0] + r*z[1]) % m
        for z1, z2 in itertools.product(window(2.5), repeat=2):
            assert f(mul(z1, z2)) == f(z1)*f(z2) % m
# every principal ideal (alpha) is c * I_{m,r} with c = gcd of coordinates, N(alpha) = c^2 m
def decompose(alpha):
    """(c, m, r) with (alpha) = c I_{m,r}."""
    c = gcd(alpha[0], alpha[1])
    m = N(alpha)//(c*c)
    beta = (alpha[0]//c, alpha[1]//c)
    # r: omega - r in (beta), i.e. omega = r mod beta; r = -x y^{-1} mod m is found by search
    rs = [r for r in range(m) if divides(beta, (-r, 1))]
    assert len(rs) == 1, (alpha, rs)
    return c, m, rs[0]
for alpha in window(30):
    if alpha == (0, 0): continue
    c, m, r = decompose(alpha)
    assert (r*r + r + 1) % m == 0 and N(alpha) == c*c*m
    for z in window(9):
        assert divides(alpha, z) == in_I(m, r, z, c), (alpha, z)
# number of ideals of norm n: sum_{c^2 | n} rho(n/c^2) = sum_{d | n} chi(d)
for n in range(1, 3001):
    lhs = sum(rho(n//(c*c)) for c in range(1, int(n**.5) + 1) if n % (c*c) == 0)
    assert lhs == sum(chi(d) for d in range(1, n + 1) if n % d == 0), n
# rho at prime powers
for p in [5, 7, 11, 13, 17, 19, 31]:
    for k in (1, 2, 3):
        assert rho(p**k) == (2 if p % 3 == 1 else 0)
assert rho(3) == 1 and rho(9) == 0 and rho(27) == 0
for m in range(2, 300):
    for c in range(2, 300):
        if gcd(m, c) == 1 and m*c < 300: assert rho(m*c) == rho(m)*rho(c)       # rho is multiplicative

# ---------------- the examples of the text
pi, pib = (1, 3), (-2, -3)
assert conj(pi) == pib and mul(pi, pib) == (7, 0)
assert decompose(pi) == (1, 7, 2) and decompose(pib) == (1, 7, 4)
assert (2*2 + 2 + 1) % 7 == 0 and (4*4 + 4 + 1) % 7 == 0 and (-1 - 2) % 7 == 4    # the roots r and r' = -1 - r
assert mul((1, 1), pi) == (-2, 1)                         # omega - 2 = (1 + omega) pi, a generator of p
lam = (1, -1)
assert decompose(lam) == (1, 3, 1)                         # (lambda) = I_{3,1} = {3 | x + y}
assert decompose((2, 0)) == (2, 1, 0)                      # (2) = 2 O
assert decompose((61, 0)) == (61, 1, 0)
assert decompose((-5, -9)) == (1, 61, 13)                 # (61, 13 - omega) of Example ex:euclid
# conjugation maps I_{m,r} to I_{m,-1-r}
for z in window(8):
    assert in_I(7, 2, z) == in_I(7, 4, conj(z))
# product p pbar from two-element generators: (7, w-2)(7, w-4) = (49, 7(w-4), 7(w-2), (w-2)(w-4)) = (7)
gens = [(49, 0), mul((7, 0), (-4, 1)), mul((7, 0), (-2, 1)), mul((-2, 1), (-4, 1))]
assert gens[3] == (7, -7)
assert sub(gens[2], gens[1]) == (14, 0) and 49 - 3*14 == 7
# p^2 = (49, 7(w-2), (w-2)^2), (w-2)^2 = 3 - 5w, p^2 = (pi^2) = (-8 - 3w) = I_{49,30}
assert mul((-2, 1), (-2, 1)) == (3, -5) and mul(pi, pi) == (-8, -3)
assert decompose((-8, -3)) == (1, 49, 30) and (30*30 + 30 + 1) % 49 == 0 and 30 % 7 == 2
assert in_I(49, 30, (3, -5))
assert mul((1, 1), (1, 1)) == W and mul(W, mul(pi, pi)) == (3, -5)    # (w-2)^2 = w pi^2 generates p^2
assert decompose(conj((-8, -3))) == (1, 49, 18) and (-1 - 30) % 49 == 18      # pbar^2 = I_{49,18}
assert mul(W, pi) == (-3, -2) and mul(W, mul(W, pi)) == (2, -1)               # the generators of p
# p^2 != p cap p: p^2 has index 7 in p; three ideals of norm 49
assert sum(1 for c in range(1, 8) if 49 % (c*c) == 0 for _ in range(rho(49//(c*c)))) == 3
# (lambda) p = (lambda) cap p, norm 21
lp = mul(lam, pi)
assert N(lp) == 21
for z in window(9):
    assert divides(lp, z) == (divides(lam, z) and divides(pi, z))
    assert divides((7, 0), z) == (divides(pi, z) and divides(pib, z))
# CRT: O/(7) -> F_7 x F_7, x + y w -> (x + 2y, x + 4y), is a bijection
assert len({((x + 2*y) % 7, (x + 4*y) % 7) for x in range(7) for y in range(7)}) == 49
# lambda^2 = -3 omega generates (3); (lambda) cap (lambda) = (lambda) != (3)
assert mul(lam, lam) == (0, -3)

# ---------------- TikZ
def P(z): x, y = xy(z); return f"({x:.4f},{y:.4f})"
GREY = r"\fill[gray!60] {} circle (1.2pt);"
def dots(pts, style, rad):
    return "\n".join(f"\\fill[{style}] {P(z)} circle ({rad});" for z in pts)
def rings(pts, style, rad):
    return "\n".join(f"\\draw[{style}] {P(z)} circle ({rad});" for z in pts)
def poly(pts, style, fill=False):
    cmd = "\\fill" if fill else "\\draw"
    return f"{cmd}[{style}] " + "--".join(P(z) for z in pts) + "--cycle;"
def label(z, text, pos, style=""):
    return f"\\node[{pos},font=\\scriptsize{(',' + style) if style else ''}] at {P(z)} {{{text}}};"
def axes(R):
    return f"\\draw[gray!35,->] ({-R-0.3:.2f},0)--({R+0.4:.2f},0); \\draw[gray!35,->] (0,{-R*0.95:.2f})--(0,{R*0.95+0.3:.2f});"
def half_open_cell(alpha, beta):
    """Points of Z[omega] in {s alpha + t beta : 0 <= s, t < 1}."""
    out = []
    det = alpha[0]*beta[1] - alpha[1]*beta[0]
    K = int(abs(alpha[0]) + abs(alpha[1]) + abs(beta[0]) + abs(beta[1])) + 2
    for a in range(-K, K + 1):
        for b in range(-K, K + 1):
            # solve (a, b) = s alpha + t beta in coordinates (exact rationals)
            s = (a*beta[1] - b*beta[0])/det; t = (alpha[0]*b - alpha[1]*a)/det
            if -1e-12 <= s < 1 - 1e-12 and -1e-12 <= t < 1 - 1e-12:
                out.append((a, b))
    return out

def fig_p7(path):
    R = 7.3
    pts = window(R)
    I = [z for z in pts if divides(pi, z)]
    gens = [mul(u, pi) for u in UNITS]
    cell = [(0, 0), pi, mul((1, 1), pi), mul(W, pi)]
    inside = half_open_cell(pi, mul(W, pi))
    assert len(inside) == 7
    res = {z: (z[0] + 2*z[1]) % 7 for z in inside}
    assert sorted(res.values()) == list(range(7))           # one point in each residue class
    out = [r"\begin{tikzpicture}[scale=0.62]", axes(R),
           poly(cell, "blue!7", fill=True), poly(cell, "blue!45"),
           poly(sorted(gens, key=lambda z: math.atan2(*xy(z)[::-1])), "red!60!black,thin"),
           "\n".join(GREY.format(P(z)) for z in pts if z not in I),
           dots([z for z in I if z not in gens], "blue!70!black", "2.6pt"),
           dots(gens, "red!70!black", "2.8pt")]
    out.append(dots([z for z in res if z != (0, 0)], "green!45!black", "1.9pt"))
    for z, v in res.items():
        if z != (0, 0):
            out.append(label(z, f"${v}$", "right=-1pt", "green!35!black"))
    out += [label((0, 0), "$0$", "below left"),
            label(pi, r"$1+3\omega$", "above right"),
            label((7, 0), "$7$", "below"),
            label(mul(W, pi), r"$\omega\pi$", "left"),
            label(mul((1, 1), pi), r"$(1+\omega)\pi=\omega-2$", "left"),
            r"\end{tikzpicture}"]
    open(path, "w").write("\n".join(out) + "\n")

def fig_ram_inert(path):
    R = 4.1
    pts = window(R)
    out = []
    for alpha, name, sub_ideal, labs in [
            (lam, r"(\lambda)", (3, 0), {(0, 0): "0", (1, 0): "1", (-1, 0): "-1"}),
            ((2, 0), "(2)", None, {(0, 0): "0", (1, 0): "1", (0, 1): r"\omega", (1, 1): r"1+\omega"})]:
        I = [z for z in pts if divides(alpha, z)]
        gens = [mul(u, alpha) for u in UNITS]
        cellpts = [(0, 0), alpha, mul((1, 1), alpha), mul(W, alpha)]
        inside = half_open_cell(alpha, mul(W, alpha))
        assert len(inside) == N(alpha)
        if alpha == lam:
            assert {(z[0] + z[1]) % 3 for z in inside} == {0, 1, 2}
            inside_lab = {z: {0: "0", 1: "1", 2: "-1"}[(z[0] + z[1]) % 3] for z in inside}
        else:
            assert {(z[0] % 2, z[1] % 2) for z in inside} == {(0, 0), (1, 0), (0, 1), (1, 1)}
            inside_lab = {z: {(0, 0): "0", (1, 0): "1", (0, 1): r"\omega", (1, 1): r"1+\omega"}[(z[0] % 2, z[1] % 2)] for z in inside}
        pic = [r"\begin{tikzpicture}[scale=0.62]", axes(R),
               poly(cellpts, "blue!7", fill=True), poly(cellpts, "blue!45"),
               poly(sorted(gens, key=lambda z: math.atan2(*xy(z)[::-1])), "red!60!black,thin"),
               "\n".join(GREY.format(P(z)) for z in pts if z not in I),
               dots([z for z in I if z not in gens], "blue!70!black", "2.6pt"),
               dots(gens, "red!70!black", "2.8pt")]
        if sub_ideal:
            sq = [z for z in pts if divides(sub_ideal, z)]
            pic.append(rings(sq, "black,thick", "4.2pt"))
        pic.append(dots([z for z in inside_lab if z != (0, 0)], "green!45!black", "1.9pt"))
        for z, v in inside_lab.items():
            if z != (0, 0):
                pic.append(label(z, f"${v}$", "right=-1pt", "green!35!black"))
        pic.append(label((0, 0), "$0$", "below left"))
        if alpha == lam:
            pic += [label(lam, r"$\lambda$", "below right"), label((3, 0), "$3$", "below right"),
                    label(mul(W, lam), r"$\omega\lambda$", "above right")]
        else:
            pic += [label((2, 0), "$2$", "below right"), label((0, 2), r"$2\omega$", "above left")]
        pic.append(r"\end{tikzpicture}")
        out.append("\n".join(pic))
    open(path, "w").write(("\n\\hfill\n").join(out) + "\n")

def fig_products(path):
    R = 7.4
    pts = window(R)
    Ip = [z for z in pts if divides(pi, z)]
    Ipb = [z for z in pts if divides(pib, z)]
    I7 = [z for z in pts if divides((7, 0), z)]
    assert set(I7) == set(Ip) & set(Ipb)
    a = [r"\begin{tikzpicture}[scale=0.4]", axes(R),
         poly(sorted([mul(u, (7, 0)) for u in UNITS], key=lambda z: math.atan2(*xy(z)[::-1])), "black!60,thin"),
         "\n".join(GREY.format(P(z)) for z in pts if z not in Ip),
         dots(Ip, "blue!70!black", "2.7pt"),
         rings(Ipb, "orange!85!black,thick", "4.2pt"),
         label(pi, r"$\pi$", "above left"), label(pib, r"$\bar\pi$", "below left"),
         label((7, 0), "$7$", "below right"), r"\end{tikzpicture}"]
    p2 = mul(pi, pi)
    Ip2 = [z for z in pts if divides(p2, z)]
    assert len(Ip2) == 7
    b = [r"\begin{tikzpicture}[scale=0.4]", axes(R),
         poly(sorted([mul(u, (7, 0)) for u in UNITS], key=lambda z: math.atan2(*xy(z)[::-1])), "black!45,thin"),
         poly(sorted([mul(u, p2) for u in UNITS], key=lambda z: math.atan2(*xy(z)[::-1])), "red!60!black,thin"),
         "\n".join(GREY.format(P(z)) for z in pts if z not in Ip),
         dots(Ip, "blue!70!black", "2.7pt"),
         rings(Ip2, "red!70!black,thick", "4.2pt"),
         label(pi, r"$\pi$", "above left"), label((7, 0), "$7$", "below right"), label(p2, r"$\pi^2$", "below left"),
         r"\end{tikzpicture}"]
    open(path, "w").write("\n".join(a) + "\n\\hfill\n" + "\n".join(b) + "\n")

out = sys.argv[1]
fig_p7(os.path.join(out, "ideal_p7.tex"))
fig_ram_inert(os.path.join(out, "ideal_ram_inert.tex"))
fig_products(os.path.join(out, "ideal_products.tex"))
print("ideal figures written; all checks passed")
