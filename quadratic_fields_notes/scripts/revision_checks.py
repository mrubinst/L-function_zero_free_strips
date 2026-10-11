"""Check the revised cusp construction and derivative conventions.

Run with: python scripts/revision_checks.py
The arithmetic helpers are loaded without running the separate, longer suite.
"""
import ast
import itertools
from pathlib import Path

import sympy as sp

source = Path(__file__).with_name('poisson_theta_checks.py').read_text()
tree = ast.parse(source)
helpers = {'mul', 'add', 'sub', 'neg', 'conj', 'N', 'divides', 'div',
           'hnf', 'canon', 'residues', 'primary', 'mmul', 'minv'}
namespace = {'gcd': __import__('math').gcd, '_hnf': {}}
nodes = [node for node in tree.body
         if isinstance(node, ast.FunctionDef) and node.name in helpers]
exec(compile(ast.Module(body=nodes, type_ignores=[]), '<arithmetic helpers>', 'exec'), namespace)
globals().update({name: namespace[name] for name in helpers})

ONE, ZERO = (1, 0), (0, 0)
THREE, LAMBDA, SQRT_MINUS_THREE = (3, 0), (1, -1), (1, 2)
UNITS = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
IDENTITY = ((ONE, ZERO), (ZERO, ONE))
S = ((ZERO, (-1, 0)), (ONE, ZERO))
GAMMA = {'0': IDENTITY,
         '+': ((ONE, ZERO), ((0, 1), ONE)),
         '-': ((ONE, ZERO), ((-1, -1), ONE))}


def quotient(a, b):
    """Euclidean division, minimizing over nearby lattice points."""
    numerator = mul(a, conj(b))
    norm = N(b)
    floors = (numerator[0] // norm, numerator[1] // norm)
    return min(((floors[0] + i, floors[1] + j)
                for i in (-1, 0, 1, 2) for j in (-1, 0, 1, 2)),
               key=lambda q: N(sub(a, mul(q, b))))


def inverse(a, modulus):
    r0, r1, s0, s1 = a, modulus, ONE, ZERO
    while r1 != ZERO:
        q = quotient(r0, r1)
        r0, r1 = r1, sub(r0, mul(q, r1))
        s0, s1 = s1, sub(s0, mul(q, s1))
    if N(r0) != 1:
        return None
    result = div(s0, r0)
    assert divides(modulus, sub(mul(a, result), ONE))
    return result


counts = {'lambda^2 divides c': 0, 'lambda exactly divides c': 0,
          'lambda coprime to c': 0}
elements = list(itertools.product(range(-3, 4), repeat=2))
for initial_a in elements:
    for initial_c in elements:
        if initial_c == ZERO or inverse(initial_a, initial_c) is None:
            continue
        normalize = initial_a if divides(LAMBDA, initial_c) else initial_c
        unit = next(u for u in UNITS if primary(mul(u, normalize)))
        a, c = mul(unit, initial_a), mul(unit, initial_c)
        if divides(LAMBDA, c):
            delta = inverse(a, mul(THREE, c))
            assert delta is not None
        else:
            delta = mul(THREE, inverse(mul(THREE, a), c))
        b = div(sub(mul(a, delta), ONE), c)
        g = ((a, b), (c, delta))
        assert sub(mul(a, delta), mul(b, c)) == ONE
        if divides(THREE, c):
            H, sigma = IDENTITY, '0'
            counts['lambda^2 divides c'] += 1
        elif divides(LAMBDA, c):
            u = next(x for x in (SQRT_MINUS_THREE, neg(SQRT_MINUS_THREE))
                     if divides(THREE, sub(x, c)))
            H = ((ONE, ZERO), (u, ONE))
            sigma = '-' if u == SQRT_MINUS_THREE else '+'
            # H differs from gamma_sigma by lower unipotents whose
            # parameters lie in Z + 3 O, where theta is invariant.
            left = mmul(H, minv(GAMMA[sigma]))
            assert left[0] == (ONE, ZERO) and left[1][1] == ONE
            assert left[1][0][1] % 3 == 0
            counts['lambda exactly divides c'] += 1
        else:
            u = (a[0] % 3, a[1] % 3)
            H = ((u, (-1, 0)), (ONE, ZERO))
            epsilon = {0: 0, 1: 1, 2: -1}[u[1]]
            sigma = {0: '0', 1: '-', -1: '+'}[epsilon]
            remainder = sub(u, (0, epsilon))
            assert remainder[1] % 3 == 0
            translation = ((ONE, remainder), (ZERO, ONE))
            rational = ((neg(ONE), neg(ONE)), (ONE, ZERO)) if epsilon == 1 else S
            assert H == mmul(translation, mmul(rational, GAMMA[sigma]))
            counts['lambda coprime to c'] += 1
        g1 = mmul(g, minv(H))
        assert all(divides(THREE, sub(g1[i][j], IDENTITY[i][j]))
                   for i in range(2) for j in range(2))
print('[PASS] Exact cusp construction and representative identities:', counts)

x, y = sp.symbols('x y', real=True)
v = sp.symbols('v', positive=True)
c, delta = sp.symbols('c delta', nonzero=True)
zprime = -delta / c - (x - sp.I*y) / (c**2 * (v**2 + x**2 + y**2))
zbarprime = -sp.conjugate(delta/c) - (x + sp.I*y) / (sp.conjugate(c)**2 * (v**2 + x**2 + y**2))
vprime = v / (c*sp.conjugate(c) * (v**2 + x**2 + y**2))
dzbar = lambda f: ((sp.diff(f, x) + sp.I*sp.diff(f, y))/2).subs({x: 0, y: 0})
assert sp.simplify(dzbar(zprime) + 1/(c**2*v**2)) == 0
assert sp.simplify(dzbar(zbarprime)) == 0
assert sp.simplify(dzbar(vprime)) == 0
print('[PASS] Symbolic cusp chain rule: derivatives of zprime, zbarprime, vprime')

q, t = sp.symbols('q t', positive=True)
s = sp.symbols('s', real=True)
old_v = 1/(q*t)
transformed = sp.expand_power_base(old_v**(2*s-3)/(q*t**2), force=True)
expected = q**(2-2*s)*t**(1-2*s)
assert sp.simplify(transformed/expected) == 1
print('[PASS] Mellin substitution factor q^(2-2s) t^(1-2s)')

lam = sp.I*sp.sqrt(3)
assert sp.simplify(lam**-3 / sp.Abs(lam**-3) - sp.I) == 0
print('[PASS] Angular normalization alpha((sqrt(-3))^-3) = i')
print('All revision checks passed.')
