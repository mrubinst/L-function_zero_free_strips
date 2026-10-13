"""Shared helpers for the t9 display tests (Appendix F of the CM guide: the plain fourth moment with an
angular factor).  Reuses, read-only, the helpers of tests_cm/phase2_2026-10-12/s5 (s5lib.py: Eisenstein
arithmetic, sextic symbols, e(z), report/falsifier with nonzero exit status on failure; hecke.py: Hecke
characters of infinity type (a/|a|)^e, conductors, Gauss sums).  No bytecode is written into that folder.

Conventions (guide (eq:d007), (eq:d008), (eq:d032)): omega = e^{2 pi i/3}; a + b omega is the pair (a, b);
primary = 1 mod 3; alpha(n) = n/|n|; e(z) = exp(4 pi i Im z / sqrt 3); lambda = 1 + 2 omega.
"""
import sys, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
S5 = os.path.normpath(os.path.join(HERE, '..', '..', 'phase2_2026-10-12', 's5'))
if S5 not in sys.path:
    sys.path.insert(0, S5)
import s5lib, hecke                                   # noqa: E402
from s5lib import (report, falsifier, summary, cx, mul, add, sub, conj, N, OM, SQ3, UNITS, chi_p,   # noqa
                   primes_upto, primary_assoc, is_primary, divides, exact_div, e_add)
from hecke import (from_primary_function, primitive, gauss, elements, residue_reps, canon_index)   # noqa

LAM = (1, 2)
TWO = (-2, 0)

def alpha_el(a):
    z = cx(a); return z / abs(z)

def cubic_mod2(a):
    """(a/2)_3: a mod 2 in F_4^x = {1, w, w^2}; a coprime to 2."""
    return {(1, 0): 1, (0, 1): OM, (1, 1): OM ** 2}[(a[0] % 2, a[1] % 2)]

def sextant_generator(a):
    """The generator of the ideal (a) with argument in [0, pi/3) (one per nonzero ideal)."""
    import cmath, math
    best = None
    for u in UNITS:
        b = mul(a, u); ang = cmath.phase(cx(b)) % (2 * math.pi)
        if ang < math.pi / 3 - 1e-12 or ang > 2 * math.pi - 1e-12:
            best = b
    assert best is not None
    return best
