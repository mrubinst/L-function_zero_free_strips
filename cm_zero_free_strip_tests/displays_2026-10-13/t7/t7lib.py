"""Shared helpers for agent t7's display tests (CM guide, Appendix C (Hecke L-functions of infinity type zbar/|z|)
and Appendix D (the zero detector for a CM target)).

Conventions (guide / [QF]): omega = e^{2 pi i/3}; a + b*omega is the integer pair (a, b); zeta = 1 + omega = e^{i pi/3};
the sixth roots of unity are indexed j = 0..5 as zeta^j (zw.UNITS); primary = 1 mod 3; lambda = 1 + 2 omega;
e(z) = exp(4 pi i Im z / sqrt 3); alpha(a) = a/|a|.

A character psi in F_CM is given by its finite part phi, a primitive character of (O/f)^x with phi(zeta) = zeta for the
six units; psi((a)) = phi(a) conj(alpha(a))  (guide eq:s2-pres).  Here phi is built as a product of
   * the unit character  u3(x) = the unit congruent to x mod 3          (a character of (O/3)^x),
   * residue characters  chi_P(x) = the root of unity = x^{(N P - 1)/d} mod P, d = gcd(6, N P - 1),
     for prime elements P not above 3 (split primary pi, inert -q, and P = 2 with d = 3).
Values are kept exactly as indices j mod 6 (zeta^j); additive characters e(hx/f) are kept exactly as indices mod Q.

zw.py is imported read-only from tests_cm/part_B_lines_0550_0997 (unchanged copy of tests_guide/zw.py).
"""
import os, sys, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(TESTS, "part_B_lines_0550_0997"))
import zw
from zw import mul, add, sub, conj, N, cx, divides, div, pw, powmod, primary, prim_assoc, UNITS, LAM, residues

ZETA6 = [cmath.exp(1j * math.pi * j / 3) for j in range(6)]
SQ3 = math.sqrt(3.0)


def unit_index(x):
    """j with x = zeta^j (x a unit)."""
    return UNITS.index(x)


def u3_index(x):
    """index j of the unit congruent to x mod 3; None if lambda | x."""
    if divides(LAM, x):
        return None
    for j, u in enumerate(UNITS):
        if divides((3, 0), sub(x, u)):
            return j
    raise ValueError(x)


def resid_index(x, P):
    """index j (zeta^j) of the residue symbol x^{(N P-1)/d} mod P, d = gcd(6, N P - 1); None if P | x."""
    if divides(P, x):
        return None
    n = N(P)
    d = math.gcd(6, n - 1)
    v = powmod(x, (n - 1) // d, P)
    for j in range(0, 6, 6 // d):
        if divides(P, sub(v, UNITS[j])):
            return j
    raise ValueError((x, P))


class Char:
    """Finite part phi of a CM character: phi(x) = zeta^{k*u3(x) + sum_P b_P chi_P(x)}; modulus f = given element."""

    def __init__(self, name, f, k3, comps):
        self.name = name
        self.f = f                      # modulus (element of O)
        self.Q = N(f)
        self.k3 = k3                    # exponent of the unit character mod 3 (0 if lambda does not divide f)
        self.comps = list(comps)        # [(P, b)]
        self.reps = residues(f)

    def idx(self, x):
        """phi(x) as an index mod 6, or None if x is not coprime to f."""
        if not coprime(x, self.f):
            return None
        j = 0
        if self.k3:
            j += self.k3 * u3_index(x)
        for P, b in self.comps:
            j += b * resid_index(x, P)
        return j % 6

    def phi(self, x):
        j = self.idx(x)
        return 0 if j is None else ZETA6[j]

    def psi_elem(self, x):
        """psi((x)) = phi(x) conj(alpha(x))  (eq:s2-pres)."""
        j = self.idx(x)
        if j is None:
            return 0
        z = cx(x)
        return ZETA6[j] * z.conjugate() / abs(z)

    def gauss_sum(self, h=(1, 0)):
        """sum_{x mod f} phi(x) e(hx/f) (complex)."""
        s = 0j
        fc = conj(self.f)
        for x in self.reps:
            j = self.idx(x)
            if j is not None:
                b = mul(mul(h, x), fc)[1] % self.Q
                s += ZETA6[j] * cmath.exp(2j * math.pi * b / self.Q)
        return s

    def gamma(self):
        return self.gauss_sum() / math.sqrt(self.Q)

    def root_number(self):
        """W(psi) = -gamma(phi) sqrt(Q)/f  (eq:d150)."""
        return -self.gamma() * math.sqrt(self.Q) / cx(self.f)


_prime_divs = {}


def prime_divisors(f):
    """Prime elements dividing f (lambda, 2, inert -q, split primary pi), from the rational norm."""
    if f in _prime_divs:
        return _prime_divs[f]
    from sympy import factorint
    out = []
    for p in factorint(N(f)):
        p = int(p)
        if p == 3:
            out.append(LAM)
        elif p % 3 == 2:
            if divides((p, 0), f):
                out.append((-p, 0) if p != 2 else (2, 0))
        else:
            r = int(math.isqrt(4 * p // 3 + 4)) + 2
            for a in range(-r, r + 1):
                for b in range(-r, r + 1):
                    if N((a, b)) == p and primary((a, b)) and divides((a, b), f):
                        out.append((a, b))
    _prime_divs[f] = out
    return out


def coprime(x, f):
    return all(not divides(P, x) for P in prime_divisors(f))


def is_primitive(ch):
    """phi primitive mod f: for each prime P | f, phi is nontrivial on {x = 1 mod f/P, x coprime to f}."""
    for P in prime_divisors(ch.f):
        g = div(ch.f, P)
        if not any(ch.idx(x) not in (None, 0) for x in ch.reps if divides(g, sub(x, (1, 0)))):
            return False
    return True


def split_primary(p):
    """the primary primes of norm p (p = 1 mod 3), sorted."""
    r = int(math.isqrt(4 * p // 3 + 4)) + 2
    return sorted((a, b) for a in range(-r, r + 1) for b in range(-r, r + 1) if N((a, b)) == p and primary((a, b)))


def test_chars():
    """Eight primitive CM characters (phi(zeta) = zeta), conductor exponents 0, 1, 2 at lambda, primes over 2, inert 5,
    split 7, 13.  A-D are the characters of tests_cm/phase2_2026-10-12/s2 (same phi, built independently here)."""
    p7 = split_primary(7)[0]
    p13 = split_primary(13)[0]
    i5 = (-5, 0)
    two = (2, 0)
    out = [
        Char("A: u3 (psi = alphabar), f = 3", (3, 0), 1, []),
        Char("B: chi_pi7, f = pi7", p7, 0, [(p7, 1)]),
        Char("C: u3^5 chi_pi13, f = 3 pi13", mul((3, 0), p13), 5, [(p13, 1)]),
        Char("D: u3^3 chi_pi7^4, f = lambda pi7", mul(LAM, p7), 3, [(p7, 4)]),
        Char("E: chi_2 chi_pi7^3, f = 2 pi7", mul(two, p7), 0, [(two, 1), (p7, 3)]),
        Char("F: u3^3 chi_5, f = lambda 5", mul(LAM, i5), 3, [(i5, 1)]),
        Char("G: u3 chi_5^3, f = 3*5", mul((3, 0), i5), 1, [(i5, 3)]),
        Char("H: chi_pi7 chi_pi13^3, f = pi7 pi13", mul(p7, p13), 0, [(p7, 1), (p13, 3)]),
    ]
    for ch in out:
        for j, u in enumerate(UNITS):
            assert ch.idx(u) == j, (ch.name, u, ch.idx(u))
        assert is_primitive(ch), ch.name
    return out


def lattice(nmax):
    """all nonzero (a, b) with N <= nmax."""
    out = []
    r = int(math.isqrt(4 * nmax // 3 + 4)) + 2
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a * a - a * b + b * b
            if 0 < n <= nmax:
                out.append((a, b))
    return out


# ------------------------------------------------------------------ exact cyclotomic arithmetic
def cyclo(M):
    """coefficients (low to high) of the M-th cyclotomic polynomial."""
    from sympy import cyclotomic_poly, Poly, symbols
    X = symbols("X")
    return [int(c) for c in reversed(Poly(cyclotomic_poly(M, X), X).all_coeffs())]


def reduce_cyclo(c, phiM):
    """remainder of sum c_k X^k modulo the monic polynomial phiM (exact integers)."""
    c = list(c)
    d = len(phiM) - 1
    for k in range(len(c) - 1, d - 1, -1):
        a = c[k]
        if a:
            base = k - d
            for i in range(d + 1):
                c[base + i] -= a * phiM[i]
    return c[:d]


# ------------------------------------------------------------------ reporting
RESULTS = []
_FORCE = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]


def report(cid, label, desc, ok, details="", falsifier_failed=None):
    """[PASS]/[FAIL] <id> <label> <desc> :: <details>.  With a falsifier, PASS needs ok and falsifier_failed."""
    if falsifier_failed is not None:
        details += " [falsifier failed as required]" if falsifier_failed else " [FALSIFIER DID NOT FAIL]"
        ok = ok and falsifier_failed
    if _FORCE[0]:
        _FORCE[0] = False
        ok = False
        details = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + details
    st = "PASS" if ok else "FAIL"
    RESULTS.append((st, cid))
    print(f"[{st}] {cid} {label} {desc} :: {details}", flush=True)
    return ok


def finish():
    nf = sum(1 for s, _ in RESULTS if s == "FAIL")
    print(f"SUMMARY: {len(RESULTS) - nf} PASS, {nf} FAIL")
    sys.stdout.flush()
    raise SystemExit(1 if nf else 0)


# ------------------------------------------------------------------ vectorized helpers (numpy)
OMEGA_C = complex(-0.5, SQ3 / 2)


def phi_table(ch):
    """phi on O/f as a Q x Q complex array indexed by the injective key (x conj f mod Q)."""
    import numpy as np
    Q = ch.Q
    T = np.zeros((Q, Q), dtype=complex)
    fc = conj(ch.f)
    for x in ch.reps:
        k = mul(x, fc)
        T[k[0] % Q, k[1] % Q] = ch.phi(x)
    return T


def keys_np(ch, A, B):
    c, d = conj(ch.f)
    return (A * c - B * d) % ch.Q, (A * d + B * c - B * d) % ch.Q


def ideal_coeffs(ch, B):
    """a(m) = sum_{N(I) = m} psi(I), m = 0..B, psi((a)) = phi(a) conj(alpha(a)); one generator per ideal
    (x + y zeta with x > 0, y >= 0)."""
    import numpy as np
    T = phi_table(ch)
    r = int(math.isqrt(B)) + 2
    xx, yy = np.meshgrid(np.arange(1, r + 1), np.arange(0, r + 1), indexing="ij")
    xx = xx.ravel(); yy = yy.ravel()
    n = xx * xx + xx * yy + yy * yy
    m = n <= B
    xx, yy, n = xx[m], yy[m], n[m]
    A = xx + yy; Bc = yy
    k0, k1 = keys_np(ch, A, Bc)
    z = A + Bc * OMEGA_C
    val = T[k0, k1] * np.conj(z) / np.abs(z)
    a = np.zeros(B + 1, dtype=complex)
    np.add.at(a, n, val)
    return a


def gp_run(script, timeout=3000, prec=38):
    """run a gp script (realprecision prec); returns the non-empty output lines; raises on a gp error."""
    import subprocess
    out = subprocess.run(["gp", "-q", "-D", f"realprecision={prec}", "-D", "parisizemax=4000000000"], input=script,
                         capture_output=True, text=True, timeout=timeout)
    errs = [l for l in (out.stdout + out.stderr).split("\n") if "***" in l and "Warning" not in l]
    if out.returncode != 0 or errs:
        raise RuntimeError(out.stdout + out.stderr)
    return [l for l in out.stdout.split("\n") if l.strip() and "***" not in l]


def gp_vec(a, start=1):
    """PARI vector literal of the complex numbers a[start:]."""
    return "[" + ",".join(f"{z.real:.17g}{z.imag:+.17g}*I" for z in a[start:]) + "]"


def gp_c(z):
    z = complex(z)
    return f"({z.real:.17g}{z.imag:+.17g}*I)"


def parse_gp_complex(sline):
    """parse a PARI complex/real printed value like '1.23 + 4.5*I' or '-0.1 - 2*I' or '3.2'."""
    t = sline.replace(" ", "").replace("E", "e")
    if "*I" not in t:
        return complex(float(t), 0)
    t = t.replace("*I", "")
    # split at the last + or - that is not part of an exponent and not leading
    for k in range(len(t) - 1, 0, -1):
        if t[k] in "+-" and t[k - 1] not in "eE":
            return complex(float(t[:k]), float(t[k:]))
    return complex(0, float(t))


def ideal_data(ch, B):
    """for each ideal I of norm <= B with psi(I) != 0: (N(I), j, A, Bc) with canonical generator g = A + Bc*omega
    (g = x + y zeta, x > 0, y >= 0) and phi(g) = zeta6^j, so that psi(I) = zeta6^j conj(g)/|g| exactly."""
    out = []
    r = int(math.isqrt(B)) + 2
    for x in range(1, r + 1):
        for y in range(0, r + 1):
            n = x * x + x * y + y * y
            if n > B:
                break
            g = (x + y, y)
            j = ch.idx(g)
            if j is not None:
                out.append((n, j, g[0], g[1]))
    return out


GP_HEAD = "w=(-1+sqrt(-3))/2;\nz6(j)=exp(Pi*I*j/3);\nal(z)=z/abs(z);"


def gp_char_setup(ch, B, name="a"):
    """gp code defining the coefficient vector <name> of L(s, psi) (exact data, working precision) and Wg = W(psi) =
    -gamma(phi) sqrt(Q)/f computed in gp from the residues mod f."""
    D = ideal_data(ch, B)
    fc = conj(ch.f)
    G = []
    for x in ch.reps:
        j = ch.idx(x)
        if j is not None:
            G.append((j, mul(x, fc)[1] % ch.Q))
    code = [GP_HEAD,
            "D=[" + ";".join(f"{n},{j},{A},{b}" for n, j, A, b in D) + "];",
            f"{name}=vector({B}); for(k=1,#D~, my(z=D[k,3]+D[k,4]*w); {name}[D[k,1]] += z6(D[k,2])*conj(z)/abs(z));",
            "GG=[" + ";".join(f"{j},{b}" for j, b in G) + "];",
            f"gam=sum(k=1,#GG~, z6(GG[k,1])*exp(2*Pi*I*GG[k,2]/{ch.Q}))/sqrt({ch.Q});",
            f"Wg=-gam*sqrt({ch.Q})/({ch.f[0]}+{ch.f[1]}*w);"]
    return "\n".join(code)
