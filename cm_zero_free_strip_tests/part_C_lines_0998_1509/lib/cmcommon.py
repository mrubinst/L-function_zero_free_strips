"""Shared arithmetic for the Part C audit.  Uses the copied helpers zw.py / common.py (conventions of the guide
and of OAI78: omega = e^{2 pi i/3}, primary = 1 mod 3, chi_n(x) = (x/n)_6 with zero values retained,
e(z) = exp(4 pi i Im z/sqrt 3), gamma_j(c) = N(c)^{-1/2} sum_{v mod c} chi_c(v)^j e(v/c), alpha(n) = n/|n|,
G(c) = conj(chi_c(4)) gamma_3(c))."""
import os, sys, math, cmath, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *            # noqa: F401,F403  (zw + chi, mu, gamma_j, G, primaries, PR_ODD, factor)
import zw

def abar(n):
    """conj(alpha(n)) = conj(n)/|n|."""
    return alpha(n).conjugate()

def pow0(c, k):
    """Power of a residue-symbol value with the zero convention: 0 stays 0 for every exponent (incl. 0)."""
    return 0 if c == 0 else c**k

def split_primes(B):
    """Primary primes of norm p = 1 mod 3, 5 <= p <= B (both conjugates)."""
    return [p for p in PR_ODD if N(p) <= B and p[1] != 0]

def inert_primes(B):
    """Primary generators -p of inert primes p = 2 mod 3, 5 <= p <= B (norm p^2)."""
    return [p for p in PR_ODD if p[1] == 0 and p[0] < 0 and -p[0] <= B and (-p[0]) % 3 == 2]

def gp(script, prec=38):
    """Run a PARI/GP script, return stdout as a list of lines."""
    out = subprocess.run(["gp", "-q", "-f", "-D", f"realprecision={prec}"], input=script,
                         capture_output=True, text=True, timeout=600)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    return [l for l in out.stdout.split("\n") if l.strip()]

def gpc(s):
    """Parse a PARI complex/real number printed by gp into a Python complex."""
    s = s.replace(" ", "").replace("E", "e")
    if "I" not in s:
        return complex(float(s))
    s = s.replace("*I", "j").replace("I", "1j")
    return complex(eval(s))
