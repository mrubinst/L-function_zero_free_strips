"""t4 / M: size checks for the two moment bounds of guide Section 14 (asymptotic bounds with unspecified constants;
proved in Appendix sec:s4-inverse, Prop. prop:s4-marked, and Appendix app:s5-plain, Prop. prop:s5-plain).
These are sanity numbers only: a finite computation cannot prove a bound "<< Z^{m+eps}".

M1  eq:cm-inverse (14.2) with the polynomials eq:d199 (14.1), WITH one prime slot (the existing check C103 of
    part_C t05 has no slot): sum_{0<N(u)<=Z^m} |M_u Q_u|^2 / Z^m for Z = 1e4, 4e4, 1.6e5, m = 1, r = 0.3, one slot
    z_1 = 0.25 (r+2z = 0.8 <= m - 0.2, 2r+8z = 2.6 <= 3m - 0.4, as required by the lemma), coefficients a_1(p) = 1,
    W = W_1 = smooth bump on [1,4]; epsilon = +1 and -1; nu = 1 and nu = chi_.(2) (a fixed ray character at 2, 3).
    M_u carries alpha(n)^eps, the slot carries no angular factor (the co-oriented CM case).  Reported: the ratios and
    the growth exponent log(ratio(4Z)/ratio(Z))/log 4 (a bound Z^{m+eps} allows only eps).
M2  eq:d084 (14.5) with eq:d083 (14.4), z = 0 (no slots), tau = 1 (psi_k = chi_.(k)), q = 0, M = m = 1, Z = 4000,
    n_1 = n_2 = 1/2 and n_1 = n_2 = 1 (no restriction on n_1, n_2 when z = 0): sum_k |S_k(n_1)S_k(n_2)|^2 / Z^M over all
    rows 0 < N(k) <= Z, for e = +1, -1.  Control e = 0 (no angular factor, OAI78's finite-order case): the principal
    rows k in K^{*6} (here k = 1, -27 = lambda^6) then carry volume-sized main terms and the full sum is ~ Z^{n_1+n_2} >> Z^M
    for n_1 = n_2 = 1, which is why OAI78 restricts to nonprincipal rows and needs the centred differences; with the
    angular factor those rows are small (the zero-frequency integral of alpha^e vanishes, guide item (c)).
Pass criteria (loose, for the size check only): M1 growth exponents <= 0.25 and ratios <= 10; M2 ratios <= 10 for
e = +-1 and the e = 0 control with principal rows >= 20 times larger for n_1 = n_2 = 1.
Exit status 1 on failure.  Run time: about 1-2 minutes.
"""
import os, sys, math, time
sys.dont_write_bytecode = True
PB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "part_B_lines_0550_0997"))
sys.path.insert(0, PB)
import numpy as np
from cmlib import chi_vec, primaries, mu, N, cx, factor, bump      # noqa: E402  (read-only helpers)

T0 = time.time()
NFAIL = [0]
def report(tid, ok, msg, status=None):
    st = status or ("PASS" if ok else "FAIL")
    if st == "FAIL":
        NFAIL[0] += 1
    print(f"[{st}] {tid} {msg}", flush=True)

def W(x):
    return bump(x, 1.0, 4.0)

def rows_upto(X):
    R = int(math.sqrt(4*X/3)) + 2
    a = np.arange(-R, R + 1)
    A, B = np.meshgrid(a, a, indexing="ij"); A = A.ravel(); B = B.ravel()
    Nk = A*A - A*B + B*B
    keep = (Nk > 0) & (Nk <= X)
    return A[keep].astype(np.int64), B[keep].astype(np.int64)

def alpha(n):
    z = complex(cx(n)); return z/abs(z)

# ------------------------------------------------------------------ M1: inverse moment with one slot
m, r, z1 = 1.0, 0.3, 0.25
res = {}
for Z in (1.0e4, 4.0e4, 1.6e5):
    A, B = rows_upto(Z**m)
    cols = [n for n in primaries(int(4*Z**r) + 1, sqfree=True) if W(N(n)/Z**r) > 0]
    slot = [p for p in primaries(int(4*Z**z1) + 1) if len(factor(p)) == 1 and W(N(p)/Z**z1) > 0]
    vec = {}
    for n in set(cols) | set(slot):
        vec[n] = chi_vec(n, A, B)                     # chi_n(u) over all rows (zero-extended)
    nu = {"nu=1": lambda n: 1.0, "nu=chi(2)": lambda n: complex(chi_vec(n, np.array([2]), np.array([0]))[0])}
    for nuname, nuf in nu.items():
        for eps in (+1, -1):
            Mu = np.zeros(len(A), complex); Qu = np.zeros(len(A), complex)
            for n in cols:
                v = vec[n] if eps == 1 else np.conj(vec[n])
                Mu += mu(n)*alpha(n)**eps*nuf(n)*v*W(N(n)/Z**r)
            for p in slot:
                v = vec[p] if eps == 1 else np.conj(vec[p])
                Qu += nuf(p)*v*W(N(p)/Z**z1)
            Mu *= Z**(-r/2); Qu *= Z**(-z1/2)
            res[(Z, nuname, eps)] = (np.sum(np.abs(Mu*Qu)**2)/Z**m, len(A), len(cols), len(slot))
ok = True; parts = []; growth = []
for nuname in ("nu=1", "nu=chi(2)"):
    for eps in (+1, -1):
        rs = [res[(Z, nuname, eps)][0] for Z in (1.0e4, 4.0e4, 1.6e5)]
        g = max(math.log(rs[1]/rs[0]), math.log(rs[2]/rs[1]))/math.log(4)
        growth.append(g); ok &= g <= 0.25 and max(rs) <= 10
        parts.append(f"{nuname}, eps={eps:+d}: " + ", ".join(f"{x:.3f}" for x in rs) + f" (growth exp {g:+.2f})")
sizes = "; ".join(f"Z={Z:.0e}: {res[(Z, 'nu=1', 1)][1]} rows, {res[(Z, 'nu=1', 1)][2]} columns, "
                  f"{res[(Z, 'nu=1', 1)][3]} slot primes" for Z in (1.0e4, 4.0e4, 1.6e5))
report("M1", ok, "eq:cm-inverse size check with one slot (m=1, r=0.3, z_1=0.25): sum|M_uQ_u|^2/Z^m at Z = 1e4, 4e4, "
       "1.6e5: " + "; ".join(parts) + f" [{sizes}]")

# ------------------------------------------------------------------ M2: plain fourth moment, z = 0
Z = 4000.0; M = 1.0
A, B = rows_upto(Z)
princ = np.zeros(len(A), bool)
for k in [(1, 0), (-27, 0)]:                          # the sixth powers a^6 with N(a^6) <= Z (a = unit, lambda)
    princ |= (A == k[0]) & (B == k[1])
assert princ.sum() == 2
out = {}
for nn in (0.5, 1.0):
    ls = [l for l in primaries(int(4*Z**nn) + 1) if W(N(l)/Z**nn) > 0]
    S = {e: np.zeros(len(A), complex) for e in (1, -1, 0)}
    for l in ls:
        v = chi_vec(l, A, B); wl = W(N(l)/Z**nn); al = alpha(l)
        for e in S:
            S[e] += al**e*v*wl
    for e in S:
        S[e] *= Z**(-nn/2)
        tot = np.abs(S[e]*S[e])**2                     # n_1 = n_2 = nn, the same test W
        out[(nn, e)] = (tot.sum()/Z**M, tot[~princ].sum()/Z**M, len(ls))
ok = all(out[(nn, e)][0] <= 10 for nn in (0.5, 1.0) for e in (1, -1)) and out[(1.0, 0)][0] >= 20*max(out[(1.0, 1)][0], 1e-9)
msg = "; ".join(f"n1=n2={nn}: e=+1 {out[(nn, 1)][0]:.3f}, e=-1 {out[(nn, -1)][0]:.3f}, control e=0 {out[(nn, 0)][0]:.3g} "
                f"(without the 2 principal rows {out[(nn, 0)][1]:.3f}) [{out[(nn, 1)][2]} plain l]" for nn in (0.5, 1.0))
report("M2", ok, f"eq:d084 size check (z=0, tau=1, M=m=1, Z=4000, {len(A)} rows): sum_k|S_k(n1)S_k(n2)|^2/Z^M: " + msg)
print(f"# t4_moments: {time.time() - T0:.1f} s, {NFAIL[0]} failure(s)")
sys.exit(1 if NFAIL[0] else 0)
