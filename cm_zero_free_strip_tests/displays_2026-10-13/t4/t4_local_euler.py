"""t4 / L1: the local content of the row formula eq:d075 (guide (12.2), Proposition prop:cm-rows).

eq:d075 says that row u of the probe is zeta_K^S(6z) L^S(w, chi.(u)) / L^S(x, psi conj chi.(u)) * H_{psi,u}, with H the
Euler product of OAI78 Lemma lem:local-euler with eta(p) replaced by psi(p).  H is defined as a quotient
(OAI78 eq:local-H-definition), so the content of eq:d075 is the closed form of the local factor P_p
(OAI78 old-eq:5.13a and the table of J_j, L4010-4041) with eta(p) -> psi(p), a_p = abar(p)^3 psi(p)^3,
b_p = abar(p)^2 psi(p)^2 conj(chi_p(4)).

Test.  The local series P_p of OAI78 eq:local-P (L3998-4006) is summed by brute force (truncated; the sums converge
geometrically at the chosen point) with the target value at p replaced by
  (a) a generic unimodular number t = exp(0.7 i) (not a root of unity: nothing like t^6 = 1 can be used),
  (b) psi(p) = conj(alpha(p)) eta(p) for the CM target psi = abar * chi.(2+omega),
  (c) eta(p) = chi_p(2+omega) itself (OAI78's finite-order case, calibration),
and compared with the closed form, for all j = v_p(u) = 0..5 and all six values of rho = chi_p(u/p^j), at 8 split and
2 inert primes, at the point x = 1.2+0.4i, w = 1.1-0.3i, z = 0.45+0.2i.  The Gauss sums G_r = g_{chi_p^r}(p,1) are
computed by brute force (helpers of tests_cm/part_B_lines_0550_0997, read only).  OAI78's evaluations of
g_{chi_p^r}(p^t, p^j') (L3878-3884) and of C_p(t,k,j') (eq:local-scalar-values) used in the series are first
re-checked by brute force from their definitions at p of norm 7 and 13 (t <= 3, k <= 2).
Falsifiers (each must fail): b_p with chi_p(4) in place of conj(chi_p(4)); a_p = abar^3 eta^3 kept with the
finite-order eta while the c-variable carries psi (i.e. psi substituted in only one place); and the closed form
evaluated at conj(t) instead of t.
Exit status 1 on any failure.  Run time: about 10 s.
"""
import os, sys, math, cmath, itertools, time
sys.dont_write_bytecode = True
PB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "part_B_lines_0550_0997"))
sys.path.insert(0, PB)
from common import chi, gamma_j, N, cx, alpha, residues, e_add, mul, pw, divides, sub, div, primaries, mu  # noqa: E402

T0 = time.time()
NFAIL = [0]
def report(tid, ok, msg):
    if not ok:
        NFAIL[0] += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {tid} {msg}", flush=True)

def zpow(c, r):
    """chi^r with the zero convention (0 stays 0 for every r, also r = 0)."""
    return 0 if c == 0 else c**r

class Local:
    def __init__(self, p):
        self.p = p; self.Q = N(p)
        self.om = chi(p, (-1, 0)).real                         # omega_p = chi_p(-1) = +-1
        rq = math.sqrt(self.Q)
        self.G = {0: -1.0}
        for r in range(1, 6):
            self.G[r] = rq*gamma_j(r, p)                       # G_r = g_{chi_p^r}(p, 1)
        self.g1 = self.G[1]/rq; self.g3 = self.G[3]/rq          # gamma_1(p), gamma_3(p)
        self.abar = alpha(p).conjugate()
        self.chi4 = chi(p, (4, 0))

    def gpp(self, r, t, jp):
        """OAI78 L3878-3884: g_{chi_p^r}(p^t, p^j') for t >= 1."""
        Q = self.Q
        if r % 6:
            return Q**(t - 1)*self.G[r % 6] if jp == t - 1 else 0
        return (Q**t if jp >= t else 0) - (Q**(t - 1) if jp >= t - 1 else 0)

    def Cp(self, t, k, jp):
        """OAI78 eq:local-scalar-values."""
        if t == 0 and k == 0:
            return 1
        if t == 0:
            return 1 if jp == 0 else 0
        if k == 0:
            return self.gpp(t, t, jp)
        if k == 1 and jp >= 1:
            return self.G[1]/self.om*self.gpp(t - 1, t, jp - 1)
        return 0

    def brute(self, x, w, z, eta_c, a_p, j, rho, L=26, K=70, M=40):
        """OAI78 eq:local-P with eta(p)^{e0} -> eta_c^{e0} and the cube phase a_p given separately."""
        Q = self.Q; om = self.om; tot = 0
        for e0 in (0, 1):
            for l in range(L + 1):
                t = e0 + 3*l
                base = ((-1)**e0*self.g1**(-e0)*eta_c**e0*(a_p*self.g3.conjugate())**l
                        * Q**(-(x + 0.5)*e0 - (1 + 3*x)*l))
                for k in range(K + 1):
                    if t > 0 and k > 1:
                        break
                    for m in range(M + 1):
                        c = self.Cp(t, k, j + 6*m)
                        if c == 0:
                            continue
                        ex = t*k - e0*l - l*(l - 1)//2
                        tot += base*c*rho**(k - t)*om**(ex % 2)*Q**(-w*k - 6*z*m)
        return tot

    def closed(self, x, w, z, eta, a_p, b_p, j, rho):
        """OAI78 old-eq:5.13a with the table of J_j (L4010-4033)."""
        Q = self.Q
        V = Q**(-6*z); R = a_p**2*Q**(4 - 6*x - 6*z); Wl = rho*Q**(-w)
        J = {0: -eta/rho*Q**(-x) + Wl*R,
             1: eta*Q**(-x - w),
             2: a_p*rho**-3*Q**(1.5 - 3*x),
             3: -eta*b_p*rho**-2*Q**(2 - 3*x - w) + b_p**2*rho**-4*Q**(2 - 4*x),
             4: -eta*a_p*rho**-3*Q**(2.5 - 4*x - w),
             5: -a_p**2*Q**(3 - 6*x)}[j]
        Ps = (1/(1 - R))*((R*(1 - 1/Q) - eta*(Q - 1)*Q**(-x - w)*V**(1 if j <= 1 else 0))/(1 - V) + J)
        return 1/(1 - V) + (Wl/(1 - Wl) if j == 0 else 0) + Ps

# ------------------------------------------------------------------ L1a: brute-force check of OAI78's local evaluations
def g_brute(p, r, t, h):
    """g_{chi_p^r}(p^t, h) = sum_{d mod p^t} chi_p(d)^r e(h d/p^t), by brute force."""
    m = pw(p, t); mm = cx(m); s = 0
    for d in residues(m):
        c = zpow(chi(p, d), r)
        if c != 0:
            s += c*e_add(cx(mul(h, d))/mm)
    return s

def Cp_brute(p, t, k, jp):
    """eq:local-scalar-definition by brute force: sum over d mod p^k with p^k | p^j' - p^t d (modulus one: d = 0)."""
    pk = pw(p, k); s = 0
    for d in (residues(pk) if k > 0 else [(0, 0)]):
        h = sub(pw(p, jp), mul(pw(p, t), d))
        if k > 0 and not divides(pk, h):
            continue
        for _ in range(k):
            h = div(h, p)
        c1 = 1 if k == 0 else zpow(chi(p, d), k)
        c2 = 1 if t == 0 else g_brute(p, t, t, h)
        s += c1*c2
    return s

err_g = err_C = 0.0; ncheck = 0
SMALL = []
for q in primaries(14):
    if N(q) in (7, 13) and mu(q) == -1:
        SMALL.append(q)
SMALL = SMALL[:3]
for p in SMALL:
    Lp = Local(p)
    for t in range(1, 4):
        for r in range(0, 6):
            for jp in range(0, 5):
                if N(p)**t > 2200:
                    continue
                gb = g_brute(p, r, t, pw(p, jp)); gf = Lp.gpp(r, t, jp)
                err_g = max(err_g, abs(gb - gf)/N(p)**t); ncheck += 1
    for t in range(0, 4):
        for k in range(0, 3):
            for jp in range(0, 5):
                if N(p)**max(t, 1)*N(p)**k > 2200 or (t == 0 and k == 0 and jp > 0):
                    continue
                cb = Cp_brute(p, t, k, jp); cf = Lp.Cp(t, k, jp)
                err_C = max(err_C, abs(cb - cf)/N(p)**(t + 1)); ncheck += 1
report("L1a", err_g < 1e-9 and err_C < 1e-9,
       f"OAI78 L3878-3884 and eq:local-scalar-values re-derived by brute force at p = {SMALL} (norms "
       f"{[N(p) for p in SMALL]}), {ncheck} cases: max rel err g {err_g:.1e}, C_p {err_C:.1e}")

# ------------------------------------------------------------------ L1b: the closed form with eta(p) -> psi(p)
x, w, z = 1.2 + 0.4j, 1.1 - 0.3j, 0.45 + 0.2j
U0 = (2, 1)                                                  # eta = chi_.(2+omega), finite order
split = [q for q in primaries(80) if mu(q) == -1 and q[1] != 0][:8]
inert = [(-5, 0), (-11, 0)]
PRIMES = split + inert
ROOTS = [cmath.exp(1j*math.pi*k/3) for k in range(6)]
worst = {"generic": 0.0, "CM": 0.0, "finite": 0.0}
fal = {"chi4 unconjugated": 0.0, "a_p with eta, c with psi": 0.0, "closed form at conj(t)": 0.0}
nid = 0
for p in PRIMES:
    Lp = Local(p)
    etap = chi(p, U0)
    cases = {"generic": cmath.exp(0.7j), "CM": Lp.abar*etap, "finite": etap}
    for name, t in cases.items():
        a_p = Lp.abar**3*t**3
        b_p = Lp.abar**2*t**2*Lp.chi4.conjugate()
        for j in range(6):
            for rho in ROOTS:
                B = Lp.brute(x, w, z, t, a_p, j, rho)
                C = Lp.closed(x, w, z, t, a_p, b_p, j, rho)
                worst[name] = max(worst[name], abs(B - C)/max(1, abs(C))); nid += 1
                if name == "CM":
                    f1 = Lp.closed(x, w, z, t, a_p, Lp.abar**2*t**2*Lp.chi4, j, rho)
                    B2 = Lp.brute(x, w, z, t, Lp.abar**3*etap**3, j, rho)
                    tc = t.conjugate()
                    f3 = Lp.closed(x, w, z, tc, Lp.abar**3*tc**3, Lp.abar**2*tc**2*Lp.chi4.conjugate(), j, rho)
                    # a falsifier "fails" if it disagrees with the brute force somewhere (max over cases)
                    fal["chi4 unconjugated"] = max(fal["chi4 unconjugated"], abs(B - f1))
                    fal["a_p with eta, c with psi"] = max(fal["a_p with eta, c with psi"], abs(B2 - C))
                    fal["closed form at conj(t)"] = max(fal["closed form at conj(t)"], abs(B - f3))
ok = all(v < 1e-11 for v in worst.values())
falok = all(v > 1e-3 for v in fal.values())
report("L1b", ok and falok,
       f"eq:d075 local factor: brute-force OAI78 eq:local-P = closed form lem:local-euler with eta(p)->psi(p), "
       f"{len(split)} split + {len(inert)} inert primes, j=0..5, 6 values of rho, 3 target values ({nid} identities): "
       f"max rel err generic unit {worst['generic']:.1e}, CM psi {worst['CM']:.1e}, finite eta {worst['finite']:.1e}; "
       f"falsifiers (max deviation, must be large): " + ", ".join(f"{k}: {v:.2f}" for k, v in fal.items()))
print(f"# t4_local_euler: {time.time() - T0:.1f} s, {NFAIL[0]} failure(s)")
sys.exit(1 if NFAIL[0] else 0)
