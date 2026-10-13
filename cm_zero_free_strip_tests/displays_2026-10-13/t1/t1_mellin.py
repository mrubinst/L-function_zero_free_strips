"""Agent t1, new test T1-M: the two Mellin displays of Lemma lem:b1-mellin that had no direct test.

eq:d030 (6.4)  int_0^oo sum_m |q_d(m) W(m/D)| D^-sigma dD/D = sum_m |q_d(m)| m^-sigma int_0^oo |W(y)| y^sigma dy/y < oo
   The left side is computed as a genuine D-integral (Gauss-Legendre in log D) of the finite sum over m <= M; the right
   side with mpmath (30 digits).  Both sides are the truncation m <= M of the display, so they must agree to quadrature
   accuracy.  Finiteness: |q_d(m)| <= r(m) = #{n primary, coprime to 6, N(n) = m} for all m <= M (exact count), and
   the partial sums of sum r(m) m^-sigma stay below (1-2^-2sigma)(1-3^-sigma) zeta(sigma) L(sigma, chi_-3).
   Falsifier: y^-sigma in place of y^sigma on the right.
eq:d031 (6.5)  F_W(s) What'(s) = What(s) What'(s) B_d(s)/L(E_d, s+1/2) = F_W'(s) What(s)   (Re s > 1)
   F_W(s) = int_0^oo S_d(D) D^-s dD/D is computed directly from S_d(D) = sum_m q_d(m) W(m/D) (eq:target), by quadrature
   in log D over 1/b < D <= DMAX (S_d(D) = 0 for D < 1/b; the tail D > DMAX is below the stated tolerance), for two
   different complex weights; What by mpmath; L(E_d, s+1/2) and B_d(s) from PARI (lfun, minimal-model local factors).
   Falsifiers: the crossed pairing F_W What = F_W' What'; B_d omitted (d = 16, where B_d != 1).
"""
import math, cmath
import numpy as np
import mpmath as mp
from t1lib import check, gp_ap, qd_array, gl_panels
from cmcore import gp, columns, N

mp.mp.dps = 30


class Weight:
    """W(y) = exp(-1/((y-a)(b-y)) + 4) * phase(y) on (a, b), zero elsewhere."""
    def __init__(self, a, b, phase_np, phase_mp, name):
        self.a, self.b, self.ph, self.phm, self.name = a, b, phase_np, phase_mp, name

    def __call__(self, y):
        y = np.asarray(y, dtype=float)
        out = np.zeros(y.shape, dtype=complex)
        m = (y > self.a) & (y < self.b)
        t = y[m]
        out[m] = np.exp(-1.0/((t - self.a)*(self.b - t)) + 4)*self.ph(t)
        return out

    def mp(self, y):
        if not (self.a < y < self.b):
            return mp.mpf(0)
        return mp.e**(-1/((y - self.a)*(self.b - y)) + 4)*self.phm(y)

    def hat(self, s):
        """What(s) = int W(y) y^s dy/y."""
        return complex(mp.quad(lambda y: self.mp(y)*y**(s - 1), [self.a, (self.a + self.b)/2, self.b]))

    def abs_hat(self, sig):
        return float(mp.quad(lambda y: abs(self.mp(y))*y**(sig - 1), [self.a, (self.a + self.b)/2, self.b]))


W1 = Weight(1.0, 2.0, lambda t: np.exp(1j*t), lambda t: mp.e**(1j*t), 'bump(1,2) e^{iy}')
W2 = Weight(0.7, 1.9, lambda t: 1 + 2j*t, lambda t: 1 + 2j*t, 'bump(0.7,1.9) (1+2iy)')


def S_on_nodes(q, W, Ds):
    """S_d(D) = sum_m q(m) W(m/D) at each D (requires len(q) - 1 >= b*max(D))."""
    M = len(q) - 1
    out = np.zeros(len(Ds), dtype=complex)
    for k, D in enumerate(Ds):
        lo = max(1, int(math.floor(W.a*D))); hi = min(M, int(math.ceil(W.b*D)))
        if hi < lo:
            continue
        m = np.arange(lo, hi + 1)
        out[k] = np.sum(q[m]*W(m/D))
    return out


def abs_sum_on_nodes(q, W, Ds):
    M = len(q) - 1
    out = np.zeros(len(Ds))
    for k, D in enumerate(Ds):
        lo = max(1, int(math.floor(W.a*D))); hi = min(M, int(math.ceil(W.b*D)))
        if hi < lo:
            continue
        m = np.arange(lo, hi + 1)
        out[k] = np.sum(np.abs(q[m]*W(m/D)))
    return out


# ------------------------------------------------------------------ eq:d030
M = 3000
cols = columns(M)
r = np.zeros(M + 1, dtype=np.int64)
for n in cols:
    r[N(n)] += 1
worst = 0.0; fal = 0.0; rbad = 0; tailbad = 0; lines = []
for d in [1, 7, -3]:
    ap = gp_ap(d, M)
    q = qd_array(d, M, ap)
    rbad += int(np.sum(np.abs(q) > r + 1e-12))
    for sig in [1.5, 2.0, 3.0]:
        T, Wt = gl_panels(math.log(1/W1.b), math.log(M/W1.a), 0.01)
        Ds = np.exp(T)
        lhs = float(np.sum(Wt*abs_sum_on_nodes(q, W1, Ds)*np.exp(-sig*T)))
        ms = np.arange(1, M + 1)
        series = float(np.sum(np.abs(q[1:])*ms**(-sig)))
        rhs = series*W1.abs_hat(sig)
        wrong = series*W1.abs_hat(-sig)
        worst = max(worst, abs(lhs - rhs)/rhs)
        fal = max(fal, abs(lhs - wrong)/rhs)
        # finiteness: partial sums of r(m) m^-sigma below the Euler product over primary n coprime to 6
        bound = float((1 - mp.mpf(2)**(-2*sig))*(1 - mp.mpf(3)**(-sig))*mp.zeta(sig)
                      * mp.nsum(lambda k: mp.mpf(1)/(3*k + 1)**sig - mp.mpf(1)/(3*k + 2)**sig, [0, mp.inf]))
        part = float(np.sum(r[1:]*ms**(-sig)))
        tailbad += int(part > bound or series > part + 1e-12)
        lines.append(f"d={d} sigma={sig}: lhs {lhs:.12f} rhs {rhs:.12f}; sum|q|m^-s {series:.6f} <= sum r m^-s "
                     f"{part:.6f} <= {bound:.6f}")
for l in lines:
    print('   ', l)
check('T1-M1', 'eq:d030', '6.4', worst < 1e-9 and rbad == 0 and tailbad == 0,
      f"d in {{1,7,-3}}, sigma in {{1.5,2,3}}, m <= {M}, weight {W1.name}: D-integral of the sum = sum times y-integral, "
      f"max rel err {worst:.1e}; |q_d(m)| <= r(m) for all m <= {M} (violations {rbad}); partial sums bounded by "
      f"(1-4^-s)(1-3^-s)zeta_K(s) (violations {tailbad}); falsifier y^-sigma: rel err {fal:.2f}",
      falsifier=fal > 1e-3)

# ------------------------------------------------------------------ eq:d031
DMAX = 6000
SS = [3.0 + 0j, 2.2 + 5j]
worst = 0.0; fal_cross = 0.0; fal_B = 0.0; lines = []
for d in [1, 7, 16]:
    sl = ",".join(f"{s.real}+{s.imag}*I" for s in SS)
    out = gp(f"""\\p 30
E=ellinit([0,0,0,0,{d}]); Nc=ellglobalred(E)[1]; L=lfuncreate(E); SV=[{sl}];
{{for(k=1,#SV,s=SV[k]; v=lfun(L,s+1/2); B=1; forprime(p=2,max(3,abs({d})),if((6*{d})%p==0, a=ellap(E,p); f=if(Nc%p, 1-a*p^(-s-1/2)+p^(-2*s), 1-a*p^(-s-1/2)); B=B/f)); print(real(v)," ",imag(v)," ",real(B)," ",imag(B)))}}""")
    LB = [tuple(float(x.replace('E', 'e')) for x in l.split()) for l in out]
    Mq = int(math.ceil(2.0*DMAX)) + 2
    ap = gp_ap(d, Mq)
    q = qd_array(d, Mq, ap)
    for k, s in enumerate(SS):
        Lv = complex(LB[k][0], LB[k][1]); Bv = complex(LB[k][2], LB[k][3])
        F = {}
        for Wt_ in (W1, W2):
            T, Wq = gl_panels(math.log(1/Wt_.b), math.log(DMAX), 0.01)
            Ds = np.exp(T)
            F[Wt_.name] = complex(np.sum(Wq*S_on_nodes(q, Wt_, Ds)*np.exp(-s*T)))
        h1, h2 = W1.hat(s), W2.hat(s)
        lhs = F[W1.name]*h2
        mid = h1*h2*Bv/Lv
        rhs = F[W2.name]*h1
        tol = 1e-9 if s.real >= 3 else 2e-6
        scale = abs(mid)
        e = max(abs(lhs - mid), abs(rhs - mid), abs(lhs - rhs))/scale
        worst = max(worst, e/tol)
        fal_cross = max(fal_cross, abs(F[W1.name]*h1 - F[W2.name]*h2)/scale/tol)
        if d == 16:
            fal_B = max(fal_B, abs(lhs - h1*h2/Lv)/scale/tol)
        lines.append(f"d={d} s={s}: F_W1 W2^ = {lhs:.10f}, W1^W2^B/L = {mid:.10f}, F_W2 W1^ = {rhs:.10f}; "
                     f"rel err {e:.1e} (tol {tol:.0e})")
for l in lines:
    print('   ', l)
check('T1-M2', 'eq:d031', '6.5', worst <= 1.0,
      f"d in {{1,7,16}}, s in {{3, 2.2+5i}}, two complex weights, F_W from direct D-quadrature of S_d(D) "
      f"(D <= {DMAX}), L and B_d from PARI: max rel err/tol {worst:.1e}; falsifiers err/tol: crossed pairing "
      f"{fal_cross:.1e}, B_d omitted (d=16) {fal_B:.1e}",
      falsifier=fal_cross > 100 and fal_B > 100)
