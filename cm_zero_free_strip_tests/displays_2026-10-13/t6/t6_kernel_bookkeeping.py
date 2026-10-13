"""t6 (2026-10-13): the analytic and bookkeeping displays of Appendix app:s1-general that are not end-to-end identities.

 K1 eq:s1-mellin-tail  the identity behind it, (-t)^N Vhat(t) = int x^t (x d/dx)^N V(x) dx/x (N = 1,2,3), and the
                       explicit bound (1+|u|)^N |Vhat(eta+iu)| <= 2^N sum_{j<=N} int x^eta |(x d/dx)^j V| dx/x on
                       eta in [-1, 1/4], |u| <= 30, for V_G (Vhat = e^{t^2}) and exp(-y-1/y) (Vhat = 2K_t(2)).
                       Falsifier: the bound without the j = N term fails for large |u|.
 K2 eq:d123            M_{A,B}(V_G) = e^{A^2} int (1+|u|)^B e^{-u^2} du for A >= 1/4 (sup at eta = -A), by a grid in eta.
 K3 eq:d128            |(x d/dx)^j V#(x)| / (min{x^{1/4},(1+x)^{-A}} M_{A,ceil(4A)+j+2}(V)) on 1e-8 <= x <= 1e6, j <= 2,
                       A = 1, 2, 3, for both tests (size check; the unnormalised ratio is t_s1_kernel.py (4c)).
 K4 eq:s1-stirling     sup_{|u| <= 1e4} |R_2(eta+iu)| (1+|u|)^{-4 eta} finite on eta in [-1.3, 6]; it grows as eta -> -4/3.
                       Falsifier: exponent 4 eta - 1/2 (ratio grows like |u|^{1/2}).
 K5 eq:d134            sup_y y|V(y)| = C_V < oo and sum_{n sf, b} |V(N(n)N(b)^3/X)|/(sqrt(N n) N b) <= C_V X sum N(n)^{-3/2}
                       sum N(b)^{-4} (truncated sums, X = 10, 1000), and the exponents (sympy).
 K6 bookkeeping        (sympy) eq:d132 exponents; eq:d133 Gamma shifts and the growth exponent 2 - 4 sigma_2 of the gamma
                       quotient; eq:s1-fe substitution v = 1/(N(c)v') giving conj(alpha(c))^4 N(c)^{2-2s} J^vee(2-s);
                       eq:s1-cm-completed: O/2 + max{M-O, 2M-O-1} <= max{M, 2M-1} for O >= 0 and Z^{max{M,2M-1}} =
                       max(Q, Q^2/Z); eq:d146: ((d-1/6)/4)_+ = 0 for 0 <= d <= 1/6 (OAI78 old-eq:5.2).
 K7 eq:d140, eq:d141   exact identities in Z[omega] (conj(alpha(x))^4 = conj(x)^4/N(x)^2, alpha(x)^2 = x/conj(x)), random
                       data; falsifiers conj(alpha(c))^2 on one side, alpha(b)^2 instead of alpha(b)^6.
 K9 eq:d122            phi from sextic symbols for three Psi_0 is 18-periodic and has the closed forms used in the tests.
 K8 eq:s1-scalar, eq:d139  |C| <= 1/81 and |C_Q| <= 1/81 over all branches (h0 for three twists modulo 18, 15 good
                       primes, all units c/p, j = 0..5), and C = -i conj(alpha(c))^2 zeta with zeta from OAI78
                       eq:reflection-row-phase (C_Q likewise with eq:marked-reflection-phase).
Usage: python3 t6_kernel_bookkeeping.py     (a few minutes; the Bessel seminorms dominate)
"""
import os, sys, time, math, cmath, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
from t6_lib import S1
from zw import mul, add, sub, conj, N, cx, divides, div, residues, sym6_prime, pw, UNITS, primary
from fastgauss import primary_primes
import numpy as np, mpmath as mp, sympy as sp
t0 = time.time()
CHECKS = []
def check(name, ok, msg):
    CHECKS.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + ": " + msg, flush=True)

# ------------------------------------------------------------------ tests V in the log variable x = e^w
w = sp.symbols('w', real=True)
VG_w = sp.exp(-w**2/4)/(2*sp.sqrt(sp.pi))           # V_G(e^w)
VK_w = sp.exp(-sp.exp(w)-sp.exp(-w))                # exp(-y-1/y) at y = e^w
DER = {}
for name, f in (('G', VG_w), ('K', VK_w)):
    DER[name] = [sp.lambdify(w, sp.diff(f, w, j), 'mpmath') for j in range(4)]   # (x d/dx)^j V = d^j/dw^j
VHAT = {'G': lambda t: mp.e**(t*t), 'K': lambda t: 2*mp.besselk(t, 2)}
WR = {'G': [-40, -10, 0, 10, 40], 'K': [-7, -2, 0, 2, 7]}
mp.mp.dps = 30
# K1
werr = 0; worst_ratio = 0
for name in ('G', 'K'):
    for t in ([mp.mpc(0.25, 1.3), mp.mpc(-0.6, 2.0), mp.mpc(-1, 3.0)] if name == 'G' else [mp.mpc(0.25, 1.3), mp.mpc(-0.6, 4.0), mp.mpc(-1, 9.5)]):
        for Nn in (1, 2, 3):
            lhs = (-t)**Nn*VHAT[name](t)
            rhs = mp.quad(lambda u: DER[name][Nn](u)*mp.e**(t*u), WR[name])
            werr = max(werr, abs(lhs-rhs)/abs(lhs))
    for eta in (-1, -0.5, 0, 0.25):
        I = [mp.quad(lambda u: abs(DER[name][j](u))*mp.e**(eta*u), WR[name]) for j in range(4)]
        for u in (0, 1, 3, 10, 30):
            vh = abs(VHAT[name](mp.mpc(eta, u)))
            for Nn in (1, 2, 3):
                worst_ratio = max(worst_ratio, (1+u)**Nn*vh/(2**Nn*sum(I[:Nn+1])))
check("K1 eq:s1-mellin-tail", werr < 1e-12 and worst_ratio <= 1,
      f"identity (-t)^N Vhat = Mellin of (x d/dx)^N V: max rel err {mp.nstr(werr, 2)}; explicit bound with constant 2^N: "
      f"max ratio {mp.nstr(worst_ratio, 3)} (<= 1) over eta in -1,-1/2,0,1/4, u in 0..30, N=1,2,3")
# K2
def MG(A, B):
    I = mp.quad(lambda u: 2*(1+u)**B*mp.e**(-u*u), [0, 5, 40])
    return max(mp.e**(eta*eta) for eta in np.linspace(-A, 0.25, 41))*I, mp.e**(max(A, 0.25)**2)*I
err2 = max(abs(a-b)/b for a, b in (MG(A, B) for A in (0.25, 1, 2, 3) for B in (0, 6, 16)))
check("K2 eq:d123 for V_G", err2 < 1e-15, f"sup over eta grid = e^(A^2) int (1+|u|)^B e^(-u^2) du for A in 1/4,1,2,3, B in 0,6,16: {mp.nstr(err2, 2)}")
print(f"[{time.time()-t0:.1f} s]", flush=True)
# K3: V# by the trapezoid rule on Re t = s0 (as in t_s1_kernel.py), normalised by M
K = float((2*mp.pi)**4/27)
R2 = lambda t: mp.gamma(mp.mpf(4)/3+t)*mp.gamma(mp.mpf(5)/3+t)/(mp.gamma(mp.mpf(4)/3-t)*mp.gamma(mp.mpf(5)/3-t))
H = 0.05; U = np.arange(-80, 80+1e-9, H)
TAB = {}
def table(kind, s0):
    if (kind, s0) not in TAB:
        TAB[(kind, s0)] = np.array([complex(VHAT[kind](-mp.mpc(s0, u))*R2(mp.mpc(s0, u))) for u in U])
    return TAB[(kind, s0)]
def Vs(x, kind, s0, j):
    t = s0+1j*U
    return complex(np.sum((-t)**j*table(kind, s0)*np.exp(-t*np.log(K*x)))*H/(2*np.pi))
def M_K(A, B):     # M_{A,B} for exp(-y-1/y): sup over an eta grid of int (1+|u|)^B |2K_{eta+iu}(2)| du (symmetric in u)
    best = 0
    uu = np.arange(0, 90+1e-9, 0.25)
    for eta in np.linspace(-A, 0.25, 5):
        vals = np.array([float(abs(2*mp.besselk(mp.mpc(eta, u), 2))) for u in uu])
        f = (1+uu)**B*vals; I = 2*(np.sum(f)-f[0]/2-f[-1]/2)*0.25
        best = max(best, I)
    return best
res = {}
for kind in ('G', 'K'):
    for A in (1, 2, 3):
        mx = 0
        for j in (0, 1, 2):
            Bo = math.ceil(4*A)+j+2
            M = float(MG(A, Bo)[1]) if kind == 'G' else M_K(A, Bo)
            for e in range(-8, 7):
                x = 10.0**e
                v = Vs(x, kind, 6 if e > 1 else 0, j)
                mx = max(mx, abs(v)/(min(x**0.25, (1+x)**(-A))*M))
        res[(kind, A)] = mx
    print(f"   K3 {kind} done [{time.time()-t0:.1f} s]", flush=True)
check("K3 eq:d128 (size check)", all(np.isfinite(v) and v < 10 for v in res.values()),
      "max over x in [1e-8,1e6], j<=2 of |(x d/dx)^j V#|/(min{x^1/4,(1+x)^-A} M_{A,ceil(4A)+j+2}(V)): "
      + ", ".join(f"{k[0]} A={k[1]}: {v:.2e}" for k, v in res.items()))
# K4
uu = np.concatenate([np.linspace(0, 10, 101), np.logspace(1, 4, 61)])
rat = {}; growth = {}
for eta in (-1.3, -1.2, -1.0, -0.5, -0.25, 0, 0.5, 2, 6):
    r = [float(abs(R2(mp.mpc(eta, u))))*(1+u)**(-4*eta) for u in uu]
    rat[eta] = max(r)
    growth[eta] = r[-1]/r[-11]                               # u = 1e4 vs u = 10^3.5
fal4 = min(float(abs(R2(mp.mpc(eta, 1e4))))*(1+1e4)**(-4*eta+0.5)/(float(abs(R2(mp.mpc(eta, 1e2))))*(1+1e2)**(-4*eta+0.5))
           for eta in (-1, 0, 2))
check("K4 eq:s1-stirling (size check)", all(np.isfinite(v) for v in rat.values()) and max(abs(g-1) for g in growth.values()) < 0.01 and fal4 > 5,
      "sup_{|u|<=1e4} |R_2(eta+iu)|(1+|u|)^{-4eta}: " + ", ".join(f"{e}: {v:.3g}" for e, v in rat.items())
      + f"; ratio at u=1e4 vs 3162: within {max(abs(g-1) for g in growth.values()):.1e} of 1; falsifier exponent 4eta-1/2 grows by {fal4:.1f} from u=1e2 to 1e4")
# K5 eq:d134
yy = np.logspace(-6, 6, 20001)
CV = {'G': float(np.max(yy*np.exp(-np.log(yy)**2/4)/(2*math.sqrt(math.pi)))), 'K': float(np.max(yy*np.exp(-yy-1/yy)))}
PR = primary_primes(20000)
sq = [(1, 0)]
for p in PR:
    if p == (-2, 0): continue
    for c in list(sq):
        cc = mul(c, p)
        if N(cc) <= 20000: sq.append(cc)
nN = np.array([N(c) for c in sq], float)
bN = np.array([N((x, y)) for x in range(-30, 31) for y in range(-30, 31) if primary((x, y)) and (x % 2 or y % 2) and N((x, y)) <= 800], float)
ok5 = True; msg5 = []
for name, V in (('G', lambda y: np.exp(-np.log(y)**2/4)/(2*math.sqrt(math.pi))), ('K', lambda y: np.exp(-y-1/y))):
    for X in (10.0, 1000.0):
        lhs = sum(np.sum(np.abs(V(nN*b**3/X))/(np.sqrt(nN)*b)) for b in bN)
        rhs = CV[name]*X*np.sum(nN**-1.5)*np.sum(bN**-4.0)
        ok5 &= lhs <= rhs; msg5.append(f"{name} X={X:g}: {lhs:.4g} <= {rhs:.4g}")
y_, Nn_, Nb_, X_ = sp.symbols('y N_n N_b X', positive=True)
expo = sp.simplify((X_/(Nn_*Nb_**3))/(sp.sqrt(Nn_)*Nb_)-X_*Nn_**sp.Rational(-3, 2)*Nb_**-4)
check("K5 eq:d134", ok5 and expo == 0, f"C_V = sup y|V(y)|: {CV['G']:.4f} (V_G), {CV['K']:.4f}; " + "; ".join(msg5) + f"; exponent identity (sympy): {expo == 0}")
# K6 bookkeeping
k, sg = sp.symbols('k sigma', real=True); s, wv, t = sp.symbols('s w t')
Nn, Nb = sp.symbols('N_n N_b', positive=True)
e1 = sp.simplify(sp.expand_log(sp.log(3**(k/6)*sp.sqrt(Nb)*(3**k*Nn*Nb**3)**(sp.Rational(1, 2)-sg)), force=True)
                 -sp.expand_log(sp.log(3**(k*(sp.Rational(2, 3)-sg))*Nn**(sp.Rational(1, 2)-sg)*Nb**(2-3*sg)), force=True))
e2 = [sp.simplify((wv+sp.Rational(1, 2))+sp.Rational(1, 3)-(wv+sp.Rational(5, 6))), sp.simplify((wv+sp.Rational(1, 2))+sp.Rational(2, 3)-(wv+sp.Rational(7, 6)))]
# growth exponent of Gamma(11/6-w)Gamma(13/6-w)/(Gamma(w+5/6)Gamma(w+7/6)) on Re w = sigma2: sum (Re arg - 1/2)
s2 = sp.symbols('sigma_2')
e3 = sp.simplify((sp.Rational(11, 6)-s2-sp.Rational(1, 2))+(sp.Rational(13, 6)-s2-sp.Rational(1, 2))-(s2+sp.Rational(5, 6)-sp.Rational(1, 2))-(s2+sp.Rational(7, 6)-sp.Rational(1, 2))-(2-4*s2))
# substitution: int (c v)^-4 G(1/(N v)) v^(2s-1) dv, v = 1/(N v'): numeric with a test G
mp.mp.dps = 30
c = mp.mpc(2.3, -1.4); Nc = abs(c)**2
G = lambda vp: vp**2*mp.e**(-vp)*(1+mp.sin(3*vp))
e4 = 0
for sv in (mp.mpf(-1), mp.mpc(-0.5, 2)):
    lhs = mp.quad(lambda v: (c*v)**-4*G(1/(Nc*v))*v**(2*sv-1), [0, 0.01, 0.1, 1, mp.inf])
    rhs = mp.conj(c/abs(c))**4*Nc**(2-2*sv)*mp.quad(lambda vp: G(vp)*vp**(2*(2-sv)-1), [0, 1, 10, mp.inf])
    e4 = max(e4, abs(lhs-rhs)/abs(rhs))
# eq:s1-cm-completed
Mg = np.linspace(0, 3, 61); Og = np.linspace(0, 6, 121)
e5 = max(o/2+max(m-o, 2*m-o-1)-max(m, 2*m-1) for m in Mg for o in Og)
Zv = 1e6; e6 = max(abs(math.log(max(Zv**m, Zv**(2*m)/Zv))/math.log(Zv)-max(m, 2*m-1)) for m in Mg)
e5f = max(1.5*o+max(m-o, 2*m-o-1)-max(m, 2*m-1) for m in Mg for o in Og)   # falsifier: 3O/2 instead of O/2
# eq:d146
dg = np.linspace(0, 1/6, 101); e7 = max(max((d-1/6)/4, 0) for d in dg); e7f = max((1/6+0.01-1/6)/4, 0)
check("K6 bookkeeping", e1 == 0 and all(x == 0 for x in e2) and e3 == 0 and e4 < 1e-15 and e5 <= 1e-12 and e6 < 1e-9 and e5f > 0.5 and e7 == 0 and e7f > 0,
      f"eq:d132 exponents (log difference) {e1}; eq:d133 Gamma shifts {e2}; growth exponent 2-4sigma_2: {e3 == 0}; eq:s1-fe substitution "
      f"(test function): {mp.nstr(e4, 2)}; eq:s1-cm-completed: max of O/2+max{{M-O,2M-O-1}}-max{{M,2M-1}} = {e5:.2g} (O in [0,6], M in [0,3]), "
      f"Q+Q^2/Z form {e6:.1e}, falsifier 3O/2 instead of O/2 exceeds by {e5f:.2f}; eq:d146: max ((d-1/6)/4)_+ on [0,1/6] = {e7}")
# K7 exact identities in Z[omega]
random.seed(7)
def rnd(): return (random.randint(-30, 30), random.randint(-30, 30))
bad = 0; badf = 0
for _ in range(300):
    cF, F, R, P = [rnd() for _ in range(4)]
    if 0 in (N(cF), N(F), N(R), N(P)): continue
    cc = mul(mul(cF, F), mul(R, P))
    lhs = mul(pw(conj(cc), 4), (N(mul(cF, F))**2*N(R)**2*N(P)**2, 0))
    rhs = mul(mul(mul(pw(conj(mul(cF, F)), 4), pw(conj(R), 4)), pw(conj(P), 4)), (N(cc)**2, 0))
    bad += lhs != rhs
    # falsifier: conj(alpha(c))^2 on the left
    badf += mul(pw(conj(cc), 2), (N(mul(cF, F))**2*N(R)**2*N(P)**2*N(cc), 0)) != mul(rhs, (1, 0))
    # eq:d141 with l' = lambda^4 l = x' n b^3, x' = u lambda^{h+4} N_F B_F^3:  alpha(l)^2 = alpha(x')^2 alpha(n)^2 alpha(b)^6
    u = random.choice(UNITS); h = random.randint(-4, 3); NF, BF, n, b = rnd(), rnd(), rnd(), rnd()
    if 0 in (N(NF), N(BF), N(n), N(b)): continue
    x = mul(mul(u, pw((1, 2), h+4)), mul(NF, pw(BF, 3))); l = mul(mul(x, n), pw(b, 3))
    # alpha(y)^2 = y / conj(y): compare l * conj(x) conj(n) conj(b)^3 with conj(l) * x n b^3
    bad += mul(l, mul(mul(conj(x), conj(n)), pw(conj(b), 3))) != mul(conj(l), mul(mul(x, n), pw(b, 3)))
lam = complex(cx((1, 2))); e8 = abs((lam/abs(lam))**8-1)
# numeric falsifier for d141: alpha(b)^2 instead of alpha(b)^6
fz = 0
for _ in range(50):
    n, b, x = [complex(cx(rnd())) or 1 for _ in range(3)]
    a2 = lambda z: (z/abs(z))**2
    fz = max(fz, abs(a2(x*n*b**3)-a2(x)*a2(n)*a2(b)))
check("K7 eq:d140, eq:d141 (exact in Z[omega])", bad == 0 and badf > 250 and e8 < 1e-14 and fz > 0.5,
      f"300 random tuples: {bad} failures; alpha(lambda)^8 = 1: {e8:.1e}; falsifiers: conj(alpha(c))^2 fails in {badf} cases, alpha(b)^2 deviates by {fz:.2f}")
print(f"[{time.time()-t0:.1f} s]", flush=True)
# K8 scalars
LAM = (1, 2)
def e_add(z): return cmath.exp(4j*math.pi*z.imag/math.sqrt(3))
RES18 = residues((18, 0))
def lam3(x): return cmath.exp(2j*math.pi/3)**((-(x[1]//3)) % 3)
def om3(x): return cmath.exp(2j*math.pi/3)**(((N(x)-1)//3) % 3)
PHI = {'1': lambda x: lam3(x), 'A': lambda x: 1, 'B': lambda x: om3(x)}
phimax = 0
for tw, f in PHI.items():
    vals = [(complex(cx(x)), f(x)) for x in RES18 if x[0] % 3 == 1 and x[1] % 3 == 0 and (x[0] % 2 or x[1] % 2)]
    for h in RES18:
        ph = sum(v*e_add(-complex(cx(h))*xc/18) for xc, v in vals)/324
        phimax = max(phimax, abs(ph))
PRS = [p for p in primary_primes(200) if p != (-2, 0) and N(p) <= 61]
locmin, locmax = 9, 0; markmax = 0
for p in PRS:
    q = N(p); RES = residues(p); P = complex(cx(p))
    def chi(x, m):
        kk = sym6_prime(x, p)
        return 0 if kk is None else cmath.exp(1j*math.pi*kk*(m % 6)/3)
    def inv(x):
        for y in RES:
            if divides(p, sub(mul(x, y), (1, 0))): return y
    def gpm(m, s_): return sum(chi(y, m)*e_add(s_*complex(cx(y))/P) for y in RES)/math.sqrt(q)
    for u in [y for y in RES if not divides(p, y)][:12]:
        sig = mul(mul(LAM, LAM), u); eps = mul((-1, 0), inv(mul(mul(pw(LAM, 3), u), sig)))
        for j in range(6):                                   # OAI78 L1832-1840
            if j == 4: varpi = gpm(4, -1)
            elif j == 0: varpi = -gpm(2, 1)*chi(eps, -2)
            else: varpi = gpm(j, -1)*gpm(j+2, 1)*chi(eps, -j-2)
            a = abs(chi(sig, -2)*varpi); locmin = min(locmin, a); locmax = max(locmax, a)
            if j == 0: markmax = max(markmax, abs(chi(sig, -2)*(-varpi)))
# C and C_Q from their factors (eq:s1-scalar, eq:d139) vs OAI78's zeta, zeta_Q: random data
rel = 0; bigC = 0
for _ in range(200):
    c = complex(random.uniform(-9, 9), random.uniform(-9, 9)); ac = np.conj(c/abs(c))
    ph = phimax*cmath.exp(2j*math.pi*random.random()); kF = cmath.exp(2j*math.pi*random.random())
    rest = np.prod([cmath.exp(2j*math.pi*random.random()) for _ in range(random.randint(0, 4))])*np.prod([1-1/random.choice([7, 13, 19, 25]) for _ in range(random.randint(0, 3))])
    C = -(ac**4)/81*ph*np.conj(kF)*rest
    zeta = -(1j/81)*ac**2*ph*np.conj(kF)*rest
    rel = max(rel, abs(C-(-1j)*ac**2*zeta)); bigC = max(bigC, 81*abs(C))
check("K8 eq:s1-scalar, eq:d139 (|C|<=1/81, C=-i conj(alpha(c))^2 zeta)",
      phimax <= 1 and abs(locmin-1) < 1e-12 and abs(locmax-1) < 1e-12 and abs(markmax-1) < 1e-12 and bigC <= 1+1e-12 and rel < 1e-15,
      f"max |phihat(h0)| over 3 twists mod 18 = {phimax:.4f}; |chi_p(sigma_p)^-2 varpi_p,j| in [{locmin:.15f}, {locmax:.15f}] and "
      f"|chi_q(sigma_q)^-2(-varpi_q,0)| <= {markmax:.15f} (15 primes, 12 units c/p, j=0..5), so |C|, |C_Q| <= 1/81; "
      f"|C - (-i conj(alpha(c))^2 zeta_OAI78)| <= {rel:.1e} (200 random data, 81|C| <= {bigC:.4f})")
# K9 eq:d122: phi(x) = chi_x(lambda)^2 Psi_0(x) (zero unless x = 1 mod 3, (x,6) = 1) computed from sextic symbols, for
# Psi_0 = 1, conj((lambda/.)_3), conj((lambda/.)_3)(omega/.)_3; it is 18-periodic and equals (lambda/x)_3, 1, (omega/x)_3.
PR9 = [p for p in primary_primes(40000) if p != (-2, 0)]
def chi6x(x, a):                     # chi_x(a) = (a/x)_6 = prod over the primes of x (with multiplicity)
    r = 1+0j; m = x
    for p in PR9:
        if N(m) == 1: break
        while divides(p, m):
            k_ = sym6_prime(a, p)
            if k_ is None: return 0
            r *= cmath.exp(1j*math.pi*k_/3); m = div(m, p)
    assert N(m) == 1
    return r
def cub(x, a):                       # (a/x)_3 = chi_x(a)^2 for primary x coprime to 6, from sextic symbols
    return chi6x(x, a)**2
def phi_def(x, tw):
    if not (x[0] % 3 == 1 and x[1] % 3 == 0) or (x[0] % 2 == 0 and x[1] % 2 == 0): return 0
    psi0 = {'1': 1, 'A': np.conj(cub(x, LAM)), 'B': np.conj(cub(x, LAM))*cub(x, (0, 1))}[tw]
    return chi6x(x, LAM)**2*psi0
random.seed(11); bad9 = 0; badf9 = 0; n9 = 0
for _ in range(150):
    x = (3*random.randint(-20, 20)+1, 3*random.randint(-20, 20))
    if x[0] % 2 == 0 and x[1] % 2 == 0: continue
    y = (random.randint(-3, 3), random.randint(-3, 3)); x2 = add(x, mul((18, 0), y)); x3 = add(x, mul((6, 0), y))
    for tw in ('1', 'A', 'B'):
        v1, v2 = phi_def(x, tw), phi_def(x2, tw)
        closed = {'1': lam3(x), 'A': 1, 'B': om3(x)}[tw]
        bad9 += abs(v1-v2) > 1e-9 or abs(v1-closed) > 1e-9; n9 += 1
        if tw == '1' and abs(phi_def(x3, tw)-v1) > 1e-9: badf9 += 1         # falsifier: period 6 for Psi_0 = 1
check("K9 eq:d122 (phi is L-periodic, L = 18)", bad9 == 0 and badf9 > 10,
      f"{n9} (x, twist) cases: phi(x) = phi(x+18y) and = closed form ((lambda/x)_3, 1, (omega/x)_3): {bad9} failures; "
      f"falsifier period 6 for Psi_0 = 1: {badf9} failures")
nf = sum(1 for _, ok in CHECKS if not ok)
print(f"{len(CHECKS)} checks, {nf} failed; time {time.time()-t0:.1f} s")
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1': nf += 1; print("CM_TESTS_FORCE_FAIL=1")
print("PASS" if nf == 0 else "FAIL")
raise SystemExit(0 if nf == 0 else 1)
