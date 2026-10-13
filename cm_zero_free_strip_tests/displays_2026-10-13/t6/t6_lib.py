"""t6 shared setup (2026-10-13): theta coefficient tables, translates and cusp data for the CM transformation formula
of Appendix app:s1-general (Proposition prop:s1-reflection), in the setting of
tests_cm/phase2_2026-10-12/s1/t_s1_endtoend_active.py:

  base twist 'A':  Psi_0(n) = conj((lambda/n)_3), so phi(x) = chi_x(lambda)^2 Psi_0(x) = 1 on primary x prime to 6
                   (period 6, used with L = 18);  one local prime pi of norm 7 with exponent j in {0..5}.

The code below is copied (unchanged in substance, reorganised into functions) from t_s1_endtoend_active.py, whose
end-to-end identity validates it.  The helpers fastgauss.py and zw.py are imported read-only from
tests_cm/phase2_2026-10-12/s1/.  Nothing in other suites is modified.
"""
import os, sys, math, cmath
HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.normpath(os.path.join(HERE, '..', '..', 'phase2_2026-10-12', 's1'))
sys.path.insert(0, S1)
sys.dont_write_bytecode = True
from fastgauss import primary_primes, g_prime, squarefree_list, cubic_prime, pairsym
from zw import mul, add, sub, conj, N, cx, divides, div, pw, primary, canon_mod, residues, reduce_mod, sym6_prime
import numpy as np

W = cmath.exp(2j*math.pi/3)
LAMp = (1, 2); LAMc = complex(cx(LAMp))
E9 = cmath.exp(2j*math.pi/9)
UE = {1: (1, 0), -1: (-1, 0), W: (0, 1), -W: (0, -1), W*W: (-1, -1), -W*W: (1, 1)}


def ce(z): return np.exp(4j*np.pi*np.real(z))
def e_add(z): return np.exp(4j*np.pi*np.imag(z)/math.sqrt(3))


def ukey(u):
    for k, v in UE.items():
        if abs(k-u) < 1e-9: return v
    raise ValueError(u)


class Setup:
    """Tables to norm B and the arithmetic data of the s1 end-to-end test."""
    def __init__(self, B, pidx=0, verbose=True):
        self.B = B
        self.plist = primary_primes(B)
        self.PI = [p for p in self.plist if N(p) == 7][pidx]; self.NPI = 7; self.PIc = complex(cx(self.PI))
        gp = {p: g_prime(p) for p in self.plist}
        self.sq = squarefree_list(B, self.plist)
        L2 = (-3, 0); wL2 = mul((0, 1), L2); w2L2 = mul((-1, -1), L2)
        cdata = []
        for c, f in self.sq:
            g = 1+0j; s1 = s2 = s3 = 1; slam = 1
            for p in f:
                g *= gp[p]; s1 *= cubic_prime(L2, p); s2 *= cubic_prime(wL2, p); s3 *= cubic_prime(w2L2, p)
                slam *= cubic_prime(LAMp, p)
            for i in range(len(f)):
                for j in range(i+1, len(f)): g *= pairsym(f[i], f[j])*pairsym(f[j], f[i])
            cdata.append((c, N(c), g, np.conj(s1)*g, np.conj(s2)*g, np.conj(s3)*g, slam, (-2, 0) in f, self.PI in f))
        self.cdata = cdata
        self.sqf = {c: f for c, f in self.sq}
        Dmax = int(round(B**(1/3)))+2
        dl = []
        for x in range(-2*Dmax, 2*Dmax+1):
            for y in range(-2*Dmax, 2*Dmax+1):
                if primary((x, y)) and N((x, y))**3 <= B: dl.append(((x, y), N((x, y)), (x % 2 == 0 and y % 2 == 0)))
        self.dl = dl
        # residue arithmetic modulo pi (Z[omega]/pi = F_7)
        P = self.PI; NP = self.NPI
        self.w0 = (-P[0]*pow(P[1], -1, NP)) % NP
        CHI = np.zeros(NP, complex)
        for r in range(1, NP): CHI[r] = cmath.exp(1j*math.pi*sym6_prime((r, 0), P)/3)
        self.CHI = CHI
        self.EPI = np.array([e_add(t/self.PIc) for t in range(NP)])
        # twist A, L = 18
        self.L18 = (18, 0); self.RES = residues(self.L18); self.RESP = residues(P)
        phivals = [(complex(cx(x)), self.phi(x)) for x in self.RES]
        self.PHIHAT = {h0: sum(v*e_add(-complex(cx(h0))*xc/18) for xc, v in phivals if v != 0)/324 for h0 in self.RES}
        self.H0S = [h0 for h0 in self.RES if abs(self.PHIHAT[h0]) > 1e-13]
        self.TR = {}
        for h0 in self.H0S:
            for hp in self.RESP:
                self.TR[(h0, hp)] = self.translate(h0, hp)
        if verbose:
            print(f"setup: B={B}, pi={P}, {len(self.cdata)} squarefree c, {len(self.TR)} translates, {len(self.H0S)} h0 classes")

    # ---------------- theta coefficient tables (identity cusp d_0 and the cusps d_+, d_-) ----------------
    def table_d0(self, Nmax, meta=False):
        L = []; Cf = []; XE = []; MK = []
        for cc_, nc, g1_, gl2, gwl2, gw2l2, _, _, _ in self.cdata:
            if nc > 27*Nmax: continue
            cc = complex(cx(cc_))
            for dd_, nd, _ in self.dl:
                base = nc*nd**3
                if base > 27*Nmax: continue
                dd = complex(cx(dd_)); absdc = abs(dd)/abs(cc); cd3 = mul(cc_, pw(dd_, 3))
                for k in range(-3, 40):
                    if base*3.0**k > Nmax*(1+1e-12): break
                    r = k % 3
                    if r == 1: continue
                    if r == 2:
                        nu = (k+4)//3
                        vals = [(1, np.conj(gl2)), (-1, np.conj(gl2)), (W, E9**-1*np.conj(gwl2)), (-W, E9**-1*np.conj(gwl2)),
                                (W*W, E9*np.conj(gw2l2)), (-W*W, E9*np.conj(gw2l2))]
                        fac = absdc*3**(nu/2+2)
                    else:
                        nu = (k+3)//3
                        vals = [(1, np.conj(g1_)), (-1, np.conj(g1_))]
                        fac = absdc*3**(nu/2+5/2)
                    lk = LAMc**k*cc*dd**3; xk = mul(pw(LAMp, k+4), cd3)
                    for u, tv in vals:
                        L.append(u*lk); Cf.append(np.conj(tv*fac)); XE.append(mul(ukey(u), xk)); MK.append((k, nc, nd))
        if meta: return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64), np.array(MK)
        return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64)

    def table_dpm(self, sign, Nmax, meta=False):
        L = []; Cf = []; XE = []; MK = []
        for cc_, nc, g1_, gl2, gwl2, gw2l2, _, _, _ in self.cdata:
            if nc > 81*Nmax: continue
            cc = complex(cx(cc_))
            for dd_, nd, _ in self.dl:
                base = nc*nd**3
                if base > 81*Nmax: continue
                dd = complex(cx(dd_)); absdc = abs(dd)/abs(cc); m0 = LAMc**-4*cc*dd**3; cd3 = mul(cc_, pw(dd_, 3))
                if sign == '+':
                    for up, t2 in [(1, 9*W*W*np.conj(gl2)), (W, 9*E9**-1*np.conj(gwl2)), (W*W, 9*W*E9*np.conj(gw2l2))]:
                        l = up*W*W*m0; L.append(l); Cf.append(W*W*np.conj(t2*absdc)*ce(l)); XE.append(mul(ukey(up*W*W), cd3))
                        MK.append((-4, nc, nd))
                else:
                    for up, t1 in [(1, 9*W*np.conj(gl2)), (W, 9*E9**-1*W*W*np.conj(gwl2)), (W*W, 9*E9*np.conj(gw2l2))]:
                        l = -up*W*m0; L.append(l); Cf.append(W*np.conj(t1*absdc)*ce(l)); XE.append(mul(ukey(-up*W), cd3))
                        MK.append((-4, nc, nd))
        if meta: return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64), np.array(MK)
        return np.array(L), np.array(Cf), np.array(XE, dtype=np.int64)

    def table(self, sigma, Nmax, meta=False):
        return self.table_d0(Nmax, meta) if sigma == '0' else self.table_dpm(sigma, Nmax, meta)

    # ---------------- arithmetic of the twist ----------------
    @staticmethod
    def lam3(x): return W**((-(x[1]//3)) % 3)          # (lambda/x)_3 for primary x
    def psi0(self, x): return np.conj(self.lam3(x))

    @staticmethod
    def phi(x):                                          # chi_x(lambda)^2 Psi_0(x) on x = 1 (3), (x,6)=1; 0 otherwise
        if not (x[0] % 3 == 1 and x[1] % 3 == 0): return 0
        if x[0] % 2 == 0 and x[1] % 2 == 0: return 0
        return 1

    def red(self, x): return (x[0]+x[1]*self.w0) % self.NPI

    def chipow(self, r, j):
        r = np.asarray(r) % self.NPI
        return np.where(r == 0, 0, self.CHI[r]**j)

    def Cpj(self, h, j):
        return sum(self.chipow(self.red(y), j)*e_add(-complex(cx(h))*complex(cx(y))/self.PIc) for y in self.RESP)/self.NPI

    def gauss_pm(self, m, sgn):
        return sum(self.chipow(self.red(y), m % 6)*e_add(sgn*complex(cx(y))/self.PIc) for y in self.RESP)/math.sqrt(self.NPI)

    # ---------------- cusp data ----------------
    @staticmethod
    def mod3eq(x, y): d = sub(x, y); return d[0] % 3 == 0 and d[1] % 3 == 0

    def vlam(self, x):
        k = 0
        while divides(LAMp, x): x = div(x, LAMp); k += 1
        return k

    def cubic_general(self, x, a):
        r = 1; m = a
        for p in self.plist:
            if N(m) == 1: break
            while divides(p, m):
                r *= cubic_prime(x, p); m = div(m, p)
        assert N(m) == 1
        return r

    @staticmethod
    def matmul(A, Bm):
        (a, b), (c, d) = A; (e, f), (g, h) = Bm
        return ((add(mul(a, e), mul(b, g)), add(mul(a, f), mul(b, h))), (add(mul(c, e), mul(d, g)), add(mul(c, f), mul(d, h))))

    @staticmethod
    def inv(A):
        (a, b), (c, d) = A; return ((d, (-b[0], -b[1])), ((-c[0], -c[1]), a))

    def cusp_data(self, a, C):
        """sigma, delta', kappa(G), and also b', H (returned for the automorphy test)."""
        v = self.vlam(C); one = (1, 0); zero = (0, 0); m3 = self.mod3eq
        if v >= 2:
            assert m3(a, one); H = ((one, zero), (zero, one)); sigma = '0'
        elif v == 1:
            assert m3(a, one)
            u0 = LAMp if m3(C, LAMp) else (-1, -2); assert m3(C, u0)
            H = ((one, zero), (u0, one)); sigma = '-' if u0 == LAMp else '+'
        else:
            assert primary(C)
            u0 = (a[0] % 3, a[1] % 3)
            eps = {0: 0, 1: 1, 2: -1}[u0[1]]
            H = ((u0, (-1, 0)), (one, zero)); sigma = {0: '0', 1: '-', -1: '+'}[eps]
        best = None
        for Rb in (14, 30):
            for x in range(-Rb, Rb+1):
                for y in range(-Rb, Rb+1):
                    dp = (x, y); t = sub(mul(a, dp), one)
                    if not divides(C, t): continue
                    bg = div(t, C); g = ((a, bg), (C, dp)); (A1, B1), (C1, D1) = self.matmul(g, self.inv(H))
                    if m3(A1, one) and m3(D1, one) and m3(B1, zero) and m3(C1, zero):
                        if best is None or N(A1) < best[0]: best = (N(A1), dp, A1, C1, bg)
            if best is not None and best[0] <= self.B: break
        _, dp, A1, C1, bg = best
        kap = 1 if C1 == (0, 0) else self.cubic_general(C1, A1)
        return sigma, dp, kap, bg, H

    def gcd_O(self, x, y):
        while y != (0, 0):
            x, y = y, reduce_mod(x, y)
        return x

    def normalize(self, a, c):
        U6 = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]
        if divides(LAMp, c):
            for u in U6:
                if self.mod3eq(mul(u, a), (1, 0)): return mul(u, a), mul(u, c)
        else:
            for u in U6:
                if self.mod3eq(mul(u, c), (1, 0)): return mul(u, a), mul(u, c)
        raise ValueError

    def translate(self, h0, hp):
        P = self.PI
        num = mul(mul(LAMp, LAMp), add(mul(h0, P), mul((18, 0), hp))); den = mul(self.L18, P)
        if num == (0, 0): a, c = (0, 0), (1, 0)
        else:
            g = self.gcd_O(num, den); a, c = div(num, g), div(den, g)
        a, c = self.normalize(a, c)
        a = add(reduce_mod(sub(a, (1, 0)), mul((3, 0), c)), (1, 0)) if divides(LAMp, c) else reduce_mod(a, mul((3, 0), c))
        sigma, dp, kap, bg, H = self.cusp_data(a, c)
        return dict(a=a, c=c, sigma=sigma, dp=dp, kap=kap, bp=bg, H=H)


def act(g, z, v):
    """Action of g = (a, b, c, d) in SL_2(C) on (z, v) in H^3 ([QF (5.1)])."""
    a, b, c, d = g
    den = abs(c*z+d)**2+abs(c)**2*v*v
    return ((a*z+b)*np.conj(c*z+d)+a*np.conj(c)*v*v)/den, v/den
