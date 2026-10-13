"""t6 (2026-10-13): finite size check (sanity only; the statements are asymptotic) of Lemma lem:s1-unmarked, eq:d137, and
of eq:s1-cm-completed, for the CM unmarked completed rows eq:s1-central

   C^(2)_V(X;m,nu) = sum_{c sf, n; (cn,6)=1} gamma_2(c) conj(alpha(c n^3))^2 nu(c n^3) chi_{cn^3}(m) V(N(c)N(n)^3/X)/(sqrt(N c) N n)

with nu = Psi_0 = conj((lambda/.)_3) (extended by zero), the log-Schwartz test V(y) = exp(-y-1/y), X = Z, and ALL nonzero
rows m in O with N(m) <= Q (units, lambda- and 2-parts and powerful parts included).

Printed: the ratios  sum_{N(m)<=Q} |C|^2 / (Q + Q^2/Z)   (eq:s1-cm-completed with Q = Z^M, M = log Q/log Z), and for the
O-dyads of eq:d137 (N(m_pow) in [P, 2P)) the ratios  sum |C|^2 / (P^{1/2} max(Q/P, Q^2/(Z P)))  = Z^{-(O/2+max{M-O,2M-O-1})}
times the dyadic sum (with Z^O = P).  The check passes if all ratios are below 10 (an arbitrary constant: the lemma has an
unspecified constant and a factor Z^eps).  No falsifier: an upper bound with an unspecified constant has none.
Usage: python3 t6_moment_sanity.py     (about 1 minute)
"""
import os, sys, time, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode = True
from t6_lib import S1
from zw import mul, N, cx, divides, div, primary, UNITS, pw
from fastgauss import primary_primes, g_prime, squarefree_list, cubic_prime, pairsym
import numpy as np
t0 = time.time()
LAMp = (1, 2)
YCUT = 30.0                                            # V(y) = exp(-y-1/y) < 1e-13 for y > 30
def V(y): return np.exp(-y-1/y)
res = {}
ok = True
for Z in (100.0, 300.0):
    B = int(YCUT*Z)
    plist = [p for p in primary_primes(B) if p != (-2, 0)]
    gp = {p: g_prime(p) for p in plist}
    sq = [(c, f) for c, f in squarefree_list(B, plist)]
    nl = [(x, y) for x in range(-40, 41) for y in range(-40, 41) if primary((x, y)) and (x % 2 or y % 2) and N((x, y))**3 <= B]
    # rows m: all nonzero elements with N(m) <= Qmax
    Qmax = 8*Z
    R = int(math.sqrt(4*Qmax/3))+2
    M = [(a, b) for a in range(-R, R+1) for b in range(-R, R+1) if 0 < N((a, b)) <= Qmax]
    MA = np.array([m[0] for m in M], np.int64); MB = np.array([m[1] for m in M], np.int64); MN = np.array([N(m) for m in M])
    # sextic symbol index (m/p)_6 for every prime p used and every row m (-1 where p | m)
    def sym6_vec(p):
        q = N(p)
        if p[1] != 0:                                   # split: O/p = F_q, omega -> w0
            w0 = (-p[0]*pow(p[1], -1, q)) % q
            r = (MA+MB*w0) % q
            e = (q-1)//6; out = np.ones_like(r); bb = r.copy()
            while e:
                if e & 1: out = (out*bb) % q
                bb = (bb*bb) % q; e >>= 1
            tab = {(u[0]+u[1]*w0) % q: j for j, u in enumerate(UNITS)}
            idx = np.array([tab.get(int(v), -2) for v in out]); idx[r == 0] = -1
            assert not np.any(idx == -2)
            return idx
        qq = -p[0]                                      # inert: O/p = F_{q^2}, elements a + b omega mod qq
        a = MA % qq; b = MB % qq; e = (q-1)//6
        ra = np.ones_like(a); rb = np.zeros_like(a); ba, bb = a.copy(), b.copy()
        while e:
            if e & 1: ra, rb = (ra*ba-rb*bb) % qq, (ra*bb+rb*ba-rb*bb) % qq
            ba, bb = (ba*ba-bb*bb) % qq, (2*ba*bb-bb*bb) % qq; e >>= 1
        tab = {(u[0] % qq, u[1] % qq): j for j, u in enumerate(UNITS)}
        idx = np.array([tab.get((int(x), int(y)), -2) for x, y in zip(ra, rb)]); idx[(a == 0) & (b == 0)] = -1
        assert not np.any(idx == -2)
        return idx
    IDX = {p: sym6_vec(p) for p in plist}
    def chi_rows(primes, mult):                        # chi_x(m) for x = prod p^mult, as a vector over the rows
        s = np.zeros(len(M), np.int64); z = np.zeros(len(M), bool)
        for p in primes:
            ip = IDX[p]; z |= ip < 0; s += mult*ip
        out = np.exp(1j*np.pi*s/3); out[z] = 0
        return out
    def lam3(x): return cmath.exp(2j*math.pi/3)**((-(x[1]//3)) % 3)
    C = np.zeros(len(M), complex)
    fac_n = {}
    for n in nl:
        f = []; r = n
        for p in plist:
            if N(r) == 1: break
            while divides(p, r): r = div(r, p); f.append(p)
        fac_n[n] = f
    for c, f in sq:
        nc = N(c)
        g = 1+0j; slam = 1
        for p in f: g *= gp[p]; slam *= cubic_prime(LAMp, p)
        for i in range(len(f)):
            for j in range(i+1, len(f)): g *= pairsym(f[i], f[j])*pairsym(f[j], f[i])
        gam2 = slam*g/math.sqrt(nc)
        chic = None
        for n in nl:
            y = nc*N(n)**3/Z
            if y > YCUT: continue
            A = mul(c, pw(n, 3)); Ac = complex(cx(A))
            coef = gam2*np.conj(Ac/abs(Ac))**2*np.conj(lam3(A))*V(y)/(math.sqrt(nc)*N(n))
            if chic is None: chic = chi_rows(f, 1)
            C += coef*chic*chi_rows(fac_n[n], 3)
    # powerful part of each row
    allp = primary_primes(int(Qmax)+1)
    def mpow(m):
        r = m; P = 1
        for p in [LAMp]+allp:
            if N(r) == 1: break
            e = 0
            while divides(p, r): r = div(r, p); e += 1
            if e >= 2: P *= N(p)**e
        return P
    MP = np.array([mpow(m) for m in M])
    for Q in (Z/8, Z/2, Z, 2*Z, 8*Z):
        sel = MN <= Q
        tot = np.sum(np.abs(C[sel])**2)
        r = tot/(Q+Q*Q/Z); res[(Z, Q)] = r; ok &= r < 10
        worstO = 0
        for k in range(0, 20):
            P0 = 2**k
            s2 = sel & (MP >= P0) & (MP < 2*P0)
            if not np.any(s2): continue
            rO = np.sum(np.abs(C[s2])**2)/(P0**0.5*max(Q/P0, Q*Q/(Z*P0)))
            worstO = max(worstO, rO)
        ok &= worstO < 10
        print(f"Z={Z:g} Q={Q:g} ({np.sum(sel)} rows): sum|C|^2 = {tot:.4g}, ratio to Q+Q^2/Z = {r:.3f}; worst O-dyad ratio {worstO:.3f}", flush=True)
    print(f"[{time.time()-t0:.1f} s]")
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1': ok = False; print("CM_TESTS_FORCE_FAIL=1")
print("PASS" if ok else "FAIL", f"time {time.time()-t0:.1f} s")
raise SystemExit(0 if ok else 1)
