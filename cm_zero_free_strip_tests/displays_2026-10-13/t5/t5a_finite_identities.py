"""t5a: finite (exact or machine-precision) identities and explicit inequalities of Appendix A.

Displays tested (labels of cm-eleven-twelfths-guide.tex, numbers of the 2026-10-13 build):
  (A.1)  eq:d096      Psi = xi_0 prod chi_pi^{j_pi} (chi^0 = 1_{pi nmid n}) is completely multiplicative on primary
                      elements and vanishes unless (n, 6) = 1 (consistency check of the definition)
  (A.3)  eq:d097      |C'| <= 1: |phihat(h0)| <= 1 (all h0 mod 72, five twists), |varpi_{pi,j}| = 1, |kappa0| = 1
  (A.12) eq:b3-coeff  coefficient of Theta_Psi at lambda^{-3} n b^3: d_0 (QF table) * phi * prod chi^j equals
                      3^{5/2}|b| gamma_2(n) Psi(n) Psi(b)^3 (gamma_2 by brute force), and the symbol identity
                      conj((lambda/n)_3) (lambda/n b^3)_3 Psi(n b^3) = Psi(n) Psi(b)^3
  (A.13)/(A.40) eq:d104, eq:d120  C_{pi,0}(0) = 1 - 1/N(pi), C_{pi,j}(0) = 0 (j != 0), so c_F restricted to a group is (A.40)
  (A.18) eq:b3-local  the local sum over h != 0 equals chi(sigma)^{-2} varpi B_{pi,j}(x), varpi independent of x, |varpi| = 1
  (A.23) eq:d108      conj(alpha(l))^2 = -conj(alpha(n))^2 conj(alpha(b))^6, N(l)^{1/2-s} = 27^{s-1/2} N(n)^{1/2-s} N(b)^{3/2-3s}
  (A.33) eq:d114      termwise bound and partial sums of sum_l |d_sigma(l)| N(l)^{1/2-s'} against 6*27 sum_m sum_n sum_b ...
                      (the termwise bound is equivalent to (A.5), eq:d099, |d_sigma| <= 27 3^{m/6}|b|, here for N(9l) <= 30000)
  (A.27) eq:d111      the second-order chain rule (sympy, generic polynomial F, at general z)
  (A.2)  count        N(72) = 5184 (the bound 5184 * 2^{|P|} on the number of groups)
Falsifiers (must fail): Psi(b) in place of Psi(b)^3; (lambda/n)_3 unconjugated; B_{pi,j} with chi^{-j}; chi(sigma h)^{+2};
exponent 3 on conj(alpha(b)); the bound with 9 in place of 27; mixed terms counted once in (A.27).
Usage: python3 -B t5a_finite_identities.py      (about 1 minute)
"""
import os, sys, time, math, cmath, random
t0 = time.time()
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t5lib import *
import cmlib                                   # read-only (part B): gamma_j_fast (vectorized brute force)
from collections import Counter
random.seed(20261013)

PI7 = [p for p in fg.primary_primes(7) if N(p) == 7][0]
PI13 = [p for p in fg.primary_primes(13) if N(p) == 13][0]
PI19 = [p for p in fg.primary_primes(19) if N(p) == 19][0]
PIIN = (-5, 0)
XI = {
    'xi0=1': lambda x: 1,
    'xi0=conj(lam/.)_3': lambda x: np.conj(lam3_euler(x)),
    'xi0=conj(lam/.)_3 (w/.)_3': lambda x: np.conj(lam3_euler(x))*om3(x),
    'xi0=psi4 (w/.)_3': lambda x: psi4(x)*om3(x),
    'xi0=(2/.)_3 psi4': lambda x: cubic_sym((-2, 0), x)*psi4(x),
}
def Psi_fn(xi0, P):
    def Psi(n):
        if not primary(n) or divides((2, 0), n): return 0
        val = xi0(n)
        for p, j in P:
            ch = sext_prime(n, p)
            if ch == 0: return 0                  # chi^j vanishes at multiples of pi, including j = 0 (1_{pi nmid n})
            val *= ch**j
        return val
    return Psi

# ------------------------------------------------------------------------------------------------ (A.2)
report("t5a.A2", "eq:reflection count", N((72, 0)) == 5184, f"N(72) = {N((72, 0))} = 5184 residues h0 mod 72 (L = 72)")

# ------------------------------------------------------------------------------------------------ (A.1)
prim = [(a, b) for a in range(-40, 41) for b in range(-40, 41) if primary((a, b)) and 1 <= N((a, b)) <= 1200]
Pset = [(PI7, 2), (PI13, 0), (PIIN, 5)]
bad_mult = 0; bad_zero = 0; cnt = 0
for name, xi0 in XI.items():
    Psi = Psi_fn(xi0, Pset)
    for _ in range(300):
        n1 = random.choice(prim); n2 = random.choice(prim)
        cnt += 1
        if abs(Psi(mul(n1, n2)) - Psi(n1)*Psi(n2)) > 1e-12: bad_mult += 1
    for n in prim:
        if divides((2, 0), n) and Psi(n) != 0: bad_zero += 1
# xi_0 must itself be multiplicative on primary elements coprime to 6 and periodic mod 72 (a character of frak G)
bad_xi = 0
for name, xi0 in XI.items():
    for _ in range(200):
        n1 = random.choice(prim); n2 = random.choice(prim)
        if divides((2, 0), n1) or divides((2, 0), n2): continue
        if abs(xi0(mul(n1, n2)) - xi0(n1)*xi0(n2)) > 1e-12: bad_xi += 1
        m = add(n1, mul((72, 0), random.choice([(1, 0), (0, 1), (3, -2)])))
        if abs(xi0(m) - xi0(n1)) > 1e-12: bad_xi += 1
report("t5a.A1", "eq:d096", bad_mult == 0 and bad_zero == 0 and bad_xi == 0,
       f"5 characters xi_0 of frak G (mod 72) x P = {{pi7: j=2, pi13: j=0, -5: j=5}}: Psi(n1 n2) = Psi(n1)Psi(n2) for {cnt} random "
       f"primary pairs (N <= 1200, even included): {bad_mult} failures; Psi(n) = 0 for even n: {bad_zero} failures; xi_0 multiplicative "
       f"and 72-periodic: {bad_xi} failures")

# ------------------------------------------------------------------------------------------------ (A.12)
TT = ThetaTables(9000, verbose=False)
tb = TT.table('0')
dmap = dict(zip(tb['m'], tb['d']))
g2cache = {}
def gamma2(n):
    if n not in g2cache: g2cache[n] = cmlib.gamma_j_fast(2, n)
    return g2cache[n]
errA12 = 0; errSym = 0; fPsi = 0; fConj = 0; cnt12 = 0; sel_bad = 0
P12 = [(PI7, 1), (PI13, 4), (PIIN, 0)]
xs = [(a, b) for a in range(-64, 65) for b in range(-64, 65) if prime_coprime6((a, b)) and N((a, b)) <= 3000]
for name, xi0 in list(XI.items())[:4]:
    Psi = Psi_fn(xi0, P12)
    for x in xs:
        m = mul(LAM, x)                                 # m = lambda^4 l with l = lambda^-3 x
        d0 = dmap.get(m, 0)
        # x = n b^3 with n squarefree, n, b primary: decomposition by factoring
        _, ps = factor_prim(x)
        cnts = Counter(ps); n = (1, 0); b = (1, 0)
        for p, e in cnts.items():
            q, r = divmod(e, 3)
            n = mul(n, pw(p, r)); b = mul(b, pw(p, q))
        u = div(x, mul(n, pw(b, 3)))
        if any(e % 3 == 2 for e in cnts.values()):
            sel_bad += abs(d0) > 1e-12                 # selection: no coefficient unless x = n b^3, n squarefree
            continue
        assert u == (1, 0)
        n = fg.prim_assoc(n); b = fg.prim_assoc(b)
        coefT = d0*lam3_euler(x)*xi0(x)*np.prod([(sext_prime(x, p)**j if sext_prime(x, p) != 0 else 0) for p, j in P12])
        rhs = 3**2.5*abs(cx(b))*gamma2(n)*Psi(n)*Psi(b)**3
        errA12 = max(errA12, abs(coefT - rhs)); cnt12 += 1
        # the symbol identity in (A.12): conj((lambda/n)_3) (lambda/n b^3)_3 Psi(n b^3) = Psi(n) Psi(b)^3
        lhs_s = np.conj(lam3_euler(n))*lam3_euler(x)*Psi(x)
        errSym = max(errSym, abs(lhs_s - Psi(n)*Psi(b)**3))
        fPsi = max(fPsi, abs(lhs_s - Psi(n)*Psi(b)))
        fConj = max(fConj, abs(lam3_euler(n)*lam3_euler(x)*Psi(x) - Psi(n)*Psi(b)**3))
report("t5a.A12", "eq:b3-coeff", errA12 < 1e-10 and errSym < 1e-12 and sel_bad == 0,
       f"{cnt12} indices l = lambda^-3 n b^3 (all primary x = nb^3 coprime to 6, N(x) <= 3000, 4 xi_0, P = {{pi7: 1, pi13: 4, -5: 0}}): "
       f"d_0 (QF table) * phi * prod chi^j vs 3^(5/2)|b| gamma_2(n) Psi(n) Psi(b)^3 (gamma_2 brute force): max abs err {errA12:.1e}; "
       f"symbol identity: {errSym:.1e}; x not of the form nb^3 with nonzero d_0: {sel_bad}")
falsifier("t5a.A12", "Psi(b) in place of Psi(b)^3", fPsi, 1e-3)
falsifier("t5a.A12", "(lambda/n)_3 unconjugated", fConj, 1e-3)

# ------------------------------------------------------------------------------------------------ (A.13), (A.40), (A.3) phihat
def residues_pi(p):
    n = N(p)
    return sorted(set(canon_mod((a, b), p) for a in range(n) for b in range(n)))
PRIMES = [PI7, PI13, PI19, PIIN]
RP = {p: residues_pi(p) for p in PRIMES}
def Ccoef(p, j, h):
    s = 0
    for y in RP[p]:
        ch = sext_prime(y, p)
        if ch == 0: continue
        s += ch**j*e_add(-complex(cx(mul(h, y)))/complex(cx(p)))
    return s/N(p)
e40 = 0
for p in PRIMES:
    e40 = max(e40, abs(Ccoef(p, 0, (0, 0)) - (1 - 1/N(p))))
    for j in range(1, 6):
        e40 = max(e40, abs(Ccoef(p, j, (0, 0))))
report("t5a.A40", "eq:d120 (and eq:d104)", e40 < 1e-12,
       f"C_pi,0(0) = 1 - 1/N(pi) and C_pi,j(0) = 0 for j = 1..5, primes of norm 7, 13, 19, 25 (inert -5): max err {e40:.1e}; "
       f"so c_F(h) on a group (h_pi = 0 off A, j_pi = 0 there) is (A.40)")
RES72 = [(a, b) for a in range(72) for b in range(72)]
XA = np.array([x[0] for x in RES72]); XB = np.array([x[1] for x in RES72])
def phihat_all(phi_vals):
    out = np.zeros(len(RES72), complex)
    for i0 in range(0, len(RES72), 400):
        ha = XA[i0:i0 + 400, None]; hb = XB[i0:i0 + 400, None]
        Q = (XA[None, :]*hb + XB[None, :]*(ha - hb)) % 72
        out[i0:i0 + 400] = (np.exp(-2j*np.pi*Q/72) @ phi_vals)/5184
    return out
phmax = 0; phinfo = []
for name, xi0 in XI.items():
    pv = np.array([(lam3_euler(x)*xi0(x) if prime_coprime6(x) else 0) for x in RES72])
    ph = phihat_all(pv)
    phmax = max(phmax, np.max(np.abs(ph))); phinfo.append(f"{np.max(np.abs(ph)):.3f}")

# ------------------------------------------------------------------------------------------------ (A.18) and |varpi| = 1
e18 = 0; f18a = 0; f18b = 0; ncase = 0; varpi_mod = 0
for p in PRIMES:
    res = RP[p]; nz = [h for h in res if h != (0, 0)]
    pc = complex(cx(p)); Np = N(p)
    inv = {h: inv_mod(h, p) for h in nz}
    chi = {y: sext_prime(y, p) for y in res}
    Cj = {(j, h): Ccoef(p, j, h) for j in range(6) for h in nz}
    def B(j, x, wrong=False):
        ch = chi[canon_mod(x, p)]
        if j not in (0, 4) or wrong:
            e = -j if wrong else -j - 2
            return ch**e if ch != 0 else 0
        if j == 0: return Np**-0.5*(ch**-2 if ch != 0 else 0)
        return Np**-0.5*(-1 + (Np if ch == 0 else 0))
    for eps in random.sample(nz, 3):
        for sig in random.sample(nz, 2):
            for j in range(6):
                ratios = []; ratios_w = []; ratios_k = []
                for x in res:
                    S = sum(Cj[(j, h)]*chi[canon_mod(mul(sig, h), p)]**-2*e_add(complex(cx(mul(mul(eps, inv[h]), x)))/pc) for h in nz)
                    Sk = sum(Cj[(j, h)]*chi[canon_mod(mul(sig, h), p)]**2*e_add(complex(cx(mul(mul(eps, inv[h]), x)))/pc) for h in nz)
                    base = chi[canon_mod(sig, p)]**-2*B(j, x)
                    ncase += 1
                    if abs(base) < 1e-12:
                        e18 = max(e18, abs(S)); continue
                    ratios.append(S/base); ratios_k.append(Sk/base)
                    bw = chi[canon_mod(sig, p)]**-2*B(j, x, wrong=True)
                    if abs(bw) > 1e-12: ratios_w.append(S/bw)
                r0 = ratios[0]
                e18 = max(e18, max(abs(r - r0) for r in ratios), abs(abs(r0) - 1))
                varpi_mod = max(varpi_mod, abs(abs(r0) - 1))
                if j not in (0, 4) and len(ratios_w) > 1: f18a = max(f18a, max(abs(r - ratios_w[0]) for r in ratios_w))
                f18b = max(f18b, max(abs(r - ratios_k[0]) for r in ratios_k) + abs(abs(ratios_k[0]) - 1))
report("t5a.A18", "eq:b3-local", e18 < 1e-10,
       f"4 primes (norms 7, 13, 19 and inert 25) x 3 eps x 2 sigma_pi x j = 0..5, all x mod pi ({ncase} cases), C_pi,j literal from (A.13): "
       f"sum_(h != 0) = chi(sigma)^-2 varpi B_pi,j(x) with varpi independent of x and |varpi| = 1: max err {e18:.1e}")
falsifier("t5a.A18", "B_pi,j with chi^-j (j not in {0,4})", f18a, 1e-3)
falsifier("t5a.A18", "chi(sigma h)^+2", f18b, 1e-3)
report("t5a.A3", "eq:d097 |C'| <= 1", phmax <= 1 + 1e-12 and varpi_mod < 1e-10,
       f"C' = phihat(h0) conj(kappa0) prod(1 - 1/N) prod chi(sigma)^-2 varpi (eq:b3-C): max |phihat(h0)| over all 5184 h0 for 5 twists "
       f"= {', '.join(phinfo)} (<= 1); |varpi| = 1 to {varpi_mod:.1e}; kappa0 is a cube root of unity (r2 T2)")

# ------------------------------------------------------------------------------------------------ (A.23)
e23 = 0; f23 = 0
ns = [n for n in prim if not divides((2, 0), n)][:120]; bs = [(1, 0), (-2, -3), (4, 3), (-5, 0), mul((-2, -3), (-2, -3)), (1, 3)]
for n in ns:
    for b in bs:
        l = LAMc**-3*complex(cx(n))*complex(cx(b))**3
        al = l/abs(l); an = complex(cx(n))/abs(cx(n)); ab = complex(cx(b))/abs(cx(b))
        e23 = max(e23, abs(np.conj(al)**2 + np.conj(an)**2*np.conj(ab)**6))
        f23 = max(f23, abs(np.conj(al)**2 + np.conj(an)**2*np.conj(ab)**3))
        for s in [0.7 + 2j, 2.5, -1.3 - 0.4j]:
            lhs = abs(l)**(2*(0.5 - s)); rhs = 27**(s - 0.5)*N(n)**(0.5 - s)*N(b)**(1.5 - 3*s)
            e23 = max(e23, abs(lhs - rhs)/abs(rhs))
report("t5a.A23", "eq:d108", e23 < 1e-12, f"{len(ns)*len(bs)} pairs (n, b), three s: max rel err {e23:.1e}")
falsifier("t5a.A23", "conj(alpha(b))^3 (zeta exponent)", f23, 1e-3)

# ------------------------------------------------------------------------------------------------ (A.33)
TT2 = ThetaTables(30000, verbose=False)
viol = 0; viol9 = 0; dup = 0; rows = []
for sp in [1.6, 2.0, 3.0]:
    lhs_tot = 0; rhs_tot = 0
    for sig in '0+-':
        t = TT2.table(sig)
        seen = set()
        for m, d, k, c, b in zip(t['m'], t['d'], t['k'], t['c'], t['b']):
            if abs(d) < 1e-12: continue
            if m in seen: dup += 1
            seen.add(m)
            Nl = N(m)/81.0
            term = abs(d)*Nl**(0.5 - sp)
            bound = 27*3**(k*(2/3 - sp))*N(c)**(0.5 - sp)*N(b)**(2 - 3*sp)
            viol += term > bound*(1 + 1e-9)
            viol9 += term > bound/3*(1 + 1e-9)
            lhs_tot += term
    # right side of (A.33) over the same range N(l) = 3^m N(n) N(b)^3 <= B/81: 6 units, m >= -4, n squarefree primary, b primary
    Bl = TT2.B/81.0
    for c, nc, *_ in TT2.cdata:
        for b in TT2.dl:
            nb = N(b)
            if nc*nb**3*3**-4 > Bl: continue
            for m in range(-4, 80):
                if 3**m*nc*nb**3 > Bl: break
                rhs_tot += 6*27*3**(m*(2/3 - sp))*nc**(0.5 - sp)*nb**(2 - 3*sp)
    # the full right side is finite for s' > 3/2: geometric m-sum times sum_n N(n)^{1/2-s'} times sum_b N(b)^{2-3s'}
    rows.append(f"s'={sp}: LHS {lhs_tot:.4g} <= RHS {rhs_tot:.4g}")
    viol += lhs_tot > rhs_tot
report("t5a.A33", "eq:d114", viol == 0 and dup == 0,
       f"all nonzero d_sigma(l), N(9l) <= 30000, three sigma: termwise |d|N(l)^(1/2-s') <= 27*3^(m(2/3-s'))N(n)^(1/2-s')N(b)^(2-3s') "
       f"(0 violations; representation unique: {dup} duplicates); partial sums over N(l) <= 30000/81: " + "; ".join(rows))
falsifier("t5a.A33", "9 in place of 27 (termwise)", viol9, 0)
# ------------------------------------------------------------------------------------------------ (A.27)
import sympy as sp
z, zb, v, c, cb, d, db = sp.symbols('z zb v c cb d db')
Xs, Ys, Vs = sp.symbols('X Y V')
Zp = -d/c - zb/(c**2*(v**2 + z*zb)); Zbp = -db/cb - z/(cb**2*(v**2 + z*zb)); Vp = v/(c*cb*(v**2 + z*zb))   # (A.25)
rng = random.Random(27)
Fpoly = sum(rng.randint(-9, 9)*Xs**i*Ys**k*Vs**m for i in range(4) for k in range(3) for m in range(3))
sub_ = {Xs: Zp, Ys: Zbp, Vs: Vp}
lhs = sp.diff(Fpoly.subs(sub_), zb, 2)
Yv = [(Xs, Zp), (Ys, Zbp), (Vs, Vp)]
rhs = sum(sp.diff(Fpoly, A_, B_).subs(sub_)*sp.diff(a_, zb)*sp.diff(b_, zb) for A_, a_ in Yv for B_, b_ in Yv) \
    + sum(sp.diff(Fpoly, A_).subs(sub_)*sp.diff(a_, zb, 2) for A_, a_ in Yv)
pt = {z: sp.Rational(3, 7) + 0*sp.I, zb: sp.Rational(-2, 5), v: sp.Rational(5, 4), c: 2 + sp.I, cb: 2 - sp.I, d: sp.Rational(1, 3), db: sp.Rational(1, 3)}
diff27 = [sp.nsimplify(sp.expand((lhs - rhs).subs(pt)))]
pt2 = {z: sp.Rational(-1, 9), zb: sp.Rational(1, 6), v: sp.Rational(2, 3), c: 3, cb: 3, d: 1, db: 1}
diff27.append(sp.nsimplify(sp.expand((lhs - rhs).subs(pt2))))
ok27 = all(x == 0 for x in diff27)
# symbolic at general (z, zbar) for a smaller generic F
F2 = sum(rng.randint(-5, 5)*Xs**i*Ys**k*Vs**m for i in range(3) for k in range(2) for m in range(2))
l2 = sp.diff(F2.subs(sub_), zb, 2)
r2_ = sum(sp.diff(F2, A_, B_).subs(sub_)*sp.diff(a_, zb)*sp.diff(b_, zb) for A_, a_ in Yv for B_, b_ in Yv) \
    + sum(sp.diff(F2, A_).subs(sub_)*sp.diff(a_, zb, 2) for A_, a_ in Yv)
ok27s = sp.simplify(l2 - r2_) == 0
# falsifier: mixed second derivatives counted once (unordered pairs)
r2_w = sum(sp.diff(Fpoly, A_, B_).subs(sub_)*sp.diff(a_, zb)*sp.diff(b_, zb) for i_, (A_, a_) in enumerate(Yv) for B_, b_ in Yv[i_:]) \
    + sum(sp.diff(Fpoly, A_).subs(sub_)*sp.diff(a_, zb, 2) for A_, a_ in Yv)
f27 = abs(complex(sp.N((lhs - r2_w).subs(pt))))
report("t5a.A27", "eq:d111", ok27 and ok27s,
       "sympy: d_zbar^2 F(z',zbar',v') = sum over ordered pairs (Y,Y') of F_YY' Y_zbar Y'_zbar + sum_Y F_Y Y_zbarzbar, with "
       "(z',zbar',v') of (A.25) and z, zbar independent: exact at two rational points for a generic F of degree (3,2,2), and "
       "identically (simplify) for a generic F of degree (2,1,1)")
falsifier("t5a.A27", "mixed terms counted once (unordered pairs)", f27, 1e-6)
finish(t0)
