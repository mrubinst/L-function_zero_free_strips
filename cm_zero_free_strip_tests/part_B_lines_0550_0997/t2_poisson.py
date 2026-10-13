"""Section 'Poisson summation reveals the cubic Gauss sums', lines 585-631: claims B09-B21.

B09  e(z) trivial on O; O self-dual for (z, w) -> e(zw).          Falsifier: exp(2 pi i Im z) is not trivial on O.
B10  Phi^ normalization (2D quadrature, 30 digits), radial.       Falsifier: without the factor 2/sqrt3.
B11  smooth mean square >= M(D,H) (CM rows A_u with alpha), Phi(t) = e^{1-t}.  Falsifier: Phi(t) = e^{-t} (< 1 on [0,1]).
B12  mu_K(n) = 0 for non-squarefree n (only squarefree columns).
B13  (z1, z2) = (g, z1 z2) = 1.                                    Falsifier: a non-squarefree column.
B14  conj(chi_{n1}(u)) chi_{n2}(u) = 1_{(u,g)=1} psi(u).            Falsifier: drop the indicator.
B15  psi primitive modulo z1 z2.
B16  |gamma(psi)| = 1.
B17  exact Poisson formula (eq:poisson-exact), Phi(t) = exp(-pi t), 7 cases; and end to end: the smoothed CM mean square
     sum_u Phi(N(u)/H)|A_u(D)|^2 with A_u = sum mu_K alpha chi_n(u) W equals its dual (Poisson) form.
     Falsifiers: drop e0 != 1; psi(h) instead of conj(psi(h)); N(e0)^2 for N(e0).
B18  psi(0) = 0 unless z1 = z2 = 1.
B19  nonzero frequencies of diagonal pairs (exist with the Gaussian; with supp Phi^ in [0, C] they need H <= C N(e0)).
B20  N(h) << N(e0)N(z1z2)/H << D^2/(H N(g)) (sympy + sampling).  Falsifier: D^2/(H N(g)^2).
B21  H > D => D^2/H < D.
Run time: about 30 seconds.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import *
import mpmath as mp, sympy as sp
from fractions import Fraction as Fr
random.seed(2026)
sq = primaries(400, sqfree=True)
allp = primaries(400)

# ---------------- B09
def e_c(z): return cmath.exp(4j*math.pi*z.imag/SQ3)
err = max(abs(e_c(cx((a, b))) - 1) for a in range(-30, 31) for b in range(-30, 31))
# self-duality: w on the grid (1/12)(Z + Z omega) with e(w) = e(omega w) = 1 iff w in O
bad = 0; cnt = 0
for i in range(-36, 37):
    for j in range(-36, 37):
        w = (i + j*W)/12
        triv = abs(e_c(w) - 1) < 1e-9 and abs(e_c(W*w) - 1) < 1e-9
        inO = (i % 12 == 0 and j % 12 == 0)
        bad += (triv != inO); cnt += 1
report("B09", err + bad, 1e-9, f"|e(a+b omega)-1| <= {err:.1e} for |a|,|b|<=30; dual lattice of O under e(zw) is O ({cnt} grid points of (1/12)O)")
falsifier("B09", "exp(2 pi i Im z) on O", max(abs(cmath.exp(2j*math.pi*cx((a, 1)).imag) - 1) for a in range(3)), 1e-3)

# ---------------- B10
mp.mp.dps = 20
def hatPhi_quad(w, norm=True):
    f = lambda x, y: mp.e**(-mp.pi*(x*x + y*y))*mp.cos(4*mp.pi*mp.im(mp.mpc(x, y)*w)/mp.sqrt(3))
    return (2/mp.sqrt(3) if norm else 1)*mp.quad(f, [-6, 6], [-6, 6])
err = 0; rad = 0
for w in [mp.mpc(0.3, 0.2), mp.mpc(-0.7, 0.5), mp.mpc(1.1, -0.4)]:
    v = hatPhi_quad(w)
    err = max(err, abs(v - 2/mp.sqrt(3)*mp.e**(-4*mp.pi*abs(w)**2/3)))
    rad = max(rad, abs(v - hatPhi_quad(w*mp.e**(1j*0.77))))
report("B10", float(err + rad), 1e-12, f"Phi(t)=exp(-pi t): Phi^(N w)=(2/sqrt3)exp(-4 pi N(w)/3), error {float(err):.1e}; rotation invariance {float(rad):.1e}")
falsifier("B10", "no factor 2/sqrt3", float(abs(hatPhi_quad(mp.mpc(0.3, 0.2), False) - 2/mp.sqrt(3)*mp.e**(-4*mp.pi*0.13/3))), 1e-3)

# ---------------- B12, B13, B14, B15, B16, B18
nonsq = [n for n in allp if mu(n) == 0]
report("B12", sum(1 for n in nonsq if len(set(factor(n))) == len(factor(n))), 0,
       f"{len(nonsq)} non-squarefree primary n (N<=400): each has a repeated prime and mu_K=0; such columns contribute 0")
bad13 = 0; bad14 = 0; err14 = 0; cnt = 0
for _ in range(400):
    n1, n2 = random.choice(sq), random.choice(sq)
    g = gcd_primary(n1, n2); z1 = div(n1, g); z2 = div(n2, g)
    bad13 += N(gcd_primary(z1, z2)) != 1 or N(gcd_primary(g, mul(z1, z2))) != 1
    for _ in range(3):
        u = (random.randint(-40, 40), random.randint(-40, 40))
        lhs = chi(n1, u).conjugate()*chi(n2, u)
        ind = 1 if all(not divides(p, u) for p in factor(g)) else 0
        err14 = max(err14, abs(lhs - ind*chi(z1, u).conjugate()*chi(z2, u))); cnt += 1
report("B13", bad13, 0, "400 random squarefree pairs")
# falsifier: non-squarefree columns
P7 = [p for p in PR_ODD if N(p) == 7][0]; P13 = [p for p in PR_ODD if N(p) == 13][0]
n1, n2 = mul(P7, P7), mul(P7, P13); g = gcd_primary(n1, n2)      # gcd of supports = P7
z1 = div(n1, g)
falsifier("B13", "non-squarefree n1 = pi^2: (g, z1) = 1 fails", 1 if N(gcd_primary(g, z1)) != 1 else 0, 0)
report("B14", err14, 1e-12, f"{cnt} (n1, n2, u) triples")
# falsifier: without the indicator 1_{(u,g)=1}
ctrl = 0
for _ in range(200):
    n1 = mul(P7, random.choice([p for p in sq if N(gcd_primary(p, P7)) == 1][:30])); n2 = mul(P7, (1, 3))
    if N(gcd_primary(n1, n2)) == 1: continue
    g = gcd_primary(n1, n2); z1 = div(n1, g); z2 = div(n2, g)
    u = mul(P7, (random.randint(1, 9), random.randint(-5, 5)))
    ctrl = max(ctrl, abs(chi(n1, u).conjugate()*chi(n2, u) - chi(z1, u).conjugate()*chi(z2, u)))
falsifier("B14", "indicator dropped (u divisible by a prime of g)", ctrl, 1e-3)
bad15 = 0; ncomp = 0; err16 = 0; pairs = 0
for _ in range(120):
    z1, z2 = random.choice(sq[:40]), random.choice(sq[:40])
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 3000: continue
    q = mul(z1, z2); pairs += 1
    if q != (1, 0): err16 = max(err16, abs(abs(gauss_psi_fast(z1, z2)) - 1))
    for p in set(factor(q)):                  # local component at p: x = 1 mod q/p, x varying mod p
        qp = div(q, p); vals = set()
        for t in residues(p):
            x = add((1, 0), mul(qp, t))
            v = chi(z1, x).conjugate()*chi(z2, x)
            if v != 0: vals.add((round(v.real, 6), round(v.imag, 6)))
        ncomp += 1; bad15 += (len(vals) <= 1)
report("B15", bad15, 0, f"{ncomp} local components of psi over {pairs} coprime pairs are all nontrivial")
report("B16", err16, 1e-12, f"|gamma(psi)|=1 for {pairs} coprime pairs (N(z1 z2)<=3000)")
bad18 = sum(1 for n in sq[:60] if (chi(n, (0, 0)) != 0) != (n == (1, 0)))
report("B18", bad18, 0, "chi_n(0)=0 for n!=1 and chi_1(0)=1, so psi(0)!=0 iff z1=z2=1 (60 squarefree n)")

# ---------------- B17: exact Poisson formula
Phi = lambda t: math.exp(-math.pi*t)
hatPhi = lambda t: 2/SQ3*math.exp(-4*math.pi*t/3)
def poisson_rhs(z1, z2, g, H, drop_e0=False, wrong_conj=False, e0sq=False):
    q = mul(z1, z2); Nq = N(q); gam = gauss_psi_fast(z1, z2)
    tot = 0
    for e0 in divisors_sqfree(g):
        if drop_e0 and e0 != (1, 0): continue
        psie = chi(z1, e0).conjugate()*chi(z2, e0)
        Ne = N(e0)**2 if e0sq else N(e0)
        inner = 0
        for h in lattice(12*Ne*Nq/H + 2):
            ps = chi(z1, h).conjugate()*chi(z2, h)
            inner += (ps if wrong_conj else ps.conjugate())*hatPhi(H*N(h)/(Ne*Nq))
        tot += mu(e0)*psie/N(e0)*inner
    return H*gam/math.sqrt(Nq)*tot
P19 = [p for p in PR_ODD if N(p) == 19][0]
# psi restricted to the units must be trivial, otherwise both sides vanish identically (radial weight, unit-invariant
# exclusion): the first eight cases are nondegenerate (psi(1+omega)=1), the last two are degenerate controls (0 = 0).
cases = [((-2, -3), (1, 3), (1, 0), 30.0), ((1, 3), (1, -6), (-5, 0), 25.0), ((1, -3), (4, 3), P19, 30.0),
         ((-5, -3), (-2, 3), mul(P7, P13), 20.0), ((4, -3), (1, 0), (1, 0), 35.0), ((1, 0), (7, 3), P13, 30.0),
         ((1, 0), (1, 0), mul(P7, P13), 20.0), ((1, 0), (1, 0), (-5, 0), 15.0),
         ((7, 9), (1, 3), (1, 0), 40.0), ((-2, -3), (4, 3), (-5, 0), 25.0)]
err = 0; f1 = f2 = f3 = 0; ndeg = 0; mags = []
for z1, z2, g, H in cases:
    lhs = 0
    for u in lattice(12*H):
        if any(divides(p, u) for p in factor(g)): continue
        lhs += Phi(N(u)/H)*chi(z1, u).conjugate()*chi(z2, u)
    unit_triv = abs(chi(z1, (1, 1)).conjugate()*chi(z2, (1, 1)) - 1) < 1e-9
    ndeg += unit_triv; mags.append(f"{abs(lhs):.3g}")
    assert unit_triv == (abs(lhs) > 1e-6), (z1, z2, g, lhs)
    sc = max(1, abs(lhs))
    err = max(err, abs(lhs - poisson_rhs(z1, z2, g, H))/sc)
    if not unit_triv: continue
    if g != (1, 0): f1 = max(f1, abs(lhs - poisson_rhs(z1, z2, g, H, drop_e0=True))/sc)
    if mul(z1, z2) != (1, 0): f2 = max(f2, abs(lhs - poisson_rhs(z1, z2, g, H, wrong_conj=True))/sc)
    if g != (1, 0): f3 = max(f3, abs(lhs - poisson_rhs(z1, z2, g, H, e0sq=True))/sc)
# end to end: the smoothed CM mean square, direct vs. Poisson-dual
def Wb(y): return math.exp(-1/((y - 1)*(2 - y)))*(1 + 0.5j*y) if 1 < y < 2 else 0.0
D = 40; H = 60.0
cols = [n for n in primaries(2*D, sqfree=True) if N(n) > D]
direct = 0
for u in lattice(14*H):
    Au = sum(mu(n)*alpha_c(n)*chi(n, u)*Wb(N(n)/D) for n in cols)
    direct += Phi(N(u)/H)*abs(Au)**2
dual = 0; dual_bad = 0     # dual_bad: alpha(n1)conj(alpha(n2)) (wrong column conjugated)
for n1 in cols:
    for n2 in cols:
        g = gcd_primary(n1, n2); z1 = div(n1, g); z2 = div(n2, g)
        coef = mu(n1)*mu(n2)*alpha_c(n1).conjugate()*alpha_c(n2)*Wb(N(n1)/D).conjugate()*Wb(N(n2)/D)
        pr = poisson_rhs(z1, z2, g, H)
        dual += coef*pr
        dual_bad += coef*(alpha_c(n1)/alpha_c(n1).conjugate())*(alpha_c(n2).conjugate()/alpha_c(n2))*pr
e2e = abs(direct - dual)/abs(direct)
report("B17", max(err, e2e), 1e-9, f"{len(cases)} cases ({ndeg} nondegenerate with psi trivial on units, incl. g with 1-2 primes, inert g, z1=1, z2=1, z1=z2=1; |LHS| = {', '.join(mags)}): max rel. err {err:.1e}; "
       f"end to end CM smoothed mean square (D={D}, H={H:g}, {len(cols)} columns, complex W): direct {direct:.6f} vs dual, rel. err {e2e:.1e}")
falsifier("B17", "end to end with alpha(n1)conj(alpha(n2))", abs(direct - dual_bad)/abs(direct), 1e-3)
falsifier("B17", "terms e0 != 1 dropped", f1, 1e-3)
falsifier("B17", "psi(h) instead of conj(psi(h))", f2, 1e-3)
falsifier("B17", "N(e0)^2 in place of N(e0) in Phi^ and weight", f3, 1e-3)

# ---------------- B11: smooth majorant >= M(D,H)
def Arow(u): return sum(mu(n)*alpha_c(n)*chi(n, u)*Wb(N(n)/D) for n in cols)
rows = {u: abs(Arow(u))**2 for u in lattice(40*H)}
M = sum(v for u, v in rows.items() if 0 < N(u) <= H)
smooth = sum(math.exp(1 - N(u)/H)*v for u, v in rows.items())
report("B11", max(0.0, M - smooth), 0, f"D={D}, H={H:g}: M(D,H)={M:.4f} <= sum_u e^(1-N(u)/H)|A_u|^2 = {smooth:.4f} (u=0 included)")
bad_smooth = sum(math.exp(-N(u)/H)*v for u, v in rows.items())
falsifier("B11", "Phi(t)=e^{-t} (not >= 1 on [0,1]): M - smooth sum", M - bad_smooth, 0)

# ---------------- B19: nonzero frequencies for diagonal pairs (z1 = z2 = 1)
# With the Gaussian Phi the h != 0 terms are present; with supp Phi^ in [0, C] they need H N(h)/N(e0) <= C.
g = mul(P7, P13); H = 20.0
nz = sum(hatPhi(H*N(h)/N(e0)) for e0 in divisors_sqfree(g) for h in lattice(10*N(e0)/H + 2) if h != (0, 0))
Dsym, th, Csym, bsym = sp.symbols('D theta C b', positive=True)
cond = sp.simplify(sp.log(Dsym**(1 + th)/(Csym*bsym*Dsym))/sp.log(Dsym))   # exponent of H/(C b D)
report("B19", 0 if nz > 1e-6 else 1, 0, f"diagonal pair g=pi7*pi13, H={H:g}: nonzero-frequency mass {nz:.3g} with the Gaussian Phi",
       note="with supp Phi^ in [0,C] and N(e0)<=N(g)<=bD, a diagonal h!=0 needs H<=C b D, i.e. D^theta<=Cb for H=D^(1+theta): "
            "these terms exist only for bounded D; the hedged 'can' is correct but vacuous asymptotically")

# ---------------- B20, B21
Ne0, Ng, Nz1, Nz2, D_, H_, b_ = sp.symbols('Ne0 Ng Nz1 Nz2 D H b', positive=True)
bound = Ne0*Nz1*Nz2/H_
worst = (Ng*(b_*D_/Ng)**2/H_)                  # N(e0)<=N(g), N(z_i)<=bD/N(g)
ok_sym = sp.simplify(worst - b_**2*D_**2/(H_*Ng)) == 0
viol = 0
for _ in range(100000):
    Ngv = 10**random.uniform(0, 4); Dv = Ngv*10**random.uniform(0, 4); Hv = Dv*10**random.uniform(0, 1); bv = random.uniform(1, 4)
    Nev = Ngv*random.random(); z1v = bv*Dv/Ngv*random.random(); z2v = bv*Dv/Ngv*random.random()
    viol += Nev*z1v*z2v/Hv > bv**2*Dv**2/(Hv*Ngv)*(1 + 1e-12)
report("B20", (0 if ok_sym else 1) + viol, 0, "N(e0)N(z1z2)/H <= b^2 D^2/(H N(g)) symbolically at the extreme and on 100000 random sets")
fal = sum(1 for _ in range(1000) if (lambda Ngv: Ngv*(2*1e4/Ngv)**2/1e5 > 2**2*1e8/(1e5*Ngv**2))(10**random.uniform(0.5, 3)))
falsifier("B20", "claimed D^2/(H N(g)^2) at N(e0)=N(g)", fal, 0)
viol21 = sum(1 for Dn in range(1, 300) for Hn in range(Dn + 1, Dn + 40) if not (Fr(Dn*Dn, Hn) < Dn))
report("B21", viol21, 0, "D^2/H < D for all integers 1<=D<300, D<H<D+40 (exact)")
falsifier("B21", "H = D (not H > D)", sum(1 for Dn in range(1, 50) if Fr(Dn*Dn, Dn) >= Dn), 0)
