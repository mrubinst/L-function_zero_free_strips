"""t0: the ten numbered displays of the Introduction of the CM notes, (0.1)-(0.10).

  T0.1  eq:intro-curve    E_d: y^2 = x^3 + d. Discriminant -432 d^2 and j = 0 (sympy, and PARI for 40 values of d).
  T0.2  eq:intro-cm-map   (x, y) -> (omega x, y) preserves the equation (sympy) and has order 3.
  T0.3  eq:intro-ap       a_p(d) = -sum_x ((x^3 + d)/p) for p not dividing 6d: brute-force point count, Legendre sum and
                          PARI ellap agree.  Falsifier: +sum.
  T0.4  eq:intro-L        The Dirichlet coefficients a_m of PARI's L(E_d, s) equal those of the product of the local
                          factors (0.5), built from a_p at good primes and 1 at bad primes.
  T0.5  eq:intro-local    Every bad prime of E_d has additive reduction (conductor exponent >= 2), and a_p = 0 there,
                          for 0 < |d| <= 200; some primes dividing 6d are primes of good reduction (e.g. 2 for d = 16).
  T0.6  eq:intro-fe       Lambda(s) = N^{s/2} (2 pi)^{-s} Gamma(s) L(E_d, s) satisfies Lambda(s) = eps Lambda(2 - s), with
                          eps = ellrootno and N = conductor, at complex points (PARI lfun).  Falsifier: Gamma(s/2).
  T0.7  eq:intro-fe-u     Lambda^u(s) = N^{s/2} (2 pi)^{-s} Gamma(s + 1/2) L(E_d, s + 1/2) = eps Lambda^u(1 - s).
                          Falsifier: Gamma(s).
  T0.8  eq:intro-psi      conj((4d/n)_6) = ((4d)^5/n)_6, so psi_d is the character of the notes; |psi_d(n)| = 1.
  T0.9  eq:intro-hecke    L(E_d, s + 1/2) = B_d(s) L^{(6d)}(s, psi_d), coefficientwise and exactly: with b the
                          coefficients of B_d (the local factors at p | 6d of the minimal model, from PARI),
                          a_m = sum_{m = m1 m2, m1 | (6d)^infty} b(m1) * sum_{N n = m2} n conj((4d/n)_6).
                          Falsifier: the unconjugated symbol.
  T0.10 eq:intro-inverse  sum_{N n = m} mu_K(n) n conj((4d/n)_6) = sqrt(m) q_d(m), the coefficient of
                          prod_{p not | 6d} (1 - a_p p^{-s} + p^{1-2s}).  Falsifier: without mu_K.

Exit status 1 if any check fails (t1lib).  Helpers of part_A_lines_0001_0549 (zw.py, cmcore.py) are imported
read-only.
"""
import os, sys, math, cmath, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 't1'))
from t1lib import check                                      # also puts part_A on sys.path; exit-status hook
from zw import mul, conj, N, primary, UNITS
from cmcore import gp, factor_primary, muK, sym6, legendre, primes_upto
import sympy as sp

DS = [1, 2, -3, 5, 16, -7, 17, -432, 6, 11]
M = 2500


def nums_of(label):
    return {'eq:intro-curve': '0.1', 'eq:intro-cm-map': '0.2', 'eq:intro-ap': '0.3', 'eq:intro-L': '0.4',
            'eq:intro-local': '0.5', 'eq:intro-fe': '0.6', 'eq:intro-fe-u': '0.7', 'eq:intro-psi': '0.8',
            'eq:intro-hecke': '0.9', 'eq:intro-inverse': '0.10'}[label]


# ------------------------------------------------------------------ T0.1, T0.2
x, y, d, w = sp.symbols('x y d omega')
disc = -16*(4*0**3 + 27*d**2)
out = gp("forstep(d=-20,20,1,if(d,E=ellinit([0,0,0,0,d]);print(d,\" \",E.disc,\" \",E.j)))")
bad = [l for l in out if (lambda t: int(t[1]) != -432*int(t[0])**2 or t[2] != '0')(l.split())]
check('T0.1', 'eq:intro-curve', nums_of('eq:intro-curve'), sp.expand(disc + 432*d**2) == 0 and not bad,
      f"discriminant -16(4a^3+27b^2) = -432 d^2 (sympy); PARI disc and j = 0 for 40 values of d; mismatches {bad[:3]}")
om = sp.exp(2*sp.pi*sp.I/3)
pres = sp.simplify((om*x)**3 + d - (x**3 + d)) == 0
order3 = sp.simplify(om**3 - 1) == 0 and sp.simplify(om - 1) != 0
check('T0.2', 'eq:intro-cm-map', nums_of('eq:intro-cm-map'), pres and order3,
      "(omega x)^3 + d = x^3 + d with omega^3 = 1 (sympy), so the map preserves y^2 = x^3 + d; it has order 3")

# ------------------------------------------------------------------ T0.3
P = [p for p in primes_upto(800) if p > 3]
ap_pari = {}
for dd in DS:
    for l in gp(f"E=ellinit([0,0,0,0,{dd}]);forprime(p=2,800,print(p,\" \",ellap(E,p)))"):
        p, a = map(int, l.split()); ap_pari[(dd, p)] = a
nbad = nfals = ncase = 0
for dd in DS:
    for p in P:
        if (6*dd) % p == 0:
            continue
        sq = {}
        for t in range(p):
            sq[t*t % p] = sq.get(t*t % p, 0) + 1
        aff = sum(sq.get((t**3 + dd) % p, 0) for t in range(p))
        brute = p + 1 - (aff + 1)
        leg = -sum(legendre(t**3 + dd, p) for t in range(p))
        ncase += 1
        nbad += not (brute == leg == ap_pari[(dd, p)])
        nfals += (-leg) != brute
check('T0.3', 'eq:intro-ap', nums_of('eq:intro-ap'), nbad == 0,
      f"{ncase} pairs (d, p), p < 800, p not dividing 6d: brute-force count, minus the Legendre sum and PARI ellap agree",
      falsifier=nfals > 0)

# ------------------------------------------------------------------ T0.4, T0.5
def ellan(dd, m):
    return [None] + [int(t) for t in gp(f"v=ellan(ellinit([0,0,0,0,{dd}]),{m});for(i=1,#v,print(v[i]))")]


def badprimes(dd):
    """{p: conductor exponent} for the bad primes of E_d."""
    out = gp(f"E=ellinit([0,0,0,0,{dd}]);g=ellglobalred(E);F=factor(g[1]);for(i=1,#F~,print(F[i,1],\" \",F[i,2]))")
    return {int(a): int(b) for a, b in (l.split() for l in out)}


nmis = 0; AP = {}
for dd in DS:
    an = ellan(dd, M)
    bp = badprimes(dd)
    ap = AP.setdefault(dd, {})
    for l in gp(f"E=ellinit([0,0,0,0,{dd}]);forprime(p=2,{M},print(p,\" \",ellap(E,p)))"):
        p, a = map(int, l.split()); ap[p] = a
    built = [0]*(M + 1); built[1] = 1
    # multiplicative build from the local factors
    pk = {}
    for p in primes_upto(M):
        seq = [1]
        if p in bp:
            seq += [0]*20                                          # additive reduction: local factor 1
        else:
            a1 = ap[p]; seq.append(a1)
            while p**len(seq) <= M:
                seq.append(a1*seq[-1] - p*seq[-2])
        pk[p] = seq
    for m in range(2, M + 1):
        v = 1
        f = sp.factorint(m)
        for p, e in f.items():
            v *= pk[p][e]
        built[m] = v
    nmis += sum(1 for m in range(1, M + 1) if built[m] != an[m])
check('T0.4', 'eq:intro-L', nums_of('eq:intro-L'), nmis == 0,
      f"PARI ellan equals the product of the local factors (good: 1 - a_p X + p X^2; bad: 1), m <= {M}, {len(DS)} values of d")
nadd = nmult = nzero_bad = 0; good_at_6d = []
for dd in [t for t in range(-200, 201) if t]:
    bp = badprimes(dd)
    nmult += sum(1 for e in bp.values() if e == 1)
    nadd += sum(1 for e in bp.values() if e >= 2)
    for p in (2, 3) + tuple(q for q in sp.primefactors(abs(dd)) if q > 3):
        if p not in bp:
            good_at_6d.append((dd, p))
apbad = []
for dd in DS:
    for p in badprimes(dd):
        if ap_pari.get((dd, p), 0) != 0:
            apbad.append((dd, p))
check('T0.5', 'eq:intro-local', nums_of('eq:intro-local'), nmult == 0 and not apbad and (16, 2) in good_at_6d,
      f"0<|d|<=200: {nadd} bad primes, all additive (conductor exponent >= 2), none multiplicative; a_p = 0 at bad p; "
      f"{len(good_at_6d)} primes dividing 6d have good reduction (e.g. (d, p) = (16, 2))")

# ------------------------------------------------------------------ T0.6, T0.7
fe = gp(r"""
default(realprecision,40);
pts=[0.7+3.1*I, 1.3-0.4*I, 0.2+10*I, 1.0+0.5*I];
{
foreach([1,2,-3,5,16,-7,17,6], d,
  my(E=ellinit([0,0,0,0,d]), Nd, ep, L, e1=0., e2=0., f1=0., f2=0., a, b);
  Nd=ellglobalred(E)[1]; ep=ellrootno(E); L=lfuncreate(E);
  foreach(pts, s,
    a=Nd^(s/2)*(2*Pi)^(-s)*gamma(s)*lfun(L,s); b=Nd^((2-s)/2)*(2*Pi)^(-(2-s))*gamma(2-s)*lfun(L,2-s);
    e1=max(e1, abs(a-ep*b)/abs(a));
    a=Nd^(s/2)*(2*Pi)^(-s)*gamma(s/2)*lfun(L,s); b=Nd^((2-s)/2)*(2*Pi)^(-(2-s))*gamma((2-s)/2)*lfun(L,2-s);
    f1=max(f1, abs(a-ep*b)/abs(a));
    a=Nd^(s/2)*(2*Pi)^(-s)*gamma(s+1/2)*lfun(L,s+1/2); b=Nd^((1-s)/2)*(2*Pi)^(-(1-s))*gamma(3/2-s)*lfun(L,3/2-s);
    e2=max(e2, abs(a-ep*b)/abs(a));
    a=Nd^(s/2)*(2*Pi)^(-s)*gamma(s)*lfun(L,s+1/2); b=Nd^((1-s)/2)*(2*Pi)^(-(1-s))*gamma(1-s)*lfun(L,3/2-s);
    f2=max(f2, abs(a-ep*b)/abs(a)));
  print(d,"|",Nd,"|",ep,"|",e1,"|",e2,"|",f1,"|",f2));
}
""")
rows = [[t.strip().replace(" ", "") for t in l.split("|")] for l in fe]
e1 = max(float(r[3]) for r in rows); e2 = max(float(r[4]) for r in rows)
f1 = min(float(r[5]) for r in rows); f2 = min(float(r[6]) for r in rows)
check('T0.6', 'eq:intro-fe', nums_of('eq:intro-fe'), e1 < 1e-25,
      f"Lambda(s) = eps Lambda(2-s) for 8 curves (N = {', '.join(r[1] for r in rows)}), 4 complex points: max rel err {e1:.1e}",
      falsifier=f1 > 1e-3)
check('T0.7', 'eq:intro-fe-u', nums_of('eq:intro-fe-u'), e2 < 1e-25,
      f"Lambda^u(s) = eps Lambda^u(1-s), same curves and points: max rel err {e2:.1e}", falsifier=f2 > 1e-3)

# ------------------------------------------------------------------ T0.8, T0.9, T0.10
def primary_upto(Mx):
    out = []
    R = int(2*math.isqrt(Mx)) + 3
    for A in range(-R, R + 1):
        for B in range(-R, R + 1):
            n = (A, B)
            if 0 < N(n) <= Mx and primary(n):
                out.append(n)
    return out


ELEMS = primary_upto(M)
ZETA = [cmath.exp(1j*math.pi*j/3) for j in range(6)]


def sym_idx(n, u):
    """index of (u/n)_6, None if not coprime."""
    j = 0
    for pi in factor_primary(n):
        s = sym6(u, pi)
        if s is None:
            return None
        j += s
    return j % 6


nok8 = nbad8 = 0
for dd in DS:
    u = (4*dd, 0); u5 = (4**5*dd**5, 0)
    for n in random.Random(dd).sample(ELEMS, 300):
        if N(n) % 2 == 0 or N(n) % 3 == 0:
            continue
        a, b = sym_idx(n, u), sym_idx(n, u5)
        if a is None:
            nok8 += b is None; nbad8 += b is not None
        else:
            nok8 += (-a) % 6 == b; nbad8 += (-a) % 6 != b
check('T0.8', 'eq:intro-psi', nums_of('eq:intro-psi'), nbad8 == 0,
      f"conj((4d/n)_6) = ((4d)^5/n)_6 for {nok8} pairs (d, n), exactly (symbol indices), zeros included; "
      "so psi_d = alpha chi_n(u_d) with u_d = (4d)^5, as in the notes")


def zw_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


n9 = bad9 = fals9 = 0; n10 = bad10 = fals10 = 0
for dd in DS:
    an = ellan(dd, M)
    ap = AP[dd]
    six = 6*abs(dd)
    def good(m): return math.gcd(m, six) == 1
    # c(m) = sum_{N n = m, (n, 6d) = 1} n conj((4d/n)_6), as Eisenstein integers; the same with mu_K, and falsifiers
    c = {}; cm = {}; cf = {}; cmf = {}
    for n in ELEMS:
        m = N(n)
        if not good(m):
            continue
        j = sym_idx(n, (4*dd, 0))
        if j is None:
            continue
        t = mul(n, UNITS[(-j) % 6]); tf = mul(n, UNITS[j % 6])
        c[m] = zw_add(c.get(m, (0, 0)), t); cf[m] = zw_add(cf.get(m, (0, 0)), tf)
        mu = muK(n)
        if mu:
            cm[m] = zw_add(cm.get(m, (0, 0)), (mu*t[0], mu*t[1]))
        cmf[m] = zw_add(cmf.get(m, (0, 0)), t)
    # b(m1): coefficients of the local factors at p | 6d (minimal model), from PARI's a_m at m1 | (6d)^infty
    def smooth(m):
        for p in sp.primefactors(m):
            if six % p:
                return False
        return True
    for m in range(1, M + 1):
        tot = (0, 0); totf = (0, 0)
        for m1 in sp.divisors(m):
            m2 = m // m1
            if not smooth(m1) or not good(m2):
                continue
            b = an[m1]
            v = c.get(m2, (0, 0)) if m2 > 1 else (1, 0)
            vf = cf.get(m2, (0, 0)) if m2 > 1 else (1, 0)
            tot = zw_add(tot, (b*v[0], b*v[1])); totf = zw_add(totf, (b*vf[0], b*vf[1]))
        n9 += 1
        bad9 += tot != (an[m], 0)
        fals9 += totf != (an[m], 0)
    # T0.10: unnormalized coefficients of prod_{p not | 6d} (1 - a_p X + p X^2)
    for m in range(1, M + 1):
        if not good(m):
            continue
        val = 1
        for p, e in sp.factorint(m).items():
            val *= {1: -ap[p], 2: p}.get(e, 0)
        got = cm.get(m, (0, 0)) if m > 1 else (1, 0)
        gotf = cmf.get(m, (0, 0)) if m > 1 else (1, 0)
        n10 += 1; bad10 += got != (val, 0); fals10 += gotf != (val, 0)
check('T0.9', 'eq:intro-hecke', nums_of('eq:intro-hecke'), bad9 == 0,
      f"a_m = sum b(m1) c(m2) exactly in Z[omega] for m <= {M}, {len(DS)} values of d (some non-minimal, e.g. 16, -432): "
      f"{n9} coefficients", falsifier=fals9 > 0)
check('T0.10', 'eq:intro-inverse', nums_of('eq:intro-inverse'), bad10 == 0,
      f"sum mu_K(n) n conj((4d/n)_6) over N n = m equals the coefficient of prod (1 - a_p X + p X^2), exactly, "
      f"{n10} pairs (d, m), m <= {M}", falsifier=fals10 > 0)
