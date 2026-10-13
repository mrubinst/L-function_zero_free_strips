"""t02: Gauss-sum identities and the formal substitution eta -> psi (guide L1107-1186, L1242-1243, L1326-1328,
L1365-1388).

Conventions (= OAI78 L566-643 and L965-972): gamma_j(c) = N(c)^{-1/2} sum_{v mod c} chi_c(v)^j e(v/c),
e(z) = exp(4 pi i Im z/sqrt 3), G(c) = conj(chi_c(4)) gamma_3(c), alpha(c) = c/|c|, c primary squarefree coprime to 6.
All Gauss sums are computed by brute force over complete residue systems (no multiplicativity assumed).

Claims tested
  C13 L1111-1113 eq:cm-weight: gamma2(c) abar(cn^3) psi(cn^3) = gamma2(c) abar(c)^2 abar(n)^6 eta(c) eta(n)^3, psi = abar*eta.
  C14 L1157-1160 eq:signal-cm: abar(c) gamma2(c) gamma1(c) = mu(c) G(c); and the OAI78 source identities
                 gamma2^3 = mu alpha, gamma1 gamma2 = mu alpha G, |G| = 1 (OAI78 eq:signal, L965-972).
  C15 L1162      abar(c)^2 gamma2(c) eta(c) gamma1(c) = mu(c) G(c) psi(c).
  C16 L1168-1170 a_p = abar^3 psi^3, b_p = abar^2 psi^2 conj chi_p(4): b_p^3 = a_p^2 because chi_p(4)^3 = (4/p)_2 = 1.
  C17 L1174      the check of eq:signal-cm at >= 60 primes has error below 4e-15.
  C18 L1181      principal row: chi_n(1) = 1 for all n, so L(w, chi(1)) = zeta_K(w) (and denominators psi conj chi(u)).
  C19 L1193/L1227 B_m weight gamma2(c) abar(cn^3)^2 nu(cn^3) equals eq:cm-weight; original-side angular factor
                 abar(n)^2 abar(b)^6 = abar(nb^3)^2.
  C20 L1242      1_{q | A} = 1 - chi_q(A)^0 with the zero convention.
  C21 L1243      abar(c)^4 is multiplicative (complete), so the extra phase separates over primes.
  C22 L1326-1328 slot factor: conj(psi(p)) * (-psi(p) conj chi_p(u)) = -conj chi_p(u) (normalization conj psi(p_{J^c})).
  C23 L1369-1370 companion: mu(n) gamma_{-1}(n) = chi_n(-1) G(n)^{-1} abar(n) gamma2(n), checked at >= 60 primes.
  C24 L1372-1373 co-oriented classes: mu*abar <-> abar^2 gamma2 (eps=-1) and conjugates (eps=+1).
  C25 L1382-1383 anti-oriented: alpha * (abar gamma2) = gamma2, i.e. mu*alpha <-> plain gamma2 (no derivative).
  C26 L1388      chi_n(u a^6) = chi_n(u) for (a,n) = 1.
  C27 L1166-1167 squarefree local term -psi(p) conj chi_p(u) N(p)^{-x} is the p-th coefficient of 1/L(x, psi conj chi(u)).
Run time: about 30 s.
"""
import os, sys, math, cmath, random, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmcommon import *
from cmreport import claim, cited, info, untested

random.seed(1)
U0 = (2, 1)                  # eta = chi_.(2+omega): a finite-order (ray class) character with non-real values

def gammas(c):
    """dict j -> gamma_j(c) for j in (1,2,3,-1), brute force over residues mod c."""
    nn = cx(c); out = {1: 0, 2: 0, 3: 0, -1: 0}
    for x in residues(c):
        ch = chi(c, x)
        if ch == 0:
            continue
        ea = e_add(cx(x)/nn)
        for j in out:
            out[j] += ch**j*ea
    r = math.sqrt(N(c))
    return {j: v/r for j, v in out.items()}

def Gfun(c, g3):
    return chi(c, (4, 0)).conjugate()*g3

primes_split = split_primes(1000)
primes_inert = inert_primes(31)
P = primes_split + primes_inert
comps = []
sp = split_primes(200)
for a, b in itertools.combinations(sp, 2):
    if N(a)*N(b) <= 2600:
        comps.append(mul(a, b))
for a in sp[:8]:
    comps.append(mul((-5, 0), a))
comps.append(mul(mul(sp[0], sp[2]), sp[4]))
comps = [c for c in comps if mu(c) != 0][:60]
C = P + comps
GAM = {c: gammas(c) for c in C}

# ---------------------------------------------------------------- C14, C17, C23, C24, C25
e_sig = e_sig_p = e_src = e_comp = e_cls = e_anti = 0.0
f_sig = f_comp = f_comp2 = 0.0
for c in C:
    g = GAM[c]; m = mu(c); G = Gfun(c, g[3]); a = alpha(c)
    lhs = a.conjugate()*g[2]*g[1]; rhs = m*G
    e_sig = max(e_sig, abs(lhs - rhs))
    if c in P:
        e_sig_p = max(e_sig_p, abs(lhs - rhs))
    f_sig = max(f_sig, abs(a*g[2]*g[1] - rhs))                       # falsifier alpha in place of abar
    e_src = max(e_src, abs(g[2]**3 - m*a), abs(g[1]*g[2] - m*a*G), abs(abs(G) - 1))
    comp_l = m*g[-1]; comp_r = chi(c, (-1, 0))*G**-1*a.conjugate()*g[2]
    e_comp = max(e_comp, abs(comp_l - comp_r))
    f_comp = max(f_comp, abs(comp_l - G**-1*a.conjugate()*g[2]))      # falsifier: drop chi_n(-1)
    f_comp2 = max(f_comp2, abs(comp_l - chi(c, (-1, 0))*G**-1*a*g[2]))  # falsifier: alpha in place of abar
    # classes: mu*abar = abar^2 gamma2 * (gamma1/G) ; conjugate for eps=+1
    e_cls = max(e_cls, abs(m*a.conjugate() - a.conjugate()**2*g[2]*g[1]/G),
                abs(m*a - (a.conjugate()**2*g[2]*g[1]/G).conjugate()))
    # anti-oriented: mu*alpha = gamma2 * gamma1/G  (plain gamma2, no angular factor)
    e_anti = max(e_anti, abs(m*a - g[2]*g[1]/G))
npr = len(P); nco = len(comps)
claim("C14", 1158, "eq:signal-cm", "E",
      f"abar(c)gamma2(c)gamma1(c) = mu(c)G(c) at {npr} primes ({len(primes_inert)} inert) and {nco} squarefree composites",
      e_sig < 1e-12, f"max err {e_sig:.1e}; OAI78 eq:signal identities (gamma2^3=mu alpha, gamma1gamma2=mu alpha G, |G|=1) "
      f"max err {e_src:.1e}; falsifier alpha(c) in place of abar(c): {f_sig:.2f}", falsifier_failed=f_sig > 0.1 and e_src < 1e-12)
P60 = sorted(P, key=N)[:60]
e60 = max(abs(alpha(c).conjugate()*GAM[c][2]*GAM[c][1] - mu(c)*Gfun(c, GAM[c][3])) for c in P60)
claim("C17", 1174, "60 primes, err<4e-15", "E", "error of eq:signal-cm at the 60 primes of smallest norm below 4e-15 "
      "(double precision; the size of the rounding error depends on the implementation)", e60 < 4e-15,
      f"max err {e60:.2e} (largest norm {N(P60[-1])}); over all {npr} primes up to norm 1000: {e_sig_p:.2e}; "
      f"falsifier: tolerance 1e-17 not met", falsifier_failed=not (e60 < 1e-17))
claim("C23", 1369, "companion", "E",
      f"mu(n)gamma_-1(n) = chi_n(-1)G(n)^-1 abar(n)gamma2(n), {npr} primes + {nco} composites",
      e_comp < 1e-12, f"max err {e_comp:.1e}; falsifiers: drop chi_n(-1): {f_comp:.2f}; alpha for abar: {f_comp2:.2f}",
      falsifier_failed=f_comp > 0.1 and f_comp2 > 0.1)
claim("C24", 1372, "co-oriented classes", "E",
      "mu abar = abar^2 gamma2 (gamma1/G) and mu alpha = conj(abar^2 gamma2 gamma1/G)", e_cls < 1e-12,
      f"max err {e_cls:.1e}; falsifier (anti-oriented pairing alpha^2... ) below in C25", falsifier_failed=f_sig > 0.1)
claim("C25", 1383, "anti-oriented", "E", "mu alpha = gamma2 (gamma1/G): alpha*abar*gamma2 = plain gamma2", e_anti < 1e-12,
      f"max err {e_anti:.1e}; falsifier mu alpha = abar^2 gamma2 gamma1/G: "
      f"{max(abs(mu(c)*alpha(c) - alpha(c).conjugate()**2*GAM[c][2]*GAM[c][1]/Gfun(c, GAM[c][3])) for c in C):.2f}",
      falsifier_failed=max(abs(mu(c)*alpha(c) - alpha(c).conjugate()**2*GAM[c][2]*GAM[c][1]/Gfun(c, GAM[c][3])) for c in C) > 0.1)

# ---------------------------------------------------------------- C15
e15 = f15 = 0.0
for c in C:
    g = GAM[c]; eta = chi(c, U0); psi = alpha(c).conjugate()*eta
    lhs = alpha(c).conjugate()**2*g[2]*eta*g[1]; rhs = mu(c)*Gfun(c, g[3])*psi
    e15 = max(e15, abs(lhs - rhs))
    f15 = max(f15, abs(lhs - mu(c)*Gfun(c, g[3])*alpha(c)*eta))     # falsifier psi = alpha*eta
claim("C15", 1162, "weight + signal", "E", "abar(c)^2 gamma2 eta gamma1 = mu G psi, psi = abar*eta, eta = chi_.(2+omega)",
      e15 < 1e-12, f"max err {e15:.1e}; falsifier psi=alpha*eta: {f15:.2f}", falsifier_failed=f15 > 0.1)

# ---------------------------------------------------------------- C13, C19
e13 = f13 = e19 = e19b = 0.0; cnt = 0
sq = [c for c in C if N(c) < 400][:25]
ns = primaries(60)
for c in sq:
    for n in ns:
        A = mul(c, pw(n, 3))
        eta = lambda x: chi(x, U0)
        psiA = alpha(A).conjugate()*eta(A)
        lhs = GAM[c][2]*alpha(A).conjugate()*psiA
        rhs = GAM[c][2]*alpha(c).conjugate()**2*alpha(n).conjugate()**6*eta(c)*eta(n)**3
        e13 = max(e13, abs(lhs - rhs))
        f13 = max(f13, abs(lhs - GAM[c][2]*alpha(c).conjugate()**2*alpha(n).conjugate()**3*eta(c)*eta(n)**3))
        e19 = max(e19, abs(GAM[c][2]*alpha(A).conjugate()**2*eta(A) - rhs))
        e19b = max(e19b, abs(alpha(c).conjugate()**2*alpha(n).conjugate()**6 - alpha(A).conjugate()**2))
        cnt += 1
claim("C13", 1112, "eq:cm-weight", "E", f"weight identity on {cnt} pairs (c sf, n), eta = chi_.(2+omega)", e13 < 1e-12,
      f"max err {e13:.1e}; falsifier abar(n)^3: {f13:.2f}", falsifier_failed=f13 > 0.1)
claim("C19", 1193, "B_m weight; table row 'angular factor, original side'", "E",
      "gamma2(c)abar(cn^3)^2 nu(cn^3) = eq:cm-weight (nu=eta); abar(n)^2abar(b)^6 = abar(nb^3)^2", max(e19, e19b) < 1e-12,
      f"max err {max(e19, e19b):.1e}; falsifier abar(nb^3)^1 (OAI weight): see C13 falsifier ({f13:.2f})",
      falsifier_failed=f13 > 0.1)

# ---------------------------------------------------------------- C16
e16 = 0.0; f16 = 0.0; e4 = 0.0
for p in P:
    c4 = chi(p, (4, 0)); e4 = max(e4, abs(c4**3 - 1))
    eta = chi(p, U0)
    if eta == 0:
        continue
    psi = alpha(p).conjugate()*eta
    ap_ = alpha(p).conjugate()**3*psi**3
    bp_ = alpha(p).conjugate()**2*psi**2*c4.conjugate()
    e16 = max(e16, abs(bp_**3 - ap_**2))
    bpf = alpha(p).conjugate()**2*psi**2*chi(p, (2, 0)).conjugate()
    f16 = max(f16, abs(bpf**3 - ap_**2))
claim("C16", 1170, "b_p^3=a_p^2", "E", f"chi_p(4)^3 = 1 and b_p^3 = a_p^2 with psi = abar*eta at {len(P)} primes",
      max(e16, e4) < 1e-12, f"max err {max(e16, e4):.1e}; falsifier chi_p(2) in b_p: {f16:.2f}", falsifier_failed=f16 > 0.1)

# ---------------------------------------------------------------- C18
e18 = max(abs(chi(n, (1, 0)) - 1) for n in primaries(2000))
f18 = max(abs(chi(n, (1, 1)) - 1) for n in primaries(2000))     # unit u = 1+omega (nontrivial Kummer character)
claim("C18", 1181, "principal row", "E", "chi_n(1) = 1 for all primary n (N<=2000): row u=1 numerator is zeta_K", e18 == 0,
      f"max err {e18}; falsifier unit row 1+omega gives |chi-1| up to {f18:.2f}", falsifier_failed=f18 > 0.5)

# ---------------------------------------------------------------- C20
e20 = 0; f20 = 0; cnt = 0
for q in split_primes(300)[:30]:
    for A in [(random.randint(-60, 60), random.randint(-60, 60)) for _ in range(40)] + [mul(q, (3, 1)), mul(q, (1, 3))]:
        if A == (0, 0):
            continue
        ind = 1 if divides(q, A) else 0
        c = chi(q, A)
        e20 += ind != 1 - pow0(c, 0)
        f20 += ind != 1 - c**0                         # Python 0**0 == 1: forgets the zero convention
        cnt += 1
claim("C20", 1242, "marks", "E", f"1_{{q|A}} = 1 - chi_q(A)^0 (zero convention), {cnt} pairs", e20 == 0,
      f"{e20} failures; falsifier without zero convention: {f20} failures", falsifier_failed=f20 > 0)

# ---------------------------------------------------------------- C21
e21 = 0.0; f21 = 0.0
for a, b in itertools.combinations(sp[:20], 2):
    ab = mul(a, b)
    e21 = max(e21, abs(alpha(ab).conjugate()**4 - alpha(a).conjugate()**4*alpha(b).conjugate()**4))
    if mu(ab) != 0 and N(ab) <= 2600:
        ga, gb, gab = gammas(a)[2], gammas(b)[2], gammas(ab)[2]
        f21 = max(f21, abs(gab - ga*gb))
claim("C21", 1243, "phases", "E", "abar(ab)^4 = abar(a)^4 abar(b)^4", e21 < 1e-12,
      f"max err {e21:.1e}; falsifier: gamma2 is only twisted-multiplicative, |gamma2(ab)-gamma2(a)gamma2(b)| up to {f21:.2f}",
      falsifier_failed=f21 > 0.1)

# ---------------------------------------------------------------- C22
e22 = f22 = 0.0
for p in split_primes(800):
    for u in [(1, 0), (2, 1), (5, -3), (7, 0), (1, 4)]:
        eta = chi(p, U0); psi = alpha(p).conjugate()*eta; cu = chi(p, u)
        if eta == 0 or cu == 0:
            continue
        e22 = max(e22, abs(psi.conjugate()*(-psi*cu.conjugate()) - (-cu.conjugate())))
        f22 = max(f22, abs(eta.conjugate()*(-psi*cu.conjugate()) - (-cu.conjugate())))   # OAI normalization conj(eta)
claim("C22", 1327, "compensated slots", "E",
      "conj(psi(p)) times the denominator's local term -psi(p)conj chi_p(u) equals -conj chi_p(u) (u=1: -1 per prime, so "
      "A_T(Z) unchanged)", e22 < 1e-12,
      f"max err {e22:.1e}; falsifier OAI normalization conj(eta(p)) leaves -abar(p)conj chi_p(u): {f22:.2f}",
      falsifier_failed=f22 > 0.1)

# ---------------------------------------------------------------- C26
e26 = f26 = 0; cnt = 0
for n in primaries(400):
    for u in [(2, 1), (5, -3), (7, 0)]:
        for a in [(2, 0), (1, 1), (4, 3), (-5, 0), (7, 3)]:
            ua6 = mul(u, pw(a, 6))
            coprime = all(not divides(p, a) for p in factor(n))
            if coprime:
                e26 += abs(chi(n, ua6) - chi(n, u)) > 1e-12; cnt += 1
            elif chi(n, u) != 0:
                f26 += abs(chi(n, ua6) - chi(n, u)) > 1e-12
claim("C26", 1388, "sixth powers", "E", f"chi_n(u a^6) = chi_n(u) for (a,n)=1 ({cnt} triples)", e26 == 0,
      f"{e26} failures; falsifier (a,n)>1: identity fails in {f26} cases", falsifier_failed=f26 > 0)

# ---------------------------------------------------------------- C27
e27 = f27 = 0.0
x = complex(1.3, 0.4)
for p in split_primes(600):
    for u in [(2, 1), (7, 0)]:
        psi = alpha(p).conjugate()*chi(p, U0); cu = chi(p, u)
        # coefficient of N(p)^{-x} in 1/L(x, psi conj chi(u)) = prod (1 - psi(p)conj chi_p(u) N(p)^{-x})
        coeff_L = -psi*cu.conjugate()
        table = -psi*cu.conjugate()                    # OAI78 J_0 with eta -> psi, rho = chi_p(u)
        e27 = max(e27, abs(coeff_L - table))
        f27 = max(f27, abs(coeff_L - (-psi*cu)))       # falsifier: rho instead of rho^{-1}
claim("C27", 1166, "local table j=0", "E", "-psi(p)conj chi_p(u) N(p)^-x is the p-coefficient of 1/L(x,psi conj chi(u))",
      e27 < 1e-12, f"max err {e27:.1e}; falsifier rho for rho^-1: {f27:.2f}; OAI78 J_0 row = -eta(p)rho^-1 Q^-x (L4021)",
      falsifier_failed=f27 > 0.1)

# ---------------------------------------------------------------- C104 (L1124-1125)
e104 = f104 = 0.0
for c in sq:
    for n in ns[:10]:
        A = mul(c, pw(n, 3)); eta = chi(A, U0)
        w_alt = GAM[c][2]*alpha(A).conjugate()*(alpha(A)*eta.conjugate())      # target alpha*etabar
        e104 = max(e104, abs(w_alt - GAM[c][2]*eta.conjugate()))
        w_cm = GAM[c][2]*alpha(A).conjugate()*(alpha(A).conjugate()*eta)       # target abar*eta
        f104 = max(f104, abs(w_cm - GAM[c][2]*eta))
claim("C104", 1125, "target alpha*etabar", "E", "with target alpha*etabar the weight is gamma2(c) etabar(cn^3) (no angular "
      "factor); with abar*eta it carries abar^2", e104 < 1e-12,
      f"max err {e104:.1e}; falsifier 'abar*eta also has no angular factor': max dev {f104:.2f}", falsifier_failed=f104 > 0.1)
