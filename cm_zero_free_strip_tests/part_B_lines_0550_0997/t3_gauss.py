"""Subsections 'The arithmetic identity and the extra angular factor' and 'The next average' (lines 633-703),
and the reverse paired identity used at lines 894-895: claims B22-B33, B73.

B22  chi_n(x)^2 = (x/n)_3, the cubic symbol computed independently with exponent (N(pi)-1)/3.
B23  eq:gauss  mu_K(n) gamma_{-1}(n) = chi_n(-1) G(n)^{-1} conj(alpha(n)) gamma_2(n), 150 squarefree n (incl. composite).
B24  |G(n)| = 1.
B25  G depends only on the class of n in frak G = (primary invertible classes mod 72): two squarefree representatives
     in every one of the 432 classes (and, finer, only on n mod 4).
B26  frak G has 432 elements; its 432 characters, built from the regular representation, are homomorphisms and
     orthonormal; c_xi = |G|^{-1} sum_t G(t^{-1}) conj(xi(t)).
B27  eq:paired with the literal sum over all 432 characters xi.  Falsifiers: zeta coefficients conj(alpha) gamma_2;
     G(t) in place of G(t^{-1}); z1, z2 swapped on the right.
B28  conj(alpha(n1)) alpha(n2) = conj(alpha(z1)) alpha(z2), and CM LHS = conj(alpha(z1))alpha(z2) * zeta LHS.
B29  conj(chi_n(e0)) chi_n(h) = chi_n(h e0^5) for (e0, n) = 1.         Falsifier: chi_n(h e0).
B30  a_xi(v0 n) = a_xi(v0) a_xi(n) chi_n(v0)^4, (v0, n) = 1, two xi.   Falsifier: chi_n(v0)^2.
B31  f = e0 v0, k = e0 h: h e0^5 v0^4 = k f^4 in O, and at the level of chi_n.  Falsifier: f = e0 v0^2.
B32  eq:initial: OpenAI's change of variables b = g/e0, f = e0 v0, k = e0 h (gv = bf, N(g)/N(e0) = N(b), Phi^ argument
     H N(k)/(N(f)^2 N(m1) N(m2))), support => N(k) <= C D^2/(H N(b)^2); Sigma = XF = D/B (sympy).
B33  bookkeeping of 'E << D^eps Sigma implies (MS)': sum_{N(b)~B} (H N(b)/D^2)(XF)^2 ~ H, so M << D^{2+theta+eps}.
B73  reverse paired identity: a_xi(u1) conj(a_xi(u2)) gamma(chi_{u1} conj chi_{u2})
       = sum_{xi1} l_{xi1} b_{xi1}(u1) conj(b_{xi1}(u2)),  b_{xi1} = mu_K conj(alpha) xi1,  (xi G)(t) = sum l_{xi1} xi1(t).
     Falsifier: b = mu_K alpha xi1.
Run time: about 30 seconds.
"""
import os, sys, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmlib import *
import sympy as sp
random.seed(72)
sq = primaries(3000, sqfree=True)
sqs = [n for n in sq if N(n) <= 600]

# ---------------- B22
err = 0
for _ in range(500):
    n = random.choice(sq[:200]); x = (random.randint(-60, 60), random.randint(-60, 60))
    err = max(err, abs(chi(n, x)**2 - cubic3(n, x)))
ctrl = max(abs(chi(n, x)**4 - cubic3(n, x)) for n in sq[1:30] for x in [(2, 0), (1, 1), (3, 5)])
report("B22", err, 1e-12, "500 random (n, x): chi_n(x)^2 = (x/n)_3 (Euler criterion with exponent (N-1)/3)")
falsifier("B22", "chi_n(x)^4", ctrl, 1e-3)

# ---------------- helpers: gamma_j via the vectorized routine (checked here against brute force)
e_chk = max(abs(gamma_j_fast(j, n) - gamma_j(j, n)) for n in sq[:40] for j in (-1, 2, 3))
assert e_chk < 1e-12, e_chk
g2c = {}; gm1c = {}; Gc = {}
def g2(n):
    if n not in g2c: g2c[n] = gamma_j_fast(2, n)
    return g2c[n]
def Gn(n):
    if n not in Gc: Gc[n] = G_fast(n)
    return Gc[n]

# ---------------- B23, B24
err = 0; e24 = 0; f1 = f2 = 0
for n in sq[:150]:
    lhs = mu(n)*gamma_j_fast(-1, n)
    rhs = complex(chi(n, (-1, 0)))/Gn(n)*alpha_c(n).conjugate()*g2(n)
    err = max(err, abs(lhs - rhs)); e24 = max(e24, abs(abs(Gn(n)) - 1))
    f1 = max(f1, abs(gamma_j_fast(-1, n) - rhs))
    f2 = max(f2, abs(lhs - complex(chi(n, (-1, 0)))/Gn(n)*alpha_c(n)*g2(n)))
ncomp = sum(1 for n in sq[:150] if len(factor(n)) > 1)
report("B23", err, 1e-11, f"150 squarefree n ({ncomp} composite), N(n)<={N(sq[149])}; vectorized Gauss sums checked against brute force ({e_chk:.1e})")
falsifier("B23", "Moebius sign dropped", f1, 1e-3)
falsifier("B23", "alpha in place of conj(alpha)", f2, 1e-3)
report("B24", e24, 1e-12, "|G(n)| = 1 for 150 squarefree n")

# ---------------- B25, B26: the group frak G of primary invertible classes mod 72
M72 = 72
els = [(a, b) for a in range(M72) for b in range(M72) if a % 3 == 1 and b % 3 == 0 and not (a % 2 == 0 and b % 2 == 0)]
idx = {e: i for i, e in enumerate(els)}
def red(x): return (x[0] % M72, x[1] % M72)
def gmul(x, y): return red(mul(x, y))
order = len(els)
mt = np.array([[idx[gmul(x, y)] for y in els] for x in els])
assert all(red(mul(x, y)) in idx for x in els for y in els[:5])
one = idx[(1, 0)]
inv = np.array([int(np.where(mt[i] == one)[0][0]) for i in range(order)])
# representatives: squarefree primary n coprime to 6 in each class (two per class)
def centered(v): return v - M72 if v > M72//2 else v
reps = {}
for e in els:
    a0, b0 = centered(e[0]), centered(e[1]); cands = []
    for i in range(-3, 4):
        for j in range(-3, 4):
            n = (a0 + M72*i, b0 + M72*j)
            if N(n) < 2:
                if n == (1, 0): cands.append(n)
                continue
            cands.append(n)
    cands.sort(key=N); got = []
    for n in cands:
        if n == (1, 0) or mu(n) != 0: got.append(n)
        if len(got) == 2: break
    reps[e] = got
assert all(len(v) == 2 for v in reps.values())
Gtab = {}; bad25 = 0; spread = 0
for e, (r1, r2) in reps.items():
    v1 = Gn(r1) if r1 != (1, 0) else 1.0 + 0j; v2 = Gn(r2)
    spread = max(spread, abs(v1 - v2)); bad25 += abs(v1 - v2) > 1e-9
    Gtab[e] = v1
# finer: G depends only on n mod 4
mod4 = {}; bad4 = 0
for e, rs in reps.items():
    for r in rs:
        k = (r[0] % 4, r[1] % 4); v = Gn(r) if r != (1, 0) else 1.0 + 0j
        if k in mod4 and abs(mod4[k] - v) > 1e-9: bad4 += 1
        mod4.setdefault(k, v)
report("B25", bad25, 0, f"two squarefree representatives in each of the {len(els)} classes; max |G(r1)-G(r2)| = {spread:.1e}; "
       f"finer: G depends only on n mod 4 ({len(mod4)} classes, {bad4} violations)")
# falsifier: gamma_2 is not a class function mod 72
sp2 = max(abs(g2(r1) - g2(r2)) for e, (r1, r2) in reps.items() if r1 != (1, 0))
falsifier("B25", "gamma_2 as a class function mod 72", sp2, 1e-3)
# characters from the regular representation: P_g e_x = e_{gx}; conj(chi) are common eigenvectors
rng = np.random.default_rng(5)
gens = [idx[e] for e in random.sample(els, 6)]
# check generation
seen = {one}; frontier = [one]
while frontier:
    x = frontier.pop()
    for g in gens:
        y = mt[g, x]
        if y not in seen: seen.add(y); frontier.append(y)
assert len(seen) == order
Mreg = np.zeros((order, order), complex)
for g in gens:
    c = rng.normal() + 1j*rng.normal()
    for x in range(order): Mreg[mt[g, x], x] += c
w, V = np.linalg.eig(Mreg)
chars = []
for k in range(order):
    v = V[:, k]; chars.append(np.conj(v/v[one]))
chars = np.array(chars)                                       # chars[k, x] = xi_k(x)
# snap to roots of unity of order dividing |G|
ang = np.angle(chars)/(2*np.pi)*order
snapped = np.exp(2j*np.pi*np.round(ang)/order)
snap_err = np.max(np.abs(snapped - chars)); chars = snapped
hom = max(np.max(np.abs(chars[:, mt[i, :]] - chars[:, [i]]*chars)) for i in range(order))
gram = chars @ chars.conj().T/order
orth = np.max(np.abs(gram - np.eye(order)))
report("B26", max(hom, orth, snap_err), 1e-8, f"|frak G| = {order}; {order} characters: homomorphism error {hom:.1e}, orthonormality {orth:.1e}")
Gvec = np.array([Gtab[e] for e in els])
c_xi = (chars.conj() @ Gvec[inv])/order                     # c_xi = |G|^{-1} sum_t G(t^{-1}) conj(xi(t))
c_bad = (chars.conj() @ Gvec)/order                          # falsifier: G(t) instead of G(t^{-1})
def cls(n): return idx[red(n)]

# ---------------- B27, B28, B73
pairs = []
cand = [n for n in sq if N(n) <= 300]
while len(pairs) < 70:
    z1, z2 = random.choice(cand), random.choice(cand)
    if N(gcd_primary(z1, z2)) > 1 or N(mul(z1, z2)) > 4000: continue
    pairs.append((z1, z2))
pairs += [((1, 0), (7, 3)), ((-2, -3), (1, 0)), ((1, 0), (1, 0))]
e27 = 0; fz = 0; fG = 0; fsw = 0; e28 = 0; e73 = 0; f73 = 0
# reverse: (xi G)(t) = sum_{xi1} l_{xi1} xi1(t) for one fixed nontrivial xi
xi0 = 7                                                       # an arbitrary character index
xiG = chars[xi0]*Gvec
l_xi1 = (chars.conj() @ xiG)/order
for z1, z2 in pairs:
    gam = gauss_psi_fast(z1, z2)
    a1 = alpha_c(z1).conjugate()**2*g2(z1) if z1 != (1, 0) else 1.0
    a2 = alpha_c(z2).conjugate()**2*g2(z2) if z2 != (1, 0) else 1.0
    m1 = mu(z1) if z1 != (1, 0) else 1; m2 = mu(z2) if z2 != (1, 0) else 1
    lhs = m1*m2*alpha_c(z1).conjugate()*alpha_c(z2)*gam
    X1 = chars[:, cls(z1)]; X2 = chars[:, cls(z2)]
    rhs = np.sum(c_xi*a1*X1*np.conj(a2*X2))
    e27 = max(e27, abs(lhs - rhs))
    za1 = alpha_c(z1).conjugate()*g2(z1) if z1 != (1, 0) else 1.0
    za2 = alpha_c(z2).conjugate()*g2(z2) if z2 != (1, 0) else 1.0
    fz = max(fz, abs(lhs - np.sum(c_xi*za1*X1*np.conj(za2*X2))))
    fG = max(fG, abs(lhs - np.sum(c_bad*a1*X1*np.conj(a2*X2))))
    fsw = max(fsw, abs(lhs - np.sum(c_xi*a2*X2*np.conj(a1*X1))))
    # B28: zeta LHS times conj(alpha(z1))alpha(z2); and conj(alpha(g z1))alpha(g z2) = conj(alpha(z1))alpha(z2)
    zeta_lhs = m1*m2*gam
    g = random.choice([p for p in sq[1:30] if N(gcd_primary(p, mul(z1, z2))) == 1])
    e28 = max(e28, abs(lhs - alpha_c(z1).conjugate()*alpha_c(z2)*zeta_lhs),
              abs(alpha_c(mul(g, z1)).conjugate()*alpha_c(mul(g, z2)) - alpha_c(z1).conjugate()*alpha_c(z2)))
    # B73 (u1, u2) = (z1, z2): a_xi(u1) conj(a_xi(u2)) gamma(chi_{u1} conj chi_{u2}) with xi = chars[xi0]
    gam_r = gauss_psi_general(z1, z2, conj_first=False)
    L73 = a1*chars[xi0, cls(z1)]*np.conj(a2*chars[xi0, cls(z2)])*gam_r
    b1 = m1*alpha_c(z1).conjugate(); b2 = m2*alpha_c(z2).conjugate()
    R73 = np.sum(l_xi1*b1*X1*np.conj(b2*X2))
    e73 = max(e73, abs(L73 - R73))
    bb1 = m1*alpha_c(z1); bb2 = m2*alpha_c(z2)
    f73 = max(f73, abs(L73 - np.sum(l_xi1*bb1*X1*np.conj(bb2*X2))))
report("B27", e27, 1e-10, f"{len(pairs)} coprime squarefree pairs (incl. z=1), literal sum over the {order} characters xi")
falsifier("B27", "zeta coefficients conj(alpha) gamma_2 xi", fz, 1e-3)
falsifier("B27", "c_xi with G(t) instead of G(t^{-1})", fG, 1e-3)
falsifier("B27", "z1, z2 swapped on the right", fsw, 1e-3)
report("B28", e28, 1e-10, "CM LHS = conj(alpha(z1))alpha(z2) x zeta LHS; conj(alpha(g z1))alpha(g z2) = conj(alpha(z1))alpha(z2)")
report("B73", e73, 1e-10, f"{len(pairs)} pairs, xi = character #{xi0}: Moebius coefficients mu_K conj(alpha) xi1 after the reverse identity")
falsifier("B73", "b = mu_K alpha xi1", f73, 1e-3)

# ---------------- B29, B30, B31
err29 = 0; f29 = 0
for _ in range(300):
    n = random.choice(sq[:150]); e0 = random.choice(sq[:40]); h = (random.randint(-30, 30), random.randint(-30, 30))
    if N(gcd_primary(n, e0)) > 1: continue
    err29 = max(err29, abs(chi(n, e0).conjugate()*chi(n, h) - chi(n, mul(h, pw(e0, 5)))))
    f29 = max(f29, abs(chi(n, e0).conjugate()*chi(n, h) - chi(n, mul(h, e0))))
report("B29", err29, 1e-12, "300 random (n, e0, h) with (e0, n) = 1")
falsifier("B29", "chi_n(h e0)", f29, 1e-3)
xi_alt = lambda n: chars[11, cls(n)]
err30 = 0; f30 = 0; cnt = 0
for _ in range(200):
    v0, n = random.choice(sq[:60]), random.choice(sq[:60])
    if v0 == (1, 0) or n == (1, 0) or N(gcd_primary(v0, n)) > 1 or N(mul(v0, n)) > 6000: continue
    for xi in (lambda m: 1.0, xi_alt):
        a = lambda m: alpha_c(m).conjugate()**2*g2(m)*xi(m)
        err30 = max(err30, abs(a(mul(v0, n)) - a(v0)*a(n)*chi(n, v0)**4))
        f30 = max(f30, abs(a(mul(v0, n)) - a(v0)*a(n)*chi(n, v0)**2))
    cnt += 1
report("B30", err30, 1e-10, f"{cnt} coprime pairs (v0, n), trivial and nontrivial xi")
falsifier("B30", "chi_n(v0)^2", f30, 1e-3)
err31 = 0; e31c = 0; f31 = 0
for _ in range(300):
    e0, v0 = random.choice(sq[:30]), random.choice(sq[:30]); h = (random.randint(-20, 20), random.randint(-20, 20))
    f = mul(e0, v0); k = mul(e0, h)
    err31 += mul(mul(h, pw(e0, 5)), pw(v0, 4)) != mul(k, pw(f, 4))
    n = random.choice(sq[:150])
    e31c = max(e31c, abs(chi(n, mul(h, pw(e0, 5)))*chi(n, v0)**4 - chi(n, k)*chi(n, f)**4))
    fb = mul(e0, pw(v0, 2))
    f31 = max(f31, abs(chi(n, mul(h, pw(e0, 5)))*chi(n, v0)**4 - chi(n, k)*chi(n, fb)**4))
report("B31", err31 + e31c, 1e-12, "300 random (e0, v0, h): exact identity in O and at the level of chi_n")
falsifier("B31", "f = e0 v0^2", f31, 1e-3)

# ---------------- B32: OpenAI's change of variables behind eq:initial, and the scale algebra
bad = 0
for _ in range(300):
    ps = random.sample(sq[1:60], 6)
    ps = [p for p in ps if len(factor(p)) == 1][:4] or [sq[1]]
    pr = sorted(set(itertools.chain.from_iterable(factor(p) for p in random.sample(sq[1:80], 4))))
    if len(pr) < 3: continue
    random.shuffle(pr)
    e0 = pr[0]; bq = pr[1]; v0 = pr[2]                        # g = b e0, f = e0 v0, (g, v0) = 1
    g = mul(bq, e0); f = mul(e0, v0); h = (random.randint(-9, 9), random.randint(-9, 9))
    if h == (0, 0): continue
    k = mul(e0, h)
    bad += mul(g, v0) != mul(bq, f)                           # g v = b f
    bad += N(g)//N(e0) != N(bq)                               # N(g)/N(e0) = N(b)
    m1, m2 = random.choice(sq[1:50]), random.choice(sq[1:50])
    lhs = N(h)*1.0/(N(e0)*N(v0)**2*N(m1)*N(m2)); rhs = N(k)*1.0/(N(f)**2*N(m1)*N(m2))
    bad += abs(lhs - rhs) > 1e-15*rhs
    bad += mu(f) != mu(e0)*mu(v0)
Dd, Hh, Bb, Ff, b0, CP, Nb, Nf, Nm1, Nm2, Nk = sp.symbols('D H B F b0 C_Phi Nb Nf Nm1 Nm2 Nk', positive=True)
# support: N(b)N(f)N(m_j) <= b0 D and H N(k)/(N(f)^2 N(m1)N(m2)) <= C_Phi  =>  N(k) <= C_Phi b0^2 D^2/(H N(b)^2)
kmax = CP*Nf**2*(b0*Dd/(Nb*Nf))**2/Hh
ok1 = sp.simplify(kmax - CP*b0**2*Dd**2/(Hh*Nb**2)) == 0
X = Dd/(Bb*Ff); ok2 = sp.simplify(X*Ff - Dd/Bb) == 0
report("B32", bad + (0 if ok1 and ok2 else 1), 0, "gv=bf, N(g)/N(e0)=N(b), mu(f)=mu(e0)mu(v0), Phi^ argument -> H N(k)/(N(f)^2N(m1)N(m2)) "
       "on random prime data; N(k) <= C b0^2 D^2/(H N(b)^2) and Sigma = XF = D/B (sympy)")
falsifier("B32", "rows << D^2/(H B) (one power of B)", 0 if sp.simplify(kmax - CP*b0**2*Dd**2/(Hh*Nb)) == 0 else 1, 0)

# ---------------- B33: bookkeeping
th = sp.symbols('theta', positive=True)
SB = Bb*(Hh*Bb/Dd**2)*(Dd/Bb)**2                              # (#b ~ B) x (H N(b)/D^2)(XF)^2 with XF = D/B
ok = sp.simplify(SB - Hh) == 0
Mexp = sp.simplify(sp.expand_log(sp.log(Dd*Hh.subs(Hh, Dd**(1 + th)))/sp.log(Dd), force=True))
report("B33", 0 if ok and sp.simplify(Mexp - (2 + th)) == 0 else 1, 0,
       "sum_{N(b)~B}(H N(b)/D^2)(XF)^2 = H, so M_D << H D^eps and M = D M_D << D^(2+theta+eps) (sympy); "
       "the smooth-kernel separation and the dyadic/log losses are not tested")
