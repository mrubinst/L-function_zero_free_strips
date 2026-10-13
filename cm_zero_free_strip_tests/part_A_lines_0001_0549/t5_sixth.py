"""Part A, Section 'Sixth powers isolate the chosen elliptic curve' (lines 457-549)."""
from cmcore import *
from fractions import Fraction as F
import sympy as sp
random.seed(5)

def Rcoef(u, m, bn):
    s = (0, 0)
    for n in bn.get(m, []):
        mu_ = muK(n)
        if mu_ == 0:
            continue
        j = chi_idx(n, u)
        if j is None:
            continue
        t = mul(UNITS[j], n)
        s = add(s, t if mu_ == 1 else neg(t))
    return s

def PY(Y, d):
    """P_Y = primary primes pi !| 6d with Y/2 < N(pi) <= Y."""
    out = []
    for p in primes_upto(int(Y)):
        if p < 5 or d % p == 0 or p % 3 == 2:
            continue                       # inert primes are added below (norm p^2)
        for pi in primes_above(p):
            if Y/2 < N(pi) <= Y:
                out.append(pi)
    q = 5
    while q*q <= Y:
        if is_prime(q) and q % 3 == 2 and d % q != 0 and Y/2 < q*q <= Y:
            out.append((-q, 0))
        q += 1
    return out

def A_rows(D, u_list, cols):
    """A_u(D) for each u in u_list (vectorised; coordinates must fit in int64)."""
    v = np.array([muK(n)*alpha(n) for n in cols])*W_test(np.array([N(n)/D for n in cols]))
    keep = np.abs(v) > 0
    cols = [c_ for c_, k in zip(cols, keep) if k]; v = v[keep]
    assert all(abs(u[0]) < 2**62 and abs(u[1]) < 2**62 for u in u_list)
    C = row_matrix(np.array([u[0] for u in u_list], dtype=np.int64), np.array([u[1] for u in u_list], dtype=np.int64), cols)
    return list(C @ v)

def run():
    # A5.01  l.459: N_0 = N((4d)^5) = |4d|^10
    ok = all(N(((4*d)**5, 0)) == abs(4*d)**10 for d in range(-200, 201) if d)
    report('A5.01', '459', '', ok, 'N((4d)^5) = |4d|^10 for 0<|d|<=200',
           falsifier=not all(N(((4*d)**5, 0)) == abs(4*d)**5 for d in range(-5, 6) if d))

    # A5.02  l.460-464: chi_n(u_0 v^6) = chi_n(u_0) 1_{(n,v)=1}, u_0 = u_d
    cols = columns(4000)
    bad = 0; fb = 0
    for _ in range(3000):
        d = random.choice([1, 2, -3, 7, 11]); u0 = ((4*d)**5, 0)
        n = random.choice(cols)
        v = (random.randint(-15, 15), random.randint(-15, 15))
        if v == (0, 0):
            continue
        cop = all(not divides(pi, v) for pi in set(factor_primary(n)))
        j0 = chi_idx(n, u0)
        expect = j0 if cop else None
        bad += chi_idx(n, mul(u0, pw(v, 6))) != expect
        fb += chi_idx(n, mul(u0, pw(v, 2))) != expect
    report('A5.02', '460-464', '', bad == 0, '3000 random (d,n,v): chi_n(u_d v^6) = chi_n(u_d) 1_(n,v)=1',
           falsifier=fb > 0)

    # A5.03  l.472-473: J = |P_Y| ~ Y/log Y (PNT) -- sanity; exact: J = 0 for Y < 7
    msg = []
    for d in [1, 7*13*19*31]:
        for Y in [10**3, 10**4, 10**5, 10**6]:
            J = len(PY(Y, d)); msg.append(f'd={d},Y={Y}: J/(Y/logY)={J/(Y/math.log(Y)):.3f}')
    J_small = [len(PY(Y, 1)) for Y in [2, 3, 4, 5, 6, 6.99]]
    info('A5.03', '472-473', '', '; '.join(msg) + f' (limit 1/2); J for Y in [2,7): {J_small}')

    # A5.04  l.473: rows u_0 pi^6 distinct with norm <= H = N_0 Y^6
    bad = 0; dup = 0
    for d in [1, 7]:
        u0 = ((4*d)**5, 0); N0 = abs(4*d)**10; Y = 3000; H = N0*Y**6
        P = PY(Y, d)
        rows = [mul(u0, pw(pi, 6)) for pi in P]
        bad += sum(N(r) > H for r in rows) + (len(set(rows)) != len(rows))
        rows_ass = [mul(u0, pw(mul(z, pi), 6)) for pi in P[:50] for z in UNITS]
        dup += len(rows_ass) - len(set(rows_ass))
    report('A5.04', '473', '', bad == 0, 'Y=3000, d=1,7: all rows u_0 pi^6 (pi in P_Y) distinct, N <= N_0 Y^6 = H',
           falsifier=dup > 0)

    # A5.05  l.474-478: A_{u0} - A_{u0 pi^6} = sum_{pi | n} mu alpha chi_n(u0) W
    D = 3000; cols = [n for n in columns(2*D) if muK(n) != 0]
    err = 0.0; fbe = 0.0
    for d in [1, 7]:
        u0 = ((4*d)**5, 0)
        P = PY(60, d)
        Au0 = A_rows(D, [u0], cols)[0]
        Ap = A_rows(D, [mul(u0, pw(pi, 6)) for pi in P], cols)
        for pi, Av in zip(P, Ap):
            rhs = sum(muK(n)*alpha(n)*chi(n, u0)*W_test(np.array([N(n)/D]))[0] for n in cols if divides(pi, n))
            wrong = sum(muK(n)*alpha(n)*chi(n, u0)*W_test(np.array([N(n)/D]))[0] for n in cols if not divides(pi, n))
            err = max(err, abs(Au0 - Av - rhs)); fbe = max(fbe, abs(Au0 - Av - wrong))
    report('A5.05', '474-478', '', err < 1e-10, f'D={D}, Y=60, d=1,7: max err {err:.1e}', falsifier=fbe > 1e-3)

    # A5.06-07  l.479-484, eq:approx: #multiples of pi = O(D/N(pi)); coefficients |mu alpha chi| <= 1;
    # |A_{u0}-A_{u0 pi^6}| <= max|W| * #multiples  (<< D/Y)
    worst_cnt = 0.0; worst_diff = 0.0; tri = 0; maxcoef = 0.0; tri_wrong = 0
    Wmax = float(np.max(np.abs(W_test(np.linspace(1, 2, 20001)))))
    for D in [1000, 3000, 8000]:
        cols = [n for n in columns(2*D) if muK(n) != 0]
        for d in [1, 7]:
            u0 = ((4*d)**5, 0)
            for Y in [20, 60, 200]:
                P = PY(Y, d)
                Au0 = A_rows(D, [u0], cols)[0]
                Ap = A_rows(D, [mul(u0, pw(pi, 6)) for pi in P], cols)
                for pi, Av in zip(P, Ap):
                    mult = sum(1 for n in cols if divides(pi, n) and abs(W_test(np.array([N(n)/D]))[0]) > 0)
                    worst_cnt = max(worst_cnt, mult/(D/N(pi)))
                    worst_diff = max(worst_diff, abs(Au0 - Av)/(D/Y))
                    tri += abs(Au0 - Av) > Wmax*mult + 1e-9
                    mult2 = sum(1 for n in cols if divides(pw(pi, 2), n))       # wrong: multiples of pi^2
                    tri_wrong += abs(Au0 - Av) > Wmax*mult2 + 1e-9
        maxcoef = max(maxcoef, max(abs(muK(n)*alpha(n)*chi(n, ((4*7)**5, 0))) for n in cols))
    report('A5.06', '479-481', '', tri == 0 and maxcoef <= 1 + 1e-12,
           f'|A_u0-A_u0pi^6| <= max|W| * #(multiples of pi) in all cases; max coefficient modulus {maxcoef:.3f}',
           falsifier=tri_wrong > 0, kind='I')
    info('A5.07', '479-484', 'eq:approx', f'max #multiples/(D/N(pi)) = {worst_cnt:.3f}; '
         f'max |A_u0-A_u0pi^6|/(D/Y) = {worst_diff:.4f} (D in 1e3..8e3, Y in 20..200)')

    # A5.08  l.485: |z|^2 <= 2|w|^2 + 2|z-w|^2  (difference = |z-2w|^2)
    a, b, c, e = sp.symbols('a b c e', real=True)
    z = a + sp.I*b; w = c + sp.I*e
    absq = lambda t: sp.expand(t*sp.conjugate(t))
    ok = sp.simplify(2*absq(w) + 2*absq(z - w) - absq(z) - absq(z - 2*w)) == 0
    rz = np.random.default_rng(1).normal(size=(10**6, 4))
    Z = rz[:, 0] + 1j*rz[:, 1]; Wc = rz[:, 2] + 1j*rz[:, 3]
    ok = ok and np.all(np.abs(Z)**2 <= 2*np.abs(Wc)**2 + 2*np.abs(Z - Wc)**2 + 1e-12)
    fb = np.any(np.abs(Z)**2 > np.abs(Wc)**2 + np.abs(Z - Wc)**2 + 1e-12)
    report('A5.08', '485', '', bool(ok), 'sympy: 2|w|^2+2|z-w|^2-|z|^2 = |z-2w|^2; 1e6 random pairs', kind='I',
           falsifier=bool(fb))

    # A5.09  l.485-492, eq:extraction: toy instance with the actual mean square M(D,H), H = N_0 Y^6
    lines = []; okall = True; noerr_fails = 0
    for (u0, N0, Y, D) in [((1, 0), 1, 8, 30), ((2, 0), 4, 8, 30), ((1, 0), 1, 8, 60)]:
        H = N0*Y**6
        P = PY(Y, 16 if u0 == (1, 0) else 1)
        cols = [n for n in columns(2*D) if muK(n) != 0]
        v = np.array([muK(n)*alpha(n) for n in cols])*W_test(np.array([N(n)/D for n in cols]))
        keep = np.abs(v) > 0; cols = [c_ for c_, k in zip(cols, keep) if k]; v = v[keep]
        xs, ys = lattice_points(H)
        M = 0.0
        for i in range(0, len(xs), 100000):
            Cm = row_matrix(xs[i:i + 100000], ys[i:i + 100000], cols)
            M += float(np.sum(np.abs(Cm @ v)**2))
        S = sum(v[k]*chi(n, u0) for k, n in enumerate(cols))
        Ap = [sum(v[k]*chi(n, mul(u0, pw(pi, 6))) for k, n in enumerate(cols)) for pi in P]
        J = len(P)
        lhs = J*abs(S)**2
        sel = sum(abs(a_)**2 for a_ in Ap)
        err2 = 2*sum(abs(S - a_)**2 for a_ in Ap)
        ok1 = lhs <= 2*sel + err2 + 1e-9 and sel <= M + 1e-9 and all(N(mul(u0, pw(pi, 6))) <= H for pi in P)
        okall &= ok1
        noerr_fails += lhs > 2*sel + 1e-9           # wrong variant without the O(J D^2/Y^2) term
        lines.append(f'u0={u0},Y={Y},H={H},D={D}: J={J}, J|S|^2={lhs:.3f} <= 2*sum_sel={2*sel:.3f} + err={err2:.3f}; '
                     f'sum_sel={sel:.3f} <= M={M:.1f}; err/(J D^2/Y^2)={err2/(J*D*D/Y**2):.4f}')
    report('A5.09', '485-492', 'eq:extraction', okall, ' | '.join(lines) + f' | variant without error term fails in '
           f'{noerr_fails} case(s)', kind='A', falsifier=noerr_fails > 0)

    # A5.10  l.494-502, eq:exponents
    D_, th, ep, N0 = sp.symbols('D theta epsilon N_0', positive=True)
    H = D_**(1 + th)
    t1 = sp.expand_log(sp.log(D_*H**sp.Rational(5, 6)*D_**ep), force=True)/sp.log(D_)
    t2 = sp.expand_log(sp.log(D_**2*H**sp.Rational(-1, 3)), force=True)/sp.log(D_)
    # derivation: 2M/J with M = D^{2+th+ep}, J ~ Y = H^{1/6}; D^2/Y^2
    t1d = sp.simplify((2 + th + ep) - (1 + th)/6); t2d = sp.simplify(2 - (1 + th)/3)
    ok = (sp.simplify(t1 - (sp.Rational(11, 6) + 5*th/6 + ep)) == 0 and sp.simplify(t2 - (sp.Rational(5, 3) - th/3)) == 0
          and sp.simplify(t1d - t1) == 0 and sp.simplify(t2d - t2) == 0)
    fal = sp.simplify(t1 - (sp.Rational(11, 6) + th/6 + ep)) != 0 and sp.simplify(t2 - (sp.Rational(5, 3) + th/3)) != 0
    report('A5.10', '494-502', 'eq:exponents', ok, 'D H^(5/6) D^eps = D^(11/6+5th/6+eps), D^2 H^(-1/3) = D^(5/3-th/3); '
           'both agree with 2M/J (M=D^(2+th+eps), J~H^(1/6)) and D^2/Y^2', kind='I', falsifier=fal)

    # A5.11  l.503-507, eq:Sbound: for every eta>0 there are 0<theta<=1/10, eps>0 with max exponent <= 11/12+eta
    bad = 0; fb = 0
    for k in range(1, 400):
        eta = F(k, 400)**3
        theta = min(F(1, 10), eta/2); eps = eta
        e1 = F(11, 12) + 5*theta/12 + eps/2; e2 = F(5, 6) - theta/6
        bad += not (max(e1, e2) <= F(11, 12) + eta and 0 < theta <= F(1, 10))
        fb += not (max(e1, e2) <= F(5, 6) + eta)
    report('A5.11', '503-507', 'eq:Sbound', bad == 0, '399 values of eta in (0,1): theta=min(1/10,eta/2), eps=eta give '
           'exponent <= 11/12+eta (exact Fractions); second term 5/6-theta/6 < 11/12', kind='I', falsifier=fb > 0)

    # A5.12  l.508-513: N_0 dependence: S << N0^{1/12} D^{11/12+5th/12+eps} + N0^{1/6} D^{5/6-th/6}; N0^{1/12}=|4d|^{5/6}
    Y = (H/N0)**sp.Rational(1, 6)
    sq1 = sp.powsimp(sp.sqrt(D_**(2 + th + ep)/Y), force=True)
    sq2 = sp.powsimp(sp.sqrt(D_**2/Y**2), force=True)
    c1 = sp.simplify(sp.log(sq1/(N0**sp.Rational(1, 12)*D_**(sp.Rational(11, 12) + 5*th/12 + ep/2))))
    c2 = sp.simplify(sp.log(sq2/(N0**sp.Rational(1, 6)*D_**(sp.Rational(5, 6) - th/6))))
    dd = sp.Symbol('d', positive=True)
    c3 = sp.simplify(((4*dd)**10)**sp.Rational(1, 12) - (4*dd)**sp.Rational(5, 6))
    ok = sp.simplify(sp.expand_log(c1, force=True)) == 0 and sp.simplify(sp.expand_log(c2, force=True)) == 0 and c3 == 0
    cw = sp.simplify(sp.expand_log(sp.log(sq1/(N0**sp.Rational(1, 6)*D_**(sp.Rational(11, 12) + 5*th/12 + ep/2))),
                                   force=True))
    report('A5.12', '508-513', '', ok, 'sqrt(D^(2+th+eps)/Y) = N0^(1/12) D^(11/12+5th/12+eps/2), sqrt(D^2/Y^2) = '
           'N0^(1/6) D^(5/6-th/6) with Y=(H/N0)^(1/6); (|4d|^10)^(1/12) = |4d|^(5/6)', kind='I', falsifier=cw != 0)

    # A5.12b remark: J = 0 when Y < 7, so the 'valid once Y >= 2' range relies on the trivial bound there.
    # |q_d(m)| <= 2^omega(m) (Hasse), i.e. Q(m)^2 <= 4^omega(m) m; and for Y < 7 the first term exceeds D^(1+th/2)/sqrt7.
    an = ellan(1, 20000)
    bad = 0
    for m in range(1, 20001):
        Q = Qd_int(1, m, lambda p: an[p])
        bad += Q*Q > 4**len(factorint(m))*m
    e_first_smallY = sp.simplify((1 + th)/12 + sp.Rational(11, 12) + 5*th/12 - (1 + th/2))
    report('A5.12b', '514', 'remark', bad == 0 and J_small == [0]*6 and e_first_smallY == 0,
           f'J=0 for 2<=Y<7 (so the extraction is empty there); but then N0>H/7^6 and the first term is '
           f'>= D^(1+th/2+eps)/sqrt7, which exceeds the trivial |S_d(D)| << D log D (|q_d(m)|<=2^omega(m) checked '
           f'm<=20000): stated bound still true', kind='I', falsifier=J_small[0] == 0)   # 'J>=1 once Y>=2' is false

    # A5.13  l.514: D^{1+theta} >= 2^6 N_0  <=>  Y >= 2
    bad = 0; fb = 0
    for Hn in range(1, 3000):
        for N0n in [1, 2, 7, 33, 64]:
            Yge2 = F(Hn, N0n) >= 2**6        # Y^6 = H/N0
            bad += Yge2 != (Hn >= 2**6*N0n)
            fb += Yge2 != (Hn >= 2**3*N0n)
    report('A5.13', '514', '', bad == 0, 'Y=(H/N0)^(1/6) >= 2 <=> H >= 2^6 N0 (integer grid)', kind='I', falsifier=fb > 0)

    # A5.14  l.514-515: nontrivial (first term < D) only if D^{1-5theta} > N_0 = |4d|^10, i.e. D larger than ~|d|^10
    lhs = sp.simplify(sp.expand_log(sp.log(N0**sp.Rational(1, 12)*D_**(sp.Rational(11, 12) + 5*th/12)/D_), force=True))
    target = sp.expand_log(sp.log(N0)/12 - (1 - 5*th)*sp.log(D_)/12, force=True)
    ok = sp.simplify(lhs - target) == 0
    report('A5.14', '514-515', '', ok, 'log(N0^(1/12)D^(11/12+5th/12)/D) = (log N0 - (1-5th) log D)/12, so the bound '
           'beats D only when D^(1-5th) > |4d|^10 (>= |d|^10)', kind='I',
           falsifier=sp.simplify(lhs - sp.expand_log(sp.log(N0)/12 - (1 + 5*th)*sp.log(D_)/12, force=True)) != 0)

    # A5.15  l.515-516: E_{dt^6} isomorphic to E_d over Q
    x, y, d_, t = sp.symbols('x y d t')
    sub_ok = sp.expand((t**3*y)**2 - (t**2*x)**3 - d_*t**6) == sp.expand(t**6*(y**2 - x**3 - d_))
    script = ("for(i=1,4,d=[1,2,-3,7][i];for(j=1,3,t=[2,3,5][j];"
              "print(ellminimalmodel(ellinit([0,0,0,0,d*t^6]))[1..5]==ellminimalmodel(ellinit([0,0,0,0,d]))[1..5],"
              "\" \",ellminimalmodel(ellinit([0,0,0,0,d*t^2]))[1..5]==ellminimalmodel(ellinit([0,0,0,0,d]))[1..5])))")
    out = [l.split() for l in gp(script)]
    ok = sub_ok and len(out) == 12 and all(o[0] == '1' for o in out)
    report('A5.15', '515-516', '', ok, '(x,y)=(t^2X,t^3Y) maps y^2=x^3+dt^6 to Y^2=X^3+d; PARI minimal models equal '
           '(d in 1,2,-3,7; t in 2,3,5)', falsifier=all(o[1] == '0' for o in out))

    # A5.16-17  l.516-519: d sixth-power-free => f_p(E_d)=2 for p>=5, p|d; |d| << N^{5/2}; |d|^10 << N^25
    ds = [dv for dv in range(-3000, 3001) if dv and all(e < 6 for e in factorint(abs(dv)).values())]
    ds += [p**5 for p in primes_upto(400) if p >= 5] + [-(p**5)*q**4 for p, q in [(5, 7), (11, 13), (17, 19)]]
    script = ("L=[" + ",".join(map(str, ds)) + "];for(i=1,#L,d=L[i];E=ellinit([0,0,0,0,d]);N=ellglobalred(E)[1];"
              "f=factor(abs(d))[,1];s=\"\";for(j=1,#f,p=f[j];if(p>=5,s=Str(s,\" \",p,\":\",valuation(N,p))));"
              "print(d,\" \",N,s))")
    bad = 0; fb = 0; worst = 0.0; worst2 = 0.0; nprimes = 0
    for line in gp(script, timeout=1200):
        parts = line.split(); dv, Nn = int(parts[0]), int(parts[1])
        for pe in parts[2:]:
            p, v = map(int, pe.split(':')); nprimes += 1
            bad += v != 2; fb += v != 1
        bad += dv*dv > 6**10*Nn**5
        worst = max(worst, abs(dv)/Nn**2.5); worst2 = max(worst2, abs(dv)/Nn**2)
    report('A5.16', '516-519', '', bad == 0,
           f'{len(ds)} sixth-power-free d (|d|<=3000, p^5, p^5q^4): f_p=2 at all {nprimes} primes p>=5 dividing d; '
           f'|d| <= 6^5 N^(5/2) (max |d|/N^2.5 = {worst:.3g}), hence |d|^10 <= 6^50 N^25', falsifier=fb > 0)
    info('A5.16b', '518', '', f'sharpness of 5/2: max |d|/N^2 = {worst2:.3g} (unbounded along d=p^5, so exponent 2 fails)')

    # A5.18  l.522-524: D^{2-1/6} = D^{11/6}; about D^{1/6} copies when H ~ D
    report('A5.18', '522-524', '', F(2) - F(1, 6) == F(11, 6) and (F(2) - F(1, 6))/2 == F(11, 12),
           '2-1/6 = 11/6, half is 11/12', kind='I', falsifier=F(2) - F(1, 3) != F(11, 6))

    # A5.19  l.528-531: chi_n(v^3) = (v/n)_2 and chi_n(v^2) = (v/n)_3
    def sym_k(v, n, k):
        """(v/n)_k for k = 2, 3 as an index of zeta (k=2: 0 or 3; k=3: 0, 2, 4); None if not coprime."""
        tot = 0
        for pi in factor_primary(n):
            if divides(pi, v):
                return None
            r = powmod(v, (N(pi) - 1)//k, pi)
            js = [j for j in range(6) if divides(pi, sub(r, UNITS[j]))]
            assert len(js) == 1 and (js[0] in ((0, 3) if k == 2 else (0, 2, 4)))
            tot += js[0]
        return tot % 6
    cols = columns(3000)
    bad = 0; fb = 0
    for _ in range(2000):
        n = random.choice(cols); v = (random.randint(-80, 80), random.randint(-80, 80))
        if v == (0, 0):
            continue
        q2, q3 = sym_k(v, n, 2), sym_k(v, n, 3)
        bad += chi_idx(n, pw(v, 3)) != q2 or chi_idx(n, pw(v, 2)) != q3
        fb += chi_idx(n, pw(v, 3)) != q3
    report('A5.19', '528-531', '', bad == 0, '2000 random (n,v): chi_n(v^3)=(v/n)_2, chi_n(v^2)=(v/n)_3',
           falsifier=fb > 0)

    # A5.20  l.532-533: rational v: (v/n)_2 = (v/N(n)) Jacobi
    cols = columns(5000)
    bad = 0; fb = 0
    for v in [x_ for x_ in range(-60, 61) if x_]:
        for n in cols:
            j = chi_idx(n, (v**3, 0))
            jac = jacobi(v, N(n))
            val = 0 if j is None else (1 if j == 0 else -1)
            assert j in (None, 0, 3)
            bad += val != jac; fb += val != jacobi(-v, N(n))
    report('A5.20', '532-533', '', bad == 0, f'v in [-60,60], {len(cols)} columns N<=5000: (v/n)_2 = Jacobi (v/N n)',
           falsifier=fb > 0)

    # A5.21  l.534-537, eq:quadratic-twists: A_{u_d v^3}(D) = sum q_d(m) (v/m) W(m/D) (coefficientwise, exact)
    X = 15000; bn = by_norm(columns(X))
    bad = 0; fb = 0; tot = 0
    for d in [1, 7, -3]:
        an = ellan(d, X); ud = ((4*d)**5, 0)
        for v in [-1, 2, -2, 3, 5, -7, 10, 6, 15]:
            for m in range(1, X + 1):
                if math.gcd(m, 6) > 1:
                    continue
                tot += 1
                target = (Qd_int(d, m, lambda p: an[p])*jacobi(v, m), 0)
                bad += Rcoef(mul(ud, (v**3, 0)), m, bn) != target
                if m < 2000:
                    fb += Rcoef(mul(ud, (v**2, 0)), m, bn) != target
    report('A5.21', '534-537', 'eq:quadratic-twists', bad == 0,
           f'{tot} coefficients (d in 1,7,-3; v in -1,+-2,3,5,-7,10,6,15; m<={X}) exact', falsifier=fb > 0)

    # A5.22  l.538-540: every nonzero element of F_p is a square in F_{p^2}
    bad = 0; fb = 0
    for p in primes_upto(500):
        if p % 3 != 2 or p == 2:
            continue
        pi = (-p, 0)
        bad += sum(powmod((x_, 0), (p*p - 1)//2, pi) != (1, 0) for x_ in range(1, p))
        fb += all(powmod((x_, 1), (p*p - 1)//2, pi) == (1, 0) for x_ in range(p))
    report('A5.22', '538-540', '', bad == 0, 'x^((p^2-1)/2)=1 in F_{p^2} for x in F_p^*, inert p<500',
           falsifier=fb == 0)

    # A5.23  l.542-543: rational slice has 2 floor(sqrt H) ~ H^{1/2} rows, lattice ~ (2pi/sqrt3) H
    ok = all(sum(1 for u in range(-int(H_**.5) - 2, int(H_**.5) + 3) if 0 < u*u <= H_) == 2*math.isqrt(H_)
             for H_ in [1, 2, 10, 99, 100, 101, 12345, 10**6])
    report('A5.23', '542-543', '', ok, '#{u in Z: 0<u^2<=H} = 2 floor(sqrt H); lattice count in A4.01',
           falsifier=not all(sum(1 for u in range(-int(H_**.5) - 2, int(H_**.5) + 3) if 0 < u*u <= H_) == math.isqrt(H_)
                             for H_ in [10]))

    # A5.24  l.544-545: sixth powers of lattice elements ~ H^{1/6}; of rational integers ~ H^{1/12}
    msg = []; okc = True; fbc = False
    for R in [50, 200, 800]:            # R = H^{1/6}
        xs, ys = lattice_points(R)
        sixth = {pw((int(a_), int(b_)), 6) for a_, b_ in zip(xs, ys)}
        okc &= len(sixth)*6 == len(xs)
        fbc |= len(sixth)*3 != len(xs)
        msg.append(f'H^(1/6)={R}: #{{v^6}}={len(sixth)}, ratio {len(sixth)/R:.4f}')
    okr = all(len({v**6 for v in range(-50, 51) if v and v**12 <= H_}) == int(round(H_**(1/12) + 1e-9)) // 1
              for H_ in [2**12, 3**12, 10**12])
    report('A5.24', '544-545', '', okc and okr, '; '.join(msg) + f' (limit (2pi/sqrt3)/6={2*math.pi/SQ3/6:.4f}); '
           'rational sixth powers: floor(H^(1/12))', falsifier=fbc, kind='E')
    untested('A5.25', '546-548', '', 'remarks (lattice primes; rational numerators still need sextic symbols): '
             'not numerical claims')

if __name__ == '__main__':
    run()
