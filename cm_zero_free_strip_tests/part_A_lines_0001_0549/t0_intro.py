"""Part A, lines 1-155: abstract, 'method in six sentences', 'higher weight', 'Why start with 11/12'."""
from cmcore import *
from fractions import Fraction as F
import sympy as sp

def run():
    # A0.01  l.33-35: E_d are exactly the j = 0 curves; CM by Z[omega]
    a, b, x, y, w = sp.symbols('a b x y w')
    j = 1728*4*a**3/(4*a**3 + 27*b**2)
    ok1 = sp.solve(sp.Eq(sp.numer(sp.together(j)), 0), a) == [0]
    out = gp("for(d=-200,200,if(d,print(ellinit([0,0,0,0,d]).j)))")
    ok2 = len(out) == 400 and all(t == '0' for t in out)
    # automorphism (x, y) -> (omega x, y) of order 3 together with -1 gives mu_6 = Z[omega]^*
    om = sp.Rational(-1, 2) + sp.sqrt(3)*sp.I/2
    ok3 = sp.simplify(sp.expand((om*x)**3 - x**3)) == 0
    fals = gp("print(ellinit([0,0,0,1,1]).j)")[0] != '0'          # y^2 = x^3 + x + 1 has j != 0
    report('A0.01', '33-35', '', ok1 and ok2 and ok3,
           'j(E_d)=0 for 0<|d|<=200 (PARI); short Weierstrass j=0 <=> a=0; (x,y)->(omega x,y) is an automorphism',
           falsifier=fals)

    untested('A0.02', '38-39,55-57', '', 'main theorems L(E_d,s)!=0 (Re s>17/12), Re s>11/8, L(s,conj(alpha)eta)!=0 '
             '(Re s>7/8): conditional analytic conclusions, not finitely testable')

    # A0.03  l.41-42, 172-174: critical strip 1/2<Re s<3/2 (weight-2 functional equation s <-> 2-s)
    out = gp("E=ellinit([0,0,0,0,1]);L=lfuncreate(E);print(lfunparams(L)[2]);print(lfuncheckfeq(L));"
             "L2=lfuncreate([n->ellan(E,n),0,[0,1],1,36,1]);print(lfuncheckfeq(L2));"
             "E7=ellinit([0,0,0,0,7]);print(lfuncheckfeq(lfuncreate(E7)))")
    k = int(out[0]); feq = int(out[1]); feq_wrong = int(out[2]); feq7 = int(out[3])
    strip = (F(k - 1, 2), F(k + 1, 2))
    ok = k == 2 and feq < -60 and feq7 < -60 and strip == (F(1, 2), F(3, 2))
    report('A0.03', '41-42,172-174', '', ok,
           f'PARI: weight k={k}, lfuncheckfeq={feq},{feq7} bits; strip ((k-1)/2,(k+1)/2)={strip}; '
           f'shifted L(E,s+1/2) has centre 1/2, unshifted centre 1',
           falsifier=feq_wrong > -20)

    # A0.04  l.43-45, 104, 149: 11/12 + 1/2 = 17/12 and 7/8 + 1/2 = 11/8
    ok = F(11, 12) + F(1, 2) == F(17, 12) and F(7, 8) + F(1, 2) == F(11, 8)
    report('A0.04', '43-45,104,149', '', ok, '11/12+1/2=17/12, 7/8+1/2=11/8 (Fraction)',
           falsifier=not (F(11, 12) + 1 == F(17, 12)), kind='I')

    # A0.05  l.45-46: Hecke character of infinite order: psi_d(pi)^k = 1 would force pi ~ conj(pi)
    bad = 0; cnt = 0
    for p in primes_upto(20000):
        if p % 3 == 1:
            pi = primes_above(p)[0]; cnt += 1
            # psi(pi)^12 = (pi/conj pi)^6 zeta^12 = 1  <=>  pi^6 = conj(pi)^6  <=>  pi = unit*conj(pi)
            if any(mul(u, conj(pi)) == pi for u in UNITS):
                bad += 1
    fals = all(any(mul(u, conj(q)) == q for u in UNITS) for q in [(-p, 0) for p in (5, 11, 17)])  # inert: finite order
    report('A0.05', '45-46', '', bad == 0,
           f'{cnt} split primes p<20000: pi not associate to conj(pi), so alpha(pi) is not a root of unity',
           falsifier=fals)

    # A0.06  l.64-65, 120-130: eta(3z)^8 (weight 4) and eta(2z)^3 eta(6z)^3 (weight 3) are CM by Q(sqrt-3) with
    # Hecke character conj(alpha)^k eta, k+1 = weight: a_p = 0 (p = 2 mod 3) and a_p = Tr(zeta^j conj(pi)^k).
    M = 3000
    def eta_prod(facs):
        c = np.zeros(M + 1, dtype=np.int64); c[0] = 1
        for (m, e) in facs:
            for _ in range(e):
                for n in range(1, M//m + 1):
                    s_ = n*m
                    c[s_:] = c[s_:] - c[:-s_].copy()
        return [int(t) for t in c]
    def check(facs, shift, k):
        c = eta_prod(facs)
        a = lambda n: c[n - shift] if n - shift >= 0 else 0
        bad = 0
        for p in primes_upto(M):
            if p in (2, 3):
                continue
            ap = a(p)
            if p % 3 == 2:
                bad += (ap != 0)
            else:
                pi = primes_above(p)[0]
                tr = [2*z[0] - z[1] for z in (mul(UNITS[j], pw(conj(pi), k)) for j in range(6))]
                bad += (ap not in tr)
                bad += (abs(ap) > 2*p**(k/2) + 1e-9)
        return bad
    b1 = check([(3, 8)], 1, 3)                # eta(3z)^8 = q prod (1-q^{3n})^8
    b2 = check([(2, 3), (6, 3)], 1, 2)        # eta(2z)^3 eta(6z)^3 = q prod (1-q^{2n})^3 (1-q^{6n})^3
    wts = (F(8, 2), F(3 + 3, 2))
    fals = check([(1, 2), (11, 2)], 1, 1) > 0  # eta(z)^2 eta(11z)^2: weight 2, not CM
    report('A0.06', '64-65,120-130', '', b1 == 0 and b2 == 0 and wts == (4, 3),
           f'q-expansions to q^{M}: weights {wts}; a_p=0 for p=2 mod 3 and a_p=Tr(zeta^j conj(pi)^k) (k=3,2), '
           f'|a_p|<=2p^(k/2); failures {b1},{b2}', falsifier=fals)

    # A0.07  l.73-74: 1/L_p(E_d,s+1/2) = 1 - lambda p^{-s} + p^{-2s} as a polynomial in X = p^{-s}: coefficients
    # q(1)=1, q(p)=-lambda, q(p^2)=1, q(p^j)=0 (j>=3)
    X, lam = sp.symbols('X lam')
    ser = sp.series(1/(1/(1 - lam*X + X**2)), X, 0, 6).removeO()
    coeffs = [sp.expand(ser).coeff(X, i) for i in range(6)]
    ok = coeffs == [1, -lam, 1, 0, 0, 0]
    report('A0.07', '72-74', '', ok, f'coefficients of 1/L_p in X=p^-s: {coeffs} (eq:reciprocal tested in A2.07)',
           falsifier=sp.expand(sp.series(1 - lam*X - X**2, X, 0, 6).removeO()).coeff(X, 2) != 1)

    # A0.08  l.76-77, 223-228: c(m) = sum_{rt=m} mu(r) mu(t) chi_{-3}(t) are the coefficients of 1/(zeta L(chi_{-3}))
    M = 20000
    def c(m):
        return sum(mobius(r)*mobius(m//r)*chi_m3(m//r) for r in divisors(m))
    def divisors(m):
        ds = [1]
        for p, e in factorint(m).items():
            ds = [d*p**i for d in ds for i in range(e + 1)]
        return ds
    cc = [0] + [c(m) for m in range(1, M + 1)]
    rr = [0] + [sum(chi_m3(e) for e in divisors(k)) for k in range(1, M + 1)]   # zeta(s) L(s,chi_{-3}) coeffs
    bad = sum(1 for m in range(1, M + 1) if sum(cc[m//k]*rr[k] for k in divisors(m)) != (1 if m == 1 else 0))
    cw = [0] + [sum(mobius(r)*mobius(m//r) for r in divisors(m)) for m in range(1, 2001)]
    fals = sum(1 for m in range(1, 2001) if sum(cw[m//k]*rr[k] for k in divisors(m)) != (1 if m == 1 else 0)) > 0
    report('A0.08', '76-77,225', '', bad == 0, f'(c * (zeta L(chi_-3)))(m) = [m=1] for m<={M}', falsifier=fals)

    # A0.09  l.82, 95: zeta sum: sum_{N(n)=m} mu_K(n) = c(m) for (m,6)=1 (primary n coprime to 6)
    X = 20000
    bn = by_norm(columns(X))
    bad = 0; fbad = 0
    for m in range(1, X + 1):
        if math.gcd(m, 6) > 1:
            continue
        s = sum(muK(n) for n in bn.get(m, []))
        bad += (s != cc[m])
        fbad += (s != mobius(m))
    report('A0.09', '82,95', '', bad == 0, f'sum_(N n=m) mu_K(n) = c(m), all m<={X} coprime to 6',
           falsifier=fbad > 0)

    # A0.10  l.86-88: D^{2+theta+eps} = D^{1+eps} H with H = D^{1+theta}
    th, ep = sp.symbols('theta epsilon', positive=True)
    ok = sp.simplify((2 + th + ep) - ((1 + ep) + (1 + th))) == 0
    report('A0.10', '87-88', '', ok, 'exponent of D^{1+eps}H equals 2+theta+eps (sympy)',
           falsifier=sp.simplify((2 + th + ep) - ((1 + ep) + 2*(1 + th))) != 0, kind='I')
    info('A0.11', '86-88,378', '', 'rows/columns counts: see A4.01 (rows ~ (2pi/sqrt3)H, columns ~ const*D)')

    # A0.12  l.91-95: H^{1/6}|S|^2 <~ DH  =>  |S| <~ D^{11/12} when H = D
    e = (1 + F(1)) - F(1, 6)    # exponent of D in DH/H^{1/6} at H = D
    ok = e/2 == F(11, 12)
    report('A0.12', '91-94', '', ok, f'|S|^2 <~ D H^(5/6) = D^{e} at H=D, so |S| <~ D^{e/2}', kind='I',
           falsifier=((1 + F(1)) - F(1, 3))/2 != F(11, 12))
    info('A0.13', '91-93', '', 'row count ~H^{1/6}/log H and O(D/N(pi)) approximation: tested in A5.03, A5.05-A5.07')

    untested('A0.14', '96-104', '', 'summary of Poisson/theta/descent steps (Sections 6-8; other chunks). Text check '
             'only: "(nbar/|n|)^2 ... (zeta: nbar/|n|)" agrees with a_xi(n)=conj(alpha(n))^2 gamma_2(n)xi(n) '
             'and "zeta version produces conj(alpha) gamma_2" at l.661-674')

    # A0.15  l.122-124 + l.57: E_d (k=1) corresponds to alpha*chi_{u_d} (Part I) and to conj(alpha)*conj(chi_{u_d})
    # (Part II convention): same L-function because a_m is real.
    d = 7; Mx = 3000
    an = ellan(d, Mx); bn = by_norm(columns(Mx))
    bad = 0; fb = 0
    for m in range(1, Mx + 1):
        if math.gcd(m, 6*d) > 1:
            continue
        s1 = (0, 0); s2 = (0, 0); s3 = (0, 0)
        for n in bn.get(m, []):
            jj = chi_idx(n, (pow(4*d, 5), 0))
            if jj is None:
                continue
            s1 = add(s1, mul(UNITS[jj], n))
            s2 = add(s2, mul(UNITS[(-jj) % 6], conj(n)))
            s3 = add(s3, mul(UNITS[(-jj) % 6], n))          # wrong: alpha * conj(chi)
        bad += (s1 != (an[m], 0)) + (s2 != (an[m], 0))
        fb += (s3 != (an[m], 0))
    report('A0.15', '57,122-124', '', bad == 0,
           f'd={d}, m<={Mx}: sum n chi_n(u_d) = sum conj(n) conj(chi_n(u_d)) = a_m (alpha vs conj(alpha) is a '
           f'convention: same L(E_d,s))', falsifier=fb > 0)

    # A0.16  l.66-67, 129-130: Sym^k L(E_d) factors into Hecke L-functions of psi^(k-2j) (and zeta/L(chi_-3) if
    # k even) -- check the local factorisation at split and inert primes (expectation itself untested).
    T = sp.symbols('T')
    p = sp.Symbol('p', positive=True); A = sp.Symbol('A')
    def symfact_bad(wrong=False):
        bad = 0
        for k in range(1, 7):
            # split p: Frobenius eigenvalues A, B = p/A
            B = p/A
            lhs = sp.prod([1 - A**(k - i)*B**i*T for i in range(k + 1)])
            rhs = 1
            for i in range((k + 1)//2):
                rhs *= (1 - A**(k - 2*i)*p**i*T)*(1 - B**(k - 2*i)*p**i*T)   # L(s-i, psi^(k-2i)) at pi, pibar
            if k % 2 == 0:
                rhs *= (1 - p**(k//2)*T)
            bad += sp.simplify(sp.expand(lhs - rhs)) != 0
            # inert p: eigenvalues +-i sqrt(p); psi((p)) = psi(-p) = -p (primary generator -p, norm p^2)
            lhs = sp.prod([1 - (sp.I*sp.sqrt(p))**(k - i)*(-sp.I*sp.sqrt(p))**i*T for i in range(k + 1)])
            rhs = 1
            for i in range((k + 1)//2):
                if wrong:
                    rhs *= (1 - (-p)**(k - 2*i)*p**(2*i)*T**2 if k % 2 else 1)*(1 if k % 2 else
                            (1 - (-1)**(k//2)*p**(k//2)*T)**2)
                else:
                    rhs *= (1 - (-p)**(k - 2*i)*p**(2*i)*T**2)
            if k % 2 == 0:
                rhs *= (1 - (-1)**(k//2)*p**(k//2)*T)             # zeta (k=0 mod 4) or L(chi_-3) (k=2 mod 4)
            bad += sp.simplify(sp.expand(lhs - rhs)) != 0
        return bad
    bad = symfact_bad()
    report('A0.16', '66-67,129-130', '', bad == 0,
           'local factorisation of Sym^k L(E_d) (k<=6) into psi^(k-2i)-factors (+ zeta/L(chi_-3) for k even) at '
           'split and inert p (the zero-free consequence itself is an expectation; untested)', kind='C',
           falsifier=symfact_bad(wrong=True) > 0)

if __name__ == '__main__':
    run()
