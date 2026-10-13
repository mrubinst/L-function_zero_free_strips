"""Part A, Section 'The triangular lattice and the curve's row' (lines 244-370)."""
from cmcore import *
import sympy as sp
random.seed(20261012)

def zsum(terms):
    s = (0, 0)
    for t in terms:
        s = add(s, t)
    return s

def row_coeff(u, m, bn, with_mu=True):
    """sqrt(m) * sum_{N(n)=m} mu_K(n) alpha(n) chi_n(u) = sum mu_K(n) n zeta^j, exactly in Z[omega]."""
    s = (0, 0)
    for n in bn.get(m, []):
        mu_ = muK(n) if with_mu else 1
        if mu_ == 0:
            continue
        j = chi_idx(n, u)
        if j is None:
            continue
        t = mul(UNITS[j], n)
        s = add(s, t if mu_ == 1 else neg(t))
    return s

def run():
    # A3.01  l.248-250: N(n) = |n|^2 = A^2 - AB + B^2
    err = 0.0; fb = 0
    for A in range(-60, 61):
        for B in range(-60, 61):
            err = max(err, abs(abs(cx((A, B)))**2 - (A*A - A*B + B*B)))
            fb += abs(abs(cx((A, B)))**2 - (A*A + A*B + B*B)) > 1e-6
    report('A3.01', '248-250', '', err < 1e-8, f'|A+B omega|^2 = A^2-AB+B^2 on |A|,|B|<=60 (max err {err:.1e})',
           falsifier=fb > 0)

    # A3.02  l.252-253: six units +-1, +-omega, +-omega^2; unique factorisation (norm-Euclidean spot check, CITED)
    units = [(A, B) for A in range(-3, 4) for B in range(-3, 4) if N((A, B)) == 1]
    expect = {(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)}       # omega^2 = -1 - omega
    ok = set(units) == expect and len(units) == 6
    eu = 0
    for _ in range(20000):
        x = (random.randint(-10**6, 10**6), random.randint(-10**6, 10**6))
        y = (random.randint(-10**4, 10**4), random.randint(-10**4, 10**4))
        if y == (0, 0):
            continue
        r = sub(x, mul((round(mul(x, conj(y))[0]/N(y)), round(mul(x, conj(y))[1]/N(y))), y))
        eu += not (4*N(r) < 4*N(y))
    report('A3.02', '252-253', '', ok and eu == 0,
           f'units = {sorted(units)}; Euclidean division by rounding gave N(r)<N(y) in 20000 random cases (UFD: CITED)',
           falsifier=len([(A, B) for A in range(-3, 4) for B in range(-3, 4) if A*A + B*B == 1]) != 6)

    # A3.03  l.253-258: unique primary associate of each element coprime to 3; primary <=> A=1, B=0 (mod 3)
    bad = 0; fb = 0
    for A in range(-40, 41):
        for B in range(-40, 41):
            n = (A, B)
            if n == (0, 0) or divides(LAM, n):
                continue
            ass = [mul(u, n) for u in UNITS]
            prim = [z for z in ass if z[0] % 3 == 1 and z[1] % 3 == 0]
            cong = [z for z in ass if divides((3, 0), sub(z, (1, 0)))]
            bad += (len(prim) != 1) or (set(prim) != set(cong))
            fb += len([z for z in ass if (z[0] % 3 in (1, 2)) and z[1] % 3 == 0]) == 1    # '= +-1 mod 3'
    report('A3.03', '253-258', '', bad == 0,
           'every n coprime to 3 in |A|,|B|<=40 has exactly one associate with A=1,B=0 (3), and it is n=1 mod 3',
           falsifier=fb == 0)

    # A3.04  l.260-261: coprime to 2 <=> A, B not both even (2 is an inert prime: no element of norm 2)
    no2 = all(N((A, B)) != 2 for A in range(-5, 6) for B in range(-5, 6))
    bad = sum(divides((2, 0), (A, B)) != (A % 2 == 0 and B % 2 == 0) for A in range(-50, 51) for B in range(-50, 51))
    fb = sum(divides((2, 0), (A, B)) != (A % 2 == 0) for A in range(-50, 51) for B in range(-50, 51))
    report('A3.04', '260-261', '', no2 and bad == 0, '2 | A+B omega <=> A,B both even; no element of norm 2',
           falsifier=fb > 0)

    # A3.05  l.264-271: mu_K(-5)=-1 (N=25), mu_K(7)=1 with 7=(1+3w)(-2-3w); mu_K(n) != mu(N(n)) in general
    f7 = factor_primary((7, 0))
    ok = (muK((-5, 0)) == -1 and N((-5, 0)) == 25 and muK((7, 0)) == 1
          and mul((1, 3), (-2, -3)) == (7, 0) and primary((1, 3)) and primary((-2, -3))
          and N((1, 3)) == 7 and N((-2, -3)) == 7 and set(f7) == {(1, 3), (-2, -3)} and mobius(25) == 0)
    report('A3.05', '264-271', '', ok,
           f'mu_K(-5)={muK((-5, 0))}, mu(25)={mobius(25)}, 7=(1+3w)(-2-3w) both primary of norm 7, mu_K(7)=1',
           falsifier=mul((1, 3), (2, 3)) != (7, 0))

    # A3.06  l.273-278: Euler criterion gives exactly one sixth root of unity; chi multiplicative in u; chi_1 = 1
    PR = [pi for p in primes_upto(400) if p > 3 for pi in primes_above(p) if N(pi) <= 400*400]
    bad = 0; fb = 0
    for _ in range(4000):
        pi = random.choice(PR)
        u1 = (random.randint(-999, 999), random.randint(-999, 999))
        u2 = (random.randint(-999, 999), random.randint(-999, 999))
        if not divides(pi, u1):
            v = powmod(u1, (N(pi) - 1)//6, pi)
            bad += sum(divides(pi, sub(v, z)) for z in UNITS) != 1
        j1, j2, j12 = sym6(u1, pi), sym6(u2, pi), sym6(mul(u1, u2), pi)
        if j1 is None or j2 is None:
            bad += j12 is not None
        else:
            bad += j12 != (j1 + j2) % 6
        js = sym6(add(u1, u2), pi)
        if j1 is not None and j2 is not None and js is not None:
            fb += js != (j1 + j2) % 6
    bad += sum(chi_idx((1, 0), (random.randint(-99, 99), random.randint(-99, 99))) != 0 for _ in range(100))
    report('A3.06', '273-278', '', bad == 0,
           '4000 random (u1,u2,pi): u^((Npi-1)/6) = exactly one unit mod pi; (u1u2/pi)_6=(u1/pi)(u2/pi); chi_1=1',
           falsifier=fb > 0)

    # columns and norm lists used below
    X = 20000
    cols = columns(X); bn = by_norm(cols)

    # A3.07  l.287-290, 300-303: conj(chi_n(4d)) = chi_n(4d)^5 = chi_n((4d)^5)  (including zero values)
    bad = 0; fb = 0; tot = 0
    for d in [1, 2, -3, 7, 13, -35, 16, 91]:
        for n in cols[:1500]:
            j1 = chi_idx(n, (4*d, 0)); j5 = chi_idx(n, ((4*d)**5, 0)); j2 = chi_idx(n, ((4*d)**2, 0))
            tot += 1
            if j1 is None:
                bad += j5 is not None
            else:
                bad += j5 != (-j1) % 6 or j5 != (5*j1) % 6
                fb += j2 != (-j1) % 6
    report('A3.07', '287-290,300-303', 'eq:ud', bad == 0, f'{tot} pairs (d,n): index((4d)^5) = 5*index(4d) = '
           f'-index(4d) mod 6, zeros preserved', falsifier=fb > 0)

    # A3.08  l.291, 339-343, eq:normgroup: q_d(m) = sum_{N(n)=m} mu alpha chi_n(u_d) (exact, scaled by sqrt m) and
    # S_d(D) = A_{u_d}(D) (numeric)
    bad = 0; fb1 = 0; fb2 = 0; tot = 0
    DL = [1, 2, -3, 7, 16, -11, 13, 5]
    for d in DL:
        an = ellan(d, X); apf = lambda p, an=an: an[p]
        ud = ((4*d)**5, 0)
        for m in range(1, X + 1):
            if math.gcd(m, 6) > 1:
                bad += m in bn
                continue
            tot += 1
            lhs = (Qd_int(d, m, apf), 0)
            bad += row_coeff(ud, m, bn) != lhs
            if m <= 3000:
                fb1 += row_coeff((4*d, 0), m, bn) != lhs             # u = 4d (not conjugated)
                fb2 += row_coeff((d**5, 0), m, bn) != lhs            # u = d^5 (no factor 4)
    # sums
    serr = 0.0
    for d in [1, 7]:
        an = ellan(d, 16000); apf = lambda p, an=an: an[p]
        ud = ((4*d)**5, 0)
        for D in [500, 2000, 8000]:
            S = sum(Qd_int(d, m, apf)/math.sqrt(m)*W_test(np.array([m/D]))[0] for m in range(1, 2*D + 1))
            A_ = sum(muK(n)*alpha(n)*chi(n, ud)*W_test(np.array([N(n)/D]))[0] for n in cols if N(n) <= 2*D)
            serr = max(serr, abs(S - A_))
    report('A3.08', '291,339-343', 'eq:normgroup', bad == 0 and serr < 1e-8,
           f'exact for d in {DL}, all m<={X} ({tot} m coprime to 6; empty for (m,6)>1); |S_d(D)-A_ud(D)| <= {serr:.1e} '
           f'(d=1,7; D=500,2000,8000); falsifiers u=4d: {fb1}, u=d^5: {fb2} mismatches', falsifier=fb1 > 0 and fb2 > 0)

    # A3.09  l.304-306: chi_n(u v^6) = chi_n(u) 1_{(n,v)=1}
    bad = 0; fb = 0
    for _ in range(3000):
        n = random.choice(cols[:3000])
        u = (random.randint(-50, 50), random.randint(-50, 50))
        v = (random.randint(-12, 12), random.randint(-12, 12))
        if u == (0, 0) or v == (0, 0):
            continue
        cop = all(not divides(pi, v) for pi in set(factor_primary(n)))
        ju = chi_idx(n, u); juv = chi_idx(n, mul(u, pw(v, 6))); juv3 = chi_idx(n, mul(u, pw(v, 3)))
        expect = ju if cop else None
        bad += juv != expect
        fb += juv3 != expect
    report('A3.09', '304-306', '', bad == 0, '3000 random (n,u,v): chi_n(uv^6) = chi_n(u) 1_(n,v)=1',
           falsifier=fb > 0)

    # A3.10  l.306-308: v supported on 2, 3 changes nothing; 16 d^5 gives the row of (4d)^5 since (4d)^5 = 16d^5 2^6
    ok = all((4*d)**5 == 16*d**5*2**6 for d in range(-50, 51))
    bad = 0; fb = 0
    for d in [1, 2, -3, 7, 11, -13]:
        for n in cols[:3000]:
            bad += chi_idx(n, (16*d**5, 0)) != chi_idx(n, ((4*d)**5, 0))
            bad += chi_idx(n, mul((5, 1), pw((3, 0), 6))) != chi_idx(n, (5, 1))          # v = 3
            bad += chi_idx(n, mul((5, 1), pw(LAM, 6))) != chi_idx(n, (5, 1))            # v = lambda
            fb += chi_idx(n, (4*d**5, 0)) != chi_idx(n, ((4*d)**5, 0))
    report('A3.10', '306-308', '', ok and bad == 0, 'rows 16d^5 and (4d)^5 agree at 3000 columns; rows times 3^6, '
           'lambda^6 unchanged', falsifier=fb > 0)

    # A3.11  l.309-310: d=16: 16*16^5 = 2^24 = (2^4)^6; y^2=x^3+16 has conductor 27 and is the row u=1;
    # y^2=x^3+1 is the row u=16
    ok = 16*16**5 == 2**24 == (2**4)**6 and conductor(16) == 27
    an16 = ellan(16, X); an1 = ellan(1, X)
    bad = 0; fb = 0
    for m in range(1, X + 1):
        if math.gcd(m, 6) > 1:
            continue
        r1 = row_coeff((1, 0), m, bn); r16 = row_coeff((16, 0), m, bn)
        bad += r1 != (Qd_int(16, m, lambda p: an16[p]), 0)
        bad += r16 != (Qd_int(1, m, lambda p: an1[p]), 0)
        fb += r1 != (Qd_int(1, m, lambda p: an1[p]), 0)
    report('A3.11', '309-310', '', ok and bad == 0,
           f'16^6=2^24, cond(E_16)={conductor(16)}; row u=1 = q_16 and row u=16 = q_1, all m<={X}', falsifier=fb > 0)

    # A3.12  l.311-312: every nonzero rational row u is the row of E_d with d = 4^5 u^5
    ok = all((4*(4**5*u**5))**5 == u*(4**5*u**4)**6 for u in range(-40, 41) if u)
    Xs = 3000
    bad = 0; fb = 0; tot = 0
    for u in [v for v in range(-30, 31) if v]:
        d = 4**5*u**5; an = ellan(d, Xs); apf = lambda p, an=an: an[p]
        dw = 4*u**5; anw = ellan(dw, 600)
        for m in range(1, Xs + 1):
            if math.gcd(m, 6) > 1:
                continue
            tot += 1
            r = row_coeff((u, 0), m, bn)
            bad += r != (Qd_int(d, m, apf), 0)
            if m <= 600:
                fb += r != (Qd_int(dw, m, lambda p: anw[p]), 0)
    report('A3.12', '311-312', '', ok and bad == 0,
           f'(4d)^5 = u (4^5u^4)^6 for d=4^5u^5; rows u in [-30,30]\\0 equal q_d (PARI) for all m<={Xs} ({tot} checks)',
           falsifier=fb > 0)

    # A3.13  l.312-315: rows not rational x sixth power have non-real coefficients in general; rational x sixth power
    # rows are real
    Xn = 3000
    nonreal = {}
    for u in [(0, 1), (1, 3), (3, 1), (5, 2), (2, 7), (-4, 1)]:
        nonreal[u] = next((m for m in range(1, Xn + 1) if math.gcd(m, 6) == 1 and row_coeff(u, m, bn)[1] != 0), None)
    realbad = 0; altered = 0
    for r in [2, -5, 7, 1]:
        for v in [(1, 3), (2, 3), (4, 1), (1, 1)]:
            u = mul((r, 0), pw(v, 6))
            for m in range(1, Xn + 1):
                if math.gcd(m, 6) > 1:
                    continue
                im = row_coeff(u, m, bn)[1] != 0
                if math.gcd(m, N(v)) == 1:
                    realbad += im                 # away from the finitely many Euler factors at primes of v
                else:
                    altered += im                 # l.304-306: only these Euler factors change
    ok = all(m is not None for m in nonreal.values()) and realbad == 0
    report('A3.13', '312-315', '', ok,
           f'first m with non-real q_u(m): {nonreal}; rows r*v^6 (r rational) are real at all m<={Xn} coprime to N(v) '
           f'({altered} non-real values only at m sharing a prime with N(v), i.e. the altered Euler factors)',
           falsifier=not any(row_coeff((7, 0), m, bn)[1] != 0 for m in range(1, 2000)))  # detector silent on a real row

    # A3.14  l.295-297: discriminant of t^2-t+a is 1-4a; J(kappa,rho) = kappa(4) J(kappa,kappa)
    t, a = sp.symbols('t a')
    okd = sp.expand(sp.discriminant(t**2 - t + a, t) - (1 - 4*a)) == 0
    bad = 0; fb1 = 0; fb2 = 0; cnt = 0
    for p in primes_upto(3000):
        if p % 3 != 1:
            continue
        g = int(sp.primitive_root(p))
        dlog = [0]*p; x = 1
        for k in range(p - 1):
            dlog[x] = k; x = x*g % p
        Wc = cmath.exp(2j*math.pi/3)
        for e in (1, 2):
            kap = lambda x_: 0 if x_ % p == 0 else Wc**(e*dlog[x_ % p] % 3)
            rho = lambda x_: 0 if x_ % p == 0 else (1 if dlog[x_ % p] % 2 == 0 else -1)
            Jkr = sum(kap(x_)*rho(1 - x_) for x_ in range(p))
            Jkk = sum(kap(x_)*kap(1 - x_) for x_ in range(p))
            cnt += 1
            bad += abs(Jkr - kap(4)*Jkk) > 1e-6
            fb1 += abs(Jkr - kap(2)*Jkk) > 1e-6
            fb2 += abs(Jkr - kap(4).conjugate()*Jkk) > 1e-6
    report('A3.14', '295-297', '', okd and bad == 0,
           f'disc(t^2-t+a)=1-4a; J(kappa,rho)=kappa(4)J(kappa,kappa) for {cnt} (p<3000, kappa cubic) cases',
           falsifier=fb1 > 0 and fb2 > 0)

    # A3.15  l.318-319: p = 1 (3) => p = pi pibar with both generators primary, N(pi) = p
    bad = 0; cnt = 0
    for p in primes_upto(100000):
        if p % 3 == 1:
            pis = primes_above(p); cnt += 1
            pi = pis[0]
            bad += not (len(pis) == 2 and mul(pi, conj(pi)) == (p, 0) and primary(pi) and primary(conj(pi))
                        and conj(pi) in pis and N(pi) == p)
    fb = any(N((A, B)) == 5 for A in range(-5, 6) for B in range(-5, 6))
    report('A3.15', '318-319', '', bad == 0, f'{cnt} primes p=1 mod 3 below 1e5: p = pi*conj(pi), both primary',
           falsifier=not fb)

    # A3.16  l.320-322, eq:point-identity: a_p(d) = pi conj(chi_pi(4d)) + conj(pi) chi_pi(4d) (exact in Z[omega])
    bad = 0; f1 = 0; f2 = 0; f3 = 0; tot = 0
    for d in [x for x in range(-60, 61) if x]:
        out = gp(f"E=ellinit([0,0,0,0,{d}]);forprime(p=7,1500,if(p%3==1&&{d}%p,print(p,\" \",ellap(E,p))))")
        for line in out:
            p, ap = map(int, line.split())
            for pi in primes_above(p):
                j = sym6((4*d, 0), pi); tot += 1
                val = add(mul(pi, UNITS[(-j) % 6]), mul(conj(pi), UNITS[j]))
                bad += val != (ap, 0)
                f1 += add(mul(pi, UNITS[j]), mul(conj(pi), UNITS[(-j) % 6])) != (ap, 0)
                jd = sym6((d, 0), pi)
                f2 += add(mul(pi, UNITS[(-jd) % 6]), mul(conj(pi), UNITS[jd])) != (ap, 0)
                f3 += add(mul(neg(pi), UNITS[(-j) % 6]), mul(neg(conj(pi)), UNITS[j])) != (ap, 0)
    report('A3.16', '320-322', 'eq:point-identity', bad == 0,
           f'{tot} (d,pi) with 0<|d|<=60, p<1500 split, p!|6d (PARI ellap), both pi above p; falsifiers: chi instead '
           f'of conj {f1}, chi_pi(d) {f2}, -pi (IR sign) {f3} mismatches', falsifier=min(f1, f2, f3) > 0)

    # A3.17  l.324-325: Ireland-Rosen Ch.18 (CITED): with pi' = -pi = 2 mod 3, #E = p+1 + conj(chi) pi' + chi conj(pi')
    bad = 0; fb = 0
    for d in [1, 2, 3, 5, -7]:
        for p in [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97]:
            if d % p == 0:
                continue
            pi = primes_above(p)[0]; pip = neg(pi); j = sym6((4*d, 0), pi)
            Np = p + 1 - ap_count(d, p)
            ir = add(mul(UNITS[(-j) % 6], pip), mul(UNITS[j], conj(pip)))      # with pi' = -pi = 2 mod 3
            bad += (Np - (p + 1), 0) != ir
            ir1 = add(mul(UNITS[(-j) % 6], pi), mul(UNITS[j], conj(pi)))       # same form with pi = 1 mod 3
            fb += (Np - (p + 1), 0) != ir1
    report('A3.17', '324-325', 'IR Ch.18', bad == 0,
           "cited form #E_d(F_p) = p+1+conj(chi)pi'+chi conj(pi') holds with pi'=-pi=2 mod 3 (sign change as stated)",
           kind='C', falsifier=fb > 0)

    untested('A3.18', '326', '', 'statement that 480 pairs were checked (Section 78) -- superseded by A3.16 '
             '(13764 (d,p) pairs = 27528 (d,pi))', kind='C')

    # A3.19  l.327: chi_{conj pi}(x) = conj(chi_pi(x)) for rational x
    bad = 0; fb = 0
    for p in primes_upto(600):
        if p % 3 != 1:
            continue
        pi = primes_above(p)[0]
        for x_ in list(range(-30, 31)) + [4*7, 4*13]:
            j1 = sym6((x_, 0), pi); j2 = sym6((x_, 0), conj(pi))
            bad += (j1 is None) != (j2 is None) or (j1 is not None and j2 != (-j1) % 6)
        for z in [(2, 1), (5, 3), (1, 4)]:
            j1 = sym6(z, pi); j2 = sym6(z, conj(pi))
            fb += j1 is not None and j2 is not None and j2 != (-j1) % 6
    report('A3.19', '327', '', bad == 0, 'chi_conj(pi)(x) = conj chi_pi(x), x in [-30,30], p<600 split',
           falsifier=fb > 0)

    # A3.20  l.328-332: (1-psi(pi)X)(1-psi(pibar)X) = 1 - lambda X + X^2
    err = 0.0; fb = 0
    for d in [1, 2, -3, 7, -11]:
        an = ellan(d, 3000)
        for p in primes_upto(3000):
            if p % 3 != 1 or (6*d) % p == 0:
                continue
            pi = primes_above(p)[0]; pb = conj(pi)
            ps1 = alpha(pi)*chi(pi, ((4*d)**5, 0)); ps2 = alpha(pb)*chi(pb, ((4*d)**5, 0))
            err = max(err, abs(ps1*ps2 - 1), abs(ps1 + ps2 - an[p]/math.sqrt(p)))
            w1 = alpha(pi)*chi(pi, (4*d, 0)); w2 = alpha(pb)*chi(pb, (4*d, 0))
            fb += abs(w1 + w2 - an[p]/math.sqrt(p)) > 1e-9
    report('A3.20', '328-332', '', err < 1e-9, f'psi_d(pi)psi_d(pibar)=1, sum=lambda_d(p): max err {err:.1e}',
           falsifier=fb > 0)

    # A3.21  l.333-337: inert p: primary prime -p of norm p^2; chi_{-p}(x)=1 for rational x prime to p;
    # (p^2-1)/6 = (p-1)(p+1)/6 with 6 | p+1; alpha(-p) = -1; local factor 1 + p^{-2s}
    bad = 0; fb = 0
    for p in primes_upto(400):
        if p % 3 != 2 or p == 2:
            continue
        pi = (-p, 0)
        bad += not (primary(pi) and N(pi) == p*p and (p + 1) % 6 == 0 and (p*p - 1)//6 == (p - 1)*((p + 1)//6))
        bad += abs(alpha(pi) + 1) > 1e-15
        for x_ in list(range(1, p)) + [-1, -2, -p - 1]:
            bad += sym6((x_, 0), pi) != 0
        fb += all(sym6((a_, 1), pi) == 0 for a_ in range(p))
        bad += sym6((3*p, 0), pi) is not None
    # local factor: (1 - psi(-p) N(-p)^-s) = 1 - (-1)(1) p^{-2s} = 1 + p^{-2s}
    report('A3.21', '333-337', '', bad == 0, 'inert 5<=p<400: -p primary of norm p^2; chi_-p(x)=1 for rational x prime '
           'to p, 0 if p|x; alpha(-p)=-1, so factor 1+p^-2s', falsifier=fb == 0)

    # A3.22  l.336-337: p | d => symbols vanish at primes above p
    bad = 0
    for d in [5, 7, 35, -91, 13*19]:
        for p in factorint(abs(d)):
            for pi in primes_above(p):
                bad += chi_idx(pi, ((4*d)**5, 0)) is not None
    report('A3.22', '336-337', '', bad == 0, 'chi_pi(u_d) = 0 for pi above p | d', falsifier=chi_idx((1, 3), (5, 0)) is not None)

    # A3.23  l.344-345: columns coprime to 6 have norm coprime to 6 (sum empty when (m,6)>1), and conversely
    cols6 = columns(100000)
    bad = sum(math.gcd(N(n), 6) != 1 for n in cols6)
    allprim = [(A, B) for A in range(-400, 401) for B in range(-400, 401)
               if A % 3 == 1 and B % 3 == 0 and 0 < N((A, B)) <= 100000]
    bad += sum((math.gcd(N(n), 6) == 1) != coprime6(n) for n in allprim)
    fb = sum(math.gcd(N(n), 6) != 1 for n in allprim)
    report('A3.23', '344-345', '', bad == 0, f'{len(cols6)} columns N<=1e5: N(n) coprime to 6 <=> n coprime to 6',
           falsifier=fb > 0)

    # A3.24  l.348-361: a_m(E_d) = sum_{N(n)=m} n chi_n(u_d) (i.e. a_m/sqrt m = sum alpha (u_d/n)_6), m coprime to 6d;
    # paper: d in {1,2,-3,7}, m<=400; here m <= 20000 and more d
    bad = 0; fb = 0; tot = 0; badp = 0
    for d in [1, 2, -3, 7, -11, 13, 16, -432]:
        an = ellan(d, X); ud = ((4*d)**5, 0)
        for m in range(1, X + 1):
            if math.gcd(m, 6*d) > 1:
                continue
            tot += 1
            r = row_coeff(ud, m, bn, with_mu=False)
            bad += r != (an[m], 0)
            if m <= 400 and d in (1, 2, -3, 7):
                badp += r != (an[m], 0)
            if m <= 3000:
                fb += row_coeff((4*d, 0), m, bn, with_mu=False) != (an[m], 0)
    report('A3.24', '348-361', '', bad == 0,
           f'{tot} coefficients (d in 1,2,-3,7,-11,13,16,-432; m<={X}, (m,6d)=1) exact vs PARI ellan '
           f'(paper\'s range: {badp} mismatches)', falsifier=fb > 0)

    # A3.25  l.358-359: Frobenius eigenvalues pi (u_d/pi)_6 and conjugate: sum a_p, product p
    bad = 0; fb = 0
    for d in [1, 2, -3, 7]:
        an = ellan(d, 5000)
        for p in primes_upto(5000):
            if p % 3 != 1 or (6*d) % p == 0:
                continue
            for pi in primes_above(p):
                bta = mul(pi, UNITS[sym6(((4*d)**5, 0), pi)])
                bad += add(bta, conj(bta)) != (an[p], 0) or mul(bta, conj(bta)) != (p, 0)
                btw = mul(pi, UNITS[sym6((4*d, 0), pi)])
                fb += add(btw, conj(btw)) != (an[p], 0)
    report('A3.25', '358-359', '', bad == 0, 'beta = pi (u_d/pi)_6: beta+conj(beta) = a_p, beta conj(beta) = p (p<5000)',
           falsifier=fb > 0)

    # A3.26  l.362-369: Kummer/Artin interpretation (CITED); spot check: for split pi | p and rational u,
    # (u/pi)_6 = 1 <=> x^6 = u has 6 roots mod p (pi splits completely in K(u^{1/6}))
    bad = 0
    for u in [2, 3, 5, -7, 12]:
        for p in primes_upto(800):
            if p % 3 != 1 or u % p == 0:
                continue
            roots = sum(1 for x_ in range(1, p) if (pow(x_, 6, p) - u) % p == 0)
            for pi in primes_above(p):
                bad += (sym6((u, 0), pi) == 0) != (roots == 6)
    report('A3.26', '362-369', '', bad == 0, 'Kummer spot check: (u/pi)_6=1 <=> x^6=u has 6 roots mod p (u rational, '
           'p<800); Artin/Dedekind factorisation itself CITED', kind='C')

if __name__ == '__main__':
    run()
