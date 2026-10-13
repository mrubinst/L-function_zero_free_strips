"""Part A, Section 'The average, in coordinates and as a double sum' (lines 372-455)."""
from cmcore import *
import sympy as sp
random.seed(4)

def Rcoef(u, m, bn):
    """sqrt(m) q_u(m) = sum_{N(n)=m} mu_K(n) n chi_n(u), exactly in Z[omega]."""
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

def M_direct(D, H):
    """eq:coordinates evaluated literally, with per-element Euler-criterion symbols (zw.sym6_prime)."""
    xs, ys = lattice_points(H)
    cols = [n for n in columns(2*D) if muK(n) != 0]
    tot = 0.0
    for x, y in zip(xs.tolist(), ys.tolist()):
        s = 0j
        for (A, B) in cols:
            Q = A*A - A*B + B*B
            w = W_test(np.array([Q/D]))[0]
            if w == 0:
                continue
            j = 0; zero = False
            for pi in factor_primary((A, B)):
                jj = sym6_prime(canon_mod((x, y), pi), pi)
                if jj is None:
                    zero = True; break
                j += jj
            if zero:
                continue
            s += muK((A, B))*(A + B*cmath.exp(2j*math.pi/3))/math.sqrt(Q)*ZETA[j % 6]*w
        tot += abs(s)**2
    return tot

def M_double(D, H, variant=0):
    """eq:double: sum_{n1,n2} mu mu conj(alpha1) alpha2 conj(W1) W2 S_H(n1,n2), S_H from the vectorised tables."""
    xs, ys = lattice_points(H)
    cols = [n for n in columns(2*D) if muK(n) != 0 and abs(W_test(np.array([N(n)/D]))[0]) > 0]
    C = row_matrix(xs, ys, cols)
    if variant == 2:
        S = C.T @ C                                   # wrong: no conjugate on chi_{n1}
    else:
        S = C.conj().T @ C                            # S_H(n1, n2) = sum_u conj(chi_n1(u)) chi_n2(u)
    mu = np.array([muK(n) for n in cols], float)
    al = np.array([alpha(n) for n in cols])
    Wv = W_test(np.array([N(n)/D for n in cols]))
    a1 = mu*(al if variant == 1 else al.conj())*Wv.conj()     # variant 1: wrong, alpha(n1) not conjugated
    a2 = mu*al*Wv
    return float(np.real(a1 @ S @ a2)), complex(a1 @ S @ a2)

def M_fast(D, H):
    xs, ys = lattice_points(H)
    cols = [n for n in columns(2*D) if muK(n) != 0]
    v = np.array([muK(n)*alpha(n) for n in cols])*W_test(np.array([N(n)/D for n in cols]))
    keep = np.abs(v) > 0
    cols = [c for c, k in zip(cols, keep) if k]; v = v[keep]
    tot = 0.0; step = 20000
    for i in range(0, len(xs), step):
        C = row_matrix(xs[i:i + step], ys[i:i + step], cols)
        tot += float(np.sum(np.abs(C @ v)**2))
    return tot, len(xs), len(cols)

def run():
    # A4.01  l.378: about H rows and about D columns per row
    msg = []
    for H in [10**3, 10**4, 10**5, 10**6]:
        xs, _ = lattice_points(H)
        msg.append(f'H={H}: rows/H={len(xs)/H:.4f}')
    for X in [10**3, 10**4, 10**5]:
        msg.append(f'X={X}: #columns(N<=X)/X={len(columns(X))/X:.4f}')
    info('A4.01', '378', '', '; '.join(msg) + f' (limits 2pi/sqrt3={2*math.pi/SQ3:.4f}, 2pi/sqrt3/12={2*math.pi/SQ3/12:.4f})')

    # A4.02  l.379-391, eq:coordinates: column set and the literal coordinate formula = eq:double
    bad = 0
    for A in range(-200, 201):
        for B in range(-200, 201):
            n = (A, B)
            if n == (0, 0):
                continue
            c1 = (A % 3 == 1 and B % 3 == 0 and not (A % 2 == 0 and B % 2 == 0))
            c2 = primary(n) and coprime6(n)
            bad += c1 != c2
    errs = []; fb1 = []; fb2 = []
    for (D, H) in [(25, 150), (40, 300), (60, 200)]:
        md = M_direct(D, H); mdb, mc = M_double(D, H)
        errs.append(abs(md - mdb)/max(1, md)); errs.append(abs(mc.imag)/max(1, md))
        fb1.append(abs(M_double(D, H, 1)[1] - md)/md); fb2.append(abs(M_double(D, H, 2)[1] - md)/md)
    ok = bad == 0 and max(errs) < 1e-10
    report('A4.02', '379-402', 'eq:coordinates, eq:double', ok,
           f'coordinate column set = primary coprime to 6 on |A|,|B|<=200; literal coordinate M = double-sum M for '
           f'(D,H)=(25,150),(40,300),(60,200): max rel err {max(errs):.1e}; falsifiers (alpha(n1) unconjugated / '
           f'chi_n1 unconjugated): rel err {min(fb1):.2f} / {min(fb2):.2f}', falsifier=min(fb1) > 1e-3 and min(fb2) > 1e-3)

    # A4.03  l.403-405: conj(alpha(n1)) alpha(n2) = conj(n1) n2 / (|n1||n2|)
    cols = columns(3000)
    err = 0.0; fb = 0
    for _ in range(2000):
        n1, n2 = random.choice(cols), random.choice(cols)
        z1, z2 = cx(n1), cx(n2)
        err = max(err, abs(alpha(n1).conjugate()*alpha(n2) - z1.conjugate()*z2/(abs(z1)*abs(z2))))
        fb += abs(alpha(n1).conjugate()*alpha(n2) - z1*z2.conjugate()/(abs(z1)*abs(z2))) > 1e-9
    report('A4.03', '403-405', '', err < 1e-12, f'max err {err:.1e} on 2000 pairs', falsifier=fb > 0)

    # A4.04  l.408-411: A_u(D) = sum_m q_u(m) W(m/D) (regrouping by norm)
    X = 6000; cols = columns(X); bn = by_norm(cols)
    err = 0.0; fbe = 0.0
    for u in [(1, 0), (5, 3), (-7, 2), (16, 0), (2, 9)]:
        for D in [300, 1500, 3000]:
            Au = sum(muK(n)*alpha(n)*chi(n, u)*W_test(np.array([N(n)/D]))[0] for n in cols if N(n) <= 2*D)
            r = 0j
            for m in range(1, 2*D + 1):
                R = Rcoef(u, m, bn)
                if R != (0, 0):
                    r += cx(R)/math.sqrt(m)*W_test(np.array([m/D]))[0]
            err = max(err, abs(Au - r))
            Au_nomu = sum(alpha(n)*chi(n, u)*W_test(np.array([N(n)/D]))[0] for n in cols if N(n) <= 2*D)
            fbe = max(fbe, abs(Au_nomu - r))
    report('A4.04', '408-411', '', err < 1e-9, f'A_u(D) = sum_m q_u(m)W(m/D), 5 rows, D<=3000: max err {err:.1e}',
           falsifier=fbe > 1e-3)

    # A4.05  l.412-425, eq:allrows: q_u multiplicative, 0 at 2|m or 3|m, local values, q_u(p^j)=0 (j>=3),
    # reduces to eq:qd at u=u_d
    X = 12000; cols = columns(X); bn = by_norm(cols)
    US = [(random.randint(-300, 300), random.randint(-300, 300)) for _ in range(25)]
    US += [(35, 0), (7, 21), mul((1, 3), (4, 1)), (11, 0), ((4*7)**5, 0)]
    bad = 0; fb1 = 0; fb2 = 0; fb3 = 0; nchk = 0
    for u in US:
        if u == (0, 0):
            continue
        R = {m: Rcoef(u, m, bn) for m in range(1, X + 1)}
        bad += sum(R[m] != (0, 0) for m in R if math.gcd(m, 6) > 1)
        for p in primes_upto(int(X**.5)):
            if p < 5:
                continue
            if p % 3 == 1:
                pi = primes_above(p)[0]; pb = conj(pi)
                j1, j2 = sym6(u, pi), sym6(u, pb)
                t1 = (0, 0) if j1 is None else mul(pi, UNITS[j1])
                t2 = (0, 0) if j2 is None else mul(pb, UNITS[j2])
                bad += R[p] != neg(add(t1, t2))
                fb1 += R[p] != add(t1, t2)
                e2 = (0, 0) if (j1 is None or j2 is None) else mul((p, 0), UNITS[(j1 + j2) % 6])
                bad += R[p*p] != e2
            else:
                j = sym6(u, (-p, 0))
                e2 = (0, 0) if j is None else mul((p, 0), UNITS[j])
                bad += R[p] != (0, 0) or R[p*p] != e2
                fb2 += R[p*p] != neg(e2)
            for k in (3, 4):
                if p**k <= X:
                    bad += R[p**k] != (0, 0)
            nchk += 1
        for _ in range(300):
            m1, m2 = random.randint(1, 110), random.randint(1, 110)
            if m1*m2 > X:
                continue
            if math.gcd(m1, m2) == 1:
                bad += R[m1*m2] != mul(R[m1], R[m2])
            elif R[m1] != (0, 0) and R[m2] != (0, 0):
                fb3 += R[m1*m2] != mul(R[m1], R[m2])
    # reduction at u = u_d: for p !| 6d, q(p^2) = 1 (split: |chi|^2; inert: chi_{-p}(rational) = 1)
    d = 7; ud = ((4*d)**5, 0); red = 0
    for p in primes_upto(100):
        if p < 5 or p == 7:
            continue
        red += Rcoef(ud, p*p, bn) != (p, 0)
    report('A4.05', '412-425', 'eq:allrows', bad == 0 and red == 0,
           f'{len(US)} rows, {nchk} (u,p) local checks (p<=109), multiplicativity, vanishing at 2,3 and p^3,p^4, '
           f'reduction to eq:qd at u_7; falsifiers (+sign at split p: {fb1}, -chi at inert p^2: {fb2}, '
           f'non-coprime product: {fb3})', falsifier=fb1 > 0 and fb2 > 0 and fb3 > 0)

    # A4.06  l.426-427: for a general row the symbols distinguish pi from pibar
    p = 7; pi = primes_above(p)[0]; pb = conj(pi)
    vals = set(); dist = 0
    for u in [(2, 1), (3, 5), (1, 4), (5, 0), (2, 0)]:
        j1, j2 = sym6(u, pi), sym6(u, pb)
        vals.add(Rcoef(u, p, by_norm(columns(50))))
        dist += (j1 != (-j2) % 6)
    dist_rat = sum(sym6((r, 0), pi) != (-sym6((r, 0), pb)) % 6 for r in [2, 3, 5, 10, -11])
    report('A4.06', '426-427', '', dist > 0 and len(vals) > 1,
           f'p=7: chi_pi(u) != conj chi_pibar(u) for {dist} of 5 rows; q_u(7) takes {len(vals)} values while mu, '
           'chi_-3 are fixed', falsifier=dist_rat == 0)   # detector silent on rational rows

    # A4.07  l.430-437, eq:MS: M(D, D^{1+theta}) << D^{2+theta+eps} -- sanity ratio only
    msg = []
    for D in [100, 200, 400, 800, 1600]:
        H = int(D**1.1)
        M, R_, K = M_fast(D, H)
        msg.append(f'D={D}: M/D^2.1={M/D**2.1:.4f} (rows {R_}, cols {K})')
    info('A4.07', '430-437', 'eq:MS', 'theta=0.1: ' + '; '.join(msg))
    untested('A4.07b', '437-439', '', 'dependence of the constant on finitely many derivatives of W; d-independence '
             '(structural remark, not finitely testable)')

    # A4.08  l.441-443: S_H(n,n) = #{u: 0<N(u)<=H, (u,n)=1} <= #{0<N(u)<=H} << H
    H = 20000; xs, ys = lattice_points(H)
    cols = columns(3000)
    sample = random.sample(cols, 150)
    C = row_matrix(xs, ys, sample)
    S = np.sum(np.abs(C)**2, axis=0)
    bad = 0; fb = 0
    for k, n in enumerate(sample):
        cop = np.ones(len(xs), bool)
        for pi in set(factor_primary(n)):
            cop &= sym_vec(pi, xs, ys) >= 0
        bad += abs(S[k] - cop.sum()) > 1e-6 or S[k] > len(xs)
        fb += abs(np.sum(np.real(C[:, k])) - cop.sum()) > 1e-6        # wrong: sum of chi instead of |chi|^2
    report('A4.08', '441-443', '', bad == 0, f'150 columns, H={H}: S_H(n,n) = #coprime u <= #rows = {len(xs)}',
           falsifier=fb > 0)

    # A4.09  l.445-451: S_H(n,n) = (2pi/sqrt3) H prod(1-1/N(pi)) + O(tau(n)(sqrt H + 1))
    worst = 0.0; worst_noprod = 0.0
    for H in [10**3, 10**4, 10**5, 4*10**5]:
        xs, ys = lattice_points(H)
        for n in random.sample(columns(5000), 60) + [mul((1, 3), mul((-2, -3), (-5, 0))), pw((1, 3), 3)]:
            cop = np.ones(len(xs), bool)
            for pi in set(factor_primary(n)):
                cop &= sym_vec(pi, xs, ys) >= 0
            main = 2*math.pi/SQ3*H*np.prod([1 - 1/N(pi) for pi in set(factor_primary(n))])
            f = factor_primary(n); tau = np.prod([f.count(pi) + 1 for pi in set(f)])
            worst = max(worst, abs(cop.sum() - main)/(tau*(math.sqrt(H) + 1)))
            worst_noprod = max(worst_noprod, abs(cop.sum() - 2*math.pi/SQ3*H)/(tau*(math.sqrt(H) + 1)))
    info('A4.09', '445-451', '', f'max |S_H(n,n) - main|/(tau(n)(sqrt H+1)) = {worst:.3f} over 62 columns x '
         f'H in 1e3..4e5 (bounded; without the Euler product the ratio reaches {worst_noprod:.1f})')

    # A4.10  l.451-452: 2pi/sqrt3 = area of disc / covolume
    cov = sp.Abs(sp.Matrix([[1, 0], [sp.Rational(-1, 2), sp.sqrt(3)/2]]).det())
    ok = sp.simplify(sp.pi/cov - 2*sp.pi/sp.sqrt(3)) == 0
    report('A4.10', '451-452', '', ok, f'covolume of Z[omega] = {cov}; pi H / covolume = 2 pi H/sqrt3', kind='I',
           falsifier=sp.simplify(sp.pi/1 - 2*sp.pi/sp.sqrt(3)) != 0)

    # A4.11  l.452-453: diagonal is O_W(HD) -- sanity
    msg = []
    for D in [100, 400, 1600]:
        H = int(D**1.1); xs, ys = lattice_points(H)
        cols = [n for n in columns(2*D) if muK(n) != 0]
        Wv = W_test(np.array([N(n)/D for n in cols]))
        C = row_matrix(xs, ys, cols)
        diag = float(np.sum(np.abs(Wv)**2*np.sum(np.abs(C)**2, axis=0)))
        msg.append(f'D={D}: diag/(HD)={diag/(H*D):.4f}')
    info('A4.11', '452-453', '', '; '.join(msg))
    untested('A4.12', '453-455', '', 'remark on what (MS) bounds (whole off-diagonal, not termwise) -- not a numerical claim')

if __name__ == '__main__':
    run()
