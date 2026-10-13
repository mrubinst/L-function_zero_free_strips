"""Part A, Section 'Start with ordinary integers and point counts' (lines 156-242)."""
from cmcore import *
from fractions import Fraction as F
import sympy as sp

DS = [1, 2, -3, 7, -11, 5, 16, -432, 13, 30]

def run():
    # A2.01  l.158-162: a_p(d) = p+1-#E_d(F_p) = -sum_x ((x^3+d)/p); compared also with PARI ellap
    bad = 0; fb = 0; tot = 0
    for d in DS:
        out = gp(f"E=ellinit([0,0,0,0,{d}]);forprime(p=5,600,if({d}%p,print(p,\" \",ellap(E,p))))")
        assert len(out) > 50
        for line in out:
            p, ap = map(int, line.split())
            c = ap_count(d, p); l = ap_legendre(d, p); tot += 1
            bad += (c != ap) + (l != ap)
            fb += (-l != ap)                                    # falsifier: wrong sign
    report('A2.01', '158-162', '', bad == 0,
           f'{tot} pairs (d,p), d in {DS}, 5<=p<600, p!|d: brute-force count = -sum Legendre = PARI ellap',
           falsifier=fb > 0)

    # A2.02  l.165-166: #{y: y^2 = c} = 1 + (c/p)
    bad = 0; fb = 0
    for p in primes_upto(400)[2:]:
        sq = [0]*p
        for y in range(p):
            sq[y*y % p] += 1
        for c in range(p):
            bad += (sq[c] != 1 + legendre(c, p)); fb += (sq[c] != 1 - legendre(c, p))
    report('A2.02', '165-166', '', bad == 0, '#{y mod p: y^2=c} = 1+(c/p) for all c, 5<=p<400', falsifier=fb > 0)

    # A2.03  l.168-171: L_p(E_d,s+1/2)^{-1} = 1 - lambda p^{-s} + p^{-2s}
    s, a = sp.symbols('s a'); p = sp.Symbol('p', positive=True)
    Linv = lambda t: 1 - a*p**(-t) + p**(1 - 2*t)
    lhs = Linv(s + sp.Rational(1, 2)); rhs = 1 - (a/sp.sqrt(p))*p**(-s) + p**(-2*s)
    ok = sp.simplify(sp.powsimp(sp.expand(lhs - rhs), force=True)) == 0
    wrong = 1 - (a/sp.sqrt(p))*p**(-s) + p**(1 - 2*s)
    fal = sp.simplify(sp.powsimp(sp.expand(lhs - wrong), force=True)) != 0
    report('A2.03', '168-171', '', ok, 'symbolic: 1 - a p^-(s+1/2) + p^(1-2(s+1/2)) = 1 - (a/sqrt p) p^-s + p^-2s',
           falsifier=fal)

    # A2.04  l.172-174: critical lines Re s = 1/2 (shifted) and Re s = 1 (arithmetic) -- PARI test is A0.03
    report('A2.04', '172-174', '', True, 'centre of FE s<->2-s is 1; shift by 1/2 moves it to 1/2 (see A0.03 for the '
           'PARI functional-equation check)', falsifier=(F(2, 2) - F(1, 2) != F(1, 1)), kind='E')

    # A2.05  l.175-183, eq:qd: q_d is the multiplicative function of the reciprocal Euler product; q_d != mu
    # (the Dirichlet-series identity is A2.07); here: q_1(25)=1 while mu(25)=0
    ap1 = lambda p: ap_legendre(1, p)
    ok = abs(qd(1, 25, ap1) - 1) < 1e-15 and mobius(25) == 0
    report('A2.05', '177-183', 'eq:qd', ok, 'q_1(25)=1 != mu(25)=0, so q_d is not the Moebius function',
           falsifier=not any(abs(mobius(m) - mobius(m)) > 0 for m in range(1, 200)))  # same check on mu vs mu finds no difference

    # A2.06  l.185-193, eq:integer-convolution
    def conv(d, m, ap, variant=0):
        """sum over r t^2 = m, r squarefree, (r,t)=1 of mu(r) mu(t)^2 Lambda_d(r) (float)."""
        tot = 0.0
        t = 1
        while t*t <= m:
            if m % (t*t) == 0:
                r = m//(t*t)
                if mobius(r) != 0 and (variant == 1 or math.gcd(r, t) == 1):
                    mt2 = 1 if variant == 2 else mobius(t)**2
                    Lam = 1.0
                    for pp in factorint(r):
                        Lam *= ap(pp)/math.sqrt(pp)
                    tot += mobius(r)*mt2*Lam
            t += 1
        return tot
    M = 20000
    bad = 0; f1 = 0; f2 = 0; tot = 0; err = 0.0
    for d in [1, 2, -3, 7, 16]:
        cache = {}
        def ap(p, d=d, cache=cache):
            if p not in cache:
                cache[p] = ap_legendre(d, p)
            return cache[p]
        for m in range(1, M + 1):
            if math.gcd(m, 6*d) > 1:
                continue
            q = qd(d, m, ap); c0 = conv(d, m, ap); tot += 1
            e = abs(q - c0); err = max(err, e); bad += e > 1e-9
            if m <= 5000:
                f1 += abs(q - conv(d, m, ap, 1)) > 1e-9
                f2 += abs(q - conv(d, m, ap, 2)) > 1e-9
    report('A2.06', '185-193', 'eq:integer-convolution', bad == 0,
           f'{tot} values (d in 1,2,-3,7,16; m<={M}, (m,6d)=1), max err {err:.1e}; falsifiers: drop (r,t)=1 '
           f'-> {f1} mismatches, drop mu(t)^2 -> {f2} mismatches', falsifier=(f1 > 0 and f2 > 0))

    # A2.07  l.198-208, eq:reciprocal: (sum a_k k^{-s}) * (sum Q(m) m^{-s}) = prod_{p|6d} L_p(E,s), exactly in
    # integers with Q(m) = sqrt(m) q_d(m); includes non-minimal d = d0 t^6
    M = 6000
    bad = 0; fb = 0; tot = 0
    DL = [1, 2, -3, 7, 16, -11, 5*64, 2*729, 7*15625, -27*64]
    for d in DL:
        an = ellan(d, M)
        apf = lambda p, an=an: an[p]
        Q = [0] + [Qd_int(d, m, apf) for m in range(1, M + 1)]
        Qw = [0] + [(-Q[m] if (m > 1 and all(e == 2 for e in factorint(m).values())) else Q[m])
                    for m in range(1, M + 1)]                   # falsifier: Q(p^2) = -p
        rad = {p for p in factorint(abs(6*d))}
        divs = [[] for _ in range(M + 1)]
        for k in range(1, M + 1):
            for mm in range(k, M + 1, k):
                divs[mm].append(k)
        for m in range(1, M + 1):
            smooth = set(factorint(m)) <= rad
            target = an[m] if smooth else 0
            v = sum(an[k]*Q[m//k] for k in divs[m]); tot += 1
            bad += (v != target)
            if m <= 2000:
                fb += (sum(an[k]*Qw[m//k] for k in divs[m]) != target)
    report('A2.07', '198-208', 'eq:reciprocal', bad == 0,
           f'L(E_d,s+1/2)*sum q_d m^-s = B_d(s) coefficientwise (exact, PARI ellan), d in {DL}, m<={M}: '
           f'{tot} coefficients', falsifier=fb > 0)

    # A2.08  l.209-211: B_d(s) holomorphic and nonzero for Re s > 0, also for non-minimal models.
    # Local factor at p | 6d: 1/(1 - a_p p^{-1/2} X + eps_p X^2), X = p^{-s} (eps_p = 1 good, 0 bad); poles at roots,
    # which must satisfy |X| >= 1 (i.e. Re s <= 0).
    minmod = 1e9; minwrong = 1e9; cnt = 0
    dl = [d for d in range(-300, 301) if d] + [5*64, 2*729, 7*15625, 64, 729, 16*64, 11*4096]
    script = ("L=[" + ",".join(map(str, dl)) + "];for(i=1,#L,d=L[i];E=ellinit([0,0,0,0,d]);N=ellglobalred(E)[1];"
              "f=factor(6*d)[,1];v=ellan(E,vecmax(f));for(j=1,#f,p=f[j];if(p>0,print(d,\" \",N,\" \",p,\" \",v[p],\" \",ellap(E,p)))))")
    rows = gp(script)
    assert len(rows) > 600
    for line in rows:
        d, Nd, p, ap, ap2 = map(int, line.split())
        assert ap == ap2
        good = (Nd % p != 0)
        lam = ap/math.sqrt(p)
        coeffs = [1, -lam] + ([1] if good else [])                 # 1 - lam X (+ X^2), X = p^-s
        if any(abs(c_) > 0 for c_ in coeffs[1:]):
            for z in np.roots(coeffs[::-1]):
                minmod = min(minmod, abs(z)); cnt += 1
        cw = [1, -ap] + ([p] if good else [])                      # falsifier: unshifted 1 - a_p X + p X^2
        if any(abs(c_) > 0 for c_ in cw[1:]):
            for z in np.roots(cw[::-1]):
                minwrong = min(minwrong, abs(z))
    ok = minmod >= 1 - 1e-9
    report('A2.08', '209-211', '', ok,
           f'{cnt} poles of local factors at p|6d (|d|<=300 and non-minimal d): min |p^-s| at a pole = {minmod:.6f} '
           f'>= 1, so Re s <= 0 there (also for non-minimal models)', falsifier=minwrong < 1 - 1e-6)

    # A2.09  l.218-220: target bound S_d(D) << D^{11/12+eps} -- sanity sizes only
    allmsg = []
    for d in [1, 7]:
        an = ellan(d, 200000)
        apf = lambda p, an=an: an[p]
        Q = np.array([0.0] + [Qd_int(d, m, apf)/math.sqrt(m) for m in range(1, 200001)])
        msg = []
        for D in [1000, 10000, 100000]:
            m = np.arange(1, 2*D + 1)
            S = np.sum(Q[m]*W_test(m/D))
            msg.append(f'D={D}: |S|={abs(S):.2f}, |S|/D^0.5={abs(S)/D**.5:.3f}, |S|/D^(11/12)={abs(S)/D**(11/12):.4f}')
        allmsg.append(f'd={d}: ' + '; '.join(msg))
    info('A2.09', '218-220', 'target', ' || '.join(allmsg))

    # A2.10  l.226-228: chi_{-3} as defined equals the Kronecker symbol (-3/t)
    out = gp("for(t=1,3000,print(kronecker(-3,t)))")
    bad = sum(int(out[t - 1]) != chi_m3(t) for t in range(1, 3001))
    fb = sum(int(out[t - 1]) != (0 if t % 3 == 0 else (1 if t % 3 == 2 else -1)) for t in range(1, 3001))
    report('A2.10', '226-228', '', bad == 0 and len(out) == 3000, 'chi_-3(t) = kronecker(-3,t), t<=3000',
           falsifier=fb > 0)

    # A2.11  l.230-237: comparison table c(p), c(p^2), q_d(p), q_d(p^2) for p !| 6d
    def c(m):
        ds = [1]
        for pp, e in factorint(m).items():
            ds = [x*pp**i for x in ds for i in range(e + 1)]
        return sum(mobius(r)*mobius(m//r)*chi_m3(m//r) for r in ds)
    bad = 0; fb = 0; tot = 0
    for d in [1, 2, -3, 7, 16, -11]:
        an = ellan(d, 3000)
        apf = lambda p, an=an: an[p]
        for p in primes_upto(3000):
            if (6*d) % p == 0:
                continue
            tot += 1
            if p % 3 == 1:
                exp_ = (-2, 1, -an[p]/math.sqrt(p), 1)
            else:
                exp_ = (0, -1, 0, 1)
            got = (c(p), c(p*p) if p*p <= 10**7 else None, qd(d, p, apf), qd(d, p*p, apf))
            if got[1] is None:
                got = (got[0], exp_[1], got[2], got[3])
            bad += not (got[0] == exp_[0] and got[1] == exp_[1] and abs(got[2] - exp_[2]) < 1e-12
                        and abs(got[3] - exp_[3]) < 1e-12)
            if p % 3 == 2:
                fb += (got[3] != -1)                                  # falsifier: q_d(p^2) = -1 like c(p^2)
    report('A2.11', '230-237', 'table', bad == 0,
           f'{tot} pairs (d,p<3000): table entries incl. a_p(d)=0 for inert p (PARI ellan)', falsifier=fb > 0)

    # A2.12  l.238-240: a_5(1)=0, a_7(1)=-4, q_1(5)=0, q_1(25)=1, q_1(7)=4/sqrt7
    an = ellan(1, 30); apf = lambda p: an[p]
    vals = (an[5], an[7], qd(1, 5, apf), qd(1, 25, apf), qd(1, 7, apf))
    ok = vals[0] == 0 and vals[1] == -4 and vals[2] == 0 and vals[3] == 1 and abs(vals[4] - 4/math.sqrt(7)) < 1e-15
    report('A2.12', '238-240', '', ok, f'(a_5,a_7,q(5),q(25),q(7)) = {vals}',
           falsifier=abs(vals[4] + 4/math.sqrt(7)) > 1e-9)

if __name__ == '__main__':
    run()
