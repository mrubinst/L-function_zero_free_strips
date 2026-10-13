"""Agent t1, new test T1-E: the five-member chain eq:b1-euler (3.17), evaluated numerically member by member.

  (1) sum_n mu_K(n) psi_d(n) N(n)^-s            n primary, coprime to 6, N(n) <= X
  (2) prod_pi (1 - psi_d(pi) N(pi)^-s)          pi primary prime not dividing 6, N(pi) <= X
  (3) prod_{p not | 6d} (1 - lambda_d(p) p^-s + p^-2s),  p <= X, lambda_d(p) = a_p/sqrt p (PARI ellap)
  (4) sum_m q_d(m) m^-s,  m <= X  (q_d from eq:qd with PARI a_p)
  (5) B_d(s)/L(E_d, s+1/2)  with L from PARI lfun (no truncation) and B_d = prod_{p | 6d} L_p(E_d, s+1/2) from the
      local data of the minimal model (PARI ellglobalred / ellap)
psi_d(n) = alpha(n) chi_n(u_d), u_d = (4d)^5 (eq:ud), computed with cmcore (Euler's criterion in Z[omega]).
Members (1)-(4) are truncated at X; the tolerance is 20 X^(1-sigma)/(sigma-1) (the tail of sum N(n)^-sigma).
Falsifiers (each must fail by more than 100 times the tolerance for some d):
  F1 u_d = 4d in (1) (sextic symbol not conjugated);  F2 B_d omitted in (5) (fails where B_d != 1, d = 16, -432);
  F3 L(E_d, s) unshifted in (5);  F4 alpha(n) replaced by conj alpha(n) in (2).
Also printed: the value of B_d(s) for each d (an example of the definition eq:d004).
"""
import math, cmath
import numpy as np
from t1lib import check
from cmcore import (columns, muK, alpha, chi, primes_upto, primes_above, gp, N, cx)

X = 20000
DS = [1, 2, -3, 7, 16, -432, 320]       # 320 = 5*2^6 is a non-minimal model of E_5
SS = [3.0 + 0j, 2.5 + 4j]


def gp_data(d, svals):
    """a_p for p <= X, and L(E_d, s+1/2), L(E_d, s), B_d(s) for each s, from PARI/GP."""
    sl = ",".join(f"{s.real}+{s.imag}*I" for s in svals)
    script = f"""\\p 40
E=ellinit([0,0,0,0,{d}]); Nc=ellglobalred(E)[1]; L=lfuncreate(E);
forprime(p=2,{X},print("AP ",p," ",ellap(E,p)));
SV=[{sl}];
{{for(k=1,#SV,s=SV[k]; v=lfun(L,s+1/2); w=lfun(L,s); B=1;
  forprime(p=2,max(3,abs({d})),if((6*{d})%p==0, a=ellap(E,p); f=if(Nc%p, 1-a*p^(-s-1/2)+p^(-2*s), 1-a*p^(-s-1/2)); B=B/f));
  print("LV ",k," ",real(v)," ",imag(v)," ",real(w)," ",imag(w)," ",real(B)," ",imag(B)))}}
print("COND ",Nc);
"""
    out = gp(script, timeout=1200)
    ap = {}; vals = {}; cond = None
    for line in out:
        t = line.split()
        if t[0] == 'AP':
            ap[int(t[1])] = int(t[2])
        elif t[0] == 'LV':
            k = int(t[1]) - 1
            f = [float(x.replace('E', 'e')) for x in t[2:]]
            vals[k] = (complex(f[0], f[1]), complex(f[2], f[3]), complex(f[4], f[5]))
        elif t[0] == 'COND':
            cond = int(t[1])
    return ap, vals, cond


cols = [n for n in columns(X) if muK(n) != 0]
plist = primes_upto(X)
pis = [pi for p in plist if p > 3 for pi in primes_above(p) if N(pi) <= X]
worst = 0.0
fals = {'F1': 0.0, 'F2': 0.0, 'F3': 0.0, 'F4': 0.0}
lines = []
for d in DS:
    ud = ((4*d)**5, 0)
    wrong_u = (4*d, 0)
    ap, vals, cond = gp_data(d, SS)
    psi = {n: alpha(n)*chi(n, ud) for n in cols}
    psi_wrong = {n: alpha(n)*chi(n, wrong_u) for n in cols}
    psi_pi = {pi: alpha(pi)*chi(pi, ud) for pi in pis}

    def qd(m):
        if math.gcd(m, 6*d) > 1:
            return 0.0
        r = 1.0
        mm = m
        for p in plist:
            if p*p > mm:
                break
            e = 0
            while mm % p == 0:
                mm //= p; e += 1
            if e == 1:
                r *= -ap[p]/math.sqrt(p)
            elif e == 2:
                r *= 1.0
            elif e >= 3:
                return 0.0
        if mm > 1:
            r *= -ap[mm]/math.sqrt(mm)
        return r
    qv = np.array([0.0] + [qd(m) for m in range(1, X + 1)])
    ms = np.arange(1, X + 1, dtype=float)
    for k, s in enumerate(SS):
        sig = s.real
        tol = 20*X**(1 - sig)/(sig - 1)
        Lsh, Lun, B = vals[k]
        m1 = sum(muK(n)*psi[n]*N(n)**(-s) for n in cols)
        m1w = sum(muK(n)*psi_wrong[n]*N(n)**(-s) for n in cols)
        m2 = 1.0 + 0j
        m2w = 1.0 + 0j
        for pi in pis:
            m2 *= 1 - psi_pi[pi]*N(pi)**(-s)
            m2w *= 1 - alpha(pi).conjugate()*chi(pi, ud)*N(pi)**(-s)          # F4: conj(alpha)
        m3 = 1.0 + 0j
        for p in plist:
            if (6*d) % p == 0:
                continue
            lam = ap[p]/math.sqrt(p)
            m3 *= 1 - lam*p**(-s) + p**(-2*s)
        m4 = complex(np.sum(qv[1:]*ms**(-s)))
        m5 = B/Lsh
        errs = [abs(m1 - m5), abs(m2 - m5), abs(m3 - m5), abs(m4 - m5), abs(m1 - m4), abs(m2 - m3)]
        e = max(errs)/tol
        worst = max(worst, e)
        fals['F1'] = max(fals['F1'], abs(m1w - m5)/tol)
        if d in (16, -432):
            fals['F2'] = max(fals['F2'], abs(1/Lsh - m4)/tol)
        fals['F3'] = max(fals['F3'], abs(B/Lun - m4)/tol)
        fals['F4'] = max(fals['F4'], abs(m2w - m5)/tol)
        lines.append(f"d={d} (N={cond}) s={s}: B_d={B:.6g}, value {m5:.10f}, max err {max(errs):.1e} (tol {tol:.1e})")
for l in lines:
    print('   ', l)
check('T1-E1', 'eq:b1-euler', '3.17', worst <= 1.0,
      f"{len(DS)} values of d (incl. non-minimal 320=5*2^6), s in {{3, 2.5+4i}}, X={X}: all five members agree; "
      f"max err/tol = {worst:.2e}; falsifiers err/tol: u_d=4d {fals['F1']:.1e}, B_d omitted {fals['F2']:.1e}, "
      f"unshifted L {fals['F3']:.1e}, conj alpha in (2) {fals['F4']:.1e}",
      falsifier=all(v > 100 for v in fals.values()))
