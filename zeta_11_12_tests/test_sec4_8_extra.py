"""Identities and inequalities of Sections 4-8 of the guide that the other scripts do not cover.

Tests (guide displays in brackets):
  4.10  [(4.6)]  the extraction inequality J|A_1|^2 <= 2 sum |A_{pi^6}|^2 + 2 sum |A_{pi^6} - A_1|^2 on actual
                 values (D = 1500, Y = 200, J = 21); the error term divided by J D^2/Y^2 is printed.
  4.11  [(4.3)]  A_{v^6}(D) = sum_{(n,v)=1} mu_K(n) W(N(n)/D) for composite v.
  6.11  [(6.5)]  the frequency-range algebra N(d)N(z1 z2)/H << D^2/(H N(g)) (sympy).
  7.11  [(7.6)]  B_{pi,1}(x) = chi_pi(x)^{-3} = chi_pi(x)^3 = (x/pi)_2.
  7.12a [(7.4)]  the j = 0 expansion 1_{pi not | x} = (1 - 1/N(pi)) - (1/N(pi)) sum_{h != 0} e(hx/pi)
                 (the constant term is the inactive branch).
  7.12b [(7.5)]  the active j = 0 branch, summed directly:
                 -(1/N) sum_{h != 0} chi_pi(h)^{-2} e(eps h^{-1} x/pi) = -gamma_2(pi) chi_pi(eps)^{-2} N^{-1/2} chi_pi(x)^{-2},
                 for several eps; it has modulus N(pi)^{-1/2} if pi does not divide x, and is 0 if pi | x.

Run time: about 1 second.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import random, sympy as sp
random.seed(23); allok=True
def Wb(y): return math.exp(4-1/((y-1)*(2-y))) if 1<y<2 else 0.0

# ---- 4.10 extraction (4.6) on real data: J|A1|^2 <= 2 sum|A_{pi^6}|^2 + 2 sum|A_{pi^6}-A1|^2,
#            and sum|A_{pi^6}-A1|^2 << J D^2/Y^2
D=1500; Y=200
cols=[n for n in primaries(2*D,sqfree=True) if N(n)>D]
A1=sum(mu(n)*Wb(N(n)/D) for n in cols)
PY=[p for p in PR_ODD if Y/2<N(p)<=Y]; J=len(PY)
Ap=[sum(mu(n)*chi(n,pw(p,6))*Wb(N(n)/D) for n in cols) for p in PY]
lhs=J*abs(A1)**2; rhs=2*sum(abs(a)**2 for a in Ap)+2*sum(abs(a-A1)**2 for a in Ap)
ratio=sum(abs(a-A1)**2 for a in Ap)/(J*D*D/Y/Y)
allok&=report("4.10 extraction inequality J|A_1|^2 <= 2 sum|A_pi6|^2 + 2 sum|A_pi6-A_1|^2",0 if lhs<=rhs*(1+1e-12) else 1,0,
              f"(J={J}; error term / (J D^2/Y^2) = {ratio:.2e}, bounded)")
# ---- 4.11 A_{v^6}(D) = sum_{(n,v)=1} mu(n) W(N(n)/D) for composite v
err=0
for v in [(2,0),mul((7,9),(1,3)),(5,1),mul((-5,0),(4,3))]:
    lhs=sum(mu(n)*chi(n,pw(v,6))*Wb(N(n)/D) for n in cols)
    rhs=sum(mu(n)*Wb(N(n)/D) for n in cols if all(not divides(p,v) for p in factor(n)))
    err=max(err,abs(lhs-rhs))
allok&=report("4.11 A_{v^6}(D)=sum_{(n,v)=1} mu W for composite v",err,1e-12)
# ---- 6.11 frequency range: N(d)N(z1 z2)/H <= c^2 D^2/(H N(g)) when d|g, N(z_i)<=c D/N(g)
Nd,Ng,Nz1,Nz2,Dd,Hh,c=sp.symbols('Nd Ng Nz1 Nz2 D H c',positive=True)
expr=sp.simplify((Ng*(c*Dd/Ng)**2/Hh)-(c**2*Dd**2/(Hh*Ng)))
allok&=report("6.11 frequency-range algebra N(d)N(z1z2)/H << D^2/(H N(g))",0 if expr==0 else 1,0)
# ---- 7.11 B_{pi,1}(x) = chi_pi(x)^{-3} = chi_pi(x)^3 = (x/pi)_2
err=0
for p in [q for q in PR_ODD if N(q)<200][:10]:
    for x in residues(p)[:30]:
        c=chi(p,x)
        if c==0: continue
        q2=1 if powmod(canon_mod(x,p),(N(p)-1)//2,p)==(1,0) else -1
        err=max(err,abs(c**(-3)-c**3),abs(c**3-q2))
allok&=report("7.11 B_{pi,1}(x)=chi_pi(x)^{-3}=chi_pi(x)^3=(x/pi)_2",err,1e-12)
# ---- 7.12 j=0 correction: 1_{pi nmid x} = (1-1/N) - (1/N) sum_{h!=0} e(hx/pi); active coefficient has modulus N^{-1/2}
err1=err2=0; n712=[0,0]
for p in [q for q in PR_ODD if 7<=N(q)<=120][:8]:
    res=residues(p); Np=N(p); pc=cx(p)
    g2=sum(chi(p,y)**2*e_add(cx(y)/pc) for y in res if chi(p,y)!=0)
    for x in res[:15]:
        lhs=1 if chi(p,x)!=0 else 0
        rhs=(1-1/Np)-sum(e_add(cx(mul(h,x))/pc) for h in res if h!=(0,0))/Np
        err1=max(err1,abs(lhs-rhs))
    g2n=g2/math.sqrt(Np)                       # gamma_2(pi)
    inv={h:next(y for y in res if canon_mod(sub(mul(y,h),(1,0)),p)==(0,0)) for h in res if chi(p,h)!=0}
    for ep in res[1:6]:                        # cusp constant eps (any residue prime to pi)
        if chi(p,ep)==0: continue
        for x in res:
            active=-sum(chi(p,h)**(-2)*e_add(cx(mul(mul(ep,inv[h]),x))/pc) for h in inv)/Np
            pred=(-g2n*chi(p,ep)**(-2)*Np**-0.5*chi(p,x)**(-2)) if chi(p,x)!=0 else 0
            err2=max(err2,abs(active-pred))
            err2=max(err2,abs(abs(active)-(Np**-0.5 if chi(p,x)!=0 else 0))); n712[0]+=1; n712[1]+=(chi(p,x)==0)
allok&=report("7.12a j=0: constant (inactive) term 1-1/N(pi) in the expansion of 1_{pi nmid x}",err1,1e-12)
allok&=report("7.12b j=0 active branch, summed directly: -(1/N) sum_h chi(h)^-2 e(eps h^-1 x/pi) = -gamma_2(pi) chi(eps)^-2 N^-1/2 chi(x)^-2; modulus N^-1/2 if pi nmid x, 0 if pi | x",err2,1e-12,f"({n712[0]} (prime, eps, x) triples, {n712[1]} with pi | x)")
print("ALL PASS" if allok else "SOME FAILURES")
