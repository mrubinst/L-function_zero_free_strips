"""Sections 7 and 8 of the guide: theta coefficients, local factors, the large sieve (sanity check),
cube inversion and the descent inequalities.

Tests (guide displays in brackets):
  7.1   [(7.1)]            d_0(lambda^{-3} n b^3) = 3^{5/2} |b| conj((lambda/n)_3) gamma_2(n), with d_0 from the
                           theta coefficients in theta.py (104 indices).
  7.2   [after (7.1)]      (lambda/n)_3 depends only on n mod 9 (Eisenstein's supplementary law).
  7.3   [before (7.1)]     among indices lambda^{-3} x with x primary, the nonzero coefficients are exactly
                           x = n b^3 with n squarefree (908 indices).
  7.4a  [(7.3), (7.4)]     the additive expansions of chi_pi^j (1 <= j <= 5) and of 1_{pi not | x}.
  7.4b  [(7.3), (7.5)]     the reflected local sums sum_{h != 0} chi_pi(h)^{-j-2} e(eps h^{-1} x/pi): a Gauss sum
                           times chi_pi(eps x)^{-j-2} for j != 4, a Ramanujan sum for j = 4; |g| = sqrt N(pi).
  7.5   [theta input]      automorphy of theta under z -> z + 1 and under the inversion.
  7.9   [after (7.6)]      lambda^4 = 9, and the least nonzero norm on lambda^{-4} Z[omega] is 1/81, attained exactly
                           at the six units divided by 9.
  7.10  [after (7.8)]      chi_n(v^2) chi_n(f)^4 = chi_n(f v^2)^4.
  7.8   [(7.7), sanity]    printed: the quadratic large sieve ratio LHS/((K+U) sum|beta|^2) for K = U = 400 and
                           six coefficient choices.  (7.7) is a bound; the ratio should stay bounded (about 0.25).
  8.0   [(7.1), helper]    the fast cubic Gauss sums of fastgauss.py agree with brute force.
  8.1   [(8.1)]            the cube inversion identity for three pairs (k, f) at X = 4000.
  8.2   [(8.2)-(8.11)]     the inequalities between scale parameters on 200 000 random admissible parameter sets,
                           and the retained-support chain (8.8)-(8.9) on a further 200 000 sets with a retained r.
                           A printed control drops condition (8.8) (keeping only F' <= H') and must violate
                           X'/N(j) >= Sigma F/H in some sets, showing that (8.8) is needed.
  8.3   [(6.11), (8.13)]   the initial scales: (D^2/(H B^2))/(D/B) = D^{-theta}/B (sympy).
  8.4   [(8.7)]            the local cancellation identity.

Run time: about 15 seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import theta as TH
TH.factor_primary=lambda n,plist=None: list(factor(n))   # fast factorization (same result)
from fastgauss import g_prime, pairsym, cubic_prime
import random, sympy as sp
random.seed(3); allok=True
sq=primaries(400,sqfree=True)
lam=LAM
def cub(n,x): return chi(n,x)**2       # (x/n)_3

# ---- 7.1 d_0(lam^-3 n b^3) = 3^{5/2}|b| conj((lam/n)_3) gamma_2(n)   [d_0(l)=conj(tau(-l)), tau from DR's table]
err=0; cnt=0
for n in [(1,0)]+sq[:25]:
    for b in [(1,0),(-2,-3),(1,3),(4,3)]:
        m=mul(n,pw(b,3))
        d0=TH.tau(neg(m)).conjugate()
        g2=gamma_j(2,n) if N(n)>1 else 1.0
        rhs=3**2.5*math.sqrt(N(b))*cub(n,lam).conjugate()*g2
        err=max(err,abs(d0-rhs)/abs(rhs)); cnt+=1
allok&=report("7.1 theta coefficient d_0(lam^-3 n b^3)",err,1e-12,f"({cnt} indices)")

# ---- 7.2 (lam/n)_3 depends only on n mod 9 (Eisenstein's supplementary law)
tab={}; bad=0
for n in primaries(1500):
    k=canon_mod(n,(9,0)); v=cub(n,lam)
    if k in tab and abs(tab[k]-v)>1e-9: bad+=1
    tab.setdefault(k,v)
allok&=report("7.2 (lam/n)_3 depends only on n mod 9",bad,0,f"({len(tab)} classes, all primary n with N(n)<=1500)")

# ---- 7.3 selection: for x primary coprime to 6, the coefficient at l=lam^-3 x is nonzero iff x = n b^3, n squarefree
bad=0; cnt=0
for x in primaries(3000):
    f=factor(x); from collections import Counter
    cntr=Counter(f); ok_form=all(e%3!=2 for e in cntr.values())
    nonzero=abs(TH.tau(neg(x)))>1e-9
    bad+=(ok_form!=nonzero); cnt+=1
allok&=report("7.3 nonzero coefficients with lam^3 l primary are exactly l=lam^-3 n b^3",bad,0,f"({cnt} indices)")

# ---- 7.4 local factors: Fourier expansion of chi_pi^j and the reflected Gauss sums S_j(x)
err1=err2=0
for p in [q for q in PR_ODD if 7<=N(q)<=100][:8]:
    res=residues(p); Np=N(p); pc=cx(p)
    def gs(k):   # g(chi^k) = sum_y chi(y)^k e(y/p)
        return sum(chi(p,y)**k*e_add(cx(y)/pc) for y in res if chi(p,y)!=0)
    for j in range(1,6):
        gj=gs(-j)
        for x in res[:12]:
            lhs=chi(p,x)**j if chi(p,x)!=0 else 0
            rhs=sum(chi(p,h)**(-j)*e_add(cx(mul(h,x))/pc) for h in res if chi(p,h)!=0)/gj
            err1=max(err1,abs(lhs-rhs))
    # j=0: 1_{p nmid x} = (1-1/N) - (1/N) sum_{h!=0} e(hx/p)
    for x in res[:12]:
        lhs=1 if chi(p,x)!=0 else 0
        rhs=(1-1/Np)-sum(e_add(cx(mul(h,x))/pc) for h in res if h!=(0,0))/Np
        err1=max(err1,abs(lhs-rhs))
    eps=(1,1)     # a unit standing for the cusp constant eps_p
    for j in range(6):
        k=j+2
        for x in res[:12]:
            S=sum(chi(p,h)**(-k)*e_add(cx(mul(mul(eps,x),inv:=next(y for y in res if canon_mod(sub(mul(y,h),(1,0)),p)==(0,0))))/pc)
                  for h in res if chi(p,h)!=0)
            if k%6!=0:
                pred=(chi(p,mul(eps,x))**(-k) if chi(p,x)!=0 else 0)*gs(k)
            else:
                pred=-1+Np*(1 if chi(p,x)==0 else 0)
            err2=max(err2,abs(S-pred))
        if k%6!=0: err2=max(err2,abs(abs(gs(k))-math.sqrt(Np)))
allok&=report("7.4a Fourier expansions of chi_pi^j (j=1..5) and of 1_{pi nmid x}",err1,1e-11)
allok&=report("7.4b reflected local sums S_j(x) = chi(eps x)^{-j-2} g(chi^{j+2}) (j!=4), Ramanujan sum (j=4); |g|=sqrt N",err2,1e-10)

# ---- 7.5 automorphy of the theta function built from these coefficients (sanity, two generators)
th=TH.Theta(R=4.5)
errA=0
for (z,v) in [(0.13+0.21j,0.93),(-0.21+0.05j,0.97)]:
    a=th(z,v); b=th(z+1,v)                     # translation by 1
    z2,v2=TH.act(((0,0),(-1,0),(1,0),(0,0)),z,v)     # inversion (0,-1;1,0)
    errA=max(errA,abs(a-b)/abs(a),abs(th(z2,v2)-a)/abs(a))
allok&=report("7.5 theta(gamma w)=theta(w) for z->z+1 and the inversion",errA,1e-8)

# ---- 7.9 nonzero reflected indices have N(l) >= 1/81.  The cusp expansions run over 0 != l in lambda^-4 O
#      (source paper2.tex, around line 1771); lambda^4 = 9, so N(l) = N(m)/81 with m in O \ {0}.
assert pw(lam,4)==(9,0)
L4=cx(lam)**4
norms=[]
for a in range(-30,31):
    for b in range(-30,31):
        if (a,b)==(0,0): continue
        l=cx((a,b))/L4; norms.append((abs(l)**2,(a,b)))
mn=min(t for t,_ in norms); attained=sorted(m for t,m in norms if abs(t-1/81)<1e-12)
allok&=report("7.9 min N(l) over 0 != l in lambda^-4 O (|a|,|b|<=30) equals 1/81, attained exactly at l = unit/9",
              abs(mn-1/81)+(0 if attained==sorted(UNITS) else 1),1e-15,f"(min = {mn:.15g}; 1/81 = {1/81:.15g})")

# ---- 7.10 chi_n(v^2) chi_n(f)^4 = chi_n(f v^2)^4
err=0
for _ in range(200):
    n=random.choice(sq); v=(random.randint(-20,20),random.randint(-20,20)); f=random.choice(sq)
    err=max(err,abs(chi(n,pw(v,2))*chi(n,f)**4-chi(n,mul(f,pw(v,2)))**4))
allok&=report("7.10 chi_n(v^2)chi_n(f)^4=chi_n(fv^2)^4",err,1e-12)

# ---- 7.8 quadratic large sieve: sanity check of the shape (K+U) sum|beta|^2
rows=[k for k in primaries(400,sqfree=True)]
colsU=[]
R=30
for a in range(-R,R+1):
    for b in range(-R,R+1):
        n=(a,b)
        if primary(n) and 0<N(n)<=400 and (not divides(LAM,n)):
            # squarefree (possibly even) column: check via factorization over primes incl. -2
            m=n; ok=True; twos=0
            while divides((2,0),m): m=div(m,(-2,0)); twos+=1
            if twos>1 or mu(m)==0: ok=False
            if ok: colsU.append(n)
K=U=400        # norm bounds in (7.7): N(k) <= K, N(n) <= U
worst=worst_cnt=0
for trial in range(6):
    beta={n:(1.0 if trial==0 else complex(random.gauss(0,1),random.gauss(0,1))) for n in colsU}
    lhs=sum(abs(sum(beta[n]*chi(k,n)**3 for n in colsU))**2 for k in rows)
    nb=sum(abs(b)**2 for b in beta.values())
    worst=max(worst,lhs/((K+U)*nb)); worst_cnt=max(worst_cnt,lhs/((len(rows)+len(colsU))*nb))
print(f"      7.8 [guide (7.7), sanity] large sieve (K=U=400, {len(rows)} rows, {len(colsU)} columns): "
      f"LHS/((K+U) sum|beta|^2) <= {worst:.3f}; LHS/((#rows+#cols) sum|beta|^2) <= {worst_cnt:.3f} over 6 coefficient choices")

# ---- 8.1 cube inversion identity (8.1), with a = conj(alpha) gamma_2 computed via fast cubic Gauss sums
def fact_all(n): return factor(n)
gp={}
def g_cub(n):    # g(n) = sum_y (y/n)_3 ce(y/n) for squarefree primary n coprime to 6
    f=factor(n); g=1+0j
    for p in f:
        if p not in gp: gp[p]=g_prime(p)
        g*=gp[p]
    for i in range(len(f)):
        for j in range(i+1,len(f)): g*=pairsym(f[i],f[j])*pairsym(f[j],f[i])
    return g
def gamma2_fast(n): return cub(n,lam)*g_cub(n)/math.sqrt(N(n)) if N(n)>1 else 1.0
e=0
for n in sq[:40]: e=max(e,abs(gamma2_fast(n)-gamma_j(2,n)))
allok&=report("8.0 gamma_2(n)=(lam/n)_3 g(n)/sqrt N(n) via fast cubic Gauss sums (helper check)",e,1e-10)
def Wb(y): return math.exp(-1/((y-1)*(2-y))) if 1<y<2 else 0.0
X=4000.0
sqX=primaries(int(2*X),sqfree=True)
bs=primaries(int((2*X)**(1/3))+1)       # includes b=1
assert bs.count((1,0))==1
for (k,f) in [((5,3),(7,9)),((2,0),(1,0)),((-1,0),mul((7,9),(1,3)))]:
    Psi=lambda n: chi(n,k)*chi(n,f)**4
    a=lambda n: alpha(n).conjugate()*gamma2_fast(n)
    lhs=sum(a(n)*Psi(n)*Wb(N(n)/X) for n in sqX)/math.sqrt(X)
    def T(Y):
        s=0
        for b in bs:
            Nb3=N(b)**3
            if Nb3>2*Y: continue
            wb=(alpha(b).conjugate()**3*Psi(b)**3*math.sqrt(N(b))) if b!=(1,0) else 1.0
            for n in sqX:
                w=Wb(N(n)*Nb3/Y)
                if w: s+=a(n)*Psi(n)*wb*w
        return s/math.sqrt(Y)
    rhs=0
    for h in bs:
        if mu(h)==0 and h!=(1,0): continue
        Nh3=N(h)**3
        if Nh3>2*X: continue
        mh=mu(h) if h!=(1,0) else 1
        ch=(alpha(h).conjugate()**3*Psi(h)**3) if h!=(1,0) else 1.0
        rhs+=mh*ch/N(h)*T(X/Nh3)
    allok&=report(f"8.1 cube inversion (8.1), k={k}, f={f}",abs(lhs-rhs)/max(abs(lhs),1e-300),1e-12,
                  f"(LHS={complex(lhs):.6g}; {sum(1 for b in bs if N(b)**3<=2*X)-1} nontrivial cubes b)")

# ---- 8.2 descent algebra (random sampling of admissible parameters)
random.seed(5); viol={}; nret=[0]; nctrl=[0]
def chk(name,cond):
    if not cond: viol[name]=viol.get(name,0)+1
for _ in range(200000):
    F=10**random.uniform(0,3); X_=10**random.uniform(0,8); Sig=X_*F
    Hh=Sig*10**random.uniform(-6,0)                 # H <= Sigma
    Hc3=min(X_,X_**2/Hh**2)
    Nh3=Hc3*random.random()                         # N(h)^3 <= H_c^3
    chk("small cubes: H + H^2 F N(h)^3/X <= 2 Sigma", Hh+Hh**2*F*Nh3/X_<=2*Sig*(1+1e-12))
    Nb3=10**random.uniform(0,math.log10(X_)) if X_>1 else 1.0; L=X_/Nb3        # N(b)^3 >= 1, so L <= X
    chk("Y_enl/Y_nat = Sigma/(LF) >= 1", Sig/(L*F)>=1-1e-12)
    R=10**random.uniform(0,3)
    Hp=Hh*L/(Sig*F*R*R); Sp=L/R
    chk("transfer: H' <= HL/(Sigma F)", Hp<=Hh*L/(Sig*F)*(1+1e-12))
    chk("transfer: H'/Sigma' <= H/Sigma", Hp/Sp<=Hh/Sig*(1+1e-12))
    chk("transfer: Sigma' <= L", Sp<=L*(1+1e-12))
    chk("H'_nat = H'*Sigma/(LF) = H/(F^2 R^2)", abs(Hp*Sig/(L*F)-Hh/(F*F*R*R))<=1e-9*Hh/(F*F*R*R))
    if Nb3>X_**2/Hh**2:
        kap=random.uniform(0.01,1); Dk=Sig/Hh       # H = Sigma D^-kappa  means D^kappa = Sigma/H
        chk("contraction: H' < H (H/Sigma)^2", Hh*L/(Sig*F)<Hh*(Hh/Sig)**2*(1+1e-12))
    # control: with only F' <= H' (no N(r)^2 support condition) the bound without the factor 2 can fail
    if Hp>=1:
        Fp=Hp; Xp=L/(R*Fp); Nj=2*R*0.999
        if Xp/Nj<Sig*F/Hh*(1-1e-9): nctrl[0]+=1
# retained divisors r (paper2.tex ~2590): f'|k' and the support of the Fourier parameter give
#   F' <= N(f') <= N(k') <= H L/(Sigma F N(r)^2),   R <= N(r) < 2R,   N(j) <= N(r) for j | r.
# Sample 200000 admissible parameter sets for which some r is retained (i.e. H L/(Sigma F) >= 1).
while nret[0]<200000:
    F=10**random.uniform(0,3); X_=10**random.uniform(0,8); Sig=X_*F
    Hh=Sig*10**random.uniform(-6,0)
    Nb3=10**random.uniform(0,math.log10(X_)) if X_>1 else 1.0; L=X_/Nb3
    top=Hh*L/(Sig*F)
    if top<1: continue
    R2=10**random.uniform(0,0.5*math.log10(top)); Nr=R2*random.uniform(1,2); cap=top/(Nr*Nr)
    if cap<1: continue
    Hp2=top/(R2*R2)
    Fp=10**random.uniform(0,math.log10(cap)); Xp=L/(R2*Fp); Nj=random.uniform(1,Nr)
    chk("H' >= F' >= 1", Hp2>=Fp*(1-1e-12) and Fp>=1)
    chk("X'/N(j) >= Sigma F N(r)^2/(H R N(j))", Xp/Nj>=Sig*F*Nr*Nr/(Hh*R2*Nj)*(1-1e-12))
    chk("X'/N(j) >= Sigma F/H >= 1", Xp/Nj>=Sig*F/Hh*(1-1e-12) and Sig*F/Hh>=1-1e-12)
    nret[0]+=1
allok&=report("8.2 descent inequalities: 200000 random admissible parameter sets, plus 200000 with a retained r",
              sum(viol.values()),0,str(viol) if viol else "")
print(f"      8.2 [guide (8.8)-(8.9)] control: dropping the retained-support condition, X'/N(j) < Sigma F/H in {nctrl[0]} sets (so that condition is needed)")
allok&=nctrl[0]>0
# initial scales: H/Sigma = D^-theta/B
D_,th_,B_=sp.symbols('D theta B',positive=True); Hs=D_**(1+th_)
allok&=report("8.3 initial scales: (D^2/(H B^2))/(D/B) = D^-theta/B",0 if sp.simplify((D_**2/(Hs*B_**2))/(D_/B_)-D_**(-th_)/B_)==0 else 1,0)
# local cancellation (1+1_{pi|k})-1_{pi|k}-1_{pi nmid k} = 1_{pi|k}
allok&=report("8.4 local cancellation identity",max(abs((1+a)-a-(1-a)-a) for a in (0,1)),0)
print("ALL PASS" if allok else "SOME FAILURES")
