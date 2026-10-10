"""Section 6 of the guide: Poisson summation in the row variable and the Moebius-Gauss identities.

Tests (guide displays in brackets):
  6.1   [(6.2)]           conj(chi_{n1}(u)) chi_{n2}(u) = 1_{(u,g)=1} psi(u), psi = conj(chi_{z1}) chi_{z2}.
  6.2   [after (6.2)]     psi is primitive modulo z1 z2 (each local component is nontrivial).
  6.3a  [(6.1)]           the normalization of Phi^: for Phi(t) = exp(-pi t), Phi^(t) = (2/sqrt 3) exp(-4 pi t/3)
                          (two-dimensional quadrature).
  6.3b  [(6.4)]           the exact Poisson formula, five cases including g != 1 and z1 = 1.
  6.5   [(6.7)]           the Moebius-Gauss identity mu_K(n) gamma_{-1}(n) = chi_n(-1) G(n)^{-1} conj(alpha(n)) gamma_2(n)
                          for 70 squarefree n, including composite n.
  6.6   [after (6.7)]     G(n) depends only on n mod 4.
  6.7   [(6.8)]           the paired identity, 22 coprime pairs; a control with the conjugation on the wrong
                          factor is printed and must fail (error of size 1).
  6.8   [before (6.10)]   twisted multiplicativity a(vm) = a(v) a(m) chi_m(v)^4.
  6.9   [(6.10)]          the character identity behind h e^5 v^4 = k f^4.
  6.10  [(3.5) + (6.4)]   end to end: the whole smoothed mean square computed directly equals its Poisson-dual
                          form (D = 40, H = 60).

Run time: about 6 seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import random, mpmath as mp
random.seed(11); allok=True
sq=primaries(400,sqfree=True)

def gcd_primary(a,b):
    fa=factor(a); fb=list(factor(b)); g=(1,0)
    for p in set(fa):
        if p in fb: g=mul(g,p)
    return g

# ---- 6.1 conj(chi_{n1}(u)) chi_{n2}(u) = 1_{(u,g)=1} psi(u), psi = conj(chi_{z1}) chi_{z2}
err=0
for _ in range(300):
    n1,n2=random.choice(sq),random.choice(sq); g=gcd_primary(n1,n2); z1=div(n1,g); z2=div(n2,g)
    u=(random.randint(-40,40),random.randint(-40,40))
    lhs=chi(n1,u).conjugate()*chi(n2,u)
    cop=all(not divides(p,u) for p in factor(g))
    rhs=(1 if cop else 0)*chi(z1,u).conjugate()*chi(z2,u)
    err=max(err,abs(lhs-rhs))
allok&=report("6.1 factorization of the paired character",err,1e-12)

# ---- 6.2 psi primitive mod z1 z2: each local component nontrivial
bad=0; cnt=0
for _ in range(60):
    z1,z2=random.choice(sq[:40]),random.choice(sq[:40])
    if N(gcd_primary(z1,z2))>1: continue
    q=mul(z1,z2)
    for p in set(factor(q)):
        # find x = 1 mod q/p, x varying mod p: CRT via search over x = 1 + (q/p) t
        qp=div(q,p); vals=set()
        for t in residues(p):
            x=add((1,0),mul(qp,t))
            v=chi(z1,x).conjugate()*chi(z2,x)
            if v!=0: vals.add((round(v.real,6),round(v.imag,6)))
        cnt+=1; bad+= (len(vals)<=1)
allok&=report("6.2 psi primitive (every local component nontrivial)",bad,0,f"({cnt} local components)")

# ---- 6.3 hat Phi normalization and the exact Poisson formula
# Phi(t)=exp(-pi t); hatPhi(N w)=(2/sqrt3) int Phi(|z|^2) e(-z w) dx dy = (2/sqrt3) exp(-4 pi N(w)/3)
mp.mp.dps=20
def hatPhi_quad(w):
    f=lambda x,y: mp.e**(-mp.pi*(x*x+y*y))*mp.cos(4*mp.pi*mp.im(mp.mpc(x,y)*w)/mp.sqrt(3))
    return 2/mp.sqrt(3)*mp.quad(f,[-6,6],[-6,6])
err=0
for w in [mp.mpc(0.3,0.2),mp.mpc(-0.7,0.5),mp.mpc(1.1,-0.4)]:
    err=max(err,abs(hatPhi_quad(w)-2/mp.sqrt(3)*mp.e**(-4*mp.pi*abs(w)**2/3)))
allok&=report("6.3a hat Phi of exp(-pi t) equals (2/sqrt3)exp(-4 pi t/3)",float(err),1e-12)
Phi=lambda t: math.exp(-math.pi*t)
hatPhi=lambda t: 2/math.sqrt(3)*math.exp(-4*math.pi*t/3)
def gauss_psi(z1,z2):
    q=mul(z1,z2); qq=cx(q); s=0
    for x in residues(q):
        v=chi(z1,x).conjugate()*chi(z2,x)
        if v!=0: s+=v*e_add(cx(x)/qq)
    return s/math.sqrt(N(q))
def lattice(Rn):
    R=int(2*math.sqrt(Rn))+2
    return [(a,b) for a in range(-R,R+1) for b in range(-R,R+1) if N((a,b))<=Rn]
def divisors_sqfree(g):
    ps=list(set(factor(g))); out=[]
    for mask in range(1<<len(ps)):
        d=(1,0)
        for i,p in enumerate(ps):
            if mask>>i&1: d=mul(d,p)
        out.append(d)
    return out
err=0; cases=0
for (z1,z2,g,H) in [((7,9),(1,3),(1,0),40.0),((-2,-3),(4,3),(-5,0),25.0),((1,0),(7,9),(1,3),30.0),
                    ((-5,0),(1,0),(-2,-3),50.0),((1,3),(1,-3),(1,0),35.0)]:
    q=mul(z1,z2); Nq=N(q)
    lhs=0
    for u in lattice(12*H):
        if any(divides(p,u) for p in factor(g)): continue
        lhs+=Phi(N(u)/H)*chi(z1,u).conjugate()*chi(z2,u)
    gam=gauss_psi(z1,z2) if Nq>1 else 1.0
    rhs=0
    for d in divisors_sqfree(g):
        psid=chi(z1,d).conjugate()*chi(z2,d)
        inner=0
        for h in lattice(12*N(d)*Nq/H+2):
            ph=(chi(z1,h).conjugate()*chi(z2,h)).conjugate() if Nq>1 else 1.0
            inner+=ph*hatPhi(H*N(h)/(N(d)*Nq))
        rhs+=mu(d)*psid/N(d)*inner
    rhs*=H*gam/math.sqrt(Nq)
    err=max(err,abs(lhs-rhs)/max(1,abs(lhs))); cases+=1
allok&=report("6.3b exact Poisson formula (6.4)",err,1e-9,f"({cases} cases incl. g>1 and z1=1)")

# ---- 6.5 Moebius-Gauss identity mu gamma_{-1} = chi_n(-1) G^{-1} conj(alpha) gamma_2
err=0
for n in sq[:70]:
    err=max(err,abs(mu(n)*gamma_j(-1,n)-chi(n,(-1,0))/G(n)*alpha(n).conjugate()*gamma_j(2,n)))
allok&=report("6.5 Moebius-Gauss identity (6.7), 70 squarefree n",err,1e-12)

# ---- 6.6 G depends only on n mod 4 (hence mod 72)
tab={}; bad=0
for n in sq[:90]:
    k=canon_mod(n,(4,0)); v=G(n)
    if k in tab and abs(tab[k]-v)>1e-9: bad+=1
    tab.setdefault(k,v)
allok&=report("6.6 G(n) depends only on n mod 4",bad,0,f"({len(tab)} classes)")
def Gclass(t):
    return tab[canon_mod(t,(4,0))]
def inv_mod4(z):
    for t in residues((4,0)):
        if canon_mod(sub(mul(z,t),(1,0)),(4,0))==(0,0): return t

# ---- 6.7 paired identity: mu mu gamma(conj chi_z1 chi_z2) = a(z1) conj(a(z2)) G(z2 z1^{-1}), a=conj(alpha)gamma_2
err=0; cnt=0
for _ in range(80):
    z1,z2=random.choice(sq[:45]),random.choice(sq[:45])
    if N(gcd_primary(z1,z2))>1 or N(mul(z1,z2))>2500: continue
    lhs=mu(z1)*mu(z2)*(gauss_psi(z1,z2) if N(mul(z1,z2))>1 else 1)
    a1=alpha(z1).conjugate()*gamma_j(2,z1) if N(z1)>1 else 1
    a2=alpha(z2).conjugate()*gamma_j(2,z2) if N(z2)>1 else 1
    t=mul(z2,inv_mod4(z1))
    rhs=a1*a2.conjugate()*Gclass(t)
    err=max(err,abs(lhs-rhs)); cnt+=1
allok&=report("6.7 paired identity (6.8)",err,1e-11,f"({cnt} coprime pairs)")
# control: conjugation on the wrong factor must fail
z1,z2=(7,9),(1,3); t=mul(z2,inv_mod4(z1))
ctrl=abs(mu(z1)*mu(z2)*gauss_psi(z1,z2)-(alpha(z1).conjugate()*gamma_j(2,z1)).conjugate()*alpha(z2).conjugate()*gamma_j(2,z2)*Gclass(t))
print(f"      6.7 [guide (6.8)] control (conjugation swapped): error {ctrl:.3f} (should be O(1))")

# ---- 6.8 twisted multiplicativity a(vm)=a(v)a(m)chi_m(v)^4
err=0; cnt=0
for _ in range(60):
    v,m=random.choice(sq[:40]),random.choice(sq[:40])
    if N(gcd_primary(v,m))>1 or N(mul(v,m))>2500 or v==(1,0) or m==(1,0): continue
    a=lambda n: alpha(n).conjugate()*gamma_j(2,n)
    err=max(err,abs(a(mul(v,m))-a(v)*a(m)*chi(m,v)**4)); cnt+=1
allok&=report("6.8 a(vm)=a(v)a(m)chi_m(v)^4",err,1e-11,f"({cnt} pairs)")

# ---- 6.9 chi_n(h e^5) chi_n(v)^4 = chi_n(k) chi_n(f)^4 with f=ev, k=eh
err=0
for _ in range(200):
    n=random.choice(sq); e=random.choice(sq[:20]); v=random.choice(sq[:20]); h=(random.randint(-20,20),random.randint(-20,20))
    f=mul(e,v); k=mul(e,h)
    err=max(err,abs(chi(n,mul(h,pw(e,5)))*chi(n,v)**4-chi(n,k)*chi(n,f)**4))
allok&=report("6.9 h e^5 v^4 = k f^4 at the level of characters",err,1e-12)

# ---- 6.10 end to end: smoothed mean square = sum over (n1,n2) of the Poisson formula
def Wb(y): return math.exp(-1/((y-1)*(2-y))) if 1<y<2 else 0.0
D=40; H=60.0
cols=[n for n in primaries(2*D,sqfree=True) if N(n)>D]
direct=0
for u in lattice(14*H):
    Au=sum(mu(n)*chi(n,u)*Wb(N(n)/D) for n in cols)
    direct+=Phi(N(u)/H)*abs(Au)**2
dual=0
gcache={}
for n1 in cols:
    for n2 in cols:
        g=gcd_primary(n1,n2); z1=div(n1,g); z2=div(n2,g); q=mul(z1,z2); Nq=N(q)
        key=(z1,z2)
        if key not in gcache: gcache[key]=gauss_psi(z1,z2) if Nq>1 else 1.0
        tot=0
        for d in divisors_sqfree(g):
            psid=chi(z1,d).conjugate()*chi(z2,d)
            inner=sum(((chi(z1,h).conjugate()*chi(z2,h)).conjugate() if Nq>1 else 1.0)*hatPhi(H*N(h)/(N(d)*Nq))
                      for h in lattice(14*N(d)*Nq/H+2))
            tot+=mu(d)*psid/N(d)*inner
        dual+=mu(n1)*mu(n2)*Wb(N(n1)/D)*Wb(N(n2)/D)*H*gcache[key]/math.sqrt(Nq)*tot
err=abs(direct-dual)/abs(direct)
allok&=report("6.10 smoothed mean square: direct vs. dual (Poisson) evaluation",err,1e-9,f"(D={D}, H={H:.0f}, {len(cols)} columns, M={direct:.6f})")
print("ALL PASS" if allok else "SOME FAILURES")
