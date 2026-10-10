"""Sections 4 and 5 of the guide: the diagonal, sixth-power rows, cube and square rows, rational rows,
and the passage from the mean square to the zero-free region.

Tests (guide displays in brackets):
  4.1   [(4.1), (4.2)]  lattice count S_H(n, n) against (2 pi/sqrt 3) H prod(1 - 1/N(pi)); the error divided
                        by tau(n)(sqrt H + 1) stays bounded (four n, H <= 10^4).
  4.2   [(4.3)]         chi_n(v^6) = 1_{(v,n)=1}.
  4.3   [(4.5)]         A_{pi^6}(D) - A_1(D) = -sum_{pi | n} mu_K(n) W(N(n)/D), and its size is << D/N(pi).
  4.4   [after (4.5)]   printed: J log(Y)/Y for Y = 10^3, 10^4, 10^5 (prime counting, tends to 1/2);
                        the rows pi^6 are distinct.
  4.5   [(4.7)]         chi_n(v^3) = (v/n)_2 and chi_n(v^2) = (v/n)_3, with the quadratic and cubic symbols
                        computed independently by Euler's criterion.
  4.6   [(4.8)]         (v/n)_2 = (v/N(n)) (Jacobi symbol) for rational v.
  4.7a  [(2.3), (3.7)]  c_1(m) = c(m) for m <= 2000.
  4.7b  [(4.9)]         c_{v^3}(m) = c(m)(v/m) for m <= 2000 and eight rational v.
  4.8   [(4.10)]        the Euler product prod_{p not | 6}(1 - chi_v(p)p^{-s})(1 - chi_v(p)chi_{-3}(p)p^{-s})
                        has coefficients c(m)(v/m).
  4.9   [after (4.10)]  printed: the counts H^{1/2} (rational rows) and H^{1/12} (rational sixth powers).
  5.1   [(5.1), (5.2)]  the exponent algebra 11/6 + 5 theta/6 and 5/3 - theta/3 (sympy).
  5.2   [(1.5), (5.4)]  sum_{(m,6)=1} c(m) m^{-s} = 1/F(s) at s = 3 and 2.5 + 4i, up to truncation error.
  5.3   [(5.5)]         the substitution int W(m/D) D^{-s} dD/D = m^{-s} W^(s) (30-digit quadrature).

Run time: a few seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import random, mpmath as mp
random.seed(7); allok=True

# ---- 4.1 diagonal count S_H(n,n) = (2pi/sqrt3) H prod(1-1/N(pi)) + O(tau(n)(sqrt H+1))
def count_coprime(H,n):
    R=int(2*math.sqrt(H))+2; c=0
    ps=set(factor(n))
    for a in range(-R,R+1):
        for b in range(-R,R+1):
            u=(a,b)
            if 0<N(u)<=H and all(not divides(p,u) for p in ps): c+=1
    return c
worst=0
for n in [(7,9),(-5,0),mul((7,9),(-5,0)),mul(mul((-2,-3),(1,3)),(-5,0))]:
    ps=set(factor(n)); tau=2**len(ps)
    for H in [100,1000,10000]:
        main=2*math.pi/math.sqrt(3)*H*np.prod([1-1/N(p) for p in ps])
        c=count_coprime(H,n); worst=max(worst,abs(c-main)/(tau*(math.sqrt(H)+1)))
allok&=report("4.1 S_H(n,n): |count-main|/(tau(n)(sqrt H+1))",worst,3.0)

# ---- 4.2 chi_n(v^6)=1_{(v,n)=1}
ns=primaries(400); err=0
for _ in range(400):
    n=random.choice(ns); v=(random.randint(-30,30),random.randint(-30,30))
    if v==(0,0): continue
    lhs=chi(n,pw(v,6)); rhs=1 if all(not divides(p,v) for p in factor(n)) else 0
    err=max(err,abs(lhs-rhs))
allok&=report("4.2 chi_n(v^6)=1_{(v,n)=1}",err,1e-12)

# ---- 4.3 A_{pi^6}(D)-A_1(D) = -sum_{pi|n} mu(n)W(N(n)/D), and size << D/N(pi)
def W(y): return math.exp(-1/((y-1)*(2-y))) if 1<y<2 else 0.0
D=1500
ns2=[n for n in primaries(2*D,sqfree=True) if N(n)>D]
def A(u): return sum(mu(n)*chi(n,u)*W(N(n)/D) for n in ns2)
A1=sum(mu(n)*W(N(n)/D) for n in ns2)
err=0; ratio=0
for p in [q for q in PR_ODD if 40<N(q)<=200][:8]:
    d=A(pw(p,6))-A1
    ident=-sum(mu(n)*W(N(n)/D) for n in ns2 if divides(p,n))
    err=max(err,abs(d-ident)); ratio=max(ratio,abs(d)/(D/N(p)))
allok&=report("4.3 A_{pi^6}-A_1 identity",err,1e-9,f"; max |A_pi6-A_1|/(D/N(pi)) = {ratio:.3f} (bounded, as claimed)")

# ---- 4.4 J = #{primary pi, Y/2<N(pi)<=Y} ~ Y/(2 log Y); rows pi^6 distinct
for Y in [10**3,10**4,10**5]:
    J=sum(1 for p in PR_ODD if Y/2<N(p)<=Y)
    print(f"      4.4 [guide text after (4.5)] Y={Y}: J={J}, J*log(Y)/Y={J*math.log(Y)/Y:.3f} (tends to 1/2)")
six=[pw(p,6) for p in PR_ODD if 1000/2<N(p)<=1000]
allok&=report("4.4 distinct rows pi^6",0 if len(set(six))==len(six) else 1,0)

# ---- 4.5 chi_n(v^3)=(v/n)_2 and chi_n(v^2)=(v/n)_3, symbols defined by Euler's criterion
def symk(v,p,k):    # (v/p)_k: the k-th root of unity congruent to v^{(N p-1)/k} mod p
    if divides(p,v): return 0
    r=powmod(canon_mod(v,p),(N(p)-1)//k,p)
    for j in range(6):
        u=UNITS[j]
        if divides(p,sub(r,u)): return ZETA[j]
    raise ValueError
err=0
for _ in range(300):
    n=random.choice(ns); v=(random.randint(-40,40),random.randint(-40,40))
    if v==(0,0): continue
    q2=np.prod([symk(v,p,2) for p in factor(n)]); q3=np.prod([symk(v,p,3) for p in factor(n)])
    err=max(err,abs(chi(n,pw(v,3))-q2),abs(chi(n,pw(v,2))-q3))
allok&=report("4.5 chi_n(v^3)=(v/n)_2, chi_n(v^2)=(v/n)_3",err,1e-12)

# ---- 4.6 rational v, n primary coprime to 6: (v/n)_2 = Jacobi (v / N(n))
def jacobi(a,m):
    assert m>0 and m%2==1
    a%=m; r=1
    while a:
        while a%2==0:
            a//=2
            if m%8 in (3,5): r=-r
        a,m=m,a
        if a%4==3 and m%4==3: r=-r
        a%=m
    return r if m==1 else 0
err=0
for n in ns[:150]:
    for v in range(-30,31):
        if v==0: continue
        err=max(err,abs(chi(n,(v,0))**3-jacobi(v,N(n))))
allok&=report("4.6 (v/n)_2=(v/N(n)) for rational v",err,1e-12)

# ---- 4.7 c_{v^3}(m)=c(m)(v/m) for (m,6)=1, with c=mu*(mu chi_{-3}) and c_u(m)=sum_{N n=m} mu chi_n(u)
MM=2000
def mobius(m):
    r=1; p=2
    while p*p<=m:
        if m%p==0:
            m//=p
            if m%p==0: return 0
            r=-r
        p+=1
    return -r if m>1 else r
def chi3(m): return 0 if m%3==0 else (1 if m%3==1 else -1)
c=[0]*(MM+1)
for d in range(1,MM+1):
    md=mobius(d)
    if md==0: continue
    for e in range(1,MM//d+1):
        c[d*e]+=md*mobius(e)*chi3(e)
byN={}
for n in primaries(MM):
    byN.setdefault(N(n),[]).append(n)
err0=0; err=0
for m in range(1,MM+1):
    if math.gcd(m,6)!=1: continue
    lst=byN.get(m,[])
    c1=sum(mu(n) for n in lst)
    err0=max(err0,abs(c1-c[m]))
    for v in [2,-1,3,5,-3,7,-6,10]:
        cu=sum(mu(n)*chi(n,(v**3,0)) for n in lst)
        err=max(err,abs(cu-c[m]*jacobi(v,m)))
allok&=report("4.7a c_1(m)=c(m) (2.3), m<=2000",err0,1e-12)
allok&=report("4.7b c_{v^3}(m)=c(m)(v/m), m<=2000, 8 values of v",err,1e-9)

# ---- 4.8 Euler product prod_{p!|6}(1-chi_v(p)p^-s)(1-chi_v(p)chi_{-3}(p)p^-s) has coefficients c(m)(v/m)
err=0
for v in [2,-1,5,-3,7]:
    coef=[0]*(MM+1); coef[1]=1
    for p in range(5,MM+1):
        if not isprime_int(p): continue
        a=jacobi(v,p); b=a*chi3(p)
        # multiply by (1-a x)(1-b x), x=p^-s
        new=coef[:]
        for k in range(MM,0,-1):
            if coef[k]==0: continue
            if k*p<=MM: new[k*p]+=-(a+b)*coef[k]
            if k*p*p<=MM: new[k*p*p]+=a*b*coef[k]
        coef=new
    for m in range(1,MM+1):
        if math.gcd(m,6)==1: err=max(err,abs(coef[m]-c[m]*jacobi(v,m)))
allok&=report("4.8 Euler product of A_{v^3} coefficients",err,0)

# ---- 4.9 counts: rational slice 2 floor(sqrt H); sixth powers
H=10**6
print(f"      4.9 [guide text after (4.10)] H={H}: rational rows {2*int(math.isqrt(H))} (~H^1/2={H**.5:.0f}); "
      f"rational sixth powers {2*int(H**(1/12))} (~H^1/12={H**(1/12):.1f})")

# ---- 5.1 exponent algebra
import sympy as sp
D_,th=sp.symbols('D theta',positive=True); Hs=D_**(1+th)
e1=sp.simplify(sp.log(D_*Hs**sp.Rational(5,6))/sp.log(D_)); e2=sp.simplify(sp.log(D_**2*Hs**sp.Rational(-1,3))/sp.log(D_))
ok=sp.simplify(sp.expand_log(e1,force=True)-(sp.Rational(11,6)+5*th/6))==0 and sp.simplify(sp.expand_log(e2,force=True)-(sp.Rational(5,3)-th/3))==0
allok&=report("5.1 exponents 11/6+5theta/6 and 5/3-theta/3",0 if ok else 1,0)

# ---- 5.2 sum_{(m,6)=1} c(m) m^-s = 1/F(s), F=(1-3^-s)(1-4^-s)zeta(s)L(s,chi_-3)
mp.mp.dps=30
def Lchi3(s): return mp.nsum(lambda k: 1/(3*k+1)**s-1/(3*k+2)**s,[0,mp.inf])
def Ftrunc(s,M):
    return mp.fsum(c[m]*mp.mpf(m)**(-s) for m in range(1,M+1) if math.gcd(m,6)==1)
for s in [mp.mpf(3),mp.mpc(2.5,4)]:
    F=(1-mp.mpf(3)**(-s))*(1-mp.mpf(4)**(-s))*mp.zeta(s)*Lchi3(s)
    val=Ftrunc(s,MM); print(f"      5.2 [guide (5.4)] s={s}: truncated series {mp.nstr(val,12)}  1/F(s) {mp.nstr(1/F,12)}  diff {mp.nstr(abs(val-1/F),3)} (tail ~ M^(1-Re s)={MM**(1-float(mp.re(s))):.1e})")
    allok&=report(f"5.2 Dirichlet series at s={s}",float(abs(val-1/F)),3*MM**(1-float(mp.re(s))))
# ---- 5.3 Mellin substitution int W(m/D) D^-s dD/D = m^-s What(s)
def Wc(y): return mp.e**(-1/((y-1)*(2-y))) if 1<y<2 else 0
s=mp.mpc(1.7,3.1)
What=mp.quad(lambda y: Wc(y)*y**(s-1),[1,1.5,2])
err=0
for m in [5,7,11,25]:
    I=mp.quad(lambda Dd: Wc(m/Dd)*Dd**(-s-1),[m/2,m/1.5,m])
    err=max(err,abs(I-mp.mpf(m)**(-s)*What))
allok&=report("5.3 Mellin substitution",float(err),1e-20)
print("ALL PASS" if allok else "SOME FAILURES")
