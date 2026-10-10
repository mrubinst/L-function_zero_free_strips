"""Sections 1 and 2 of the guide: the integer coefficients c(m), the Eisenstein integers, primary
generators, the splitting of rational primes, and the norm-grouping identity (2.3).

Tests (guide displays in brackets):
  1.1   [(1.1)]  chi_{-3} as defined equals the Kronecker symbol (-3/m) and is completely multiplicative.
  1.2a  [(1.2)]  coefficientwise c * (1 * chi_{-3}) = delta for all m <= 10^5, i.e. the Dirichlet series
                 of c times zeta(s) L(s, chi_{-3}) is 1.
  1.2b  [(1.2)]  the truncated series sum_{m <= 10^5} c(m) m^{-s} against 1/(zeta(s) L(s, chi_{-3})) at
                 s = 3, 3+4i, 4 (mpmath; L(s, chi_{-3}) via Hurwitz zeta).
  1.3a-c         the local table (1.3) for every prime power <= 10^5, multiplicativity of c for all
                 m <= 10^5, and the examples c(7) = -2, c(49) = 1, c(25) = -1.
  1.5a-c [(1.5)] removing 2 and 3: coefficientwise, the Euler factors 1 - 4^{-s} at 2 and 1 - 3^{-s} at 3,
                 and numerically at s = 3+2i.
  2.0-2.2        omega^2 + omega + 1 = 0; N(a + b omega) = a^2 - ab + b^2; exactly six units; "primary"
                 (defined as 3 | n - 1 in Z[omega]) is a = 1, b = 0 mod 3, and each n coprime to 3 has
                 exactly one primary associate.
  2.3a           p = 1 mod 3 splits into two conjugate primary primes of norm p; p = 2 mod 3 has no element
                 of norm p; 3 = -lambda^2 (all p < 3000).
  2.3b-d [(2.3)] #{primary n : N(n) = m} = sum_{d | m} chi_{-3}(d) (unique factorization);
                 sum_{N(n)=m} mu_K(n) = c(m) for all m <= 10^5 coprime to 6; the local factors.
  2.4   [(2.4)]  A_1(D) computed from (1.4) and from (2.4) agree (complex smooth weight).

Independence: c(m) is computed from its definition (1.2) with sympy's ordinary Moebius function;
mu_K(n) is computed from the factorization of N(n) in Z (not from a list of Eisenstein primes); and
"primary" is tested from its definition.  So the two sides of (2.3) are computed independently.

Run time: about 5 seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import random, numpy as np, mpmath as mp
from sympy import mobius as smob, jacobi_symbol, factorint, isprime, primerange
random.seed(11); allok=True

M=100000
# ---------- smallest-prime-factor sieve for fast integer factorizations
spf=np.zeros(M+1,dtype=np.int64)
for i in range(2,M+1):
    if spf[i]==0: spf[i:M+1:i][spf[i:M+1:i]==0]=i
def ifac(m):
    out={}
    while m>1:
        p=int(spf[m]); out[p]=out.get(p,0)+1; m//=p
    return out

# ---------- (1.1) chi_{-3}
def chi3(m): return 0 if m%3==0 else (1 if m%3==1 else -1)
bad=0
for m in range(1,3001):
    k=jacobi_symbol(m,3) if m%3 else 0          # (m/3) = (-3/m) for m>0 by quadratic reciprocity
    bad+=(chi3(m)!=k)
for _ in range(2000):
    a,b=random.randint(1,10**6),random.randint(1,10**6); bad+=(chi3(a*b)!=chi3(a)*chi3(b))
allok&=report("1.1 chi_{-3}(m) as defined = Kronecker (-3/m); completely multiplicative",bad,0)

# ---------- (1.2) c(m) = sum_{de=m} mu(d) mu(e) chi_{-3}(e), by direct convolution from the definition
mu_int=np.array([0]+[int(smob(m)) for m in range(1,M+1)],dtype=np.int64)
ch=np.array([0]+[chi3(m) for m in range(1,M+1)],dtype=np.int64)
c=np.zeros(M+1,dtype=np.int64)
for d in range(1,M+1):
    if mu_int[d]==0: continue
    e=np.arange(1,M//d+1); c[d*e]+=mu_int[d]*mu_int[e]*ch[e]
# formal Dirichlet-series identity: c * (1 * chi) = delta, i.e. (sum c m^-s) * zeta(s) L(s,chi) = 1
one_chi=np.zeros(M+1,dtype=np.int64)
for d in range(1,M+1):
    e=np.arange(1,M//d+1); one_chi[d*e]+=ch[e]
conv=np.zeros(M+1,dtype=np.int64)
for d in range(1,M+1):
    if c[d]==0: continue
    e=np.arange(1,M//d+1); conv[d*e]+=c[d]*one_chi[e]
delta=np.zeros(M+1,dtype=np.int64); delta[1]=1
allok&=report("1.2a coefficientwise: c * (1*chi_{-3}) = delta, i.e. 1/(zeta L) has coefficients c",
              int(np.abs(conv[1:]-delta[1:]).max()),0,f"(all m <= {M})")
mp.mp.dps=30
def L3(s): return mp.power(3,-s)*(mp.zeta(s,mp.mpf(1)/3)-mp.zeta(s,mp.mpf(2)/3))
errs=[]
for s in [mp.mpf(3),mp.mpc(3,4),mp.mpf(4)]:
    lhs=mp.fsum(int(c[m])*mp.power(m,-s) for m in range(1,M+1) if c[m])
    rhs=1/(mp.zeta(s)*L3(s)); errs.append(abs(lhs-rhs))
allok&=report("1.2b sum_{m<=1e5} c(m) m^-s against 1/(zeta(s)L(s,chi_-3)) at s=3, 3+4i, 4",float(max(errs)),1e-9,
              "(truncation error only)")

# ---------- (1.3) local table and multiplicativity
bad=0; cnt=0
for p in primerange(2,M+1):
    if p==3: continue
    j=1
    while p**j<=M:
        pred={1:-2,2:1}.get(j,0) if p%3==1 else {1:0,2:-1}.get(j,0)
        bad+=(c[p**j]!=pred); cnt+=1; j+=1
bad+=(c[3]!=-1)+sum(c[3**j]!=0 for j in range(2,11))          # 3: factor 1-3^-s
allok&=report("1.3a local table of c(p^j) for all primes p != 3 (and p = 3: c(3)=-1, c(3^j)=0, j>=2)",bad,0,
              f"({cnt} prime powers <= {M})")
bad=0
for m in range(2,M+1):
    prod=1
    for p,e in ifac(m).items(): prod*=int(c[p**e])
    bad+=(prod!=c[m])
allok&=report("1.3b c is multiplicative: c(m) = prod c(p^e) for all m <= 1e5",bad,0)
allok&=report("1.3c examples c(7)=-2, c(49)=1, c(25)=-1 (and c != mu)",
              abs(c[7]+2)+abs(c[49]-1)+abs(c[25]+1)+int(c[49]==mu_int[49]),0)

# ---------- (1.5) removing 2 and 3; Euler factors at 2 and 3
cp=c.copy(); cp[np.gcd(np.arange(M+1),6)!=1]=0
rhs=cp.copy()
for d,sg in [(3,-1),(4,-1),(12,1)]:          # (1-3^-s)(1-4^-s) = 1 - 3^-s - 4^-s + 12^-s
    rhs[d::d]+=sg*cp[1:M//d+1][:len(rhs[d::d])]
allok&=report("1.5a coefficientwise: (sum_{(m,6)=1} c m^-s)(1-3^-s)(1-4^-s) = sum c m^-s",
              int(np.abs(rhs[1:]-c[1:]).max()),0,f"(all m <= {M})")
loc2=[int(c[2**j]) for j in range(0,17)]; loc3=[int(c[3**j]) for j in range(0,11)]
allok&=report("1.5b local factors: sum_j c(2^j)2^-js = 1-4^-s, sum_j c(3^j)3^-js = 1-3^-s",
              int(loc2!=[1,0,-1]+[0]*14)+int(loc3!=[1,-1]+[0]*9),0)
s=mp.mpc(3,2)
lhs=mp.fsum(int(cp[m])*mp.power(m,-s) for m in range(1,M+1) if cp[m])
rhs_=1/((1-mp.power(3,-s))*(1-mp.power(4,-s))*mp.zeta(s)*L3(s))
allok&=report("1.5c numerically at s=3+2i",float(abs(lhs-rhs_)),1e-9,"(truncation error only)")

# ---------- Section 2
w=cmath.exp(2j*math.pi/3)
allok&=report("2.0 omega^2+omega+1=0",abs(w*w+w+1),1e-15)
err=0
for a in range(-40,41):
    for b in range(-40,41):
        err=max(err,abs(abs(a+b*w)**2-(a*a-a*b+b*b)))
allok&=report("2.1 N(a+b omega) = |a+b omega|^2 = a^2-ab+b^2",err,1e-10,"(|a|,|b|<=40)")
units=[(a,b) for a in range(-3,4) for b in range(-3,4) if N((a,b))==1]
pred=sorted([(1,0),(-1,0),(0,1),(0,-1),pw((0,1),2),neg(pw((0,1),2))])
allok&=report("2.1b exactly six units, namely +-1, +-omega, +-omega^2",int(sorted(units)!=pred),0)
def is_primary_def(n): return divides((3,0),sub(n,(1,0)))        # n = 1 mod 3 in Z[omega]
bad=0; cnt=0
for a in range(-60,61):
    for b in range(-60,61):
        n=(a,b)
        if n==(0,0): continue
        bad+=(is_primary_def(n)!=(a%3==1 and b%3==0))
        k=sum(1 for u in UNITS if is_primary_def(mul(u,n)))
        bad+=(k!=(1 if N(n)%3 else 0)); cnt+=1      # lambda | n  iff 3 | N(n)
allok&=report("2.2 n = 1 mod 3 iff a=1, b=0 mod 3; each n coprime to 3 has exactly one primary associate",bad,0,
              f"({cnt} elements)")

# splitting of rational primes
bad=0; nsplit=ninert=0
for p in primerange(2,3000):
    els=[(a,b) for a in range(-2*int(p**.5)-2,2*int(p**.5)+3) for b in range(-2*int(p**.5)-2,2*int(p**.5)+3)
         if N((a,b))==p]
    if p==3:
        bad+=(mul((-1,0),pw(LAM,2))!=(3,0))    # 3 = -lambda^2 (lambda = i sqrt3), ramified
    elif p%3==1:
        prim=[n for n in els if is_primary_def(n)]
        bad+=(len(els)!=12)+(len(prim)!=2)
        if len(prim)==2:
            x,y=prim; bad+=(conj(x)!=y)+(mul(x,y)!=(p,0)); nsplit+=1
    else:
        bad+=(len(els)!=0)+(not is_primary_def((-p,0)))+(N((-p,0))!=p*p); ninert+=1
allok&=report("2.3a p=1 mod 3 splits into two conjugate primary primes of norm p; p=2 mod 3 has no element of norm p "
              "(so stays prime, norm p^2); 3 = -lambda^2",bad,0,f"({nsplit} split, {ninert} inert)")

# mu_K from the factorization of N(n) in Z (independent of the list of Eisenstein primes)
def muK(n):
    a,b=n; r=0
    for p,e in ifac(N(n)).items():
        if p%3==2:
            if e!=2: return 0
            r+=1
        else:                                  # p = 1 mod 3
            if e==1: r+=1
            elif e==2 and a%p==0 and b%p==0: r+=2    # n divisible by p = pi*conj(pi)
            else: return 0
    return (-1)**r
# all primary elements with N(n) <= M, bucketed by norm
B=int(2*math.sqrt(M/3))+3
cnt_by=np.zeros(M+1,dtype=np.int64); mu_by=np.zeros(M+1,dtype=np.int64)
for b in range(-B,B+1):
    if b%3: continue
    # a^2 - ab + b^2 <= M  <=>  |2a-b| <= sqrt(4M-3b^2)
    r=4*M-3*b*b
    if r<0: continue
    s=math.isqrt(r)
    for a in range((b-s)//2-1,(b+s)//2+2):
        if a%3!=1: continue
        m=a*a-a*b+b*b
        if 0<m<=M:
            cnt_by[m]+=1
            if m%2 and m%3: mu_by[m]+=muK((a,b))
# cross-check muK against the Eisenstein-prime factorization on a sample
bad=0
for n in random.sample(primaries(20000),400): bad+=(muK(n)!=mu(n))
allok&=report("2.3b mu_K from N(n) agrees with mu_K from Eisenstein factorization",bad,0,"(400 random n)")
# count of primary elements of norm m (unique factorization + splitting) = sum_{d|m} chi_{-3}(d), (m,3)=1
r1=np.zeros(M+1,dtype=np.int64)
for d in range(1,M+1):
    e=np.arange(1,M//d+1); r1[d*e]+=ch[d]
ms=np.array([m for m in range(1,M+1) if m%3])
allok&=report("2.3c #{n primary: N(n)=m} = sum_{d|m} chi_{-3}(d) for all m <= 1e5, 3 nmid m",
              int(np.abs(cnt_by[ms]-r1[ms]).max()),0)
ms6=np.array([m for m in range(1,M+1) if m%2 and m%3])
allok&=report("2.3 (2.3) sum_{N(n)=m} mu_K(n) = c(m) for all (m,6)=1, m <= 1e5",
              int(np.abs(mu_by[ms6]-c[ms6]).max()),0,f"({len(ms6)} values of m)")
# local factors (1-p^-s)^2 and 1-p^-2s
bad=0
for p in primerange(5,400):
    j=0; ser=[]
    while p**j<=M: ser.append(int(mu_by[p**j]) if j else 1); j+=1
    pred=([1,-2,1]+[0]*20)[:len(ser)] if p%3==1 else ([1,0,-1]+[0]*20)[:len(ser)]
    bad+=(ser!=pred)
allok&=report("2.3d local factors: (1-p^-s)^2 for split p, 1-p^-2s for inert p (p < 400)",bad,0)

# (2.4) A_1(D) via integers (1.4) equals the Eisenstein sum (2.4), complex smooth weight
def Wc(y):
    if not 0.5<y<2.3: return 0
    return math.exp(-1/((y-0.5)*(2.3-y)))*cmath.exp(0.7j*y)
err=0
for D in [300.0,3000.0,40000.0]:
    A_int=sum(int(c[m])*Wc(m/D) for m in range(1,M+1) if m%2 and m%3 and c[m])
    A_eis=0
    for b in range(-B,B+1):
        if b%3: continue
        r=4*M-3*b*b
        if r<0: continue
        s=math.isqrt(r)
        for a in range((b-s)//2-1,(b+s)//2+2):
            if a%3!=1: continue
            m=a*a-a*b+b*b
            if 0<m<=M and m%2:
                w_=Wc(m/D)
                if w_: A_eis+=muK((a,b))*w_
    err=max(err,abs(A_int-A_eis)/max(1e-300,abs(A_eis)))
    print(f"      D={D:g}: A_1 = {complex(A_eis):.10g}")
allok&=report("2.4 (2.4) A_1(D) via integers (1.4) = sum over Eisenstein integers",err,1e-12,"(D=300, 3000, 40000)")
print("ALL PASS" if allok else "SOME FAILURES")
