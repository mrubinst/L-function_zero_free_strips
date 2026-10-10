"""Section 3 of the guide: the sextic symbol, the coordinate formula and the double sum.

Tests (guide displays in brackets):
  3.0a  [text before (3.1)]  6 | N(pi) - 1 for the first 3000 primary primes pi coprime to 6.
  3.0b,c [after (3.1)]       chi_n is multiplicative in the numerator and in the denominator.
  3.1   [(3.4)]              "primary and coprime to 6" <=> a = 1, b = 0 mod 3, (a, b) not both even.
  3.2a  [(3.4)]              the integer-coordinate formula equals M(D, H) of (3.2)-(3.3) (D = 60, H = 90,
                             complex weight W).
  3.2b  [(3.5), (3.6)]       the double sum over columns with S_H computed directly equals M(D, H).
  3.2c  [(3.7)]              A_u(D) = sum_{(m,6)=1} c_u(m) W(m/D) for 40 rows u.
  3.2d  [after (3.7)]        printed example: c_u(p) depends on the residue symbols, not only on mu, chi_{-3}.
  3.3   [(3.8), sanity]      printed: M(D, D^1.1)/(D H) for D = 100, 200, 400.  (3.8) is a bound, not an
                             identity; the ratio should stay bounded (it is about 0.3).

Run time: a few seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import random
random.seed(17); allok=True

# ---- 3.0 facts in the definition of the sextic symbol
bad=sum(1 for p in PR_ODD[:3000] if (N(p)-1)%6!=0)
allok&=report("3.0a 6 | N(pi)-1 for primary primes pi coprime to 6 (first 3000)",bad,0)
ns=primaries(300); err=0
for _ in range(400):
    n=random.choice(ns); u=(random.randint(-30,30),random.randint(-30,30)); v=(random.randint(-30,30),random.randint(-30,30))
    err=max(err,abs(chi(n,mul(u,v))-chi(n,u)*chi(n,v)))
allok&=report("3.0b chi_n is multiplicative in the numerator",err,1e-12)
err=0
for _ in range(300):
    n1,n2=random.choice(ns),random.choice(ns); u=(random.randint(-30,30),random.randint(-30,30))
    err=max(err,abs(chi(mul(n1,n2),u)-chi(n1,u)*chi(n2,u)))
allok&=report("3.0c chi_n is multiplicative in the denominator",err,1e-12)

# ---- 3.1 coordinate description: primary and coprime to 6  <=>  a=1 (3), b=0 (3), a,b not both even
bad=0
for a in range(-40,41):
    for b in range(-40,41):
        n=(a,b)
        if N(n)==0: continue
        lhs=primary(n) and not divides((2,0),n) and not divides(LAM,n)
        rhs=(a%3==1 and b%3==0 and not (a%2==0 and b%2==0))
        bad+=(lhs!=rhs)
allok&=report("3.1 primary & coprime to 6 <=> a=1 mod 3, b=0 mod 3, (a,b) not both even",bad,0)

# ---- 3.2 (3.4) = definition of M; (3.5) with S_H of (3.6); grouping by norms
def Wc(y): return complex(math.exp(4-1/((y-1)*(2-y))),0.3*math.sin(5*y)) if 1<y<2 else 0   # complex weight
D=60; H=90
cols=[n for n in primaries(2*D) if N(n)>D]              # all primary n coprime to 6 (mu kills non-squarefree)
rows=[(x,y) for x in range(-12,13) for y in range(-12,13) if 0<N((x,y))<=H]
A={u:sum(mu(n)*chi(n,u)*Wc(N(n)/D) for n in cols) for u in rows}
M_def=sum(abs(v)**2 for v in A.values())
M_coord=0
for (x,y) in rows:
    s=0
    for a in range(-20,21):
        for b in range(-20,21):
            if a%3==1 and b%3==0 and not (a%2==0 and b%2==0) and D<a*a-a*b+b*b<2*D:
                s+=mu((a,b))*chi((a,b),(x,y))*Wc((a*a-a*b+b*b)/D)
    M_coord+=abs(s)**2
allok&=report("3.2a integer-coordinate sum = M(D,H)",abs(M_def-M_coord)/M_def,1e-12,f"(M={M_def:.4f})")
def S_H(n1,n2): return sum(chi(n1,u).conjugate()*chi(n2,u) for u in rows)
M_double=sum(mu(n1)*mu(n2)*Wc(N(n1)/D).conjugate()*Wc(N(n2)/D)*S_H(n1,n2) for n1 in cols for n2 in cols if mu(n1) and mu(n2))
allok&=report("3.2b (3.5) with S_H of (3.6)",abs(M_def-M_double)/M_def,1e-12)
err=0
for u in rows[:40]:
    byN={}
    for n in cols: byN.setdefault(N(n),0); byN[N(n)]+=mu(n)*chi(n,u)
    grouped=sum(cm*Wc(m/D) for m,cm in byN.items() if math.gcd(m,6)==1)
    err=max(err,abs(grouped-A[u]))
allok&=report("3.2c A_u(D)=sum_{(m,6)=1} c_u(m)W(m/D)",err,1e-12)
# c_u(m) is not determined by mu and chi_{-3}: two rows with different c_u(p) at the same split p
p=(7,9); P=N(p); u1,u2=(2,0),(5,1)
c1=sum(mu(n)*chi(n,u1) for n in [p,conj(p)]); c2=sum(mu(n)*chi(n,u2) for n in [p,conj(p)])
print(f"      3.2d [guide text after (3.7)] c_u({P}) for u=2: {complex(c1):.3f};  for u=5+omega: {complex(c2):.3f}  (depend on the residue symbols)")

# ---- 3.3 sanity check of the size in (3.8): M(D,D^{1+theta}) / (D^{1+theta} D)
for D in [100,200,400]:
    H=int(D**1.1)
    cols=[n for n in primaries(2*D,sqfree=True) if N(n)>D]
    R=int(2*math.sqrt(H))+2
    rws=[(x,y) for x in range(-R,R+1) for y in range(-R,R+1) if 0<N((x,y))<=H]
    M=sum(abs(sum(mu(n)*chi(n,u)*Wc(N(n)/D) for n in cols))**2 for u in rws)
    print(f"      3.3 [guide (3.8), sanity] D={D}, H=D^1.1={H}: M/(D H) = {M/(D*H):.3f}  ((3.8) predicts at most D^eps)")
print("ALL PASS" if allok else "SOME FAILURES")
