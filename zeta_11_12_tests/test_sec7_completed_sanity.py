"""Sanity check (not a proof) of the completed estimate, guide (7.8):

    sum_{0 < N(k) <= H} |T(X; k, f)|^2  <<  D^eps (H + H^2 N(f)/X).

T(X; k, f) is the completed sum (7.2), computed directly from the cubic Gauss coefficients and the cube
indices b (X = 1500).  For f = 1 and f = 7 + 9 omega and H = 30, 100, 300, 1000 the script prints the
partial sums and their ratio to H + H^2 N(f)/X.  The estimate is a bound, so there is no pass/fail
criterion; the ratios should stay bounded (they are below 0.3).

Run time: about 12 seconds.
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from fastgauss import g_prime, pairsym
gp={}
def g_cub(n):
    f=factor(n); g=1+0j
    for p in f:
        if p not in gp: gp[p]=g_prime(p)
        g*=gp[p]
    for i in range(len(f)):
        for j in range(i+1,len(f)): g*=pairsym(f[i],f[j])*pairsym(f[j],f[i])
    return g
def a(n): return alpha(n).conjugate()*chi(n,LAM)**2*g_cub(n)/math.sqrt(N(n)) if N(n)>1 else 1.0
def Wb(y): return math.exp(4-1/((y-1)*(2-y))) if 1<y<2 else 0.0   # bump normalized to max 1
X=1500.0
terms=[]
for b in primaries(int((2*X)**(1/3))+1):
    Nb3=N(b)**3
    for n in primaries(int(2*X/Nb3),sqfree=True):
        w=Wb(N(n)*Nb3/X)
        if w: terms.append((n,b,a(n)*(alpha(b).conjugate()**3 if b!=(1,0) else 1)*math.sqrt(N(b))*w))
print("terms in T:",len(terms))
def T(k,f):
    s=0
    for n,b,c in terms:
        Psi=lambda m: chi(m,k)*chi(m,f)**4
        s+=c*Psi(n)*(Psi(b)**3 if b!=(1,0) else 1)
    return s/math.sqrt(X)
R=40; rows=sorted([(x,y) for x in range(-R,R+1) for y in range(-R,R+1) if 0<N((x,y))<=1000],key=N)
for f in [(1,0),(7,9)]:
    vals=[]; tot=0; i=0
    for H in [30,100,300,1000]:
        while i<len(rows) and N(rows[i])<=H:
            tot+=abs(T(rows[i],f))**2; i+=1
        bound=H+H*H*N(f)/X
        vals.append((H,tot,tot/bound))
    print(f"[guide (7.8), sanity] X={X:g}, f={f}: "+"; ".join(f"H={H}: sum={s:.2f}, ratio to H+H^2N(f)/X = {r:.3f}" for H,s,r in vals))
