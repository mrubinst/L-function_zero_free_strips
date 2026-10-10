"""Arithmetic in Z[omega], omega = e^{2 pi i/3}; elements are pairs (a,b) = a + b*omega.
Generates the LaTeX tables for the notes and checks every numerical claim made in them.

Usage: python3 zomega.py <output directory>
Writes rational_primes.tex, primary_elements.tex, cubic_table.tex and the four sextic_*.tex tables."""
import cmath, math, itertools, sys, os
from fractions import Fraction
W = cmath.exp(2j*math.pi/3)
def mul(x,y): a,b=x; c,d=y; return (a*c-b*d, a*d+b*c-b*d)
def conj(x): a,b=x; return (a-b,-b)
def N(x): a,b=x; return a*a-a*b+b*b
def tr(x): a,b=x; return 2*a-b
def neg(x): return (-x[0],-x[1])
def cx(x): return x[0]+x[1]*W
def divides(d,x):
    n=N(d); q=mul(x,conj(d)); return q[0]%n==0 and q[1]%n==0
def div(x,d):
    n=N(d); q=mul(x,conj(d)); assert q[0]%n==0 and q[1]%n==0; return (q[0]//n,q[1]//n)
UNITS=[(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]   # 1,-1,w,-w,1+w=-w^2,-1-w=w^2
def primary(x): return x[0]%3==1 and x[1]%3==0
def primary_assoc(x):
    c=[mul(u,x) for u in UNITS if primary(mul(u,x))]; assert len(c)==1,x; return c[0]
def isprime_int(n): return n>1 and all(n%k for k in range(2,int(n**.5)+1))
def fmt(x):
    a,b=x
    if b==0: return f"{a}"
    bs = "" if abs(b)==1 else str(abs(b))
    if a==0: return ("-" if b<0 else "")+bs+r"\omega"
    return f"{a}{'+' if b>0 else '-'}{bs}\\omega"
# ---------- checks: basic identities
assert sorted(UNITS)==sorted([x for x in itertools.product(range(-2,3),repeat=2) if N(x)==1])
lam=(1,-1)
assert mul(lam,lam)==(0,-3)                      # lambda^2 = -3 omega
w2=mul((0,1),(0,1)); assert w2==(-1,-1)
assert mul(neg(w2),mul(lam,lam))==(3,0)          # 3 = -omega^2 lambda^2
assert conj(lam)==mul(neg(w2),lam)               # conj(lambda) = -omega^2 lambda
# primary elements up to norm X, with factorization into primary primes
X=100
prim=[(a,b) for a in range(-12,13) for b in range(-12,13) if primary((a,b)) and 0<N((a,b))<=X]
pp=sorted([x for x in [(a,b) for a in range(-20,21) for b in range(-20,21)] if primary(x) and N(x)>1 and N(x)<=400 and
          ((isprime_int(N(x)) and N(x)%3==1) or (N(x)==int(round(N(x)**.5))**2 and isprime_int(int(round(N(x)**.5))) and int(round(N(x)**.5))%3==2))],
          key=lambda x:(N(x),-x[1],x[0]))
def factor(x):
    out=[]
    for p in pp:
        while N(x)>1 and divides(p,x): out.append(p); x=div(x,p)
    assert x==(1,0), ("unit left",x)   # primary => no unit factor
    return out
def sortkey(x): return (N(x),-x[1],x[0])
prim.sort(key=sortkey)
# ---------- primary-element table
rows=[]
for x in prim:
    f=factor(x)
    if x==(1,0): fs=r"\text{unit}"
    elif len(f)==1: fs=r"\text{prime}"
    else:
        from collections import Counter
        c=Counter(f); fs=r"\cdot".join(("("+fmt(p)+")" if p[1]!=0 or p[0]<0 else fmt(p))+(f"^{c[p]}" if c[p]>1 else "") for p in sorted(c,key=sortkey))
    rows.append((fmt(x),N(x),fs))
# check: count of primary elements of norm n = number of ideals of norm n prime to 3 = sum_{d|n} chi(d) for 3 not dividing n
chi=lambda n: 0 if n%3==0 else (1 if n%3==1 else -1)
for n in range(1,X+1):
    cnt=sum(1 for x in prim if N(x)==n)
    exp=0 if n%3==0 else sum(chi(d) for d in range(1,n+1) if n%d==0)
    assert cnt==exp,(n,cnt,exp)
    # all elements of norm n: 6 * sum chi(d)
    tot=sum(1 for a in range(-15,16) for b in range(-15,16) if N((a,b))==n)
    assert tot==6*sum(chi(d) for d in range(1,n+1) if n%d==0),(n,tot)
print("primary elements of norm <=",X,":",len(prim))
half=(len(rows)+1)//2
with open(os.path.join(sys.argv[1],"primary_elements.tex"),"w") as f:
    for i in range(half):
        L=rows[i]; R=rows[i+half] if i+half<len(rows) else None
        f.write(f"${L[0]}$ & {L[1]} & ${L[2]}$ & " + (f"${R[0]}$ & {R[1]} & ${R[2]}$" if R else "& &") + r"\\"+"\n")
# ---------- rational primes table
def two_sq(p):   # p = x^2+3y^2, x,y>=0
    for y in range(0,int(p**.5)+1):
        r=p-3*y*y
        if r>=0 and int(round(r**.5))**2==r: return (int(round(r**.5)),y)
prow=[]
for p in [q for q in range(2,100) if isprime_int(q)]:
    if p==3:
        prow.append((p,0,"ram.","","","","")); continue
    if p%3==2:
        prow.append((p,2,"inert",f"{p}=-1\\cdot({-p})","--","--","--")); continue
    pis=sorted([x for x in pp if N(x)==p],key=lambda x:(-x[1],x[0])); assert len(pis)==2 and pis[1]==conj(pis[0])
    assert mul(pis[0],pis[1])==(p,0)               # p = pi * conj(pi), no unit
    a,b=pis[0]; L=tr(pis[0]); M=b//3
    assert 4*p==L*L+27*M*M and L%3==2
    A=-L                                            # Gauss's A = 1 mod 3
    # Gauss: #{x^3+y^3=1 in F_p} = p-2+A
    cubes={}
    for t in range(p): cubes[t**3%p]=cubes.get(t**3%p,0)+1
    cnt=sum(cubes.get(u,0)*cubes.get((1-u)%p,0) for u in range(p))
    assert cnt==p-2+A,(p,cnt,A)
    # Gauss: 2 is a cube mod p iff L, M even iff p = x^2+27y^2
    two_cube=pow(2,(p-1)//3,p)==1
    x27=any(int(round((p-27*y*y)**.5))**2==p-27*y*y for y in range(0,int((p/27)**.5)+1))
    assert two_cube==(L%2==0)==x27,(p)
    x,y=two_sq(p)
    # a^2-ab+b^2 representation with a,b>=0
    rep=[(c,d) for c in range(0,p) for d in range(0,c+1) if c*c-c*d+d*d==p][0]
    prow.append((p,1,"split",f"{p}=({fmt(pis[0])})({fmt(pis[1])})",f"{rep[0]}^2-{rep[0]}\\cdot{rep[1]}+{rep[1]}^2",f"{x}^2+3\\cdot{y}^2",(f"({A})^2" if A<0 else f"{A}^2")+f"+27\\cdot{abs(M)}^2"))
with open(os.path.join(sys.argv[1],"rational_primes.tex"),"w") as f:
    for r in prow:
        p,m,t,fa,r1,r2,r3=r
        if p==3: f.write(r"$3$ & $0$ & ramified & $3=-\omega^2\lambda^2$ & $1^2-1\cdot2+2^2$ & $0^2+3\cdot1^2$ & --\\"+"\n"); continue
        ch = "1" if m==1 else "-1"
        f.write(f"${p}$ & ${ch}$ & {t} & ${fa}$ & " + (f"${r1}$ & ${r2}$ & ${r3}$" if m==1 else "-- & -- & --") + r"\\"+"\n")
assert 1-2+4==3
# ---------- residue symbols
def powmod_pi(x,e,pi):
    """x^e reduced mod the prime pi; returns representative pair. Works for any prime pi."""
    n=N(pi)
    def red(y):
        # reduce y modulo pi: subtract nearest multiple
        q=mul(y,conj(pi)); qa=round(q[0]/n); qb=round(q[1]/n)
        return (y[0]-mul((qa,qb),pi)[0], y[1]-mul((qa,qb),pi)[1])
    r=(1,0); b=red(x)
    while e:
        if e&1: r=red(mul(r,b))
        b=red(mul(b,b)); e>>=1
    return r
SIX=[(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1)]   # zeta^j, zeta = 1+omega = e^{i pi/3}: 1, zeta, zeta^2=omega, -1, zeta^4=omega^2, zeta^5
assert all(abs(cx(SIX[j])-cmath.exp(1j*math.pi*j/3))<1e-12 for j in range(6))
def sym(x,pi,k):
    """(x/pi)_k for k in {2,3,6} as exponent j: value = zeta^(j*6/k)... return the root of unity as pair"""
    n=N(pi); assert (n-1)%k==0
    if divides(pi,x): return None
    v=powmod_pi(x,(n-1)//k,pi)
    roots={6:SIX,3:[SIX[0],SIX[2],SIX[4]],2:[SIX[0],SIX[3]]}[k]
    for u in roots:
        if divides(pi,(v[0]-u[0],v[1]-u[1])): return u
    raise Exception
def cube_exp(x,pi):
    u=sym(x,pi,3); return {(1,0):0,(0,1):1,(-1,-1):2}[u]
# cubic reciprocity for primary primes (both split and inert), norms <= 400
test=[p for p in pp if N(p)>4]   # exclude 2 (fine for cubic, but keep reciprocity test to odd)
for p1,p2 in itertools.combinations(pp,2):
    if N(p1)==N(p2) and p1!=conj(p2) and False: pass
    if math.gcd(N(p1),N(p2))==1 or (N(p1)==N(p2) and p1!=p2):
        if N(p1)==N(p2): continue
        assert cube_exp(p1,p2)==cube_exp(p2,p1),(p1,p2)
# supplements in our convention: pi primary = (3m+1) + 3n omega
for pi in pp:
    if pi==(-2,0): pass
    a,b=pi; m=(a-1)//3; n=b//3
    assert cube_exp((0,1),pi)==(-(m+n))%3                        # (omega/pi)_3
    assert cube_exp((0,1),pi)==((N(pi)-1)//3)%3
    assert cube_exp(lam,pi)==m%3,(pi,cube_exp(lam,pi),m)          # (lambda/pi)_3 = omega^m
    assert cube_exp((3,0),pi)==n%3,(pi,)                          # (3/pi)_3 = omega^n
    if N(pi)%6==1:
        assert sym((-1,0),pi,6)==((1,0) if ((N(pi)-1)//6)%2==0 else (-1,0))
        # rational integers coprime to p: (k/pi)_6 * (k/conj pi)_6 = 1 ; inert: (k/q)_6=1
        for k in range(2,12):
            if math.gcd(k,N(pi))==1:
                s1=sym((k,0),pi,6); s2=sym((k,0),conj(pi),6)
                assert mul(s1,s2)==(1,0)
                if pi[1]==0: assert s1==(1,0)
print("cubic reciprocity and supplements verified for",len(pp),"primary primes")
# Gauss and Jacobi sums: our convention pi = 1 mod 3, N(pi)=p
for pi in [x for x in pp if isprime_int(N(x))]:
    p=N(pi); a,b=pi
    r=(-a*pow(b,-1,p))%p                     # omega = r mod pi
    assert (r*r+r+1)%p==0 and divides(pi,( -r,1))
    def chi3(t):
        t%=p
        if t==0: return 0
        v=pow(t,(p-1)//3,p)
        return {1:1, r:W, (r*r)%p:W*W}[v]
    g=sum(chi3(t)*cmath.exp(2j*math.pi*t/p) for t in range(1,p))
    assert abs(abs(g)**2-p)<1e-6
    assert abs(g**3-(-p*cx(pi)))<1e-6*p**1.5, (pi,g**3/(p*cx(pi)))      # g^3 = -p*pi
    J=sum(chi3(t)*chi3(1-t) for t in range(p))
    assert abs(J-(-cx(pi)))<1e-9                                         # J(chi,chi) = -pi
    # chi3(t) agrees with cube_exp on rational t
    for t in range(1,min(p,30)): assert abs(chi3(t)-[1,W,W*W][cube_exp((t,0),pi)])<1e-12
print("Gauss/Jacobi sums: g^3 = -p*pi and J = -pi for primary pi (=1 mod 3)")
# cubic residue table among small primary primes
small=[x for x in pp if N(x)<=31]
with open(os.path.join(sys.argv[1],"cubic_table.tex"),"w") as f:
    f.write("&"+"&".join(f"${fmt(q)}$" for q in small)+r"\\ \midrule"+"\n")
    for p1 in small:
        cells=[]
        for p2 in small:
            if p1==p2: cells.append(r"\cdot")
            else:
                e=cube_exp(p1,p2); cells.append(["1",r"\omega",r"\omega^2"][e])
        f.write(f"${fmt(p1)}$&"+"&".join(f"${c}$" for c in cells)+r"\\"+"\n")
# Euclid example: gcd(61, 13-omega)
x,y=(61,0),(13,-1); steps=[]
while N(y)>0:
    q=mul(x,conj(y)); n=N(y); qq=(round(q[0]/n),round(q[1]/n))
    rr=(x[0]-mul(qq,y)[0],x[1]-mul(qq,y)[1]); steps.append((x,qq,y,rr,N(rr))); x,y=y,rr
print("Euclid steps gcd(61,13-w):",[(fmt(s[0]),fmt(s[1]),fmt(s[2]),fmt(s[3]),s[4]) for s in steps])
print("primary associate:",fmt(primary_assoc(steps[-2][3] if N(steps[-1][3])==0 else steps[-1][3])))
# L(1,chi_-3) and class number formula
L1=sum(chi(n)/n for n in range(1,2000001))
print("L(1,chi_-3) ~",L1," pi/(3 sqrt3) =",math.pi/(3*math.sqrt(3)))
h=-(6/(2*3))*sum(a*chi(a) for a in range(1,3)); print("h from finite formula:",h)
print("small primary primes:",[ (fmt(x),N(x)) for x in pp if N(x)<=100])

# conjugate pairs: (pi/conj pi)_3 vs (conj pi/pi)_3
for x in [y for y in pp if isprime_int(N(y))]:
    e1=cube_exp(x,conj(x)); e2=cube_exp(conj(x),x)
    print("pair",fmt(x),N(x),"(pi/pibar)=w^",e1,"(pibar/pi)=w^",e2)

# ---------- sextic residue symbol: checks and tables
SIXFMT={(1,0):"1",(1,1):"-\\omega^2",(0,1):"\\omega",(-1,0):"-1",(-1,-1):"\\omega^2",(0,-1):"-\\omega"}
def s6(x,pi):
    """(x/pi)_6 as a pair (a unit), or None if pi | x."""
    return sym(x,pi,6)
def f6(x,pi):
    u=s6(x,pi); return r"\cdot" if u is None else SIXFMT[u]
def unit_pow(u,e):
    r=(1,0)
    for _ in range(e%6): r=mul(r,u)
    return r
def uorder(u):
    for k in range(1,7):
        if unit_pow(u,k)==(1,0): return k
odd_pp=[p for p in pp if N(p)%2==1]                  # primary primes prime to 6, norm <= 400
# units: (u/pi)_6 = u^{(N pi-1)/6} exactly, so it depends only on N pi mod 36
for pi in odd_pp:
    for u in UNITS:
        assert s6(u,pi)==unit_pow(u,(N(pi)-1)//6)
    assert s6((0,1),pi)==unit_pow((0,1),(N(pi)-1)//6)              # (omega/pi)_6 = omega^{(N pi-1)/6}
    assert s6((-1,0),pi)==((1,0) if N(pi)%12==1 else (-1,0))        # -1 iff N pi = 7 mod 12
for p1 in odd_pp:
    for p2 in odd_pp:
        if N(p1)%36==N(p2)%36:
            assert all(s6(u,p1)==s6(u,p2) for u in UNITS)
# rational integers k, p the rational prime below pi
for pi in odd_pp:
    n=N(pi)
    if pi[1]==0:                                      # inert: pi = -q, N pi = q^2
        q=-pi[0]; assert q%6==5 and n==q*q
        for k in range(1,60):
            if k%q: assert s6((k,0),pi)==(1,0)
            # Fermat route: k^{(q^2-1)/6} = (k^{q-1})^{(q+1)/6} = 1 mod q
            if k%q: assert pow(k,(n-1)//6,q)==1
    else:                                             # split: N pi = p
        p=n
        for k in range(1,60):
            if k%p==0: continue
            u=s6((k,0),pi)
            assert mul(u,s6((k,0),conj(pi)))==(1,0)                   # (k/pibar)_6 = conj (k/pi)_6
            leg=1 if pow(k,(p-1)//2,p)==1 else -1
            assert unit_pow(u,3)==(leg,0)                              # (k/pi)_2 = Legendre (k/p)
            is_cube=any(pow(t,3,p)==k%p for t in range(1,p))
            assert (unit_pow(u,2)==(1,0))==is_cube                      # (k/pi)_3 = 1 iff k is a cube mod p
            # O/pi = F_p via omega -> r; the symbol is the unit congruent to k^{(p-1)/6} mod p
            a,b=pi; r=(-a*pow(b,-1,p))%p
            v=pow(k,(p-1)//6,p); assert (u[0]+u[1]*r-v)%p==0
# (k/m)_6 = 1 for coprime rational integers k, m prime to 6 (m composite included)
for m in [5,7,11,13,25,35,49,55,65,77,91,175,245,1001]:
    fs=[]
    for q in [5,7,11,13]:
        e=0
        while m%(q**(e+1))==0: e+=1
        for _ in range(e): fs+= [x for x in odd_pp if (N(x)==q and q%3==1) or (x==(-q,0))]
    for k in range(1,40):
        if math.gcd(k,m)==1:
            v=(1,0)
            for x in fs: v=mul(v,s6((k,0),x))
            assert v==(1,0),(k,m)
# Table: the sextic character modulo the primes of norm 7 and 13 (rational representatives x)
def ordmod(x,p):
    k=1; y=x%p
    while y!=1: y=y*x%p; k+=1
    return k
with open(os.path.join(sys.argv[1],"sextic_small.tex"),"w") as f:
    for p in (7,13):
        pi=[x for x in odd_pp if N(x)==p][0]; pib=conj(pi)
        assert primary(pib) and pib in odd_pp
        vals=[s6((x,0),pi) for x in range(1,p)]
        assert sorted(set(vals))==sorted(UNITS)                        # onto the sixth roots of unity
        assert all(vals.count(u)==(p-1)//6 for u in UNITS)
        for x in range(1,p):                                           # kernel = sixth powers mod p
            assert (vals[x-1]==(1,0))==any(pow(t,6,p)==x for t in range(1,p))
        if p==7:
            assert all(uorder(vals[x-1])==ordmod(x,7) for x in range(1,7))  # an isomorphism F_7^* -> mu_6
        rows=[("x",[str(x) for x in range(1,p)]),
              (f"\\text{{order of }}x\\bmod{p}",[str(ordmod(x,p)) for x in range(1,p)]),
              (f"(x/\\pi)_6",[f6((x,0),pi) for x in range(1,p)]),
              (f"(x/\\bar\\pi)_6",[f6((x,0),pib) for x in range(1,p)]),
              (f"(x/\\pi)_3",[SIXFMT[unit_pow(s6((x,0),pi),2)] for x in range(1,p)]),
              (f"(x/\\pi)_2",[SIXFMT[unit_pow(s6((x,0),pi),3)] for x in range(1,p)]),
              (f"\\bigl(\\frac x{{{p}}}\\bigr)",[("1" if pow(x,(p-1)//2,p)==1 else "-1") for x in range(1,p)])]
        f.write(f"\\multicolumn{{13}}{{@{{}}l}}{{$p={p}$, \\ $\\pi={fmt(pi)}$, \\ $\\bar\\pi={fmt(pib)}$}}\\\\\\midrule\n")
        for name,cells in rows:
            f.write(f"${name}$&"+"&".join(f"${c}$" for c in cells)+"&"*(12-len(cells))+r"\\"+"\n")
        if p==7: f.write(r"\midrule"+"\n")
# Table: ((a+b omega)/5)_6 on O/5 = F_25
five=(-5,0)
cnt={}
for a in range(5):
    for b in range(5):
        if (a,b)==(0,0): continue
        u=s6((a,b),five); cnt[u]=cnt.get(u,0)+1
        assert s6(conj((a,b)),five)==conj(u)                           # conj(alpha) -> conj symbol (5 = conj 5)
        assert s6(mul((0,1),(a,b)),five)==mul((0,1),u)                 # times omega -> times (omega/5)_6 = omega
        if b==0: assert u==(1,0)
assert all(cnt[u]==4 for u in UNITS)
with open(os.path.join(sys.argv[1],"sextic_five.tex"),"w") as f:
    f.write("$b\\backslash a$&"+"&".join(f"${a}$" for a in range(5))+r"\\ \midrule"+"\n")
    for b in range(5):
        f.write(f"${b}$&"+"&".join(f"${f6((a,b),five)}$" for a in range(5))+r"\\"+"\n")
# Table: units, lambda and small rational integers modulo the primary primes of norm <= 121 prime to 6
cols=[(-1,0),(0,1),(1,1),lam,(2,0),(3,0),(5,0),(7,0)]
with open(os.path.join(sys.argv[1],"sextic_units.tex"),"w") as f:
    for pi in [x for x in odd_pp if N(x)<=121]:
        f.write(f"${fmt(pi)}$&${N(pi)}$&${N(pi)%36}$&"+"&".join(f"${f6(x,pi)}$" for x in cols)+r"\\"+"\n")
# 2 is a cube mod p (Gauss) iff p = x^2+27y^2: the 2-column entry is +-1
for pi in [x for x in odd_pp if isprime_int(N(x))]:
    p=N(pi); x27=any(int(round((p-27*y*y)**.5))**2==p-27*y*y for y in range(0,int((p/27)**.5)+1))
    assert (unit_pow(s6((2,0),pi),2)==(1,0))==x27
# 3 = -omega^2 lambda^2: (3/pi)_6 = (-omega^2/pi)_6 (lambda/pi)_6^2
for pi in odd_pp: assert s6((3,0),pi)==mul(s6((1,1),pi),unit_pow(s6(lam,pi),2))
# Table: (pi/rho)_6 for the primary primes of norm <= 31 prime to 6
smallo=[x for x in odd_pp if N(x)<=31]
with open(os.path.join(sys.argv[1],"sextic_recip.tex"),"w") as f:
    f.write("&"+"&".join(f"${fmt(q)}$" for q in smallo)+r"\\ \midrule"+"\n")
    for p1 in smallo:
        f.write(f"${fmt(p1)}$&"+"&".join("$\\cdot$" if p1==p2 else f"${f6(p1,p2)}$" for p2 in smallo)+r"\\"+"\n")
# reciprocity: (pi/rho)_6 = +-(rho/pi)_6; the sign (pi/rho)_2 (rho/pi)_2 depends only on pi, rho mod 4,
# and it is 1 if pi or rho is = 1, omega or omega^2 mod 4.  Checked for all pairs of norm < 1000.
big=sorted([x for x in [(a,b) for a in range(-40,41) for b in range(-40,41)] if primary(x) and N(x)%2==1 and 1<N(x)<1000 and
            ((isprime_int(N(x)) and N(x)%3==1) or (x[1]==0 and isprime_int(-x[0]) and (-x[0])%3==2))],key=sortkey)
signs={}; cube_roots_mod4=[(1,0),(0,1),(3,3)]
for p1,p2 in itertools.combinations(big,2):
    u,v=s6(p1,p2),s6(p2,p1)
    assert unit_pow(u,2)==unit_pow(v,2)                                 # cubic reciprocity
    sg=mul(unit_pow(u,3),unit_pow(v,3))[0]
    c1=(p1[0]%4,p1[1]%4); c2=(p2[0]%4,p2[1]%4); key=tuple(sorted([c1,c2]))
    assert signs.setdefault(key,sg)==sg                                 # depends only on classes mod 4
    if c1 in cube_roots_mod4 or c2 in cube_roots_mod4: assert sg==1
print("sextic symbol: tables written;",len(big),"primary primes of norm < 1000 used for the reciprocity sign;",
      sum(1 for s in signs.values() if s==-1),"of",len(signs),"class pairs mod 4 have sign -1")
