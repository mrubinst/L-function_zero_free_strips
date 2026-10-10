"""End-to-end check of the theta transformation of guide Section 7 ("The transformation and its
purpose") at a single cusp, with one zbar-derivative (the case needed for zeta).

For a cusp a/C (a, C in Z[omega] coprime) the script compares, for a smooth weight V,

  LHS(X) = sum_{0 != l in lambda^{-3} O} d_0(l) ce(l a/C) conj(alpha(l)) N(l)^{-1/2} V(N(l)/X)

with the reflected sum

  RHS(X) = - conj(kappa(g1)) conj(alpha(C))^2
           * sum_{0 != l in lambda^{-4} O} d_sigma(l) alpha(l) ce(-delta' l/C) N(l)^{-1/2} V#((2 pi)^4 N(l) X/N(C)^2),

  V#(y) = (1/2 pi i) int_{Re t = 0} V^(-t) Q(t) y^{-t} dt,    Q(t) = Gamma(7/6+t) Gamma(5/6+t) / (Gamma(7/6-t) Gamma(5/6-t)),

which is the kernel of OpenAI's 11/12 paper (its equation (theta-weight)) and the kernel quoted in the
guide's table in Section 9.  Here
  * ce(z) = exp(2 pi i (z + conj z)), alpha(l) = l/|l|;
  * d_sigma, sigma in {0, +, -}, are the Fourier coefficients of conj(theta) at the three cusps, from
    Patterson's formulas as given by Dunn-Radziwill [DR, (5.7)-(5.8), (5.13)-(5.15)] and restated in
    OpenAI's 11/12 paper, Appendix A, equations (cusp-fourier-coefficients) and
    (conjugate-cusp-coefficients); d_0(l) = conj(tau(-l));
  * g = (a b; C delta') in SL_2(O), H is the cusp normalizer chosen as in that appendix (cases 3 | C,
    v_lambda(C) = 1, and (C, lambda) = 1), g1 = g H^{-1} = (A1 B1; C1 D1) = I mod 3, and
    kappa(g1) = (C1/A1)_3 is the cubic Kubota character (kappa = 1 if C1 = 0);
  * the sign -1 = (-1)^j, the factor alpha(l)^j and the power conj(alpha(C))^{2j} are those for j = 1
    derivative.
The weight is V(y) = exp(-(log y)^2/(2 s^2)) with s = 0.5, so V^ is a Gaussian and V# is computed by
quadrature along Re t = 0.

The cusps cover all three types (3 | C, v_lambda(C) = 1, C coprime to lambda) and, for each cusp, every
value of kappa that occurs among the matrices g found.  At a few cusps the l and -l terms cancel and
LHS vanishes identically; there the script checks that RHS vanishes too ("degenerate cusps").

Controls, each of which must fail by an amount of size 1:
  * the sign flipped (+ instead of -);
  * the kernel Gamma(4/3 +- t) Gamma(5/3 +- t) of two zbar-derivatives in place of Q;
  * conj(alpha(C))^4 in place of conj(alpha(C))^2.

PASS criterion: relative error < 1e-6 at every nondegenerate cusp, |RHS| < 1e-8 at degenerate cusps,
and every control error > 1e-2.  Last run: worst relative error 5.5e-9 over 39 cusp configurations;
controls 2.0, 5.3, 2.0.

Usage: python3 test_sec7_reflection_j1.py [B]   (B = norm bound for the coefficient tables, default 80000).
Run time: about 75 seconds.

[DR] A. Dunn and M. Radziwill, Bias in cubic Gauss sums: Patterson's conjecture, Ann. of Math. 200 (2024).
"""
import os, sys, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastgauss import *
from scipy.special import loggamma
t0=time.time()
LAMp=(1,2); LAMc=cx(LAMp)
s_w=0.5; K=6.5; YCUT=4e5        # weight width; LHS cutoff exp(K s_w); RHS cutoff for V# (rapid decay)
def V(y):                       # smooth weight V(y) = exp(-(log y)^2/(2 s_w^2))
    return np.exp(-np.log(y)**2/(2*s_w*s_w))
def Q(t):                       # kernel for one zbar-derivative
    return np.exp(loggamma(7/6+t)+loggamma(5/6+t)-loggamma(7/6-t)-loggamma(5/6-t))
def Q2(t):                      # kernel for two zbar-derivatives (control)
    return np.exp(loggamma(4/3+t)+loggamma(5/3+t)-loggamma(4/3-t)-loggamma(5/3-t))
# quadrature nodes t = i*tau on Re t = 0; FV, FV2 hold V^(-t) Q(t) dt/(2 pi i) for the two kernels
TAU=np.arange(-14/s_w,14/s_w+1e-9,0.05); TT=1j*TAU
FV=np.sqrt(2*np.pi)*s_w*np.exp(s_w*s_w*TT*TT/2)*Q(TT)*0.05/(2*np.pi)
FV2=np.sqrt(2*np.pi)*s_w*np.exp(s_w*s_w*TT*TT/2)*Q2(TT)*0.05/(2*np.pi)
def Vnat(y,FVv=None):
    """V#(y) for an array y (kernel Q, or the control kernel if FVv=FV2)."""
    FVv=FV if FVv is None else FVv
    y=np.asarray(y,float); out=np.empty(len(y),complex)
    for i in range(0,len(y),4000):
        ly=np.log(y[i:i+4000]); out[i:i+4000]=np.exp(-np.outer(ly,TT))@FVv
    return out
def ce(z): return np.exp(4j*np.pi*np.real(z))
# ---------- cusp data ----------
def mod3eq(x,y): d=sub(x,y); return d[0]%3==0 and d[1]%3==0
def vlam(x):
    k=0
    while divides(LAMp,x): x=div(x,LAMp); k+=1
    return k
def cubic_general(x,a,plist):
    """(x/a)_3 for primary a (coprime to 3; may be even / non-squarefree)"""
    assert primary(a)
    r=1; m=a
    for p in plist:
        if N(m)==1: break
        while divides(p,m):
            r*=cubic_prime(x,p); m=div(m,p)
    assert N(m)==1,(a,m)
    return r
def matmul(A,B):
    (a,b),(c,d)=A; (e,f),(g,h)=B
    return ((add(mul(a,e),mul(b,g)),add(mul(a,f),mul(b,h))),(add(mul(c,e),mul(d,g)),add(mul(c,f),mul(d,h))))
def inv(A):
    (a,b),(c,d)=A; return ((d,(-b[0],-b[1])),((-c[0],-c[1]),a))
def cusp_data(a,C,plist):
    """For the cusp a/C: choose H as in OpenAI's 11/12 paper, Appendix A, and list (sigma, delta', kappa(g1)) for
    matrices g = (a b; C delta') in SL_2(O) with g H^{-1} = I mod 3, one for each value of kappa that occurs."""
    v=vlam(C); one=(1,0); zero=(0,0)
    if v>=2:
        assert mod3eq(a,one); H=((one,zero),(zero,one)); sigma='0'
    elif v==1:
        assert mod3eq(a,one)
        u0=LAMp if mod3eq(C,LAMp) else (-1,-2); assert mod3eq(C,u0)
        H=((one,zero),(u0,one)); sigma='-' if u0==LAMp else '+'
    else:
        assert primary(C)
        u0=(a[0]%3,a[1]%3)                  # fixed representatives of O/3
        eps={0:0,1:1,2:-1}[u0[1]]
        H=((u0,(-1,0)),(one,zero)); sigma={0:'0',1:'-',-1:'+'}[eps]
    found={}
    R=14
    for x in range(-R,R+1):
        for y in range(-R,R+1):
            dp=(x,y)
            t=sub(mul(a,dp),one)
            if not divides(C,t): continue
            bg=div(t,C); g=((a,bg),(C,dp))
            g1=matmul(g,inv(H))
            (A1,B1),(C1,D1)=g1
            if mod3eq(A1,one) and mod3eq(D1,one) and mod3eq(B1,zero) and mod3eq(C1,zero):
                if N(A1)>4000: continue
                kap=1 if C1==(0,0) else cubic_general(C1,A1,plist)
                key=(round(kap.real,6),round(kap.imag,6)) if kap!=1 else (1.0,0.0)
                if key not in found or N(dp)<N(found[key][0]): found[key]=(dp,kap)
    assert found,(a,C)
    return [dict(sigma=sigma,delta=dp,kappa=kap,v=v) for dp,kap in found.values()]
# ---------- coefficient tables: Gauss sums of all squarefree c with N(c) <= B, and the cube factors d ----------
B=int(sys.argv[1]) if len(sys.argv)>1 else 80000
plist=primary_primes(B)
print("primes:",len(plist),"time",round(time.time()-t0,1)); 
gp={}
for p in plist: gp[p]=g_prime(p)
print("prime Gauss sums done",round(time.time()-t0,1))
sq=squarefree_list(B,plist)
L2=(-3,0); wL2=mul((0,1),L2); w2L2=mul((-1,-1),L2)
cdata=[]
for c,f in sq:
    g=1+0j; s1=s2=s3=1
    for p in f:
        g*=gp[p]; s1*=cubic_prime(L2,p); s2*=cubic_prime(wL2,p); s3*=cubic_prime(w2L2,p)
    for i in range(len(f)):
        for j in range(i+1,len(f)): g*=pairsym(f[i],f[j])*pairsym(f[j],f[i])
    # g(r,c)=conj((r/c)_3) g(c)
    cdata.append((cx(c),N(c),g,np.conj(s1)*g,np.conj(s2)*g,np.conj(s3)*g))
print("squarefree c:",len(cdata),"time",round(time.time()-t0,1))
Dmax=int(round(B**(1/3)))+2
dl=[]
for x in range(-2*Dmax,2*Dmax+1):
    for y in range(-2*Dmax,2*Dmax+1):
        if primary((x,y)) and N((x,y))**3<=B: dl.append((cx((x,y)),N((x,y))))
E9=cmath.exp(2j*math.pi/9)
# d_0 on lam^-3 O : list of (l, coef)
def table_d0(Nmax):
    """Indices l in lambda^{-3} O with N(l) <= Nmax and the coefficients d_0(l) = conj(tau(-l))."""
    L=[];Cf=[]
    for cc,nc,g1_,gl2,gwl2,gw2l2 in cdata:
        if nc>27*Nmax: continue
        for dd,nd in dl:
            base=nc*nd**3
            if base>27*Nmax: continue
            absdc=abs(dd)/abs(cc)
            for k in range(-3,40):
                if base*3.0**k>Nmax*(1+1e-12): break
                r=k%3
                if r==1: continue
                if r==2:
                    nu=(k+4)//3
                    vals=[(1,np.conj(gl2)),(-1,np.conj(gl2)),(W,E9**-1*np.conj(gwl2)),(-W,E9**-1*np.conj(gwl2)),(W*W,E9*np.conj(gw2l2)),(-W*W,E9*np.conj(gw2l2))]
                    fac=absdc*3**(nu/2+2)
                else:
                    nu=(k+3)//3
                    vals=[(1,np.conj(g1_)),(-1,np.conj(g1_))]
                    fac=absdc*3**(nu/2+5/2)
                lk=LAMc**k*cc*dd**3
                for u,tv in vals:
                    L.append(u*lk); Cf.append(np.conj(tv*fac))      # d_0(l)=conj(tau(l)), tau even
    return np.array(L),np.array(Cf)
def table_dpm(sign,Nmax):
    """Indices l in lambda^{-4} O with N(l) <= Nmax and the coefficients d_+(l) or d_-(l) at the other two cusps."""
    L=[];Cf=[]
    for cc,nc,g1_,gl2,gwl2,gw2l2 in cdata:
        if nc>81*Nmax: continue
        for dd,nd in dl:
            base=nc*nd**3
            if base>81*Nmax: continue
            absdc=abs(dd)/abs(cc); m0=LAMc**-4*cc*dd**3
            if sign=='+':
                # tau_2 at mu'=-u' lam^-4 c d^3 ; l=u' w^2 lam^-4 c d^3 ; d_+(l)=w^2 conj(tau_2(mu')) ce(l)
                for up,t2 in [(1,9*W*W*np.conj(gl2)),(W,9*E9**-1*np.conj(gwl2)),(W*W,9*W*E9*np.conj(gw2l2))]:
                    l=up*W*W*m0; L.append(l); Cf.append(W*W*np.conj(t2*absdc)*ce(l))
            else:
                # tau_1 at mu''=u' lam^-4 c d^3 ; l=-u' w lam^-4 c d^3 ; d_-(l)=w conj(tau_1(mu'')) ce(l)
                for up,t1 in [(1,9*W*np.conj(gl2)),(W,9*E9**-1*W*W*np.conj(gwl2)),(W*W,9*E9*np.conj(gw2l2))]:
                    l=-up*W*m0; L.append(l); Cf.append(W*np.conj(t1*absdc)*ce(l))
    return np.array(L),np.array(Cf)
print("tables ready",round(time.time()-t0,1))
cusps=[((1,0),(-3,0)),((4,0),(-3,0)),((1,0),pw(LAMp,3)),((7,0),(-6,0)),((-2,-3),(-6,0)),
       ((1,0),LAMp),((1,0),(-1,-2)),((1,0),mul((-2,0),LAMp)),((7,0),mul((2,0),LAMp)),((1,0),mul(LAMp,(-2,-3))),
       ((0,0),(1,0)),((0,1),(1,0)),((0,-1),(1,0)),((1,0),(-2,0)),((0,1),(-2,0)),((0,-1),(-2,0)),((1,1),(-2,-3)),((2,0),(-2,-3)),((0,2),(-2,-3))]
worst=0; cmax=[0,0,0]; zero_cases=[]
for a,C in cusps:
  for cd in cusp_data(a,C,plist):
      NC=N(C); Cc=cx(C); ac=cx(a)/Cc
      if cd['sigma']=='0': X=3.16*NC
      else: X=5.46*NC
      NL=X*math.exp(K*s_w); NR=YCUT*NC**2/((2*math.pi)**4*X)
      if 27*NL>B or (27 if cd['sigma']=='0' else 81)*NR>B:
          print("skip (B too small)",a,C); continue
      Ll,Cl=table_d0(NL)
      lhs=np.sum(Cl*ce(Ll*ac)*np.conj(Ll/np.abs(Ll))**1*np.abs(Ll)**-1*V(np.abs(Ll)**2/X))
      if cd['sigma']=='0': Lr,Cr=table_d0(NR)
      else: Lr,Cr=table_dpm(cd['sigma'],NR)
      dl_=cx(cd['delta'])
      nl=np.abs(Lr)**2
      core=np.conj(cd['kappa'])*np.conj(Cc/abs(Cc))**2*np.sum(Cr*(Lr/np.abs(Lr))**1*ce(-dl_*Lr/Cc)*nl**-0.5*Vnat((2*np.pi)**4*nl*X/NC**2))
      rhs=-core
      if abs(lhs)<1e-8:
          zero_cases.append(abs(rhs)); print(f'   degenerate cusp a={a} C={C}: LHS={abs(lhs):.1e}, RHS={abs(rhs):.1e} (both vanish: l and -l cancel)'); continue
      err=abs(lhs-rhs)/abs(lhs); worst=max(worst,err)
      # controls: wrong sign, kernel of two derivatives, wrong power of conj(alpha(C))
      c1=abs(lhs-core)/abs(lhs)
      c2=abs(lhs+np.conj(cd['kappa'])*np.conj(Cc/abs(Cc))**2*np.sum(Cr*(Lr/np.abs(Lr))**1*ce(-dl_*Lr/Cc)*nl**-0.5*Vnat((2*np.pi)**4*nl*X/NC**2,FV2)))/abs(lhs)
      c3=abs(lhs+np.conj(cd['kappa'])*np.conj(Cc/abs(Cc))**4*np.sum(Cr*(Lr/np.abs(Lr))**1*ce(-dl_*Lr/Cc)*nl**-0.5*Vnat((2*np.pi)**4*nl*X/NC**2)))/abs(lhs)
      cmax[0]=max(cmax[0],c1); cmax[1]=max(cmax[1],c2); cmax[2]=max(cmax[2],c3)
      print(f"a={a} C={C} v_lam={cd['v']} sigma={cd['sigma']} delta'={cd['delta']} kappa={complex(cd['kappa']):.3f} X={X:.1f} #L={len(Ll)} #R={len(Lr)}  LHS={complex(lhs):.10f}  RHS={complex(rhs):.10f}  rel.err={err:.1e}")
print("worst rel. error (nondegenerate cusps)",worst,"; degenerate cusps max |RHS|",max(zero_cases) if zero_cases else 0)
print("controls (largest error over cusps): sign flipped %.3f, two-derivative kernel %.3f, conj(alpha(C))^4 %.3f"%tuple(cmax),"time",round(time.time()-t0,1))
print("[PASS]" if worst<1e-6 and (not zero_cases or max(zero_cases)<1e-8) and min(cmax)>1e-2 else "[FAIL]","7.6 [guide Sec. 7, theta transformation] single-cusp reflection with one zbar-derivative (zeta case)")
