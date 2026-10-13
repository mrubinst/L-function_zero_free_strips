"""Claims B53 (eq:weight) and B55 (dual angular factor alpha(l)^2, cusp prefactor -conj(alpha(c))^4/81):
end-to-end check of the theta transformation at single cusps with TWO zbar-derivatives (the CM case).

Adapted from tests_guide/test_sec7_reflection_j1.py (one derivative, zeta case; copied, not edited).  For a cusp a/C
(g = (a b; C delta') in SL_2(O), cusp normalizer and Kubota character kappa as in that script) it compares

  LHS = sum_{0 != l in lambda^-3 O} d_0(l) ce(l a/C) conj(alpha(l))^2 N(l)^{-1/2} V(27 N(l)/X)

with the reflected sum

  RHS = (+1) conj(kappa) conj(alpha(C))^4 sum_{0 != l in lambda^-4 O} d_sigma(l) alpha(l)^2 ce(-delta' l/C) N(l)^{-1/2}
        V#(N(l) X / N(C)^2),
  V#(x) = (1/2 pi i) int_{Re t=0} V^(-t) Gamma(4/3+t)Gamma(5/3+t)/(Gamma(4/3-t)Gamma(5/3-t)) ((2 pi)^4 x/27)^{-t} dt,

which is eq:weight verbatim.  The completed sum T of eq:T equals -(1/81) times LHS for the corresponding twist
(t4_theta.py, B47: -1 = conj(alpha(lambda^-3))^2 and 81 = 3^(5/2) 27^(1/2)), so RHS = LHS is the statement that the
fixed cusp prefactor is -conj(alpha(c))^4/81 and the dual angular factor is alpha(l)^2 (B55), with the kernel and the
scaling 27 of eq:weight (B53).  The sign +1 = (-(cv)^-2)^2 comes from the chain rule eq:chain.

Controls (each must fail at size >= 1e-2): sign flipped; one-derivative kernel Gamma(7/6+-t)Gamma(5/6+-t);
conj(alpha(C))^2 in place of ^4; alpha(l)^1 in place of alpha(l)^2; scaling (2 pi)^4 x without /27.
Weight V(y) = exp(-(log y)^2/(2 s^2)), s = 0.5 (V^ Gaussian).  X = 27 X~ with X~ the j=1 script's choice.
Usage: python3 t5_reflection_j2.py [B]   (B = norm bound of the coefficient tables, default 80000).  About 90 s.
"""
import os, sys, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastgauss import *
from scipy.special import loggamma
t0=time.time()
LAMp=(1,2); LAMc=cx(LAMp)
s_w=0.5; K=6.5; YCUT=4e5        # weight width; LHS cutoff exp(K s_w); RHS cutoff for V# (rapid decay)
def V(y):                       # smooth weight V(y) = exp(-(log y)^2/(2 s_w^2))
    return np.exp(-np.log(y)**2/(2*s_w*s_w))
def Q(t):                       # kernel for one zbar-derivative (zeta; used as a control)
    return np.exp(loggamma(7/6+t)+loggamma(5/6+t)-loggamma(7/6-t)-loggamma(5/6-t))
def Q2(t):                      # kernel for two zbar-derivatives (the CM kernel of eq:weight)
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

worst=0; cmax=[0,0,0,0,0]; zero_cases=[]; ncases=0; lines=[]
for a,C in cusps:
  for cd in cusp_data(a,C,plist):
      NC=N(C); Cc=cx(C); ac=cx(a)/Cc
      Xt=3.16*NC if cd['sigma']=='0' else 5.46*NC          # X~ of the j=1 script; the CM scale is X = 27 X~
      X=27*Xt
      NL=Xt*math.exp(K*s_w); NR=YCUT*NC**2/((2*math.pi)**4*Xt)
      if 27*NL>B or (27 if cd['sigma']=='0' else 81)*NR>B:
          print("skip (B too small)",a,C); continue
      Ll,Cl=table_d0(NL)
      lhs=np.sum(Cl*ce(Ll*ac)*np.conj(Ll/np.abs(Ll))**2*np.abs(Ll)**-1*V(27*np.abs(Ll)**2/X))
      if cd['sigma']=='0': Lr,Cr=table_d0(NR)
      else: Lr,Cr=table_dpm(cd['sigma'],NR)
      dl_=cx(cd['delta']); nl=np.abs(Lr)**2
      kap=np.conj(cd['kappa']); aC=np.conj(Cc/abs(Cc))
      arg=(2*np.pi)**4*nl*X/(27*NC**2)                     # ((2 pi)^4 x/27) with x = N(l) X/N(C)^2
      base=Cr*ce(-dl_*Lr/Cc)*nl**-0.5
      core=kap*aC**4*np.sum(base*(Lr/np.abs(Lr))**2*Vnat(arg,FV2))
      rhs=+core
      if abs(lhs)<1e-8:
          zero_cases.append(abs(rhs)); print(f'   degenerate cusp a={a} C={C}: LHS={abs(lhs):.1e}, RHS={abs(rhs):.1e}'); continue
      err=abs(lhs-rhs)/abs(lhs); worst=max(worst,err); ncases+=1
      c=[abs(lhs+core)/abs(lhs),
         abs(lhs-kap*aC**4*np.sum(base*(Lr/np.abs(Lr))**2*Vnat(arg,FV)))/abs(lhs),
         abs(lhs-kap*aC**2*np.sum(base*(Lr/np.abs(Lr))**2*Vnat(arg,FV2)))/abs(lhs),
         abs(lhs-kap*aC**4*np.sum(base*(Lr/np.abs(Lr))**1*Vnat(arg,FV2)))/abs(lhs),
         abs(lhs-kap*aC**4*np.sum(base*(Lr/np.abs(Lr))**2*Vnat(arg*27,FV2)))/abs(lhs)]
      cmax=[max(u,w) for u,w in zip(cmax,c)]
      cmin_case=min(c)
      print(f"a={a} C={C} v_lam={cd['v']} sigma={cd['sigma']} kappa={complex(cd['kappa']):.3f} X={X:.1f} #L={len(Ll)} #R={len(Lr)}  LHS={complex(lhs):.10f}  RHS={complex(rhs):.10f}  rel.err={err:.1e}  min control err={cmin_case:.2f}")
      lines.append(cmin_case)
from cmlib import report, falsifier
deg=max(zero_cases) if zero_cases else 0
detail=(f"{ncases} nondegenerate cusp configurations (all cusp types, all kappa values found), worst rel. error {worst:.1e}; "
        f"{len(zero_cases)} degenerate (LHS=0) with max |RHS| {deg:.1e}")
report("B53", max(worst/1e-6, deg/1e-8), 1.0, detail+"; kernel Gamma(4/3+-t)Gamma(5/3+-t) and scaling (2pi)^4 x/27 of eq:weight")
report("B55", max(worst/1e-6, deg/1e-8), 1.0, detail+"; dual angular factor alpha(l)^2 and prefactor +conj(kappa)conj(alpha(C))^4 on LHS, i.e. -conj(alpha(c))^4/81 on T")
falsifier("B55","sign flipped (prefactor +conj(alpha(c))^4/81 on T)",cmax[0],1e-2)
falsifier("B53","one-derivative kernel Gamma(7/6+-t)Gamma(5/6+-t)",cmax[1],1e-2)
falsifier("B55","conj(alpha(C))^2 instead of ^4",cmax[2],1e-2)
falsifier("B55","alpha(l)^1 instead of alpha(l)^2",cmax[3],1e-2)
falsifier("B53","(2 pi)^4 x without the factor 1/27",cmax[4],1e-2)
print("note: a control can coincide with the claim at special cusps (e.g. alpha(C) = +-1 makes conj(alpha(C))^2 = conj(alpha(C))^4);"
      " each falsifier above is required to fail at some cusp, and does (max error over cusps).")
print("time",round(time.time()-t0,1))
