"""t08: sanity check of the conclusion of Theorem thm:78 (guide L1013-1017) for a few curves, at finite height.

  C98 L1015  L(E_d, s) has no zero in 11/8 < Re s (checked on 1.38 <= Re s <= 2.65, 0 <= Im s <= 40, by the argument
             principle with PARI's lfun), for d in {1, 2, -3, 5, 16}.  Falsifier: the same winding count on the box
             0.9 <= Re s <= 1.1, 0.5 <= Im s <= 30 around the central line Re s = 1 must find the zeros listed by lfunzeros.
This is only a finite-height numerical check of a statement that GRH would imply; it says nothing about the proof.
Run time: about 30 s.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmcommon import gp
from cmreport import claim

GPCODE = r"""
wind(L,a,b,c,e,n)={my(P=[],s,v,tot=0.,prev,cur,d);
  P=concat(P,vector(n,k,a+c*I+(b-a)*(k-1)/n));
  P=concat(P,vector(4*n,k,b+c*I+(e-c)*I*(k-1)/(4*n)));
  P=concat(P,vector(n,k,b+e*I-(b-a)*(k-1)/n));
  P=concat(P,vector(4*n,k,a+e*I-(e-c)*I*(k-1)/(4*n)));
  P=concat(P,[P[1]]);
  prev=arg(lfun(L,P[1])); my(mn=abs(lfun(L,P[1])));
  for(k=2,#P, v=lfun(L,P[k]); mn=min(mn,abs(v)); cur=arg(v); d=cur-prev;
      while(d>Pi,d-=2*Pi); while(d<-Pi,d+=2*Pi); tot+=d; prev=cur);
  [tot/(2*Pi), mn]};
"""
rows = []; ok = True; fal_ok = True
for d in [1, 2, -3, 5, 16]:
    out = gp(GPCODE + f"""
E=ellinit([0,0,0,0,{d}]); L=lfuninit(E,[1.6,1.1,42]);
r=wind(L,1.38,2.65,0.,40.,1200); print(r[1]); print(r[2]);
r2=wind(L,0.9,1.1,0.5,30.,1200); print(r2[1]); print(#select(x->x>0.5,lfunzeros(L,30)));
""", prec=19)
    fl = lambda z: float(z.replace(" ", ""))
    w1, mn, w2, nz = fl(out[0]), fl(out[1]), fl(out[2]), int(out[3])
    ok &= abs(w1) < 0.01 and mn > 1e-6
    fal_ok &= abs(w2 - nz) < 0.01 and nz > 0
    rows.append(f"d={d}: zeros in box = {w1:.3f} (min |L| on boundary {mn:.3g}); critical box winding {w2:.3f} vs lfunzeros {nz}")
claim("C98", 1015, "thm:78 conclusion (finite height)", "A",
      "L(E_d,s) != 0 on 1.38<=Re s<=2.65, 0<=Im s<=40 for d in {1,2,-3,5,16} (argument principle)", ok,
      "; ".join(rows) + "; falsifier: the same count finds the central zeros", falsifier_failed=fal_ok)
