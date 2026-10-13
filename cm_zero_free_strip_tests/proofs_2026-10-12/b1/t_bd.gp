\\ Block b1: local factors at p | 6d (Lemma b1-euler(c)).
bad=0; good23=0; goodp=0; nchecked=0; falsif=0; cmdone=0;
{
L=concat(vector(6001,i,i-3001), concat(vector(41,k,5^6*(k-21)), concat(vector(21,k,7^6*(k-11)), vector(11,k,2^6*3^6*5^6*(k-6))))); for(ii=1,#L, d=L[ii]; if(d==0,next);
  E=ellinit([0,0,0,0,d]); M=ellminimalmodel(E); gr=ellglobalred(E); N=gr[1];
  fa=factor(abs(6*d))[,1];
  for(i=1,#fa, p=fa[i]; nchecked++;
    isgood = (N%p!=0);
    ap = ellap(M,p);
    if(isgood,
       if(ap^2>4*p, bad++; print("Hasse fails ",d," ",p));
       if(p<=3, good23++; if(abs(ap)>p, bad++; print("|a_p|>p ",d," ",p)));
       if(p>=5, goodp++; v=valuation(d,p); if(v%6!=0, bad++; print("good but 6 !| v ",d," ",p));
          d1=d/p^v; if(ellap(ellinit([0,0,0,0,d1]),p)!=ap, bad++; print("ap mismatch ",d," ",p))),
       \\ bad reduction: local factor 1 - a_p T with a_p in {0,+-1}; for j=0 expect a_p=0
       if(p>=5 && valuation(d,p)%6==0, bad++; print("bad but 6 | v ",d," ",p));
       if(ap!=0, bad++; print("bad prime with a_p!=0 ",d," ",p)));
    \\ roots of L_p(E,s+1/2)^{-1} in X=p^{-s}: all |X|>=1 (no zeros in Re s>0)
    if(isgood, P=1-ap/sqrt(p)*x+x^2, P=1-ap/sqrt(p)*x);
    if(poldegree(P)>0, r=polroots(P); for(j=1,#r, if(abs(r[j])<1-1e-12, bad++; print("root inside ",d," ",p))));
    \\ falsifier: unshifted factor 1-a_p X + p X^2 has roots of modulus p^{-1/2}<1 at good p
    if(isgood, r=polroots(1-ap*x+p*x^2); if(vecmin(apply(abs,r))<1, falsif++))
  ));
print("checked (d,p|6d) pairs: ",nchecked,"  good p in {2,3}: ",good23,"  good p>=5 dividing d: ",goodp,"  failures: ",bad,"  falsifier detections: ",falsif);
cmdone=1;
}
\\ Exit status (added 2026-10-12): 0 only if the block above ran to the end (no GP error), found no failure, and the
\\ falsifier fired at every good prime.  CM_TESTS_FORCE_FAIL=1 forces status 1 (self-test of the wiring).
cmforce=(getenv("CM_TESTS_FORCE_FAIL")=="1");
cmst=if(cmdone && bad==0 && falsif==good23+goodp && !cmforce, 0, 1);
print("EXIT STATUS ",cmst,if(cmst,": a failure, a falsifier miss, a GP error, or CM_TESTS_FORCE_FAIL=1",": all checks passed"));
quit(cmst)
