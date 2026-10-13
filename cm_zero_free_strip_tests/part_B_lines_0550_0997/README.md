# Part B audit: lines 550-997 of the CM guide

This folder tests every equality, inequality, bound and table entry in lines 550-997 of
`cm-eleven-twelfths-guide.tex`, which is identical to the public `cm_zero_free_strip.tex`. Those lines are the sections
"Why the bound excludes zeros", "Poisson summation reveals the cubic Gauss sums", "The completed sums and the second
theta derivative", "The quadratic average and the completed bound" and "Remove the cubes, then reduce the number of rows".

## Result

There are **88 claims**: 61 exact identities (E), 16 explicit inequalities (I), 7 asymptotic or analytic statements (A)
and 4 cited results (C).

- **86 PASS** and **0 FAIL**.
- **2 UNTESTED**, both cited: B05 (L(E_d,s) is entire, by modularity) and B71 (the cube-removal inequality, OpenAI's
  Lemma 5.3).
- Every exact or inequality test that admits one has a falsifier: a deliberately wrong variant that is required to
  fail. All 53 falsifiers fail as required.
- The 7 class-A passes are bookkeeping or sanity checks only: B01, B02, B19, B33, B61, B64 and B85. So are the passes
  of the cited B37 and B62. None of them proves an asymptotic bound.

Exit status (2026-10-12): `cmlib.report` and `cmlib.falsifier` record failures, and an atexit hook makes a script
exit 1 if any check failed or any falsifier did not fail; `CM_TESTS_FORCE_FAIL=1` makes the first check fail, as a
self-test. `run_all.py` exits 1 if any script exits nonzero. Caveat: the random induction chains of B79 and B83 test
scale bookkeeping under the stated transition formulas, not the analytic transfer inequality.

**No discrepancy was found in lines 550-997.** The CM-specific normalizations agree with exact computation.

- **conj(alpha)^2 gamma_2 (B27).** The paired identity eq:paired holds with these coefficients, summed over all 432
  characters of the classes mod 72. The zeta coefficients conj(alpha) gamma_2 fail.
- **conj(alpha(lambda^-3))^2 = -1 and conj(alpha(b))^6 (B46, B68).** The cube inversion with conj(alpha(h))^3 fails.
- **The -1/81 normalization (B47).** Each term of T equals -(1/81) times the theta-coefficient form.
- **The two-derivative reflection (B53, B55).** The kernel Gamma(4/3+-t)Gamma(5/3+-t), the scaling (2 pi)^4 x/27, the
  dual angular factor alpha(l)^2 and the cusp prefactor -conj(alpha(c))^4/81 are all confirmed end to end at 41
  single-cusp configurations, with relative error at most 2.3e-8. Each of the following controls fails at size 1 or
  more: the sign flipped, the one-derivative kernel, conj(alpha(c))^2, alpha(l)^1, and the scaling without /27.

### Observations (not errors)

1. **B19, lines 623-624.** The paper says "there can also be nonzero frequencies from these diagonal pairs". When
   Phi^ is supported in [0, C], a diagonal pair (z1 = z2 = 1) has a frequency h != 0 only if H <= C N(e0) <= C b D. With
   H = D^(1+theta), that means D^theta <= C b. So these terms exist only for bounded D. The hedged "can" is correct
   but is vacuous asymptotically.
2. **Proposition, lines 553-554.** (Withdrawn 2026-10-12.) This observation said one theta suffices. That is wrong: a
   fixed theta gives S_d(D) << D^{11/12+5theta/12+eps}, i.e. only Re s > 17/12 + 5theta/12, so the hypothesis "every
   fixed 0 < theta <= 1/10" is needed, as the guide states.
3. **eq:local, lines 788-795.** The j = 0 entry is the active branch only. The inactive branch contributes the scalar
   1 - 1/N(pi). This is consistent with the appendix's active-set notation, and the finite sums are checked in B57.
4. **The helper test suite, outside this chunk.** Test 6.3b of `tests_guide/test_sec6.py` checks the exact Poisson
   formula, but all five of its cases are degenerate. In each one psi is nontrivial on the units, which makes both
   sides vanish identically, so the test only checks 0 = 0. Here B17 uses eight nondegenerate cases, with psi(1+omega)
   = 1, together with two degenerate controls. The helper's end-to-end test 6.10 is nondegenerate.

### What is not tested

- The analytic estimates eq:MS, eq:qls, eq:completed, eq:remove and the transfer estimate, and the final descent
  bound. Where it is possible, their sizes are checked on finite data as sanity checks (B62, B64, B85), with
  arbitrary loose criteria.
- The theta transformation itself (Kubota, Patterson, Dunn-Radziwill). It is tested only end to end at single cusps
  (B53, B55) and through the automorphy spot check (B37).
- The Kubota shift chi(h)^-2 in the local factors (B57).
- The statement that only the identity cusp has a constant term (B52).
- The separation of smooth kernels in "E << D^eps Sigma implies (MS)" (B33). Only the exponent bookkeeping is
  checked.

## Running

```sh
python3 run_all.py            # about 4 minutes; one line per claim; outputs saved to outputs/*.txt
python3 run_all.py --quick    # skips t5 (the 2-minute end-to-end reflection)
python3 t4_theta.py           # any script also runs on its own and prints RESULT|id|PASS/FAIL|... lines
```

The scripts need Python 3 with numpy, scipy, sympy and mpmath, plus PARI/GP (`gp`) for t1. The modules `zw.py`,
`common.py`, `fastgauss.py` and `theta.py` are unmodified copies from `tests_guide/`. `cmlib.py` adds vectorized sextic
symbols and Gauss sums, which are checked against the brute-force versions. `claims.py` is the registry of claims.
`t5_reflection_j2.py` is adapted from `tests_guide/test_sec7_reflection_j1.py`, changed from one derivative to two.

The conventions are the paper's. Elements a + b omega are stored as pairs (a, b), and primary means 1 mod 3. lambda =
1 + 2 omega = i sqrt3, and chi_n(x) = (x/n)_6 is extended by zero, for every power. The additive characters are
e(z) = exp(4 pi i Im z/sqrt3) and ce(z) = exp(2 pi i (z + conj z)), and alpha(n) = n/|n|. The theta coefficients are
transcribed independently from the paper's appendix tables. B37 checks them against `theta.py`.

## Claim table

The class column uses E = exact, I = inequality, A = asymptotic or analytic, and C = cited. The test column names the
script (t1 to t7). The full detail of each result, including errors and tolerances, is in `outputs/*.txt` and in the
output of `run_all.py`.

| id | line(s) | label | claim | class | test | range | result |
|---|---|---|---|---|---|---|---|
| B01 | 552-555 | Prop. | (MS) for all complex W and 0<theta<=1/10 implies L(E_d,s)!=0 for Re s>17/12 (logical chain; exponent 11/12+1/2=17/12) | A | t1: exponent arithmetic (Fraction) + chain of B02-B08 | logical; ingredients only | PASS |
| B02 | 563 | - | F_W(s)=int S_d(D)D^{-s}dD/D is holomorphic for Re s>11/12 (from S_d(D)<<D^{11/12+eta}); sanity: actual growth of S_d(D) | A | t1: integrability exponents (Fraction); PARI sanity of S_d(D) | 199 eta; d in {1,2,-3,7,16,-432,5}, D<=3e4 | PASS |
| B03 | 564 | - | S_d(D)=0 for D<1/b (W supported in [a,b], m>=1) | E | t1: direct evaluation of W(m/D) | D<1/b, m<1e4 | PASS |
| B04 | 565-568 | eq:mellin-zero | L(E_d,s+1/2)F_W(s)=W^(s)B_d(s) for Re s>1 (substitution + Dirichlet series vs PARI lfun and minimal-model local factors) | E | t1: mpmath substitution identity; PARI lfun x truncated series vs minimal-model B_d | 7 d x 3 s (Re s=3,2.5,1.7), M=6e4 | PASS (2 falsifiers fail) |
| B05 | 569-570 | - | L(E_d,s) entire, so the identity extends to Re s>11/12 (modularity; cited) | C | not tested: cited (modularity of E_d; L(E_d,s) entire); the identity-theorem step is standard complex analysis | - | UNTESTED |
| B06 | 573-577 | - | W(y)=y^{-s0}phi(y) gives W^(s0)=int_1^2 phi(y)dy/y>0 | E | t1: mpmath quadrature (30 digits) | 3 complex s0 | PASS (1 falsifier fail) |
| B07 | 578-579 | - | right side W^(s0)B_d(s0) nonzero: B_d holomorphic and nonzero for Re s>0 (local roots) | E | t1: PARI roots of local polynomials at p|6d | 7 values of d | PASS |
| B08 | 572-573,582 | - | Re rho>17/12 and s0=rho-1/2 give Re s0>11/12 | I | t1: Fraction | 1999 rationals | PASS (1 falsifier fail) |
| B09 | 592 | - | e(z)=exp(4 pi i Im z/sqrt3) is trivial on O and O is self-dual for (z,w)->e(zw) (needed for e(x/q) and Poisson) | E | t2: e on O; dual lattice of O under e(zw) | |a|,|b|<=30; 5329 points of (1/12)O | PASS (1 falsifier fail) |
| B10 | 593-594 | - | Phi^(N(w)) normalization: Phi(t)=exp(-pi t) gives (2/sqrt3)exp(-4 pi N(w)/3), radial | E | t2: 2D quadrature (20 digits) | 3 w + rotations | PASS (1 falsifier fail) |
| B11 | 596-597 | - | smooth mean square over all u (incl. 0) >= M(D,H) when Phi>=0 and Phi>=1 on [0,1] | I | t2: CM rows A_u with alpha, Phi=e^{1-t} | D=40, H=60 | PASS (1 falsifier fail) |
| B12 | 599 | - | only squarefree columns occur (mu_K(n)=0 otherwise) | E | t2: factorization | N<=400 | PASS |
| B13 | 599-600 | - | n_i=g z_i with g=gcd gives (z1,z2)=(g,z1 z2)=1 | E | t2: random squarefree pairs | 400 pairs | PASS (1 falsifier fail) |
| B14 | 601-605 | - | conj(chi_{n1}(u))chi_{n2}(u)=1_{(u,g)=1}psi(u), psi=conj(chi_{z1})chi_{z2} | E | t2: random (n1,n2,u) | 1200 triples | PASS (1 falsifier fail) |
| B15 | 606 | - | psi is primitive modulo z1 z2 | E | t2: local components of psi | 46 pairs, 99 components | PASS |
| B16 | 607-609 | - | normalized Gauss sum: \|gamma(psi)\|=1 for primitive psi | E | t2: vectorized Gauss sums | 46 pairs, N(q)<=3000 | PASS |
| B17 | 611-621 | eq:poisson-exact | exact Poisson formula with the exclusion divisor e0 (incl. g!=1, z1=1, z2=1, z1=z2=1); end-to-end CM mean square | E | t2: lattice sum vs dual sum, Phi=exp(-pi t); end-to-end CM mean square | 10 cases (8 nondegenerate); D=40,H=60 | PASS (4 falsifiers fail) |
| B18 | 622-623 | - | zero frequency contributes only when z1=z2=1 | E | t2: chi_n(0) | 60 n | PASS |
| B19 | 623-624 | - | nonzero frequencies can occur for diagonal pairs (requires H N(h)/N(e0)<=C_Phi) | A | t2: Gaussian Phi; support bookkeeping | g=pi7 pi13, H=20 | PASS |
| B20 | 625-629 | - | N(h)<<N(e0)N(z1z2)/H<<D^2/(H N(g)) | I | t2: sympy + 1e5 random sets | - | PASS (1 falsifier fail) |
| B21 | 630-631 | - | H>D gives frequency norms D^2/H<D | I | t2: Fraction | D<300 | PASS (1 falsifier fail) |
| B22 | 640 | - | gamma_2 is a cubic Gauss sum: chi_n(x)^2=(x/n)_3 | E | t3: Euler criterion (N-1)/3 | 500 (n,x) | PASS (1 falsifier fail) |
| B23 | 642-645 | eq:gauss | mu_K(n)gamma_{-1}(n)=chi_n(-1)G(n)^{-1}conj(alpha(n))gamma_2(n) | E | t3: vectorized Gauss sums (checked vs brute force) | 150 squarefree n, N<=529 | PASS (2 falsifiers fail) |
| B24 | 647-648 | - | \|G(n)\|=1 | E | t3 | 150 n | PASS |
| B25 | 648-649 | - | G depends only on the primary invertible class mod 72 | E | t3: two squarefree reps per class | all 432 classes mod 72 | PASS (1 falsifier fail) |
| B26 | 651-658 | - | the class group frak G mod 72 (\|G\|=432) and its characters xi; c_xi well defined | E | t3: characters from the regular representation | |frak G|=432 | PASS |
| B27 | 659-666 | eq:paired | mu(z1)mu(z2)conj(alpha(z1))alpha(z2)gamma(conj chi_{z1} chi_{z2})=sum_xi c_xi a_xi(z1)conj(a_xi(z2)), a_xi=conj(alpha)^2 gamma_2 xi | E | t3: literal sum over all 432 xi | 73 coprime pairs | PASS (3 falsifiers fail) |
| B28 | 667-670 | - | the second conj(alpha) comes from conjugating the first column: conj(alpha(n1))alpha(n2)=conj(alpha(z1))alpha(z2); CM LHS = conj(alpha(z1))alpha(z2) x zeta LHS | E | t3 | 73 pairs | PASS |
| B29 | 684-685 | - | e0 contributes chi_n(h e0^5): conj(chi_n(e0))chi_n(h)=chi_n(h e0^5) | E | t3: random | 300 triples | PASS (1 falsifier fail) |
| B30 | 686-690 | - | a_xi(v0 n)=a_xi(v0)a_xi(n)chi_n(v0)^4 for (v0,n)=1 | E | t3: random, two xi | 71 pairs | PASS (1 falsifier fail) |
| B31 | 691-692 | - | f=e0 v0, k=e0 h give h e0^5 v0^4=k f^4 (ring and character level) | E | t3: exact in O + characters | 300 triples | PASS (1 falsifier fail) |
| B32 | 694-698 | eq:initial | X=D/(BF), Sigma=XF=D/B, H_rows<<D^2/(H B^2) (from N(bf m_i)<=b0 D and supp Phi^) | I | t3: OAI change of variables on prime data; sympy | 300 samples | PASS (1 falsifier fail) |
| B33 | 699-701 | - | E<<D^eps Sigma in these ranges implies (MS): exponent bookkeeping sum_{N(b)~B}(H N(b)/D^2)(XF)^2 ~ H | A | t3: sympy bookkeeping | - | PASS |
| B34 | 717-719 | eq:T | the b=1 term of T(X;k,f) is X^{-1/2} times the column sum of (eq:energy) | E | t4: direct | X=1400, xi=(omega/n)_3 | PASS |
| B35 | 718 | eq:T | all sums in T are finite (N(n)N(b)^3<=bX) | E | t4: enumeration | X=1400 | PASS |
| B36 | 723 | - | lambda=1+2 omega=i sqrt3 | E | t4 | - | PASS |
| B37 | 724-727 | - | conj(theta) has terms d_0(l)vK_{1/3}(4 pi\|l\|v)ce(lz) [DR]; spot check: appendix table = helper transcription of DR (which passes automorphy) | C | t4: own transcription of appendix table vs theta.py; automorphy | 528 indices, |m|<=11 | PASS |
| B38 | 728-732 | eq:theta-coeff | d_0(lambda^{-3}n b^3)=3^{5/2}\|b\|conj((lambda/n)_3)gamma_2(n) (from the appendix table) | E | t4: appendix table vs RHS | 30 n x 7 b | PASS (1 falsifier fail) |
| B39 | 733-735 | - | restricting lambda^3 l primary coprime to 6 selects exactly l=lambda^{-3}nb^3, n squarefree | E | t4: appendix table | 363 x, N<=1200 | PASS (1 falsifier fail) |
| B40 | 735-737 | - | reflected indices lie in l=u lambda^m n b^3, u unit, m>=-4 | E | t4: d_0, d_+, d_- from the appendix | lambda^-4 m, |m|<=9 | PASS |
| B41 | 738-739 | - | d_0(lambda^{-4})=0 | E | t4: appendix table | - | PASS (1 falsifier fail) |
| B42 | 739-740 | - | (lambda/n)_3 depends only on n mod 9 | E | t4: Euler criterion | all primary n, N<=3000 | PASS (1 falsifier fail) |
| B43 | 740-741 | - | inserting (lambda/lambda^3 l)_3 cancels (lambda/n)_3: (lambda/n b^3)_3=(lambda/n)_3 | E | t4 | 60 n x 4 b | PASS |
| B44 | 743-744 | - | d^2/dzbar^2 ce(lz)=(2 pi i conj l)^2 ce(lz) | E | t4: sympy | - | PASS (1 falsifier fail) |
| B45 | 744-745 | - | Mellin removes the radial power leaving conj(alpha(l))^2 (Bessel-Mellin, two derivatives) | E | t4: mpmath quadrature (25 digits) | 3 l x 3 s | PASS (1 falsifier fail) |
| B46 | 746-748 | - | conj(alpha(l))^2=conj(alpha(lambda^{-3}))^2 conj(alpha(n))^2 conj(alpha(b))^6 and conj(alpha(lambda^{-3}))^2=-1 | E | t4 | 280 (n,b) | PASS (1 falsifier fail) |
| B47 | 749 | - | every angular power in (eq:T) is explained: T summand = -(1/81) d_0(l)(lambda/lambda^3 l)_3 conj(alpha(l))^2 N(l)^{-1/2} Psi(lambda^3 l) V(27N(l)/X) | E | t4: termwise, T vs theta form | 40 n x 7 b | PASS (1 falsifier fail) |
| B48 | 752-756 | - | z', v' of g^{-1}(a/c+z,v) | E | t4: sympy (ad-bc=1) + numeric SL_2(Z[omega]) | - | PASS (1 falsifier fail) |
| B49 | 757-759 | - | at z=zbar=0: pure zbar-derivatives of zbar', v' vanish, of z' of order>=2 vanish, d_zbar z'=-(cv)^{-2} | E | t4: sympy | orders 1-5 | PASS |
| B50 | 760-764 | eq:chain | d_zbar^2 F(z',zbar',v')\|_0=(cv)^{-4}(d_{z'}^2 F)(-delta/c,1/(N(c)v)) | E | t4: sympy, generic F | degree (4,2,2) | PASS (2 falsifiers fail) |
| B51 | 765 | - | no lower derivative survives (coefficients of F_{z'}, F_{zbar'}, F_v, mixed terms vanish) | E | t4: sympy, Faa di Bruno coefficients | - | PASS |
| B52 | 765-768 | - | the derivative kills the constant Fourier term (only the identity cusp has one: cited DR) | E | t4: sympy, generic phi(v') | - | PASS |
| B53 | 770-775 | eq:weight | reflected weight V#(x) with Gamma(4/3+-t)Gamma(5/3+-t) and ((2 pi)^4 x/27)^{-t}: end-to-end single-cusp reflection with two derivatives | E | t5: end-to-end reflection, two derivatives, kernel of eq:weight | 41 cusp configs, tables N<=8e4 | PASS (2 falsifiers fail) |
| B54 | 776 | - | CM gamma quotient = zeta quotient shifted by 1/2 (7/6,5/6 -> 5/3,4/3) | E | t4: Fraction + sympy | - | PASS |
| B55 | 777-779 | - | reflected angular factor alpha(l)^2 and fixed cusp prefactor -conj(alpha(c))^4/81 | E | t5: end-to-end reflection, alpha(l)^2 and conj(alpha(C))^4 | 41 cusp configs | PASS (3 falsifiers fail) |
| B56 | 785-786 | - | Psi(n)=xi(n)chi_n(k)chi_n(f)^4 = xi_0(n) prod_pi chi_pi(n)^{j_pi} with xi_0 a class function (tested mod 72) | E | t4: ratio as a function of n mod 72 | 6 (k,f), N(n)<=2500 | PASS |
| B57 | 787-795 | eq:local | local reflected factors B_{pi,j} (finite sums; the Kubota shift -2 is the cited input) | E | t4: reflected local sums (Kubota shift -2 cited) | 7 primes, j=0..5, 1854 cases | PASS (1 falsifier fail) |
| B58 | 796-799 | - | B_{pi,1}(x)=chi_pi(x)^3=(x/pi)_2 | E | t4: Euler criterion (N-1)/2 | 10 primes | PASS |
| B59 | 800-801 | - | squarefree row primes away from 6f have exponent j=1 | E | t4 | 6 (k,f) | PASS |
| B60 | 804-806 | - | effective dual range ~Q^2/X is shorter than X iff Q<X | I | t4: Fraction | Q<120, X<200 | PASS (1 falsifier fail) |
| B61 | 807-808 | - | V# decays rapidly (negligible tail when X/Q^2 is large) | A | t4: contour-shift bound |V#(x)|<=C_c (Kx)^-c | c=1,3,5; x<=1e4 | PASS |
| B62 | 812-822 | eq:qls | quadratic large sieve [GL]; sanity ratio LHS/((K+U)sum\|beta\|^2) | C | t6: sanity ratio | K=U=400 | PASS |
| B63 | 819-821 | - | squarefree restriction matters: square columns give repeated values 1 (no cancellation) | E | t6: beta on squares vs squarefree | K=U=400,1600,6400 | PASS |
| B64 | 824-835 | eq:completed | completed estimate sum_k\|T(X;k,f)\|^2<<D^eps(H+H^2N(f)/X); sanity with CM T | A | t6: sanity ratio with CM T | X=1400, H<=600, 3 f | PASS |
| B65 | 837-838 | - | every k!=0 is u s v0^2 with u unit, s squarefree | E | t6: factorization incl. lambda, 2 | all k, N(k)<=3000 | PASS |
| B66 | 839-841 | - | chi_n(v0^2)chi_n(f)^4=chi_n(f v0^2)^4 | E | t6: random + T identity | 400 + 8 cases | PASS (1 falsifier fail) |
| B67 | 842-844 | - | sum_{v0}N(v0)^{-2}<oo and the bookkeeping (H/N(v0)^2, f v0^2) -> N(v0)^{-2}(H+H^2N(f)/X) | E | t6: partial sums; sympy | R<=1e4 | PASS |
| B68 | 853-859 | eq:inversion | exact cube inversion with mu_K(h)conj(alpha(h))^6 Psi(h)^3/N(h) and T(X/N(h)^3) | E | t7: LHS vs sum over h of T | X=4000, 3 (k,f,xi) | PASS (2 falsifiers fail) |
| B69 | 860-863 | - | sum_{h\|b'}mu_K(h)=1_{b'=1}; alpha and Psi completely multiplicative | E | t7: divisor sums; random pairs | N(b')<=3000 | PASS |
| B70 | 865-876 | - | H<=Sigma, N(h)<=H_c: H+H^2 F N(h)^3/X<<Sigma (and X/N(h)^3>=1) | I | t7: exact rationals | 77237 cases | PASS (1 falsifier fail) |
| B71 | 877-883 | eq:remove | cube-removal inequality (OAI Lemma 5.3) | C | not tested: cited analytic estimate with unspecified constants (OAI Lemma 5.3); no finite test of its truth; its range conditions are tested in B70, B72 | - | UNTESTED |
| B72 | 884-885,956-957 | - | H^2<=X makes the supremum empty; a nonempty supremum has N(b)^3>X^2/H^2 | I | t7: exact rationals | 455878 cases | PASS (1 falsifier fail) |
| B73 | 894-895 | - | reverse paired identity: a(u1)conj(a(u2))gamma(chi_{u1}conj chi_{u2})=mu mu G(u1/u2) conj(alpha(u1))alpha(u2), i.e. Moebius coefficients mu_K conj(alpha) xi_1 | E | t3: reverse identity with xi-expansion | 73 pairs | PASS (1 falsifier fail) |
| B74 | 897-914 | eq:P, eq:Q | Y_{C,e0} with N(e0)/N(C)^2 suppressed is Y_enlarged; natural y-range <= Y_{C,e0} (c_*>=4C_Phi v^2); same normalizations as OAI | E | t7: sympy; comparison with OAI source | - | PASS |
| B75 | 916-918 | - | first Poisson frequency ~L^2/H; y=h f^2 gives ~L^2F^2/H | I | t7: sympy | - | PASS |
| B76 | 920-925 | eq:enlarge | Y_enl/Y_nat=Sigma/(LF)>=1 | I | t7: sympy + sampling | 5e4 samples | PASS |
| B77 | 926-928 | - | without enlargement H'_nat=H/(F^2R^2) (no reduction at F=R=1) | E | t7: sympy | - | PASS |
| B78 | 932-936 | - | local cancellation (1+1_{pi\|k'})-1_{pi\|k'}-1_{pi nmid k'}=1_{pi\|k'} | E | t7: both values | - | PASS (1 falsifier fail) |
| B79 | 937-940,952 | - | H'>=F'>=1 and child column scales >=1 (retained-support chain) | I | t7: sampling of the retained-support chain | 1e5 sets | PASS |
| B80 | 942-950 | eq:transfer | transfer ranges H'<=HL/(Sigma F), H'/Sigma'<=H/Sigma, Sigma'<=L (consistent with H'=HL/(Sigma F R^2), Sigma'=L/R); the estimate itself is cited | I | t7: sympy + sampling | 1e5 sets | PASS |
| B81 | 955-962 | eq:contract | H'<=H/(N(b)^3F^2)<H(H/XF)^2<=H D^{-2 kappa} | I | t7: sampling + sympy | 79906 sets | PASS (1 falsifier fail) |
| B82 | 963-964 | - | H'/Sigma'<=D^{-kappa} is preserved | I | t7: sampling | 79906 sets | PASS |
| B83 | 965-966 | - | finitely many steps reach H<=1 (random chains) | I | t7: random induction chains | 5e4 chains | PASS |
| B84 | 966-967 | - | H<=1 leaves at most six rows (the units); counting gives E<<X<=Sigma | E | t7: unit count; computed base-case averages | 3 (X,F) | PASS |
| B85 | 968-973 | - | E(H,X,F)<<D^eps Sigma for H<=Sigma D^{-kappa}; sanity with CM coefficients | A | t7: sanity E/Sigma with CM coefficients | X<=1200, H<=160 | PASS |
| B86 | 974-975 | - | initial ranges: H_rows/Sigma = D^{-theta}/B (up to constants) | E | t7: sympy | - | PASS (1 falsifier fail) |
| B87 | 976 | - | kappa=theta/2 works for D large (C D^{-theta}/B<=D^{-theta/2}) | I | t7: logarithms | 2e4 samples | PASS |
| B88 | 980-991 | table | table of zeta->CM changes (mu alpha; conj(alpha)^2 gamma_2; d_zbar^2; conj(alpha(b))^6; alpha(l)^2), each row tested with a zeta-form control that must fail | E | run_all: B17, B27, B46, B50, B53, B55, B68 and their zeta-form controls | - | PASS |
