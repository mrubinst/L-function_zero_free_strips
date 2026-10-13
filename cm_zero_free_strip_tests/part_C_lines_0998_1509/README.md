# Audit of guide lines 998–1509 (Part II: the 7/8 boundary for $y^2=x^3+d$)

Target: `cm-eleven-twelfths-guide.tex` (= public `cm_zero_free_strip.tex`), lines 998–1509: sections *Part II*,
*The probe for a CM target*, *Rows, zeros and prime slots*, *Counting bad rows with an angular factor*,
*The final inequality*, *What has been checked, and what remains*. Source of the adapted argument (OAI78):
`sources/openai/qrh-september-30-2026-paper.tex` (a copy of OpenAI's `build/paper.tex`) (16,677 lines; "L1234"
means that file). G78 = `zeta_seven_eighths/seven-eighths-guide.tex`. The paper was not modified.

## Result
*Re-run 2026-10-12 by the main agent:* `run_all.py` reports **123** claims (65 PASS, 35 CITED, 4 INFO, 19 UNTESTED,
0 FAIL). The counts below are from the agent's first run.


* **122 claims** recorded (every equality, inequality, bound, exponent, table entry, numerical value and quoted
  line number in the chunk; some statements appear twice, once as a restatement check against OAI78 and once as an
  independent computation). By class: **E 46, I 16, A 18, C 42**.
* **64 PASS, 0 FAIL, 35 CITED (restatement found in OAI78/G78 and consistent), 0 CITED-X, 4 INFO (sanity numbers for
  bounds with unspecified constants), 19 UNTESTED** (analytic statements and proof-inspection claims; reasons in the table).
* Every PASS has a falsifier (a deliberately wrong variant) that failed; `run_all.py` reports FAIL if a falsifier passes.
  Since 2026-10-12 each script also exits 1 if it reports a FAIL or CITED-X line (atexit hook of `lib/cmreport.py`),
  and `run_all.py` exits 1 if any script exits nonzero; `CM_TESTS_FORCE_FAIL=1` makes the first line of each script
  a FAIL, as a self-test. On 2026-10-13 `lib/cmreport.py` was pointed at the frozen copy
  `../cm-eleven-twelfths-guide_2026-10-09.tex` (the version whose lines these tests cite): C105 in `t07` counts the
  items of the old "Six statements must still be written out" list, which the current guide no longer has, so `t07`
  crashed when it read the live guide. With the frozen copy its output is the same as before.
  The source-text checks have their own falsifier (C62: wrong needles must not be found).
* **No genuine error was found.** The final exponent follows exactly: with $\sigma_0=7/8$, $C(s)=s-11/16$,
  $C(7/8)=3/16$ (the compensated low exponent), the certificate $-E(h)\ge 49/440640$ (re-derived symbolically from
  OAI78 eq:common-high-exponent; true minimum $2.2815\cdot10^{-4}$ at $x=1/2,\ \delta\approx0.3865$),
  $E_{\rm actual}(h)-\Delta\le-\tfrac{49}{440640}-\tfrac{51}{64}\Delta+\epsilon_{\rm total}<0$, and the removal of
  $h\Delta/4$ when $\kappa=3/4$ (the loss is exactly OAI78's capacity comparison, whose baseline $2/9=1/(6\cdot\tfrac34)$ is
  $\kappa=3/4$). Hence $\beta_*\le 7/8$ and, via the verified identity $L(E_d,s+\tfrac12)=L(s,\psi_d^*)$ up to
  nonvanishing local factors, $L(E_d,s)\ne0$ for $\Re s>7/8+1/2=11/8$ — *conditional on the six statements of
  L1484–1496*, as the guide itself says.

### Observations (not errors)

1. **L1022–1023, scope.** "$\psi_d^*(-p)=-1$ for $p\equiv2\pmod 3$" needs $p\nmid d$ (the split clause says $p\nmid6d$,
   the inert clause does not repeat it). For $p\mid d$ the value is $0$, e.g. $d=5,p=5$: $\psi_5^*(-5)=0$ (C03b). Harmless:
   those factors are among the deleted ones.
2. **L1174, "error below $4\cdot10^{-15}$".** Reproduced at the 60 primes of smallest norm (max $3.3\cdot10^{-15}$);
   over all 165 primes of norm $\le1000$ the double-precision error reaches $7.6\cdot10^{-15}$. Implementation detail only.
3. **L1202, "exactly as in lem:unmarked-completed-moment".** OAI78 states the bound per powerful-part dyad,
   $Z^{O/2+\max\{M-O,2M-O-N\}+\epsilon}$; the guide's $Z^\varepsilon(Q+Q^2/Z)$ is the sum over dyads at $N=1$ (C48).
   A simplification, not a change.
4. **L1234.** The guide writes $\min(x^{1/4},x^{-A})$ for OAI78's $\min\{x^{1/4},(1+x)^{-A}\}$; equivalent.
5. **L1389–1390.** (MS) is a special case of eq:cm-inverse only with the margin $c_1\le\vartheta$ (fixed for each
   fixed $\vartheta\in(0,1/10]$) and $c_2\le1+3\vartheta$; fine, since OAI78 allows any fixed positive $c_1,c_2$ (C50).
6. **L1396.** The CM plain moment is stated for both $e=\pm1$ with the fixed orientation $\chi_\bullet(k)$, i.e. it also
   covers anti-oriented products; that is consistent with the guide's remark that orientation matters only for the
   theta step of the inverse moment.
7. **L1126–1129** (constant term only at the identity cusp) could not be checked without Dunn–Radziwiłł App. A (C68).

## What was tested, briefly

* *Curves (t01):* point counts vs $\psi_d^*$ for 2352 split and 854 inert pairs; $a_m(E_d)/\sqrt m=\sum_{N(n)=m}\psi_d^*(n)$
  against PARI for 8 curves, $m\le600$; local factors at $p\mid6d$; $j$-invariants via `qfbclassno`/`polclass`.
* *Gauss sums (t02):* eq:signal-cm, its companion, OAI78 eq:signal, the co-oriented and anti-oriented classes, the weight
  identity eq:cm-weight, $b_p^3=a_p^2$, marks, sixth powers — all by brute-force Gauss sums at 165 primes and 60
  squarefree composites.
* *Reflection table (t03):* Bessel–Mellin for one and two $\bar z$-derivatives (rel. err $10^{-30}$), the chain rule for
  $\partial_{\bar z}^2$ through the cusp map (sympy), and a full numerical reflected/direct ratio reproducing OAI78's
  $-\tfrac{i}{81}\bar\alpha(c)^2R(t)\alpha(\mu)$ (calibration) and the CM $-\tfrac1{81}\bar\alpha(c)^4R_{\rm CM}(t)\alpha(\mu)^2$;
  pole sets, Stirling exponents, and the $V^\sharp$ bound.
* *Growth (t04):* conductors of $\psi_d^*$ by brute force give $3N(\mathfrak f)=N(E_d)$ for 13 curves; the functional
  equation with $\Gamma(s+\tfrac12)$ and $(3Q)^{s/2}$ against PARI; exponent $6/5$; $Q^{3/5}$; the theta Poisson step;
  $|1/L^S|\le\zeta_K(2)$; dual profile; conductor change only at $\lambda$; closure under finite-order twists.
* *Rows/moments (t05):* label set $\mathcal X_u$ (model), co-orientation, vanishing lattice main terms, exponent algebra of
  Prop. cm-completed, (MS) as a special case; sanity numbers for $\sum|B_m|^2$ and the inverse moment.
* *Final inequality (t06):* all exponent bookkeeping in exact arithmetic (Fractions/sympy), including the OAI78 certificate
  identity (completed square) and its positivity argument, on top of an exact rational grid.
* *Sources (t07):* every quoted OAI78 formula, label and line number (L2331, L3871, L4056, L2003–2279, the range table)
  and every G78 section reference.
* *Theorem sanity (t08):* argument principle: no zeros of $L(E_d,s)$ in $1.38\le\Re s\le2.65$, $0\le\Im s\le40$ for
  $d\in\{1,2,-3,5,16\}$ (falsifier: the same count recovers the central zeros listed by `lfunzeros`).

## Running

```sh
python3 run_all.py          # all nine scripts, about 2 minutes (t04 ~80 s); needs numpy, sympy, mpmath and PARI/GP `gp`
python3 run_all.py --quick  # skips t04 and t08
```
One line per claim: `[STATUS] <id> L<guide line> <label> {class} <claim> :: <details>`. Outputs of the last run are in
`outputs/` (`outputs/run_all.txt` is the full log). `lib/` holds copies of `zw.py`, `common.py`, `fastgauss.py` from
`tests_guide/` (unchanged) plus `cmcommon.py` and `cmreport.py`. Conventions are OAI78's: $\omega=e^{2\pi i/3}$,
primary $\equiv1\bmod3$, $\chi_n(x)=(x/n)_6$ with zero values retained, $e(z)=\exp(4\pi i\,\Im z/\sqrt3)$,
$\gamma_j(c)=N(c)^{-1/2}\sum_{v}\chi_c(v)^je(v/c)$, $G(c)=\overline{\chi_c(4)}\gamma_3(c)$, $\alpha(n)=n/|n|$.

Classes: **E** exact identity, **I** explicit inequality/exponent, **A** asymptotic/analytic, **C** cited from OAI78/G78.
Statuses: PASS, FAIL, CITED (restatement found in the source and consistent), CITED-X, INFO, UNTESTED.

## Claim → test table

| line | label | claim | class | test | range | result |
|---|---|---|---|---|---|---|
| 1011 | thm:78 | L(s, abar*eta) != 0 for Re s > 7/8, every finite-order eta (hence L(E_d,s) != 0 for Re s > 11/8) | A | t09_untested.py U01 | - | UNTESTED: the theorem is labelled 'claimed'; it is conditional on the six statements of L1484-1496. Only finite-height numerics (C98) and the exponent bookkeeping (C51-C6 |
| 1015 | thm:78 conclusion (finite height) | L(E_d,s) != 0 on 1.38<=Re s<=2.65, 0<=Im s<=40 for d in {1,2,-3,5,16} (argument principle) | A | t08_theorem_sanity.py C98 | 5 curves, Im s<=40 | PASS |
| 1015 | thm:78 | 7/8+1/2 = 11/8 and 1/2 < 11/8 < 3/2 | I | t01_curves.py C01 | exact Fractions | PASS |
| 1022 | psi*_d | sqrt(p)(psi*(pi)+psi*(pibar)) = a_p(d), 2352 split pairs (d,p), p<700, 20 d | E | t01_curves.py C02 | 20 d, split p<700, both primes (2352 pairs) | PASS |
| 1022 | psi*_d(-p) | psi*_d(-p) = -1 and a_p(d)=0 for 854 inert pairs (p<500) | E | t01_curves.py C03 | 20 d, inert p<500 (854 pairs) | PASS |
| 1023 | L(E_d,s+1/2)=L(s,psi*) | a_m/sqrt(m) = sum_{N(n)=m} psi*_d(n), (m,6d)=1, m<=600, d in [1, 2, -3, 7, 16, 5, -432, -1] (1531 coefficients) | E | t01_curves.py C04 | 8 curves, m<=600 coprime to 6d (1531 coeffs), PARI ellan | PASS |
| 1023 | scope of psi*_d(-p) = -1 | the inert statement needs p not dividing d (as in the split case): for p \| d the value is 0 | E | t01_curves.py C03b | p\|d, 6 cases | INFO |
| 1024 | Euler factors | omitted local factors at p\|6d have no zeros in Re s>0 (PARI a_p, conductor) | E | t01_curves.py C05 | 18 local factors (PARI) | PASS |
| 1025 | counts | numbers of checked pairs at least 480 split / 112 inert | E | t01_curves.py C06 | counts | PASS |
| 1028 | psi_d | psi*_d = conj(psi_d) and alpha conj(chi_n(4d)) = alpha chi_n((4d)^5) | E | t01_curves.py C07 | 6 d, 677 primary n | PASS |
| 1029 | same L | psi_d and psi*_d have equal Euler factors at split and inert primes | E | t01_curves.py C09 | 6 d, split p<600, inert p<200 | PASS |
| 1030 | psi*(pibar) | psi*_d(pibar) = conj psi*_d(pi) (4d rational) | E | t01_curves.py C08 | 2352 split pairs | PASS |
| 1032 | classification | every unitary Hecke character of K with infinity type (zbar/\|z\|)^{+-1} is abar*eta or alpha*etabar; covers all weight-2 CM-by-K newforms | C | t01_curves.py C12 | - | UNTESTED: standard theory (quotient of two characters with the same infinity type has trivial infinity type, hence finite order); no finite test. Consistent with C04 (E_d |
| 1036 | j-invariants | orders Z+fO_K with h=1 are f=1,2,3 with j = 0, 54000, -12288000 | E | t01_curves.py C10 | f<=300 (PARI) | PASS |
| 1040 | 16,700 lines | OAI78 source has about 16,700 lines | C | t07_source_crossrefs.py C63 | 16677 lines | CITED |
| 1046 | section-by-section check | each proof of OAI78 uses only properties CM characters share; no new idea needed | C | t09_untested.py U02 | - | UNTESTED: requires re-reading 16,677 lines of OAI78; not redone here (prior audits: notes78/audit_*.md). |
| 1065 | closure under finite-order twists | psi*_d * chi_.(lambda) has infinity type zbar/\|z\| (a modulus g with value abar(a) on a = 1 mod g exists), so its primitive inducer lies in F_CM | E | t04_growth_fe.py C99 | 3 d | PASS |
| 1068 | conjugation | L(s,conj psi) = conj L(conj s, psi) for psi = abar*chi(2+omega) (truncated, s=3+1.7i) | E | t01_curves.py C11 | N(n)<=3000, s=3+1.7i | PASS |
| 1071 | no poles | every L(s,psi), psi in F_CM, is entire | A | t09_untested.py U03 | - | UNTESTED: Hecke's theorem; the mechanism (no constant term in the theta transform) is tested in C41 and C101. |
| 1075 | thm:main | finite-order Hecke L-functions over Q(sqrt-3) have no zero in Re s>7/8 | C | t07_source_crossrefs.py C64 | OAI78 L106, L107, L108 | CITED |
| 1082 | continuation criterion | OAI78 Prop. lem:continuation-criterion holds verbatim with eta -> psi | C | t09_untested.py U04 | - | UNTESTED: proof inspection; its three inputs are tested: \|1/L^S\| <= zeta_K(2) (C42), nonvanishing local factors (C05), beta_* (def.). |
| 1082 | continuation criterion | Prop. lem:continuation-criterion; reciprocal Euler product bounded on Re s>=2; sigma_0=7/8, C(s)=s-11/16 | C | t07_source_crossrefs.py C65 | OAI78 L400, L450, L514 | CITED |
| 1083 | G78 references | G78 Prop. 1, Sections 3.1, 3.2, 3.3 (bin ceiling, witnesses, sextic sieve), 4.2, 4.4 (Q_i(u)), 5.1, 6 are the cited items | C | t07_source_crossrefs.py C97 | G78 Prop. 1: L137 | CITED |
| 1087 | \|1/L^S\| <= zeta_K(2) | sup over unimodular psi of \|prod(1-psi(p)N(p)^-s)\| on Re s>=2 is prod(1+N(p)^-2) = zeta_K(2)/zeta_K(4) < zeta_K(2) | I | t04_growth_fe.py C42 | Re s=2, t<=60; primes N<=2e5 | PASS |
| 1092 | sigma_0, C(s), 3/16 | sigma_0=7/8, C(s)=s-11/16, C(7/8)=3/16 = l_x/2+b/12; design b=1/8, ell=1/6, h=13/16, l_x=17/48, l_y=23/48 consistent | I | t06_final_inequality.py C51 | exact | PASS |
| 1107 | probe weight | eta enters the probe only through gamma_2(c) abar(A) Psi(A), Psi(a) = eta(a) x fixed phases | C | t07_source_crossrefs.py C66 | OAI78 L3386, L3363, L3397 | CITED |
| 1112 | eq:cm-weight | weight identity on 425 pairs (c sf, n), eta = chi_.(2+omega) | E | t02_gauss.py C13 | 425 pairs (c,n) | PASS |
| 1125 | target alpha*etabar | with target alpha*etabar the weight is gamma2(c) etabar(cn^3) (no angular factor); with abar*eta it carries abar^2 | E | t02_gauss.py C104 | 250 pairs | PASS |
| 1126 | cusp constant terms | theta has a nonzero constant term at the identity cusp and none at the other two cusps | C | t07_source_crossrefs.py C68 | - | UNTESTED: needs the cusp expansions of Dunn-Radziwill App. A (F_1, F_19, F_10); OAI78 L1936-1940 only defines C_H(v) |
| 1129 | L2331 | OAI78 L2331: the zbar-derivative removes the constant term C_H(v) | C | t07_source_crossrefs.py C67 | OAI78 L2331 | CITED |
| 1140 | rows | row u = zeta^S(6z) L^S(w,chi(u)) / L^S(x, eta conj chi(u)) H_{eta,u} (eta -> psi in the guide) | C | t07_source_crossrefs.py C69 | OAI78 L4045, L4046, L3831 | CITED |
| 1146 | Prop. cm-rows properties | H_{psi,u} converges in OAI78's regions, H << N(u)^eps, H_{psi,1} = 1+O(P_0^{-4/5}) | C | t09_untested.py U05 | - | UNTESTED: needs the full local Euler factor of OAI78 lem:local-euler; restatement checked (C70), formal substitution checked (C15, C16, C27). |
| 1148 | L4056, H = 1 + O(P_0^-4/5) | second region (L4056); H << q_u^eps; H_{eta,1} = 1+O(P_0^{-c_H}), c_H = 4/5 | C | t07_source_crossrefs.py C70 | OAI78 L4056, L4060, L4064 | CITED |
| 1151 | local table formal in psi(p) | every entry of OAI78's local table holds with eta(p) -> psi(p) | C | t09_untested.py U06 | - | UNTESTED: only the entries quoted by the guide are tested (C16, C27, C73); the full six-row table J_0..J_5 is not re-derived. |
| 1154 | L3871 | OAI78 L3871: no ray-class restriction on p | C | t07_source_crossrefs.py C71 | OAI78 L3871 | CITED |
| 1158 | eq:signal-cm | abar(c)gamma2(c)gamma1(c) = mu(c)G(c) at 165 primes (5 inert) and 60 squarefree composites | E | t02_gauss.py C14 | 165 primes (5 inert), 60 sf composites | PASS |
| 1160 | eq:signal | gamma1 gamma2 = mu alpha G for squarefree c (rearranged in eq:signal-cm) | C | t07_source_crossrefs.py C72 | OAI78 L967, L964 | CITED |
| 1162 | weight + signal | abar(c)^2 gamma2 eta gamma1 = mu G psi, psi = abar*eta, eta = chi_.(2+omega) | E | t02_gauss.py C15 | 225 sf c | PASS |
| 1166 | local table j=0 | -psi(p)conj chi_p(u) N(p)^-x is the p-coefficient of 1/L(x,psi conj chi(u)) | E | t02_gauss.py C27 | split p<600 | PASS |
| 1168 | a_p, b_p | a_p = abar^3 eta^3, b_p = abar^2 eta^2 conj chi_p(4), b_p^3 = a_p^2; J_0 = -eta(p) rho^-1 Q^-x + ... | C | t07_source_crossrefs.py C73 | OAI78 L3863, L3864, L3865, L4025 | CITED |
| 1170 | b_p^3=a_p^2 | chi_p(4)^3 = 1 and b_p^3 = a_p^2 with psi = abar*eta at 165 primes | E | t02_gauss.py C16 | 165 primes | PASS |
| 1174 | 60 primes, err<4e-15 | error of eq:signal-cm at the 60 primes of smallest norm below 4e-15 (double precision; the size of the rounding error depends on the implementation) | E | t02_gauss.py C17 | 60 smallest primes | PASS |
| 1181 | principal row | chi_n(1) = 1 for all primary n (N<=2000): row u=1 numerator is zeta_K | E | t02_gauss.py C18 | N(n)<=2000 | PASS |
| 1182 | c_S | principal data: H_eta = H_{eta,1}(s,1,1/6), c_S = W1^(1) M(1/6) (Res zeta^S)^2/6 (positive iff M(1/6)>0) | C | t07_source_crossrefs.py C74 | OAI78 L5520, L5521 | CITED |
| 1192 | B_m = C_V with abar -> abar^2 | OAI78's completed row C_V(X;m,nu) | C | t07_source_crossrefs.py C75 | OAI78 L2822 | CITED |
| 1193 | B_m weight; table row 'angular factor, original side' | gamma2(c)abar(cn^3)^2 nu(cn^3) = eq:cm-weight (nu=eta); abar(n)^2abar(b)^6 = abar(nb^3)^2 | E | t02_gauss.py C19 | 425 pairs | PASS |
| 1202 | Prop. cm-completed (sanity) | sum_{N(m)<=Q}\|B_m\|^2/(Q+Q^2/Z), Z=1500, 428 completed terms, abar^2 (CM) vs abar (OAI78) | A | t05_rows_moments.py C49 | Z=1500, Q<=3000 | INFO |
| 1202 | Prop. cm-completed exponent | O/2 + max(M-O, 2M-O-N) <= max(M, 2M-N) at N=1 for 0<=O<=2M (summing O-dyads gives Z^eps(Q+Q^2/Z)) | I | t05_rows_moments.py C48 | M in [0,3] grid | PASS |
| 1202 | Prop. cm-completed | sum_{N(m)<=Q}\|B_m\|^2 << Z^eps(Q+Q^2/Z) for the CM weight abar^2 | A | t09_untested.py U07 | - | UNTESTED: an analytic estimate (needs the CM transformation formula for d_zbar^2 theta at general moduli: item 1 of L1486); exponent algebra C48, sanity numbers C49, refl |
| 1202 | lem:unmarked-completed-moment | bound Z^{O/2+max(M-O,2M-O-N)+eps} (=> Z^eps(Q+Q^2/Z), see C48) | C | t07_source_crossrefs.py C76 | OAI78 L3139, L3159 | CITED |
| 1205 | marked lemmas (CM) | CM forms of lem:hybrid-energy, lem:reflected-energy, lem:probe-row-norm | A | t09_untested.py U08 | - | UNTESTED: stated by the guide as not written out (item 1). |
| 1207 | marked lemmas | lem:hybrid-energy, lem:reflected-energy, lem:probe-row-norm, Cor. cor:marked-reflection exist | C | t07_source_crossrefs.py C77 | OAI78 L7510, L7747, L8111, L7361 | CITED |
| 1215 | two derivatives through the cusp map | dbar of the composite at z=0 is -(cv)^-2 F_w (OAI78 L2324) and dbar^2 is +(cv)^-4 F_ww (generic non-holomorphic F) | E | t03_reflection.py C30 | sympy, exact point | PASS |
| 1221 | table, OAI78 column | dbar theta; Gamma(s+1/3)Gamma(s+2/3); R(t)=prod Gamma(1+t+-1/6)/Gamma(1-t+-1/6); -i/81 abar(c)^2; abar(n)abar(b)^3; alpha(mu) | C | t07_source_crossrefs.py C78 | OAI78 L2296, L2299, L1741, L1833, L1729, L1768 | CITED |
| 1223 | table: direct Mellin factor | Bessel-Mellin for one and two zbar-derivatives; CM factor Gamma(w+5/6)Gamma(w+7/6) in OAI78's variable w=s-1/2 | E | t03_reflection.py C28 | 3 (mu,s) points, 30 digits | PASS |
| 1224 | table: reflected kernel and cusp scalar | reflected/direct = -(1/81)abar(c)^4 q^-1/2 (Kq/q_c^2)^-t R_CM(t) alpha(mu)^2; OAI78 j=1 analogue -(i/81)abar(c)^2 R_1(t) reproduced as calibration | E | t03_reflection.py C31 | 2 mu x 2 s, 30 digits | PASS |
| 1228 | table: dual angular factor | dual side carries alpha(mu)^2 (CM) vs alpha(mu) (OAI78) | E | t03_reflection.py C32 | as C31 | PASS |
| 1232 | kernel bounds | \|R(sigma+iu)\| << (1+\|u\|)^{4 sigma} on -1/4<=sigma<=A; first pole -5/6; min(x^{1/4},(1+x)^{-A}) | C | t07_source_crossrefs.py C79 | OAI78 L2414, L2412, L2409, L1862 | CITED |
| 1232 | kernel growth | \|R(sigma+iu)\|/(1+\|u\|)^{4 sigma} bounded for sigma in {-1/4,0,0.3,1,2,3}, u in [0,1e6], both kernels | A | t03_reflection.py C34 | sigma in [-1/4,3], u<=1e6 | PASS |
| 1233 | holomorphy of kernels | CM kernel holomorphic for Re t > -4/3, OAI78's for Re t > -5/6; strip -1/4 <= Re t contains no pole | I | t03_reflection.py C33h | exact | PASS |
| 1234 | reflected weight bound | V#(x) << min(x^{1/4}, x^{-A}) for both kernels (Gaussian test, A=2, x in 1e-8..1e4) | A | t03_reflection.py C35 | x in [1e-8,1e4] | PASS |
| 1236 | dual lengths, sieve | dual lengths unchanged; alpha(mu)^2 absorbed by the quadratic large sieve | A | t09_untested.py U09 | - | UNTESTED: argument only; unimodularity of alpha(mu)^2 is trivial; no finite test of the sieve. |
| 1242 | marks | 1_{q\|A} = 1 - chi_q(A)^0 (zero convention), 1260 pairs | E | t02_gauss.py C20 | 1260 pairs | PASS |
| 1242 | marks | chi_p(x)^0 = 1_{p not \| x} (zero convention) | C | t07_source_crossrefs.py C80 | OAI78 L1674 | CITED |
| 1243 | phases | abar(ab)^4 = abar(a)^4 abar(b)^4 | E | t02_gauss.py C21 | 190 pairs | PASS |
| 1246 | balanced low estimate Z^{1/4} | C_I(11/12) = 11/12 - 2/3 = 1/4 = (1/2)(1/2) (X^{1/2}, X=Z^{1/2}) | I | t06_final_inequality.py C52 | exact | PASS |
| 1246 | low estimates | \|I_psi(Z^1/2,Z^1/2,Z)\| << Z^{1/4+eps} and \|I_mod(Z)\| << Z^{3/16+eps} for the CM target | A | t09_untested.py U10 | - | UNTESTED: analytic estimates; exponents checked in C51, C52, C81. |
| 1247 | Z^{3/16} | compensated low estimate Z^{l_x/2+b/12+eps} = Z^{3/16+eps} | C | t07_source_crossrefs.py C81 | OAI78 L8570, L509 | CITED |
| 1251 | rel. err 3e-26 | Bessel-Mellin identities agree to ~1e-26 or better at 30 digits | E | t03_reflection.py C29 | 30 digits | PASS |
| 1254 | Stirling exponents | gamma-quotient exponents: kernels 4 sigma; reflected/direct products 2-4 sigma (OAI78) and 4-4 sigma (CM, J-variable) | E | t03_reflection.py C36 | u=1e5..1e6 | PASS |
| 1263 | 3Q in Lambda | 3 N(f_psi*) = N(E_d) for 13 curves (conductor of psi* by brute force over moduli 2^i lambda^j prod p^k) | E | t04_growth_fe.py C37 | 13 curves | PASS |
| 1265 | functional equation | Lambda(s)=(3Q)^{s/2}(2pi)^{-s}Gamma(s+1/2)L(s,psi) = eps Lambda(1-s), eps=+-1, 5 curves x 4 points (PARI values parsed in double precision) | E | t04_growth_fe.py C38 | 5 curves x 4 points | PASS |
| 1266 | lem:hecke-strip-growth | \|L(sigma+it,psi)\| << Q^{3/5}(3+\|t\|)^2 on [-1/10,11/10]; FE with Gamma(s); \|L(-1/10+it)\| << Q^{3/5}(3+\|t\|)^{6/5} | C | t07_source_crossrefs.py C82 | OAI78 L1429, L1430, L1440, L1450 | CITED |
| 1266 | strip growth (PL bookkeeping) | 6/5 <= 2 and \|s+2\| <= (31/30)(3+\|t\|) on the strip, so PL on Q^{-3/5}L/(s+2)^2 gives Q^{3/5}(3+\|t\|)^2 | I | t04_growth_fe.py C40 | t<=1000 (exact part); PARI t<=100 | PASS |
| 1274 | theta series of psi | phi(zeta a)zetabar = phi(a) for phi = psi/abar; sum_{a in O} phi(a) abar(a) N(a)^-s = 6 sum_n psi(n) N(n)^{1/2-s} (N<=1500) | E | t04_growth_fe.py C101 | N(a)<=1500 | PASS |
| 1276 | theta: no constant term | FT of zbar e^{-pi a\|z\|^2} = -i a^{-2} xibar e^{-pi\|xi\|^2/a} (a=1/2,1,2), zero at xi=0 | E | t04_growth_fe.py C41 | a=1/2,1,2; 4 xi | PASS |
| 1280 | exponent 6/5 and Q^{3/5} | \|Gamma(3/2-s)/Gamma(1/2+s)\| and \|Gamma(1-s)/Gamma(s)\| ~ (3+\|t\|)^{6/5} on Re s=-1/10; (1-2sigma)/2 = 3/5 | E | t04_growth_fe.py C39 | t<=1e6 | PASS |
| 1286 | lem:logarithmic-control | logarithmic control lemma, with principal clauses | C | t07_source_crossrefs.py C83 | OAI78 L1531, L1533 | CITED |
| 1286 | logarithmic control for F_CM | OAI78 lem:logarithmic-control holds for F_CM with the same exponents | A | t09_untested.py U11 | - | UNTESTED: item 2 of L1491; its growth input is tested (C37-C40). |
| 1292 | X_u and Theta | Theta = <eta, T^>, X_u = {nu chi(u)^sigma}, Q_psi << q_u; slot characters lie in X_u | C | t07_source_crossrefs.py C84 | OAI78 L4224, L4233, L4248, L15029 | CITED |
| 1299 | label set X_u | X_u = {theta chi(u)^+-1} u {abar theta conj chi(u), alpha theta chi(u)}: finite, conjugation-closed, contains numerator, denominator and conjugates (model Theta_0 = Z/6 x Z/2) | E | t05_rows_moments.py C45 | model group Z/6 x Z/2 | PASS |
| 1306 | detector with X_u | X_u, Lemma cm-growth, Q_psi << N(u) and non-principality are all the detector uses | A | t09_untested.py U12 | - | UNTESTED: item 3 of L1492; X_u properties tested (C45, C46), conductor sizes (C102). |
| 1308 | Q_psi << N(u) | conductor norm of psi*_d (row u = (4d)^5 ~ 4d) relative to N(4d) | A | t04_growth_fe.py C102 | 13 d | INFO |
| 1320 | co-orientation | all CM members of X_u and both witnesses have k*sigma = -1 (co-oriented) | E | t05_rows_moments.py C46 | model | PASS |
| 1327 | compensated slots | conj(psi(p)) times the denominator's local term -psi(p)conj chi_p(u) equals -conj chi_p(u) (u=1: -1 per prime, so A_T(Z) unchanged) | E | t02_gauss.py C22 | split p<800, 5 rows u | PASS |
| 1338 | P_i^{delta/2} | OAI78 lem:bin-prime-bound exponent a-1/2 equals delta/2; q<=delta/2 <=> x in [0,1/2] | I | t06_final_inequality.py C54 | sympy | PASS |
| 1338 | lem:bin-prime-bound | \|Q_i(u;z)\| << U^eps P_i^{a-1/2+O(e)} | C | t07_source_crossrefs.py C85 | OAI78 L15015, L15021 | CITED |
| 1340 | old-eq:3.5 | uniform slot bound Z^{kappa z+eps} under beta_* <= (1+kappa)/2, used only for slot characters (only zero-free use in lem:plain: L12833-12898) | C | t07_source_crossrefs.py C86 | OAI78 L12873, L12833, L12896 | CITED |
| 1342 | kappa = 3/4 | kappa=2(7/8)-1=3/4 in [3/4,1]; (1+kappa)/2 = 7/8 (OpenAI's theorem); OAI's kappa = 3/4+2Delta > 3/4 under beta_*>7/8 | I | t06_final_inequality.py C55 | exact/sympy | PASS |
| 1350 | lem:sextic-large-sieve | arbitrary coefficients c_n independent of the row; bound (KD)^eps{K+D+(KD)^{2/3}} | C | t07_source_crossrefs.py C87 | OAI78 L4708, L4717 | CITED |
| 1351 | only row count of Part I; numerator U^{delta/2} | G78: Part I uses only M_r (sextic large sieve); G78 6: the size U^{delta/2} of the numerator on the reflected line | C | t07_source_crossrefs.py C106 | G78 L338, L581 | CITED |
| 1354 | lem:marked | conditions r+2z<=m-c_1, 2r+8z<=3m-c_2; M_u, Q_u with psi_u = nu chi(u)^eps on n and p; bound Z^{m+eps}; either sign | C | t07_source_crossrefs.py C88 | OAI78 L9254, L9233, L9227, L9242, L9259, L9264 | CITED |
| 1362 | eq:cm-inverse | sum_{N(u)<<Z^m}\|M_u Q_u\|^2 << Z^{m+eps} for co-oriented coefficients | A | t09_untested.py U13 | - | UNTESTED: item 4 of L1493; class identities tested (C23-C25), sanity numbers without slots (C103). |
| 1363 | eq:cm-inverse (sanity, no slots) | sum_{N(u)<=H}\|M_u\|^2/H with D=400, H=D^1.1=728, 119 columns, 2646 rows | A | t05_rows_moments.py C103 | D=400, H=728 | INFO |
| 1369 | companion | mu(n)gamma_-1(n) = chi_n(-1)G(n)^-1 abar(n)gamma2(n), 165 primes + 60 composites | E | t02_gauss.py C23 | 165 primes + 60 composites | PASS |
| 1372 | co-oriented classes | mu abar = abar^2 gamma2 (gamma1/G) and mu alpha = conj(abar^2 gamma2 gamma1/G) | E | t02_gauss.py C24 | 225 sf c | PASS |
| 1379 | descent never sees alpha | the descent, masks, punctures and seminorm propagation never see alpha | C | t09_untested.py U14 | - | UNTESTED: proof inspection of OAI78 L9209-12477; not redone. |
| 1383 | anti-oriented | mu alpha = gamma2 (gamma1/G): alpha*abar*gamma2 = plain gamma2 | E | t02_gauss.py C25 | 225 sf c | PASS |
| 1387 | lem:inverse-amplification | input is lem:marked with no prime slots | C | t07_source_crossrefs.py C89 | OAI78 L12364 | CITED |
| 1388 | sixth powers | chi_n(u a^6) = chi_n(u) for (a,n)=1 (1737 triples) | E | t02_gauss.py C26 | 1737 triples | PASS |
| 1389 | (MS) as a special case | (MS) [theta in (0,1/10]] = eq:cm-inverse with eps=1, nu=1, no slots, r=1, m=1+theta, c_1=theta, c_2=1; m <= 11/10 | I | t05_rows_moments.py C50 | theta in (0,1/10] | PASS |
| 1393 | lem:plain | 3/4<=kappa<=1; z=0 no restriction; n1+n2+6 kappa z<=M; beta_*<=(1+kappa)/2; M=m+q; uniform in moving moduli, masks, frozen labels | C | t07_source_crossrefs.py C90 | OAI78 L12532, L12561, L12565, L12592, L12577 | CITED |
| 1402 | CM plain fourth moment | sum_k \|S_k(n1)S_k(n2)Q\|^2 << Z^{M+eps} with alpha(l)^e on the plain variables | A | t09_untested.py U15 | - | UNTESTED: item 5 of L1494; ingredients tested: dual profile (C43), conductor change only at lambda (C44), vanishing lattice main terms (C47), kappa = 3/4 (C55, C86). |
| 1411 | dual profile | L(s)/L(1-s) = eps(3Q)^{1/2-s}(2pi)^{2s-1}Gamma(3/2-s)/Gamma(s+1/2), so MW#(s)=MW(1-s)Gamma(s+1/2)/Gamma(3/2-s); both quotients ~\|t\|^{2 sigma-1} | E | t04_growth_fe.py C43 | sigma in [-1/2,2], t<=1e6 | PASS |
| 1413 | conductor change only at lambda | conductors of eta=chi_.(4d), abar*eta and alpha*eta agree away from lambda | E | t04_growth_fe.py C44 | 7 d | PASS |
| 1420 | lattice main terms vanish | int alpha^e f(\|z\|^2) = 0 for e=1,2,-1,3; coset sums sum_{z in a+mO} alpha(z)^e f(N z/X) have no volume term | E | t05_rows_moments.py C47 | 4 cosets, e=1,2,-1, X<=3e5 | PASS |
| 1431 | E(d) two forms | OAI78 eq:common-high-exponent: first form = second form with C_0 = -1/48 | E | t06_final_inequality.py C56 | sympy | PASS |
| 1431 | eq:common-high-exponent | E(d) = C_0 + 2delta/3 + q/6 - h(1-R) + (d-h)(R+delta/2-17/50), C_0 = -1/48 | C | t07_source_crossrefs.py C91 | OAI78 L15732, L15734 | CITED |
| 1437 | numerator U^{delta/2} | FE on Re w = 1-a: conductor power U^{(1-2(1-a))/2} = U^{delta/2} | I | t06_final_inequality.py C57 | sympy | PASS |
| 1439 | no CM gamma factor in E(d) | the four sources of E(d) do not involve the gamma factor of a CM L-function | A | t09_untested.py U16 | - | UNTESTED: argument only (the gamma factor enters through fixed powers of the height); E(d) itself checked (C56, C91). |
| 1444 | endpoint certificate | -E(h) >= 49/440640 for 0<=delta<=5/6, 0<=x<=1/2 (q=x delta, R=R_*): derivation, OAI78 identity, positivity, grid | I | t06_final_inequality.py C58 | sympy + 15351-point exact grid + 161k float grid | PASS |
| 1444 | lem:balanced-endpoint | -E_* >= 49/440640 for 0<=delta<=5/6, 0<=x<=1/2, q = x delta, R = R_* | C | t07_source_crossrefs.py C92 | OAI78 L15994, L15995, L16000 | CITED |
| 1450 | eq:adaptive-endpoint | 1-h/4 = 51/64 (h=13/16); E_actual-Delta <= -49/440640-(1-h/4)Delta+eps_total < 0 for 0<Delta<=1/24, eps_total<49/440640 | I | t06_final_inequality.py C59 | Delta in (0,1/24] | PASS |
| 1450 | eq:adaptive-endpoint | E_actual(h) - Delta <= -49/440640 - (1-h/4)Delta + eps_total = ... - (51/64)Delta | C | t07_source_crossrefs.py C93 | OAI78 L16098, L16099, L16101, L16104 | CITED |
| 1452 | h Delta/4 from kappa | kappa = 3/4+2Delta; capacity comparison raises the row exponent by at most Delta/4; coefficient of R is d | C | t07_source_crossrefs.py C94 | OAI78 L15926, L15960, L16090 | CITED |
| 1452 | kappa=3/4 removes h Delta/4 | capacity loss 2q(1-2m)(2/9 - 1/(9/2+12 Delta)) <= Delta/4 is the only Delta-dependent row loss; the baseline 2/9 = 1/(6*3/4) is kappa = 3/4; dropping h Delta/4 lowers the bound | I | t06_final_inequality.py C60 | sympy | PASS |
| 1454 | bin ceiling delta <= 5/6 | a <= beta_* <= 11/12 implies delta = 2a-1 <= 5/6 | I | t06_final_inequality.py C53 | exact | PASS |
| 1455 | beta_* <= 11/12 for F_CM | the balanced probe with the sextic row count gives beta_* <= 11/12 for F_CM | A | t09_untested.py U17 | - | UNTESTED: item 6 of L1495; not proved in the notes (Part I treats one curve at a time). |
| 1460 | contradiction logic | for every Delta in (0,1/24]: omega=eps<Delta and sigma = 49/440640+(51/64)Delta-eps > 0; kappa=3/4 gives larger sigma | I | t06_final_inequality.py C61 | Delta in {1e-6..1/24} | PASS |
| 1470 | table of OAI78 line ranges | eight consecutive parts 1-1648,...,14985-16677 start at section/part boundaries and contain the cited lemmas | C | t07_source_crossrefs.py C95 | boundaries [1648, 3313, 4687, 6807, 9209, 12477, 14985] are \section/\part lines: True | CITED |
| 1484 | six statements; eight parts | the status list has 6 items and the table 8 consecutive parts | E | t07_source_crossrefs.py C105 | guide text | PASS |
| 1488 | L2003-2279 | cross-prime phase bookkeeping lies in OAI78 L2003-2279 | C | t07_source_crossrefs.py C96 | supplementary-law periodicity at L2003, eq:reflection-cross-prime-phase at L2271 | CITED |
