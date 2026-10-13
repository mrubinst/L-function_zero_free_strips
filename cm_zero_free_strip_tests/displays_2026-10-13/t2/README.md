# t2: tests of the displays (7.1)-(9.4) of the CM guide (2026-10-13)

Guide: `sources/cm_zero_free_strip.tex` (equation numbers from `sources/cm_zero_free_strip.aux`
of 2026-10-13). Inventory: `sources/inventory/inventory_t2.json`, 48 displays in Sections 7-9 (Poisson summation,
the Moebius-Gauss and paired identities, Lemma b2-paired, Proposition 3 (b2-reduction), the completed sums and the second
theta derivative, the quadratic average and the completed estimate). One row per display is in `rows.json`.

```sh
./run_all.sh          # new tests t2a-t2d, then reruns of the existing tests cited in rows.json; 311 s on 2026-10-13
./run_all.sh --quick  # skips the end-to-end reflection test (part B t5, 130 s)
```
Outputs: `outputs/*.txt` (new tests), `outputs/reruns/*.txt` (existing tests, unchanged scripts run read-only with
`python3 -B`), `outputs/run_all.txt` (summary). Every new script prints `RESULT|id|PASS/FAIL|...` and falsifier lines, and
exits 1 if a check fails or a falsifier does not fail (`cmlib.report`/`cmlib.falsifier`; `CM_TESTS_FORCE_FAIL=1`
forces a failure, checked). `run_all.sh` exits 1 if any script exits nonzero or prints a failure line.

## Result (run of 2026-10-13: overall status 0)

| kind | count | outcome |
|---|---|---|
| E exact identity | 30 | all tested and PASS; each has a falsifier that fails, except eq:d045 (equality part), eq:d050, eq:d054 |
| I explicit inequality | 3 | all tested and PASS |
| A asymptotic/analytic bound | 4 | eq:b2-hyp and eq:completed: size checks (ratios at most 0.17 and 0.30); eq:d043, eq:d052: bookkeeping only, not testable by computation |
| D definition | 9 | consistency checked (the tests implement them; examples evaluated) |
| C cited | 2 | eq:theta-coeff = [QF, Cor. 24(2)] verbatim (and B38 PASS); eq:qls = [OAI, Lemma 6.5] / [GL, Thm 1.1] |
| R restatement | 1 | eq:d056 restates eq:completed (9.2) |

**33 of 33 E/I displays tested, 33 PASS. No failure. No display found wrong.**

Runtime: 311 s in total (new tests 22 s: t2a 1 s, t2b 17 s, t2c 1 s, t2d 3 s; reruns 289 s, of which part B t2 73 s and t5 130 s).

### New tests (written here)

- `t2a_gauss_identities.py`: eq:b2-g1g2, eq:b2-crt (both directions), eq:b2-two, eq:d040, eq:d041 (each previously only
  implicit in the class-form tests or "covered" by scripts outside tests_cm), eq:d042 with 6 characters xi, eq:d039 (all
  432 x 432 coefficients), eq:d037, eq:b2-forward/eq:b2-reverse and eq:paired as literal sums over 432 characters,
  eq:b2-twisted / eq:d038 with |a_xi| = 1, eq:d044, eq:d046 with g != 1. The character table of frak G
  (`t2lib.py`) is built from an explicit basis (orders 12, 6, 6), independently of the eigenvector table of part B.
- `t2b_sxi_endtoend.py`: eq:b2-Sxi end to end in its own (g, e0, h, z1, z2) form (before the change of variables), with
  the literal xi-sum: M_D = Z + sum_xi c_xi S_xi to relative error 7e-17, 3e-16, 4e-16 at (D, H) = (120, 170), (120, 40),
  (300, 360); falsifiers (zeta coefficients, h e0 for h e0^5, G(t) for G(t^-1), N(g) for N(g)/N(e0), mu(e0) dropped) fail
  by 1e-3 to 1.3. M_D agrees to all printed digits with the independent b2 script (final form eq:b2-output), which was also
  run at the new size D = 300, H = 360 (rel. err 1.3e-14). Also eq:d045 (Z equals the h = 0 terms, 4e-15; size check of
  Z << H||W||^2: ratio 0.08-0.095 for D = 30..5000) and eq:d047 numerically.
- `t2c_change_of_variables.py`: eq:d048 as a bijection, enumerated exactly in Z[omega] (45150 tuples (g, e0, v0, h) with
  (v0, h) = 1 onto 45150 triples (r0, f, k)), with the stated inverse; eq:d049, eq:d050 exactly on every tuple; eq:d047
  in sympy. Previously only the identities on random data (B31, B32) were tested, not the bijection.
- `t2d_sizecheck_hyp.py`: size check of the hypothesis eq:b2-hyp, E computed literally from eq:energy in the ranges
  eq:b2-ranges (theta = 1/10): E/(XF max|W|^2) <= 0.17 over 24 cases, decreasing in D (0.16, 0.13, 0.12 at B = F = 1
  for D = 500, 1500, 4000). Sanity only.

### Existing tests rerun (all PASS on 2026-10-13)

Part B `t2_poisson.py` (B09, B10, B14, B16, B17, B20), `t3_gauss.py` (B22, B23, B26-B28, B30-B33), `t4_theta.py`
(B34, B35, B38, B47, B48, B50, B51, B57, B58), `t5_reflection_j2.py` (B53), `t6_sieve.py` (B62, B64, B66);
`proofs_2026-10-12/b2/b2_reduction_check.py` (D = 120 with H = 170, 40; and D = 300, H = 360), `b2_normalization_check.py`,
`b3/b3_symbolic_checks.py`, `b4/b4_checks.py` (C1, C2); `phase2_2026-10-12/r1/test_classgroup.py` (C2, C5-C7),
`test_step45_bookkeeping.py` (S1-S4); `r2/t2_grouping_cusp_data.py`, `r2/t3_gl_family.py`;
`section10_descent/test_sec10_descent.py` (S10.11, S10.11a, S10.11b, S10.18).

### Observations (none is an error)

1. eq:d039: "all these coefficients have absolute value at most one" holds; the maximum is 1/2.
2. eq:qls: [OAI, Lemma 6.5] has N(n) asymptotic to U, the guide has N(n) <= U as in [GL, Thm 1.1]; the dyadic splitting costs
   a log, absorbed in (KU)^eps. The symbol is chi_k(n)^3 = (n/k)_2 in both.
3. eq:local: the j = 0 entry is the active branch only (part B observation 3); it is literally [OAI, (6.5)].
4. eq:theta-coeff is classified C (it is [QF, Corollary 24(2)], eq. (d0), with sqrt(-3) = lambda); it is also tested against the
   appendix tables (B38).
5. eq:b2-hyp, eq:d043, eq:d052, eq:completed have unspecified constants; the size checks above are sanity checks only.
   eq:b2-hyp follows from (10.25)-(10.26); eq:d043 and eq:d052 are proved in Proposition 3; eq:completed in Section 9.1.

## Table (same content as rows.json)

| number | label | kind | what was done | script (relative to tests_cm/) | result |
|---|---|---|---|---|---|
| 7.1 | `eq:d032` | D | Definitions of $e$ and $\widehat\Phi$: $e$ trivial on $\mathcal O$, $\mathcal O$ self-dual; Gaussian gives $(2/\sqrt3)e^{-4\pi t/3}$ | part_B_lines_0550_0997/t2_poisson.py B09, B10 | definition (B10 PASS, err 9e-22) |
| 7.2 | `eq:d033` | E | 1200 random triples $(n_1,n_2,u)$, including $u$ sharing primes with $g$; dropping the indicator fails | part_B_lines_0550_0997/t2_poisson.py B14 | PASS (max err 2e-16) |
| 7.3 | `eq:d034` | D | Implemented by brute force modulo $z_1z_2$; $\|\gamma(\psi)\|=1$ for 46 pairs; change of generator checked | part_B_lines_0550_0997/t2_poisson.py B16; proofs_2026-10-12/b2/b2_normalization_check.py | definition (B16 PASS, err 9e-16) |
| 7.4 | `eq:poisson-exact` | E | Lattice sum against dual sum, Gaussian $\Phi$, 10 cases (8 nondegenerate, $g$ with 1-2 primes); 4 falsifiers fail | part_B_lines_0550_0997/t2_poisson.py B17 | PASS (max rel err 5e-15) |
| 7.5 | `eq:d035` | I | $N(e_0)N(z_1z_2)/H\le b^2D^2/(HN(g))$ symbolically and on $10^5$ random sets; falsifier fails | part_B_lines_0550_0997/t2_poisson.py B20 | PASS (exact) |
| 7.6 | `eq:d036` | D | Vectorized $\gamma_j$ equals brute force; $\gamma_2$ is cubic (Euler criterion); $G$ equals OAI's closed form | displays_2026-10-13/t2/t2a_gauss_identities.py helper; part_B_lines_0550_0997/t3_gauss.py B22; phase2_2026-10-12/r1/test_classgroup.py C2 | definition (PASS, err 2e-14) |
| 7.7 | `eq:gauss` | E | 150 squarefree $n$ (53 composite) with brute-force-checked Gauss sums; Moebius-sign and $\alpha$ falsifiers fail | part_B_lines_0550_0997/t3_gauss.py B23 | PASS (max err 2e-15) |
| 7.8 | `eq:d037` | D | $c_\xi$ from its definition, independent table of 432 characters: $\sum_\xi c_\xi\xi(t)=G(t^{-1})$, $\sum c_\xi=1$, $\|c_\xi\|\le1/2$ | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d037; part_B_lines_0550_0997/t3_gauss.py B26 | definition (PASS, err 1e-14) |
| 7.9 | `eq:paired` | E | Literal sum over all 432 characters, 73 coprime pairs; again with an independent table (60 pairs); zeta coefficients fail | part_B_lines_0550_0997/t3_gauss.py B27; displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-forward | PASS (max err 2e-14) |
| 7.10 | `eq:energy` | D | Definition; implemented literally in the computed averages of Section 10 tests and in t2d | section10_descent/test_sec10_descent.py S10.11b, S10.18; displays_2026-10-13/t2/t2d_sizecheck_hyp.py | definition |
| 7.11 | `eq:d038` | E | 120 coprime squarefree pairs times 6 characters, plus 71 and 300 pairs in B30, S10.11a; $\chi_n(v_0)^2$ fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d038; part_B_lines_0550_0997/t3_gauss.py B30 | PASS (max err 2e-14) |
| 7.12 | `eq:initial` | I | Support bounds on prime data and sympy: $N(k)\le CD^2/(HN(r_0)^2)$, $\Sigma=XF=D/B$; one-power-of-$B$ falsifier fails | part_B_lines_0550_0997/t3_gauss.py B32; phase2_2026-10-12/r1/test_step45_bookkeeping.py S1 | PASS (exact) |
| 7.13 | `eq:b2-g1g2` | E | New: 150 squarefree $n$ (53 composite), brute-force-checked Gauss sums; $\alpha$ and Moebius-sign falsifiers fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-g1g2 | PASS (max err 2e-15) |
| 7.14 | `eq:b2-quotient` | E | All 144 pairs of classes mod 4, $G$ and $\mathcal R$ from squarefree representatives; dropping $\chi_x(-1)$ fails | phase2_2026-10-12/r1/test_classgroup.py C5 | PASS (max err 4e-15) |
| 7.15 | `eq:b2-remove` | I | Both sides computed from (7.10), 3 cases; the exact identity (10.14) behind it holds to 3e-15 | section10_descent/test_sec10_descent.py S10.11b, S10.11 | PASS (ratio at most 0.31) |
| 7.16 | `eq:d039` | D | All $432^2$ coefficients $c_\xi(\xi_1)$, $\ell_{\xi_1}(\xi)$: expansions reproduce; $c_\xi(\mathbf1)=c_\xi$; moduli at most 1/2 | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d039 | definition; claims PASS (err 3e-14) |
| 7.17 | `eq:b2-forward` | E | Literal sum over the 432 characters xi of an independent table: 60 pairs times 4 characters xi_1; the zeta coefficients fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-forward; proofs_2026-10-12/b4/b4_checks.py C1 | PASS (max err 2e-14) |
| 7.18 | `eq:b2-reverse` | E | Literal sum over the 432 characters xi_1 with l_xi1(xi): 60 pairs times 4 characters xi; beta without conj(alpha) (zeta form) fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-reverse; proofs_2026-10-12/b4/b4_checks.py C2 | PASS (max err 1.3e-14) |
| 7.19 | `eq:b2-twisted` | E | 120 coprime squarefree pairs times 6 characters: twisted multiplicativity and $\|a_\xi\|=1$; $\chi_b(a)^2$ fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-twisted; part_B_lines_0550_0997/t3_gauss.py B30 | PASS (max err 2e-14) |
| 7.20 | `eq:b2-crt` | E | New: brute-force Gauss sums mod $z_1z_2$, 83 coprime pairs (38 with $\mathcal R=-1$), both directions; dropping $\mathcal R$ fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-crt | PASS (max err 2e-15) |
| 7.21 | `eq:b2-two` | E | New: both identities for 150 squarefree $n$; $G$ for $\overline G$, and unconjugated $a^\zeta$, fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:b2-two | PASS (max err 2e-15) |
| 7.22 | `eq:d040` | E | New: both equalities, 83 coprime pairs; $G(z_1z_2^{-1})$ on the right, and dropping $\chi_{z_1}(-1)$, fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d040 | PASS (max err 2e-15) |
| 7.23 | `eq:d041` | E | New: both identities for 149 squarefree $n\ne1$; the two conjugation variants fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d041 | PASS (max err 2e-15) |
| 7.24 | `eq:d042` | E | New: 83 coprime pairs times 6 characters $\xi$; $(\xi G)(z_2z_1^{-1})$ on the right fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d042; proofs_2026-10-12/b2/b2_reduction_check.py | PASS (max err 1e-14) |
| 7.25 | `eq:b2-hyp` | A | Hypothesis; it follows from (10.25), (10.26). Size check: $E/(XF\max\|W\|^2)\le0.17$, 24 cases, $D\le4000$ | displays_2026-10-13/t2/t2d_sizecheck_hyp.py | size check: ratio at most 0.17 |
| 7.26 | `eq:b2-ranges` | D | Ranges; the Step 4 scales satisfy them at 3777 support points and after the swaps $j\mid r_0$ | phase2_2026-10-12/r1/test_step45_bookkeeping.py S1-S3 | definition |
| 7.27 | `eq:d043` | A | Conclusion of Proposition 3; only the bookkeeping $\sum(HN(r_0)/D^2)(D/B)^2\ll H$ is checked | part_B_lines_0550_0997/t3_gauss.py B33; phase2_2026-10-12/r1/test_step45_bookkeeping.py S4 | not testable by computation; proved in Proposition 3 |
| 7.28 | `eq:d044` | E | 72 triples $(g,z_1,z_2)$, $g$ with 0-2 primes; the unconjugated variant fails | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d044; part_B_lines_0550_0997/t3_gauss.py B28 | PASS (max err 4e-16) |
| 7.29 | `eq:d045` | E | $Z$ equals the $h=0$ Poisson terms and enters the end-to-end identity; size check $Z/(H\max\|W\|^2)\le0.095$ | displays_2026-10-13/t2/t2b_sxi_endtoend.py eq:d045, eq:d045-bound | PASS (rel err 4e-15) |
| 7.30 | `eq:d046` | E | 72 triples (69 with $g\ne1$), literal sum over 432 characters; $z_1=z_2=1$ gives $\sum c_\xi=1$; two falsifiers fail | displays_2026-10-13/t2/t2a_gauss_identities.py eq:d046 | PASS (max err 1e-14) |
| 7.31 | `eq:d047` | E | sympy with positive norms and complex values of $W$; numerically on every term of the end-to-end sum | displays_2026-10-13/t2/t2c_change_of_variables.py eq:d047; displays_2026-10-13/t2/t2b_sxi_endtoend.py eq:d047 | holds (symbolic) |
| 7.32 | `eq:b2-Sxi` | E | New end to end: $M_D=Z+\sum_\xi c_\xi S_\xi$ in this form, $D=120$ ($H=40,170$), $D=300$; five falsifiers fail | displays_2026-10-13/t2/t2b_sxi_endtoend.py eq:b2-Sxi | PASS (max rel err 4e-16) |
| 7.33 | `eq:d048` | E | New: bijection enumerated exactly in $\mathbb Z[\omega]$, 45150 tuples onto 45150 triples, inverse checked; three falsifiers fail | displays_2026-10-13/t2/t2c_change_of_variables.py eq:d048 | PASS (exact) |
| 7.34 | `eq:d049` | E | Exact in $\mathbb Z[\omega]$ on all 45150 tuples of the bijection; $e_0$ to the first power fails | displays_2026-10-13/t2/t2c_change_of_variables.py eq:d049; part_B_lines_0550_0997/t3_gauss.py B31, B32 | PASS (exact) |
| 7.35 | `eq:d050` | E | Exact rationals on all 45150 tuples of the bijection | displays_2026-10-13/t2/t2c_change_of_variables.py eq:d050 | PASS (exact) |
| 7.36 | `eq:b2-output` | E | End to end: $M_D=Z+\sum_\xi c_\xi S_\xi$ in this final form, $D=120$ ($H=40,170$), $D=300$; four falsifiers fail | proofs_2026-10-12/b2/b2_reduction_check.py (D 120, 300) | PASS (max rel err 1e-14) |
| 7.37 | `eq:d051` | D | Scales; at 3777 support points $N(k)\le\mathcal H$, $\varrho\in[1,4)$, $\varkappa\le C_1$, support in $[a_0/4,b_0]$ | phase2_2026-10-12/r1/test_step45_bookkeeping.py S1, S2 | definition |
| 7.38 | `eq:d052` | A | Only the dyadic bookkeeping is checked: exact sums give ratio at most 1 for every dyadic $B\le2D$ | phase2_2026-10-12/r1/test_step45_bookkeeping.py S4 | not testable by computation; proved in Proposition 3, Step 5 |
| 8.1 | `eq:T` | D | Definition; $b=1$ term is $X^{-1/2}$ times the column sum; sums finite; termwise theta form with $-1/81$ | part_B_lines_0550_0997/t4_theta.py B34, B35, B47 | definition |
| 8.2 | `eq:theta-coeff` | C | Identical to [QF, Corollary 24(2)] with $\lambda=\sqrt{-3}$; also checked against the appendix tables, 210 indices | part_B_lines_0550_0997/t4_theta.py B38 | checked against source |
| 8.3 | `eq:d053` | E | sympy with $ad-bc=1$ and a numeric $\mathrm{SL}_2(\mathcal O)$ example; $g$ for $g^{-1}$ fails | part_B_lines_0550_0997/t4_theta.py B48 | PASS (max err 4e-16) |
| 8.4 | `eq:chain` | E | sympy, generic $F$ of degree (4,2,2), no lower derivatives survive; also in b3; two falsifiers fail | part_B_lines_0550_0997/t4_theta.py B50, B51; proofs_2026-10-12/b3/b3_symbolic_checks.py | holds (symbolic) |
| 8.5 | `eq:weight` | E | End-to-end single-cusp reflection with two derivatives, 41 cusp configurations; one-derivative kernel and missing 1/27 fail | part_B_lines_0550_0997/t5_reflection_j2.py B53 | PASS (max rel err 2e-8) |
| 8.6 | `eq:local` | E | Reflected local sums, 7 primes, $j=0,\dots,5$, 1854 cases; grouped with actual cusp data; equals [OAI, (6.5)] | part_B_lines_0550_0997/t4_theta.py B57; phase2_2026-10-12/r2/t2_grouping_cusp_data.py | PASS (max err 7e-14) |
| 8.7 | `eq:d054` | E | Euler criterion, 10 primes times 40 residues | part_B_lines_0550_0997/t4_theta.py B58 | PASS (max err 3e-15) |
| 9.1 | `eq:qls` | C | Matches [OAI, Lemma 6.5] ($\chi_k(n)^3=(n/k)_2$) and [GL, Theorem 1.1]; family hypotheses and a sanity ratio 0.25 checked | part_B_lines_0550_0997/t6_sieve.py B62; phase2_2026-10-12/r2/t3_gl_family.py | checked against source |
| 9.2 | `eq:completed` | A | Sanity ratio with the CM $T$ at most 0.30 ($X=1400$, $H\le600$, 3 values of $f$); proved in Section 9.1 | part_B_lines_0550_0997/t6_sieve.py B64 | size check: ratio at most 0.30 |
| 9.3 | `eq:d055` | E | 400 random $(n,v,f)$ including non-coprime cases, and the $T$ identity in 8 cases; $\chi_n(fv)^4$ fails | part_B_lines_0550_0997/t6_sieve.py B66 | PASS (max err 4e-15) |
| 9.4 | `eq:d056` | R | Precise form of (9.2) (Proposition 4); quantifiers agree with [OAI, Proposition 5.2] | - | see eq:completed |
