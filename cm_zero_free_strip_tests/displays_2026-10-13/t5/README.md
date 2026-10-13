# t5: numbered displays (A.1)-(A.41) of the CM guide (Appendix A, the CM transformation formula)

Agent t5, 2026-10-13. Inventory: `sources/inventory/inventory_t5.json` (41 displays: the statement of
Proposition 8 (prop:b3-reflection), the theta coefficients, Lemmas 12-13 (uniformity, decay) and the six-step proof).
Equation numbers are those of the 2026-10-13 build (`sources/cm_zero_free_strip.aux`; all 41 checked
against the inventory). Every display was read in context in `cm-eleven-twelfths-guide.tex` (lines 3239-3776).

## Counts

| kind | count | displays |
|---|---|---|
| E (exact identity) | 27 | A.2, A.10, A.12, A.14-A.18, A.20, A.22-A.32, A.34-A.36, A.38-A.41 |
| I (explicit inequality) | 4 | A.3, A.5, A.9, A.33 |
| A (asymptotic, unspecified constants) | 3 | A.7, A.8, A.37 |
| D (definition) | 7 | A.1, A.4, A.6, A.11, A.13, A.19, A.21 |
| C (quoted) / R (restatement) | 0 | (displays taken from [QF] or [OAI1112], e.g. A.5, A.10, A.18, A.25, were tested directly and are counted as E/I; the source match is noted in rows.json) |

**All 31 E/I displays were tested and pass**: 23 by new tests in this folder (several of them also by existing
tests), and 8 (A.2, A.17, A.22, A.26, A.34, A.38, A.39, A.41) only by existing tests that really check that display
and were rerun today. The 3 A displays have size checks; the 7 D displays were
checked for consistency (each is implemented literally in the tests). Every new test has falsifiers, and every
falsifier fails as required.

**Failures: none. Displays suspected wrong in the guide: none.**

## New tests (this folder)

Run `python3 -B run_all.py` (about 80 s wall; `--reruns` also reruns the 11 existing scripts below in parallel, writing
only into `outputs/reruns/`). Each script exits 1 on any failed check or on a falsifier that does not fail;
`CM_TESTS_FORCE_FAIL=1` forces a failure (self-test, verified). Helpers: `t5lib.py` (imports `zw.py` and
`fastgauss.py` read-only from `part_B_lines_0550_0997`; `t5a` and `t5c` also use `cmlib.gamma_j_fast` from there).

| script | time | displays | method |
|---|---|---|---|
| `t5a_finite_identities.py` | 13 s | A.1, A.2, A.3, A.12, A.18, A.23, A.27, A.33 (= A.5 termwise), A.40 | exact / machine-precision finite identities in Z[omega]; brute-force gamma_2; literal C_(pi,j); sympy |
| `t5b_theta_pointwise.py` | 43 s | A.4/A.29, A.14, A.15, A.16, A.25, A.28, A.30, A.31, A.35 | the actual theta function of [QF, Thm 23] (tables equal part D's literal transcription to 1.3e-13), own cusp data per [QF (5.18)] for 61 cusps of all types and for all translates of 8 twists (L = 72, P empty or one prime of norm 7 with j = 0, 1, 2, 4) |
| `t5c_mellin.py` | 20 s | A.10, A.20 (with A.19), A.24, A.32, A.36 | mpmath (30 digits); Mellin integrals of the truncated twisted theta series by quadrature against the Dirichlet series built from alpha, gamma_2, Psi |
| `t5d_decay.py` | 11 s | A.7 (size check), A.9 | V# of a compactly supported W = y^(-1/2) (log y log(2/y))^18 by contour integrals on four lines |

## Existing tests rerun today (all PASS; outputs in `outputs/reruns/`)

| script | wall time | displays it checks |
|---|---|---|
| `proofs_2026-10-12/b3/b3_symbolic_checks.py` | 6 s | A.8 (6), A.22 (1), A.23 (7), A.26 and A.27 at z = 0 (part 2), A.34 (4), A.36 (5), A.37 (8), A.38 (3) |
| `proofs_2026-10-12/b3/t_b3_prop_endtoend.py` | 118 s | A.2 and A.39 (P empty, before grouping): rel err 1.1e-12 |
| `phase2_2026-10-12/r2/t1_prop_grouped_endtoend.py` | 186 s | A.2 grouped with A.41: rel err <= 2.5e-12 for j = 0..5, 7 falsifiers fail |
| `phase2_2026-10-12/r2/t2_grouping_cusp_data.py` | 17 s | A.15 (g1 in Gamma(3)), A.17 (exact), A.18 with actual cusp data, kappa0 = [OAI1112 (A.10)] |
| `part_B_lines_0550_0997/t4_theta.py` | 17 s | B44 (A.30), B45 (A.22), B46 (A.23), B48 (A.25), B49 (A.26), B51 (A.27), B57 (A.18), B59 (A.6) |
| `part_B_lines_0550_0997/t5_reflection_j2.py` | 133 s | single-cusp reflection (B53, B55): structure of A.2 |
| `part_D_lines_1510_end/t02_chain_rule.py`, `t03_mellin.py`, `t04_kernels.py` | 43, 30, 48 s | chain rule (A.27/A.28 at z = 0), mode Mellin, kernel growth D04.3 (A.8), V# size D04.4 (A.7) |
| `part_D_lines_1510_end/t07_theta_appendix.py` | 54 s | D07.2-D07.3 (A.4, A.29), D07.4 (A.5, sharp) |
| `part_D_lines_1510_end/t08_reflection_j.py` | 209 s | D08.2: single-cusp reflection with j = 2 (A.2) |

**Total runtime of the tests: new tests 82 s wall (74 s CPU); reruns 861 s summed over the scripts (the times above,
first run, outputs in `outputs/reruns_first/`). The final `python3 -B run_all.py --reruns` (outputs in `outputs/reruns/`
and `outputs/run_all.txt`) took 285 s wall with all 11 reruns in parallel. ALL PASS.**

## Remarks (no errors found)

1. Lemma 13 (A.7): the small-x exponent 1/4 is valid but not sharp. The fitted slope of log|V#| at x -> 0 is 4/3
   (the first pole of the CM kernel is at t = -4/3), for j = 0, 1, 2. This matches the proof, which only needs the
   line Re t = -1/4.
2. Lemma 13 is stated for W in C_c^infinity, but its proof uses only ||W||_(C^J), J = ceil(4A)+j+2. The size check uses a
   C^17 bump, which is admissible for A <= 3, j <= 2. The observed large-x decay is about x^(-4.4), which is consistent
   with "4A derivatives buy x^(-A)". The implied constants are large (up to 6e15 for A = 3, j = 2), as expected for a
   narrow bump with large ||W||_(C^J).
3. The bound (A.5) is attained (|d_0| = 27 3^(m/6) |b| on row 4 of the tau table: D07.4, and again here up to
   N(9l) <= 30000). The representation l = u lambda^m n b^3 is unique (no duplicates among about 75000 coefficients).
4. Cusp data: the own implementation of [QF (5.18)] (delta, H, sigma, kappa(g1) = (c1/a1)_3) agrees with
   kappa = (b1/d1)_3 for every cusp used. For each of the 12 matrices H, H gamma_sigma^(-1) = S U with S in SL_2(Z)
   and U in Gamma(3), kappa(U) = 1, so theta(HP) = theta(gamma_sigma P) by (T1).
5. In `t_b3_prop_endtoend.py` the control conj(alpha(c))^2 coincides with the claim for twists '1' and 'A', whose cusps
   are real. It fails only for twist 'B', at rel err 1.0. The new tests t5b.A35 and r2 T1 have non-real c = c0 pi, and
   there the conj(alpha(c))^2 variant fails.
6. Two bugs in the first drafts of the new tests were fixed before the final run. Neither was a problem in the
   guide. (a) The coefficient check of (A.14) used ce(x z_h/lambda) instead of ce(x z_h/lambda^3). The pointwise
   check of the same display passed throughout. (b) The mpmath quadrature in (A.10) was inaccurate at s = -0.3; the
   substitution v = u^15 fixed it.

## Rows

`rows.json` has one row per display (label, printed number, kind, what was done, script, result):

| no. | label | kind | what was done | script (relative to tests_cm/) | result |
|---|---|---|---|---|---|
| (A.1) | `eq:d096` | D | Definition of $\Psi$; checked complete multiplicativity and vanishing at even $n$: 5 characters, 1500 random pairs | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A1 | definition; consistency PASS (exact) |
| (A.2) | `eq:reflection` | E | End to end: $T(X;\Psi)$ from (8.1) vs the grouped sum, 72 groups, $P=\{\pi\}$, $j=0,\dots,5$; also $P$ empty | phase2_2026-10-12/r2/t1_prop_grouped_endtoend.py T1; proofs_2026-10-12/b3/t_b3_prop_endtoend.py | PASS (max rel err 2.5e-12; 7 falsifiers fail) |
| (A.3) | `eq:d097` | I | Factors of $C'$ in (A.41): $|\hat\phi(h_0)|\le1$ for all $h_0$ (5 twists), $|\varpi|=1$, $|\kappa_0|=1$ | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A3 | PASS (max $|\hat\phi|$ = 0.083; $|\varpi|=1$ to 3e-14) |
| (A.4) | `eq:d098` | D | Same as [QF, Thm 23 (T3)]; these coefficients reproduce $\overline\theta$ at the three cusps numerically | part_D_lines_1510_end/t07_theta_appendix.py D07.2, D07.3; displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A15 | definition; checked against source (rel err 4e-15) |
| (A.5) | `eq:d099` | I | Bound checked for every nonzero $d_\sigma(\ell)$ with $N(9\ell)\le30000$; sharp; same as [QF Cor. 24] | part_D_lines_1510_end/t07_theta_appendix.py D07.4; displays_2026-10-13/t5/t5a_finite_identities.py t5a.A33 | PASS (0 violations; equality attained) |
| (A.6) | `eq:d100` | D | Definition: (A.1) with $j_\pi=1$ for $\pi\mid k_0$; row primes have $j=1$ (B59), hence are active | part_B_lines_0550_0997/t4_theta.py B59; displays_2026-10-13/t5/t5a_finite_identities.py t5a.A40 | definition |
| (A.7) | `eq:d101` | A | Size check with compact $W\in C^{17}$: ratio to $\min(x^{1/4},x^{-A})$ bounded, $A\le3$, $j\le2$; slope 4/3 at 0 | displays_2026-10-13/t5/t5d_decay.py t5d.A7; part_D_lines_1510_end/t04_kernels.py D04.4 | size check: bounded on $[10^{-8},10^8]$; proved in Lemma 13 |
| (A.8) | `eq:b3-stirling` | A | Fitted growth exponents of $\|R(\eta+iu)\|$ equal $4\eta$ for five $\eta$; the ratio stays bounded on the sampled grid | R(\eta+iu) | $ equal $4\eta$ for five $\eta$; finite sup of the ratio | proofs_2026-10-12/b3/b3_symbolic_checks.py (6); part_D_lines_1510_end/t04_kernels.py D04.3 | size check: exponents $4\eta$ within 2e-5; proved by Stirling |
| (A.9) | `eq:d102` | I | Inequality at 33 values $x\in[10^{-8},10^8]$, $j\le2$, $\eta=-1/4,1,2$; left side from the line $\Re t=0$ | displays_2026-10-13/t5/t5d_decay.py t5d.A9 | PASS (0 violations) |
| (A.10) | `eq:b3-bessel` | E | mpmath, 30 digits, 8 cases including $s=-0.3$ and complex $s$; same as [QF (5.21)] | displays_2026-10-13/t5/t5c_mellin.py t5c.A10 | PASS (max rel err 1e-29) |
| (A.11) | `eq:d103` | D | Definition of $\Theta_\Psi$; implemented in the pointwise tests of (A.14) and (A.31) | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A14, t5b.A31 | definition |
| (A.12) | `eq:b3-coeff` | E | 3424 indices: $d_0$ from the QF table times the twist vs the right side with brute-force $\gamma_2$ | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A12 | PASS (max abs err 2e-13; 2 falsifiers fail) |
| (A.13) | `eq:d104` | D | Definitions of $\hat\phi$ and $C_{\pi,j}$; implemented literally with $L=72$ | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A40; displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A14 | definition |
| (A.14) | `eq:b3-translates` | E | 8 twists, $L=72$: coefficient identity for all $N(x)\le1500$ and pointwise with the actual $\theta$; $\sum c_F=0$ | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A14 | PASS (max err 9e-14) |
| (A.15) | `eq:d105` | E | $g_1\in\Gamma(3)$ in exact arithmetic (528 translates); $\theta(HP)=\theta(\gamma_\sigma P)$ numerically, 12 matrices $H$ | phase2_2026-10-12/r2/t2_grouping_cusp_data.py (b); displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A15 | PASS (exact for $g_1$; $\theta$ to rel err 4e-15) |
| (A.16) | `eq:b3-automorphy` | E | Pointwise at 61 cusps of all types, $\theta$ from [QF Thm 23], own cusp data | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A16 | PASS (max rel err 4e-14) |
| (A.17) | `eq:b3-crt` | E | Exact in $\mathbb Z[\omega]$ for 46 groups (528 translates), $L=72$, $M=\lambda^{12}72^4$ | phase2_2026-10-12/r2/t2_grouping_cusp_data.py (d) | PASS (exact) |
| (A.18) | `eq:b3-local` | E | Literal $C_{\pi,j}$: 4 primes (one inert), all $j$, all $x$; also with actual cusp data | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A18; phase2_2026-10-12/r2/t2_grouping_cusp_data.py (e) | PASS (max err 6e-14) |
| (A.19) | `eq:d106` | D | Definition of $\mathcal T(w)$; used as the product of the two series in the test of (A.20) | displays_2026-10-13/t5/t5c_mellin.py t5c.A20 | definition |
| (A.20) | `eq:b3-inversion` | E | $T(X;\Psi)$ from (8.1) vs the $w$-integral on $\Re w=1.5,2$, $X=60,300$ | displays_2026-10-13/t5/t5c_mellin.py t5c.A20 | PASS (max rel err 8e-13; 2 falsifiers fail) |
| (A.21) | `eq:d107` | D | Definitions of $f$ and $\mathcal J$; implemented in the tests of (A.24) and (A.35) | displays_2026-10-13/t5/t5c_mellin.py t5c.A24; displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A35 | definition |
| (A.22) | `eq:b3-mode` | E | mpmath, 40 digits, 9 pairs $(\ell,s)$ | proofs_2026-10-12/b3/b3_symbolic_checks.py (1); part_B_lines_0550_0997/t4_theta.py B45 | PASS (max rel err 1e-40) |
| (A.23) | `eq:d108` | E | 720 pairs $(n,b)$, three values of $s$; $\alpha(\lambda^{-3})=i$ | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A23; proofs_2026-10-12/b3/b3_symbolic_checks.py (7) | PASS (max rel err 4e-15) |
| (A.24) | `eq:b3-J` | E | Quadrature of $\int f(v)v^{2s-1}dv$ from the twisted theta series vs the Dirichlet series, $s=2.5$, $3+1.2i$ | displays_2026-10-13/t5/t5c_mellin.py t5c.A24 | PASS (max rel err 2e-14; 3 falsifiers fail) |
| (A.25) | `eq:d109` | E | Formula vs the action of $g^{-1}$ at 61 cusps; sympy with $ad-bc=1$ | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A25; part_B_lines_0550_0997/t4_theta.py B48 | PASS (max err 3e-15) |
| (A.26) | `eq:d110` | E | sympy: the six first and second derivatives | proofs_2026-10-12/b3/b3_symbolic_checks.py part 2; part_B_lines_0550_0997/t4_theta.py B49 | PASS (exact) |
| (A.27) | `eq:d111` | E | sympy, generic polynomial $F$, general $z$; at $z=0$ only $Y=Y'=z'$ survives | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A27; part_B_lines_0550_0997/t4_theta.py B51 | PASS (exact) |
| (A.28) | `eq:b3-cusp` | E | Pointwise at 61 cusps: second $\bar z$-derivative from the expansion at infinity vs the cusp side | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A28 | PASS (max rel err 2e-14; 3 falsifiers fail) |
| (A.29) | `eq:d112` | E | Expansions from $d_\sigma$ equal $\overline{\theta(\gamma_\sigma P)}$; constant term only for $\sigma=0$ | part_D_lines_1510_end/t07_theta_appendix.py D07.2, D07.3; displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A15, t5b.A16 | PASS (rel err 4e-15) |
| (A.30) | `eq:b3-G` | E | Used in the pointwise test of (A.28); the variant with $\overline{\alpha(\ell)}^2$ fails | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A28; part_B_lines_0550_0997/t4_theta.py B44 | PASS (max rel err 2e-14) |
| (A.31) | `eq:d113` | E | $f(v)$ from the twisted series vs the sum over 36 to 252 translates, 8 twists | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A31 | PASS (max rel err 5e-14) |
| (A.32) | `eq:b3-dual` | E | Single mode (mpmath) and truncated series at three cusps by quadrature | displays_2026-10-13/t5/t5c_mellin.py t5c.A32 | PASS (max rel err 2e-14) |
| (A.33) | `eq:d114` | I | Termwise bound for all $d_\sigma(\ell)$, $N(9\ell)\le30000$; partial sums for $s'=1.6,2,3$ | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A33 | PASS (0 violations) |
| (A.34) | `eq:d115` | E | Substitution by mpmath quadrature at two $s$; $c^{-4}N(c)^2=\overline{\alpha(c)}^4$ | proofs_2026-10-12/b3/b3_symbolic_checks.py (4) | PASS (max rel err 1e-41) |
| (A.35) | `eq:b3-fe` | E | $\mathcal J(-3+0.5i)$ by quadrature in $v$ vs the right side, 8 twists | displays_2026-10-13/t5/t5b_theta_pointwise.py t5b.A35 | PASS (max rel err 2e-10; 2 falsifiers fail) |
| (A.36) | `eq:d116` | E | sympy: (A.24) at $s=w+1/2$ | displays_2026-10-13/t5/t5c_mellin.py t5c.A36; proofs_2026-10-12/b3/b3_symbolic_checks.py (5) | holds (symbolic) |
| (A.37) | `eq:d117` | A | Fitted exponents 6.0, 4.0, 2.4 at $\sigma_2=-1,-0.5,-0.1$ equal $2-4\sigma_2$ | proofs_2026-10-12/b3/b3_symbolic_checks.py (8) | size check: exponents to 1e-6; proved by Stirling |
| (A.38) | `eq:d118` | E | sympy, and numerically at three values of $t$ | proofs_2026-10-12/b3/b3_symbolic_checks.py (3) | holds (symbolic) |
| (A.39) | `eq:d119` | E | End to end, $P$ empty, three twists, all cusp types: (8.1) vs this formula at $\Re t=0$ | proofs_2026-10-12/b3/t_b3_prop_endtoend.py | PASS (max rel err 1.1e-12; 3 falsifiers fail) |
| (A.40) | `eq:d120` | E | $C_{\pi,0}(0)=1-1/N(\pi)$ and $C_{\pi,j}(0)=0$ for $j\ne0$, 4 primes | displays_2026-10-13/t5/t5a_finite_identities.py t5a.A40 | PASS (max err 4e-16) |
| (A.41) | `eq:b3-C` | E | Grouped end to end with this $C$, $j=0,\dots,5$; $\kappa_0$ agrees with [OAI1112 (A.10)] | phase2_2026-10-12/r2/t1_prop_grouped_endtoend.py T1; phase2_2026-10-12/r2/t2_grouping_cusp_data.py (c) | PASS (max rel err 2.5e-12) |
