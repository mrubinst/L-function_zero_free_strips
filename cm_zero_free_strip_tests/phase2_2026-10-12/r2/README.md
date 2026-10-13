# r2 tests: Sections 8, 9, 10 and Appendix A of the CM guide

Run from this folder (Python 3 with numpy, scipy, sympy; PARI/GP `gp` on the PATH for T2):

```sh
python3 -B t2_grouping_cusp_data.py      # ~10 s   exact arithmetic; output outputs/t2_grouping_cusp_data.txt
python3 -B t1_prop_grouped_endtoend.py   # ~90 s   output outputs/t1_prop_grouped_endtoend.txt
python3 -B t3_gl_family.py               # ~2 s    output outputs/t3_gl_family.txt
```
`r2lib.py` holds the exact Z[omega] arithmetic, the residue symbols (norms factored with PARI/GP) and the cusp data of
[OAI1112, App. A.2] with L = 72 (a verbatim copy of the helper part of `t2_grouping_cusp_data.py`). T1 also imports
`fastgauss.py` and `zw.py` read-only from `tests_cm/proofs_2026-10-12/b3/`.

**Last run (2026-10-12): T1, T2, T3 PASS; every falsifier fails.** The existing tests that cover these sections were
rerun unchanged (outputs in `outputs/reruns/`, nothing written into `tests_cm/`): all PASS.

## New tests

| id | claim (guide location) | method | result |
|---|---|---|---|
| T1 | Proposition `prop:b3-reflection` (A.1) **in its grouped form**, with one prime in P: local factors `B_{pi,j}` (8.5), constant (A.16) `C = -conj(alpha(c))^4 phihat(h0) conj(kappa0) prod(1-1/N) prod chi(sigma_pi)^-2 varpi`, psi of (A.7), cusp data from Step 2 with L = 72 | Psi = conj((lambda/.)_3) chi_pi^j, N(pi) = 7, j = 0..5 (j = 0 with active and inactive branch, j = 4 Ramanujan case); left side = (8.1) summed directly (log-normal V); right side = sum over the 72 groups (h0 in O/72, A) of C times the dual series with the CM kernel (8.4) | **PASS**: rel. error 2.5e-12, 8.6e-13, 3.1e-13, 2.6e-13, 6.8e-13, 1.4e-13 for j = 0..5; groups cover sigma in {0,+,-} and N(c) in {1,4,7,28} (c = c0*pi not real). Falsifiers (max rel. error): alpha(c)^4 unconjugated 1.38; conj(alpha(c))^2 1.96; zeta prefactor -i/81 1.41; kappa0 unconjugated 3.49; chi^-j instead of chi^(-j-2) 4.02; one-derivative kernel 5.52; conj(alpha(b))^3 on the left 1.52 |
| T2 | Appendix A, Step 2 (a)-(e) and the Step 6 grouping, with **L = 72, M = lambda^12 72^4** (E12, E13); multiplier formula [OAI1112 (A.10)]; (A.7) = [OAI (A.11)]; (A.8) = [OAI (A.12)] with the actual kappa(g1) and delta' | exact arithmetic in Z[omega] (big integers); h_pi lifted = 0 mod M^2; delta' by the CRT congruences of [OAI App. A.2]; kappa(g1) = (b1/d1)_3 with N(d1) ~ 1e24 factored by PARI; 46 groups (40 with one split prime of norm 7, 3 with the inert prime -5, 3 with two primes of norms 7 and 13), 528 translates | **PASS**: c constant; c0 divides 8 lambda^2; a = sigma_pi h_pi mod pi; g in SL_2(O), g1 in Gamma(3); sigma and delta' mod lambda^3 c0 constant on each group; kappa(g1) prod chi_pi(a)^-2 constant and equal to (A.10) in all three cases; (A.7) exact in O; grouped local identity for every j-vector at 48 points x per vector, max abs. error 1.1e-14. 13 cusp types (sigma, v_lambda(c), N(c0)) occur, covering every possible pair (sigma, v_lambda(c)): sigma in {0,+,-} for v = 0, sigma = +- for v = 1, sigma = 0 for v >= 2. Falsifiers: kappa instead of conj(kappa) 7.55; chi^-j 5.76; psi dropped 12.6 |
| T3 | Section 9.1, the large sieve (E7): chi_(k) = (./k)_2 kappa_lambda^{e_k} is trivial on units and equals (n/k)_2 on primary n; reciprocity factor depends only on classes mod 4 (GL Def. 1 (2)); k1 = k2 mod 4 gives e equal and product induced by (./k1k2)_2, primitive (GL Def. 1 (3)); (-2n'/k)_2 = (-2/k)_2 (n'/k)_2 | 421 squarefree primary k (N <= 1500), 5712 coprime pairs (all 144 class pairs mod 4), 298 pairs in one class | **PASS**, 0 violations. Falsifiers: classes mod 2 (2748 violations); (./k)_2 without kappa^e is not trivial on units (636) |

## Existing tests reused (rerun 2026-10-12, unchanged; outputs in `outputs/reruns/`)

| claim | test | rerun |
|---|---|---|
| (8.3) chain rule, six derivatives, (A.10) mode Mellin, (A.11) constant, (A.12)-(A.14), c^-4 N(c)^2 = conj(alpha(c))^4, (A.2) growth 4 eta, PL exponent 2-4 sigma_2, Step 6 constant -1/81 | `tests_cm/proofs_2026-10-12/b3/b3_symbolic_checks.py` | PASS |
| Proposition 7 with P empty, end to end (three twists, all cusp types) | `.../b3/t_b3_prop_endtoend.py` | PASS, 1.1e-12 |
| single-cusp reflection with two derivatives, 41 cusp configurations incl. non-real c (C = 4 - omega, -2 - 3 omega) | `tests_cm/part_B_lines_0550_0997/t5_reflection_j2.py` (B53, B55) | PASS, 2.3e-8 |
| (10.5)-(10.13), (10.15)-(10.20), (10.23), (10.1), (10.11), citations of Section 10 | `tests_cm/section10_descent/test_sec10_descent.py` | ALL PASS |
| Lemmas 7, 8 (first and second Poisson identities, end to end), regrouping tau(r) 1_{f'\|k'}, (b4-long) with conj(alpha)^6, Hecke-family facts, scale identities | `.../b4/b4_checks.py` (C1-C8, C10, C11) | 14 PASS |
| support and bound abs(d_sigma) <= 27 3^{m/6} abs(b) (sharp) | `tests_cm/part_D_lines_1510_end/t07_theta_appendix.py` (D07.4) | not rerun (unchanged) |
| local factors B_{pi,j} as finite sums | part B `t4_theta.py` (B57) | not rerun; superseded by T2 (e) with actual cusp data |

Not testable: the analytic inequalities themselves (9.1), (9.2), (10.4), Proposition 4 (transfer) and the induction,
beyond the size checks already in `tests_cm`; Lemma 11 (uniformity) for two rows in one class modulo M^2 (N(M^2) is
about 1.5e41); T2 checks the constancy within a group, which is the part of Lemma 11 used in the grouping.
